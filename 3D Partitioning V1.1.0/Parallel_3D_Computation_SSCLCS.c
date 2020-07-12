#include "Parallel_3D_SSCLCS.h" 

/******* Travaux préliminaires (Calcul des valeurs de sigma) ********/
void allocations(){
    X = (char*) malloc(sizeof(char)*m);
    Y = (char*) malloc(sizeof(char)*n);
    listparts = (ListPartitions*) malloc(sizeof(ListPartitions)*depth);
    for (int i = 0; i < depth; i++){
        listparts[i].partition = (char*) malloc(sizeof(char)*contraintLength);
        if (listparts[i].partition == NULL) {
            printf("error: Partition allocation failed");
            exit(0);
        }
    }
    listpartsTmp = (ListPartitions*) malloc(sizeof(ListPartitions)*l);
    for (int i = 0; i < l; i++){
        listpartsTmp[i].partition = (char*) malloc(sizeof(char)*contraintLength);
        if (listpartsTmp[i].partition == NULL) {
            printf("error: Partition allocation failed");
            exit(0);
        }
    }
    partitionTmp = (char*) malloc(sizeof(char)*contraintLength);
    prevMatchTab_X = (long*) malloc(sizeof(long)*m*s);
    prevMatchTab_Y = (long*) malloc(sizeof(long)*n*s);
    f1 = (long*) calloc((m)*(n)*(depth) + 2, (sizeof(long)));
    f2 = (long*) calloc((m)*(n)*(depth) + 2, (sizeof(long)));
    f1Tmp = (long*) calloc((m)*(n) + 2, (sizeof(long)));
    f2Tmp = (long*) calloc((m)*(n) + 2, (sizeof(long)));
    listStartPos_X = (ListstartPostTab_X*) malloc(sizeof(ListstartPostTab_X)*depth);
    listStartPos_Y = (ListstartPostTab_Y*) malloc(sizeof(ListstartPostTab_Y)*depth);
    checkAllocationMemory1D(prevMatchTab_X,"prevMatch");
    checkAllocationMemory1D(prevMatchTab_Y,"prevMatch");

    if ((f1 == NULL) || (f2 == NULL)){
        printf("Allocation DP tables has Error");
        exit(0);
    }
    
}

void getData(){
    FILE* DNAConstraint = NULL;
    FILE* DNASequence_X = NULL;
    FILE* DNASequence_Y = NULL;
    DNAConstraint = fopen("Resources/constraint.fasta","r");
    DNASequence_X = fopen("Resources/sequence_1.fasta","r");
    DNASequence_Y = fopen("Resources/sequence_2.fasta","r");
    if ((DNASequence_X == NULL) || (DNASequence_Y == NULL) || (DNAConstraint == NULL)) {
		printf("error: missing file(s) DNA resources");
		exit(0);
	}

    fgets(X, m+1, DNASequence_X);
    fgets(Y, n+1, DNASequence_Y);
    
    /** Partitions allocation **/
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < l; i++){
        fgets(listpartsTmp[i].partition, contraintLength+1, DNAConstraint); 
        fseek(DNAConstraint, contraintLength+10, SEEK_SET);
    }

    int strid = 0;
    int beg =  (rank == 0) ? maxkrecv : maxkrecv + 1;
    for(int i = beg; i <= maxk; i++){
        listparts[strid].partition = listpartsTmp[i].partition;
        strid++;
        //free(listpartsTmp[i].partition);        
    }

    fclose(DNASequence_X);
    fclose(DNASequence_Y);
    fclose(DNAConstraint);
    free(listpartsTmp);
}

void loadSequence_X(char* tmp){
    for(int i = 0; i < m; i++){
        X[i] = tmp[i+1];
    }
}
void loadSequence_Y(char* tmp){
    for(int i = 0; i < n; i++){
        Y[i] = tmp[i+1];
    }
}

void preliminaryWork(){
	Sigma = malloc(4*sizeof(char));
	if (Sigma == NULL){
		printf("Erreur allocation mémoire pour Sigma\n");
		exit(0);
	}
    
    //Communication round to send Last partition to other processors
    MPI_Type_contiguous(contraintLength, MPI_CHAR, &linePart);
	MPI_Type_commit(&linePart);
    if (rank != 0){
        MPI_Recv(&partitionTmp[0],1,linePart,rank-1,rank-1,MPI_COMM_WORLD,&status);
    }
    if (rank != (maxProcess - 1)){
        MPI_Send(&listparts[depth-1].partition[0],1,linePart,rank + 1,rank,MPI_COMM_WORLD);
    }
    MPI_Type_free(&linePart);

    prevMatchTab_X = prevMatch(X);
    prevMatchTab_Y = prevMatch(Y);

    //Compute StartPosTable for all partitions of string X and Y
    for (long i = 0; i < depth; i++){
        listStartPos_X[i].startPost = startPosCom(X,listparts[i].partition,prevMatchTab_X);
        listStartPos_Y[i].startPost = startPosCom(Y,listparts[i].partition,prevMatchTab_Y);
    }
    free(prevMatchTab_X);
    free(prevMatchTab_Y);

    f1[(n-1)*n*depth + ((n-1)*depth) + (depth - 1) + 1] = 0;
    f2[(n-1)*n*depth + ((n-1)*depth) + (depth - 1) + 1] = 0;
    f1[i_j_k2(n-1, n-1, depth - 1) + 2] = LONG_MIN;
    f2[i_j_k2(n-1, n-1, depth - 1) + 2] = LONG_MIN;

    f1Tmp[(n-1)*n + ((n-1)) + 1] = 0;
    f2Tmp[(n-1)*n + ((n-1)) + 1] = 0;
    f1Tmp[(n-1)*n + ((n-1)) + 2] = LONG_MIN;
    f2Tmp[(n-1)*n + ((n-1)) + 2] = LONG_MIN;
}

void blockCompute(BlockData blockData){
	int lig, col;
	Boundary bound = computeLimitBlock(blockData);
	lig = bound.lig;
	col = bound.col;
	Coord ci = blockData.coord;

	int count, blockLengths, stride, strideCol;
	count = blockData.dim.nbrColumn;
	blockLengths = 1;
	stride = depth;
	strideCol = depth*m;

	MPI_Type_vector(count, blockLengths, stride, MPI_LONG, &line);
	MPI_Type_commit(&line);

    MPI_Type_vector(count, blockLengths, 1, MPI_LONG, &vectorLin);
	MPI_Type_commit(&vectorLin);

    tmpDouble = MPI_Wtime();
    if((rank != 0)){
		    ///***************** Receive calculated table cells ****************/
        for (int i = 0; i < blockData.dim.nbrLine; i++){
            MPI_Recv(&f1Tmp[i_j_k3(ci.i+i,ci.j,0)], 1, vectorLin, rank - 1, i + blockData.id, MPI_COMM_WORLD, &status);
		    MPI_Recv(&f2Tmp[i_j_k3(ci.i+i,ci.j,0)], 1, vectorLin, rank - 1, i + m*m + blockData.id, MPI_COMM_WORLD, &status);
        }
	}

    comTime += MPI_Wtime() - tmpDouble;

	 //***************** Calcul du bloc de données ... *****************//
	tmpDouble = MPI_Wtime();
	computeBlock(lig, col, blockData);
    calculTime += MPI_Wtime() - tmpDouble;

	if((rank != maxProcess - 1)){
	    tmpDouble = MPI_Wtime();
			/****** sends calculated table cells ****************/ 
        for (int i = 0; i < blockData.dim.nbrLine; i++){
            MPI_Isend(&f1[i_j_k2(ci.i+i,ci.j,depth-1)], 1, line, rank + 1, i + blockData.id, MPI_COMM_WORLD, &request);
            MPI_Isend(&f2[i_j_k2(ci.i+i,ci.j,depth-1)], 1, line, rank + 1, i + m*m + blockData.id, MPI_COMM_WORLD, &request);
        }
        comTime = comTime +  MPI_Wtime() - tmpDouble;
	}

    MPI_Type_free(&line);
    MPI_Type_free(&vectorLin);
}

Boundary computeLimitBlock(BlockData blockData){
	Boundary bound;
	bound.lig = blockData.coord.i + blockData.dim.nbrLine - 1;
	bound.col = blockData.coord.j + blockData.dim.nbrColumn - 1;
	return bound;
}

void allBlockCompute(){
	int i = 0;
	while(i < pow(maxProcess,2)){
		blockCompute(arrayData[i]);
	    i++;
	}
}

void computeBlock(int lig, int col, BlockData blockData){
	for (long i = blockData.coord.i; i <= lig; i++){
        for (long j = blockData.coord.j; j <= col; j++){
            for (long k = 0; k < depth; k++){
                long res1 = 0;
                long res2 = 0;
                long res3 = 0;
                long value1 = 0;
                long value2 = 0;

                if((X[i] == Y[j]) && (Y[j] == searchPartitionInList(k)) && (k >= 0)){
                        long w;
                        if (k == 0){
                            w = (rank == 0) ? 0 : computeWeightPartition(partitionTmp, listparts[k].partition);;
                        }else{
                            w = computeWeightPartition(listparts[k-1].partition, listparts[k].partition);
                        }

                        long l_k = strlen(listparts[k].partition);
                        long startPosX_i = listStartPos_X[k].startPost[i_j(i,w,m)];
                        long startPosY_j = listStartPos_Y[k].startPost[i_j(j,w,n)];
                        long startPosX_0 = listStartPos_X[k].startPost[i_j(i,l_k,m)];
                        long startPosY_0 = listStartPos_Y[k].startPost[i_j(j,l_k,n)];

                        long f1Value;
                        long f2Value;
                        if((k == 0) && (rank != 0)){
                            f1Value = f1Tmp[i_j_k3(startPosX_0,startPosY_0,0)];
                            f2Value = f2Tmp[i_j_k3(startPosX_i,startPosY_j,0)];                            
                        }else{
                            f1Value = f1[i_j_k2(startPosX_0,startPosY_0,k-1)];
                            f2Value = f2[i_j_k2(startPosX_i,startPosY_j,k-1)];
                        }
                        f2[i_j_k2(i,j,k)] = ((f2Value + l_k - w) > (f1Value + l_k)) ? (f2Value + l_k - w) : (f1Value + l_k);
                        value1 = f2[i_j_k2(i,j,k)];
                }else{
                    if((X[i] == Y[j])){
                        res1 = f1[i_j_k2(i-1,j-1,k)] + 1;
                    }else{
                        res3 = (f1[i_j_k2(i-1,j,k)] > f1[i_j_k2(i,j-1,k)]) ? f1[i_j_k2(i-1,j,k)] 
                                            : f1[i_j_k2(i,j-1,k)];
                        f2[i_j_k2(i,j,k)] = (f2[i_j_k2(i-1,j,k)] > f2[i_j_k2(i,j-1,k)]) ? f2[i_j_k2(i-1,j,k)] 
                                            : f2[i_j_k2(i,j-1,k)];
                        value2 = f2[i_j_k2(i,j,k)];
                    }
                }

                res2 =  (value1 >= value2) ? value1 : value2;
                f1[i_j_k2(i,j,k)] = maxNumber(res1, res2, res3);
                if((rank == maxProcess - 1) && (i == m-1) && (j == n-1) && (k == depth-1)) printf("P_%d.\t f1[%ld] = %ld \n computation time = %f \n",
                                    rank, i_j_k2(m-1,n-1,depth-1), f1[i_j_k2(m-1,n-1,depth-1)], comTime);
			}
		}
	}
}

int isSamePlanScale(int k){
    return (k <= maxkrecv) ? 0 : 1;
}

void findAlphabet(){
	char *res = NULL;
	int i, l;
	int j = 0;

		for(i = 0; i < n; i++){
			res = strrchr(Sigma, *(X+i));
			if(res == NULL)
			{
				*(Sigma+j) = *(X+i);
				j++;
			}
		}

		for(i = 0; i < m; i++){
			res = strrchr(Sigma, *(Y+i));
			if(res == NULL)
			{
				*(Sigma+j) = *(Y+i);
				j++;
			}
		}
		*(Sigma+j) = '\0';
}

void findAlphabet_standard(){
	Sigma = "ACGT";
}


void freeMemory(){
    free(Sigma);
    free(X);
    free(Y);
    for (int i = 0; i < depth; i++){
        free(listparts[i].partition);
    }
    free(listparts);
    free(listStartPos_X);
    free(listStartPos_Y);
    free(arrayData);
    free(f1Tmp);
    free(f2Tmp);
    free(f2);
    free(f1);
}

//Function to compute Previous Match table
long* prevMatch(char* stringChar){
    long* prevMatchTab = NULL;
    prevMatchTab = (long*) malloc(sizeof(long)*strlen(stringChar)*s);
    long position = -1;
    long newPosition = -1;
    long i,j;
    for ( i = 0; i < s; i++){
        for (j = 0; j < strlen(stringChar); j++){
            if ((Sigma[i] == stringChar[j])){
                newPosition = j;
            }      
            prevMatchTab[(i*strlen(stringChar))+j] = position;
            position = newPosition;        
        }
        position = -1;
        newPosition = -1;     
    }
    
    return prevMatchTab;
}

void printPrevMatch(char* stringChar, long* prevM){
    printf("#### \t X = %s \t",stringChar);
    printf("#### \t Y = %s \t",Y);
    printf("#### \t C1 = %s \t",partition);
    printf("#### \t Sigma = %s \n\n",Sigma);
    printf("\t\t*************** PREVIOUS MATCH TABLE ***************\n\n");
    printf("\t");
    for (long i = 0; i < strlen(stringChar); i++){
        if(i == 0){
            printf("\t----");
        }
        printf("\t%c(%ld)",stringChar[i],i);
    }
    long k = 0;
    long w = 0;
    for (long j = 0; j < (s*strlen(stringChar)); j++){
        w++;
        if((j == 0) || (w == strlen(stringChar))){
            w = 0;
            printf("\n\n\t\t(%ld)%c",k,Sigma[k]);
            k++;
        }
        printf("\t %ld",prevM[j]);
    }
    printf("\n\n");
}

//Function to compute StartPos Table
long* startPosCom(char* stringChar, char* stringPart, long* prevM){
    long prevTabVal = -1;
    long* startPosTab = NULL;
    long mc = strlen(stringChar);
    long lp = strlen(stringPart);
    startPosTab = (long*) malloc(sizeof(long)*mc*(lp+1));
    for (long i = 0; i < mc; i++){
        for (long j = lp-1; j >= 0; j--){
            if(stringPart[lp-1] == stringChar[i]){
                if((lp-1 == j)){
                    startPosTab[j*mc+i] = i;
                }else{
                    long sj = searchIndiceInSigma(stringPart[j]);
                    long pos = (sj*mc) + i;
                    prevTabVal = prevM[pos];
                    long startPosVal = startPosTab[(j+1)*mc+i];
                    if(startPosVal != -1){
                        while((startPosVal < prevTabVal) && (prevTabVal != -1)){
                            pos--;
                            prevTabVal = prevM[pos];                         
                        }
                        startPosTab[j*mc+i] = prevTabVal;
                    }else{
                        startPosTab[j*mc+i] = -1;
                    }
                }
            }else{
                startPosTab[j*mc+i] = -1; 
            }
        }
        startPosTab[lp*mc+i] = startPosTab[i] - 1; ; 
    }
    return startPosTab;
};

//Search previous position of character of X in prevMatch table
long searchPrevMatchTab(long indiceSigma, long indiceX, long* prevM){
    long pos = (indiceX * m) + indiceSigma;
    return prevM[pos];
}

void printStartPos(char* stringChar, char* stringPart, long* startP){
    printf("\t\t*************** Start Positions TABLE ***************\n\n");
    printf("\t");
    for (long i = 0; i < strlen(stringChar); i++){
        if(i == 0){
            printf("\t----");
        }
        printf("\t%c(%ld)",stringChar[i],i);
    }
    long k = 0;
    long w = 0;
    for (long j = 0; j < (strlen(stringChar)*(strlen(stringPart)+1)); j++){
        w++;
        if((j == 0) || (w == strlen(stringChar))){
            w = 0;
            char ch = '0';
            printf("\n\n\t\t(%ld)%c",k, (j<=((strlen(stringChar)*strlen(stringPart))-1)) 
                                ? stringPart[k] : ch);
            k++;
        }
        printf("\t %ld",startP[j]);
    }
    printf("\n");     
};

//Function Search incide of one chararacter in Alphabet Sigma
long searchIndiceInSigma(char ch){
    long indice = 0;
	while(((Sigma[indice]) != ch) && (indice <= strlen(Sigma))){
		indice++;
	}
	if(indice > strlen(Sigma))
		return -1;
	return indice;	
}

//Function to compute set of weight between two successifs partitions
long computeWeightPartition(char* c1, char* c2){
    //long* W;
    long p = 0;
    for (long i = 0; i < strlen(c1); i++){
        for (long j = p; j < strlen(c2); j++){
            if(c1[i] == c2[j]){
                p++;
                //W[p] == p;
            }
        }        
    }
    return p;
}

//Function to get 1D array address for a 3D array
long i_j_k2(long i, long j, long k){
    if(((k==-1) && (i==-1) && (j>=-1)) || ((k==-1) && (i>=-1) && (j==-1)))
        return (n-1)*n*depth + ((n-1)*depth) + (depth - 1) + 1;
    if(((i < -1) || (j < -1)) || ((k >= 0) && ((i <= -1) || (j <= -1))))
        return (n-1)*n*depth + ((n-1)*depth) + (depth - 1) + 2;
        else
            return (i*(n)*(depth))+ (j*(depth)) + k;
}

long i_j_k3(long i, long j, long k){
    if(((k==-1) && (i==-1) && (j>=-1)) || ((k==-1) && (i>=-1) && (j==-1)))
        return (n-1)*n + n;
    if(((i < -1) || (j < -1)) || ((k >= 0) && ((i <= -1) || (j <= -1))))
        return (n-1)*n + n + 1;
        else
            return (i*n)+ j + k;
}

long i_j(long i, long j, long stringLength){
    return j*stringLength+i;
}

//Function to compute partition indice
char searchPartitionInList(long k){
    long length = strlen(listparts[k].partition);
    return listparts[k].partition[length - 1];
}

//Compute max number between many
long maxNumber(long nb1, long nb2, long nb3){
    long max = (nb1 >= nb2) ? nb1 : nb2;
    max = (max >= nb3) ? max : nb3;
    return max;
} 

//Functions check memory for allocation 1D, 2D and 3D (D = dimension)
void checkAllocationMemory1D(long* tab, char* text){
    if(tab == NULL){
        printf("Allocated memory failed for table: %s\n",text);
        exit(0);
    }
};
void checkAllocationMemory2D(long** tab, char* text){
    if(tab == NULL){
        printf("Allocated memory failed for table: %s\n",text);
        exit(0);
    }
};
void checkAllocationMemory3D(long*** tab, char* text){
    if(tab == NULL){
        printf("Allocated memory failed for table: %s\n",text);
        exit(0);
    }
};