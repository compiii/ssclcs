/*
  * Created by Hermann (19/01/2020), Sequential "C" algorithm to solved SSCLCS problem playground82gitup
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>


/**
 * Data structure definition, Functions prototypes and globals declarations
*/

typedef struct ListPartitions_L ListPartitions;
struct ListPartitions_L{
    char* partition;
    //ListPartitions* next;
};

typedef struct ListstartPostTab_LX ListstartPostTab_X;
struct ListstartPostTab_LX{
    long* startPost;
};

typedef struct ListstartPostTab_LY ListstartPostTab_Y;
struct ListstartPostTab_LY{
    long* startPost;
};


long m, s, l, n, r_l;
char* X = NULL;
char* Y = NULL;
char* sigma = NULL;
char* partition = NULL;
char* partition2 = NULL;
long* prevMatchTab_X = NULL;
long* prevMatchTab_Y = NULL;
short* f1 = NULL;
short* f2 = NULL;
short* f1Tmp = NULL;
short* f2Tmp = NULL;
ListPartitions* listparts = NULL;
ListstartPostTab_X* listStartPos_X = NULL;
ListstartPostTab_Y* listStartPos_Y = NULL;
double startTime, endTime, totalTime;

long* prevMatch(char* stringChar);
void printPrevMatch(char* stringChar, long* prevM);
long* startPosCom(char* stringChar, char* stringPart, long* prevM);
void printStartPos(char* stringChar, char* stringPart, long* startP);
long searchPrevMatchTab(long indiceSigma, long indiceX, long* prevM);
long searchIndiceInSigma(char ch);
long computeWeightPartition(char* c1, char* c2);
long i_j_k(long i, long j, long k);
long i_j_k2(long i, long j, long k);
long i_j(long i, long j, long stringLength);
char searchPartitionInList(long k);
long maxNumber(long nb1, long nb2, long nb3);
void getData();
void loadSequence_X(char* tmp);
void loadSequence_Y(char* tmp);


/**
 * Main function
*/
void checkAllocationMemory1D(long* tab, char* text);
void checkAllocationMemory2D(long** tab, char* text);
void checkAllocationMemory3D(long*** tab, char* text);

long main(long argc, char const *argv[]){
    long clk_tck = CLOCKS_PER_SEC;

    startTime  = endTime = totalTime = 0;
    int TAILLE_MAX; 
    int maxNumberConstraint;
    if (argc == 4){
	    TAILLE_MAX = atol(argv[1]);
        r_l = atol(argv[2]);
        maxNumberConstraint = atol(argv[3]);
        // offsetSequence = atol(argv[3]);
        // offsetConstraint = atol(argv[4]);
    }else{
        printf("Arguments incorrect !\n");
        exit(1);
    }
    m=TAILLE_MAX, s=4, l=maxNumberConstraint, n=TAILLE_MAX;
    s = 4;

    X = (char*) malloc(sizeof(char)*m);
    Y = (char*) malloc(sizeof(char)*n);
    sigma = (char*) malloc(sizeof(char)*s);
    listparts = (ListPartitions*) malloc(sizeof(ListPartitions)*l);
    for (int i = 0; i < l; i++){
        listparts[i].partition = (char*) malloc(sizeof(char)*r_l);
        if (listparts[i].partition == NULL) {
            printf("error: Partition allocation failed");
            exit(0);
        }
    }
    prevMatchTab_X = (long*) malloc(sizeof(long)*m*s);
    prevMatchTab_Y = (long*) malloc(sizeof(long)*n*s);
    f1 = (short*) calloc((m)*(n)*(l) + 2, (sizeof(short))); // The last  +1 is for case who have i, j or k > -1
    f2 = (short*) calloc((m)*(n)*(l) + 2, (sizeof(short))); // The last  +1 is for case who have i, j or k > -1
    listStartPos_X = (ListstartPostTab_X*) malloc(sizeof(ListstartPostTab_X)*l);
    listStartPos_Y = (ListstartPostTab_Y*) malloc(sizeof(ListstartPostTab_Y)*l);
    checkAllocationMemory1D(prevMatchTab_X,"prevMatch");
    checkAllocationMemory1D(prevMatchTab_Y,"prevMatch");
    checkAllocationMemory1D(f1,"first array function");
    checkAllocationMemory1D(f2,"second array function");

    startTime = (double)clock();
    getData();
    prevMatchTab_X = prevMatch(X);
    prevMatchTab_Y = prevMatch(Y);

    //Compute StartPosTable for all partitions of string X and Y
    for (long i = 0; i < l; i++){
        listStartPos_X[i].startPost = startPosCom(X,listparts[i].partition,prevMatchTab_X);
        listStartPos_Y[i].startPost = startPosCom(Y,listparts[i].partition,prevMatchTab_Y);
    }

    //First and second recursive function (f1, f2) to compute SSCLCS
    f1[(n-1)*n*l + ((n-1)*l) + (l - 1) + 1] = 0;
    f2[(n-1)*n*l + ((n-1)*l) + (l - 1) + 1] = 0;
    f1[(n-1)*n*l + ((n-1)*l) + (l - 1) + 2] = SHRT_MIN;
    f2[(n-1)*n*l + ((n-1)*l) + (l - 1) + 2] = SHRT_MIN;
    for (long i = 0; i < m; i++){
        for (long j = 0; j < n; j++){
            for (long k = 0; k < l; k++){
                long res1 = 0;
                long res2 = 0;
                long res3 = 0;
                long value1 = 0;
                long value2 = 0;
                
                if((X[i] == Y[j]) && (Y[j] == searchPartitionInList(k)) && (k >= 0)){
                    long w = (k == 0) ? 0 : computeWeightPartition(listparts[k-1].partition, listparts[k].partition);
                    long l_k = strlen(listparts[k].partition);
                    long startPosX_i = listStartPos_X[k].startPost[i_j(i,w,m)];
                    long startPosY_j = listStartPos_Y[k].startPost[i_j(j,w,n)];
                    long startPosX_0 = listStartPos_X[k].startPost[i_j(i,l_k,m)];
                    long startPosY_0 = listStartPos_Y[k].startPost[i_j(j,l_k,n)];

                    f2[i_j_k(i,j,k)] = ((f2[i_j_k(startPosX_i,startPosY_j,k-1)] + l_k - w) > 
                                        (f1[i_j_k(startPosX_0,startPosY_0,k-1)] + l_k))
                                        ? (f2[i_j_k(startPosX_i,startPosY_j,k-1)] + l_k - w) : 
                                        (f1[i_j_k(startPosX_0,startPosY_0,k-1)] + l_k);
                    value1 = f2[i_j_k(i,j,k)];    
                        // if((i==m-1) && (n==n-1) && (k==l-1)) printf("value1 = %ld", value1);
                }else{
                    if((X[i] == Y[j])){
                        res1 = f1[i_j_k(i-1,j-1,k)] + 1;
                        // if((i==m-1) && (n==n-1) && (k==l-1)) printf("res1 = %ld", res1);
                    }else{
                        res3 = (f1[i_j_k(i-1,j,k)] > f1[i_j_k(i,j-1,k)]) ? f1[i_j_k(i-1,j,k)] 
                                            : f1[i_j_k(i,j-1,k)];
                        f2[i_j_k(i,j,k)] = (f2[i_j_k(i-1,j,k)] > f2[i_j_k(i,j-1,k)]) ? f2[i_j_k(i-1,j,k)] 
                                            : f2[i_j_k(i,j-1,k)];
                        value2 = f2[i_j_k(i,j,k)];                   
                        // if((i==m-1) && (n==n-1) && (k==l-1)) printf("res3 = %ld && value2 = %ld", res3, value2);
                    }
                }
                    
                res2 =  (value1 >= value2) ? value1 : value2;
                f1[i_j_k(i,j,k)] = maxNumber(res1, res2, res3);
                // printf("2.(%ld,%ld,%ld) f1[%ld] = %ld \n",i,j,k, i_j_k(i,j,k), f1[i_j_k(i,j,k)]);
            }            
        }     
    }

    printf("f2(%ld) = %ld \n", i_j_k(m-1,n-1,l-1),f1[i_j_k(m-1,n-1,l-1)]);
    free(prevMatchTab_X);
    free(prevMatchTab_Y);
    endTime = (double)clock();
    
    totalTime = (endTime - startTime) / (double)clk_tck;

    FILE *fip = NULL;
    fip = fopen("SSCLCS_Result_Sequentiel.csv","a");
    fseek(fip, 0L, SEEK_END);
	    if (ftell(fip) == 0) 
            fprintf(fip, "Size_X;Size_Y;Size_Contraint;Contraints_Number;Total_Time;SSCLCS\n");
    fprintf(fip, "\n");
    fprintf(fip, "%d;", TAILLE_MAX);
	fprintf(fip, "%d;", TAILLE_MAX);
	fprintf(fip, "%ld;", r_l);
	fprintf(fip, "%d;", maxNumberConstraint);
	fprintf(fip, "%lf;", totalTime);
	fprintf(fip, "%ld;", f1[i_j_k(m-1,n-1,l-1)]);

    printf("Compute SUCCESS! \nComputed Time = %lf\n", totalTime);
    return 1;
}

void getData(){
    sigma = "ACGT";
    FILE* DNAConstraint = NULL;
    FILE* DNASequence_X = NULL;
    FILE* DNASequence_Y = NULL;
    DNAConstraint = fopen("Resources/constraint.fasta","r");
    DNASequence_X = fopen("Resources/sequence_1.fasta","r");
    DNASequence_Y = fopen("Resources/sequence_2.fasta","r");
    if ((DNASequence_X == NULL) || (DNAConstraint == NULL) || (DNASequence_Y == NULL)) {
		printf("error: missing file(s) DNA resources");
		exit(0);
	}
    fgets(X, m+1, DNASequence_X);
    fgets(Y, n+1, DNASequence_Y);
    
    /** Partitions allocation **/
    fseek(DNAConstraint, (0)*(r_l), SEEK_SET);
    for(int i = 0; i < l; i++){
        fgets(listparts[i].partition, r_l + 1, DNAConstraint); 
        fseek(DNAConstraint, 10 + r_l, SEEK_SET);
        // printf("c_%d = %s \t ConstraintLength = %ld\n and size = %ld\n", i, listparts[i].partition, r_l, strlen(listparts[i].partition));
    }
    // printf("\t****** X = %s \t m = %ld = %ld \n\t****** Y = %s\t n = %ld = %ld\n", X, m, strlen(X), Y, n, strlen(Y)) ;
    fclose(DNASequence_X);
    fclose(DNASequence_Y);
    fclose(DNAConstraint);
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


/**
 * Helper functions (preprocessing)
*/

//Function to compute Previous Match table
long* prevMatch(char* stringChar){
    long* prevMatchTab = NULL;
    prevMatchTab = (long*) malloc(sizeof(long)*strlen(stringChar)*s);
    long position = -1;
    long newPosition = -1;
    long i,j;
    for ( i = 0; i < s; i++){
        for (j = 0; j < strlen(stringChar); j++){
            if ((sigma[i] == stringChar[j])){
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
    printf("#### \t Sigma = %s \n\n",sigma);
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
            printf("\n\n\t\t(%ld)%c",k,sigma[k]);
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
	while(((sigma[indice]) != ch) && (indice <= strlen(sigma))){
		indice++;
	}
	if(indice > strlen(sigma))
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
    i = (i == -1) ? m : i;
    j = (j == -1) ? n : j;
    k = (k == -1) ? l : k;
	return (i*(m)*(l))+ (j*(l)) + k;
}
long i_j_k(long i, long j, long k){
    if(((k==-1) && (i==-1) && (j>=-1)) || ((k==-1) && (i>=-1) && (j==-1)))
        return (n-1)*n*l + ((n-1)*l) + (l - 1) + 1;
    if(((i < -1) || (j < -1)) || ((k >= 0) && ((i <= -1) || (j <= -1))))
        return (n-1)*n*l + ((n-1)*l) + (l - 1) + 2;
        else
            return (i*(n)*(l))+ (j*(l)) + k;
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
        exit(EXIT_FAILURE);
    }
};
void checkAllocationMemory2D(long** tab, char* text){
    if(tab == NULL){
        printf("Allocated memory failed for table: %s\n",text);
        exit(EXIT_FAILURE);
    }
};
void checkAllocationMemory3D(long*** tab, char* text){
    if(tab == NULL){
        printf("Allocated memory failed for table: %s\n",text);
        exit(EXIT_FAILURE);
    }
};