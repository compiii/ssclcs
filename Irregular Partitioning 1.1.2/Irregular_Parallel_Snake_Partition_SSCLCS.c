#include "Irregular_Parallel_Snake_SSCLCS.h"
 
void init(){ 
    maxBlock = pow(maxProcess,2);
    maxDiag = (2*maxProcess) - 1;
    sizePartA = (m / maxProcess);
    sizePartB = (n / maxProcess);

    depth = l;
    maxk = l-1;
    maxkrecv = 0;
    id_Tmp = 0;
    
    arrayData = malloc(maxBlock * sizeof(BlockData));
    arrayDataSplit = malloc(((24*k_frag) + (maxBlock)) * sizeof(BlockData));//maxBlock*(k_frag+maxProcess)
    if((arrayData == NULL)){
        printf("Echec d'allocation de la mémoire pour la table de données (arrayData et tabBlock) !\n");
        exit(1);
    }

    // if(rank == 0){
    //     printEnterInfo();
    // } 
} 
 
void constructBlockData(){
    int n1 = 0;
    int m1 = 0;
    int level = 1;
    int level2 = 1;
    int id = 0;
    int idSnake = 0;
    int invDiag = maxProcess, tmpInit;
    BlockData blockData;
    for (int i = 0; i < maxProcess; ++i){
        for (int j = 0; j < maxProcess; ++j){
            blockData.diag = level;
            blockData.id = id;
            blockData.i = m1;
            blockData.j = n1;
            blockData.k = maxk;
            n1 = sizeBlock(n1, sizePartB, n);
            blockData.dim = computeDimBlock(depth, n1, sizeBlock(m1, sizePartA, m));

            /** Begin **/
            if (j == 0){
                idSnake = idSnake + i;
                tmpInit = idSnake;
            }else{
                if(level <= maxProcess)
                    idSnake = idSnake + level;
                    else{
                        if(level == maxProcess + 1)
                            idSnake = idSnake + maxProcess;
                            else{
                                invDiag = invDiag - 1;
                                idSnake = idSnake + invDiag;
                            }
                    }
            }
            blockData.id = idSnake;
            blockData.rank = (idSnake < maxProcess) ? idSnake : (idSnake % maxProcess);
            level++;

            arrayData[id] = blockData;
            id++;
            // if(rank == 0) printBlockData(blockData);
        }
        idSnake = tmpInit;
        invDiag = maxProcess;
        /** End **/

        m1 = sizeBlock(m1, sizePartA, m);
        n1 = 0;
        level2++;
        level = level2;
    }
}

void loadBlock_Array(){
    snake_rank = 0;
    is_begin = 0;
    indice_Itteration = 0;
    diagonal_level = 1;
    int level = 0;
    BlockData* bd1 = NULL;
    BlockData* bd2 = NULL;
    bd1 = malloc(4 * sizeof(BlockData));
    bd2 = malloc(4 * sizeof(BlockData));
    for (int i = 0; i < maxBlock; i++){
        if((arrayData[i].diag == 1) || (arrayData[i].diag == (2*maxProcess - 1)) || (((arrayData[i].diag == 2)  || (arrayData[i].diag == (2*maxProcess - 2))) && maxProcess > 2)){
            if((arrayData[i].diag == (2*maxProcess - 2)) && (maxProcess > 2))
                diagonal_level = arrayData[i].diag + 2;

            if(arrayData[i-1].id != 0)
                level = ((arrayData[i].id - arrayData[i-1].id) == 2) ? 0 : 2;

            bd1 = fragmentationBlock(arrayData[i], level);
            if((is_begin == 0) && (maxProcess > 2) && (k_frag > 1)) {
                diagonal_level = 1;
                is_begin = 1;
            }
            int ii = 0;
            while (ii < 4){
                if((k_frag == 2) && ((arrayData[i].diag == 1))){
                    for(int j = 0; j < 3; j++){
                        if(bd1[j-1].id != 0)
                            level = ((bd1[j].id - bd1[j-1].id) == 2) ? 0 : 2;
                        bd2 = fragmentationBlock(bd1[j], level);
                        
                        for(int jj = 0; jj < 4; jj++){
                            bd2[jj].id = indice_Itteration;
                            bd2[jj].rank = snake_rank;
                            arrayDataSplit[indice_Itteration] = bd2[jj];
                            snake_rank = (snake_rank < (maxProcess - 1)) ? snake_rank + 1 : 0;
                            indice_Itteration++;
                            // if(rank == 0) printBlockData(arrayDataSplit[indice_Itteration - 1]);
                        }
                        ii++;
                    }
                    is_begin = 2;
                }
                bd1[ii].diag += ((ii == 3) && (arrayData[i].diag == 1) && (maxProcess > 2) && (k_frag > 1)) ? 3 : 0;
                bd1[ii].id = indice_Itteration;
                bd1[ii].rank = snake_rank;
                arrayDataSplit[indice_Itteration] = bd1[ii];
                indice_Itteration++;
                snake_rank = (snake_rank < (maxProcess - 1)) ? snake_rank + 1 : 0;
                // if(rank == 0) printBlockData(arrayDataSplit[indice_Itteration - 1]);
                ii++;

                if((k_frag == 2) && ((arrayData[i].diag == (2*maxProcess - 1)))){
                    is_begin = 4;
                    for(int j = 1; j < 4; j++){                        
                        diagonal_level = diagonal_level + ((ii == 1) ? 2 : 0);
                        level = (ii == 2) ? 2 : 0;
                        bd2 = fragmentationBlock(bd1[j], level);
                        for(int jj = 0; jj < 4; jj++){
                            bd2[jj].id = indice_Itteration;
                            bd2[jj].rank = snake_rank;
                            arrayDataSplit[indice_Itteration] = bd2[jj];
                            snake_rank = (snake_rank < (maxProcess - 1)) ? snake_rank + 1 : 0;
                            indice_Itteration++;
                            // if(rank == 0) printBlockData(arrayDataSplit[indice_Itteration - 1]);
                        }
                        ii++;
                    }
                }
            }
        }else{
            arrayData[i].diag = arrayData[i].diag + 3  + (((k_frag > 1) && (maxProcess > 2)) ? 3 : 0);
            arrayData[i].id = indice_Itteration;
            arrayData[i].rank = snake_rank;
            arrayDataSplit[indice_Itteration] = arrayData[i];
            indice_Itteration++;
            snake_rank = (snake_rank < (maxProcess - 1)) ? snake_rank + 1 : 0;
            // if(rank == 0) printBlockData(arrayDataSplit[indice_Itteration - 1]);
        }
    }
    tri_insertion(indice_Itteration);
    // for (int i = 0; i < indice_Itteration; i++){
    //    if(rank == 0) printBlockData(arrayDataSplit[i]);
    // }
    free(arrayData);    
}

void tri_insertion(int sizeTab){
    int i, j, idNew = 0;
    BlockData en_cours;
 
    for (i = 1; i < sizeTab; i++) {
        en_cours = arrayDataSplit[i];
        for (j = i; j > 0 && arrayDataSplit[j - 1].diag > en_cours.diag; j--) {
            arrayDataSplit[j - 1].id = j;
            arrayDataSplit[j] = arrayDataSplit[j - 1];
        }
        en_cours.id = j;
        arrayDataSplit[j] = en_cours;
    }
}

int sizeBlock(int sizeString, int sizePartString, int mn){
    return ((sizeString+sizePartString) >= (maxProcess*sizePartString)) ? 
            ((sizeString + sizePartString) + (mn%maxProcess)) : (sizeString + sizePartString);
}

Dim computeDimBlock(int plan, int n1, int m1){
    Dim dim;
    dim.nbrLine = (m1 >= m) ? sizePartA + (m % maxProcess) : sizePartA;
    dim.nbrColumn = (n1 >= n) ? sizePartB + (n % maxProcess) : sizePartB;
    dim.plan = plan;
    return dim;
}


BlockData* fragmentationBlock(BlockData block, int level){
    static BlockData* bd = NULL;
    bd = malloc(4 * sizeof(BlockData));
    int sizeB_i = block.dim.nbrLine / 2;
    int sizeB_j = block.dim.nbrColumn / 2;
    int m_i = sizeB_i;
    int n_j = sizeB_j;
    int ii = 0;
    int iid = 0;
    int line = block.i;
    int column = block.j;
    diagonal_level = ((block.diag >= (2*maxProcess - 3)) && (block.diag != 2*maxProcess - 1)  && (is_begin != 4) && (maxProcess > 2)) ? diagonal_level - level + 3 : diagonal_level - level;
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            bd[ii].dim.nbrLine = m_i + ((j == 1) ? (block.dim.nbrLine % 2) : 0);
            bd[ii].dim.nbrColumn = n_j + ((i == 1) ? (block.dim.nbrColumn % 2) : 0);
            bd[ii].diag = diagonal_level + (((k_frag > 1) && (is_begin == 2) && (maxProcess > 2)) ? 3 : 0);
            bd[ii].i = line;
            bd[ii].j = column;
            bd[ii].k = block.k;
            bd[ii].id = iid;

            column += sizeB_j + ((j == 1) ? (block.dim.nbrColumn % 2) : 0);
            ii++;
            iid++;
            diagonal_level++;
        }
        diagonal_level--;
        column = block.j;
        line += sizeB_i + ((i == 1) ? (block.dim.nbrLine % 2) : 0);
    }

    return bd;
}

int profondeur(int numproc){
    int rest;
    rest = l % maxProcess;
    int result, difference;
    if(rest == 0){
        result = (l / maxProcess);
    }else{
        difference = maxProcess - rest;
        if (numproc > (difference - 1)){
            result = (l / maxProcess) + 1;
        }else{
            result = (l / maxProcess);
        }
    }
    return result;
}


////***************************** Functions for output Data **********************//
void printEnterInfo(){
    printf("********* Valeurs des chaines : X=%ld, Y=%ld et C=%ld\n", n,m,l);
    printf("max processor: %d", maxProcess);
    printf("\nmax diagonal: %d", maxDiag);
    printf("\nrank processor: %d", rank);
    printf("\nsize string X: %ld", n);
    printf("\nsize string Y: %ld", m);
    printf("\nsize string C for each processor: %d", depth);
    printf("\nmax block: %d", maxBlock);
    printf("\nmax block per processor: %d\n", maxBlockPerProcess);
}
void printBlockData(BlockData bd){
    printf("\nid block : %d", bd.id);
    printf("\n\tRun processor : %d , depth=%d, maxk=%d", bd.rank,depth, maxk);
    printf("\n\tcoordinates : ");
    printf("\n\t\ti -> %d", bd.i);
    printf("\n\t\tj -> %d", bd.j);
    printf("\n\t\tk -> %d", bd.k);
    printf("\n\tdimensions : ");
    printf("\n\t\tnbLine -> %d", bd.dim.nbrLine);
    printf("\n\t\tnbColumn -> %d", bd.dim.nbrColumn);
    printf("\n\tdiagonale : %d", bd.diag);
    printf("\n\n***************************************************\n \n \n");
}

void printTabBlockData(){   
    int i = 0;
    while(i < (pow(maxProcess,2))){
        printBlockData(arrayData[i]);
        i++;
    }    
}

void easyPrint(){    
    char filname[10];
    filname[0]='\0';
    strcat(filname,"P");
    gcvt((double) rank,5,convertstr);
    strcat(filname,convertstr);
    FILE *fip = NULL;
    fip = fopen(filname,"a+");
}

