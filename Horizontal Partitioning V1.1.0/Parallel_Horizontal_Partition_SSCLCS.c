#include "Parallel_Horizontal_SSCLCS.h"
 
void init(){ 
    maxBlock = pow(maxProcess,2);
    maxDiag = (2*maxProcess) - 1;
    sizePartA = (m / maxProcess);
    sizePartB = (n / maxProcess);

    depth = l;
    maxk = l-1;
    maxkrecv = 0;
    
    arrayData = malloc(maxBlock * sizeof(BlockData));
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
    BlockData blockData;
        for (int i = 0; i < maxProcess; ++i){
            for (int j = 0; j < maxProcess; ++j){
                blockData.diag = level;
                blockData.id = id;
                level++;
                blockData.i = m1;
                blockData.j = n1;

                n1 = sizeBlock(n1, sizePartB, n);
                blockData.dim = computeDimBlock(depth, n1, sizeBlock(m1, sizePartA, m));

                blockData.k = maxk;
                blockData.rank = i;
                arrayData[id] = blockData;
                id++;
            }
            m1 = sizeBlock(m1, sizePartA, m);
            n1 = 0;
            level2++;
            level = level2;
        }
        // for (int i = 0; i < maxBlock; i++){        
        //     if(rank == 0) printBlockData(arrayData[i]);
        // }
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


////***************************** Fonctions pour la sortie des données **********************//
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

