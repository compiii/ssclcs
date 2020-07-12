/**
 * @file datasets.c
 * @author Jerry Lacmou (jerrylacmou@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2019-01-11
 * 
 * @copyright Compiii Thesis (c) 2019 
 * 
 */
#include "datasets.h"
#include "config.h"
#include "clogger.h"

void initDataSetsMod(){
    maxSizeData = TAILLE_MAX;
    maxSizeConstraint = sizeMaxC;
    //maxNumberConstraint = 1;
    listFistSequenceData = createListSequenceData();
    listSecondSequenceData = createListSequenceData();
    listConstraintSequenceData = createListSequenceData();
    
    listFistSequenceData = readFileData(1);
    listSecondSequenceData = readFileData(2);
    listConstraintSequenceData = readFileData(3);
}

void createDatasets(char* datasets){
    int i, d, offset = 5000, nbLine;
    char* tmp, *sequence, *pathSequence;
    FILE* f;

    nbLine = (int) ceil(MAX_SIZE_SEQUENCE / 70.0); 
    maxSizeData = 70 * nbLine;
    
    tmp = malloc(70 * sizeof(char));
    sequence = malloc(maxSizeData * sizeof(char));
    pathSequence = malloc(1000 * sizeof(char));
    
    file = fopen(datasets, "r");
    if (file == NULL) {
			printf("error: missing file(s) %s", datasets);
			exit (EXIT_FAILURE);
	}

    fseek(file, 72 * offset, SEEK_SET);

    for(i = 0; i < nbLine; i++){
        fscanf(file, "%70s\n", tmp);
        sequence = concat(3, sequence, tmp, "\n");
    }

    pathSequence = concat(1, "../resources/datasets/dna/NC_002516.2/NC_002516.2_3/sequence_1.fasta");
    f = fopen(pathSequence, "w");
    fprintf(f, "%s", sequence);
    fclose(f);
    free(sequence);
    
    rewind(file);
    sprintf(sequence, "%s", "");
    offset = 10000;
    
    fseek(file, 72 * offset, SEEK_SET);

    for(i = 0; i < nbLine; i++){
        fscanf(file, "%70s\n", tmp);
        sequence = concat(3, sequence, tmp, "\n");
    }

    pathSequence = concat(1, "../resources/datasets/dna/NC_002516.2/NC_002516.2_3/sequence_2.fasta");
    f = fopen(pathSequence, "w");
    fprintf(f, "%s", sequence);
    fclose(f);
    free(sequence);

    nbLine = (int) ceil(MAX_SIZE_CONSTRAINT / 70.0);
    maxSizeData = 70 * nbLine;
    
    sequence = malloc(maxSizeData * sizeof(char));
    srand(time(NULL));

    for(d = 0; d < MAX_NUMBER_CONSTRAINT; d++){
        rewind(file);
        sprintf(sequence, "%s", "");
        
        offset = 15000 + (rand() % 40000);

        fseek(file, 72 * offset, SEEK_SET);

        for(i = 0; i < nbLine; i++){
            fscanf(file, "%70s\n", tmp);
            sequence = concat(3, sequence, tmp, (i == nbLine - 1 ? "\n\n" : "\n"));
        }

        pathSequence = concat(1, "../resources/datasets/dna/NC_002516.2/NC_002516.2_3/constraint.fasta");
        f = fopen(pathSequence, (d == 0? "w" : "a"));
        fprintf(f, "%s", sequence);
        fclose(f);
        free(sequence);
    }

    fclose(file);

}


ListSequenceData* readFileData(int numSequence){
    int i, d, nbLine;
    char* tmp, *sequence, *pathSequence;
    
    ListSequenceData* listSequenceData = createListSequenceData();

    tmp = malloc(70 * sizeof(char));
    pathSequence = malloc(1000 * sizeof(char));
    
    if(numSequence == 1)
        pathSequence = concat(1, getFilenameFirstSequenceDatasets());
    else if (numSequence == 2)
        pathSequence = concat(1, getFilenameSecondSequenceDatasets());
    else
        pathSequence = concat(1, getFilenameConstraintDatasets());
        
    file = fopen(pathSequence, "r");
    if (file == NULL) {
			logE("Error: missing file(s) %s", pathSequence);
			exit (EXIT_FAILURE);
	}    

    if(numSequence == 1 || numSequence == 2){
        nbLine = (int) ceil(MAX_SIZE_SEQUENCE / 70.0);
        sequence = malloc(MAX_SIZE_SEQUENCE * sizeof(char));
        int j = 0;
        for(i = 0; i < nbLine; i++){
            fscanf(file, "%70s\n", tmp);
            for(j = 0; j < 70; j++)
                sequence[j+i*70] = tmp[j];
        }

        free(tmp);
        tmp = malloc(MAX_SIZE_SEQUENCE * sizeof(char));
        sprintf(tmp, "%s", sequence + offsetSequence);
        free(sequence);
        sequence = malloc((maxSizeData+1) * sizeof(char));
        for(i = 0; i < maxSizeData; i++)
                sequence[i] = tmp[i];
        sequence[maxSizeData] = '\0';
        addSequenceToListSequenceData(listSequenceData, sequence);
    }
    else{
        nbLine = (int) ceil(MAX_SIZE_CONSTRAINT / 70.0);
        sequence = malloc(MAX_SIZE_CONSTRAINT * sizeof(char));
        long sizeConstraint;
        for(d = 0; d < maxNumberConstraint; d++){
            sprintf(sequence, "%s", "");
            for(i = 0; i < nbLine; i++){
                fscanf(file, "%70s\n", tmp);
                sequence = concat(2, sequence, tmp);
            }
            free(tmp);
            tmp = malloc(MAX_SIZE_CONSTRAINT * sizeof(char));
            sprintf(tmp, "%s", sequence + offsetConstraint);
            free(sequence);
            sizeConstraint = (long) floor(maxSizeConstraint / maxNumberConstraint);
            if((sizeConstraint * maxNumberConstraint + d + 1) <= maxSizeConstraint)
                sizeConstraint++;
            sequence = malloc((sizeConstraint+1) * sizeof(char));
            for(i = 0; i < sizeConstraint; i++)
                sequence[i] = tmp[i];
            sequence[sizeConstraint] = '\0'; 
            addSequenceToListSequenceData(listSequenceData, sequence);
            free(tmp);
            free(sequence);
            tmp = malloc(70 * sizeof(char));
            sequence = malloc(MAX_SIZE_CONSTRAINT * sizeof(char));
        }
    }

    fclose(file);

    return listSequenceData;
}

char* getAlphabet(){
    char* pathSequence, *alpha, *symbol;
    int i = 0;
    alpha = malloc(24 * sizeof(char));
    symbol = malloc(sizeof(char));
    pathSequence = malloc(1000 * sizeof(char));
    pathSequence = concat(1, getFilenameSymbolsDatasets());
    
    file = fopen(pathSequence, "r");
    if (file == NULL) {
		logE("error: missing file(s) %s\n", pathSequence);
		exit (EXIT_FAILURE);
	}

    while(!feof(file)){
        fscanf(file, "%s\n", symbol);
        alpha[i] = symbol[0];
        i++;
    }
    alpha[i] = '\0';
    free(symbol);
    fclose(file);
    return alpha;
}

long getNewNumberConstraint(){
    return getNumberSequenceOfListSequenceData(listConstraintSequenceData);
}

long getNewSizeConstraint(){
    return getSizeAllSequenceOfListSequenceData(listConstraintSequenceData);
}

ListSequenceData* getListFistSequenceData(){
    return listFistSequenceData;
}

ListSequenceData* getListSecondSequenceData(){
    return listSecondSequenceData;
}

ListSequenceData* getListConstraintSequenceData(){
    return listConstraintSequenceData;
}

ListSequenceData* createListSequenceData(){
    ListSequenceData* list = (ListSequenceData*) malloc(sizeof(ListSequenceData));
    if(list == NULL){
        exit(EXIT_FAILURE);
    }
    list->sequence = (char*) malloc(MAX_SIZE_SEQUENCE * sizeof(char));
    sprintf(list->sequence, "%s", "");
    list->next = NULL;
    return list;
}

int isEmptyListSequenceData(ListSequenceData* listSequenceData){
    return strcmp(listSequenceData->sequence, "") == 0 || listSequenceData == NULL;
}

int isSubStringOfAnySequence(ListSequenceData* listSequenceData, char* sequence){
    if(isEmptyListSequenceData(listSequenceData))
        return 0;
    else{
        ListSequenceData* path = listSequenceData;
        int isPresent = 0;
        while(path != NULL && isPresent == 0){
            if(strstr(path->sequence, sequence) != NULL)
                isPresent = 1;
            path = path->next;
        }
        return isPresent;
    }
}

void addSequenceToListSequenceData(ListSequenceData* listSequenceData, char* sequence){
    if(isEmptyListSequenceData(listSequenceData)){
        if(listSequenceData == NULL)
            listSequenceData = createListSequenceData();
        strcpy(listSequenceData->sequence, sequence);
    }
    else if(!isSubStringOfAnySequence(listSequenceData, sequence)){
        ListSequenceData* path = listSequenceData;
        while(path->next != NULL)
            path = path->next;
        ListSequenceData* new = createListSequenceData();
        strcpy(new->sequence, sequence);
        path->next = new;
    }
}

long getNumberSequenceOfListSequenceData(ListSequenceData* listSequenceData){
    int count = 0;
    if(isEmptyListSequenceData(listSequenceData))
        return 0;
    else{
        ListSequenceData* path = listSequenceData;
        while(path != NULL){
            ++count;
            path = path->next;
        }
        return count;
    }
}

long getSizeAllSequenceOfListSequenceData(ListSequenceData* listSequenceData){
    int count = 0;
    if(isEmptyListSequenceData(listSequenceData))
        return 0;
    else{
        ListSequenceData* path = listSequenceData;
        while(path != NULL){
            count += (long) strlen(path->sequence);
            path = path->next;
        }
        return count;
    }
}

void freeListSequenceData(ListSequenceData* listSequenceData){
    free(listSequenceData);
}
