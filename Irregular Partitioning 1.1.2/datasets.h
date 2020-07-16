/**
 * @file datasets.h
 * @author Jerry Lacmou (jerrylacmou@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2019-01-11
 * 
 * @copyright Compiii Thesis (c) 2019 
 * 
 */
#ifndef DATASETS_H_
#define DATASETS_H_

//#include "threeDimensionalPartitioning.h"
#include "utils.h"

#define MAX_SIZE_SEQUENCE getMaxSizeSequenceDatasets()
#define MAX_SIZE_CONSTRAINT getMaxSizeConstraintDatasets()
#define MAX_NUMBER_CONSTRAINT getMaxNumberConstraintDatasets()

typedef struct ListSequenceData_T ListSequenceData;
struct ListSequenceData_T{
    char* sequence;
    ListSequenceData *next;
};

extern int offsetSequence, offsetConstraint, maxNumberConstraint;
int maxSizeData, maxSizeConstraint;//, maxNumberConstraint, offsetSequence, offsetConstraint;
FILE* file;

ListSequenceData* listFistSequenceData, *listSecondSequenceData, *listConstraintSequenceData;

void initDataSetsMod();
void createDatasets(char* datasets);
ListSequenceData* readFileData(int numSequence);
char* getAlphabet();
long getNewNumberConstraint();
long getNewSizeConstraint();

ListSequenceData* getListFistSequenceData();
ListSequenceData* getListSecondSequenceData();
ListSequenceData* getListConstraintSequenceData();

ListSequenceData* createListSequenceData();
int isEmptyListSequenceData(ListSequenceData* listSequenceData);
void addSequenceToListSequenceData(ListSequenceData* listSequenceData, char* sequence);
void freeListSequenceData(ListSequenceData* listSequenceData);
int isSubStringOfAnySequence(ListSequenceData* listSequenceData, char* sequence);
long getNumberSequenceOfListSequenceData(ListSequenceData* listSequenceData);
long getSizeAllSequenceOfListSequenceData(ListSequenceData* listSequenceData);

#endif
