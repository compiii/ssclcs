/**
 * @file config.h
 * @author Jerry Lacmou (jerrylacmou@gmail.com)
 * @brief 
 * @version 4.1
 * @date 2019-02-11
 * 
 * @copyright Compiii Thesis (c) 2019 
 * 
 */
#ifndef CONFIG_H_
#define CONFIG_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "jsmn.h"

#define JSON_FILE_PATH "Resources/config/config.json"
#define BUFFER_SIZE 5000
#define MAX_TOKEN_COUNT 128

void readfile(char* fileContent);
char* getValueByKey(char* key);
int jsoneq(const char *json, jsmntok_t *tok, const char *s);

char* getSoftName();
char* getSoftVersion();
char* getInputDataSetId();
char* getInputDataSetOffset();

char* getFileNamePrimarySurvey(int number);
char* getFileNameSecondarySurvey(int number);

char* getFilenameDatasets();
char* getFilenameSymbolsDatasets();
char* getFilenameFirstSequenceDatasets();
char* getFilenameSecondSequenceDatasets();
char* getFilenameConstraintDatasets();
char* getDirnameCreateDatasets();

long getMaxSizeSequenceDatasets();
long getMaxSizeConstraintDatasets();
long getMaxNumberConstraintDatasets();
long getMaxNumberFragmentation();
long getMaxSizeSequenceInputDatasets();
long getOffsetSequenceInputDatasets();
long getMaxSizeConstraintInputDatasets();
long getMaxNumberConstraintInputDatasets();
long getOffsetConstraintInputDatasets();

int getAlgorithmInputDatasets();
#endif