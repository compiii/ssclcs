/**
 * @file output.h
 * @author Jerry Lacmou (jerrylacmou@gmail.com)
 * @brief 
 * @version 4.1
 * @date 2019-02-11
 * 
 * @copyright Compiii Thesis (c) 2019 
 * 
 */
#ifndef OUTPUT_H_
#define OUTPUT_H_

#define MAX_SIZE_STRING 1024
#define CONVERT_TO_MEGABYTE 1048576.0

typedef struct GeneralInfo_T GeneralInfo;
struct GeneralInfo_T{
    char softName[MAX_SIZE_STRING]; 
    char softVersion[MAX_SIZE_STRING];
    char dataset[MAX_SIZE_STRING];
    char algorithm[MAX_SIZE_STRING];
    char offset[MAX_SIZE_STRING];
    char processorName[MAX_SIZE_STRING];
    int maxProcess;
    int rank;  
};

typedef struct InputInfo_T InputInfo;
struct InputInfo_T{
    long maxSizeSequence;
    long maxSizeConstraint;
    long maxNumberConstraint;
    long maxFragmentation;
};

typedef struct OutputInfo_T OutputInfo;
struct OutputInfo_T{
    unsigned short ssclcs;
    double totalTime;
    double calculTime;
    double comTime;
    double idlenessTime;
    long nbMessageExchange;
    long sizeDataCom;
};

void saveStatistic(char* filename, GeneralInfo genInfo, InputInfo inInfo, OutputInfo outInfo);

#endif
