/**
 * @file output.c
 * @author Jerry Lacmou (jerrylacmou@gmail.com)
 * @brief 
 * @version 4.1
 * @date 2019-02-11
 * 
 * @copyright Compiii Thesis (c) 2019 
 * 
 */

#include "output.h"
#include "clogger.h"

void saveStatistic(char* filename, GeneralInfo genInfo, InputInfo inInfo, OutputInfo outInfo){
    FILE* file = fopen(filename, "a");
    if (file == NULL){
		logW("Error: no such file %s from process %d on %s\r\n", filename, genInfo.rank, genInfo.processorName);
        return;
    }
    fseek(file, 0L, SEEK_END);
	if (ftell(file) == 0) {
		//fprintf(file, "%s;", "soft-name");
		//fprintf(file, "%s;", "version");
		//fprintf(file, "%s;", "algorithm");
		//fprintf(file, "%s;", "processor-name");
		fprintf(file, "%s;", "max-process");
		fprintf(file, "%s;", "rank");
		//fprintf(file, "%s;", "dataset");
		//fprintf(file, "%s;", "offset");

		fprintf(file, "%s;", "max-size-sequence");
		fprintf(file, "%s;", "max-size-constraint");
		fprintf(file, "%s;", "max-number-constraint");
		fprintf(file, "%s;", "max-fragmentation");

		fprintf(file, "%s;", "total-time");
		fprintf(file, "%s;", "calcul-time");
		fprintf(file, "%s;", "com-time");
        fprintf(file, "%s;", "idleness-time");
		fprintf(file, "%s;", "calcul-rate");
		fprintf(file, "%s;", "com-rate");
        fprintf(file, "%s;", "idleness-rate");
        fprintf(file, "%s;", "nb-message-exchange");
		fprintf(file, "%s;", "size-data-com");
        fprintf(file, "%s;", "size-data-com-in-mb");
        fprintf(file, "%s;", "ssclcs");
		fprintf(file, "\n");
	}

    //fprintf(file, "%s;", genInfo.softName);
	//fprintf(file, "%s;", genInfo.softVersion);
	//fprintf(file, "%s;", genInfo.algorithm);
	//fprintf(file, "%s;", genInfo.processorName);
	fprintf(file, "%d;", genInfo.maxProcess);
	fprintf(file, "%d;", genInfo.rank);
	//fprintf(file, "%s;", genInfo.dataset);
	//fprintf(file, "%s;", genInfo.offset);

	fprintf(file, "%ld;", inInfo.maxSizeSequence);
	fprintf(file, "%ld;", inInfo.maxSizeConstraint);
	fprintf(file, "%ld;", inInfo.maxNumberConstraint);
	fprintf(file, "%ld;", inInfo.maxFragmentation);
	
    fprintf(file, "%f;", outInfo.totalTime);
	fprintf(file, "%f;", outInfo.calculTime);
	fprintf(file, "%f;", outInfo.comTime);
    fprintf(file, "%f;", outInfo.idlenessTime);
	fprintf(file, "%f;", outInfo.calculTime / outInfo.totalTime * 100);
	fprintf(file, "%f;", outInfo.comTime / outInfo.totalTime * 100);
    fprintf(file, "%f;", outInfo.idlenessTime / outInfo.totalTime * 100);
    fprintf(file, "%ld;", outInfo.nbMessageExchange);
	fprintf(file, "%ld;", outInfo.sizeDataCom);
    fprintf(file, "%f;", outInfo.sizeDataCom / CONVERT_TO_MEGABYTE);
    fprintf(file, "%ld;", (long) outInfo.ssclcs);
	fprintf(file, "\n");
    
    fclose(file);
}