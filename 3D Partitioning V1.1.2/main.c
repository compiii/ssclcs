#include "Parallel_3D_SSCLCS.h"
#include "Parallel_3D_Partition_SSCLCS.c"
#include "Parallel_3D_Computation_SSCLCS.c"
#include "clogger.h"
#include "clogger.c"
#include "config.c"
#include "config.h"
#include "datasets.c"
#include "datasets.h"
#include "output.c"
#include "output.h"
#include "utils.c"
#include "utils.h"
#include "jsmn.c"
#include "jsmn.h"

int main(int argc, char **argv)
{   

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &maxProcess);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(processorName, &processorNameSize);
    MPI_Status status;
    
    if (argc == 4){
	    TAILLE_MAX = atol(argv[1]);
	    contraintLength = atol(argv[2]);
        maxNumberConstraint = atol(argv[3]);
        // offsetSequence = atol(argv[3]);
        // offsetConstraint = atol(argv[4]);
    }else{
        printf("Arguments incorrect !\n");
        exit(1);
    }
    GeneralInfo generalInfo;
    InputInfo inputInfo;
    OutputInfo outputInfo;
    char fileName[1000];
    char* print = malloc(sizeof(char));
    char tmp[1024];
    generalInfo.rank = rank;
    generalInfo.maxProcess = maxProcess;

    inputInfo.maxSizeSequence = TAILLE_MAX;
    inputInfo.maxSizeConstraint = contraintLength;
    inputInfo.maxNumberConstraint = maxNumberConstraint;
    inputInfo.maxFragmentation = 0;

    m=TAILLE_MAX, s=4, l=maxNumberConstraint, n=TAILLE_MAX;
    // contraintLength = (m/4)+(m%4);
    startTime = tmpDouble = endTime = calculTime = comTime = comRate = calculRate = 
    totalTime = idlenessTime = idlenessRate = sizeDataCom = nbMessageExchange = 0;
    init();
    MPI_Barrier(MPI_COMM_WORLD);
    constructBlockData();
    MPI_Barrier(MPI_COMM_WORLD);
    allocations();
    MPI_Barrier(MPI_COMM_WORLD);
    getData();
    MPI_Barrier(MPI_COMM_WORLD);
    preliminaryWork();
    MPI_Barrier(MPI_COMM_WORLD);
    startTime = MPI_Wtime();
    allBlockCompute();
    MPI_Barrier(MPI_COMM_WORLD);
    endTime = MPI_Wtime(); 
    idlenessTime = endTime - startTime - (calculTime + comTime);
        totalTime = calculTime + comTime + idlenessTime;

        calculRate = calculTime / totalTime * 100;
        comRate = comTime / totalTime * 100;
        idlenessRate = idlenessTime / totalTime * 100;

        
        outputInfo.totalTime = totalTime;
        outputInfo.calculTime = calculTime;
        outputInfo.comTime = comTime;
        outputInfo.idlenessTime = idlenessTime;
        outputInfo.nbMessageExchange = nbMessageExchange;
        outputInfo.sizeDataCom = sizeDataCom * sizeof(unsigned short);
	    outputInfo.ssclcs = 0;


    if(rank == maxProcess - 1){
	outputInfo.ssclcs = f1[i_j_k2(m-1,n-1,depth-1)];
            strcpy(fileName, getFileNamePrimarySurvey(rank));
            saveStatistic(fileName, generalInfo, inputInfo, outputInfo);
               print = concat(1, "");

                sprintf(tmp, "Output datas of %s (process %d) for Hermann", processorName, rank);
                print = concat(2, tmp, "\n------------------------------\n");
                sprintf(tmp, "ssclcs: %d\n", outputInfo.ssclcs);
                print = concat(2, print, tmp);
                sprintf(tmp, "total time: %f\n", totalTime);
                print = concat(2, print, tmp);
                sprintf(tmp, "calcul time: %f\n", calculTime);
                print = concat(2, print, tmp);
                sprintf(tmp, "comunication time: %f\n", comTime);
                print = concat(2, print, tmp);
                sprintf(tmp, "idleness time: %f\n", idlenessTime);
                print = concat(2, print, tmp);
                sprintf(tmp, "calcul rate: %f\n", calculRate);
                print = concat(2, print, tmp);
                sprintf(tmp, "comunication rate: %f\n", comRate);
                print = concat(2, print, tmp);
                sprintf(tmp, "idleness rate: %f\n", idlenessRate);
                print = concat(2, print, tmp);
                sprintf(tmp, "number of message exchange: %ld\n", nbMessageExchange);
                print = concat(2, print, tmp);
                sprintf(tmp, "size of data communicate: %f MB\n", sizeDataCom * sizeof(unsigned short) / 1048576.0);
                print = concat(2, print, tmp);
                logI(print); 
       
    }
    else{
        strcpy(fileName, getFileNameSecondarySurvey(rank));
        saveStatistic(fileName, generalInfo, inputInfo, outputInfo);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    freeMemory();
    MPI_Finalize();
    return 0;
}
