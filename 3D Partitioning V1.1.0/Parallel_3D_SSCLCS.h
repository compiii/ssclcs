#ifndef THREEDIMENSIONALPARTITIONING_H_ 
#define THREEDIMENSIONALPARTITIONING_H_
#include <stdlib.h>
#include <stdio.h>
#include <string.h> 
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <mpi.h>

#define controlFile 0
#define LIMIT 4

/**
 * @brief Definition of the Coord_T data structure that refers to the coordinates of a block
 * 
 */
typedef struct Coord_T Coord;
struct Coord_T{
    int i;
    int j;
    int k;
};

typedef struct Boundary_T Boundary;
struct Boundary_T{
    int lig;
    int col;
};

/**
 * @brief Definition of the Dimension_T data structure that refers to the dimensions of a block
 * 
 */
typedef struct Dimension_T Dim;
struct Dimension_T{
    int nbrLine;
    int nbrColumn;
    int plan;
};

/**
 * @brief Definition of the BlockCoreData_T data structure that refers to the core data of 
 * a block
 * 
 */
typedef struct BlockData_T BlockData;
struct BlockData_T{
    int id;
    Coord coord;
    Dim dim;
    int diag;
    int rank;
};

typedef struct ListPartitions_L ListPartitions;
struct ListPartitions_L{
    char* partition;
};

typedef struct ListstartPostTab_LX ListstartPostTab_X;
struct ListstartPostTab_LX{
    long* startPost;
};

typedef struct ListstartPostTab_LY ListstartPostTab_Y;
struct ListstartPostTab_LY{
    long* startPost;
};

int maxNumberConstraint, offsetSequence, offsetConstraint;
int maxDiag, maxProcess, ch1, rank, maxBlock, maxBlockDiag, processorNameSize,
    sizePartA, sizePartB, depth, maxBlockPerProcess, maxk, maxkrecv;
double startTime, tmpDouble, endTime, calculTime, comTime, comRate, calculRate, totalTime, idlenessTime, idlenessRate;
long sizeDataCom, nbMessageExchange, ssclcs;
BlockData* arrayData = NULL;        //Contains the address where the BlockData values are store
char *X = NULL;
char *Y = NULL;
char *Sigma = NULL;
int t,nb_ch, TAILLE_MAX, sizeMaxC, contraintLength;
long r, m, s, l, n;
char* partition = NULL;
char* partition2 = NULL;
char* partitionTmp = NULL;
long* prevMatchTab_X = NULL;
long* prevMatchTab_Y = NULL;
long* f1 = NULL;
long* f2 = NULL;
long* f1Tmp = NULL;
long* f2Tmp = NULL;
ListPartitions* listparts = NULL;
ListPartitions* listpartsTmp = NULL;
ListstartPostTab_X* listStartPos_X = NULL;
ListstartPostTab_Y* listStartPos_Y = NULL;

MPI_Request request;
MPI_Status  statut;
MPI_Request request1;
MPI_Status  status1;
MPI_Status status;

char filename[14];                    /* Name of file */
char convertstr[5];
FILE *fp = NULL;

char processorName[MPI_MAX_PROCESSOR_NAME];
MPI_Datatype vectorLin, vectorCol, line, column, tmpcolumn, linePart;

/**
 * Functions prototypes
 * 
 */
void init();
void constructBlockData();
Dim computeDimBlock(int plan, int n1, int m1);
int sizeBlock(int sizeString, int sizePartString, int mn);
int profondeur(int numproc);

void printEnterInfo();
void printBlockData(BlockData blockData);
void printTabBlockData();

FILE* createFile();
void fileBlockData(FILE *fp, BlockData bd);

//******* Fonction pour le calcul des blocs ***************//
void initArrayDP();
void preliminaryWork();
void blockCompute(BlockData blockData);
Boundary computeLimitBlock(BlockData blockData);
void allBlockCompute();
void computeBlock(int lig, int col, BlockData blockData);
int sigma(int k, char ch);
void findAlphabet();
void findAlphabet_standard();
void allocations();
void getData();
void loadSequence_X(char* tmp);
void loadSequence_Y(char* tmp);

void freeMemory();

long* prevMatch(char* stringChar);
void printPrevMatch(char* stringChar, long* prevM);
long* startPosCom(char* stringChar, char* stringPart, long* prevM);
void printStartPos(char* stringChar, char* stringPart, long* startP);
long searchPrevMatchTab(long indiceSigma, long indiceX, long* prevM);
long searchIndiceInSigma(char ch);
long computeWeightPartition(char* c1, char* c2);
long i_j_k2(long i, long j, long k);
long i_j_k3(long i, long j, long k);
long i_j(long i, long j, long stringLength);
char searchPartitionInList(long k);
long maxNumber(long nb1, long nb2, long nb3);
int isSamePlanScale(int k);
void checkAllocationMemory1D(long* tab, char* text);
void checkAllocationMemory2D(long** tab, char* text);
void checkAllocationMemory3D(long*** tab, char* text);

#endif
