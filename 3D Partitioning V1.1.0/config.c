/**
 * @file config.c
 * @author Jerry Lacmou (jerrylacmou@gmail.com)
 * @brief 
 * @version 4.1
 * @date 2019-02-11
 * 
 * @copyright Compiii Thesis (c) 2019 
 * 
 */
#include "config.h"
#include "utils.h"
#include "clogger.h"

void readfile(char* fileContent){
    FILE *f;
    char c;
    int index = 0;
    f = fopen(JSON_FILE_PATH, "rt");
    if(f == NULL){
        printf("Impossible");
        logE("Failed to open %s", JSON_FILE_PATH);
        exit(EXIT_FAILURE);
    }
    while((c = fgetc(f)) != EOF){
        fileContent[index] = c;
        index++;
    }
    fileContent[index] = '\0';
    fclose(f);
}

char* getValueByKey(char* key){
    char JSON_STRING[BUFFER_SIZE];

    char* value = malloc(1024 * sizeof(char));
    readfile(JSON_STRING);

    int i;
    int r;

    jsmn_parser p;
    jsmntok_t t[MAX_TOKEN_COUNT];

    jsmn_init(&p);

    r = jsmn_parse(&p, JSON_STRING, strlen(JSON_STRING), t, sizeof(t)/(sizeof(t[0])));

    if (r < 0) {
        logE("Failed to parse JSON: %d", r);
        exit(EXIT_FAILURE);
    }

    /* Assume the top-level element is an object */
    if (r < 1 || t[0].type != JSMN_OBJECT) {
        logE("Object expected");
        exit(EXIT_FAILURE);
    }

    for (i = 1; i < r; i++){
        if (jsoneq(JSON_STRING, &t[i], key) == 0){
            jsmntok_t json_value = t[i+1];
            int string_length = json_value.end - json_value.start;
            int idx;
            for (idx = 0; idx < string_length; idx++){
                value[idx] = JSON_STRING[json_value.start + idx ];
            }
            value[string_length] = '\0';
            return value;
        }
        i++;
    }
    logW("The key %s is not found in config.json", key);
    return NULL;
}

int jsoneq(const char *json, jsmntok_t *tok, const char *s){
    if (tok->type == JSMN_STRING && (int) strlen(s) == tok->end - tok->start &&
			strncmp(json + tok->start, s, tok->end - tok->start) == 0) {
		return 0;
	}
	return -1;
}

char* getFilenameDatasets(){
    return concat(5, getValueByKey("datasets-dir"), getValueByKey("datasets-type"), getValueByKey("datasets-dirname"), getValueByKey("datasets-filename"), getValueByKey("datasets-extension"));
}

char* getFilenameSymbolsDatasets(){
    return concat(4, getValueByKey("datasets-dir"), getValueByKey("datasets-type"), getValueByKey("datasets-symbols-filename"), getValueByKey("datasets-symbols-extension"));
}

char* getFilenameFirstSequenceDatasets(){
    return concat(7, getValueByKey("datasets-dir"), getValueByKey("datasets-type"), getValueByKey("datasets-dirname"), getValueByKey("datasets-filename"), getValueByKey("datasets-input-id"), "sequence_1", getValueByKey("datasets-extension"));
}

char* getFilenameSecondSequenceDatasets(){
    return concat(7, getValueByKey("datasets-dir"), getValueByKey("datasets-type"), getValueByKey("datasets-dirname"), getValueByKey("datasets-filename"), getValueByKey("datasets-input-id"), "sequence_2", getValueByKey("datasets-extension"));
}

char* getFilenameConstraintDatasets(){
    return concat(7, getValueByKey("datasets-dir"), getValueByKey("datasets-type"), getValueByKey("datasets-dirname"), getValueByKey("datasets-filename"), getValueByKey("datasets-input-id"), "constraint", getValueByKey("datasets-extension"));
}

char* getDirnameCreateDatasets(){
    return concat(5, getValueByKey("datasets-dir"), getValueByKey("datasets-type"), getValueByKey("datasets-dirname"), getValueByKey("datasets-filename"), getValueByKey("datasets-create-id"));
}

long getMaxSizeSequenceDatasets(){
    return atol(getValueByKey("datasets-max-size-sequence"));
}

long getMaxSizeConstraintDatasets(){
    return atol(getValueByKey("datasets-max-size-constraint"));
}

long getMaxNumberConstraintDatasets(){
    return atol(getValueByKey("datasets-max-number-constraint"));
}

long getMaxNumberFragmentation(){
    return atol(getValueByKey("datasets-input-max-fragmentation"));
}

long getMaxSizeSequenceInputDatasets(){
    return atol(getValueByKey("datasets-input-max-size-sequence"));
}

long getOffsetSequenceInputDatasets(){
    return atol(getValueByKey("datasets-input-offset-sequence"));
}

long getMaxSizeConstraintInputDatasets(){
    return atol(getValueByKey("datasets-input-max-size-constraint"));
}

long getMaxNumberConstraintInputDatasets(){
    return atol(getValueByKey("datasets-input-max-number-constraint"));
}

long getOffsetConstraintInputDatasets(){
    return atol(getValueByKey("datasets-input-offset-constraint"));
}

int getAlgorithmInputDatasets(){
    return atoi(getValueByKey("datasets-input-algorithm"));
}

char* getSoftName(){
    return getValueByKey("soft-name");
}

char* getSoftVersion(){
    return getValueByKey("soft-version");
}

char* getInputDataSetId(){
    return concat(2,getValueByKey("datasets-filename"),getValueByKey("datasets-input-id"));
}

char* getInputDataSetOffset(){
    return concat(3,getValueByKey("datasets-input-offset-sequence"),"-",getValueByKey("datasets-input-offset-constraint"));
}

char* getFileNamePrimarySurvey(int number){
    char* num;
    num = (char *) malloc(3 * sizeof(char));
    sprintf(num, "%d", number);
    return concat(4, getValueByKey("output-dirname"), getValueByKey("output-filename"), num, getValueByKey("output-extension"));
}

char* getFileNameSecondarySurvey(int number){
    char* num;
    num = (char *) malloc(3 * sizeof(char));
    sprintf(num, "%d", number);
    return concat(5, getValueByKey("output-dirname"), getValueByKey("output-dirname-secondary"), getValueByKey("output-filename"), num, getValueByKey("output-extension"));
}