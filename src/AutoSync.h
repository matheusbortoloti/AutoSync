#ifndef __AUTO_SYNC_H__
#define __AUTO_SYNC_H__

#include <pthread.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>

/* Comment to suppress AutoSync messages */
#define AUTO_SYNC_VERBOSE

/* Definiton of custom errors */
#define AUTO_SYNC_OK 0

/* Definition of constants */
#define MAX_DEPENDENCIES 10 /* Increase it if more is needed */

/* EXTERNAL VARIABLES */

/* DATA STRUCTURES */
typedef struct xAutoSyncIntentionsStruct 
{  
  void* pvDependsOn[MAX_DEPENDENCIES];  
  bool bConstantInitByMain;  
  bool bSlicedArray;
  uint64_t uiFirstAccess;
  uint64_t uiLastAccess;  
} const xAutoSyncIntentions;  

typedef int8_t xAutoSyncEvent;


/*******************************************************************************
*                              INTERFACE DEFINTION
*******************************************************************************/
int8_t iAutoSyncCreate(void);
int8_t iAutoSyncDestroy(void);
int8_t iAutoSyncRead(void* pvValue, void* pvSharedVar, size_t xSizeData, xAutoSyncIntentions xIntention);
int8_t iAutoSyncWrite(void* pvSharedVar, void* pvValue, size_t xSizeData, xAutoSyncIntentions xIntention);

int8_t iAutoSyncReadToUpdate(void* pvValue, void* pvSharedVar, size_t xSizeData, xAutoSyncIntentions xIntention);
int8_t iAutoSyncUpdate(void* pvSharedVar, void* pvValue, size_t xSizeData, xAutoSyncIntentions xIntention);

int8_t iAutoSyncSharedVarAsArg(void* pvSharedVar);

int8_t iAutoSyncProceedOnEvent(xAutoSyncEvent xEvent, uint8_t uiNoOfThreads); 
#endif
