#include <stdio.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include "../../src/AutoSync.h"

#define ARRAY_SIZE 600
#define NO_OF_THREADS 4
#define PATTERN 3

/**************************** Thread's prototypes *****************************/
void* SearchThread(void* args);

/****************************** Shared Variables ******************************/
uint32_t uiCountOccurrences = 0;
uint32_t uiSearchArray[ARRAY_SIZE];

void vFillArray(uint32_t* puiArray, size_t xSize);

/************************************* MAIN ***********************************/
int main(int argc, char const *argv[])
{  
  int16_t iRetVal;
  double ElapsedTime;
  pthread_t xThreadHandle[NO_OF_THREADS]; 

  /* Initialize the random number generator */
  srand(time(0));

  vFillArray(uiSearchArray, sizeof(uiSearchArray)/sizeof(uint32_t));
  
  /* To calculate the execution time */
  clock_t t;
  t = clock();

  iAutoSyncCreate();
  
  printf("We must have exactly %d occurrences!\n", ARRAY_SIZE/10);
  
  printf("Starting search threads...\n\n");
  /* Create POSIX threads */
  for (uint8_t i = 0; i < NO_OF_THREADS; i++)
  {
    uint8_t* a = malloc(sizeof(uint8_t));
    *a = i;

    iRetVal = pthread_create(&xThreadHandle[i],
                             NULL,
                             &SearchThread,
                             a);
  }
    
  /* Wait Search threads finished */
  for (uint8_t i = 0; i < NO_OF_THREADS; i++)
  {
    iRetVal = pthread_join(xThreadHandle[i], NULL);   
  }

  iAutoSyncDestroy();

  printf("Total No Of Occurrences: %d\n", uiCountOccurrences);

  t = clock() - t;
  ElapsedTime = ((double)t)/CLOCKS_PER_SEC; // in seconds
  printf("\n >>>>> Execution Time: %f seconds \n", ElapsedTime);

  return 0;
}

void vFillArray(uint32_t* puiArray, size_t xSize)
{    
    uint8_t num;
    
    /* Generate random number between 5 and 10                */
    /* For the indexes multiple of 10, insert pattern         */
    /* That way we have a deterministic number of occurrences */
    for (uint32_t i = 0; i < xSize; i++)
    {
        puiArray[i] = (rand() % (6)) + 5;

        if (i % 10 == 0)puiArray[i] = PATTERN;
    
    }      
}

/****************************** Thread's Bodies *******************************/
void* SearchThread(void* args)
{
    uint32_t uiIndex = *(uint8_t*)args;
    uint32_t uiStart = uiIndex * (ARRAY_SIZE/NO_OF_THREADS);
    uint32_t uiEnd = uiStart + ARRAY_SIZE/NO_OF_THREADS -1;
    uint32_t uiCountLocal;

    for (uint32_t i = uiStart; i < uiEnd; i++)
    {
        if(uiSearchArray[i] == PATTERN)
        {            
            iAutoSyncReadToUpdate(&uiCountLocal, &uiCountOccurrences, sizeof(uiCountLocal));
            uiCountLocal++;
            iAutoSyncUpdate(&uiCountOccurrences, &uiCountLocal, sizeof(uiCountLocal));
        }        
    } 
}





