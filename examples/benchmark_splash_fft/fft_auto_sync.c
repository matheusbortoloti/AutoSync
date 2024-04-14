/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

/*************************************************************************/
/*                                                                       */
/*  Perform 1D fast Fourier transform using six-step FFT method          */
/*                                                                       */
/*  1) Performs staggered, blocked transposes for cache-line reuse       */
/*  2) Roots of unity rearranged and distributed for only local          */
/*     accesses during application of roots of unity                     */
/*  3) Small set of roots of unity elements replicated locally for       */
/*     1D FFTs (less than root N elements replicated at each node)       */
/*  4) Matrix data structures are padded to reduce cache mapping         */
/*     conflicts                                                         */
/*                                                                       */
/*  Command line options:                                                */
/*                                                                       */
/*  -mM : M = even integer; 2**M total complex data points transformed.  */
/*  -pP : P = number of processors; Must be a power of 2.                */
/*  -nN : N = number of cache lines.                                     */
/*  -lL : L = Log base 2 of cache line length in bytes.                  */
/*  -s  : Print individual processor timing statistics.                  */
/*  -t  : Perform FFT and inverse FFT.  Test output by comparing the     */
/*        integral of the original data to the integral of the data      */
/*        that results from performing the FFT and inverse FFT.          */
/*  -o  : Print out complex data points.                                 */
/*  -h  : Print out command line options.                                */
/*                                                                       */
/*  Note: This version works under both the FORK and SPROC models        */
/*                                                                       */
/*************************************************************************/

#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "../00_AutoSync/AutoSync.h" 
#define PAGE_SIZE               4096
#define NUM_CACHE_LINES        65536
#define LOG2_LINE_SIZE             4
#define PI                         3.1416
#define DEFAULT_M                 10
#define DEFAULT_P                  1

MAIN_ENV

#define SWAP_VALS(a,b) {double tmp; tmp=a; a=b; b=tmp;}

struct GlobalMemory {
  long id;
  long *transtimes;
  long *totaltimes;
  unsigned long starttime;
  unsigned long finishtime;
  unsigned long initdonetime;
} *Global;

int  is_output = 1;
long P = DEFAULT_P;
long M = DEFAULT_M;
long N;                  /* N = 2^M                                */
long rootN;              /* rootN = N^1/2                          */
double *x;              /* x is the original time-domain data     */
double *trans;          /* trans is used as scratch space         */
double *umain;          /* umain is roots of unity for 1D FFTs    */
double *umain2;         /* umain2 is entire roots of unity matrix */
long test_result = 0;
long doprint = 0;
long dostats = 0;
/* long transtime = 0; AutoSync: there is no reason for this variable be global! */
/* long transtime2 = 0; AutoSync: there is no reason for this variable be global! */
/* long avgtranstime = 0; AutoSync: there is no reason for this variable be global! */
/* long avgcomptime = 0; AutoSync: there is no reason for this variable be global! */
/* unsigned long transstart = 0; AutoSync: variable not used! */
/* unsigned long transend = 0; AutoSync: variable not used! */
/* long maxtotal=0; AutoSync: there is no reason for this variable be global! */
/* long mintotal=0; AutoSync: there is no reason for this variable be global! */
/* double maxfrac=0; AutoSync: there is no reason for this variable be global! */
/* double minfrac=0; AutoSync: there is no reason for this variable be global! */
/* double avgfractime=0; AutoSync: there is no reason for this variable be global! */
/* long orig_num_lines = NUM_CACHE_LINES; AutoSync: there is no reason for this variable be global! */
long num_cache_lines = NUM_CACHE_LINES;    /* number of cache lines */
/* long log2_line_size = LOG2_LINE_SIZE; AutoSync: there is no reason for this variable be global! */
/* long line_size; AutoSync: there is no reason for this variable be global! */
/* long rowsperproc; AutoSync: there is no reason for this variable be global! */
/* double ck1; AutoSync: there is no reason for this variable be global! */
/* double ck3; AutoSync: there is no reason for this variable be global! */
long pad_length; 

xAutoSyncIntentions xNoSpecialIntention;  
xAutoSyncIntentions xConstantInitByMain = {.bConstantInitByMain = true};
xAutoSyncIntentions xIntentionTransTimes = {.pvDependsOn[0] = &P}; /* Sliced Array!*/
xAutoSyncIntentions xIntentionTotalTimes = {.pvDependsOn[0] = &P}; /* Sliced Array!*/
xAutoSyncIntentions xIntentionN = {.bConstantInitByMain = true, .pvDependsOn[0] = &M};
xAutoSyncIntentions xIntentionX = {.pvDependsOn[0] = &N, .pvDependsOn[1] = &rootN, .pvDependsOn[2] = &pad_length};
xAutoSyncIntentions xIntentionTrans = {.pvDependsOn[0] = &N, .pvDependsOn[1] = &rootN, .pvDependsOn[2] = &pad_length};
xAutoSyncIntentions xIntentionUmain2 = {.bConstantInitByMain = true, .pvDependsOn[0] = &N, .pvDependsOn[1] = &rootN, .pvDependsOn[2] = &pad_length};
xAutoSyncIntentions xIntentionUmain = {.bConstantInitByMain = true, .pvDependsOn[0] = &rootN};
xAutoSyncIntentions xIntentionRootN = {.bConstantInitByMain = true, .pvDependsOn[0] = &M };

xAutoSyncEvent xTwiddleDone;
xAutoSyncEvent xFFT1DDone;
xAutoSyncEvent xFFTDone;
xAutoSyncEvent xTransposeDone;

xAutoSyncIntentions xIntentionSlicedArray = {.bSlicedArray = true, 
                                            uiFirstAccess = 0,
                                            uiLastAccess  = 0 };

iAutoSyncRead(&uiLocalVar, &uiSlicedArray, xIntentionSlicedArray);

void SlaveStart(void);
double TouchArray(double *x, double *scratch, double *u, double *upriv, long MyFirst, long MyLast);
double CheckSum(double *x);
void InitX(double *x);
void InitU(long N, double *u);
void InitU2(long N, double *u, long n1);
long BitReverse(long M, long k);
void FFT1D(long direction, long M, long N, double *x, double *scratch, double *upriv, double *umain2,
	   long MyNum, long *l_transtime, long MyFirst, long MyLast, long pad_length, long test_result, long dostats);
void TwiddleOneCol(long direction, long n1, long j, double *u, double *x, long pad_length);
void Scale(long n1, long N, double *x);
void Transpose(long n1, double *src, double *dest, long MyNum, long MyFirst, long MyLast, long pad_length);
void CopyColumn(long n1, double *src, double *dest);
void Reverse(long N, long M, double *x);
void FFT1DOnce(long direction, long M, long N, double *u, double *x);
void PrintArray(long N, double *x);
void printerr(char *s);
long log_2(long number);

void srand48(long int seedval);
double drand48(void);

int main(int argc, char *argv[])
{
  long i;
  long c;
  extern char *optarg;
  long m1;
  long factor;
  long pages;
  unsigned long start;

  long localP = DEFAULT_P;
  long localM = DEFAULT_M;
  long localN;                 
  long localrootN;   
  long local_num_cache_lines = NUM_CACHE_LINES;    
  long local_test_result = 0;
  long local_doprint = 0;
  long local_dostats = 0;
  long local_pad_length;
  long orig_num_lines = NUM_CACHE_LINES;
  long log2_line_size = LOG2_LINE_SIZE;
  long rowsperproc;
  long line_size;
  long transtime = 0;
  long totaltime = 0;
  long transtime2 = 0;
  long avgtranstime = 0;
  long avgcomptime = 0;
  long maxtotal=0;
  long mintotal=0;
  double maxfrac=0;
  double minfrac=0;
  double avgfractime=0;
  double ck1;
  double ck3;                        /* checksums for testing answer */
  long local_global_id = 0;  

  CLOCK(start);

  iAutoSyncCreate();

  while ((c = getopt(argc, argv, "p:m:n:l:stoh")) != -1) {
    switch(c) {
      case 'p': localP = atoi(optarg);                
                if (localP < 1) {
                  printerr("P must be >= 1\n");
                  exit(-1);
                }
                if (log_2(localP) == -1) {
                  printerr("P must be a power of 2\n");
                  exit(-1);
                }
                /* AutoSync: Since P is plausible, write it to the shared-variable */
                iAutoSyncWrite(&P, &localP, sizeof(P), xConstantInitByMain);
	        break;
      case 'm': localM = atoi(optarg);
                m1 = localM/2;
                if (2*m1 != localM) {
                  printerr("M must be even\n");
                  exit(-1);
                }
                /* AutoSync: Since M is plausible, write it to the shared-variable */
                iAutoSyncWrite(&M, &localM, sizeof(M), xConstantInitByMain);
	        break;
      case 'n': local_num_cache_lines = atoi(optarg);
                orig_num_lines = local_num_cache_lines;
                if (local_num_cache_lines < 1) {
                  printerr("Number of cache lines must be >= 1\n");
                  exit(-1);
                }
                /* AutoSync: Since num_cache_lines is plausible, write it to the shared-variable */
                iAutoSyncWrite(&num_cache_lines, &local_num_cache_lines, sizeof(num_cache_lines), xConstantInitByMain);
	        break;
      case 'l': log2_line_size = atoi(optarg);
                if (log2_line_size < 0) {
                  printerr("Log base 2 of cache line length in bytes must be >= 0\n");
                  exit(-1);
                }                
	        break;
      case 's': iAutoSyncReadToUpdate(&local_dostats, &dostats, sizeof(dostats), xConstantInitByMain);
                local_dostats = !local_dostats;
                iAutoSyncUpdate(&dostats, &local_dostats, sizeof(dostats), xConstantInitByMain);
	        break;
      case 't': iAutoSyncReadToUpdate(&local_test_result, &test_result, sizeof(test_result), xConstantInitByMain);
                local_test_result = !local_test_result;
                iAutoSyncUpdate(&test_result, &local_test_result, sizeof(test_result), xConstantInitByMain);               
	        break;
      case 'o': iAutoSyncReadToUpdate(&local_doprint, &doprint, sizeof(doprint), xConstantInitByMain);
                local_doprint = !local_doprint;
                iAutoSyncUpdate(&doprint, &local_doprint, sizeof(doprint), xConstantInitByMain);
	        break;
      case 'h': printf("Usage: FFT <options>\n\n");
                printf("options:\n");
                printf("  -mM : M = even integer; 2**M total complex data points transformed.\n");
                printf("  -pP : P = number of processors; Must be a power of 2.\n");
                printf("  -nN : N = number of cache lines.\n");
                printf("  -lL : L = Log base 2 of cache line length in bytes.\n");
                printf("  -s  : Print individual processor timing statistics.\n");
                printf("  -t  : Perform FFT and inverse FFT.  Test output by comparing the\n");
                printf("        integral of the original data to the integral of the data that\n");
                printf("        results from performing the FFT and inverse FFT.\n");
                printf("  -o  : Print out complex data points.\n");
                printf("  -h  : Print out command line options.\n\n");
                printf("Default: FFT -m%1d -p%1d -n%1d -l%1d\n",
                       DEFAULT_M,DEFAULT_P,NUM_CACHE_LINES,LOG2_LINE_SIZE);
		exit(0);
	        break;
    }
  }

  MAIN_INITENV(,80000000);  

  iAutoSyncRead(&localN, &N, sizeof(N), xIntentionN);
  iAutoSyncRead(&localM, &M, sizeof(M), xConstantInitByMain); 
  localN = 1<<localM;
  iAutoSyncWrite(&N, &localN, sizeof(N), xIntentionN);

  iAutoSyncRead(&localrootN, &rootN, sizeof(rootN), xIntentionRootN);
  iAutoSyncRead(&localM, &M, sizeof(M), xConstantInitByMain); 
  localrootN = 1<<(localM/2);
  iAutoSyncWrite(&rootN, &localrootN, sizeof(rootN), xIntentionRootN);

  iAutoSyncRead(&localrootN, &rootN, sizeof(rootN), xIntentionRootN);
  iAutoSyncRead(&localP, &P, sizeof(P), xConstantInitByMain);
  rowsperproc = rootN/P;
  if (rowsperproc == 0) {
    printerr("Matrix not large enough. 2**(M/2) must be >= P\n");
    exit(-1);
  }

  line_size = 1 << log2_line_size;
  if (line_size < 2*sizeof(double)) {
    printf("WARNING: Each element is a complex double (%ld bytes)\n",2*sizeof(double));
    printf("  => Less than one element per cache line\n");
    printf("     Computing transpose blocking factor\n");
    factor = (2*sizeof(double)) / line_size;    
    local_num_cache_lines = orig_num_lines / factor;
    iAutoSyncWrite(&num_cache_lines, &local_num_cache_lines, sizeof(num_cache_lines), xConstantInitByMain);
  }
  if (line_size <= 2*sizeof(double)) {
    local_pad_length = 1;
  } else {
    local_pad_length = line_size / (2*sizeof(double));
  }
  iAutoSyncWrite(&pad_length, &local_pad_length, sizeof(pad_length), xConstantInitByMain);

  iAutoSyncRead(&local_pad_length, &pad_length, sizeof(pad_length), xConstantInitByMain);
  iAutoSyncRead(&localrootN, &rootN, sizeof(rootN), xIntentionRootN);
  if (rowsperproc * localrootN * 2 * sizeof(double) >= PAGE_SIZE) {
    pages = (2 * pad_length * sizeof(double) * rowsperproc) / PAGE_SIZE;
    if (pages * PAGE_SIZE != 2 * local_pad_length * sizeof(double) * rowsperproc) {
      pages ++;
    }
    
    local_pad_length = (pages * PAGE_SIZE) / (2 * sizeof(double) * rowsperproc);    
  } else {
    local_pad_length = (PAGE_SIZE - (rowsperproc * rootN * 2 * sizeof(double))) /
                 (2 * sizeof(double) * rowsperproc);
    
    if (local_pad_length * (2 * sizeof(double) * rowsperproc) !=
        (PAGE_SIZE - (rowsperproc * localrootN * 2 * sizeof(double)))) {
      printerr("Padding algorithm unsuccessful\n");
      exit(-1);
    }
  }
  iAutoSyncWrite(&pad_length, &local_pad_length,  sizeof(pad_length), xConstantInitByMain);

  struct GlobalMemory* local_pGlobal;
  double* local_px;
  double* local_ptrans;
  double* local_pumain;
  double* local_pumain2;

  local_pGlobal = (struct GlobalMemory *) G_MALLOC(sizeof(struct GlobalMemory));
  iAutoSyncWrite(&Global, &local_pGlobal, sizeof(Global), xNoSpecialIntention);

  iAutoSyncRead(&localN, &N, sizeof(N), xIntentionN);
  iAutoSyncRead(&localrootN, &rootN, sizeof(rootN), xIntentionRootN);
  iAutoSyncRead(&local_pad_length, &pad_length, sizeof(pad_length), xConstantInitByMain);
  local_px = (double *) G_MALLOC(2*(localN+localrootN*local_pad_length)*sizeof(double)+PAGE_SIZE);
  iAutoSyncWrite(&x, &local_px, sizeof(x), xIntentionX);

  iAutoSyncRead(&localN, &N, sizeof(N), xIntentionN);
  iAutoSyncRead(&localrootN, &rootN, sizeof(rootN), xIntentionRootN);
  iAutoSyncRead(&local_pad_length, &pad_length, sizeof(pad_length), xConstantInitByMain);
  local_ptrans = (double *) G_MALLOC(2*(localN+localrootN*local_pad_length)*sizeof(double)+PAGE_SIZE);
  iAutoSyncWrite(&trans, &local_ptrans, sizeof(trans), xIntentionTrans);

  iAutoSyncRead(&localrootN, &rootN, sizeof(rootN), xIntentionRootN);
  local_pumain = (double *) G_MALLOC(2*localrootN*sizeof(double));
  iAutoSyncWrite(&umain, &local_pumain, sizeof(umain), xIntentionUmain);

  iAutoSyncRead(&localN, &N, sizeof(N), xIntentionN);
  iAutoSyncRead(&localrootN, &rootN, sizeof(rootN), xIntentionRootN);
  iAutoSyncRead(&local_pad_length, &pad_length, sizeof(pad_length), xConstantInitByMain);
  local_pumain2 = (double *) G_MALLOC(2*(localN+localrootN*local_pad_length)*sizeof(double)+PAGE_SIZE);
  iAutoSyncWrite(&umain2, &local_pumain2, sizeof(umain2), xIntentionUmain2);

  long* local_ptranstimes;
  long* local_ptotaltimes;

  local_ptranstimes = (long *) G_MALLOC(P*sizeof(long));
  iAutoSyncWrite(&Global->transtimes, &local_ptranstimes, sizeof(local_ptranstimes), xIntentionTransTimes);

  local_ptotaltimes = (long *) G_MALLOC(P*sizeof(long));
  iAutoSyncWrite(&Global->totaltimes, &local_ptotaltimes, sizeof(local_ptotaltimes), xIntentionTotalTimes);

  iAutoSyncRead(&local_pGlobal, &Global, sizeof(Global), xNoSpecialIntention);
  iAutoSyncRead(&local_px, &x, sizeof(x), xIntentionX);
  iAutoSyncRead(&local_ptrans, &trans, sizeof(trans), xIntentionTrans);
  iAutoSyncRead(&local_pumain, &umain, sizeof(umain), xIntentionUmain);
  iAutoSyncRead(&local_pumain2, &umain2, sizeof(umain2), xIntentionUmain2);
  if (local_pGlobal == NULL) {
    printerr("Could not malloc memory for Global\n");
    exit(-1);
  } else if (local_px == NULL) {
    printerr("Could not malloc memory for x\n");
    exit(-1);
  } else if (local_ptrans == NULL) {
    printerr("Could not malloc memory for trans\n");
    exit(-1);
  } else if (local_pumain == NULL) {
    printerr("Could not malloc memory for umain\n");
    exit(-1);
  } else if (local_pumain2 == NULL) {
    printerr("Could not malloc memory for umain2\n");
    exit(-1);
  }

  iAutoSyncReadToUpdate(&local_px, &x, sizeof(x), xIntentionX);
  local_px = (double *) (((unsigned long) local_px) + PAGE_SIZE - ((unsigned long) local_px) % PAGE_SIZE);
  iAutoSyncUpdate(&x, &local_px, sizeof(x), xIntentionX);

  iAutoSyncReadToUpdate(&local_ptrans, &trans, sizeof(trans), xIntentionTrans);
  local_ptrans = (double *) (((unsigned long) local_ptrans) + PAGE_SIZE - ((unsigned long) local_ptrans) % PAGE_SIZE);
  iAutoSyncUpdate(&trans, &local_ptrans, sizeof(trans), xIntentionTrans);

  iAutoSyncReadToUpdate(&local_pumain2, &umain2, sizeof(umain2), xIntentionUmain2);
  umain2 = (double *) (((unsigned long) umain2) + PAGE_SIZE - ((unsigned long) umain2) % PAGE_SIZE);
  iAutoSyncUpdate(&umain2, &local_pumain2, sizeof(umain2), xIntentionUmain2);

/* In order to optimize data distribution, the data structures x, trans,
   and umain2 have been aligned so that each begins on a page boundary.
   This ensures that the amount of padding calculated by the program is
   such that each processor's partition ends on a page boundary, thus
   ensuring that all data from these structures that are needed by a
   processor can be allocated to its local memory */

/* POSSIBLE ENHANCEMENT:  Here is where one might distribute the x,
   trans, and umain2 data structures across physically distributed
   memories as desired.

   One way to place data is as follows:

   double *base;
   long i;

   i = ((N/P)+(rootN/P)*pad_length)*2;
   base = &(x[0]);
   for (j=0;j<P;j++) {
    Place all addresses x such that (base <= x < base+i) on node j
    base += i;
   }

   The trans and umain2 data structures can be placed in a similar manner.

   */

  iAutoSyncRead(&localN, &N, sizeof(N), xIntentionN);
  iAutoSyncRead(&localP, &P, sizeof(P), xConstantInitByMain);
  iAutoSyncRead(&local_num_cache_lines, &num_cache_lines, sizeof(num_cache_lines), xConstantInitByMain);
  printf("\n");
  printf("FFT with Blocking Transpose\n");
  printf("   %ld Complex Doubles\n",localN);
  printf("   %ld Processors\n",localP);
  if (local_num_cache_lines != orig_num_lines) {
    printf("   %ld Cache lines\n",orig_num_lines);
    printf("   %ld Cache lines for blocking transpose\n",local_num_cache_lines);
  } else {
    printf("   %ld Cache lines\n",local_num_cache_lines);
  }
  printf("   %d Byte line size\n",(1 << log2_line_size)); 
  printf("   %d Bytes per page\n",PAGE_SIZE);
  printf("\n");

  local_global_id = 0;
  iAutoSyncWrite(&Global->id, &local_global_id, sizeof(local_global_id), xNoSpecialIntention);

  iAutoSyncSharedVarAsArg(&x);
  iAutoSyncReadToUpdate(&local_px, &x, sizeof(x), xIntentionX);
  InitX(local_px);                  /* place random values in x */
  iAutoSyncUpdate(&x, &local_px, sizeof(x), xIntentionX);

  iAutoSyncRead(&local_test_result, &test_result, sizeof(test_result), xConstantInitByMain);
  if (local_test_result) {
    iAutoSyncRead(&local_px, &x, sizeof(x), xIntentionX);
    iAutoSyncSharedVarAsArg(&x);
    ck1 = CheckSum(local_px);
  }

  iAutoSyncRead(&local_doprint, &doprint, sizeof(doprint), xConstantInitByMain);
  if (local_doprint) {
    printf("Original data values:\n");
    iAutoSyncRead(&local_px, &x, sizeof(x), xIntentionX);
    iAutoSyncRead(&localN, &N, sizeof(N), xIntentionN);
    iAutoSyncSharedVarAsArg(&x);
    iAutoSyncSharedVarAsArg(&N);
    PrintArray(localN, local_px);
  }

  iAutoSyncRead(&localN, &N, sizeof(N), xIntentionN);  
  iAutoSyncReadToUpdate(&local_pumain, &umain, sizeof(umain), xIntentionUmain);
  iAutoSyncSharedVarAsArg(&N);
  iAutoSyncSharedVarAsArg(&umain);
  InitU(localN, local_pumain);               /* initialize u arrays*/
  iAutoSyncUpdate(&umain, &local_pumain, sizeof(umain), xIntentionUmain);

  iAutoSyncRead(&localrootN, &rootN, sizeof(rootN), xIntentionRootN);
  iAutoSyncRead(&localN, &N, sizeof(N), xIntentionN);
  iAutoSyncReadToUpdate(&local_pumain2, &umain2, sizeof(umain2), xIntentionUmain2);
  iAutoSyncSharedVarAsArg(&N);
  iAutoSyncSharedVarAsArg(&umain2);
  iAutoSyncSharedVarAsArg(&rootN);
  InitU2(localN, local_pumain2, localrootN);
  iAutoSyncUpdate(&umain2, &local_pumain2, sizeof(umain2), xIntentionUmain2);

  /* fire off P processes */

  SPLASH3_ROI_BEGIN();

  CREATE(SlaveStart, P);
  WAIT_FOR_END(P);

  SPLASH3_ROI_END();

  iAutoSyncRead(&local_doprint, &doprint, sizeof(doprint), xConstantInitByMain);
  if (local_doprint) {
    iAutoSyncRead(&local_test_result, &test_result, sizeof(test_result), xConstantInitByMain);
    if (local_test_result) {
      printf("Data values after inverse FFT:\n");
    } else {
      printf("Data values after FFT:\n");
    }
    iAutoSyncRead(&localN, &N, sizeof(N), xIntentionN);
    iAutoSyncRead(&local_px, &x, sizeof(x), xIntentionX);
    iAutoSyncSharedVarAsArg(&N);
    iAutoSyncSharedVarAsArg(&x);
    PrintArray(localN, local_px);
  }


  iAutoSyncRead(&Global->transtimes[0], &transtime, sizeof(transtime), xNoSpecialIntention);
  iAutoSyncRead(&Global->totaltimes[0], &totaltime, sizeof(totaltime), xNoSpecialIntention);
  printf("\n");
  printf("                 PROCESS STATISTICS\n");
  printf("            Computation      Transpose     Transpose\n");
  printf(" Proc          Time            Time        Fraction\n");
  printf("    0        %10ld     %10ld      %8.5f\n",
         totaltime, transtime,
         ((double)transtime)/totaltime);
  iAutoSyncRead(&local_dostats, &dostats, sizeof(dostats), xConstantInitByMain);
  if (local_dostats) {
    transtime2 = transtime;
    avgtranstime = transtime;
    avgcomptime = totaltime;
    maxtotal = totaltime;
    mintotal = totaltime;
    maxfrac = ((double)transtime)/totaltime;
    minfrac = ((double)transtime)/totaltime;
    avgfractime = ((double)transtime)/totaltime;
    
    iAutoSyncRead(&localP, &P, sizeof(P), xConstantInitByMain);
    iAutoSyncRead(&local_pGlobal, &Global, sizeof(Global), xNoSpecialIntention);
    for (i=1;i<localP;i++) {
      if (local_pGlobal->transtimes[i] > transtime) {
        transtime = local_pGlobal->transtimes[i];
      }
      if (local_pGlobal->transtimes[i] < transtime2) {
        transtime2 = local_pGlobal->transtimes[i];
      }
      if (local_pGlobal->totaltimes[i] > maxtotal) {
        maxtotal = local_pGlobal->totaltimes[i];
      }
      if (local_pGlobal->totaltimes[i] < mintotal) {
        mintotal = local_pGlobal->totaltimes[i];
      }
      if (((double)local_pGlobal->transtimes[i])/local_pGlobal->totaltimes[i] > maxfrac) {
        maxfrac = ((double)local_pGlobal->transtimes[i])/local_pGlobal->totaltimes[i];
      }
      if (((double)local_pGlobal->transtimes[i])/local_pGlobal->totaltimes[i] < minfrac) {
        minfrac = ((double)local_pGlobal->transtimes[i])/local_pGlobal->totaltimes[i];
      }
      printf("  %3ld        %10ld     %10ld      %8.5f\n",
             i,local_pGlobal->totaltimes[i],local_pGlobal->transtimes[i],
             ((double)local_pGlobal->transtimes[i])/local_pGlobal->totaltimes[i]);
      avgtranstime += local_pGlobal->transtimes[i];
      avgcomptime += local_pGlobal->totaltimes[i];
      avgfractime += ((double)local_pGlobal->transtimes[i])/local_pGlobal->totaltimes[i];
    }
    printf("  Avg        %10.0f     %10.0f      %8.5f\n",
           ((double) avgcomptime)/P,((double) avgtranstime)/localP,avgfractime/localP);
    printf("  Max        %10ld     %10ld      %8.5f\n",
	   maxtotal,transtime,maxfrac);
    printf("  Min        %10ld     %10ld      %8.5f\n",
	   mintotal,transtime2,minfrac);
  }
  local_pGlobal->starttime = start;
  printf("\n");
  printf("                 TIMING INFORMATION\n");
  printf("Start time                        : %16lu\n",
	  local_pGlobal->starttime);
  printf("Initialization finish time        : %16lu\n",
	  local_pGlobal->initdonetime);
  printf("Overall finish time               : %16lu\n",
	  local_pGlobal->finishtime);
  printf("Total time with initialization    : %16lu\n",
	  local_pGlobal->finishtime-Global->starttime);
  printf("Total time without initialization : %16lu\n",
	  local_pGlobal->finishtime-Global->initdonetime);
  printf("Overall transpose time            : %16ld\n",
         transtime);
  printf("Overall transpose fraction        : %16.5f\n",
         ((double) transtime)/(local_pGlobal->finishtime-local_pGlobal->initdonetime));
  printf("\n");
  
  iAutoSyncRead(&local_test_result, &test_result, sizeof(test_result), xConstantInitByMain);
  if (local_test_result) {    
    iAutoSyncRead(&local_px, &x, sizeof(x), xIntentionX);
    iAutoSyncSharedVarAsArg(&x);
    ck3 = CheckSum(local_px);
    printf("              INVERSE FFT TEST RESULTS\n");
    printf("Checksum difference is %.3f (%.3f, %.3f)\n",
	   ck1-ck3, ck1, ck3);
    if (fabs(ck1-ck3) < 0.001) {
      printf("TEST PASSED\n");
    } else {
      printf("TEST FAILED\n");
    }
  }

  iAutoSyncDestroy();

  MAIN_END;
}


void SlaveStart()
{
  long i;
  long MyNum;
  long NewId;
  double *upriv;
  long initdone;
  long finish;
  long l_transtime=0;
  long MyFirst;
  long MyLast;
  long localrootN;
  double *local_umain;
  long localP;
  long local_dostats;
  long local_test_result;
  long local_timediff;

  /* AutoSync: this use seems ok. Only the granularity of the lock could be adjusted */
  iAutoSyncReadToUpdate(&NewId, &Global->id, sizeof(NewId), xNoSpecialIntention);  
  MyNum = NewId;
  NewId++;
  iAutoSyncUpdate(&Global->id, &NewId, sizeof(NewId), xNoSpecialIntention);

/* POSSIBLE ENHANCEMENT:  Here is where one might pin processes to
   processors to avoid migration */  

  iAutoSyncRead(&localrootN, &rootN, sizeof(rootN), xIntentionRootN);
  upriv = (double *) malloc(2*(localrootN-1)*sizeof(double));
  if (upriv == NULL) {
    fprintf(stderr,"Proc %ld could not malloc memory for upriv\n",MyNum);
    exit(-1);
  }
  
  iAutoSyncRead(&localrootN, &rootN, sizeof(rootN), xIntentionRootN);
  iAutoSyncRead(&local_umain, &umain, sizeof(umain), xIntentionUmain);  
  for (i=0;i<2*(localrootN-1);i++) {
    upriv[i] = local_umain[i];
  } 

  iAutoSyncRead(&localrootN, &rootN, sizeof(rootN), xIntentionRootN);
  iAutoSyncRead(&localP, &P, sizeof(P), xConstantInitByMain);
  MyFirst = localrootN*MyNum/localP;
  MyLast = localrootN*(MyNum+1)/localP;

  iAutoSyncSharedVarAsArg(&x);
  iAutoSyncSharedVarAsArg(&trans);
  iAutoSyncSharedVarAsArg(&umain2);
  TouchArray(x, trans, umain2, upriv, MyFirst, MyLast);

/* POSSIBLE ENHANCEMENT:  Here is where one might reset the
   statistics that one is measuring about the parallel execution */

  iAutoSyncRead(&local_dostats, &dostats, sizeof(dostats), xConstantInitByMain);
  if ((MyNum == 0) || (local_dostats)) {
    CLOCK(initdone);
  }

  /* perform forward FFT */
  iAutoSyncSharedVarAsArg(&M);
  iAutoSyncSharedVarAsArg(&N);
  iAutoSyncSharedVarAsArg(&x);
  iAutoSyncSharedVarAsArg(&trans);
  iAutoSyncSharedVarAsArg(&umain2);
  iAutoSyncSharedVarAsArg(&pad_length);
  iAutoSyncSharedVarAsArg(&test_result);
  iAutoSyncSharedVarAsArg(&dostats);
  FFT1D(1, M, N, x, trans, upriv, umain2, MyNum, &l_transtime, MyFirst,
	MyLast, pad_length, test_result, dostats);  

  iAutoSyncRead(&local_test_result, &test_result, sizeof(test_result), xConstantInitByMain);
  /* perform backward FFT */
  if (local_test_result) {
    iAutoSyncSharedVarAsArg(&M);
    iAutoSyncSharedVarAsArg(&N);
    iAutoSyncSharedVarAsArg(&x);
    iAutoSyncSharedVarAsArg(&trans);
    iAutoSyncSharedVarAsArg(&umain2);
    iAutoSyncSharedVarAsArg(&pad_length);
    iAutoSyncSharedVarAsArg(&test_result);
    iAutoSyncSharedVarAsArg(&dostats);
    FFT1D(-1, M, N, x, trans, upriv, umain2, MyNum, &l_transtime, MyFirst,
	  MyLast, pad_length, test_result, dostats);
  }

  iAutoSyncRead(&local_dostats, &dostats, sizeof(dostats), xConstantInitByMain);
  if ((MyNum == 0) || (local_dostats)) {
    CLOCK(finish);
    iAutoSyncWrite(&Global->transtimes[MyNum], &l_transtime, sizeof(l_transtime), xIntentionTransTimes);

    local_timediff= finish-initdone;
    iAutoSyncWrite(&Global->totaltimes[MyNum], &local_timediff, sizeof(local_timediff), xIntentionTotalTimes);    
  }
  if (MyNum == 0) {
    iAutoSyncWrite(&Global->finishtime, &finish, sizeof(finish), xNoSpecialIntention);      /*ToDo MATHEUS: is this intention correct? */
    iAutoSyncWrite(&Global->initdonetime, &initdone, sizeof(initdone), xNoSpecialIntention); /*ToDo MATHEUS: is this intention correct? */
  }
}


double TouchArray(double *x, double *scratch, double *u, double *upriv, long MyFirst, long MyLast)
{
  long i,j,k;
  double tot = 0.0;
  
  /* touch my data */
  for (j=0;j<2*(rootN-1);j++) {
    tot += upriv[j];
  }
  
  for (j=MyFirst; j<MyLast; j++) {
    k = j * (rootN + pad_length);
    for (i=0;i<rootN;i++) {
      tot += x[2*(k+i)] + x[2*(k+i)+1] +
             scratch[2*(k+i)] + scratch[2*(k+i)+1] +
	     u[2*(k+i)] + u[2*(k+i)+1];
    }
  }
  return tot;
}


double CheckSum(double *x)
{
  long i,j,k;
  double cks;

  cks = 0.0;
  for (j=0; j<rootN; j++) {
    k = j * (rootN + pad_length);
    for (i=0;i<rootN;i++) {
      cks += x[2*(k+i)] + x[2*(k+i)+1];
    }
  }

  return(cks);
}


void InitX(double *x)
{
  long i,j,k;

  srand48(0);
  for (j=0; j<rootN; j++) {
    k = j * (rootN + pad_length);
    for (i=0;i<rootN;i++) {
      x[2*(k+i)] = drand48();
      x[2*(k+i)+1] = drand48();
    }
  }
}


void InitU(long N, double *u)
{
  long q;
  long j;
  long base;
  long n1;

  for (q=0; 1<<q<N; q++) {
    n1 = 1<<q;
    base = n1-1;
    for (j=0; j<n1; j++) {
      if (base+j > rootN-1) {
	return;
      }
      u[2*(base+j)] = cos(2.0*PI*j/(2*n1));
      u[2*(base+j)+1] = -sin(2.0*PI*j/(2*n1));
    }
  }
}


void InitU2(long N, double *u, long n1)
{
  long i,j,k;

  for (j=0; j<n1; j++) {
    k = j*(rootN+pad_length);
    for (i=0; i<n1; i++) {
      u[2*(k+i)] = cos(2.0*PI*i*j/(N));
      u[2*(k+i)+1] = -sin(2.0*PI*i*j/(N));
    }
  }
}


long BitReverse(long M, long k)
{
  long i;
  long j;
  long tmp;

  j = 0;
  tmp = k;
  for (i=0; i<M; i++) {
    j = 2*j + (tmp&0x1);
    tmp = tmp>>1;
  }
  return(j);
}


void FFT1D(long direction, long M, long N, double *x, double *scratch, double *upriv, double *umain2,
           long MyNum, long *l_transtime, long MyFirst, long MyLast, long pad_length, long test_result, long dostats)
{
  long j;
  long m1;
  long n1;
  unsigned long clocktime1;
  unsigned long clocktime2;

  m1 = M/2;
  n1 = 1<<m1;

  printf("iter_num = %lu\n", MyLast - MyFirst);

  if ((MyNum == 0) || (dostats)) {
    CLOCK(clocktime1);
  }

  /* transpose from x into scratch */
  Transpose(n1, x, scratch, MyNum, MyFirst, MyLast, pad_length);

  if ((MyNum == 0) || (dostats)) {
    CLOCK(clocktime2);
    *l_transtime += (clocktime2-clocktime1);
    printf("Step 1: %8lu\n", clocktime2-clocktime1);
  }

  /* do n1 1D FFTs on columns */
  for (j=MyFirst; j<MyLast; j++) {
    FFT1DOnce(direction, m1, n1, upriv, &scratch[2*j*(n1+pad_length)]);
    TwiddleOneCol(direction, n1, j, umain2, &scratch[2*j*(n1+pad_length)], pad_length);
  }

  iAutoSyncProceedOnEvent(xTwiddleDone, P);  

  if ((MyNum == 0) || (dostats)) {
    CLOCK(clocktime1);
    printf("Step 2: %8lu\n", clocktime1-clocktime2);
  }
  /* transpose */
  Transpose(n1, scratch, x, MyNum, MyFirst, MyLast, pad_length);

  if ((MyNum == 0) || (dostats)) {
    CLOCK(clocktime2);
    *l_transtime += (clocktime2-clocktime1);
    printf("Step 3: %8lu\n", clocktime2-clocktime1);
  }

  /* do n1 1D FFTs on columns again */
  for (j=MyFirst; j<MyLast; j++) {
    FFT1DOnce(direction, m1, n1, upriv, &x[2*j*(n1+pad_length)]);
    if (direction == -1)
      Scale(n1, N, &x[2*j*(n1+pad_length)]);
  }

  iAutoSyncProceedOnEvent(xTwiddleDone, P); /* xFFT1Done */  

  if ((MyNum == 0) || (dostats)) {
    CLOCK(clocktime1);
    printf("Step 4: %8lu\n", clocktime1-clocktime2);
  }

  /* transpose back */
  Transpose(n1, x, scratch, MyNum, MyFirst, MyLast, pad_length);

  if ((MyNum == 0) || (dostats)) {
    CLOCK(clocktime2);
    *l_transtime += (clocktime2-clocktime1);
    printf("Step 5: %8lu\n", clocktime2-clocktime1);
  }


  /* copy columns from scratch to x */
  if ((test_result) || (doprint)) 
  {
    /* BARRIER(Global->start, P); */
    iAutoSyncProceedOnEvent(xTwiddleDone, P);  /* xTrasnposeDone */

    for (j=MyFirst; j<MyLast; j++) {
      CopyColumn(n1, &scratch[2*j*(n1+pad_length)], &x[2*j*(n1+pad_length)]);
    }
  }
  /* BARRIER(Global->start, P); */
  iAutoSyncProceedOnEvent(xTwiddleDone, P); /* xFFTDone */
}


void TwiddleOneCol(long direction, long n1, long j, double *u, double *x, long pad_length)
{
  long i;
  double omega_r;
  double omega_c;
  double x_r;
  double x_c;

  for (i=0; i<n1; i++) {
    omega_r = u[2*(j*(n1+pad_length)+i)];
    omega_c = direction*u[2*(j*(n1+pad_length)+i)+1];
    x_r = x[2*i];
    x_c = x[2*i+1];
    x[2*i] = omega_r*x_r - omega_c*x_c;
    x[2*i+1] = omega_r*x_c + omega_c*x_r;
  }
}


void Scale(long n1, long N, double *x)
{
  long i;

  for (i=0; i<n1; i++) {
    x[2*i] /= N;
    x[2*i+1] /= N;
  }
}


void Transpose(long n1, double *src, double *dest, long MyNum, long MyFirst, long MyLast, long pad_length)
{
  long i;
  long j;
  long k;
  long l;
  long m;
  long blksize;
  long numblks;
  long firstfirst;
  long h_off;
  long v_off;
  long v;
  long h;
  long n1p;
  long row_count;
  long iter_num = 0;

  blksize = MyLast-MyFirst;
  numblks = (2*blksize)/num_cache_lines;
  if (numblks * num_cache_lines != 2 * blksize) {
    numblks ++;
  }
  blksize = blksize / numblks;
  firstfirst = MyFirst;
  row_count = n1/P;
  n1p = n1+pad_length;
  for (l=MyNum+1;l<P;l++) {
    v_off = l*row_count;
    for (k=0; k<numblks; k++) {
      h_off = firstfirst;
      for (m=0; m<numblks; m++) {
        for (i=0; i<blksize; i++) {
	  v = v_off + i;
          for (j=0; j<blksize; j++) {
            iter_num ++;
	    h = h_off + j;
            dest[2*(h*n1p+v)] = src[2*(v*n1p+h)];
            dest[2*(h*n1p+v)+1] = src[2*(v*n1p+h)+1];
          }
        }
	h_off += blksize;
      }
      v_off+=blksize;
    }
  }
  printf("Transpose: iter_num = %lu\n", iter_num);

  for (l=0;l<MyNum;l++) {
    v_off = l*row_count;
    for (k=0; k<numblks; k++) {
      h_off = firstfirst;
      for (m=0; m<numblks; m++) {
        for (i=0; i<blksize; i++) {
	  v = v_off + i;
          for (j=0; j<blksize; j++) {
            h = h_off + j;
            dest[2*(h*n1p+v)] = src[2*(v*n1p+h)];
            dest[2*(h*n1p+v)+1] = src[2*(v*n1p+h)+1];
          }
        }
	h_off += blksize;
      }
      v_off+=blksize;
    }
  }

  v_off = MyNum*row_count;
  for (k=0; k<numblks; k++) {
    h_off = firstfirst;
    for (m=0; m<numblks; m++) {
      for (i=0; i<blksize; i++) {
        v = v_off + i;
        for (j=0; j<blksize; j++) {
          h = h_off + j;
          dest[2*(h*n1p+v)] = src[2*(v*n1p+h)];
          dest[2*(h*n1p+v)+1] = src[2*(v*n1p+h)+1];
	}
      }
      h_off += blksize;
    }
    v_off+=blksize;
  }
}


void CopyColumn(long n1, double *src, double *dest)
{
  long i;

  for (i=0; i<n1; i++) {
    dest[2*i] = src[2*i];
    dest[2*i+1] = src[2*i+1];
  }
}


void Reverse(long N, long M, double *x)
{
  long j, k;

  for (k=0; k<N; k++) {
    j = BitReverse(M, k);
    if (j > k) {
      SWAP_VALS(x[2*j], x[2*k]);
      SWAP_VALS(x[2*j+1], x[2*k+1]);
    }
  }
}


void FFT1DOnce(long direction, long M, long N, double *u, double *x)
{
  long j;
  long k;
  long q;
  long L;
  long r;
  long Lstar;
  double *u1;
  double *x1;
  double *x2;
  double omega_r;
  double omega_c;
  double tau_r;
  double tau_c;
  double x_r;
  double x_c;
  long iter_num = 0;

  Reverse(N, M, x);

  for (q=1; q<=M; q++) {
    L = 1<<q; r = N/L; Lstar = L/2;
    u1 = &u[2*(Lstar-1)];
    for (k=0; k<r; k++) {
      x1 = &x[2*(k*L)];
      x2 = &x[2*(k*L+Lstar)];
      for (j=0; j<Lstar; j++) {
        iter_num ++;
	omega_r = u1[2*j];
        omega_c = direction*u1[2*j+1];
	x_r = x2[2*j];
        x_c = x2[2*j+1];
	tau_r = omega_r*x_r - omega_c*x_c;
	tau_c = omega_r*x_c + omega_c*x_r;
	x_r = x1[2*j];
        x_c = x1[2*j+1];
	x2[2*j] = x_r - tau_r;
	x2[2*j+1] = x_c - tau_c;
	x1[2*j] = x_r + tau_r;
	x1[2*j+1] = x_c + tau_c;
      }
    }
  }
  if(is_output == 1){
      printf("FFt1DOnce: iter_num = %lu\n", iter_num);
      is_output = 0;
  }
}


void PrintArray(long N, double *x)
{
  long i, j, k;

  for (i=0; i<rootN; i++) {
    k = i*(rootN+pad_length);
    for (j=0; j<rootN; j++) {
      printf(" %4.2f %4.2f", x[2*(k+j)], x[2*(k+j)+1]);
      if (i*rootN+j != N-1) {
        printf(",");
      }
      if ((i*rootN+j+1) % 8 == 0) {
        printf("\n");
      }
    }
  }
  printf("\n");
  printf("\n");
}


void printerr(char *s)
{
  fprintf(stderr,"ERROR: %s\n",s);
}


long log_2(long number)
{
  long cumulative = 1, out = 0, done = 0;

  while ((cumulative < number) && (!done) && (out < 50)) {
    if (cumulative == number) {
      done = 1;
    } else {
      cumulative = cumulative * 2;
      out ++;
    }
  }

  if (cumulative == number) {
    return(out);
  } else {
    return(-1);
  }
}
