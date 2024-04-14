# Cross-compile changes by Thomas E. Hansen (CodingCellist) 2021-03-08
TARGET = FFT
OBJS = fft_auto_sync2.o

CC := $(TOOLCHAIN_PREFIX)gcc
CFLAGS := -O2 -pthread -D_XOPEN_SOURCE=500 -D_POSIX_C_SOURCE=200112 -std=c11 -g -fno-strict-aliasing
LDFLAGS := -lm

EXECUTE_FFT := ./FFT -m28 -p4 -n512 -l6 -t

# BASEDIR needs to be set to the same directory as this Makefile
# BASEDIR := $(HOME)/Splash-3/codes
BASEDIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
MACROS := pthread.m4.stougie
M4 := m4 -Ulen -Uindex

ifndef BASEDIR
$(error Please define the BASEDIR variable)
endif

.SECONDARY: fft.c

x = *

$(TARGET): $(OBJS) $(EXTRADEPS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

auto_sync: 
	$(M4) $(MACROS) fft_auto_sync.c > fft_auto_sync2.c
	echo "M4 Macros applied"

build: auto_sync
	# Copy AutoSync files
	cd .. && make parse C_FILE=../10_SPLASH_FFT/fft_auto_sync2.c && make generate_code C_FILE=../10_SPLASH_FFT/fft_auto_sync2.c
	cp ../05_Workspace/temp.c ../../SPLASH_3/codes/kernels/fft_auto_sync/fft_auto_sync2.c
	cp ../05_Workspace/_AutoSync.h ../../SPLASH_3/codes/kernels/fft_auto_sync/_AutoSync.h
	cp ../05_Workspace/_AutoSync.c ../../SPLASH_3/codes/kernels/fft_auto_sync/_AutoSync.c
	cd ../../SPLASH_3/codes/kernels/fft_auto_sync && gcc -c _AutoSync.c -o AutoSync.o && ar rcs AutoSyncLib.a AutoSync.o
	cd ../../SPLASH_3/codes/kernels/fft_auto_sync && make

run: build 	
	cd ../../SPLASH_3/codes/kernels/fft_auto_sync && $(EXECUTE_FFT)
	@echo "How to run FFT:" $(EXECUTE_FFT)
	
valgrind: build
	cd ../../SPLASH_3/codes/kernels/fft_auto_sync && valgrind $(EXECUTE_FFT)

clean:
	rm -rf *.c *.h *.o $(TARGET)

%.h: %.h.in
	$(M4) $(MACROS) $< > $@

%.c: %.c.in
	$(M4) $(MACROS) $< > $@

%.o: %.c $(EXTRADEPS)
	$(CC) -c $(CFLAGS) $<