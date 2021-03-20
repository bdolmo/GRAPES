
CC=g++
CFLAGS= -std=c++11 -pthread -g
RM=rm

SEQLIB= ../SeqLib
INC= -I ../htslib -I ../SeqLib/htslib/htslib -I ../SeqLib -L ../SeqLib/htslib
LIB= -L ../SeqLib/htslib -L ../SeqLib/bin -Lhtslib


%.o: src/%.c*
	@echo ' Building file: $<'
	@echo ' Invoking: GCC C++ Compiler'
	$(CC) -I${SEQLIB} $(INC) $(LIB) -O2 -lhts -lseqlib -lbwa -lfml -lboost_system -lboost_filesystem  -lboost_timer -fopenmp -c -std=c++11 -o "$@" "$<"
	@echo ' Finished building: $<'
	@echo ' '

all: GRAPES
install: GRAPES


OBJS += \
grapes.o \
Assembler.o \
Aligner.o \
alignContig.o \
callSV.o\
cnv.o \
clusterSC.o \
clusterDiscordant.o \
dna2bit.o \
excludeRegion.o \
largeSV.o\
Trimmer.o \
Utils.o \
smallSV.o\
SV.o \
ssw_cpp.o \
ssw.o\
varCall.o\
hashAligner.o\


GRAPES: $(OBJS)
	@echo ' Building target: $@'
	@echo ' Invoking: GCC C++ Linker' 
	$(CC) -o $@ $(CFLAGS) $(OBJS) $(LIB)  $(INC) -I${SEQLIB} -O2 -lhts -lseqlib -lbwa -lfml -lboost_program_options -lboost_timer  -lboost_system -lboost_filesystem -fopenmp -lz 
	export LD_LIBRARY_PATH=../SeqLib/htslib
	@echo ' Finished building target: $@'
	@echo ' '

clean:
	@echo ' Removing object files: $(OBJS)'
	$(RM) *.o


