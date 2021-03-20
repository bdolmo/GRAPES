
CC=g++
CFLAGS= -std=c++11 -pthread
RM=rm

SEQLIB= ../SeqLib
INC= -I ../htslib -I ../SeqLib/htslib/htslib -I ../SeqLib -L ../SeqLib/htslib
LIB= -L ../SeqLib/htslib -L ../SeqLib/bin -Lhtslib
LDFLAGS= ../SeqLib/bin/libseqlib.a


%.o: src/%.c*
	@echo ' Building file: $<'
	@echo ' Invoking: GCC C++ Compiler'
	$(CC) -I${SEQLIB} $(INC) $(LIB) -O2 -lhts -lboost_system -lboost_filesystem -c -fopenmp -std=c++11 -o "$@" "$<"
	@echo ' Finished building: $<'
	@echo ' '

all: AnnFeat
install: AnnFeat

OBJS += \
AnnFeat.o \
Utils.o \
varCall.o \


AnnFeat: $(OBJS)
	@echo ' Building target: $@'
	@echo ' Invoking: GCC C++ Linker' 
	$(CC) -o $@ $(CFLAGS) $(OBJS) $(LDFLAGS) $(INC) -O2 -lhts -lboost_program_options -lboost_system -lboost_filesystem -lz -fopenmp
	export LD_LIBRARY_PATH=../SeqLib/htslib
	@echo ' Finished building target: $@'
	@echo ' '

clean:
	@echo ' Removing object files: $(OBJS)'
	$(RM) *.o
	
