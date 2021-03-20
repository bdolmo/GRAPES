
CC=g++
CFLAGS= -std=c++11 -pthread
RM=rm

SEQLIB= ../SeqLib
INC= -I ../htslib -I ../SeqLib/htslib/htslib -I ../SeqLib -L ../SeqLib/htslib
LIB= -L ../SeqLib/htslib -L ../SeqLib/bin -Lhtslib


%.o: src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CC) -I${SEQLIB} $(INC) $(LIB) -O2 -lhts -lseqlib -c -std=c++11 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

all: TargetDepth
install: TargetDepth

OBJS += \
TargetDepth.o \
ExtractDepth.o \

TargetDepth: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker' 
	$(CC) -o $@ $(CFLAGS) $(OBJS) $(LIB) $(INC) -I${SEQLIB} -O2 -lhts -lseqlib -lboost_program_options
	export LD_LIBRARY_PATH=../SeqLib/htslib
	@echo 'Finished building target: $@'
	@echo ' '

clean:
	@echo 'Removing object files: $(OBJS)'
	$(RM) *.o
