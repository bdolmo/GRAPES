
CC=g++
CFLAGS= -std=c++11 -pthread
RM=rm

SEQLIB= ../SeqLib
INC= -I ../htslib -I ../SeqLib/htslib/htslib -I ../SeqLib -L ../SeqLib/htslib
LIB= -L ../SeqLib/htslib -L ../SeqLib/bin -Lhtslib


%.o: src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CC) -I${SEQLIB} $(INC) $(LIB) -lhts -lseqlib -c -std=c++11 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

all: extractOfftarget
install: extractOfftarget

OBJS += \
extractOfftarget.o \

extractOfftarget: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker' 
	$(CC) -o $@ $(CFLAGS) $(OBJS) $(LIB)  $(INC) -I${SEQLIB} -lhts -lseqlib -lboost_program_options
	export LD_LIBRARY_PATH=../SeqLib/htslib
	@echo 'Finished building target: $@'
	@echo ' '

clean:
	@echo 'Removing object files: $(OBJS)'
	$(RM) *.o
