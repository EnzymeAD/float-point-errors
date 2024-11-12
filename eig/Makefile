HOME ?= /home/sbrantq
ENZYME_PATH ?= $(HOME)/sync/Enzyme/build/Enzyme/ClangEnzyme-15.so
CLANG_PATH ?= $(HOME)/llvms/llvm15/build/bin

SHELL = /bin/sh

.SUFFIXES: .cpp .o

CXX = $(CLANG_PATH)/clang++
LLC = $(CLANG_PATH)/llc

SOURCES = eig.cpp
OBJECTS = $(SOURCES:.cpp=.o)

CXXFLAGS = -I$(HOME)/include -L$(HOME)/lib \
           -I/usr/include/c++/13 \
           -I/usr/include/x86_64-linux-gnu/c++/13 \
           -L/usr/lib/gcc/x86_64-linux-gnu/13 \
           -O3 -I. -Wall \
           -fno-exceptions \
           -mllvm -enzyme-loose-types \
           -fplugin=$(ENZYME_PATH) \
           -mllvm -enzyme-inline \
           -ffast-math -fno-finite-math-only

all: eig.exe eig-logged.exe eig-fpopt.exe

eig.exe: $(SOURCES)
	@echo "Building original executable: $@"
	time $(CXX) $(SOURCES) $(CXXFLAGS) -lm -o $@

eig-logged.exe: $(SOURCES) fp-logger.o
	@echo "Building logged executable: $@"
	time $(CXX) $(SOURCES) fp-logger.o $(CXXFLAGS) -lm -o $@ -DLOGGED

fp-logger.o: fp-logger.cpp
	@echo "Compiling fp-logger.o"
	time $(CXX) -c fp-logger.cpp $(CXXFLAGS) -o $@

eig-fpopt.exe: eig-fpopt.o
	@echo "Building optimized executable: $@"
	time $(CXX) eig-fpopt.o $(CXXFLAGS) -lm -o $@

eig.txt: eig-logged.exe
	@echo "Running logged executable: $<"
	./$< > $@

eig-fpopt.o: eig.cpp eig.txt
	@echo "Compiling optimized object: $@"
	time $(CXX) -c $(SOURCES) $(CXXFLAGS) \
		-mllvm -enzyme-enable-fpopt \
		-mllvm -fpopt-log-path=eig.txt \
		-mllvm -fpopt-enable-solver \
		-mllvm -fpopt-enable-pt \
		-mllvm -fpopt-target-func-regex=eig \
		-mllvm -fpopt-comp-cost-budget=0 \
		-mllvm -fpopt-num-samples=1024 \
		-mllvm -fpopt-cost-dom-thres=0.0 \
		-mllvm -fpopt-acc-dom-thres=0.0 \
		-mllvm -enzyme-print-fpopt \
		-mllvm -fpopt-cache-path=cache \
		-mllvm -herbie-timeout=3000 \
		-mllvm -herbie-num-iters=4 \
        -mllvm -herbie-num-threads=12 \
		-mllvm -fpopt-cost-model-path=cm.csv \
		-o $@

clean:
	@echo "Cleaning up..."
	/bin/rm -f *.o *~ *.exe
	/bin/rm -rf *.dSYM
	/bin/rm -rf herbie_input_* herbie_output_*

.PHONY: all clean
