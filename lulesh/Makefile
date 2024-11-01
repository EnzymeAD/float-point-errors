#default build suggestion of MPI + OPENMP with gcc on Livermore machines you might have to change the compiler name

# ENZYME_PATH ?= /home/ubuntu/Enzyme/enzyme/build/Enzyme/ClangEnzyme-14.so
# CLANG_PATH ?= /home/ubuntu/llvm-project/build/bin

# ENZYME_PATH ?= /home/wmoses/git/Enzyme/enzyme/build13D/Enzyme/ClangEnzyme-13.so
# CLANG_PATH ?= /mnt/sabrent/wmoses/llvm13/buildall/bin

HOME ?= /home/sbrantq
ENZYME_PATH ?= $(HOME)/sync/Enzyme/build/Enzyme/ClangEnzyme-15.so
CLANG_PATH ?= $(HOME)/llvms/llvm15/build/bin

#ENZYME_PATH ?= /home/ubuntu/Enzyme/enzyme/build/Enzyme/ClangEnzyme-14.so
#CLANG_PATH ?= /home/ubuntu/llvm-project/build/bin

OPENMP_PATH ?= $(CLANG_PATH)/../projects/openmp/runtime/src
MPI_PATH ?= /usr/lib/x86_64-linux-gnu/openmpi/lib
OPENMP_LIB ?= $(CLANG_PATH)/../lib/libomp.so

SHELL = /bin/sh
.SUFFIXES: .cc .o

CXX = $(CLANG_PATH)/clang++
LLC = $(CLANG_PATH)/llc

SOURCES2.0 = \
	lulesh.cc \
	lulesh-viz.cc \
	lulesh-util.cc \
	lulesh-init.cc
OBJECTS2.0 = $(SOURCES2.0:.cc=.o)

#Default build suggestions with OpenMP for g++
CXXFLAGS = -I$(HOME)/include -L$(HOME)/lib -I/usr/include/c++/13 -I/usr/include/x86_64-linux-gnu/c++/13 -L/usr/lib/gcc/x86_64-linux-gnu/13  
CXXFLAGS += -O3 -I. -Wall -I $(OPENMP_PATH) -L $(MPI_PATH) -fno-exceptions -mllvm -enzyme-loose-types -fplugin=$(ENZYME_PATH) -lmpi -mllvm -enzyme-inline -ffast-math -fno-finite-math-only 
# CXXFLAGS += -mllvm -attributor-max-iterations=128 -mllvm -capture-tracking-max-uses-to-explore=256 -mllvm -memdep-block-scan-limit=70000 -mllvm -dse-memoryssa-walklimit=70000 
# -mllvm -enzyme-print -Rpass=enzyme
# enzyme-print-perf
# CXXFLAGS += -mllvm -debug-only=openmp-opt -mllvm -openmp-opt-print-module
# -mllvm -enzyme-attributor=0
# -mllvm -enzyme-print -mllvm -enzyme-print-type


# all: ser-mpi-gradient.exe

all: ser-single-gradient.exe

# all: ser-single-forward.exe omp-single-forward.exe ompM-single-forward.exe ser-mpi-forward.exe omp-mpi-forward.exe ompM-mpi-forward.exe ser-single-gradient.exe omp-single-gradient.exe ompM-single-gradient.exe ser-mpi-gradient.exe omp-mpi-gradient.exe ompM-mpi-gradient.exe ompOpt-single-forward.exe ompOpt-mpi-forward.exe ompOpt-single-gradient.exe ompOpt-mpi-gradient.exe ser-single-gradientnomat.exe

# all: omp-single-forward.exe omp-single-gradient.exe

simple.ll: simple.cpp
	time $(CXX) simple.cpp $(CXXFLAGS) -S -emit-llvm -o $@ 

ser-single-forward.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ -DFORWARD=1

omp-single-forward.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -DFORWARD=1 -mllvm -enzyme-omp-opt=0

ompM-single-forward.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=1 -DUSE_MPI=0 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -DFORWARD=1 -mllvm -enzyme-omp-opt=0

ompOpt-single-forward.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=1 -DUSE_MPI=0 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -DFORWARD=1 -mllvm -enzyme-omp-opt=1

ser-mpi-forward.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=1 $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ -DFORWARD=1

omp-mpi-forward.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=1 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -DFORWARD=1 -mllvm -enzyme-omp-opt=0

ompM-mpi-forward.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=1 -DUSE_MPI=1 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -DFORWARD=1 -mllvm -enzyme-omp-opt=0

ompOpt-mpi-forward.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=1 -DUSE_MPI=1 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -DFORWARD=1 -mllvm -enzyme-omp-opt=1

ser-single-gradient.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 $(SOURCES2.0) fp-logger.cpp $(CXXFLAGS) -lm -o $@

ser-single-gradientnomat.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ -mllvm -enzyme-rematerialize=0

omp-single-gradient.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -mllvm -enzyme-omp-opt=0

ompM-single-gradient.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=1 -DUSE_MPI=0 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -mllvm -enzyme-omp-opt=0

ompOpt-single-gradient.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=1 -DUSE_MPI=0 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -mllvm -enzyme-omp-opt=1

ser-mpi-gradient.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=1 $(SOURCES2.0) $(CXXFLAGS) -lm -o $@

omp-mpi-gradient.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=1 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -mllvm -enzyme-omp-opt=0

ompM-mpi-gradient.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=1 -DUSE_MPI=1 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -mllvm -enzyme-omp-opt=0

ompOpt-mpi-gradient.exe: $(SOURCES2.0)
	time $(CXX) -DOMP_MERGE=1 -DUSE_MPI=1 -fopenmp $(SOURCES2.0) $(CXXFLAGS) -lm -o $@ $(OPENMP_LIB) -mllvm -enzyme-omp-opt=1

clean:
	/bin/rm -f *.o *~ *.exe
	/bin/rm -rf *.dSYM
	/bin/rm -rf herbie_input_* herbie_output_*

lulesh-fpopt.o:
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -DFORWARD=1 -c lulesh.cc $(CXXFLAGS) \
		-mllvm -enzyme-enable-fpopt \
		-mllvm -fpopt-log-path=lulesh.txt \
		-mllvm -fpopt-enable-solver \
		-mllvm -fpopt-target-func-regex=LagrangeLeapFrog \
		-mllvm --fpopt-enable-pt \
		-mllvm -fpopt-comp-cost-budget=0 \
		-mllvm --fpopt-num-samples=1000 \
		-mllvm --fpopt-cost-model-path=cm.csv -o $@
	# time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -DFORWARD=1 lulesh.cc $(CXXFLAGS) -mllvm -enzyme-enable-fpopt -mllvm -fpopt-log-path=lulesh.txt -mllvm -enzyme-print-herbie -mllvm -enzyme-print-fpopt -S -emit-llvm -o lulesh-fpopt.ll # Doesn't work
	# time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -DFORWARD=1 lulesh.cc $(CXXFLAGS) -mllvm -enzyme-enable-fpopt -mllvm -fpopt-log-path=lulesh.txt -mllvm -enzyme-print-herbie -mllvm -enzyme-print-fpopt -mllvm -herbie-disable-taylor -mllvm -fpopt-enable-solver -mllvm -fpopt-solver-type=dp -S -emit-llvm -o lulesh-fpopt.ll  # Works
	# time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -DFORWARD=1 lulesh.cc $(CXXFLAGS) -mllvm -enzyme-enable-fpopt -mllvm -fpopt-log-path=lulesh.txt -mllvm -enzyme-print-herbie -mllvm -enzyme-print-fpopt -mllvm -herbie-disable-taylor -mllvm -herbie-disable-branch-expr -mllvm -fpopt-enable-solver -mllvm -fpopt-solver-type=dp -S -emit-llvm -o lulesh-fpopt.ll  # Works
	# time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -DFORWARD=1 lulesh.cc $(CXXFLAGS) -mllvm -enzyme-enable-fpopt -mllvm -fpopt-log-path=lulesh.txt -mllvm -enzyme-print-herbie -mllvm -enzyme-print-fpopt -mllvm -herbie-disable-taylor -mllvm -herbie-disable-regime -mllvm -herbie-disable-branch-expr -mllvm -fpopt-enable-solver -mllvm -fpopt-solver-type=dp -S -emit-llvm -o lulesh-fpopt.ll  # Works
	# time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -DFORWARD=1 lulesh.cc $(CXXFLAGS) -mllvm -enzyme-enable-fpopt -mllvm -fpopt-log-path=lulesh-fake.txt -mllvm -enzyme-print-herbie -mllvm -enzyme-print-fpopt -S -emit-llvm -o lulesh-fpopt.ll
	# $(LLC) lulesh-fpopt.ll --filetype=obj -O3 -o lulesh-fpopt.o

ser-single-forward-fpopt.exe: lulesh-fpopt.o
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -c lulesh-viz.cc lulesh-util.cc lulesh-init.cc $(CXXFLAGS)
	time $(CXX) lulesh-fpopt.o lulesh-viz.o lulesh-util.o lulesh-init.o $(CXXFLAGS) -lm -o $@ -DFORWARD=1

lulesh-others.ll:
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -DFORWARD=1 -S -emit-llvm lulesh-viz.cc $(CXXFLAGS) -o lulesh-viz.ll
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -DFORWARD=1 -S -emit-llvm lulesh-util.cc $(CXXFLAGS) -o lulesh-util.ll
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 -DFORWARD=1 -S -emit-llvm lulesh-init.cc $(CXXFLAGS) -o lulesh-init.ll

