#default build suggestion of MPI + OPENMP with gcc on Livermore machines you might have to change the compiler name

# ENZYME_PATH ?= /home/ubuntu/Enzyme/enzyme/build/Enzyme/ClangEnzyme-14.so
# CLANG_PATH ?= /home/ubuntu/llvm-project/build/bin

# ENZYME_PATH ?= /home/wmoses/git/Enzyme/enzyme/build13D/Enzyme/ClangEnzyme-13.so
# CLANG_PATH ?= /mnt/sabrent/wmoses/llvm13/buildall/bin

ENZYME_PATH ?= /home/wmoses/git/Enzyme/enzyme/build12/Enzyme/ClangEnzyme-12.so
CLANG_PATH ?= /home/wmoses/llvms/llvm12/buildD/bin

#ENZYME_PATH ?= /home/ubuntu/Enzyme/enzyme/build/Enzyme/ClangEnzyme-14.so
#CLANG_PATH ?= /home/ubuntu/llvm-project/build/bin

OPENMP_PATH ?= $(CLANG_PATH)/../projects/openmp/runtime/src
MPI_PATH ?= /usr/lib/x86_64-linux-gnu/openmpi/include
OPENMP_LIB ?= $(CLANG_PATH)/../lib/libomp.so

SHELL = /bin/sh
.SUFFIXES: .cc .o

MPI_INC = /opt/local/include/openmpi
MPI_LIB = /opt/local/lib

CXX = $(CLANG_PATH)/clang++

SOURCES2.0 = \
	lulesh.cc \
	lulesh-viz.cc \
	lulesh-util.cc \
	lulesh-init.cc
OBJECTS2.0 = $(SOURCES2.0:.cc=.o)

#Default build suggestions with OpenMP for g++
CXXFLAGS = -I/usr/include/c++/11 -I/usr/include/x86_64-linux-gnu/c++/11 -L/usr/lib/gcc/x86_64-linux-gnu/11 -O3 -I. -Wall -I $(OPENMP_PATH) -I $(MPI_PATH) -fno-exceptions -mllvm -enzyme-loose-types -fplugin=$(ENZYME_PATH) -lmpi -mllvm -attributor-max-iterations=128 -mllvm -capture-tracking-max-uses-to-explore=256 -ffast-math -mllvm -memdep-block-scan-limit=70000 -mllvm -dse-memoryssa-walklimit=70000 -mllvm -enzyme-inline 
# -mllvm -enzyme-print -Rpass=enzyme
# enzyme-print-perf
# CXXFLAGS += -mllvm -debug-only=openmp-opt -mllvm -openmp-opt-print-module
# -mllvm -enzyme-attributor=0
# -mllvm -enzyme-print -mllvm -enzyme-print-type


all: ser-single-forward.exe omp-single-forward.exe ompM-single-forward.exe ser-mpi-forward.exe omp-mpi-forward.exe ompM-mpi-forward.exe ser-single-gradient.exe omp-single-gradient.exe ompM-single-gradient.exe ser-mpi-gradient.exe omp-mpi-gradient.exe ompM-mpi-gradient.exe ompOpt-single-forward.exe ompOpt-mpi-forward.exe ompOpt-single-gradient.exe ompOpt-mpi-gradient.exe ser-single-gradientnomat.exe

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
	time $(CXX) -DOMP_MERGE=0 -DUSE_MPI=0 $(SOURCES2.0) $(CXXFLAGS) -lm -o $@

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

