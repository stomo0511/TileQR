UNAME = $(shell uname)
MYNAME = $(shell whoami)
ifeq ($(UNAME),Linux)
	ifeq ($(MYNAME),v74001)
# Wisteria (Tokyo)
		CXX = FCCpx
		CXXFLAGS = -Nclang -Nlibomp -Kfast -Kopenmp
		BLAS_LIBS = -SSL2BLAMP
	else
# AOBA
		ifeq ($(MYNAME),y01505)
			MKL_LIB_DIR = $(MKLROOT)/lib/intel64
			MKL_INC_DIR = $(MKLROOT)/include
			MKL_LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

			PLASMA_DIR = /uhome/y01505/plasma
			PLASMA_LIB_DIR = $(PLASMA_DIR)/lib
			PLASMA_INC_DIR = $(PLASMA_DIR)/include
			PLASMA_LIBS = -L$(PLASMA_LIB_DIR) -lplasma_core_blas

			CXX = clang++
			CXXFLAGS =  -DMKL -I$(MKL_INC_DIR) -I$(PLASMA_INC_DIR) -fopenmp
			BLAS_LIBS = -L$(MKL_LIB_DIR) $(MKL_LIBS) -lpthread -lm -ldl
		else
# Linux server
			MKL_LIB_DIR = $(MKLROOT)/lib/intel64
			MKL_INC_DIR = $(MKLROOT)/include
			MKL_LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
		
			PLASMA_DIR = /opt/plasma
			PLASMA_LIB_DIR = $(PLASMA_DIR)/lib
			PLASMA_INC_DIR = $(PLASMA_DIR)/include
			PLASMA_LIBS = -L$(PLASMA_LIB_DIR) -lplasma_core_blas

			CXX = mpic++
			CXXFLAGS = -DMKL -I$(MKL_INC_DIR) -I$(PLASMA_INC_DIR) -fopenmp
			BLAS_LIBS = -L$(MKL_LIB_DIR) $(MKL_LIBS) -lpthread -lm -ldl
		endif
	endif
endif
ifeq ($(UNAME),Darwin)
# M1 MacBook Pro
	PLASMA_DIR = /opt/plasma
	PLASMA_INC_DIR = $(PLASMA_DIR)/include
	PLASMA_LIB_DIR = $(PLASMA_DIR)/lib
	PLASMA_LIBS = -L$(PLASMA_LIB_DIR) -lcoreblas

	OPENBLAS_DIR = /opt/homebrew/opt/openblas
	BLAS_INC_DIR = $(OPENBLAS_DIR)/include
	BLAS_LIB_DIR = $(OPENBLAS_DIR)/lib
	BLAS_LIBS = -L$(BLAS_LIB_DIR) -lopenblas -lgomp

	CXX = g++-13
	CC = gcc-13
	CXXFLAGS = -fopenmp -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR)
endif

OBJS = TileQR.o CoreBlas.o detach.o

# for Performance evaluation
# CXXFLAGS += -O3

# for Debug
CXXFLAGS += -DDEBUG -g

# for Trace
# CXXFLAGS += -DTRACE
# OBJS += trace.o

all: TileQR

TileQR : $(OBJS)
	$(CXX) $(CFLAGS) -o $@ $(OBJS) $(PLASMA_LIBS) $(BLAS_LIBS) 

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

detach.o: detach.cpp mpi-detach.h
	$(CXX) $(CXX_FLAGS) -DOMPI_SKIP_MPICXX=1 -c detach.cpp -g

trace.o: trace.c
	$(CXX) -O3 -c -o $@ $<

clean:
	rm -f *.o TileQR
