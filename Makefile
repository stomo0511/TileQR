UNAME = $(shell uname)
MYNAME = $(shell whoami)
ifeq ($(UNAME),Linux)
	ifeq ($(MYNAME),v74001)
# Wisteria (Tokyo)
		CXX = FCCpx
		CXXFLAGS = -Nclang -Nlibomp -Kfast -Kopenmp
		BLAS_LIBS = -SSL2BLAMP
	else
# Linux server
		MKL_LIB_DIR = $(MKLROOT)/lib/intel64
		MKL_INC_DIR = $(MKLROOT)/include
		MKL_LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
		
		PLASMA_DIR = /opt/plasma
		PLASMA_LIB_DIR = $(PLASMA_DIR)/lib
		PLASMA_INC_DIR = $(PLASMA_DIR)/include
		PLASMA_LIBS = -lplasma_core_blas

		CXX = g++-9
		CXXFLAGS = -DMKL -I$(MKL_INC_DIR) -I$(PLASMA_INC_DIR) -fopenmp
		BLAS_LIBS = -L$(MKL_LIB_DIR) $(MKL_LIBS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -lpthread -lm -ldl
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

OBJS = SuperM.o Matrix.o Tile.o TMatrix.o CoreBlas.o Check_Accuracy.o TileQR.o
RL_OBJS = $(OBJS) Right_Looking.o
LL_OBJS = $(OBJS) Left_Looking.o

# for Performance evaluation
CXXFLAGS += -O3

# for Debug
CXXFLAGS += -DDEBUG -g

# for Trace
# CXXFLAGS += -DTRACE
# OBJS += trace.o

all: RL LL

RL : $(RL_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(RL_OBJS) $(PLASMA_LIBS) $(BLAS_LIBS) 

LL : $(LL_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(LL_OBJS) $(PLASMA_LIBS) $(BLAS_LIBS)

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

trace.o: trace.c
	$(CXX) -O3 -c -o $@ $<

clean:
	rm -f *.o RL LL
