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
		PLASMA_LIBS = -lcoreblas

		CXX = g++-8
		CXXFLAGS = -DMKL -I$(MKL_INC_DIR) -I$(PLASMA_INC_DIR) -fopenmp -g
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
	BLAS_LIBS = -L$(BLAS_LIB_DIR) -lopenblas 

	CXX = g++-12
	CC = gcc-12
	CXXFLAGS = -fopenmp -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR)
endif

RT_OBJS = SuperM.o Matrix.o Tile.o TMatrix.o CoreBlas.o Check_Accuracy.o TileQR.o Right_Looking_Task.o

# for Performance evaluation
CXXFLAGS += -O3

# for Debug
# CXXFLAGS += -DDEBUG -g

# for Trace
# CXXFLAGS += -DTRACE
# RT_OBJS += trace.o

all: RT

RT : $(RT_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(RT_OBJS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(BLAS_LIBS) -lgomp

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

trace.o: trace.c
	$(CXX) -O3 -c -o $@ $<

clean:
	rm -f *.o RT
