#
BLAS_INC_DIR = $(MKLROOT)/include
BLAS_LIB_DIR = $(MKLROOT)/lib/intel64
SBLAS_LIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
#
PLASMA_ROOT = /opt/plasma
PLASMA_INC_DIR = $(PLASMA_ROOT)/include
PLASMA_LIB_DIR = $(PLASMA_ROOT)/lib
PLASMA_LIBS = -lplasma -lcoreblas
#
CXX = g++-7 -fopenmp

# for Debug
#CXXFLAGS = -DDEBUG -g -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR)

# for Performace Ecaluation
CXXFLAGS = -O3 -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR)

RT_OBJS = SuperM.o Matrix.o Tile.o TMatrix.o CoreBlas.o Check_Accuracy.o TileQR.o Right_Looking_Task.o

all: RT

RT : $(RT_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(RT_OBJS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(SBLAS_LIBS) $(OTHER_LIBS)

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<
clean:
	rm -f Makefile~ *.cpp~ *.h~ *.o
