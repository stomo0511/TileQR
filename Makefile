#
BLAS_ROOT = /opt/OpenBLAS
BLAS_INC_DIR = $(BLAS_ROOT)/include
BLAS_LIB_DIR = $(BLAS_ROOT)/lib
BLAS_LIBS = -lopenblas
SBLAS_LIBS = -lopenblas_seq
#
LAPACK_ROOT = /opt/LAPACK_3.6.0
LAPACK_INC_DIR = $(LAPACK_ROOT)/include
LAPACK_LIB_DIR = $(LAPACK_ROOT)/lib
LAPACK_LIBS = -llapacke -llapack
#
PLASMA_ROOT = /opt/PLASMA
PLASMA_INC_DIR = $(PLASMA_ROOT)/include
PLASMA_LIB_DIR = $(PLASMA_ROOT)/lib
PLASMA_LIBS = -lplasma -lcoreblas -lquark 
#
TMATRIX_ROOT = /Users/stomo/WorkSpace/TileAlgorithm/TileMatrix
TMATRIX_INC_DIR = $(TMATRIX_ROOT)
TMATRIX_LIB_DIR = $(TMATRIX_ROOT)
TMATRIX_LIBS = -lTileMatrix
#
COREBLAS_ROOT = /Users/stomo/WorkSpace/TileAlgorithm/CoreBlas
COREBLAS_INC_DIR = $(COREBLAS_ROOT)
COREBLAS_LIB_DIR = $(COREBLAS_ROOT)
COREBLAS_LIBS = -lCoreBlasTile
#
CXX =	/usr/local/bin/g++ -fopenmp
# for DEBUG
#CXXFLAGS =	-DDEBUG -g -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR) -I$(TMATRIX_INC_DIR) -I$(COREBLAS_INC_DIR)
# for Performance evaluation
CXXFLAGS =	-O2 -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR) -I$(TMATRIX_INC_DIR) -I$(COREBLAS_INC_DIR)

LLOBJS =		TileQR.o Check_Accuracy.o Progress.o LeftLooking.o 
SPOBJS =		TileQR.o Check_Accuracy.o Progress.o StaticPipeline.o
DSOBJS =		TileQR.o Check_Accuracy.o Progress.o DynamicSched.o
RLOBJS =		TileQR.o Check_Accuracy.o RightLooking.o
RTOBJS =		TileQR.o Check_Accuracy.o RightLooking_Task.o
QROBJS =		geqrf.o Check_Accuracy.o

all:	RL RT LL SP DS

geqrf: $(QROBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(QROBJS) \
				-L$(BLAS_LIB_DIR) $(BLAS_LIBS) \
				-L$(LAPACK_LIB_DIR) $(LAPACK_LIBS)

LL:	$(LLOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(LLOBJS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(SBLAS_LIBS)

SP:	$(SPOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(SPOBJS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(SBLAS_LIBS)

DS:	$(DSOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(DSOBJS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(SBLAS_LIBS)

RL:	$(RLOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(RLOBJS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(SBLAS_LIBS)

RT:	$(RTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(RTOBJS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(SBLAS_LIBS)

clean:
	rm -f $(RTOBJS) $(RLOBJS) $(LLOBJS) $(SPOBJS) $(DSOBJS) $(QROBJS)
