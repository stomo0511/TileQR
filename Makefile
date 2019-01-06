#
UNAME = $(shell uname)
ifeq ($(UNAME),Linux)
  BLAS_ROOT = /opt/intel/compilers_and_libraries/linux/mkl
  BLAS_INC_DIR = $(BLAS_ROOT)/include
  BLAS_LIB_DIR = $(BLAS_ROOT)/lib/intel64
  SBLAS_LIBS = -Wl,--start-group $(BLAS_LIB_DIR)/libmkl_intel_lp64.a $(BLAS_LIB_DIR)/libmkl_sequential.a $(BLAS_LIB_DIR)/libmkl_core.a -Wl,--end-group -lpthread -ldl -lm
  TMATRIX_ROOT = /home/stomo/WorkSpace/TileMatrix
  COREBLAS_ROOT = /home/stomo/WorkSpace/CoreBlas
  CXX =	eztrace_cc g++
endif
ifeq ($(UNAME),Darwin)
  BLAS_ROOT = /opt/intel/compilers_and_libraries/mac/mkl
  BLAS_INC_DIR = $(BLAS_ROOT)/include
  BLAS_LIB_DIR = $(BLAS_ROOT)/lib
  SBLAS_LIBS = $(BLAS_LIB_DIR)/libmkl_intel_lp64.a $(BLAS_LIB_DIR)/libmkl_sequential.a $(BLAS_LIB_DIR)/libmkl_core.a -lpthread -ldl -lm
  TMATRIX_ROOT = /Users/stomo/WorkSpace/TileAlgorithm/TileMatrix
  COREBLAS_ROOT = /Users/stomo/WorkSpace/TileAlgorithm/CoreBlas
  CXX =	/usr/local/bin/g++ 
endif

#
PLASMA_ROOT = /opt/plasma-17.1
PLASMA_INC_DIR = $(PLASMA_ROOT)/include
PLASMA_LIB_DIR = $(PLASMA_ROOT)/lib
PLASMA_LIBS = -lcoreblas -lplasma
#
TMATRIX_INC_DIR = $(TMATRIX_ROOT)
TMATRIX_LIB_DIR = $(TMATRIX_ROOT)
TMATRIX_LIBS = -lTileMatrix
#
COREBLAS_INC_DIR = $(COREBLAS_ROOT)
COREBLAS_LIB_DIR = $(COREBLAS_ROOT)
COREBLAS_LIBS = -lCoreBlasTile
#
CXXFLAGS =	-fopenmp -m64 -O2 -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR) -I$(TMATRIX_INC_DIR) -I$(COREBLAS_INC_DIR)
#
LLOBJS =		TileQR.o Check_Accuracy.o Progress.o LeftLooking.o 
LTOBJS =		TileQR.o Check_Accuracy.o LeftLooking_Task.o 
SPOBJS =		TileQR.o Check_Accuracy.o Progress.o StaticPipeline.o
DSOBJS =		TileQR.o Check_Accuracy.o Progress.o DynamicSched.o
RLOBJS =		TileQR.o Check_Accuracy.o RightLooking.o
RTOBJS =		TileQR.o Check_Accuracy.o RightLooking_Task.o
QROBJS =		geqrf.o Check_Accuracy.o

all: RL RT DS SP LL LT

RT:	$(RTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(RTOBJS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(SBLAS_LIBS)

LT:	$(LTOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(LTOBJS) \
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

SP:	$(SPOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(SPOBJS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(SBLAS_LIBS)

LL:	$(LLOBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(LLOBJS) \
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

geqrf: $(QROBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(QROBJS) \
				-L$(BLAS_LIB_DIR) $(BLAS_LIBS) 

clean:
	rm -f $(RTOBJS) $(LTOBJS) $(RLOBJS) $(LLOBJS) $(SPOBJS) $(DSOBJS) RL RT DS SP LL LT
