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
#VT_INC = /opt/openmpi/include/vampirtrace
VT_INC = /usr/include/openmpi-x86_64
#
CXX = g++-7 -fopenmp -DDEBUG
# for Debug
#CXXFLAGS = -DDEBUG -g -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR)

# for Performace Ecaluation
CXXFLAGS = -O3 -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR)

MR_OBJS = SuperM.o Matrix.o Tile.o TMatrix.o CoreBlas.o Progress.o Check_Accuracy.o TileQR.o MRight_Looking.o

TS_OBJS = SuperM.o Matrix.o Tile.o TMatrix.o CoreBlas.o Check_Accuracy.o TSQR.o 

DS_OBJS = SuperM.o Matrix.o Tile.o TMatrix.o CoreBlas.o Progress.o Check_Accuracy.o TileQR.o Dynamic_Sched.o

LL_OBJS = SuperM.o Matrix.o Tile.o TMatrix.o CoreBlas.o Progress.o Check_Accuracy.o TileQR.o Left_Looking.o

SP_OBJS = SuperM.o Matrix.o Tile.o TMatrix.o CoreBlas.o Progress.o Check_Accuracy.o TileQR.o Static_Pipeline.o

ML_OBJS = SuperM.o Matrix.o Tile.o TMatrix.o CoreBlas.o Progress.o Check_Accuracy.o TileQR.o MLeft_Looking.o

RL_OBJS = SuperM.o Matrix.o Tile.o TMatrix.o CoreBlas.o Check_Accuracy.o TileQR.o Right_Looking.o

RT_OBJS = SuperM.o Matrix.o Tile.o TMatrix.o CoreBlas.o Check_Accuracy.o TileQR.o Right_Looking_Task.o

all: MR TS RL ML DS LL SP RT

MR : $(MR_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(MR_OBJS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(SBLAS_LIBS) $(OTHER_LIBS)


SS : $(SS_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(SS_OBJS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(SBLAS_LIBS) $(OTHER_LIBS)

TS : $(TS_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(TS_OBJS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(SBLAS_LIBS) $(OTHER_LIBS)

DS : $(DS_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(DS_OBJS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(SBLAS_LIBS) $(OTHER_LIBS)

LL : $(LL_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(LL_OBJS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(SBLAS_LIBS) $(OTHER_LIBS)

SP : $(SP_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(SP_OBJS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(SBLAS_LIBS) $(OTHER_LIBS)

ML : $(ML_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(ML_OBJS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(SBLAS_LIBS) $(OTHER_LIBS)

RL : $(RL_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(RL_OBJS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(SBLAS_LIBS) $(OTHER_LIBS)

RT : $(RT_OBJS)
	$(CXX) $(CFLAGS) -o $@ $(RT_OBJS) -L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) -L$(BLAS_LIB_DIR) $(SBLAS_LIBS) $(OTHER_LIBS)

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<
clean:
	rm -f Makefile~ *.cpp~ *.h~ *.o
