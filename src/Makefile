CXX = g++


INC = -I/home/anandps/local/adolc/include -I/home/anandps/local/armadillo/usr/include -I$(PWD) -I/home/anandps/local/superlu/usr/local/include
INC += -I/home/anandps/local/petsc_3_4_5_without_mpi/include

LIB = -L/home/anandps/local/adolc/lib64 -ladolc -lgfortran  -L/home/anandps/local/superlu/usr/local/lib -lsuperlu -llapack -lblas
LIB += -L/home/anandps/local/armadillo/usr/lib -larmadillo
LIB += -I/home/anandps/local/petsc_3_4_5_without_mpi/lib

CFLAGS = --std=c++11 -O3

all:
	$(CXX) main.cpp -o main_eigen $(CFLAGS) $(INC) $(LIB) -DENABLE_EIGEN
	$(CXX) main.cpp -o main_arma $(CFLAGS) $(INC) $(LIB)	-DENABLE_ARMA

clean:
	rm main
