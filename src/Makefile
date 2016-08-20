CXX = g++


INC = -I/home/anandps/local/adolc/include -I/home/anandps/local/armadillo/usr/include -I$(PWD) -I/home/anandps/local/superlu/usr/local/include
INC += -I/home/anandps/local/petsc_3_4_5_without_mpi/include

LIB = -L/home/anandps/local/adolc/lib64 -ladolc -lgfortran
LIB += -L/home/anandps/local/armadillo/usr/lib -larmadillo
LIB += -L/home/anandps/local/petsc_3_4_5_without_mpi/lib
LIB += -L/home/anandps/local/superlu/usr/local/lib/ -lsuperlu -llapack -lblas
CFLAGS = --std=c++11 -O3

make: arma eigen
eigen:
	$(CXX) main.cpp -o main_eigen $(CFLAGS) $(INC) $(LIB) -DENABLE_EIGEN
arma:
	$(CXX) main.cpp -o main_arma $(CFLAGS) $(INC) $(LIB)	-DENABLE_ARMA

clean:
	rm main
