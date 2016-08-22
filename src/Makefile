CXX = g++
MPICXX = mpicxx
CFLAGS = --std=c++11 -O3 -fopenmp

INC = -I$(PWD) -I$(PWD)/model -I$(PWD)/solver -I$(PWD)/linearsolver -I$(PWD)/utils
INC += -I$(PWD)/thirdparty
LIB = -L$(PWD)/thirdparty -lcnpy -lz

INC_ADOLC = -I/home/anandps/local/adolc_new/include
LIB_ADOLC = -L/home/anandps/local/adolc_new/lib64 -ladolc

INC_ARMA = -I/home/anandps/local/armadillo/usr/include -I/home/anandps/local/superlu/usr/local/include
LIB_ARMA = -L/home/anandps/local/armadillo/usr/lib -larmadillo -L/home/anandps/local/superlu/usr/local/lib/ -lsuperlu -llapack -lblas -lgfortran

INC_PETSC = -I/home/anandps/local/petsc_3_7_3/include
LIB_PETSC = -L/home/anandps/local/petsc_3_7_3/lib -lpetsc -lgfortran

make: arma eigen petsc explicit

explicit:
	$(CXX) main.cpp -o ../bin/main_explicit $(CFLAGS) $(INC) $(LIB)
eigen:
	$(CXX) main.cpp -o ../bin/main_eigen $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC)  -DENABLE_EIGEN
arma:
	$(CXX) main.cpp -o ../bin/main_arma $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC) -DENABLE_ARMA $(INC_ARMA) $(LIB_ARMA)
petsc:
	$(MPICXX) main.cpp -o ../bin/main_petsc $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC) -DENABLE_PETSC $(INC_PETSC) $(LIB_PETSC)
clean:
	rm main
