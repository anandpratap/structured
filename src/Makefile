CXX = g++
MPICXX = mpicxx


INC = -I/home/anandps/local/adolc_new/include  -I$(PWD) -I$(PWD)/model -I$(PWD)/solver -I$(PWD)/linearsolver -I$(PWD)/utils
LIB = -L/home/anandps/local/adolc_new/lib64 -ladolc


INC += -I$(PWD)/thirdparty
LIB += -L$(PWD)/thirdparty -lcnpy -lz
CFLAGS = --std=c++11 -O3 -fopenmp


INC_ARMA = -I/home/anandps/local/armadillo/usr/include -I/home/anandps/local/superlu/usr/local/include
LIB_ARMA = -L/home/anandps/local/armadillo/usr/lib -larmadillo -L/home/anandps/local/superlu/usr/local/lib/ -lsuperlu -llapack -lblas -lgfortran

INC_PETSC = -I/home/anandps/local/petsc_3_7_3/include
LIB_PETSC = -L/home/anandps/local/petsc_3_7_3/lib -lpetsc -lgfortran

make: arma eigen petsc

eigen:
	$(CXX) main.cpp -o ../bin/main_eigen $(CFLAGS) $(INC) $(LIB) -DENABLE_EIGEN
arma:
	$(CXX) main.cpp -o ../bin/main_arma $(CFLAGS) $(INC) $(LIB) -DENABLE_ARMA $(INC_ARMA) $(LIB_ARMA)
petsc:
	$(MPICXX) main.cpp -o ../bin/main_petsc $(CFLAGS) $(INC) $(LIB) -DENABLE_PETSC $(INC_PETSC) $(LIB_PETSC)
clean:
	rm main
