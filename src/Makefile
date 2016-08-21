CXX = g++
MPICXX = mpicxx


INC = -I/home/anandps/local/adolc/include -I/home/anandps/local/armadillo/usr/include -I$(PWD) -I/home/anandps/local/superlu/usr/local/include
INC += -I/home/anandps/local/petsc_3_7_3/include
INC += -I/home/anandps/Dropbox/research/fall_2016/structured/src/thirdparty

LIB = -L/home/anandps/local/adolc/lib64 -ladolc -lgfortran
LIB += -L/home/anandps/local/armadillo/usr/lib -larmadillo
LIB += -L/home/anandps/local/petsc_3_7_3/lib -lpetsc
LIB += -L/home/anandps/local/superlu/usr/local/lib/ -lsuperlu -llapack -lblas
LIB += -L$(PWD)/thirdparty -lcnpy -lz
CFLAGS = --std=c++11 -O3 -fopenmp

make: arma eigen petsc
eigen:
	$(CXX) main.cpp -o main_eigen $(CFLAGS) $(INC) $(LIB) -DENABLE_EIGEN
arma:
	$(CXX) main.cpp -o main_arma $(CFLAGS) $(INC) $(LIB) -DENABLE_ARMA
petsc:
	$(MPICXX) main.cpp -o main_petsc $(CFLAGS) $(INC) $(LIB) -DENABLE_PETSC

clean:
	rm main
