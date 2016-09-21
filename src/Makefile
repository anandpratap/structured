CXX = g++
1;95;0c1;95;0cMPICXX = mpicxx
CFLAGS = -std=c++14 -O3 -fopenmp  -fstrict-aliasing -g -DRA_BOUNDSCHECK -g #-ftree-vectorize -ftree-vectorizer-verbose=1 -fopt-info-vec-all

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

SRC = main.cpp linearsolver/ls_eigen.cpp utils/io.cpp utils/config.cpp model/bc.cpp model/flux.cpp model/fluid.cpp model/eulerequation.cpp solver/solution.cpp solver/solver.cpp model/reconstruction.cpp utils/mesh.cpp
OBJ = $(subst .cpp,.o,$(SRC))
test:$(OBJ)
	$(CXX) -o test $(OBJ) $(LIB) $(LIB_ADOLC) $(CFLAGS)
	mv mujo $(BIN).

%.o:%.cpp
	$(CXX) $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC) -DENABLE_EIGEN  -o $@ -c $<
explicit:
	$(CXX) main.cpp -o ../bin/main_explicit $(CFLAGS) $(INC) $(LIB)

check:
	cppcheck -j 4 --enable=all main.cpp -I$(PWD) -I$(PWD)/model -I$(PWD)/solver -I$(PWD)/linearsolver -I$(PWD)/utils

eigen:
	$(CXX) main.cpp -o ../bin/main_eigen $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC)  -DENABLE_EIGEN

cpp:
	$(CXX) -c linearsolver/ls_eigen.cpp -o ls_eigen.o $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC)  -DENABLE_EIGEN
	$(CXX) -c utils/config.cpp -o config.o $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC)  -DENABLE_EIGEN
	$(CXX) -c model/bc.cpp -o bc.o $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC)  -DENABLE_EIGEN
	$(CXX) -c model/flux.cpp -o flux.o $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC)  -DENABLE_EIGEN
	$(CXX) -c model/fluid.cpp -o fluid.o $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC)  -DENABLE_EIGEN
	$(CXX) -c model/eulerequation.cpp -o eulerequation.o $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC)  -DENABLE_EIGEN
	$(CXX) -c solver/solution.cpp -o solution.o $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC)  -DENABLE_EIGEN
	$(CXX) -c solver/solver.cpp -o solver.o $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC)  -DENABLE_EIGEN
	$(CXX) -c utils/io.cpp -o io.o $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC)  -DENABLE_EIGEN
arma:
	$(CXX) main.cpp -o ../bin/main_arma $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC) -DENABLE_ARMA $(INC_ARMA) $(LIB_ARMA)
petsc:
	$(MPICXX) main.cpp -o ../bin/main_petsc $(CFLAGS) $(INC) $(LIB) $(INC_ADOLC) $(LIB_ADOLC) -DENABLE_PETSC $(INC_PETSC) $(LIB_PETSC)
clean:
	rm main $(OBJ)
