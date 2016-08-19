CXX = g++

INC = -I/home/anandps/local/adolc/include -I/home/anandps/local/armadillo/usr/include -I$(PWD) -I/home/anandps/local/superlu/usr/local/include

LIB = -L/home/anandps/local/adolc/lib64 -ladolc -lgfortran  -L/home/anandps/local/superlu/usr/local/lib -lsuperlu -llapack -lblas
#LIB += -L/home/anandps/local/armadillo/usr/lib -larmadillo

CFLAGS = --std=c++11 -O3

all:
	$(CXX) main.cpp -o main $(CFLAGS) $(INC) $(LIB)

clean:
	rm main
