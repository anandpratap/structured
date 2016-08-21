CXX = g++
all:
	$(CXX) -fPIC -c cnpy.cpp -o cnpy.o
	$(CXX) -shared  -Wl,-soname,libcnpy.so -o libcnpy.so cnpy.o
