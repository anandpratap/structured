#ifndef __ARRAY_H
#define __ARRAY_H

template<class T>
class Array3D{
public:
	T *data;
	uint ni = 0U;
	uint nj = 0U;
	uint nk = 0U;
	Array3D(){data = nullptr;};
	Array3D(uint val_ni, uint val_nj, uint val_nk){
		ni = val_ni; nj = val_nj; nk = val_nk;
		data  = new T[ni*nj*nk]();
	};
	
	void increment(T value, uint i, uint j, uint k){
		data[i*nj*nk + j*nk + k] += value;
	}
	
	void set(T value, uint i, uint j, uint k){
		data[i*nj*nk + j*nk + k] = value;
	}
	
	inline T &operator()(uint i, uint j, uint k){
		//		spdlog::get("console")->info("{} {} {}", i, j, k);
		//		assert(i>=0 && i<ni);
		//assert(j>=0 && j<nj);
		//assert(k>=0 && k<nk);
		//if(data != nullptr){
			return data[i*nj*nk + j*nk + k];
			//	}
	};
	//const T operator()(uint i, uint j, uint k) const {return data[i*nj*nk + j*nk + k];};
	
	inline T& operator()(uint i){return data[i];};
	//const T operator()(uint i) const {return data[i];};
	
	uint get_size(){return ni*nj*nk;};

	void print(){
		for(uint i=0U; i<get_size(); i++){
			std::cout<<i<<" "<<data[i]<<std::endl;
		}
	}
	
	~Array3D(){
		delete[] data;
	};
};


template<class T>
class Array2D{
public:
	T *data;
	uint ni = 0U;
	uint nj = 0U;
	Array2D(){data = nullptr;};
	Array2D(uint val_ni, uint val_nj){
		ni = val_ni; nj = val_nj;
		data  = new T[ni*nj]();
	};
	inline T &operator()(uint i, uint j){
		//	assert(i>=0 && i<ni);
		//assert(j>=0 && j<nj);
		//if( data != nullptr){
			return data[i*nj + j];
			//}
	};
	//const T operator()(uint i, uint j) const {return data[i*nj + j];};
	void increment(T value, uint i, uint j){
		data[i*nj + j] += value;
	}
	
	void set(T value, uint i, uint j){
		data[i*nj + j] = value;
	}

	//	T &operator()(uint i){return data[i];};
	//const T operator()(uint i) const {return data[i];};
	inline T& operator()(uint i){return data[i];};
	uint get_size(){return ni*nj;};
	void print(){
		for(uint i=0U; i<get_size(); i++){
			std::cout<<i<<" "<<data[i]<<std::endl;
		}
	}
	~Array2D(){
		delete[] data;
	};
};


template<class T>
class Array1D{
public:
	T *data;
	uint ni = 0U;
	Array1D(){data=nullptr;};
	Array1D(uint val_ni){
		ni = val_ni;
		data  = new T[ni]();
	};
	inline T &operator()(uint i){
		//assert(i>=0 && i<ni);
		//if( data != nullptr){
			return data[i];
			//}
	};
	//const T operator()(uint i) const {return data[i];};

	void increment(T value, uint i){
		data[i] += value;
	}
	
	void set(T value, uint i){
		data[i] = value;
	}
	uint get_size(){return ni;};
	void print(){
		for(uint i=0U; i<get_size(); i++){
			std::cout<<i<<" "<<data[i]<<std::endl;
		}
	}
	~Array1D(){
		delete[] data;
	};
};

#endif
