
#ifndef MATRIXTOBJPTR_H
#define MATRIXTOBJPTR_H
#include <string>
#include <iostream>

template<typename T>
class matrixPtrHolder {
		public : matrixPtrHolder(int n, int m): nrow(n), ncol(m){ref.resize(n*m);}
				 virtual ~matrixPtrHolder(){}
				 int flatten(int n, int m){
						if(n > nrow-1 || m > ncol -1 ) {
								std::cout<<" ==== ERROR: index exceeds the range!!! current shape ["<<nrow<<", "<<ncol<<"] ===="<<std::endl;
								return 0;
						}
						 // the start index of row/column is 0 rather than 1
						 return n+m*nrow;}
				 T* at(int n, int m){return ref[flatten(n, m)];}
				 T* operator()(int n, int m){ return at(n, m);}
				 void add(T* ptr, int i, int j){
						// ref.push_back(ptr);
						 ref[flatten(i,j)] = ptr;
				 }
				 void cleanAll(){for(auto & it: ref) delete it;}
				 void transpose(){
						 std::vector<T*> tmp(ref);
						 for(int i=0; i<nrow; i++){
								 for(int j=0; j<ncol; j++){
										 //cout<<"i "<<i<<", j "<<j<<endl;
										 ref[j+i*ncol] = tmp[flatten(i,j)];
										 //cout<<tmp[flatten(i,j)]->GetName()<<endl;
								 }
						 }
						 tmp.clear();
						 int trans = nrow;
						 nrow = ncol;
						 ncol = trans;
				 }
		public :
				 std::vector<T*> ref;
				 int nrow, ncol;
};

template<typename T> 
class matrixTObjPtr : public matrixPtrHolder<T>{
		public : matrixTObjPtr(const char * _name, int n, int m): matrixPtrHolder<T>(n,m) {
						 name = _name;
				 };
				 virtual ~matrixTObjPtr(){};
				 void autoLoad(TFile* f){
						 for(int j=0; j<matrixPtrHolder<T>::ncol; ++j){
								 for(int i=0; i<matrixPtrHolder<T>::nrow; i++){
										 auto hn = name+"_"+i+"_"+j;
										 auto ptr = (T*)f->Get(hn);
										 TString hname = ptr->GetName();
										 std::cout<<"adding "+hname<<std::endl;
										 this->add(ptr, i,j);
								 }
						 }
				 };
				 //				 void write(){ for(auto & it:ref) it->Write();}

				 std::string name;
};

#endif
