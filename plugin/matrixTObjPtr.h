
#ifndef MATRIXTOBJPTR_H
#define MATRIXTOBJPTR_H
#include <string>
#include <iostream>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"

template<typename T>
class matrixPtrHolder{
		public : matrixPtrHolder(){};
				 matrixPtrHolder(int n, int m){ setup(n,m);}
				 virtual ~matrixPtrHolder(){};
				 void setup(int n, int m){nrow = n, ncol=m, ref.resize(n*m);};
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
				 int size(){return nrow*ncol;};
				 bool isValid(int i, int j){// return 1 means i, j is valid
						 return flatten(i,j) < size(); }
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
				 int Ncol(){ return ncol; }
				 int Nrow(){ return nrow; }
		public :
				 std::vector<T*> ref;
				 int nrow, ncol;
};

template<typename T> 
class matrixTObjPtr : public matrixPtrHolder<T>{
		public :
				matrixTObjPtr(): matrixPtrHolder<T>(){};	
				matrixTObjPtr(const char * _name, int n, int m): matrixPtrHolder<T>(n,m) {
						name = _name;
				};
				virtual ~matrixTObjPtr(){
					//	if(doFree){
					//			std::cout<<"deleting "<<matrixPtrHolder<T>::ref.size()<<" objects in "<<name<<std::endl;
					//		   	for(auto & it : matrixPtrHolder<T>::ref){
					//					std::cout<<"deleting "<<it->GetName()<<std::endl;
					//				   	delete it;
					//			}
					//	}
						doFree = 0;
				};
				matrixTObjPtr * deep_clone(const char * name_){
						auto m2 = new matrixTObjPtr<T>();
						m2->setup(name_, matrixPtrHolder<T>::nrow, matrixPtrHolder<T>::ncol);
						return m2;
				}
				void setup(const char* _name, int n, int m){name = _name, matrixPtrHolder<T>::setup(n,m);};
				void autoLoad(TFile* f){
						for(int j=0; j<matrixPtrHolder<T>::ncol; ++j){
								for(int i=0; i<matrixPtrHolder<T>::nrow; i++){
										auto hn = name+"_"+i+"_"+j;
										//std::cout<<"adding "+hn<<std::endl;
										auto ptr = (T*) f->Get(hn);
										TString hname = ptr->GetName();
										std::cout<<"adding "+hname<<std::endl;
										this->add(ptr, i,j);
								}
						}
				};
				matrixTObjPtr<T>* load_m2TObj(const char * name, int n, int m, TFile* file){
						auto newm2 = matrixTObjPtr<T>(name, n, m);
						newm2->autoLoad(file);
						return newm2;
				}
				void write(){ for(auto & it:matrixPtrHolder<T>::ref) it->Write();}
				const char * getName(){return name.c_str();}

				std::string name;
				bool doFree = 0;
};
/*
class matrixTH1Ptr : public matrixTObjPtr<TH1>{
		public : matrixTH1Ptr(): matrixTObjPtr<TH1>(){};
				 virtual ~matrixTH1Ptr(){ };
				 matrixTH1Ptr(const char * _name, int n, int m): matrixTObjPtr<TH1>(_name, n,m) {}
				 matrixTH1Ptr(matrixTH1Ptr& m2): matrixTObjPtr<TH1>(TString("cl_"+m2.name), m2.Nrow(), m2.Ncol()) {}
				 void setName(const char* _name) {
						 name = _name;
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 auto hn = name+"_"+i+"_"+j;
										 at(i,j)->SetName(hn);
								 }
						 }
				 }
				 void scale(float s) {
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 at(i,j)->Scale(s);
								 }
						 }
				 }
				 matrixTH1Ptr * clone(const char * newname){
						 auto m2 = new matrixTH1Ptr(newname, matrixTObjPtr<TH1>::nrow, matrixTObjPtr<TH1>::ncol);
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 auto hn = std::string(newname)+"_"+i+"_"+j;
										 auto hh = (TH1*) at(i,j)->Clone(hn);
										 m2->add(hh, i, j);
								 }
						 }
						 m2->doFree = 1;
						 return m2;
				 }
				 matrixTH1Ptr * operator+( matrixTH1Ptr & rhs){
						 //check if the matrix shapes are the same
						 if( matrixTObjPtr<TH1>::ncol != rhs.ncol || matrixTObjPtr<TH1>::nrow != rhs.nrow) return 0;
						 auto m2 = clone((std::string("sum_")+name).c_str());
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 m2->at(i,j)->Add(rhs(i,j));
								 }
						 }
						 return m2;
				 }
				 matrixTH1Ptr * operator-( matrixTH1Ptr & rhs){
						 //check if the matrix shapes are the same
						 if( matrixTObjPtr<TH1>::ncol != rhs.ncol || matrixTObjPtr<TH1>::nrow != rhs.nrow) return 0;
						 auto m2 = clone((std::string("sum_")+name).c_str());
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 m2->at(i,j)->Add(m2->at(i,j), rhs(i,j), 1, -1);
								 }
						 }
						 return m2;
				 }
				 matrixTH1Ptr * operator/( matrixTH1Ptr & rhs){
						 if( matrixTObjPtr<TH1>::ncol != rhs.ncol || matrixTObjPtr<TH1>::nrow != rhs.nrow) return 0;
						 auto m2 = clone((std::string("division_")+name).c_str());
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 m2->at(i,j)->Divide(rhs(i,j));
								 }
						 }
						 return m2;
				 }

				 matrixTH1Ptr * operator%( matrixTH1Ptr & rhs){
						 if( matrixTObjPtr<TH1>::ncol != rhs.ncol || matrixTObjPtr<TH1>::nrow != rhs.nrow) return 0;
						 auto m2 = clone((std::string("bdivided_")+name).c_str());
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 m2->at(i,j)->Divide(m2->at(i,j), rhs(i,j), 1, 1, "B");
								 }
						 }
						 return m2;
				 }
				 void normalized_by_area(float x = 1, float y = -1){
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 auto h = at(i,j);
										 float area;
										 if(x< y){
												 int n1 = h->GetXaxis()->FindBin(x);
												 int n2 = h->GetXaxis()->FindBin(y);
												 area = h->Integral(n1, n2, "width");
										 }else area = h->Integral("width");
										 h->Scale(1.0/area);
								 }
						 }
						 return;
				 }
};*/



#endif
