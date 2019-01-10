
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
		public :
				matrixTObjPtr(): matrixPtrHolder<T>(){};	
				matrixTObjPtr(const char * _name, int n, int m): matrixPtrHolder<T>(n,m) {
						name = _name;
				};
				virtual ~matrixTObjPtr(){
						if(doFree) for(auto & it : matrixPtrHolder<T>::ref) delete it;
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
				//				 void write(){ for(auto & it:ref) it->Write();}

				std::string name;
				bool doFree = 0;
};

class matrixTH1Ptr : public matrixTObjPtr<TH1>{
		public : matrixTH1Ptr(): matrixTObjPtr<TH1>(){};
				 matrixTH1Ptr(const char * _name, int n, int m): matrixTObjPtr<TH1>(_name, n,m) {}
				 void setName(const char* _name) {
						 name = _name;
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 auto hn = name+"_"+i+"_"+j;
										 at(i,j)->SetName(hn);
								 }
						 }
				 }
				 matrixTH1Ptr * clone(const char * newname){
						 auto m2 = new matrixTH1Ptr(newname, matrixTObjPtr<TH1>::nrow, matrixTObjPtr<TH1>::ncol);
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 auto hn = string(newname)+"_"+i+"_"+j;
										 auto hh = (TH1*) at(i,j)->Clone(hn);
										 m2->add(hh, i, j);
								 }
						 }
						 return m2;
				 }
				 matrixTH1Ptr * operator+( matrixTH1Ptr & rhs){
						 //check if the matrix shapes are the same
						 if( matrixTObjPtr<TH1>::ncol != rhs.ncol || matrixTObjPtr<TH1>::nrow != rhs.nrow) return 0;
						 auto m2 = clone((string("sum_")+name).c_str());
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 m2->at(i,j)->Add(rhs(i,j));
								 }
						 }
						 return m2;
				 }
				 matrixTH1Ptr * operator/( matrixTH1Ptr & rhs){
						 if( matrixTObjPtr<TH1>::ncol != rhs.ncol || matrixTObjPtr<TH1>::nrow != rhs.nrow) return 0;
						 auto m2 = clone((string("division_")+name).c_str());
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 m2->at(i,j)->Divide(rhs(i,j));
								 }
						 }
						 return m2;
				 }
				 matrixTH1Ptr * operator%( matrixTH1Ptr & rhs){
						 if( matrixTObjPtr<TH1>::ncol != rhs.ncol || matrixTObjPtr<TH1>::nrow != rhs.nrow) return 0;
						 auto m2 = clone((string("bdivided_")+name).c_str());
						 for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
								 for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
										 m2->at(i,j)->Divide(m2->at(i,j), rhs(i,j), 1, 1, "B");
								 }
						 }
						 return m2;
				 }

};

template<typename T> 
class multi_canvas : public TCanvas{
		public: multi_canvas(matrixTObjPtr<T> & p): TCanvas(TString(p.name+"_canvas"), TString(p.name), int(p.ncol*350), int(p.nrow*325))
				{
						nrow= p.nrow;
						ncol= p.ncol;
						m2ptr = &p;
						Divide(ncol, nrow);
				}
				multi_canvas(TString name, TString title, int mrow, int mcol): TCanvas(name, title, mcol*325, mrow*325)
		{
				nrow= mrow;
				ncol= mcol;
				Divide(ncol, nrow);
		}
				int flatten(int n, int m){
						if(n > nrow -1 || m > ncol -1 ) {
								std::cout<<" ==== ERROR: index ("<<n<<", "<<m<<") exceeds the range!!! current shape ["<<nrow<<", "<<ncol<<"] ===="<<std::endl;
								return 0;
						}
						//slightly different from the matrixPtrHolder, the row and column in canvas start from 1 and up, 0 stands for the whole pad.
						return n*ncol+m+1;}
				void CD(int n, int m){ this->cd(flatten(n,m));}
				void draw(){
						this->SetMargin(0.5, 0.5, 0.4, 0.40);
						//gStyle->SetOptStat(0);
						for(int i=0; i< nrow; ++i){
								for(int j=0; j< ncol; ++j){
										cd(i,j);
										m2ptr(i,j)->Draw();
								}
						}
				}

		public:
				int nrow, ncol;
				matrixTObjPtr<T> *m2ptr;
};


#endif
