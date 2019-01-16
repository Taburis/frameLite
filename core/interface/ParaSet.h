
#ifndef ParaSet_H
#define ParaSet_H
#include "any.h"
#include <iostream>
#include <unordered_map>
#include <vector>
//create or set a parameter named pt with initial value 10 : setPara("pt", 10)
//get a exist parameter named pt : getPara("pt")
//create a new parameter named pt with inital value 10 : addPara("pt", 10)
class ParaSetBase{
		public: 
				ParaSetBase(){};
				ParaSetBase(const char * name) : set_name(name){};
				virtual ~ParaSetBase(){
						table.clear();
				}
				bool exists(std::string name){
						return  table.find(name) == table.end() ? 0 : 1;
				}
				template <typename T>
						bool addPara(const char* pname, const T & val){
								if(exists(pname)){
										std::cout<<"parameter: "<<pname<<" exists in the Set"<<std::endl;
										return 1;
								}
								table[pname] = any(T(val));
								return 0;
						}
				template <typename T>
						bool setPara(const char* pname, const T & val){
								table[pname] = any(val);
								return 0;
						}
				template <typename T>
						T getPara(const char* pname){
								if(!exists(pname)) {
										std::cout<<"parameter '"<<pname<<"' haven't been  defined in the ParaSet:"<<set_name<<std::endl;
										return T(0);
								}
								return any_cast<T>(table[pname]);
						}

				template <typename T>
						T securePara(const char* pname, T def){
								if(!exists(pname)) {
										std::cout<<"parameter '"<<pname<<"' not defined! using defalut value."<<std::endl;
										setPara<T>(pname, def);
										return T(def);
								}
								return any_cast<T>(table[pname]);
						}
				ParaSetBase* clone(){
						auto cps = new ParaSetBase();
						cps->table = table;
						return cps;
				}
				ParaSetBase& operator=(const ParaSetBase & rhs){
						set_name = rhs.set_name;
						table = rhs.table;
						return (*this);
				};
		public: 
				std::string set_name="";
				std::unordered_map<std::string, any> table;
};


class ParaSet : public ParaSetBase{
		public: ParaSet(): ParaSetBase(){};
				ParaSet(const char * name) : ParaSetBase(name){};
				virtual ~ParaSet(){
						for(auto it = vtable.begin(); it!= vtable.end(); ++it){
								//(*it).second->clear();
								//swap with an empty vector will free the memory of the vector
								std::vector<any>().swap(*((*it).second));
								delete (*it).second;
						}
						vtable.clear();
				}
				bool exists(std::string name){
						if(!(table.find(name) == table.end())) return 1;
						else if(!(vtable.find(name) == vtable.end())) return 1;
						else return 0;
				}
				ParaSet& operator=(const ParaSet& rhs){
						set_name = rhs.set_name;
						table = rhs.table;
						vtable = rhs.vtable;
						return (*this);
				};
				template <typename T>
						T getVPara(const char * pname, int i){
								if(!existsV(pname)) {
										std::cout<<"parameter vector '"<<pname<<"' haven't been  defined in the ParaSet:"<<set_name<<std::endl;
										return NAN;
								}
								if( vtable[pname]->size()<=size_t(i)) {
										std::cout<<"exceeds the length of parameter vector '"<<pname<<"'->("<<vtable[pname]->size()<<") in "<<set_name<<std::endl;
										return NAN;
								}
								return any_cast<T>(vtable[pname]->at(i));
						}
				std::unordered_map<std::string, std::vector<any>*> clone_vtable(){
						std::unordered_map<std::string, std::vector<any>*> vtable1;
						for( auto &it : vtable){
								auto vec = new std::vector<any>(*(it.second));
								vtable1[it.first] = vec;
						}
						return vtable1;
				}
				ParaSet* clone(){
						auto cps = new ParaSet();
						cps->table = table;
						if(vtable.size()!=0){
								cps->vtable = clone_vtable();
						}
						return cps;
				}
				size_t length(const char* pname){
						if(existsV(pname)) return vtable[pname]->size();
						else return 0;
				}
				std::vector<any> * pullVector(const char* pname){
								if(!existsV(pname)) {
										std::cout<<"the parameter vector '"<<pname<<"' doesnt exists in set:"<<set_name<<std::endl;
										return NULL;
								}
								return vtable[pname];
				}
				template <typename T>
						std::vector<T> * castParaVector(const char* pname){
								auto p = pullVector(pname);
								if( p==NULL) return NULL;
								std::vector<T> *vec = new std::vector<T>();
								vec->reserve(p->size());
								for(auto it = p->begin(); it<p->end(); ++it){
										//cout<<any_cast<T>(*it)<<endl;
										vec->push_back(any_cast<T>(*it));
								}
								return vec;
						}
				bool existsV(std::string name){
						return  vtable.find(name) == vtable.end() ? 0 : 1;
				}
				template <typename T>
						bool insertVPara(const char* pname, const T & val){
								if(!existsV(pname)){
										vtable[pname] = new std::vector<any>();
								}
								vtable[pname]->push_back(any(T(val)));
								return 0;
						}
				template <typename T>
						bool setParaVector(const char* pname, int n, const T* arr){
								if(existsV(pname)){
										std::cout<<"parameter vector '"<<pname<<"' exists in the ParaSet:"<<set_name<<std::endl;
										return 1;
								}
								vtable[pname] = new std::vector<any>();
								vtable[pname]->reserve(size_t(n));
								for(size_t jj = 0; jj<size_t(n); ++jj){
										vtable[pname]->push_back(any(arr[jj]));
								}
								return 0;
						}
				template <typename T>
						T* getVectorAsArray(const char * pname){
								auto vec = pullVector(pname);
								T* p_ = new T[vec->size()];
//								cout<<"vector size is "<<vec->size()<<endl;
								for(auto it = vec->begin(); it<vec->end(); ++it){
										p_[it-vec->begin()] = any_cast<T>(*(it));
								}	
								return p_;
						}
				bool checkConflict(ParaSet &b){
						// return 1 if the key exists in both this and a;
						for(auto & it: b.table )  if(exists(it.first)) return 1;
						for(auto & it: b.vtable ) if(exists(it.first)) return 1;
						return 0;
				}

				void clear(){
						table.clear();
						for(auto it = vtable.begin(); it!= vtable.end(); ++it){
								//(*it).second->clear();
								//swap with an empty vector will free the memory of the vector
								std::vector<any>().swap(*((*it).second));
						}
						vtable.clear();
				}

				ParaSet operator+(ParaSet & rhs){
						ParaSet sum = rhs;
						
						if(checkConflict(rhs)){
								sum.clear();
							   	return sum;
						}
						sum.table.insert(table.begin(), table.end());
						sum.vtable = rhs.clone_vtable();
						auto vtable2 = clone_vtable();
						sum.vtable.insert(vtable2.begin(), vtable2.end());
						return sum;
				}

		public: 
				std::unordered_map<std::string, std::vector<any>*> vtable;
};


#endif
