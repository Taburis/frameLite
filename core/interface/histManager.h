
#ifndef EMD_HISTMANAGER_H
#define EMD_HISTMANAGER_H
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <unordered_map>

class histManager {
		public: 
				histManager(){};
				~histManager(){
//						for(auto & it: hkey){
//								delete it.second;
//						}
						for(auto & it: th1){
								delete it;
						}
						hkey.clear();
				};
				template <typename T> T* regHist (const char* name, const char* title, int nbin, double x, double y);
				template <typename T> T* regHist (const char* name, const char* title, int nbin, const float *bins);
				template <typename T> T* regHist (const char* name, const char* title, int nbin, const Double_t *bins);
				template <typename T> T* regHist (const char* name, const char* title, int nbin, const Double_t *bins, Int_t nbin2, Double_t ylow, Double_t yup);
				template <typename T> T* regHist (const char* name, const char* title, int nx, double x1, double x2, 
								int ny, double y1, double y2);
				template <typename T> T* get(const char* name);
				template <typename T> int fill(const char * name, double x, double w=1){
						return get<T>(name)->Fill(x, w);
				}
				void sumw2();
				void write();
				TH1* operator[](const char* name){
						return hkey[name];
				}
		public:
				std::vector<TH1*> th1;
				std::unordered_map<const char*, TH1*> hkey;
};

template <typename T> T* histManager::regHist (const char* name, const char* title, int nbin, double x, double y){
		T *h = new T(name, title, nbin, x, y);
		th1.push_back(dynamic_cast<TH1*>(h));
		hkey[name] = th1.back();
		return h;
}

template <typename T> T* histManager::regHist (const char* name, const char* title, int nbin, const float* bins){
		T *h = new T(name, title, nbin, bins);
		th1.push_back(dynamic_cast<TH1*>(h));
		hkey[name] = th1.back();
		return h;
}
template <typename T> T* histManager::regHist (const char* name, const char* title, int nbin, const Double_t* bins){
		T *h = new T(name, title, nbin, bins);
		th1.push_back(dynamic_cast<TH1*>(h));
		hkey[name] = th1.back();
		return h;
}

template <typename T> T* histManager::regHist (const char* name, const char* title, int nbin, const Double_t* bins, Int_t nbin2, Double_t ylow, Double_t yup){
		T *h = new T(name, title, nbin, bins, nbin2, ylow, yup);
		th1.push_back(dynamic_cast<TH1*>(h));
		hkey[name] = th1.back();
		return h;
}

template <typename T> T* histManager::regHist (const char* name, const char* title, int nx, double x1, double x2,
				int ny, double y1, double y2){
		T *h = new T(name, title, nx, x1, x2, ny, y1, y2);
		th1.push_back(dynamic_cast<TH1*>(h));
		hkey[name] = th1.back();
		return h;
}

template <typename T> T* histManager::get (const char* name){
		return dynamic_cast<T*>(hkey[name]);
}

void histManager::sumw2(){
		for(auto x=th1.begin(); x!= th1.end(); ++x){
				(*x)->Sumw2();
		}
}
void histManager::write(){
		for(auto x=th1.begin(); x!= th1.end(); ++x) (*x)->Write();
}

#endif
