
#ifndef JTCQA_H
#define JTCQA_H

#include "matrixTObjPtr.h"
#include "edmJtcUtility.h"
#include "ParaSet.h"

Color_t color_vec[6]={kBlue+1, kRed+1, kGreen+2, kAzure+7, kMagenta+2, kBlack};

class jtcQaMonitor{
		public : jtcQaMonitor(){
				 };
				 multi_canvas<TH1>* overlay(TString savename);
				 void flash(){
						 vm2th1.clear();
				 }
				 void addm2TH1(matrixTH1Ptr* m2){
						 vm2th1.push_back(m2);
				 }
				 void setTitlePosition(float a, float b){xtitle = a; ytitle = b;}

		public :// ParaSet * ps;
				 jtc_utility::index2d (*pad_map) (int, int);
				 TString (*pad_title) (int, int) = 0;

				 // config to plot
				 float x1, x2;
				 float y1, y2;
				 float xtitle = 0.5, ytitle =0.85;
				 bool fixYrange = 0;
				 //pad config
				 bool doSave = 0, makeTitle = 0;
				 int ncol = 3, nrow =2;
				 int npt, ncent;
				 vector<matrixTH1Ptr*> vm2th1; 
};

multi_canvas<TH1>* jtcQaMonitor::overlay(TString savename){
		va_list ap;
		auto cm = new multi_canvas<TH1>("c_"+savename, "", nrow, ncol);
        auto tx = new TLatex();  tx->SetTextSize(.06);
        auto tl = new TLine(); tl->SetLineStyle(2);
        TString tmp;
		int npt   = vm2th1[0]->nrow;
		int ncent = vm2th1[0]->ncol;
        float min[npt*ncent];
        float max[npt*ncent];
        for(int i=0; i<npt ; ++i){
                for(int j=0; j<ncent; ++j){
                        min[i+npt*j] = 0;
                        max[i+npt*j] = 0;
                }
        }
		auto m2th = vm2th1[0];
        for(int i =0; i<npt ; ++i){
                for(int j=0; j<ncent; ++j){
                        max[i+j*npt] = m2th->at(i,j)->GetMaximum();
                        min[i+j*npt] = m2th->at(i,j)->GetMinimum();
                }
        }
		for(int k=1; k<vm2th1.size(); ++k){
				m2th = vm2th1[k];
                for(int i=0; i<npt ; ++i){
                        for(int j=0; j<ncent; ++j){
                                float holder = m2th->at(i,j)->GetMaximum();
                                if(max[i+j*npt]< holder)  max[i+j*npt] = holder ;
                                holder = m2th->at(i,j)->GetMinimum();
                                if(min[i+j*npt]> holder)  min[i+j*npt] = holder ;
                        }
                }
        }
		if(makeTitle && pad_title != 0){
				std::cout<<"ERROR: please specify the pad_title function first!"<<std::endl;
		}
		for(int k=0; k<vm2th1.size(); ++k){
				m2th = vm2th1[k];
                for(int i=0; i<npt ; ++i){
                        for(int j=0; j<ncent; ++j){
                                float grid = (max[i+j*npt]-min[i+j*npt])/16;

                                m2th->at(i,j)->SetLineColor(color_vec[k]);
                                m2th->at(i,j)->SetMarkerStyle(20);
                                m2th->at(i,j)->SetMarkerSize(0.8);
                                m2th->at(i,j)->SetMarkerColor(color_vec[k]);
                                m2th->at(i,j)->GetXaxis()->SetNdivisions(505);
                                m2th->at(i,j)->SetAxisRange(x1, x2,"X");
                                if(fixYrange) m2th->at(i,j)->SetAxisRange(y1, y2,"Y");
                                else m2th->at(i,j)->SetAxisRange(min[i+j*npt]-1.5*grid, max[i+j*npt]+1.5*grid,"Y");
								auto index = pad_map(i,j);
								cout<<index.i1<<", "<<index.i2<<endl;
								cm->CD(index.i1, index.i2);
                                m2th->at(i,j)->Draw("same");
								if(makeTitle)  tx->DrawLatexNDC(xtitle, ytitle, pad_title(i,j));
                                //m2th->at(i,j)->DrawNormalized("same");
						}
				}
		}
		if(doSave) cm->SaveAs(savename);
		return cm;
}

#endif
