
#ifndef JTCQA_H
#define JTCQA_H

#include "matrixTObjPtr.h"
#include "edmJtcUtility.h"
#include "ParaSet.h"
#include "TLatex.h"
#include "TLine.h"

Color_t color_vec[6]={kBlue+1, kRed+1, kGreen+2, kAzure+7, kMagenta+2, kBlack};

class jtcQaMonitor{
		public : jtcQaMonitor(){
				 };
				 multi_canvas<TH1>* overlay(TString savename = "");
				 void flash(){
						 vm2th1.clear();
				 }
				 void addm2TH1(matrixTH1Ptr* m2){
						 vm2th1.push_back(m2);
				 }
				 void setTitlePosition(float a, float b){xtitle = a; ytitle = b;}
				 void setYrange(float a, float b){y1 = a; y2 = b; fixYrange = 1;}
				 void addhLine(float a){yline = a; drawLine = 1;}
				 void errorDrivenRange(TH1* h, float x01, float x02){// the error range
						 int n1 = h->GetXaxis()->FindBin(x01);
						 int n2 = h->GetXaxis()->FindBin(x02);
						 float mean = h->Integral(n1, n2)/(n2-n1+1);
						 int nsum = fabs(n2-n1)+1;
						 float err =0;
						 for(int i=n1; i< n2+1; i++){
								 err += h->GetBinError(i);
						 }
						 float dvt = err/nsum;
						 //cout<<"mean: "<<mean<<", dvt: "<<dvt<<endl;
						 h->SetAxisRange(mean-5*dvt, mean+10*dvt, "Y");
				 }
				 void setXndivision(int x){xndivision = x;}
				 void setYndivision(int y){yndivision = y;}

		public :// ParaSet * ps;
				 jtc_utility::index2d (*pad_map) (int, int);
				 TString (*pad_title) (int, int) = 0;

				 // config to plot
				 float x1, x2;
				 float y1, y2;
				 float xtitle = 0.5, ytitle =0.85;
				 bool fixYrange = 0, autoYrange = 1;
				 //pad config
				 bool doSave = 0, makeTitle = 0;
				 int ncol = 3, nrow =2;
				 int npt, ncent;
				 float xline , yline;
				 bool drawLine = 0;
				 std::vector<matrixTH1Ptr*> vm2th1; 
				 // pad style config
				 int xndivision = 505, yndivision = 505;
};

multi_canvas<TH1>* jtcQaMonitor::overlay(TString savename){
		va_list ap;
		auto cm = new multi_canvas<TH1>("c_"+savename, "", nrow, ncol);
		TLatex tx;  tx.SetTextSize(.06);
		TLine tl; tl.SetLineStyle(2);
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
		if(autoYrange){
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
		}
		if(makeTitle && pad_title == 0){
				std::cout<<"ERROR: please specify the pad_title function first!"<<std::endl;
				return nullptr;
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
								m2th->at(i,j)->GetXaxis()->SetNdivisions(xndivision);
								m2th->at(i,j)->GetYaxis()->SetNdivisions(yndivision);
								m2th->at(i,j)->SetAxisRange(x1, x2,"X");
								if(autoYrange) m2th->at(i,j)->SetAxisRange(min[i+j*npt]-1.5*grid, max[i+j*npt]+1.5*grid,"Y");
								if(fixYrange) m2th->at(i,j)->SetAxisRange(y1, y2, "Y");
								auto index = pad_map(i,j);
								//std::cout<<index.i1<<", "<<index.i2<<std::endl;
								cm->CD(index.i1, index.i2);
								m2th->at(i,j)->GetXaxis()->CenterTitle();
								m2th->at(i,j)->Draw("same");
								if(makeTitle)  tx.DrawLatexNDC(xtitle, ytitle, pad_title(i,j));
								if(drawLine ){
										if(autoYrange && min[i+j*npt]-1.5*grid < yline && max[i+j*npt]+1.5*grid > yline) tl.DrawLine(x1, yline, x2, yline);
										else tl.DrawLine(x1, yline, x2, yline);
								}
								//m2th->at(i,j)->DrawNormalized("same");
						}
				}
		}
		if(doSave) cm->SaveAs(savename);
		return cm;
}

#endif
