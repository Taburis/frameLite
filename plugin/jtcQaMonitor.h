
#ifndef JTCQA_H
#define JTCQA_H

#include "edmJtcUtility.h"
#include "jtcTH1Player.h"
#include "ParaSet.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"

Color_t color_vec[6]={kBlue+1, kRed+1, kGreen+2, kAzure+7, kMagenta+2, kBlack};

class jtcQaMonitor{
		public : jtcQaMonitor(){
						 tl.SetLineStyle(2);
						 tx.SetTextSize(.06);
				 };
				 ~jtcQaMonitor(){
						 for(auto &it : vm2trash) delete it;
						 //						 for(auto &it : vcanvas_trash){
						 //it->Close();
						 //								 delete it;
						 //						 }
				 }
				 multi_canvas<TH1>* overlay(TString savename = "", bool drawShape = 0);
				 multi_canvas<TH1>* overlayR(TString savename = "Ra", TString opt = "");
				 multi_canvas<TH1>* drawBkgErrorCheck(TString savenname = "");
				 void flash(){
						 if(needDelete)for(auto &it : vm2pair){
								 //								 it.first->cleanAll();
								 //								 delete it.second;
						 }
						 if(needDelete){
								 //								 for(auto &it : subpad) delete it;
								 for(auto it : vm2th1) vm2trash.push_back(it);
						 }
						 if(tleg!=nullptr) delete tleg;
						 tleg=nullptr;
						 vm2th1.clear();
						 vm2th1e.clear();
						 vm2pair.clear();
						 vm2Epair.clear();
						 needDelete = 0;
						 doSystError = 0;
				 }
				 void drawLegend(){
						 if(makeLegend && tleg!=nullptr ) tleg->Draw();
						 else if(makeLegend && tleg == nullptr){
								 std::cout<<"ERROR: please call bookLegend() before drawing"<<std::endl;
						 }
				 }
				 void bookLegend(float x1=0.65, float y1= 0.55, float x2=.96, float y2=.75){
						 tleg = new TLegend(x1,y1, x2, y2);
						 tleg->SetLineColor(0);
				 }
				 void addm2TH1(matrixTH1Ptr* m2){
						 vm2th1.push_back(m2);
				 }
				 TPad* padcd(int i, int j, int isdown=0){
						 return subpad[2*i*ncol+2*j+isdown];
				 }
				 void style0(TH1*h, Color_t color){
						 // used by overlay and overlayR(upper panel);
						 h->SetLineColor(color);
						 h->SetMarkerStyle(20);
						 h->SetMarkerSize(0.8);
						 h->SetMarkerColor(color);
				 }
				 void style0_for_systError(TH1*h, Color_t color){
						 h->SetFillStyle(1001);
						 h->SetFillColorAlpha(color, 0.4);
						 h->SetMarkerStyle(20);
						 h->SetMarkerSize(0.8);
						 h->SetMarkerColor(color);
				 }
				 void upper_pad_cfg(TH1* h){
						 h->GetYaxis()->SetTitleSize(0.08);
						 h->GetYaxis()->SetLabelSize(0.08);
						 h->GetYaxis()->SetNdivisions(yndivision);
				 }
				 void lower_pad_cfg(TH1* h){
						 h->GetYaxis()->SetTitleSize(0.14);
						 h->GetYaxis()->SetLabelSize(0.15);
						 h->GetXaxis()->SetLabelSize(0.16);
						 h->GetXaxis()->SetTitleSize(0.18);
						 h->GetXaxis()->SetNdivisions(xndivision);
						 h->GetXaxis()->SetTickSize(0.08);
						 h->GetYaxis()->SetTickSize(0.06);
						 h->SetMarkerStyle(20);
						 h->SetMarkerSize(0.7);
				 }
				 void addm2TH1Pair(matrixTH1Ptr *m2num, matrixTH1Ptr *m2den);
				 void addm2TH1ErrorPair(matrixTH1Ptr *m2num, matrixTH1Ptr *m2den);
				 void setTitlePosition(float a, float b){xtitle = a; ytitle = b;}
				 void setYrange(float a, float b){y1 = a; y2 = b; fixYrange = 1;}
				 void setLowerYrange(float a, float b){yR1 = a; yR2 = b; fixRatioRange = 1;}
				 void setXrange(float a, float b){x1 = a; x2 = b; fixXrange = 1;}
				 void addhLine(float a){yline = a; drawLine = 1;}
				 void addRLine(float a){yratioLine = a; ratioLine = 1;}
				 void errorDrivenRange(TH1* h, float x01, float x02, int upscale=10, int lowscale=8){// the error range
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
						 h->SetAxisRange(mean-lowscale*dvt, mean+upscale*dvt, "Y");
				 }
				 void setXndivision(int x){xndivision = x;}
				 void setYndivision(int y){yndivision = y;}
				 void addLegendEntry(TString s, int n, TString opt = "pl"){
						 tleg->AddEntry((vm2th1[n])->at(0,0), s, opt);
						 makeLegend =1;
				 }
				 void addLegendPair(TString s1, TString s2, int n, TString opt = "pl"){
						 tleg->AddEntry(((vm2pair[n]).first)->at(0,0), s1, opt);
						 tleg->AddEntry(((vm2pair[n]).second)->at(0,0), s2, opt);
						 makeLegend =1;
				 }

		public :// ParaSet * ps;
				 jtc_utility::index2d (*pad_map) (int, int) = nullptr;
				 TString (*pad_title) (int, int) = nullptr;

				 // config to plot
				 float x1 = 1, x2 = -1;
				 float y1, y2, yR1, yR2;
				 float xtitle = 0.5, ytitle =0.85;
				 bool fixXrange = 0, fixYrange = 0, autoYrange = 1, fixRatioRange = 0;
				 //pad config
				 bool doSave = 0, makeTitle = 0, needDelete = 0, makeLegend = 0, doSystError = 0;
				 int ncol = 1, nrow = 1;
				 int npt, ncent;
				 float xline , yline, yratioLine;
				 bool drawLine = 0 , ratioLine = 0;
				 std::vector<matrixTH1Ptr*> vm2th1;
				 std::vector<matrixTH1Ptr*> vm2th1e;
				 std::vector<matrixTH1Ptr*> vm2trash;
				 std::vector<TCanvas*> vcanvas_trash;
				 std::vector<std::pair<matrixTH1Ptr*, matrixTH1Ptr*>> vm2pair;
				 std::vector<std::pair<matrixTH1Ptr*, matrixTH1Ptr*>> vm2Epair; // pair of error
				 std::vector<TPad*> subpad; 
				 // pad style config
				 int xndivision = 505, yndivision = 505;
				 TLine tl;
				 TBox box;
				 TLatex tx;  
				 TLegend* tleg = nullptr;  
};

multi_canvas<TH1>* jtcQaMonitor::overlay(TString savename, bool drawShape){
		va_list ap;
		//if( vm2th1.size() != 0) ncol = vm2th1[0]->Ncol(); nrow = vm2th1[0]->Nrow();
		auto cm = new multi_canvas<TH1>("c_"+savename, "", nrow, ncol);
		TString tmp;
		int npt   = vm2th1[0]->nrow;
		int ncent = vm2th1[0]->ncol;
		auto m2th = vm2th1[0];
		float min[npt*ncent];
		float max[npt*ncent];
		if(autoYrange){
				std::cout<<"setting auto scale for Y axis ... "<<std::endl;
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
										if(!(m2th->isValid(i,j))) continue;
										float holder = m2th->at(i,j)->GetMaximum();
										if(max[i+j*npt]< holder)  max[i+j*npt] = holder ;
										holder = m2th->at(i,j)->GetMinimum();
										if(min[i+j*npt]> holder)  min[i+j*npt] = holder ;
								}
						}
				}
		}
		if(makeTitle && pad_title == nullptr){
				std::cout<<"ERROR: please specify the pad_title function first!"<<std::endl;
				return nullptr;
		}
		//		std::cout<<"plot here"<<std::endl;
		for(int k=0; k<vm2th1.size(); ++k){
				m2th = vm2th1[k];
				for(int i=0; i<npt ; ++i){
						for(int j=0; j<ncent; ++j){
								if(!(m2th->isValid(i,j))) continue;
								float grid = (max[i+j*npt]-min[i+j*npt])/16;
								m2th->at(i,j)->SetLineColor(color_vec[k]);
								m2th->at(i,j)->SetMarkerStyle(20);
								m2th->at(i,j)->SetMarkerSize(0.8);
								m2th->at(i,j)->SetMarkerColor(color_vec[k]);
								m2th->at(i,j)->GetXaxis()->SetNdivisions(xndivision);
								m2th->at(i,j)->GetYaxis()->SetNdivisions(yndivision);
								if(fixXrange) m2th->at(i,j)->SetAxisRange(x1, x2,"X");
								if(autoYrange) {
										m2th->at(i,j)->SetAxisRange(min[i+j*npt]-1.5*grid, max[i+j*npt]+1.5*grid,"Y");
								}
								if(fixYrange) m2th->at(i,j)->SetAxisRange(y1, y2, "Y");
								if(pad_map != nullptr){ 
										auto index = pad_map(i,j);
										cm->CD(index.i1, index.i2);
								} else cm->cd(j+1);
								m2th->at(i,j)->GetXaxis()->CenterTitle();
								if(drawShape)m2th->at(i,j)->DrawNormalized("same");
								else m2th->at(i,j)->Draw("same");
								if(makeTitle)  tx.DrawLatexNDC(xtitle, ytitle, pad_title(i,j));
								if(i+j==0) drawLegend();
								if(drawLine ){
										if(autoYrange && min[i+j*npt]-1.5*grid < yline && max[i+j*npt]+1.5*grid > yline) tl.DrawLine(x1, yline, x2, yline);
										else tl.DrawLine(x1, yline, x2, yline);
								}
								//m2th->at(i,j)->DrawNormalized("same");
						}
				}
		}
		//vcanvas_trash.push_back(cm);
		if(doSave) cm->SaveAs(savename);
		//		std::cout<<"plot here"<<std::endl;
		return cm;
}

void jtcQaMonitor::addm2TH1Pair(matrixTH1Ptr *m2num, matrixTH1Ptr *m2den){
		vm2pair.push_back(std::make_pair(m2num, m2den));
		return;
}

void jtcQaMonitor::addm2TH1ErrorPair(matrixTH1Ptr *m2num, matrixTH1Ptr *m2den){
		vm2Epair.push_back(std::make_pair(m2num, m2den));
		doSystError =1;
		return;
}

multi_canvas<TH1>* jtcQaMonitor::overlayR(TString savename, TString opt){
		// the ratio will be store in vm2th1
		TString cmname = "c_"+savename;
		auto cm = new multi_canvas<TH1>(cmname, "", nrow, ncol, 325, 400);
		//		cm->Divide(ncol, nrow);
		float r = 0.325;
		for(auto &it : vm2pair ){
				matrixTH1Ptr *h;
				if(opt == "B") h = (*(it.first))%(*(it.second));
				else h = (*(it.first))/(*(it.second));
				vm2th1.push_back(h);
		}
		if(doSystError){
				for(auto &it : vm2Epair ){
						matrixTH1Ptr *h;
						if(opt == "B") h = (*(it.first))%(*(it.second));
						else h = (*(it.first))/(*(it.second));
						vm2th1e.push_back(h);
				}
		}
		npt   = vm2th1[0]->nrow;
		ncent = vm2th1[0]->ncol;
		float min[npt*ncent];
		float max[npt*ncent];
		subpad.resize(2*npt*ncent);

		for(int i=0; i<npt; i++){
				for(int j=0; j<ncent; j++){
						int i1, i2;
						if(pad_map!= nullptr){
								auto indx = pad_map(i,j);
								i1=indx.i1, i2 = indx.i2;
								cm->CD(i1, i2);
						} else {
								i1 = i; i2 = j;
								cm->cd(i1+i2+1);
						}
						//cout<<i<<","<<j<<endl;
						TString padname = Form("subpad_%d_%d_0",i,j);
						subpad[2*i1*ncol+2*i2] = new TPad(padname, "", 0.0, r, 1, 1);
						padname = Form("subpad_%d_%d_1",i,j);  
						subpad[2*i1*ncol+2*i2+1] = new TPad(padname, "", 0.0, 0, 1, r);
						subpad[2*i1*ncol+2*i2]  ->SetTopMargin(0.05);
						subpad[2*i1*ncol+2*i2]  ->SetLeftMargin(0.14);
						subpad[2*i1*ncol+2*i2+1]->SetLeftMargin(0.14);
						subpad[2*i1*ncol+2*i2]  ->SetBottomMargin(0);
						subpad[2*i1*ncol+2*i2+1]->SetTopMargin(0);
						subpad[2*i1*ncol+2*i2]  ->SetRightMargin(0.02); 
						subpad[2*i1*ncol+2*i2]  ->SetTopMargin(0.03);
						subpad[2*i1*ncol+2*i2+1]->SetRightMargin(0.02);
						subpad[2*i1*ncol+2*i2+1]->SetBottomMargin(0.5);
						padcd(i1,i2,0)->Draw(); 
						padcd(i1,i2,1)->Draw();
				}
		}

		matrixTH1Ptr* m2th1, *m2th2, *m2th1E, *m2th2E, *m2ratE;
		if(autoYrange){
				auto m2th1 = vm2pair[0].first;
				auto m2th2 = vm2pair[0].second;
				for(int i =0; i<npt ; ++i){
						for(int j=0; j<ncent; ++j){
								max[i+j*npt] = m2th1->at(i,j)->GetMaximum();
								min[i+j*npt] = m2th1->at(i,j)->GetMinimum();
						}
				}
				for(int k=0; k<vm2th1.size(); ++k){
						m2th1 = vm2pair[k].first;
						m2th2 = vm2pair[k].second;
						for(int i=0; i<npt ; ++i){
								for(int j=0; j<ncent; ++j){
										float holder = m2th1->at(i,j)->GetMaximum();
										if(max[i+j*npt]< holder)  max[i+j*npt] = holder ;
										holder = m2th2->at(i,j)->GetMaximum();
										if(max[i+j*npt]< holder)  max[i+j*npt] = holder ;
										holder = m2th1->at(i,j)->GetMinimum();
										if(min[i+j*npt]> holder)  min[i+j*npt] = holder ;
										holder = m2th2->at(i,j)->GetMinimum();
										if(min[i+j*npt]> holder)  min[i+j*npt] = holder ;
								}
						}
				}
		}

		if(makeTitle && pad_title == nullptr){
				std::cout<<"ERROR: please specify the pad_title function first!"<<std::endl;
				return nullptr;
		}

		for(int k=0; k<vm2pair.size(); ++k){
				m2th1 = vm2pair[k].first;
				m2th2 = vm2pair[k].second;
				auto m2rat = vm2th1[k];
				if(doSystError){
						m2th1E= vm2Epair[k].first;
						m2th2E= vm2Epair[k].second;
						m2ratE = vm2th1e[k];
				}
				int color_index = 0;
				for(int i=0; i<npt ; ++i){
						for(int j=0; j<ncent; ++j){
								//cout<<i<<", "<<j<<endl;
								jtc_utility::index2d index;
								float grid = (max[i+j*npt]-min[i+j*npt])/16;
								style0(m2th1->at(i,j), color_vec[color_index]);
								style0(m2th2->at(i,j), color_vec[color_index+1]);
								if(fixXrange) m2th1->at(i,j)->SetAxisRange(x1, x2,"X");
								if(autoYrange) m2th1->at(i,j)->SetAxisRange(min[i+j*npt]-1.5*grid, max[i+j*npt]+1.5*grid,"Y");
								if(fixYrange) m2th1->at(i,j)->SetAxisRange(y1, y2, "Y");
								if(pad_map != nullptr){ 
										index = pad_map(i,j);
										padcd(index.i1, index.i2, 0)->cd();
								} else padcd(i, j ,0)->cd();
								upper_pad_cfg(m2th1->at(i,j));
								m2th1->at(i,j)->Draw();
								m2th2->at(i,j)->Draw("same");

								if(doSystError){
										style0_for_systError(m2th1E->at(i,j), color_vec[color_index]);
										style0_for_systError(m2th2E->at(i,j), color_vec[color_index+1]);
										m2th1E->at(i,j)->Draw("same e2");
										m2th2E->at(i,j)->Draw("same e2");
								}

								if(i+j==0) drawLegend();
								if(makeTitle)  tx.DrawLatexNDC(xtitle, ytitle, pad_title(i,j));
								if(drawLine ){
										if(autoYrange && min[i+j*npt]-1.5*grid < yline && max[i+j*npt]+1.5*grid > yline) tl.DrawLine(x1, yline, x2, yline);
										else tl.DrawLine(x1, yline, x2, yline);
								}
								if(pad_map != nullptr){ padcd(index.i1, index.i2, 1)->cd();
								} else padcd(i, j ,1)->cd();
								lower_pad_cfg(m2rat->at(i,j));
								if(fixXrange) m2rat->at(i,j)->SetAxisRange(x1, x2,"X");
								m2rat->at(i,j)->GetXaxis()->SetTitle(m2th1->at(i,j)->GetXaxis()->GetTitle());
								if(fixRatioRange) m2rat->at(i,j)->SetAxisRange(yR1, yR2, "Y");
								m2rat->at(i,j)->Draw();
								if(doSystError){
										style0_for_systError(m2ratE->at(i,j), kGray+4);
										m2ratE->at(i,j)->Draw("same e2");
								}
								if(ratioLine){
										if(!fixXrange){
												x1=m2rat->at(i,j)->GetXaxis()->GetXmin();
												x2=m2rat->at(i,j)->GetXaxis()->GetXmax();
										}
										tl.DrawLine(x1, yratioLine, x2, yratioLine);
								}
						}
				}
				color_index = color_index+2;
		}

		cm->Update();
		needDelete = 1;
		vcanvas_trash.push_back(cm);
		return cm;
}

multi_canvas<TH1>* jtcQaMonitor::drawBkgErrorCheck(TString savename){
		// should add the projX jtcTH1Player first and this jtcTH1Player should got the bkg error already.
		setXrange(-3, 2.99);
		int npt   = vm2th1[0]->nrow;
		int ncent = vm2th1[0]->ncol;
		jtcTH1Player* m2th =(jtcTH1Player*) vm2th1[0];
		for(int i=0; i<npt ; ++i){
				for(int j=0; j<ncent; ++j){
						auto h = m2th->at(i,j);
						float bg_err = m2th->bkgError(i,j);
						float me_err = m2th->meError(i,j);
						float mean  = m2th->bkgLevel(i,j);
						float err = pow(bg_err*bg_err+me_err*me_err, .5);
						h->SetAxisRange(mean-4*err, mean+4*err, "Y");
				}
		}
		fixYrange = 0;
		autoYrange = 0;
		auto cm = overlay(savename);
		for(int i=0; i<npt ; ++i){
				for(int j=0; j<ncent; ++j){
						auto indx = pad_map(i,j);
						int i1=indx.i1, i2 = indx.i2;
						cm->CD(i1, i2);
						box.SetFillColorAlpha(kAzure+7, 0.4);
						float bg_err = m2th->bkgError(i,j);
						float me_err = m2th->meError(i,j);
						float mean  = m2th->bkgLevel(i,j);
						float err = pow(bg_err*bg_err+me_err*me_err, .5);
						auto h = m2th->at(i,j);
						float center = h->GetBinContent(h->FindBin(0));
						box.DrawBox(-2.5, mean+bg_err, 2.499, mean-bg_err);
						box.SetFillColorAlpha(kRed+1, 0.4);
						box.DrawBox(-2.5, mean+me_err, -1.5, mean-me_err);
						box.DrawBox(1.5, mean+me_err, 2.5, mean-me_err);
						tl.DrawLine(-2.5, mean, 2.5, mean);
						float percen = err/center;
						tx.DrawLatexNDC(0.3, 0.2, Form("error/center: %f",percen));
				}
		}
		return cm;
}
#endif
