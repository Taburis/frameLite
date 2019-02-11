
#ifndef jtc_utility
#define jtc_utility
#include <iostream>
#include "TH2.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TString.h"
#include "TString.h"
//#include "matrixTObjPtr.h"
#include "matrixTH1Ptr.h"
//#include "jtcTH1Player.h"


namespace jtc_utility {

		struct index2d{
				int i1;
				int i2;
		};

		typedef index2d (*mapper_func)(int, int) ; 
		typedef TString (*title_func)(int, int) ; 

		Double_t etabin[24] ={-3.5, -3, -2.5,-2.,-1.5, -1., -0.8, -0.6, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 1.5,2.,2.5, 3, 3.5};

		Double_t phibin[18] ={-1.50796, -1.00531,-0.879646, -.75398, -0.628319,-0.502655, -0.376991, -0.251327, -0.125664, 0.125664, 0.251327, 0.376991, 0.502655, 0.628319,.75398, 0.879646, 1.00531,1.50796};

		void invariant_TH2(TH2* h){ h->Scale(1.0/h->GetXaxis()->GetBinWidth(1)/h->GetYaxis()->GetBinWidth(1));
		}

		Double_t seagull_pol2_par3(Double_t *x, Double_t *par){
				//fitting function is: a0+a1*|x-x0|+a2*x+a3*x^2;
				// total 4 parameters: npar = 4 
				float x0 = 0.28;
				if(x[0] < x0 && x[0] > -x0) return par[0];
				else if( x[0]<-x0 ) return par[1]*(x[0]+x0)+par[0]+par[3]*pow(x[0]-x0, 2);
				else return par[2]*(x[0]-x0)+par[0]+par[3]*pow(x[0]-x0, 2);
				//else if( x[0]<-x0 ) return -par[1]*(x[0]+x0)+par[0]+par[3]*pow(x[0]-x0, 2);
				//else return par[1]*(x[0]-x0)+par[0]+par[3]*pow(x[0]-x0, 2);
				//else if( x[0]<-x0 ) return -par[1]*(x[0]+x0)+par[0]+par[2]*pow(x[0], 2);
				//else return par[1]*(x[0]-x0)+par[0]+par[2]*pow(x[0], 2);
		}

		TString histNameScheme(const char * caption, bool isGenJet, bool isGenTrack, bool isPTweighted = 0, bool isMix = 0){
				TString histType = caption;
				if(isGenJet){
						histType = histType+"GenJet_";
				} else {
						histType = histType+"RecoJet_";
				}
				if(isGenTrack){
						histType = histType+"GenTrack";
				}else {
						histType = histType+"RecoTrack";
				}
				if( isPTweighted) histType = histType + "_pTweighted";
				if( isMix) histType = histType + "_mixing";
				return histType;
		}

		TH2D* sideBandMixingTableMaker(TH2D* h2, float sidemin, float sidemax){
				int s1 = h2->GetYaxis()->FindBin(sidemin);
				int s2 = h2->GetYaxis()->FindBin(sidemax-0.001);
				int dbin = s2-s1;
				TH1D* temp = (TH1D*) h2->ProjectionX("_sideMix_deta", s1, s2);
				TString name = h2->GetName();
				TH2D* ME = (TH2D*) h2->Clone(name);
				float mean=0;
				// the middle position left/right edge
				float midLeft = -0.15;
				float midRight = 0.15;
				int binLeft = temp->FindBin(midLeft);
				int binRight= temp->FindBin(midRight)+1;
				for(int i=binLeft;i<binRight; i++){
						mean += temp->GetBinContent(i);
				}
				mean = mean/(temp->FindBin(midRight)-temp->FindBin(midLeft)+1)/h2->GetNbinsY();
				temp->Scale(1.0/mean);
				for(int ix=1; ix<h2->GetNbinsX()+1; ix++){
						for(int iy=1; iy<h2->GetNbinsY()+1; iy++){
								if( ix< binRight && ix>= binLeft){
										ME->SetBinContent(ix, iy, 1);
										ME->SetBinError(ix, iy, 0);
								}
								else{
										ME->SetBinContent(ix, iy, temp->GetBinContent(ix)/dbin);
										ME->SetBinError(ix, iy, temp->GetBinError(ix)/sqrt(dbin));
								}
						}
				}
				temp->Delete();
				return ME;
		}
		TH2D* mixingTableMaker(TH2D* mix, bool doSmooth){
				float midLeft = -0.15;
				float midRight = 0.15;
				//make the mixing table invariant
				mix->Scale(1.0/mix->Integral()/mix->GetXaxis()->GetBinWidth(1)/mix->GetYaxis()->GetBinWidth(1));
				TH1D* temp = (TH1D*)mix->ProjectionX("_eta");
				TString name = mix->GetName();
				if( doSmooth) name = name+"_smoothed";
				else name = name+"_noEmptyBin";
				TH2D* ME = (TH2D*) mix->Clone(name);
				float mean=0;
				int binLeft = temp->FindBin(midLeft);
				int binRight= temp->FindBin(midRight)+1;
				for(int i=binLeft;i<binRight; i++){
						mean += temp->GetBinContent(i);
				}
				mean = mean /(temp->FindBin(midRight)-temp->FindBin(midLeft)+1)/mix->GetNbinsY();
				temp->Scale(1.0/mean);
				if(doSmooth){
						for(int ix=1; ix<mix->GetNbinsX()+1; ix++){
								for(int iy=1; iy<mix->GetNbinsY()+1; iy++){
										if( ix< binRight && ix>= binLeft){
												ME->SetBinContent(ix, iy, 1);
												ME->SetBinError(ix, iy, 0);
										}
										else{
												ME->SetBinContent(ix, iy, temp->GetBinContent(ix)/mix->GetNbinsY());
												ME->SetBinError(ix, iy, temp->GetBinError(ix)/sqrt(mix->GetNbinsY()));
										}
								}
						}
				}
				else { //only fill the empty bin of the mixing by the content we have
						ME->Scale(1.0/mean);
						for(int ix=1; ix<mix->GetNbinsX()+1; ix++){
								for(int iy=1; iy<mix->GetNbinsY()+1; iy++){
										if( ix< binRight && ix>= binLeft){
												ME->SetBinContent(ix, iy, 1);
												ME->SetBinError(ix, iy, 0);
										}
										else if(ME->GetBinContent(ix,iy)==0){
												ME->SetBinContent(ix, iy, temp->GetBinContent(ix)/mix->GetNbinsY());
												ME->SetBinError(ix, iy, temp->GetBinError(ix)/sqrt(mix->GetNbinsY()));
										}
								}
						}
				}
				delete temp;
				return ME;
		}
		TH2D* testf(){
				std::cout<<"bkg"<<std::endl;
				return 0;
		}
		TH2D* getV2Bkg(TH2D* signal, float sideMin,float sideMax){
				TString stemp = signal->GetName();
				TH2D* bkg = (TH2D*)signal->Clone("bkg_"+stemp);
				int outterRight= signal->GetXaxis()->FindBin(sideMax);
				int innerRight = signal->GetXaxis()->FindBin(sideMin);
				int outterLeft = signal->GetXaxis()->FindBin(-signal->GetXaxis()->GetBinCenter(outterRight));
				int innerLeft  = signal->GetXaxis()->FindBin(-signal->GetXaxis()->GetBinCenter(innerRight));
				TH1D* temp = (TH1D*)bkg->ProjectionY("bkg_sideBandProj", innerRight, outterRight, "e");
				temp->Add(bkg->ProjectionY("", outterLeft, innerLeft, "e"));

				TH2D* aux_bkg = (TH2D*)signal->Clone("aux_bkg");
				for(int iy=1; iy<bkg->GetNbinsY()+1; iy++){
						for(int ix=1; ix<bkg->GetNbinsX()+1; ix++){
								aux_bkg->SetBinContent(ix, iy, 1);
						}
				}
				TH1D* aux_proj = (TH1D*)aux_bkg->ProjectionY("aux_proj", innerRight, outterRight);
				aux_proj->Add(aux_bkg->ProjectionY("", outterLeft, innerLeft));
				for(int iy=1; iy<bkg->GetNbinsY()+1; iy++){
						int nbin = int(aux_proj->GetBinContent(iy));
						for(int ix=1; ix<bkg->GetNbinsX()+1; ix++){
								bkg->SetBinContent(ix, iy, temp->GetBinContent(iy)/nbin);
								bkg->SetBinError(ix, iy, temp->GetBinError(iy)/sqrt(nbin));
						}
				}
				temp->Delete();
				aux_bkg->Delete();
				aux_proj->Delete();
				return bkg;
		}
		void drIntegral(TH2D* signal, TH1D* drDist, bool isStatError){
				// this should return the invariant dr distribution (dr bin width will be divided out )
				float content;
				float error;
				float width;
				float dr;
				float xwidth, ywidth;
				drDist->Sumw2();
				for(int i=1; i<drDist->GetNbinsX()+1;i++){
						drDist->SetBinContent(i, 0);
						drDist->SetBinError(i, 0);
				}
				for(int jx=0; jx<signal->GetNbinsX(); jx++){
						for(int jy=0; jy<signal->GetNbinsY(); jy++){
								dr = sqrt( pow(signal->GetXaxis()->GetBinCenter(jx),2) +\
												pow(signal->GetYaxis()->GetBinCenter(jy),2));
								xwidth = signal->GetXaxis()->GetBinWidth(jx);
								ywidth = signal->GetYaxis()->GetBinWidth(jy);
								// integrand f(x,y)dxdy
								content = signal->GetBinContent(jx,jy)*xwidth*ywidth;
								if( content ) {
										if(isStatError){
												error = sqrt(pow(drDist->GetBinError(drDist->FindBin(dr)),2)+\
																pow(signal->GetBinError(jx,jy)*xwidth*ywidth,2));
										}
										else {
												float err = signal->GetBinError(jx,jy)*xwidth*ywidth/content;
												error = sqrt(pow(drDist->GetBinError(drDist->FindBin(dr)),2)+pow(err,2));
										}
										drDist->Fill(dr, content);
										drDist->SetBinError(drDist->FindBin(dr), error);
								}
						}
				}
				// make the histogram invariant
				for(int i=1; i<drDist->GetNbinsX()+1;i++){
						content = drDist->GetBinContent(i);
						error= drDist->GetBinError(i);
						width= drDist->GetBinWidth(i);
						drDist->SetBinContent(i, content/width);
						if( isStatError) drDist->SetBinError(i, error/width);
						else drDist->SetBinError(i, content*error/width);
				}
				return;
		}
		TH1D* drDistMaker(TH2D* signal, TString name, TString title, int nbin, const float * bins){
				TH1D* drDist = new TH1D(name, title, nbin, bins);
				drIntegral(signal, drDist, 1);
				return drDist;
		}
		TH1D* doDrIntegral(TString name , TH2D* sig, int ndr, float * drbin){
				TString title = sig->GetTitle();
				TString hname = "dr_"+name;
				auto dr_integral = drDistMaker(sig, hname, title, ndr, drbin );
				dr_integral->GetXaxis()->SetTitle("#Deltar");
				return dr_integral;
		}


		TH1* invariantRebin(TH1* h1, TString name , int n, Double_t * bins){
				// rebin the histogram based on the bins given in the parameter
				if(n == h1->GetNbinsX() ) return h1;
				TH1* h=(TH1*) h1->Clone("tmp");
				//input h needs to be invariant
				for(int i=1; i<h->GetNbinsX()+1; ++i){
						Double_t wd = h->GetBinWidth(i);
						h->SetBinContent(i, h->GetBinContent(i)*wd);
						h->SetBinError  (i, h->GetBinError(i)*wd);
				}
				TH1* hh = h->Rebin(n, name, bins);
				for(int i=1; i<hh->GetNbinsX()+1; ++i){
						Double_t wd = hh->GetBinWidth(i);
						if(!hh->GetBinContent(i)) continue; //skip empty bin
						hh->SetBinContent(i, hh->GetBinContent(i)/wd);
						hh->SetBinError  (i, hh->GetBinError(i)/wd);
				}
				delete h;
				return hh;
		}

		float range_based_on_error(TH1 &h, float x1, float x2){
				int n1 = h.GetXaxis()->FindBin(x1);
				int n2 = h.GetXaxis()->FindBin(x2);
				int nsum = fabs(n2-n1)+1;
				float err =0;
				for(int i=n1; i< n2+1; i++){
						err += h.GetBinError(i);
				}
				return err/nsum;
		}
		TH1* projectionX(TH2D* h, float x, float y, TString opt){
				int xbin = h->GetYaxis()->FindBin(x);
				int ybin = h->GetYaxis()->FindBin(y);
				ybin = h->GetYaxis()->FindBin(y-h->GetYaxis()->GetBinWidth(ybin));
				TString name = h->GetName(); name = "projX_"+name;
				return h->ProjectionX(name, xbin, ybin , opt);
		}

		TH1* projectionY(TH2D* h, float x, float y, TString opt){
				int xbin = h->GetXaxis()->FindBin(x);
				int ybin = h->GetXaxis()->FindBin(y)-1;
				ybin = h->GetXaxis()->FindBin(y-h->GetXaxis()->GetBinWidth(ybin));
				TString name = h->GetName(); name = "projY_"+name;
				return h->ProjectionY(name, xbin, ybin ,opt);
		}

		TH1* projX(bool doRebin, TH2D *h2, float x, float y, TString opt){
				// here h2 needs to be invariant
				TH1* h=jtc_utility::projectionX(h2, x, y, opt);
				h->Scale(h2->GetYaxis()->GetBinWidth(1));
				if(doRebin){
						TString name = h->GetName();
						name = "rebined_"+name;
						TH1* hh=h;
						h=jtc_utility::invariantRebin(h,name, 23, etabin);
						delete hh;
				}
				return h;
		}

		TH1* projY(bool doRebin, TH2D*h2, float x, float y, TString opt){
				// here h2 needs to be invariant
				TH1* h=projectionY(h2, x, y, opt);
				h->Scale(h2->GetXaxis()->GetBinWidth(1));
				if(doRebin){
						TString name = h->GetName();
						name = "rebined_"+name;
						TH1* hh=h;
						h=invariantRebin(h,name, 17, phibin);
						delete hh;
				}
				return h;
		}

		void ring_corr( TH2D* h, TH1D* corr, float drmax = 1){
				for(int k=1; k<h->GetNbinsX()+1; ++k){
						for(int l=1; l<h->GetNbinsY()+1; ++l){
								float deta = h->GetXaxis()->GetBinCenter(k);
								float dphi = h->GetYaxis()->GetBinCenter(l);
								float dr = pow(deta*deta+dphi*dphi, 0.5);
								if(dr > drmax) continue;
								int nn = corr->GetXaxis()->FindBin(dr);
								h->SetBinContent(k,l, h->GetBinContent(k,l)/corr->GetBinContent(nn));
						}
				}
				return ;
		}

		matrixTH1Ptr * getDr(const char * name, matrixTH1Ptr * m2h, int ndr, float * drbin){
				auto m2 = new matrixTH1Ptr(name, int(m2h->Nrow()), int(m2h->Ncol()));
				for(int i=0; i<m2h->Nrow(); i++){
						for(int j=0; j<m2h->Ncol(); j++){
								auto h = doDrIntegral(TString(Form("%s_%d_%d", name, i,j)), (TH2D*)m2h->at(i,j),ndr, drbin );
								m2->add(h, i,j );
						}
				}
				m2->doFree = 1;
				return m2;
		}
}


#endif
