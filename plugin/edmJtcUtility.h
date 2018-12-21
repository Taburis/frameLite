
namespace utility {

		void invariant_TH2(TH2* h){
				h->Scale(1.0/h->GetXaxis()->GetBinWidth(1)/h->GetYaxis()->GetBinWidth(1));
		}

		TH2D* sideBandMixingTableMaker(TH2D* h2, float sidemin, float sidemax){
				int s1 = h2->GetYaxis()->FindBin(sidemin);
				int s2 = h2->GetYaxis()->FindBin(sidemax-0.001);
				int dbin = s2-s1;
				TH1D* temp = (TH1D*) h2->ProjectionX("_sideMix_deta", s1, s2);
				TString name = h2->GetName();
				TH2D* ME = (TH2D*) h2->Clone(name);
				float mean=0;
				int binLeft = temp->FindBin(midLeft);
				int binRight= temp->FindBin(midRight)+1;
				for(int i=binLeft;i<binRight; i++){
						mean += temp->GetBinContent(i);
				}
				mean = mean /(temp->FindBin(midRight)-temp->FindBin(midLeft)+1)/h2->GetNbinsY();
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
		TH2D* getV2Bkg(TH2D* signal, float sideMin,float sideMax){
				TString stemp = signal->GetName();
				TH2D* bkg = (TH2D*)signal->Clone(stemp+"_bkg");
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
		TH1D* doDrIntegral(TString name , TH2D* sig){
				TString title = sig->GetTitle();
				TString hname = "dr_"+name;
				dr_integral = drDistMaker(sig, hname, title, ndr, drbin );
				dr_integral->GetXaxis()->SetTitle("#Deltar");
				return dr_integral;
		}
		TH1D* drDistMaker(TH2D* signal, TString name, TString title, int nbin, const float * bins){
				TH1D* drDist = new TH1D(name, title, nbin, bins);
				drIntegral(signal, drDist);
				return drDist;
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

		Double_t seagull_pol2_par3(Double_t *x, Double_t *par){
				//fitting function is: a0+a1*|x-x0|+a2*x+a3*x^2;
				float x0 = 0.5;
				if(x[0] < x0 && x[0] > -x0) return par[0];
				else if( x[0]<-x0 ) return -par[1]*(x[0]+x0)+par[0]+par[2]*pow(x[0]+x0, 2);
				else return par[1]*(x[0]-x0)+par[0]+par[2]*pow(x[0]-x0, 2);
		}
}

