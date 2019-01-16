
#ifndef EDMJTCSIGNALPRODUCER_H
#define EDMJTCSIGNALPRODUCER_H
#include "edmJtcUtility.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"

class edmJtcSignalProducer {
		public : edmJtcSignalProducer(){};
				 virtual ~edmJtcSignalProducer(){
						 func_seagull->Delete();
						 //if( ndrbin != 0) delete drbins;
				 }
				 void addSig(TH2D* h1, TH2D* h2){
						 raw_sig = h1; mix = h2;
						 name = h1->GetName();
				 }
				 TH2D* getSignal(TString);
				 float getMean(TH1* h, float x1, float x2){
						 int n1 = h->GetXaxis()->FindBin(x1);
						 int n2 = h->GetXaxis()->FindBin(x2);
						 return h->Integral(n1, n2)/fabs(n2-n1+1);
				 }
				 void fixSeagull(TH2D* hsig);
				 void setFitFunction(int npar = 2, Double_t (*fcn)(Double_t *, Double_t*) = 0){
						 if( fcn == 0) fcn = jtc_utility::seagull_pol2_par3;
						 if( func_seagull !=0 ) func_seagull->Delete();
						 func_seagull = new TF1("func_"+name, fcn, -3., 3., npar);
				 }
				 void check_seagull();

				 TString name;
				 void write();
				 void produce();
				 TH2D* raw_sig = 0;	 
				 TH2D* sig = 0;	 
				 TH2D* sig_step2 = 0;	 //raw signal after acceptance correction 
				 TH2D* bkg = 0;	 
				 TH2D* mix = 0;	 
				 TH2D* mix_norm = 0;	 
				 TH1D* dr_shape = 0;
				 TH1D* dr_shape_step2 = 0; // dr integral over sig_step2
				 TH1D* dr0 = 0; // dr integral over raw_sig
				 // for seagull fitting
				 TF1 * func_seagull = 0;
				 TH1D* hsideband = nullptr;
				 //				 Double_t  (*func) (Double_t *, Double_t *) = jtc_utility::seagull_pol2_par3;
				 float sideMin = 1.5, sideMax = 2.5;
				 bool doSmoothME = true;
				 int ndrbin = 0;
				 float* drbins;
				 bool doSeagullCorr = 0;
				 bool dodrIntegral=1;
};

TH2D* edmJtcSignalProducer::getSignal(TString name){
		//		raw_sig=(TH2D*) raw_sig->Clone("raw_"+name);
		raw_sig->Scale(1.0/raw_sig->GetXaxis()->GetBinWidth(1)/raw_sig->GetYaxis()->GetBinWidth(1)); //make the h2 invariant
		sig = (TH2D*) raw_sig->Clone("signal_"+name);
		sig->GetXaxis()->SetTitle("d#eta");
		sig->GetYaxis()->SetTitle("d#phi");
		sig->GetXaxis()->CenterTitle();
		sig->GetYaxis()->CenterTitle();
		if(mix == 0) {
				std::cout<<"no mixing table has been assigned!"<<std::endl;
				return 0;
		}
		else {
				mix_norm= jtc_utility::mixingTableMaker(mix, doSmoothME);
		}
		mix_norm->SetName("smoothed_mixing_"+name);
		sig->Divide(mix_norm);
		if(doSeagullCorr)   fixSeagull(sig);
		sig_step2 = (TH2D*) sig->Clone("sig_mix_corrected_"+name);
		bkg = (TH2D*) jtc_utility::getV2Bkg(sig,sideMin , sideMax );
		bkg->SetName("bkg_"+name);
		sig->Add(sig, bkg, 1, -1);
		//		if(shiftSignal){
		//				getSignal_phiSideBand("sideBand");
		//				float mean = getMean(side_deta, -0.5, 0.5);
		//				std::cout<<"mean = "<<mean<<std::endl;
		//				std::cout<<"name : "<<name<<std::endl;
		//				shiftTH2D(sig, -mean);
		//		}
		return sig;
}

void edmJtcSignalProducer::produce(){
		TString info1 = "pull signal from input: '"+name+"'";
		TString info2;
		if(	mix ==0 ) info2 = " withOUT mix...";
		else info2 = " with mix: '"+TString(mix->GetName())+"'";
		std::cout<<info1+info2<<std::endl;
		getSignal(name);	
		if(dodrIntegral){
				dr0 = jtc_utility::doDrIntegral("raw_"+name, raw_sig, ndrbin, drbins);
				dr_shape_step2 = jtc_utility::doDrIntegral("mix_correctd_"+name, sig_step2, ndrbin, drbins);
				dr_shape = jtc_utility::doDrIntegral("signal_"+name, sig, ndrbin, drbins);
		}
}

void edmJtcSignalProducer::write(){
		if(sig!=0) sig->Write();
		if(sig_step2 != 0) sig_step2->Write();
		if(mix_norm!=0) mix_norm->Write();
		if(dr0 !=0) dr0->Write();
		if(dr_shape !=0) dr_shape->Write();
		if(dr_shape_step2 !=0) dr_shape_step2->Write();
}

void edmJtcSignalProducer::fixSeagull(TH2D* hsig){
		int n1 = hsig->GetYaxis()->FindBin(1.4);
		int n2 = hsig->GetYaxis()->FindBin(1.8);
		hsideband = (TH1D*) hsig->ProjectionX("side_band_"+name, n1, n2);
		hsideband->Scale(1.0/(n2-n1+1));
		hsideband->Rebin(6); hsideband->Scale(1.0/6);
		float mean = getMean(hsideband, -0.5, 0.5);
		hsideband->SetAxisRange(-3, 3, "X");
		hsideband->Draw();
		TLine line; line.DrawLine(-3, mean, 3, mean);
		std::cout<<"fitting function: "<<func_seagull->GetName()<<std::endl;
		std::cout<<func_seagull<<std::endl;
		auto ptr = hsideband->Fit(func_seagull, "", "", -3., 2.99);
		auto res = ptr.Get();
		auto err = res->Errors();
		cout<<"ERROR: "<<err[0]<<endl;
		//mean = func_seagull->GetParameter(0);
		//htm1->Fit(func);
		for(int i=1; i<hsig->GetNbinsX()+1; ++i){
				float x = hsig->GetXaxis()->GetBinCenter(i);
				//              if(fabs(x)<0.3) continue;
				float corr = mean/func_seagull->Eval(x);
				for(int j=1; j<hsig->GetNbinsY()+1; ++j){
						if(hsig->GetBinContent(i,j) == 0) continue;
						float cont = hsig->GetBinContent(i,j);
						hsig->SetBinContent(i,j, cont*corr);
				}
		}
		/*
		*/
		return;
}

void edmJtcSignalProducer::check_seagull(){
		float sidebandmin = 1.4, sidebandmax = 1.8;
		TH1 *h1, *h2;
		h1 = (TH1D*) jtc_utility::projX(1, sig, -1, 1, "e"); 
		h2 = jtc_utility::projX(1, sig, sidebandmin, sidebandmax, "e"); 
		std::cout<<h2->GetName()<<std::endl;
		//histStyle(h2);
		float mean = h2->GetBinContent(h1->FindBin(0));
		float dvt = jtc_utility::range_based_on_error(*h2, -2.5, -2);
		h2->SetAxisRange(mean-6*dvt, mean+10*dvt, "Y");
		h1->SetLineWidth(1);
		h2->SetLineWidth(1);
		h1->SetLineColor(kBlack);
		h2->SetLineColor(kOrange+7);
		h2->SetAxisRange(-2.5, 2.49, "X");
		//        h2->Draw();
		h1->Draw("same");
		TLine* l = new TLine(); l->SetLineStyle(2);
		l->DrawLine(-2.5, 0, 2.5, 0);
		delete h1; 
		delete h2;
}

#endif
