
#ifndef JTCERROR_H
#define JTCERROR_H

#include "edmJtcUtility.h"
#include "TMath.h"

class jtcError {
		public : jtcError(){}
				 ~jtcError(){};
				 void getBkgError(TH1D* h);
				 void makeDrError(TH1* hdr); 
				 void drawBkgErrorCheck();
				 jtcError& operator = (jtcError & je){
						 kEmpty = je.kEmpty;
						 bg_err = je.bg_err;
						 me_err = je.me_err;
						 bkg_level = je.bkg_level;
						 return *this;
				 }
				 
		public : bool kEmpty = 1;
				 float bg_err = 0, me_err = 0, bkg_level = 0;
				 float chi2 = 0;
				 
};

void jtcError::getBkgError(TH1D* sig_step2_proX){
        TH1D* h1= sig_step2_proX;
        int l1p5 = h1->FindBin(-1.4999)-1;
        int l2   = h1->FindBin(-2.0001);
        int r2   = h1->FindBin( 2.0);
        int r1p5 = h1->FindBin( 1.5);

        float mean = h1->GetBinContent(l1p5)+h1->GetBinContent(l2)+h1->GetBinContent(r2)+h1->GetBinContent(r1p5);
		mean = mean /4;
        float left_ave, right_ave, in_ave, out_ave;
		bkg_level = mean;
        left_ave  = (h1->GetBinContent(l1p5)+h1->GetBinContent(l2))/2;
        right_ave = (h1->GetBinContent(r1p5)+h1->GetBinContent(r2))/2;
        in_ave = (h1->GetBinContent(l1p5)+h1->GetBinContent(r1p5))/2;
        out_ave= (h1->GetBinContent(l2)+h1->GetBinContent(r2))/2;
        me_err = fmax(fabs(left_ave-mean), fabs(right_ave-mean));
        bg_err = fmax(fabs(in_ave-mean), fabs(out_ave-mean));
		kEmpty= 0;
		return;
}

void jtcError::makeDrError(TH1* hdr){
		// calculate and substitute the error by the syst error in dr, but leave the bin contant unchanged.
		if( kEmpty ){
			   	std::cout<<"please load or get the error info. first!"<<std::endl;
				return;
		}
		int ndr = hdr->GetNbinsX();
		for(int k=1; k< ndr+1; ++k){
                float ring =TMath::Pi()*(pow(hdr->GetBinLowEdge(k)+hdr->GetBinWidth(k),2)-
                                pow(hdr->GetBinLowEdge(k),2))/2/hdr->GetBinWidth(k);
                float err = ring*pow(bg_err*bg_err+me_err*me_err, 0.5);
                //  if(withbg) err = ring*pow(bg_err*bg_err+me_err*me_err, 0.5);
                //  else err    = ring*me_err;
                hdr->SetBinError(k, err);
//                hdr->SetBinContent(k, dr->GetBinContent(k));
        }
}


#endif
