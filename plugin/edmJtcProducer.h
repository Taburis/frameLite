
#include "histManager.h"
#include "edmMCBase.h"
#include "ParaSet.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include <algorithm>
#include "jetSet.h"
#include "xAxis.h"

//using for producing the inlcusive jet-track correlation on CMS detector MC
//can be customized to adapt to different type of jet-track correlation demands.


class edmJtcProducer: public edmMCBase{
		public : 
				struct histCase{
						jetSet * set;
						double weight = 0;
						// current jet axis 
						float jtpt, jteta, jtphi;
						// current particle axis 
						float ppt, peta, pphi;
						TH2D** sig;
						TH2D** sig_pTweighted;
						TH1D** jet_pt;
						TH1D** jet_eta;
						TH1D** jet_phi;
				};

				edmJtcProducer(ParaSet & s);
				virtual ~edmJtcProducer(){
						std::cout<<"deleting edmJtcProducer"<<std::endl;
						delete hm;
				};
				void regHist(std::string cap, histCase &hc, jetSet *jts = 0);
				void removeHC(histCase &hc){
						delete hc.sig;
						delete hc.sig_pTweighted;
						delete hc.jet_pt;
						delete hc.jet_eta;
						delete hc.jet_phi;
				}
				virtual void beginJob(); virtual void endJob();
				virtual void analyzer();
				void resetHC( histCase & hc){
						hc.weight = 0; 
				}
				void fillJetInfo(histCase &hc);
				void loadJet(histCase &hc, int j){
						hc.weight = weight; 
						hc.jtpt = hc.set->jt_pt ->at(j);
						hc.jteta= hc.set->jt_eta->at(j);
						hc.jtphi= hc.set->jt_phi->at(j);
				}
				void loadParticle(histCase &hc, unsigned int j){
						hc.ppt = pt->at(j);
						hc.peta= eta->at(j);
						hc.pphi= phi->at(j);
				}
				double findDr(double eta1, double phi1, double eta2, double phi2){
						double dphi = phi1 - phi2;
						while (dphi > M_PI) dphi-= 2*M_PI;
						while (dphi <= -M_PI) dphi += 2*M_PI;
						return sqrt(pow(dphi,2)+pow(eta1-eta2,2));
				}
				bool jetAcceptanceCut(histCase &hc){
						if(hc.jtpt<120) return 1;
						if(fabs(hc.jteta)>1.6) return 1;
						return 0;
				}
				bool particleAcceptanceCut(histCase &hc){
						if(hc.ppt<1) return 1;
						if(fabs(hc.peta)>2.4) return 1;
						return 0;
				}

				void jtc_wf001(histCase &hc);
				void fillCase(histCase & hc){
						float deta;
						float dphi;
						deta = hc.jteta - hc.peta;
						dphi = hc.jtphi - hc.pphi;
						while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
						while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}
						if(hc.ppt>400 ) hc.ppt=399;
						int i = ptax.findBin(hc.ppt);
						for(int j=0; j<nCent; ++j){
								//cout<<"index: "<<(i+j*nPt)<<" which jetpt: "<<hc.ppt<<endl;
								hc.sig[i+j*nPt]->Fill(deta, dphi, hc.weight);
								hc.sig_pTweighted[i+j*nPt]->Fill(deta, dphi, hc.weight*hc.ppt);
						}
				}

		public : 
				// variables to setup the histograms--------------
				Double_t *jet_pt_bin;
				int n_jet_pt_bin;
				std::string * ptlabel;
				int nPt, nCent;
				TString output_name;
				//------------------------------------------------
				xAxis ptax;
				histManager *hm;
				// define the histCase for each type of jet you need to run
				histCase hc; 
};

edmJtcProducer::edmJtcProducer(ParaSet & s) :
		jet_pt_bin(s.getVectorAsArray<Double_t>("jet_pt_bin")),
		n_jet_pt_bin(s.length("jet_pt_bin")-1),
		nPt(s.getPara<int>("npt")),
		nCent(s.getPara<int>("ncent"))
{
		hm = new histManager();
		float ptbins[7] = {1, 2, 3, 4, 8, 12, 400};
		ptax.setAxis(nPt, ptbins);
		output_name = s.getPara<TString>("output_name");
}


void edmJtcProducer::regHist(std::string cap, histCase &hc, jetSet *jts){
		hc.set = jts;
		std::cout<<"register the histCase: "<<cap<<std::endl;
		float ptbins     [7] = {1, 2, 3, 4, 8, 12, 400};
		TString ptLabel  [7] = {"Trk1", "Trk2", "Trk3", "Trk4", "Trk8", "trk12", "Trk400"};
		float centbins   [3] = {0, 60, 200};
		TString centLabel[3] = {"Cent0", "Cent30", "Cent100"};
		int nHistoBinsX = 500;
		int nHistoBinsY = 200;
		hc.sig= new TH2D*[nPt*nCent];
		hc.sig_pTweighted= new TH2D*[nPt*nCent];
		TString tmp, name;
		name = cap;
		//		std::cout<<"creating hists"<<std::endl;
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						tmp = centLabel[j]+"_"+centLabel[j+1]+"_"+ptLabel[i]+"_"+ptLabel[i+1];
						hc.sig[i+j*nPt] = hm->regHist<TH2D>(name+Form("_%d_%d",i, j), tmp,
										nHistoBinsX,-5,5,nHistoBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
						hc.sig_pTweighted[i+j*nPt] = hm->regHist<TH2D>(name+Form("_pTweighted_%d_%d",i, j), tmp,
										nHistoBinsX,-5,5,nHistoBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
				}
		}
		//		const float newbin [21] = {110, 120, 136, 152, 168, 184, 200, 216, 232, 248, 264, 280, 296, 312, 328, 344, 360,
		//				380, 400, 432, 500};
		//		int nbin = 20;
		hc.jet_pt = new TH1D*[nCent];
		hc.jet_eta = new TH1D*[nCent];
		hc.jet_phi = new TH1D*[nCent];

		for(int j=0; j<nCent; ++j){
				tmp = centLabel[j]+" to "+centLabel[j+1];
				hc.jet_pt[j] = hm->regHist<TH1D>(name+Form("_pt_%d",j), tmp, n_jet_pt_bin, jet_pt_bin);
				hc.jet_eta[j] = hm->regHist<TH1D>(name+Form("_eta_%d",j), tmp, 100, -2.0, 2.0);
				hc.jet_phi[j] = hm->regHist<TH1D>(name+Form("_phi_%d",j), tmp, 72, -TMath::Pi(), TMath::Pi());
		}
		return ;
}
void edmJtcProducer::beginJob(){
		//handle the Gen particles and jets
		handleGenParticle("particle_tree");
		auto jes4 = handleJetSet("ak4ES", 1);
		regHist("ak4Escheme", hc, jes4);
		hm->sumw2();
}

void edmJtcProducer::endJob(){
		auto wf = TFile::Open(output_name, "recreate");
		wf->cd();
		hm->write();
}

void edmJtcProducer::fillJetInfo(histCase &hc){
		for(int j=0; j<nCent; ++j){
				hc.jet_pt[j] ->Fill(hc.jtpt, hc.weight);
				hc.jet_eta[j]->Fill(hc.jteta,hc.weight);
				hc.jet_phi[j]->Fill(hc.jtphi,hc.weight);
		}
}

void edmJtcProducer::analyzer(){
		if(pthat < 80) return;
		jtc_wf001(hc);
}

void edmJtcProducer::jtc_wf001(histCase &hc){
		//the standared workflow to produce the jtc signal
		for(unsigned int j=0; j<hc.set->jt_pt->size(); ++j){
				loadJet(hc, j);
				if(jetAcceptanceCut(hc)) continue;
				fillJetInfo(hc);
				for(int k=0; k<pt->size(); ++k){
						if(chg->at(k) == 0) continue;
						loadParticle(hc, k);
						if(particleAcceptanceCut(hc)) continue;
						fillCase(hc);
				}
		}

}
