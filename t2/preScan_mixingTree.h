
#ifndef PRESCAN_MIXINGTREE_H
#define  PRESCAN_MIXINGTREE_H
#include "rootEDM.cc"
#include "histManager.h"
#include "/home/wang4187/frameLite/crab/lib/configReader.h"

class preScan_mixingTree : public rootEDMProducer{
	public :
		preScan_mixingTree(){}
		preScan_mixingTree(TString inputf){ infstr = inputf;}
		void beginJob();
		void endJob();
		void analyze();
		void loadCfg(ParaSet &);
	public :
		TString infstr, outfstr;	
		TFile *wf;
		TTree *t;
		std::vector<float> *gen_jtpt=0, *gen_jteta=0, *gen_jtphi=0, *reco_jtpt= 0, *reco_jteta=0, *reco_jtphi=0;
		std::vector<float> *esm_jteta= 0, *esm_jtphi = 0;
		TH2F* h2eta, *h2phi;
		float vz, pthat;
		TH1F* hvz, *hpthat, *hpthatStat, *hrecojtpt;
		histManager* hm;
		//		TArrayF *pThatInterval, *pThatEntries;
		ParaSet *ps;
		bool (*evtCut)(float) = 0;
		bool testWeight=0, isdata=0;
		int HBHE, pvfilter;
		configReader* cr= 0 ;
		int DBCode = -1;
};

void preScan_mixingTree::analyze(){
	float weight = 1;
	if( testWeight && !isdata) weight = (cr->pThatWeight(pthat))*(cr->vzWeight(vz));	
	if(evtCut(vz)) return 0;
	if(isdata){
		if(HBHE == 0 || pvfilter == 0) return 0;
	}
	hvz->Fill(vz, weight);
	for(int i =0 ;i<reco_jtpt->size(); ++i){	
		if(reco_jtpt->at(i) < 80 || fabs(reco_jteta->at(i)) > 1.6) continue;
		hrecojtpt->Fill(reco_jtpt->at(i), weight);
		h2eta->Fill(reco_jteta->at(i), esm_jteta->at(i), weight);
		h2phi->Fill(reco_jtphi->at(i), esm_jtphi->at(i), weight);
	}
	if( DBCode == 4 || DBCode == 5) return 0;
	hpthat->Fill(pthat, weight);
	hpthatStat->Fill(pthat);
}

void preScan_mixingTree::beginJob(){

	DBCode = ps->getPara<int>("DBCode");
	auto cr = new configReader("config_bjetMc.root");
	t = handle("mixing_tree");
	if(testWeight){
		if(cr == 0)  std::cout<<"Please load the config.root for testing!"<<std::endl;	
	}
	if(isdata){
		t->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHE);
		t->SetBranchAddress("pPAprimaryVertexFilter", &pvfilter);
	}else if(DBCode != 4 && DBCode != 5){
		t->SetBranchAddress("pthat", &pthat);
	}
	t->SetBranchAddress("vz", &vz);
	//t->SetBranchAddress("genpt", &gen_jtpt);
	t->SetBranchAddress("pf_jtpt", &reco_jtpt);
	t->SetBranchAddress("pf_jteta", &esm_jteta);
	t->SetBranchAddress("pf_jtphi", &esm_jtphi);
	t->SetBranchAddress("pf_wta_eta", &reco_jteta);
	t->SetBranchAddress("pf_wta_phi", &reco_jtphi);
	hm = new histManager();
	hvz = hm->regHist<TH1F>("vz", "vz distribution", 200, -20, 20);
	if(DBCode != 4 && DBCode != 5){
		hpthat = hm->regHist<TH1F>("pthat", "pthat distribution", 200, 0, 400);
		Int_t nbin = ps->getPara<Int_t>("npThat");
		const Float_t* ptedge= ps->getPara<const Float_t*>("pThatInterval");
		hpthatStat = hm->regHist<TH1F>("pthatStat", "pthat Stats distribution", nbin, ptedge);
	}
	Int_t njtptbin = ps->getPara<Int_t>("njtpTbin");
	const Float_t* jtptedge= ps->getPara<const Float_t*>("jtpTbin");
	hrecojtpt = hm->regHist<TH1F>("reco_jtpt", "reco jet pT distribution", njtptbin, jtptedge);
	h2eta = hm->regHist<TH2F>("h2eta", "esm-eta vs wta-eta distribution", 100, -2, 2, 100, -2, 2);
	h2phi = hm->regHist<TH2F>("h2phi", "esm-phi vs wta-phi distribution", 100, -TMath::Pi(), TMath::Pi(), 100, -TMath::Pi(), TMath::Pi());
	hm->sumw2();
	evtCut = ps->getPara<bool (*)(float)>("evtCut");
}

void preScan_mixingTree::endJob(){
	wf = TFile::Open(outfstr, "recreate");
	wf->cd();
	//pThatInterval->Write();
	//pThatEntries->Write();
	hm->write();
	wf->Write();
}

void preScan_mixingTree::loadCfg(ParaSet &ps0){
	ps = &ps0;
}

#endif
