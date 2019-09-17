
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
		float vz, pthat;
		TH1F* hvz, *hpthat, *hpthatStat, *hrecojtpt;
		histManager* hm;
//		TArrayF *pThatInterval, *pThatEntries;
		ParaSet *ps;
		bool (*evtCut)(float) = 0;
		bool testWeight=0, isdata=0;
		int HBHE, pvfilter;
		configReader* cr= 0 ;
};

void preScan_mixingTree::analyze(){
	if(evtCut(vz)) return 0;
	float weight = 1;
	if(isdata){
		if(HBHE == 0 || pvfilter == 0) return 0;
	} else if( testWeight){
		weight = (cr->pThatWeight(pthat))*(cr->vzWeight(vz));	
	}
	hvz->Fill(vz, weight);
	if(!isdata){
		hpthat->Fill(pthat, weight);
		hpthatStat->Fill(pthat, weight);
	}
	for(int i =0 ;i<reco_jtpt->size(); ++i){	
		hrecojtpt->Fill(reco_jtpt->at(i), weight);
	}	
}

void preScan_mixingTree::beginJob(){
	t = handle("mixing_tree");
	if(testWeight){
		if(cr == 0)  std::cout<<"Please load the config.root for testing!"<<std::endl;	
	}
	if(isdata){
		t->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHE);
		t->SetBranchAddress("pPAprimaryVertexFilter", &pvfilter);
	}else {
		t->SetBranchAddress("pthat", &pthat);
	}
	t->SetBranchAddress("vz", &vz);
	//t->SetBranchAddress("genpt", &gen_jtpt);
	t->SetBranchAddress("pf_jtpt", &reco_jtpt);
	hm = new histManager();
	hvz = hm->regHist<TH1F>("vz", "vz distribution", 200, -20, 20);
	hpthat = hm->regHist<TH1F>("pthat", "pthat distribution", 200, 0, 400);
	Int_t nbin = ps->getPara<Int_t>("npThat");
	const Float_t* ptedge= ps->getPara<const Float_t*>("pThatInterval");
	hpthatStat = hm->regHist<TH1F>("pthatStat", "pthat Stats distribution", nbin, ptedge);
	Int_t njtptbin = ps->getPara<Int_t>("njtpTbin");
	const Float_t* jtptedge= ps->getPara<const Float_t*>("jtpTbin");
	hrecojtpt = hm->regHist<TH1F>("reco_jtpt", "reco jet pT distribution", njtptbin, jtptedge);
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
