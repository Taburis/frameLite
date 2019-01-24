
#include "edmJtcAnalyzer.h" 

class bjtcAnalyzer_Step3 : public rootEDMAnalyzer {
		public : bjtcAnalyzer_Step3(ParaSet & pset) : rootEDMAnalyzer(){
						 _nodeClassName_ = "bjtcAnalyzer_Step3";
						 ps = &pset;
						 loadStep2OutputFiles();
				 }
				 void get_res_corr();
				 void loadStep2OutputFiles();
				 matrixTH1Ptr* getPhiProjection(matrixTH1Ptr &m, float x, float y, TString opt = "e"){
						 TString name = m.name;
						 auto m2 = new matrixTH1Ptr(TString("phisb_"+name), m.Nrow(), m.Ncol());
						 for(auto i=0; i<m.ref.size() ; ++i) {
								 auto h = jtc_utility::projX(1, (TH2D*)m.ref[i], x, y, opt);
								 (m2->ref).at(i) = (TH1*)h;
						 }
						 m2->doFree =1;
						 return m2;
				 }
				 void recoJetCheck();
				 void recoTrackCheck();
				 void bjtc_jff();
				 int analyze(){return 0;}
		public: 
				 int npt= 6, ncent=1;
				 ParaSet *ps;
				 TString inputPath;
				 TFile *fbMC_gg, *fbMC_rg, *fbMC_rr, *fdMC_gg, *fdMC_rg, *fdMC_rr, *fcorr;
				 jtc_utility::index2d (*pad_map) (int, int);
};

void bjtcAnalyzer_Step3::loadStep2OutputFiles(){
		TString path = ps->getPara<TString>("step2output_folder");
		TString cap = ps->getPara<TString>("bjetMC_step2output_name");
		TString file = jtc_utility::histNameScheme("_", 0, 1)+".root";
		fbMC_rg = TFile::Open(path+cap+file);
		file = jtc_utility::histNameScheme("_", 1, 1)+".root";
		fbMC_gg = TFile::Open(path+cap+file);
		file = jtc_utility::histNameScheme("_", 0, 0)+".root";
		fbMC_rr = TFile::Open(path+cap+file);
		//dijet MC loading section
		cap = ps->getPara<TString>("dijetMC_step2output_name");
		file = jtc_utility::histNameScheme("_", 0, 1)+".root";
		fdMC_rg = TFile::Open(path+cap+file);
		file = jtc_utility::histNameScheme("_", 1, 1)+".root";
		fdMC_gg = TFile::Open(path+cap+file);
		file = jtc_utility::histNameScheme("_", 0, 0)+".root";
		fdMC_rr = TFile::Open(path+cap+file);
}

void bjtcAnalyzer_Step3::recoJetCheck(){
		gQA->pad_map = ps->getPara<mapper_func>("pad_map");
		matrixTH1Ptr m2in_rg("signal_inclJet_RecoJet_GenTrack_noCorr", npt, ncent);
		matrixTH1Ptr m2in_gg("signal_inclJet_GenJet_GenTrack_noCorr", npt, ncent);
		m2in_rg.autoLoad(fdMC_rg);
		m2in_gg.autoLoad(fdMC_gg);
		auto m2side_rg = getPhiProjection(m2in_rg, -1, 1);
		auto m2side_gg = getPhiProjection(m2in_gg, -1, 1);
		gQA->addm2TH1Pair(m2side_rg, m2side_gg);
		gQA->bookLegend();
		gQA->addLegendPair("Reco-Gen", "Gen-Gen", 0);
		gQA->setTitlePosition(0.625,0.85);
		gQA->addRLine(1);
		gQA->addhLine(0);
		gQA->setXrange(-2, 1.99);
		gQA->makeTitle = 1;
		gQA->setLowerYrange(0.6, 1.4);
		gAnaIO.saveCanvas(gQA->overlayR(), "inclJetReco_signal");
		delete m2side_rg;
		delete m2side_gg;

		matrixTH1Ptr *m2in_dr_rg= new matrixTH1Ptr("dr_raw_inclJet_RecoJet_GenTrack_noCorr", npt, ncent);
		matrixTH1Ptr *m2in_dr_gg= new matrixTH1Ptr("dr_raw_inclJet_GenJet_GenTrack_noCorr", npt, ncent);
		m2in_dr_rg->autoLoad(fdMC_rg);
		m2in_dr_gg->autoLoad(fdMC_gg);
		auto m2 = (*m2in_dr_rg)/(*m2in_dr_gg);
		gQA->addm2TH1Pair(m2in_dr_rg, m2in_dr_gg);
		gQA->bookLegend();
		gQA->addLegendPair("Reco-Gen", "Gen-Gen", 0);
		gQA->setXrange(0, 0.99);
		gQA->setLowerYrange(0.9, 1.1);
		gAnaIO.saveCanvas(gQA->overlayR(), "inclJetReco_effect_rawSig_ratio");

		matrixTH1Ptr *m2raw_dr_rg= new matrixTH1Ptr("dr_signal_inclJet_RecoJet_GenTrack_noCorr", npt, ncent);
		matrixTH1Ptr *m2raw_dr_gg= new matrixTH1Ptr("dr_signal_inclJet_GenJet_GenTrack_noCorr", npt, ncent);
		m2raw_dr_rg->autoLoad(fdMC_rg);
		m2raw_dr_gg->autoLoad(fdMC_gg);
		gQA->addm2TH1Pair(m2raw_dr_rg, m2raw_dr_gg);
		gQA->bookLegend();
		gQA->addLegendPair("Reco-Gen", "Gen-Gen", 0);
		gQA->setXrange(0, 0.99);
		gQA->setLowerYrange(0.9, 1.1);
		gAnaIO.saveCanvas(gQA->overlayR(), "inclJetReco_effect_signal_ratio");
}

void bjtcAnalyzer_Step3::recoTrackCheck(){
		gQA->pad_map = ps->getPara<mapper_func>("pad_map");
		matrixTH1Ptr* m2tt_rr = new matrixTH1Ptr("dr_raw_taggedTrueBJet_RecoJet_RecoTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2tt_rg = new matrixTH1Ptr("dr_raw_taggedTrueBJet_RecoJet_GenTrack_noCorr", npt, ncent);
		m2tt_rr->autoLoad(fbMC_rr);
		m2tt_rg->autoLoad(fbMC_rg);
		gQA->bookLegend();
		gQA->addm2TH1Pair(m2tt_rr, m2tt_rg);
		gQA->addLegendPair("Reco-Reco", "Reco-Gen", 0);
		gQA->setTitlePosition(0.625,0.85);
		gQA->addRLine(1);
		gQA->addhLine(0);
		gQA->setXrange(0, 0.99);
		gQA->setLowerYrange(0.6, 1.05);
		gQA->makeTitle = 1;
//		gQA->overlayR();
		gAnaIO.saveCanvas(gQA->overlayR(), "taggedTrueBJetRatio_trackRecoCorr");
}

void bjtcAnalyzer_Step3::bjtc_jff(){
		gQA->pad_map = ps->getPara<mapper_func>("pad_map");
		matrixTH1Ptr* m2tt_rg = new matrixTH1Ptr("dr_raw_taggedTrueBJet_RecoJet_GenTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2tb_gg = new matrixTH1Ptr("dr_raw_trueBJet_GenJet_GenTrack_noCorr", npt, ncent);
		m2tt_rg->autoLoad(fbMC_rg);
		m2tb_gg->autoLoad(fbMC_gg);
		gQA->bookLegend();
		gQA->addm2TH1Pair(m2tt_rg, m2tb_gg);
		gQA->addLegendPair("Reco-Gen", "Gen-Gen", 0);
		gQA->setTitlePosition(0.625,0.85);
		gQA->addRLine(1);
		gQA->addhLine(0);
		gQA->setXrange(0, 0.99);
		gQA->setLowerYrange(0.8, 1.2);
		gQA->makeTitle = 1;
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_pp_jffCorr");
}
