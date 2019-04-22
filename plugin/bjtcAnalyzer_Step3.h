
#include "edmJtcAnalyzer.h" 
#include "jtcTH1Player.h"

class bjtcAnalyzer_Step3 : public rootEDMAnalyzer {
		public : bjtcAnalyzer_Step3(ParaSet & pset) : rootEDMAnalyzer(){
						 _nodeClassName_ = "bjtcAnalyzer_Step3";
						 ps = &pset;
						 loadStep2OutputFiles();
						 gAnaIO.output_plot_path = ps->getPara<TString>("step3output_folder");
						 ndrbin = ps->getPara<int>("ndr");
						 drbins = ps->getVectorAsArray<float>("drbins");
				 }
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
				 void testRun();
				 void JER();
				 void get_res_corr();
				 void check_dijtc_result();
				 void recoJetCheck();
				 void forBjet_recoTrackCheck();
				 void fordijet_recoTrackCheck();
				 void check_bjtc_jff();
				 void check_dijtc_jff();
				 void get_dijtc_correction();
				 int analyze(){return 0;}
				 void check_dijtc_WTAvsEsch_result();
				 void check_WTAvsEsch_result();
				 void produce_bjtc();
				 void produce_djtc();
				 void get_bjtc_correction();
				 void data_overlay_injtc_vs_bjtc();
				 void bjtc_validation();
				 void AN_plot();
				 void ending_plots();
		public: 
				 int npt= 6, ncent=1;
				 int ndrbin;
				 float *drbins;
				 ParaSet *ps;
				 TString inputPath;
				 TFile *fbMC_gg, *fbMC_rg, *fbMC_rr, *fdMC_gg, *fdMC_rg, *fdMC_rr, *fcorr, *fdata;
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
		file = ps->getPara<TString>("pp5TeVData_step2output_name");
		fdata = TFile::Open(path+file+".root");
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

void bjtcAnalyzer_Step3::forBjet_recoTrackCheck(){
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

void bjtcAnalyzer_Step3::fordijet_recoTrackCheck(){
		gQA->pad_map = ps->getPara<mapper_func>("pad_map");
		matrixTH1Ptr* m2tt_rr = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_RecoTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2tt_rg = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_GenTrack_noCorr", npt, ncent);
		m2tt_rr->autoLoad(fdMC_rr);
		m2tt_rg->autoLoad(fdMC_rg);
		gQA->bookLegend();
		gQA->addm2TH1Pair(m2tt_rr, m2tt_rg);
		gQA->addLegendPair("Reco-Reco", "Reco-Gen", 0);
		gQA->setTitlePosition(0.625,0.85);
		gQA->addRLine(1);
		gQA->addhLine(0);
		gQA->setXrange(0, 0.99);
		gQA->setLowerYrange(0.6, 1.05);
		gQA->makeTitle = 1;
		gAnaIO.saveCanvas(gQA->overlayR(), "inclJetRatio_trackRecoCorr");
}

void bjtcAnalyzer_Step3::check_dijtc_jff(){
		gQA->pad_map = ps->getPara<mapper_func>("pad_map");
		matrixTH1Ptr* m2tt_rg = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_GenTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2tb_gg = new matrixTH1Ptr("dr_raw_inclJet_GenJet_GenTrack_noCorr", npt, ncent);
		m2tt_rg->autoLoad(fdMC_rg);
		m2tb_gg->autoLoad(fdMC_gg);
		gQA->bookLegend();
		gQA->addm2TH1Pair(m2tt_rg, m2tb_gg);
		gQA->addLegendPair("Reco-Gen", "Gen-Gen", 0);
		gQA->setTitlePosition(0.625,0.85);
		gQA->addRLine(1);
		gQA->addhLine(0);
		gQA->setXrange(0, 0.99);
		gQA->setLowerYrange(0.8, 1.2);
		gQA->makeTitle = 1;
		gAnaIO.saveCanvas(gQA->overlayR(), "dijtc_pp_jffCorr");
}

void bjtcAnalyzer_Step3::check_bjtc_jff(){
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

void bjtcAnalyzer_Step3::get_dijtc_correction(){
		gQA->pad_map = ps->getPara<mapper_func>("pad_map");
		matrixTH1Ptr* m2y_rg = new matrixTH1Ptr("dr_signal_inclJet_RecoJet_GenTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2y_gg = new matrixTH1Ptr("dr_signal_inclJet_GenJet_GenTrack_noCorr", npt, ncent);
		//matrixTH1Ptr* m2y_rg = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_GenTrack_noCorr", npt, ncent);
		//matrixTH1Ptr* m2y_gg = new matrixTH1Ptr("dr_raw_inclJet_GenJet_GenTrack_noCorr", npt, ncent);
		m2y_rg->autoLoad(fdMC_rg);
		m2y_gg->autoLoad(fdMC_gg);
		auto m2y_jff = (*m2y_rg)%(*m2y_gg);
		auto m2y_jff2= m2y_jff->clone("dijtc_yield_jff_noSmooth");
		m2y_jff->setName("dijtc_yield_jff");
		//matrixTH1Ptr* m2j_rg = new matrixTH1Ptr("dr_signal_inclJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		//matrixTH1Ptr* m2j_gg = new matrixTH1Ptr("dr_signal_inclJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		//matrixTH1Ptr* m2j_rg = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		//matrixTH1Ptr* m2j_gg = new matrixTH1Ptr("dr_raw_inclJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_rg = new matrixTH1Ptr("dr_mix_seagull_corrected_inclJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_gg = new matrixTH1Ptr("dr_mix_seagull_corrected_inclJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		m2j_rg->autoLoad(fdMC_rg);
		m2j_gg->autoLoad(fdMC_gg);
		auto m2j_jff = (*m2j_rg)%(*m2j_gg);
		//auto m2j_jff = (*m2j_rg)-(*m2j_gg);
		auto m2j_jff2= m2j_jff->clone("dijtc_jetShape_jff_noSmooth");
		m2j_jff->smooth();
		m2j_jff->setName("dijtc_jetShape_jff");
		gQA->bookLegend(.5, .2, .9, .4);
		gQA->addm2TH1(m2j_jff);
		gQA->addm2TH1(m2j_jff2);
		gQA->addLegendEntry("jff corr.", 0);
		gQA->addLegendEntry("original", 1);
		gQA->setXrange(0. , 2.499);
	//	gQA->setYrange(0.5, 1.2);
		gQA->addhLine(1);
		gAnaIO.saveCanvas(gQA->overlay(), "dijtc_jffCorr");
		

		//matrixTH1Ptr* m2y_rr = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_RecoTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2y_rr = new matrixTH1Ptr("dr_signal_inclJet_RecoJet_RecoTrack_noCorr", npt, ncent);
		m2y_rr->autoLoad(fdMC_rr);
		auto m2y_trk = (*m2y_rr)%(*m2y_rg);
		auto m2y_trk2 = m2y_trk->clone("dijtc_yield_trkCorr_noSmooth");
		m2y_trk->setName("dijtc_yield_trkCorr");

		//matrixTH1Ptr* m2j_rr = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_rr = new matrixTH1Ptr("dr_mix_seagull_corrected_inclJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		//matrixTH1Ptr* m2j_rr = new matrixTH1Ptr("dr_signal_inclJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		m2j_rr->autoLoad(fdMC_rr);
		auto m2j_trk = (*m2j_rr)%(*m2j_rg);
		auto m2j_trk2 = m2j_trk->clone("dijtc_jetShape_trkCorr_noSmooth");
		m2j_trk->setName("dijtc_jetShape_trkCorr");
		m2j_trk->smooth();

		gQA->bookLegend(.5, .2, .9, .5);
		gQA->addm2TH1(m2j_trk);
		gQA->addm2TH1(m2j_trk2);
		gQA->addLegendEntry("trk corr.", 0);
		gQA->addLegendEntry("original", 1);
		gQA->addhLine(1);
		gAnaIO.saveCanvas(gQA->overlay(), "dijtc_trkCorr");

		TString folder = ps->getPara<TString>("step3output_folder");
		auto wf = TFile::Open(folder+"inclJet_correction.root", "recreate");
		m2y_trk->smooth();
		m2y_jff->smooth();
		m2j_jff->write();
		m2y_jff->write();
		m2j_trk->write();
		m2y_trk->write();
		m2j_jff2->write();
		m2y_jff2->write();
		m2j_trk2->write();
		m2y_trk2->write();
		wf->Close();
		/*
		   gQA->bookLegend();
		   gQA->addm2TH1Pair(m2y_rg, m2y_gg);
		   gQA->addLegendPair("Reco-Gen", "Gen-Gen", 0);
		   gQA->setTitlePosition(0.625,0.85);
		   gQA->addRLine(1);
		   gQA->addhLine(0);
		   gQA->setXrange(0, 0.99);
		   gQA->setLowerYrange(0.8, 1.2);
		   gQA->makeTitle = 1;
		   gAnaIO.saveCanvas(gQA->overlayR(), "dijtc_pp_jffCorr");
		   */
}

void bjtcAnalyzer_Step3::produce_djtc(){
		TString folder = ps->getPara<TString>("step3output_folder");

		matrixTH1Ptr* m2j_data = new matrixTH1Ptr("sig_mix_seagull_corrected_inclJet_Data_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2y_data = new matrixTH1Ptr("sig_mix_seagull_corrected_inclJet_Data_noCorr", npt, ncent);
		m2y_data->autoLoad(fdata);
		m2j_data->autoLoad(fdata);
		std::cout<<"here"<<std::endl;
		matrixTH1Ptr* m2j_trk = new matrixTH1Ptr("dijtc_jetShape_trkCorr", npt, ncent);
		matrixTH1Ptr* m2y_trk = new matrixTH1Ptr("dijtc_yield_trkCorr", npt, ncent);
		matrixTH1Ptr* m2j_jff = new matrixTH1Ptr("dijtc_jetShape_jff", npt, ncent);
		matrixTH1Ptr* m2y_jff = new matrixTH1Ptr("dijtc_yield_jff", npt, ncent);
		folder = ps->getPara<TString>("step3output_folder");
		auto f0 = TFile::Open(folder+"inclJet_correction.root");
		std::cout<<folder+"inclJet_correction.root"<<std::endl;
		m2j_trk->autoLoad(f0);
		m2y_trk->autoLoad(f0);
		m2j_jff->autoLoad(f0);
		m2y_jff->autoLoad(f0);

		auto m2j_step1 = (jtcTH1Player*) m2j_data->clone("dijtc_jetShape_step1");
		auto m2y_step1 = (jtcTH1Player*) m2y_data->clone("dijtc_yield_step1");
		m2j_step1-> ring_corr(m2j_trk, 2.5);
		m2y_step1-> ring_corr(m2y_trk, 2.5);
		//step2: apply jff correction
		auto m2j_step2 = (jtcTH1Player*) m2j_step1->clone("dijtc_jetShape_step2");
		auto m2y_step2 = (jtcTH1Player*) m2y_step1->clone("dijtc_yield_step2");
		m2j_step2 ->ring_corr(m2j_jff, 2.5);
		m2y_step2 ->ring_corr(m2y_jff, 2.5);

		auto m2j_step2_prox = m2j_step2->getBkgError();
		auto m2y_step2_prox = m2y_step2->getBkgError();
		gQA->addm2TH1(m2j_step2_prox);
		gAnaIO.saveCanvas(gQA->drawBkgErrorCheck(), "djtc_data_js_bkgSystCheck");
		gQA->addm2TH1(m2y_step2_prox);
		gAnaIO.saveCanvas(gQA->drawBkgErrorCheck(), "djtc_data_py_bkgSystCheck");

		//step3: subtract the bkg
		auto m2j_step3 = m2j_step2->bkgSub("dijtc_jetShape_step3", 1.5, 2.5);
		auto m2y_step3 = m2y_step2->bkgSub("dijtc_yield_step3", 1.5, 2.5);
		
		auto dr_m2j_step3 = m2j_step3->drIntegral("dr_signal_dijtc_jetShape_step3", ndrbin, drbins);
		auto dr_m2y_step3 = m2y_step3->drIntegral("dr_signal_dijtc_yield_step3", ndrbin, drbins);
		auto dr_m2j_err = (jtcTH1Player*) dr_m2j_step3->clone("dr_dijtc_jetShape_systError");
		auto dr_m2y_err = (jtcTH1Player*) dr_m2y_step3->clone("dr_dijtc_yield_systError");
		dr_m2j_err->setDrError(m2j_step2);
		dr_m2y_err->setDrError(m2y_step2);
		
		folder = ps->getPara<TString>("step3output_folder");
		auto wf = TFile::Open(folder+"inclJet_data_final.root", "recreate");
		dr_m2j_step3->write();
		dr_m2y_step3->write();
		m2j_step2->write();
		m2y_step2->write();
		m2j_step3->write();
		m2y_step3->write();
		dr_m2j_err->write();
		dr_m2y_err->write();
		wf->Close();
}

void bjtcAnalyzer_Step3::check_dijtc_result(){
		gQA->pad_map = ps->getPara<mapper_func>("pad_map");
		matrixTH1Ptr* m2j_data = new matrixTH1Ptr("dr_signal_dijtc_jetShape_step2", npt, ncent);
		matrixTH1Ptr* m2y_data = new matrixTH1Ptr("dr_signal_dijtc_yield_step2", npt, ncent);
		TString folder = ps->getPara<TString>("step3output_folder");
		auto f = TFile::Open(folder+"inclJet_data_final.root");
		m2j_data->autoLoad(f);
		m2y_data->autoLoad(f);
		gQA->addm2TH1(m2j_data);
		gAnaIO.saveCanvas(gQA->overlay(), "djtc_pp_Rho_result");
		gQA->addm2TH1(m2y_data);
		gAnaIO.saveCanvas(gQA->overlay(), "djtc_pp_yield_result");
		f->Close();
}

void bjtcAnalyzer_Step3::check_dijtc_WTAvsEsch_result(){
		gQA->pad_map = ps->getPara<mapper_func>("pad_map");
		matrixTH1Ptr* m2j_data = new matrixTH1Ptr("dr_signal_dijtc_jetShape_step2", npt, ncent);
		matrixTH1Ptr* m2y_data = new matrixTH1Ptr("dr_signal_dijtc_yield_step2", npt, ncent);
		TString folder = ps->getPara<TString>("step3output_folder");
		auto f = TFile::Open(folder+"inclJet_data_final.root");
		m2j_data->autoLoad(f);
		m2y_data->autoLoad(f);

		matrixTH1Ptr* m2j_calo = new matrixTH1Ptr("inclJTC_pp_Data_pTweighted", npt, ncent);
		matrixTH1Ptr* m2y_calo = new matrixTH1Ptr("inclJTC_pp_Data", npt, ncent);
		auto fcalo = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/myCorrelation/trunk/v3/incl_pp_referenceForbJTC2.root");
		m2j_calo->autoLoad(fcalo);
		m2y_calo->autoLoad(fcalo);

		gQA->addm2TH1Pair(m2j_data, m2j_calo);
		gQA->makeTitle = 1;
		gQA->setXrange(0, 0.99);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addRLine(1);
		gQA->addLegendPair("PF WTA(rec)", "PF EScheme", 0);
		gQA->setTitlePosition(0.625,0.85);
		gAnaIO.saveCanvas(gQA->overlayR(), "djtc_pp_Rho_WTAvsEsc_result");
		gQA->addm2TH1Pair(m2y_data, m2y_calo);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair("PF WTA(rec)", "PF EScheme", 0);
		gQA->addRLine(1);
		gAnaIO.saveCanvas(gQA->overlayR(), "djtc_pp_yield_WTAvsEsc_result");
		f->Close();
		fcalo->Close();
}

void bjtcAnalyzer_Step3::get_bjtc_correction(){
		gQA->pad_map = ps->getPara<mapper_func>("pad_map");
		matrixTH1Ptr* m2y_tt_rg = new matrixTH1Ptr("dr_mix_seagull_corrected_taggedTrueBJet_RecoJet_GenTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2y_tb_gg = new matrixTH1Ptr("dr_mix_seagull_corrected_trueBJet_GenJet_GenTrack_noCorr", npt, ncent);
		m2y_tt_rg->autoLoad(fbMC_rg);
		m2y_tb_gg->autoLoad(fbMC_gg);

		gQA->addm2TH1Pair(m2y_tt_rg, m2y_tb_gg);
		gQA->makeTitle = 1;
		gQA->setXrange(0, 0.99);
		gQA->addRLine(1);
		gQA->setTitlePosition(0.625,0.85);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair("reco tag&true b-jet", "gen true b-jet", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_yield_pp_jff");

		auto m2y_jff = (*m2y_tt_rg)%(*m2y_tb_gg);
		m2y_jff->setName("bjtc_yield_jff");
		matrixTH1Ptr* m2j_tt_rg = new matrixTH1Ptr("dr_mix_seagull_corrected_taggedTrueBJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_tb_gg = new matrixTH1Ptr("dr_mix_seagull_corrected_trueBJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		m2j_tt_rg->autoLoad(fbMC_rg);
		m2j_tb_gg->autoLoad(fbMC_gg);
		auto m2j_jff = (*m2j_tt_rg)%(*m2j_tb_gg);
		m2j_jff->setName("bjtc_jetShape_jff");

		m2j_jff->smooth();
		m2y_jff->smooth();
		gQA->addm2TH1(m2j_jff);
		gQA->setTitlePosition(0.425,0.85);
		gQA->setYrange(0.6, 1.4);
		gQA->addhLine(1);
		gAnaIO.saveCanvas(gQA->overlay(), "bjtc_JS_pp_jff");

		gQA->addm2TH1Pair(m2j_tt_rg, m2j_tb_gg);
		gQA->setTitlePosition(0.625,0.85);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair("reco tag&true b-jet", "gen true b-jet", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_JS_pp_ratio_jff");

		matrixTH1Ptr* m2y_tt_rr = new matrixTH1Ptr("dr_mix_seagull_corrected_taggedTrueBJet_RecoJet_RecoTrack_noCorr", npt, ncent);
		m2y_tt_rr->autoLoad(fbMC_rr);
		//matrixTH1Ptr* m2y_tt_gr = new matrixTH1Ptr("dr_raw_taggedTrueBJet_GenJet_RecoTrack_noCorr", npt, ncent);
		//m2y_tt_gr->autoLoad(fbMC_gr);
		matrixTH1Ptr* m2y_tt_gg = new matrixTH1Ptr("dr_mix_seagull_corrected_taggedTrueBJet_GenJet_GenTrack_noCorr", npt, ncent);
		m2y_tt_gg->autoLoad(fbMC_gg);
		auto m2y_trk = (*m2y_tt_rr)%(*m2y_tt_rg);
//		auto m2y_trk_au = (*m2y_tt_gr)%(*m2y_tt_gg);
		auto m2y_trk_or = m2y_trk->clone("bjtc_yield_trkCorr_au");
		m2y_trk->setName("bjtc_yield_trkCorr");
		m2y_trk->smooth();

		// AN plots
		gQA->addm2TH1(m2y_trk);
		gQA->addm2TH1(m2y_trk_or);
		gQA->bookLegend(.5, .3, .9, .6);
		gQA->setYrange(0.5, 1.1);
		gQA->addhLine(1);
		gQA->setTitlePosition(0.425,0.85);
		gQA->addLegendEntry("trkCorr.", 0);
		gQA->addLegendEntry("Reco-jet.", 1);
		gAnaIO.saveCanvas(gQA->overlay(), "bjtc_yield_AN_pp_trkCorr");
		gQA->setTitlePosition(0.525,0.85);

		gQA->addm2TH1Pair(m2y_tt_rr, m2y_tt_rg);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair("t&t b-jet(reco trk)", "t&t b-jet(gen trk)", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_yield_pp_trkCorr");

		//matrixTH1Ptr* m2j_tt_rr = new matrixTH1Ptr("dr_raw_taggedTrueBJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_tt_rr = new matrixTH1Ptr("dr_mix_seagull_corrected_taggedTrueBJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		m2j_tt_rr->autoLoad(fbMC_rr);
		auto m2j_trk = (*m2j_tt_rr)%(*m2j_tt_rg);
		m2j_trk->setName("bjtc_jetShape_trkCorr");

		gQA->addm2TH1Pair(m2j_tt_rr, m2j_tt_rg);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair("t&t b-jet(reco trk)", "t&t b-jet(gen trk)", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_JS_pp_trkCorr");

		matrixTH1Ptr* m2y_co_gg = new matrixTH1Ptr("dr_mix_seagull_corrected_contJet_GenJet_GenTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_co_gg = new matrixTH1Ptr("dr_mix_seagull_corrected_contJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2y_in_gg = new matrixTH1Ptr("dr_mix_seagull_corrected_inclJet_GenJet_GenTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_in_gg = new matrixTH1Ptr("dr_mix_seagull_corrected_inclJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		m2y_co_gg->autoLoad(fdMC_gg);
		m2j_co_gg->autoLoad(fdMC_gg);
		m2y_in_gg->autoLoad(fdMC_gg);
		m2j_in_gg->autoLoad(fdMC_gg);

		matrixTH1Ptr* m2j_tg_gg = new matrixTH1Ptr("dr_mix_seagull_corrected_taggedBJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		m2j_tg_gg->autoLoad(fdMC_gg);

		auto contCorrJ = (*m2j_in_gg)/(*m2j_co_gg);
		auto contCorrY = (*m2y_in_gg)/(*m2y_co_gg);
		contCorrJ->setName("bjtc_jetShape_contCorr");
		contCorrY->setName("bjtc_yield_contCorr");
		contCorrJ->setAxisRange(0.101, 0.9, "X");	
		contCorrY->setAxisRange(0.101, 0.9, "X");	
		contCorrJ->smooth(1, "R");
		contCorrY->smooth(1, "R");
		gQA->addm2TH1Pair(m2j_co_gg, m2j_in_gg);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair("cont. jet", "incl. jet", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_JS_pp_contamination");
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->setYrange(0.8, 1.5);
		gQA->addm2TH1(contCorrJ);
		gQA->addm2TH1(contCorrY);
		gQA->addLegendEntry("JS", 0);
		gQA->addLegendEntry("yield", 1);

	//	auto test = (*m2j_tg_rg)/(*m2j_in_rg);
	//	test->smooth();
	//	gQA->addm2TH1(test);

		gQA->addhLine(1);
		gAnaIO.saveCanvas(gQA->overlay(), "bjtc_JS_pp_contCorr");

		TString folder = ps->getPara<TString>("step3output_folder");
		auto wf = TFile::Open(folder+"bjtc_correction.root", "recreate");
		m2j_trk->smooth();
		m2y_trk->smooth();
		m2j_jff->write();
		m2y_jff->write();
		m2j_trk->write();
		m2y_trk->write();
		contCorrJ->write();
		contCorrY->write();
		wf->Close();
}

void bjtcAnalyzer_Step3::testRun(){
		TString folder = ps->getPara<TString>("step3output_folder");
		jtcTH1Player* m2j_in_data = new jtcTH1Player("sig_mix_corrected_inclJet_Data_pTweighted_noCorr", npt, ncent);
		jtcTH1Player* m2j_in_data2= new jtcTH1Player("sig_mix_seagull_corrected_inclJet_Data_pTweighted_noCorr", npt, ncent);
		m2j_in_data->autoLoad(fdata);
		m2j_in_data2->autoLoad(fdata);
		auto m2p = m2j_in_data->projX("tes", -1, 1);
		auto m2err = m2j_in_data->getBkgError();
		auto m2err2= m2j_in_data2->getBkgError();
		auto m2dr = m2j_in_data->drIntegral("dr_integrated", ndrbin, drbins);
		gQA->addm2TH1(m2err);
		gQA->setXrange(-3, 2.99);
		gAnaIO.output_plot_path = "/Users/tabris/frameLite/output/testQA/";
	//	gAnaIO.saveCanvas(gQA->overlay(), "output_test");
		gAnaIO.saveCanvas(gQA->drawBkgErrorCheck(), "inclJet_bkgSystError");
		gQA->addm2TH1(m2err2);
		gAnaIO.saveCanvas(gQA->drawBkgErrorCheck(), "inclJet_sgFixed_bkgSystError");
}

void bjtcAnalyzer_Step3::produce_bjtc(){
		float purity = 0.66;
		TString folder = ps->getPara<TString>("step3output_folder");
		auto f0 = TFile::Open(folder+"bjtc_correction.root");
		auto m2j_trk = new jtcTH1Player("bjtc_jetShape_trkCorr", npt, ncent);
		auto m2y_trk = new jtcTH1Player("bjtc_yield_trkCorr", npt, ncent);
		auto m2j_jff = new jtcTH1Player("bjtc_jetShape_jff", npt, ncent);
		auto m2y_jff = new jtcTH1Player("bjtc_yield_jff", npt, ncent);
		auto m2j_cont = new jtcTH1Player("bjtc_jetShape_contCorr", npt, ncent);
		auto m2y_cont = new jtcTH1Player("bjtc_yield_contCorr", npt, ncent);
		m2j_cont->autoLoad(f0);
		m2y_cont->autoLoad(f0);
		m2j_trk->autoLoad(f0);
		m2y_trk->autoLoad(f0);
		m2j_jff->autoLoad(f0);
		m2y_jff->autoLoad(f0);

		auto m2j_in_data = new jtcTH1Player("sig_mix_seagull_corrected_inclJet_Data_pTweighted_noCorr", npt, ncent);
		auto m2y_in_data = new jtcTH1Player("sig_mix_seagull_corrected_inclJet_Data_noCorr", npt, ncent);
		auto m2j_tg_data = new jtcTH1Player("sig_mix_seagull_corrected_taggedBJet_Data_pTweighted_noCorr", npt, ncent);
		auto m2y_tg_data = new jtcTH1Player("sig_mix_seagull_corrected_taggedBJet_Data_noCorr", npt, ncent);
		m2y_tg_data->autoLoad(fdata);
		m2j_tg_data->autoLoad(fdata);
		m2y_in_data->autoLoad(fdata);
		m2j_in_data->autoLoad(fdata);

		auto m2j_co_data = (jtcTH1Player*) m2j_in_data->clone("m2j_co_data");
		auto m2y_co_data = (jtcTH1Player*) m2y_in_data->clone("m2y_co_data");
		m2j_co_data->ring_corr(m2j_cont, 0.4, 1);
		m2y_co_data->ring_corr(m2y_cont, 0.4, 1);

		auto dr_co_data = m2j_co_data->drIntegral("dr_contamination", ndrbin, drbins);
		gQA->addm2TH1(dr_co_data);
		gAnaIO.saveCanvas(gQA->overlay(), "contamination_template");

		m2j_co_data->scale(1-purity);
		m2y_co_data->scale(1-purity);

		auto m2j_tt_data = (*m2j_tg_data)-(*m2j_co_data);
		auto m2y_tt_data = (*m2y_tg_data)-(*m2y_co_data);

		m2j_tt_data->scale(1.0/purity);
		m2y_tt_data->scale(1.0/purity);

//		std::cout<<folder+"bjtc_correction.root"<<std::endl;


		auto m2j_step1 = (jtcTH1Player*) m2j_tt_data->clone("bjtc_jetShape_step1");
		auto m2y_step1 = (jtcTH1Player*) m2y_tt_data->clone("bjtc_yield_step1");
		m2j_step1->ring_corr(m2j_trk, 2.5);
		m2y_step1->ring_corr(m2y_trk, 2.5);
		//step2: apply jff correction
		auto m2j_step2 =(jtcTH1Player*) m2j_step1->clone("signal_bjtc_jetShape_step2");
		auto m2y_step2 =(jtcTH1Player*) m2y_step1->clone("signal_bjtc_yield_step2");
		m2j_step2->ring_corr(m2j_jff, 2.5);
		m2y_step2->ring_corr(m2y_jff, 2.5);


		auto m2j_step2_prox = m2j_step2->getBkgError();
		auto m2y_step2_prox = m2y_step2->getBkgError();
		gQA->addm2TH1(m2j_step2_prox);
		gAnaIO.saveCanvas(gQA->drawBkgErrorCheck(), "bjtc_data_js_bkgSystCheck");
		gQA->addm2TH1(m2y_step2_prox);
		gAnaIO.saveCanvas(gQA->drawBkgErrorCheck(), "bjtc_data_py_bkgSystCheck");
		//step3: subtract the bkg
		auto m2j_step3 = m2j_step2->bkgSub("signal_bjtc_jetShape_step3", 1.5, 2.5);
		auto m2y_step3 = m2y_step2->bkgSub("signal_bjtc_yield_step3", 1.5, 2.5);
		auto dr_m2j_step3 = m2j_step3->drIntegral("dr_signal_bjtc_jetShape_step3", ndrbin, drbins);
		auto dr_m2y_step3 = m2y_step3->drIntegral("dr_signal_bjtc_yield_step3", ndrbin, drbins);

		auto dr_m2j_err = (jtcTH1Player*) dr_m2j_step3->clone("dr_bjtc_jetShape_systError");
		auto dr_m2y_err = (jtcTH1Player*) dr_m2y_step3->clone("dr_bjtc_yield_systError");
		dr_m2j_err->setDrError(m2j_step2);
		dr_m2y_err->setDrError(m2y_step2);

		folder = ps->getPara<TString>("step3output_folder");
		auto wf = TFile::Open(folder+"bjtc_data_final.root", "recreate");
		m2j_step2->write();
		m2y_step2->write();
		m2j_step3->write();
		m2y_step3->write();
		dr_m2j_err->write();
		dr_m2j_step3->write();
		dr_m2y_step3->write();
		dr_m2j_err->write();
		dr_m2y_err->write();

		wf->Close();
}

void bjtcAnalyzer_Step3::data_overlay_injtc_vs_bjtc(){
		matrixTH1Ptr* m2j_data = new matrixTH1Ptr("dr_signal_dijtc_jetShape_step2", npt, ncent);
		matrixTH1Ptr* m2y_data = new matrixTH1Ptr("dr_signal_dijtc_yield_step2", npt, ncent);
		TString folder = ps->getPara<TString>("step3output_folder");
		auto f = TFile::Open(folder+"inclJet_data_final.root");
		m2j_data->autoLoad(f);
		m2y_data->autoLoad(f);

		matrixTH1Ptr* m2j_bjtc = new matrixTH1Ptr("dr_signal_bjtc_jetShape_step2", npt, ncent);
		matrixTH1Ptr* m2y_bjtc = new matrixTH1Ptr("dr_signal_bjtc_yield_step2", npt, ncent); 
		auto fb = TFile::Open(folder+"bjtc_data_final.root");
		m2j_bjtc->autoLoad(fb);
		m2y_bjtc->autoLoad(fb);

		gQA->addm2TH1Pair(m2j_bjtc, m2j_data);
		gQA->setLowerYrange(0.5 , 2);
		gQA->setXrange(0, 0.99);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addRLine(1);
		gQA->addLegendPair( "b-jet", "incl-jet", 0);
		gQA->setTitlePosition(0.625,0.85);
		gAnaIO.saveCanvas(gQA->overlayR(), "pp_Rho_WTA_incl_vs_bjtc_result");
		gQA->addm2TH1Pair(m2y_bjtc, m2y_data);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair( "b-jet", "incl-jet", 0);
		gQA->addRLine(1);
		gAnaIO.saveCanvas(gQA->overlayR(), "pp_yield_WTA_incl_vs_bjtc_result");
		f->Close();
		fb->Close();
}

void bjtcAnalyzer_Step3::bjtc_validation(){
		float purity = 0.7;
		gQA->setXrange(0, 0.99);
		gQA->addRLine(1);
		gQA->setTitlePosition(0.625,0.85);
		auto m2j_tg_rr = new jtcTH1Player("sig_mix_seagull_corrected_taggedBJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		auto m2j_tt_rr = new jtcTH1Player("sig_mix_seagull_corrected_taggedTrueBJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		auto m2j_in_rr = new jtcTH1Player("sig_mix_seagull_corrected_inclJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		auto m2j_coCorr = new jtcTH1Player("bjtc_jetShape_contCorr", npt, ncent);
		auto folder = ps->getPara<TString>("step3output_folder");
		auto f0 = TFile::Open(folder+"bjtc_correction.root");
		m2j_coCorr->autoLoad(f0);
		//matrixTH1Ptr* m2j_tg_rr = new matrixTH1Ptr("sig_mix_seagull_corrected_taggedBJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		//matrixTH1Ptr* m2j_tt_rr = new matrixTH1Ptr("sig_mix_seagull_corrected_taggedTrueBJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		//matrixTH1Ptr* m2j_in_rr = new matrixTH1Ptr("sig_mix_seagull_corrected_inclJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		m2j_tg_rr->autoLoad(fdMC_rg);
		m2j_tt_rr->autoLoad(fdMC_rg);
		m2j_in_rr->autoLoad(fdMC_rg);

		m2j_in_rr->ring_corr(m2j_coCorr, 0.4, 1);
		m2j_in_rr->scale(1-purity);
		auto m2j_pu_rr = (*m2j_tg_rr)-(*m2j_in_rr);
		m2j_pu_rr->scale(1.0/purity);
		auto m2j_pu_rr_dr = ((jtcTH1Player*)m2j_pu_rr)->drIntegral("purified_rr", ndrbin,  drbins);
		auto m2j_tg_rr_dr = m2j_tg_rr->drIntegral("tg_rr", ndrbin,  drbins);
		auto m2j_tt_rr_dr = m2j_tt_rr->drIntegral("tt_rr", ndrbin,  drbins);
		gQA->addm2TH1Pair(m2j_tg_rr_dr, m2j_tt_rr_dr);
//		gQA->addm2TH1Pair(m2j_tg_rr, m2j_tt_rr);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair( "tagged b-jet", "t&t b-jet", 0);
		gQA->setLowerYrange(0.8 , 1.2);
		gAnaIO.saveCanvas(gQA->overlayR(), "validation_before_purification");
		gQA->addm2TH1Pair(m2j_pu_rr_dr, m2j_tt_rr_dr);
//		gQA->addm2TH1Pair(m2j_pu_rr, m2j_tt_rr);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->box.SetFillColorAlpha(kGray+2, 0.5);
		gQA->addBox_ratioPad(0, .95, 1,  1.05);
		gQA->addLegendPair( "purified b-jet", "t&t b-jet", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "validation_purification_withContCorr");

		auto m2j_in_gg_sig = new jtcTH1Player("signal_inclJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		auto m2j_tb_gg_sig = new jtcTH1Player("signal_trueBJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		auto m2j_in_gg = new jtcTH1Player("dr_signal_inclJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		auto m2j_tb_gg = new jtcTH1Player("dr_signal_trueBJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		m2j_in_gg->autoLoad(fdMC_gg);
		m2j_tb_gg->autoLoad(fbMC_gg);
		m2j_in_gg_sig->autoLoad(fdMC_gg);
		m2j_tb_gg_sig->autoLoad(fbMC_gg);
		auto m2j_in_gg_x = m2j_in_gg_sig->getBkgError();
		auto m2j_tb_gg_x = m2j_tb_gg_sig->getBkgError();

		auto m2j_in_gg_err = (jtcTH1Player*) m2j_in_gg->clone("m2j_in_gg_err");
		auto m2j_tb_gg_err = (jtcTH1Player*) m2j_tb_gg->clone("m2j_tb_gg_err");
		m2j_in_gg_err->loadBkgError(m2j_in_gg_sig);
		m2j_tb_gg_err->loadBkgError(m2j_tb_gg_sig);
		m2j_in_gg_err->setDrError();
		m2j_tb_gg_err->setDrError();
		gQA->autoYrange = 1;
		gQA->addm2TH1Pair(m2j_tb_gg, m2j_in_gg);
		gQA->addm2TH1ErrorPair(m2j_tb_gg_err, m2j_in_gg_err);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair( "incl b-jet (GenGen)", "t&t b-jet (GenGen)", 0);
		gQA->setLowerYrange(0.5 , 2);
		gAnaIO.saveCanvas(gQA->overlayR(), "validation_GenGen_inclOverTrueB");
}

void bjtcAnalyzer_Step3::check_WTAvsEsch_result(){
		TString folder = ps->getPara<TString>("step3output_folder");
		auto fwta_in = TFile::Open(folder+"inclJet_data_final.root");
		auto fesm_in = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/myCorrelation/trunk/v3/incl_pp_referenceForbJTC2.root");
		auto fwta_b = TFile::Open(folder+"bjtc_data_final.root");
		auto fesm_b = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/myCorrelation/trunk/v3/pp_data_csv_v2_purity66_v3sequence.root");
		matrixTH1Ptr* m2j_esm_in_data = new matrixTH1Ptr("inclJTC_pp_Data_pTweighted", npt, ncent);
		matrixTH1Ptr* m2j_esm_b_data = new matrixTH1Ptr("dr_pt_final_data", npt, ncent);
		matrixTH1Ptr* m2j_wta_in_data = new matrixTH1Ptr("dr_signal_dijtc_jetShape_step2", npt, ncent);
		matrixTH1Ptr* m2j_wta_b_data = new matrixTH1Ptr("dr_signal_bjtc_jetShape_step2", npt, ncent);
		m2j_wta_in_data->autoLoad(fwta_in);
		m2j_wta_b_data->autoLoad(fwta_b);
		m2j_esm_in_data->autoLoad(fesm_in);
		m2j_esm_b_data->autoLoad(fesm_b);
		auto m2j_esm_in_data_rebin = m2j_esm_in_data->rebin("rebined_esm_in", ndrbin, drbins);

		auto ratio_esm = (*m2j_esm_b_data)/(*m2j_esm_in_data);
		auto ratio_wta = (*m2j_wta_b_data)/(*m2j_wta_in_data);
		gQA->addm2TH1(ratio_esm);
		gQA->addm2TH1(ratio_wta);
		gQA->makeTitle = 1;
		gQA->setXrange(0, 0.99);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addRLine(1);
		gQA->addLegendEntry("EScheme", 0);
		gQA->addLegendEntry("WTA", 0);
		gAnaIO.saveCanvas(gQA->overlay(), "bjtc_wta_vs_esm");
}

void bjtcAnalyzer_Step3::JER(){
		TString folder = ps->getPara<TString>("step2output_folder");
		auto fd_rg_gA = TFile::Open(folder+"Signal_PYTHIA6_dijetSample_allJets_genAxis_RecoJet_GenTrack.root");
		auto fb_rg_gA = TFile::Open(folder+"Signal_PYTHIA6_bjetSample_allJets_genAxis_RecoJet_GenTrack.root");
		auto djtc_rg_ga = new jtcTH1Player("dr_signal_inclJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		auto djtc_gg = new jtcTH1Player("dr_signal_inclJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		auto djtc_rg = new jtcTH1Player("dr_signal_inclJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		djtc_rg_ga->autoLoad(fd_rg_gA);
		djtc_rg->autoLoad(fdMC_rg);
		djtc_gg->autoLoad(fdMC_gg);
		gQA->addm2TH1Pair(djtc_rg, djtc_rg_ga);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addRLine(1);
		gQA->addLegendPair( "reco-jet-axis", "gen-jet-axis", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "djtc_jer_syst");

		gQA->addm2TH1Pair(djtc_rg_ga, djtc_gg);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair( "reco-jet (gen-axis)", "gen-jet", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "djtc_jer_2_syst");

		auto bjtc_rg_ga = new jtcTH1Player("dr_signal_trueBJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		auto bjtc_gg = new jtcTH1Player("dr_signal_trueBJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		auto bjtc_rg = new jtcTH1Player("dr_signal_trueBJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		bjtc_rg_ga->autoLoad(fb_rg_gA);
		bjtc_rg->autoLoad(fbMC_rg);
		bjtc_gg->autoLoad(fbMC_gg);

		gQA->addm2TH1Pair(bjtc_rg, bjtc_rg_ga);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair( "reco-jet-axis", "gen-jet-axis", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_jer_syst");

		gQA->addm2TH1Pair(bjtc_rg_ga, bjtc_gg);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair( "reco-jet (gen-axis)", "gen-jet", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_jer_2_syst");
}


void bjtcAnalyzer_Step3::ending_plots(){
		int ndr = 21;
		float drs[] = {0.,.025, 0.05, 0.075, .1, .125, .15, .175, .2, 0.25,.3,0.35,0.4,0.45,.5, 0.6,0.8,1.,1.2, 1.5, 2, 2.5};
		TString folder = ps->getPara<TString>("step3output_folder");
		auto fd = TFile::Open(folder+"inclJet_data_final.root");
		auto fb = TFile::Open(folder+"bjtc_data_final.root");
		//auto sigj_di_data = new jtcTH1Player("signal_dijtc_jetShape_step2", npt, ncent);
		//auto sigj_b_data = new jtcTH1Player("signal_bjtc_jetShape_step2", npt, ncent);
		auto h2data_d= new jtcTH1Player("dijtc_jetShape_step3", npt, ncent);
		auto data_di = new jtcTH1Player("dr_signal_dijtc_jetShape_step3");
		auto data_b = new jtcTH1Player("dr_signal_bjtc_jetShape_step3");
		auto h2data_b= new jtcTH1Player("signal_bjtc_jetShape_step3", npt, ncent);
		auto data_di_err= new jtcTH1Player("dr_dijtc_jetShape_systError", npt, ncent);
		auto data_b_err= new jtcTH1Player("dr_bjtc_jetShape_systError", npt, ncent);
		h2data_d->autoLoad(fd);
		h2data_b->autoLoad(fb);
		data_di->autoLoad(fd);
		data_b->autoLoad(fb);
		data_b_err->autoLoad(fb);
		data_di_err->autoLoad(fd);
		data_b_err->add_frac_error(0.11);
		data_di_err->add_frac_error(0.07);
		auto dr_data_d = h2data_d->drIntegral("dr_djtc_js", ndr, drs);
		auto dr_data_b = h2data_b->drIntegral("dr_bjtc_js", ndr, drs);
		dr_data_d->doGeoCorr(h2data_d);
		dr_data_b->doGeoCorr(h2data_b);
		auto xdata_b = h2data_b->projX("prox_data_bjtc", -1, 1, "e", 0);
		auto xdata_d = h2data_d->projX("prox_data_incl", -1, 1, "e", 0);
		xdata_b->rebin();
		xdata_d->rebin();
		gQA->bookLegend(.7, .6, .95, .8);
		gQA->addm2TH1Pair(xdata_b, xdata_d);
		gQA->setTitlePosition(0.6, 0.85);
		gQA->addLegendPair("b-jet", "incl.", 0);
		gQA->addRLine(1);
		gQA->addhLine(0);
		gQA->setXrange(-1, .99);
		gAnaIO.saveCanvas(gQA->overlayR(), "overlay_ratio_data_deta");
		gQA->addm2TH1Pair(dr_data_b, dr_data_d);
		gQA->bookLegend(.7, .6, .95, .8);
		gQA->setXrange(0, .99);
		gQA->addLegendPair("b-jet", "incl.", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "overlay_ratio_data_dr");
		//sigj_di_data->autoLoad(fd);
		//sigj_b_data->autoLoad(fb);

		auto mc_di_sig= new jtcTH1Player("signal_inclJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		auto mc_b_sig= new jtcTH1Player("signal_trueBJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		mc_di_sig->autoLoad(fdMC_gg);
		mc_b_sig->autoLoad(fbMC_gg);
		auto mc_di_prox = mc_di_sig->getBkgError();
		auto mc_b_prox = mc_b_sig->getBkgError();
		auto dr_mc_di = mc_di_sig->drIntegral("dr_dijtc_jetShape_gg",ndrbin, drbins);
		auto x_mc_di = mc_di_sig->projX("prox_gg_incl", -1, 1, "e", 0);
		auto x_mc_b  = mc_b_sig->projX("prox_gg_b", -1, 1, "e", 0);
		auto dr_mc_b = mc_b_sig->drIntegral("dr_bjtc_jetShape_gg", ndrbin, drbins);
		auto mc_di_err= (jtcTH1Player*)dr_mc_di->clone("dr_dijtc_jetShape_systError");
		auto mc_b_err = (jtcTH1Player*)dr_mc_b->clone("dr_bjtc_jetShape_systError");
		mc_di_err->setDrError(mc_di_sig);
		mc_b_err->setDrError(mc_b_sig);

		auto m2_mc_b0 = new jtcTH1Player("signal_trueBJet_GenJet_GenTrack_pTweighted_noCorr");
		auto dr_mc_b0 = new jtcTH1Player("dr_signal_trueBJet_GenJet_GenTrack_pTweighted_noCorr");
		auto fmc0 = TFile::Open("/Users/tabris/frameLite/output/step2/Signal_PYTHIA6_bjetSample_noGPS_GenJet_GenTrack.root_GenJet_GenTrack.root");
		dr_mc_b0->autoLoad(fmc0);
		m2_mc_b0->autoLoad(fmc0);
		m2_mc_b0->getBkgError();
		auto dr_mc_b0_err = (jtcTH1Player*)dr_mc_b0->clone("dr_gg_mc0_err");
		dr_mc_b0_err->setDrError(m2_mc_b0);

		auto ratio_data = (*data_b)/(*data_di);
		auto ratio_data_err = (*data_b_err)/(*data_di_err);
		auto ratio_mc = (*dr_mc_b)/(*dr_mc_di);
		auto ratio_mc_err = (*mc_b_err)/(*mc_di_err);
		auto ratio_mc0 = (*dr_mc_b0)/(*dr_mc_di);
		auto ratio_mc0_err = (*dr_mc_b0_err)/(*mc_di_err);

		gQA->addm2TH1(ratio_data);
		gQA->addm2TH1(ratio_mc);
		gQA->addm2TH1(ratio_mc0);
		gQA->bookLegend(.3, .6, .9, .8);
		gQA->setTitlePosition(0.325,0.85);
		gQA->addLegendEntry("b-jet/incl-jet (data)", 0);
		gQA->addLegendEntry("b-jet/incl-jet (GSP weighted)", 1);
		gQA->addLegendEntry("b-jet/incl-jet (PYTHIA6)", 2);
		gQA->addm2TH1Error(ratio_data_err);
		gQA->addm2TH1Error(ratio_mc_err);
		gQA->addm2TH1Error(ratio_mc0_err);
		gQA->setXrange(0, 0.99);
		gQA->addRLine(1);
		gQA->setYrange(0.5, 3);
		gAnaIO.saveCanvas(gQA->overlay(), "overlay_ratio_data_vs_mc");

		auto dr_data_b_all = data_b->contractX("dr_data_bjtc_js_allpt");
		auto dr_data_d_all = data_di->contractX("dr_data_incl_js_allpt");
		auto dr_data_b_err = data_b_err->contractX("dr_data_bjtc_allpt_syt_err");
		auto dr_data_d_err = data_di_err->contractX("dr_data_incl_allpt_syst_err");
		auto dr_mc_d_all= dr_mc_di->contractX("dr_mc_incl_js_allpt");
		auto dr_mc_b_all= dr_mc_b->contractX("dr_mc_bjtc_js_allpt");
		auto dr_mc_b_err= mc_b_err->contractX("dr_mc_bjtc_allpt_syst_err");
		auto dr_mc_d_err= mc_di_err->contractX("dr_mc_incl_allpt_syst_err");
		auto dr_mc0_b_all= dr_mc_b0->contractX("dr_mc_gsp_unweighted_bjtc_syst_allpt");
		auto dr_mc0_b_err= dr_mc_b0_err->contractX("dr_mc_gsp_unweighted_bjtc_allpt_syst_err");

		auto ratio_data_all = (*dr_data_b_all)/(*dr_data_d_all);
		auto ratio_data_all_err = (*dr_data_b_err)/(*dr_data_d_err);
		auto ratio_mc_all = (*dr_mc_b_all)/(*dr_mc_d_all);
		auto ratio_mc_all_err = (*dr_mc_b_err)/(*dr_mc_d_err);
		auto ratio_mc0_all = (*dr_mc0_b_all)/(*dr_mc_d_all);
		auto ratio_mc0_all_err = (*dr_mc0_b_err)/(*dr_mc_d_err);

		gQA->addm2TH1(ratio_data_all);
		gQA->addm2TH1Error(ratio_data_all_err);
		gQA->addm2TH1(ratio_mc_all);
		gQA->addm2TH1Error(ratio_mc_all_err);
//		gQA->addm2TH1(ratio_mc_all_err, "ce5");
		gQA->addm2TH1(ratio_mc0_all);
		gQA->addm2TH1Error(ratio_mc0_all_err);
		gQA->bookLegend(0.55, 0.55, 0.9, 0.83);
		gQA->addLegendEntry("Data", 0);
		gQA->addLegendEntry("PYTHIA6(GSP weighted)", 1);
		gQA->addLegendEntry("PYTHIA6", 2);
		gQA->ncol = 1;
		gQA->nrow = 1;
		gQA->setXrange(.0, .99);
		gQA->setLowerYrange(.5, 2.);
		bjtc_pp_config::ptString[0]= "p^{track}_{T} > 1 GeV";
		gQA->Ytitle = "R(#Delta r)";
		gAnaIO.saveCanvas(gQA->overlay(), "integratedJSRatio");


		gQA->autoYrange = 1;
		gQA->fixYrange = 0;
		gQA->addm2TH1(dr_data_b_all);
		gQA->addm2TH1Error(dr_data_b_err);
		gQA->addm2TH1(dr_data_d_all);
		gQA->addm2TH1Error(dr_data_d_err);
		gQA->bookLegend(0.6, 0.6, 0.9, 0.8);
		gQA->addLegendEntry("b-jet", 0);
		gQA->addLegendEntry("incl. jet", 1);
		gAnaIO.saveCanvas(gQA->overlay(), "integratedJS");
	
		auto ym2_mc_d= new jtcTH1Player("signal_inclJet_GenJet_GenTrack_noCorr", npt, ncent);
		auto ym2_mc_b= new jtcTH1Player("signal_trueBJet_GenJet_GenTrack_noCorr", npt, ncent);
		ym2_mc_b->autoLoad(fbMC_gg);
		ym2_mc_d->autoLoad(fdMC_gg);
		ym2_mc_b->getBkgError();
		ym2_mc_d->getBkgError();
		auto ymc_d = ym2_mc_b->drIntegral("dr_signal_incl_yield_mc", ndrbin, drbins);
		auto ymc_b = ym2_mc_b->drIntegral("dr_signal_bjtc_yield_mc", ndrbin, drbins);
		auto ymc_d_err = (jtcTH1Player*) ymc_d->clone("dr_signal_incl_yield_mc_systError");
		auto ymc_b_err = (jtcTH1Player*) ymc_d->clone("dr_signal_bjtc_yield_mc_systError");
		ymc_d_err->setDrError(ym2_mc_d);
		ymc_b_err->setDrError(ym2_mc_b);

		auto ym2_mc0_d= new jtcTH1Player("signal_inclJet_GenJet_GenTrack_noCorr", npt, ncent);
		auto ym2_mc0_b= new jtcTH1Player("signal_trueBJet_GenJet_GenTrack_noCorr", npt, ncent);
		ym2_mc0_b->autoLoad(fmc0);
		ym2_mc0_d->autoLoad(fmc0);
		ym2_mc0_b->getBkgError();
		ym2_mc0_d->getBkgError();
		auto ymc0_d = ym2_mc0_b->drIntegral("dr_signal_incl_yield_mc0", ndrbin, drbins);
		auto ymc0_b = ym2_mc0_b->drIntegral("dr_signal_bjtc_yield_mc0", ndrbin, drbins);
		auto ymc0_d_err = (jtcTH1Player*) ymc0_d->clone("dr_signal_incl_yield_mc0_systError");
		auto ymc0_b_err = (jtcTH1Player*) ymc0_d->clone("dr_signal_bjtc_yield_mc0_systError");
		ymc0_d_err->setDrError(ym2_mc0_d);
		ymc0_b_err->setDrError(ym2_mc0_b);

		auto ydata_d = new jtcTH1Player("dr_signal_dijtc_yield_step3");
		auto ydata_d_err = new jtcTH1Player("dr_dijtc_yield_systError");
		auto ydata_b = new jtcTH1Player("dr_signal_bjtc_yield_step3");
		auto ydata_b_err = new jtcTH1Player("dr_bjtc_yield_systError");
		ydata_d->autoLoad(fd);
		ydata_d_err->autoLoad(fd);
		ydata_b->autoLoad(fb);
		ydata_b_err->autoLoad(fb);
		ydata_d_err->add_frac_error(0.07);
		ydata_b_err->add_frac_error(0.11);

		auto fw = TFile::Open("AN17337Result.root", "recreate");
		data_di->write();
		data_b->write();
		data_b_err->write();
		data_di_err->write();
		dr_data_b_all->write();
		dr_data_d_all->write();
		dr_data_b_err->write();
		dr_data_d_err->write();
		ydata_d->write();
		ydata_b->write();
		ydata_d_err->write();
		ydata_b_err->write();
		ymc_b->write();
		ymc_d->write();
		ymc_b_err->write();
		ymc_d_err->write();
		ymc0_b->write();
		ymc0_d->write();
		ymc0_b_err->write();
		ymc0_d_err->write();
		//dr_mc_d_all->write();
		//dr_mc_b_all->write();
		//dr_mc_b_err->write();
		//dr_mc_d_err->write();
		ratio_data_all    ->write();
		ratio_data_all_err->write();
		ratio_mc_all      ->write();
		ratio_mc_all_err  ->write();
		ratio_mc0_all     ->write();
		ratio_mc0_all_err ->write();

		/*
		gQA->fixYrange = 0;
		gQA->autoYrange = 1;
		gQA->addm2TH1Pair(data_di, dr_mc_di);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair( "inclJet Data", "inclJet MC", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "overlay_inclJet_data_vs_mc");
		gQA->addm2TH1Pair(data_b, dr_mc_b);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair( "b-jet Data", "b-jet MC", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "overlay_bJet_data_vs_mc");
		*/
}

