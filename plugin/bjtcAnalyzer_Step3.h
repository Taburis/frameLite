
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
				 void get_res_corr();
				 void check_dijtc_result();
				 void recoJetCheck();
				 void forBjet_recoTrackCheck();
				 void fordijet_recoTrackCheck();
				 void check_bjtc_jff();
				 void check_dijtc_jff();
				 void get_dijtc_correction();
				 int analyze(){return 0;}
				 void apply_dijtc_correction();
				 void check_dijtc_WTAvsEsch_result();
				 void check_WTAvsEsch_result();
				 void produce_bjtc();
				 void get_bjtc_correction();
				 void overlay_injtc_vs_bjtc();
				 void bjtc_validation();
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
		matrixTH1Ptr* m2y_rg = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_GenTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2y_gg = new matrixTH1Ptr("dr_raw_inclJet_GenJet_GenTrack_noCorr", npt, ncent);
		m2y_rg->autoLoad(fdMC_rg);
		m2y_gg->autoLoad(fdMC_gg);
		auto m2y_jff = (*m2y_rg)%(*m2y_gg);
		m2y_jff->setName("dijtc_yield_jff");
		matrixTH1Ptr* m2j_rg = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_gg = new matrixTH1Ptr("dr_raw_inclJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		m2j_rg->autoLoad(fdMC_rg);
		m2j_gg->autoLoad(fdMC_gg);
		auto m2j_jff = (*m2j_rg)%(*m2j_gg);
		m2j_jff->setName("dijtc_jetShape_jff");

		matrixTH1Ptr* m2y_rr = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_RecoTrack_noCorr", npt, ncent);
		m2y_rr->autoLoad(fdMC_rr);
		auto m2y_trk = (*m2y_rr)%(*m2y_rg);
		m2y_trk->setName("dijtc_yield_trkCorr");

		matrixTH1Ptr* m2j_rr = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		m2j_rr->autoLoad(fdMC_rr);
		auto m2j_trk = (*m2j_rr)%(*m2j_rg);
		m2j_trk->setName("dijtc_jetShape_trkCorr");

		TString folder = ps->getPara<TString>("step3output_folder");
		auto wf = TFile::Open(folder+"inclJet_correction.root", "recreate");
		m2j_jff->write();
		m2y_jff->write();
		m2j_trk->write();
		m2y_trk->write();
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

void bjtcAnalyzer_Step3::apply_dijtc_correction(){
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

		//step1: apply track correction
		edmJtcAnalyzer p2yield(*ps), p2shape(*ps);
		p2yield.initialization();
		p2shape.initialization();

		auto m2j_step1 = p2shape.m2_ring_corr("dijtc_jetShape_step1", m2j_data, m2j_trk);
		auto m2y_step1 = p2yield.m2_ring_corr("dijtc_yield_step1", m2y_data, m2y_trk);
		//step2: apply jff correction
		auto m2j_step2 = p2shape.m2_ring_corr("dijtc_jetShape_step2", m2j_step1, m2j_jff);
		auto m2y_step2 = p2yield.m2_ring_corr("dijtc_yield_step2", m2y_step1, m2y_jff);
		//step3: subtract the bkg
		p2shape.add2Step3( m2j_step2);
		p2yield.add2Step3( m2y_step2);
		p2shape.runProduce();
		p2yield.runProduce();
		folder = ps->getPara<TString>("step3output_folder");
		auto wf = TFile::Open(folder+"inclJet_data_final.root", "recreate");
		p2shape.write();
		p2yield.write();
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
		matrixTH1Ptr* m2y_tt_rg = new matrixTH1Ptr("dr_raw_taggedTrueBJet_RecoJet_GenTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2y_tb_gg = new matrixTH1Ptr("dr_raw_trueBJet_GenJet_GenTrack_noCorr", npt, ncent);
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
		matrixTH1Ptr* m2j_tt_rg = new matrixTH1Ptr("dr_raw_taggedTrueBJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_tb_gg = new matrixTH1Ptr("dr_raw_trueBJet_GenJet_GenTrack_pTweighted_noCorr", npt, ncent);
		m2j_tt_rg->autoLoad(fbMC_rg);
		m2j_tb_gg->autoLoad(fbMC_gg);
		auto m2j_jff = (*m2j_tt_rg)%(*m2j_tb_gg);
		m2j_jff->setName("bjtc_jetShape_jff");

		gQA->addm2TH1Pair(m2j_tt_rg, m2j_tb_gg);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair("reco tag&true b-jet", "gen true b-jet", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_JS_pp_jff");

		matrixTH1Ptr* m2y_tt_rr = new matrixTH1Ptr("dr_raw_taggedTrueBJet_RecoJet_RecoTrack_noCorr", npt, ncent);
		m2y_tt_rr->autoLoad(fbMC_rr);
		auto m2y_trk = (*m2y_tt_rr)%(*m2y_tt_rg);
		m2y_trk->setName("bjtc_yield_trkCorr");

		gQA->addm2TH1Pair(m2y_tt_rr, m2y_tt_rg);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair("t&t b-jet(reco trk)", "t&t b-jet(gen trk)", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_yield_pp_trkCorr");

		matrixTH1Ptr* m2j_tt_rr = new matrixTH1Ptr("dr_raw_taggedTrueBJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		m2j_tt_rr->autoLoad(fbMC_rr);
		auto m2j_trk = (*m2j_tt_rr)%(*m2j_tt_rg);
		m2j_trk->setName("bjtc_jetShape_trkCorr");

		gQA->addm2TH1Pair(m2j_tt_rr, m2j_tt_rg);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair("t&t b-jet(reco trk)", "t&t b-jet(gen trk)", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_JS_pp_trkCorr");

		matrixTH1Ptr* m2y_co_rr = new matrixTH1Ptr("dr_raw_contJet_RecoJet_RecoTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_co_rr = new matrixTH1Ptr("dr_raw_contJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2y_in_rr = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_RecoTrack_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_in_rr = new matrixTH1Ptr("dr_raw_inclJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		m2y_co_rr->autoLoad(fdMC_rr);
		m2j_co_rr->autoLoad(fdMC_rr);
		m2y_in_rr->autoLoad(fdMC_rr);
		m2j_in_rr->autoLoad(fdMC_rr);
		auto contCorrJ = (*m2j_co_rr)/(*m2j_in_rr);
		auto contCorrY = (*m2y_co_rr)/(*m2y_in_rr);
		contCorrJ->smooth();
		contCorrY->smooth();
		gQA->addm2TH1Pair(m2j_co_rr, m2j_in_rr);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair("cont. jet", "incl. jet", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "bjtc_JS_pp_contamination");
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->setYrange(0.8, 1.2);
		gQA->addm2TH1(contCorrJ);
		gQA->addm2TH1(contCorrY);
		gQA->addLegendEntry("JS", 0);
		gQA->addLegendEntry("yield", 1);
		gAnaIO.saveCanvas(gQA->overlay(), "bjtc_JS_pp_contCorr");

		TString folder = ps->getPara<TString>("step3output_folder");
		auto wf = TFile::Open(folder+"bjtc_correction.root", "recreate");
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
		jtcTH1Player* m2j_in_data = new jtcTH1Player("signal_inclJet_Data_pTweighted_noCorr", npt, ncent);
		m2j_in_data->autoLoad(fdata);
		auto m2p = m2j_in_data->projX("tes", -1, 1);
		auto m2err = m2j_in_data->getBkgError();
		auto m2dr = m2j_in_data->drIntegral("dr_integrated", ndrbin, drbins);
		gQA->addm2TH1(m2err);
		gQA->setXrange(-3, 2.99);
		gAnaIO.output_plot_path = "./";
	//	gAnaIO.saveCanvas(gQA->overlay(), "output_test");
		gAnaIO.saveCanvas(gQA->drawBkgErrorCheck(), "output_test");
}

void bjtcAnalyzer_Step3::produce_bjtc(){
		float purity = 0.7;
		TString folder = ps->getPara<TString>("step3output_folder");

		matrixTH1Ptr* m2j_in_data = new matrixTH1Ptr("sig_mix_seagull_corrected_inclJet_Data_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2y_in_data = new matrixTH1Ptr("sig_mix_seagull_corrected_inclJet_Data_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_tg_data = new matrixTH1Ptr("sig_mix_seagull_corrected_taggedBJet_Data_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2y_tg_data = new matrixTH1Ptr("sig_mix_seagull_corrected_taggedBJet_Data_noCorr", npt, ncent);
		m2y_tg_data->autoLoad(fdata);
		m2j_tg_data->autoLoad(fdata);
		m2y_in_data->autoLoad(fdata);
		m2j_in_data->autoLoad(fdata);
		m2j_in_data->scale(1-purity);
		m2y_in_data->scale(1-purity);
		auto m2j_tt_data = (*m2j_tg_data)-(*m2j_in_data);
		auto m2y_tt_data = (*m2y_tg_data)-(*m2y_in_data);

		m2j_tt_data->scale(1.0/purity);
		m2y_tt_data->scale(1.0/purity);

		matrixTH1Ptr* m2j_trk = new matrixTH1Ptr("bjtc_jetShape_trkCorr", npt, ncent);
		matrixTH1Ptr* m2y_trk = new matrixTH1Ptr("bjtc_yield_trkCorr", npt, ncent);
		matrixTH1Ptr* m2j_jff = new matrixTH1Ptr("bjtc_jetShape_jff", npt, ncent);
		matrixTH1Ptr* m2y_jff = new matrixTH1Ptr("bjtc_yield_jff", npt, ncent);
		folder = ps->getPara<TString>("step3output_folder");
		auto f0 = TFile::Open(folder+"bjtc_correction.root");
//		std::cout<<folder+"bjtc_correction.root"<<std::endl;
		m2j_trk->autoLoad(f0);
		m2y_trk->autoLoad(f0);
		m2j_jff->autoLoad(f0);
		m2y_jff->autoLoad(f0);

		//step1: apply track correction
		edmJtcAnalyzer p2yield(*ps), p2shape(*ps);
		p2yield.initialization();
		p2shape.initialization();

		auto m2j_step1 = p2shape.m2_ring_corr("bjtc_jetShape_step1", m2j_tt_data, m2j_trk);
		auto m2y_step1 = p2yield.m2_ring_corr("bjtc_yield_step1", m2y_tt_data, m2y_trk);
		//step2: apply jff correction
		auto m2j_step2 = p2shape.m2_ring_corr("bjtc_jetShape_step2", m2j_step1, m2j_jff);
		auto m2y_step2 = p2yield.m2_ring_corr("bjtc_yield_step2", m2y_step1, m2y_jff);
		//step3: subtract the bkg
		p2shape.add2Step3( m2j_step2);
		p2yield.add2Step3( m2y_step2);
		p2shape.runProduce();
		p2yield.runProduce();
/*
		gQA->bookLegend();
		gQA->addLegendEntry("trk corrected", 0);
		gQA->addLegendEntry("jff corrected", 1);
		gAnaIO.saveCanvas(gQA->overlay(), "bjtc_yield_corrected");
		*/

		folder = ps->getPara<TString>("step3output_folder");
		auto wf = TFile::Open(folder+"bjtc_data_final.root", "recreate");
		p2shape.write();
		p2yield.write();
		wf->Close();
}

void bjtcAnalyzer_Step3::overlay_injtc_vs_bjtc(){
		matrixTH1Ptr* m2j_data = new matrixTH1Ptr("dr_signal_dijtc_jetShape_step2", npt, ncent);
		matrixTH1Ptr* m2y_data = new matrixTH1Ptr("dr_signal_dijtc_yield_step2", npt, ncent);
		TString folder = ps->getPara<TString>("step3output_folder");
		auto f = TFile::Open(folder+"inclJet_data_final.root");
		m2j_data->autoLoad(f);
		m2y_data->autoLoad(f);

		matrixTH1Ptr* m2j_bjtc = new matrixTH1Ptr("dr_signal_bjtc_jetShape_step2", npt, ncent);
		matrixTH1Ptr* m2y_bjtc = new matrixTH1Ptr("dr_signal_bjtc_yield_step2", npt, ncent); auto fb = TFile::Open(folder+"bjtc_data_final.root");
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
		float purity = 0.66;
		gQA->setXrange(0, 0.99);
		gQA->addRLine(1);
		gQA->setTitlePosition(0.625,0.85);
		matrixTH1Ptr* m2j_tg_rr = new matrixTH1Ptr("sig_mix_seagull_corrected_taggedBJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_tt_rr = new matrixTH1Ptr("sig_mix_seagull_corrected_taggedTrueBJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		matrixTH1Ptr* m2j_in_rr = new matrixTH1Ptr("sig_mix_seagull_corrected_inclJet_RecoJet_GenTrack_pTweighted_noCorr", npt, ncent);
		//matrixTH1Ptr* m2j_tg_rr = new matrixTH1Ptr("sig_mix_seagull_corrected_taggedBJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		//matrixTH1Ptr* m2j_tt_rr = new matrixTH1Ptr("sig_mix_seagull_corrected_taggedTrueBJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		//matrixTH1Ptr* m2j_in_rr = new matrixTH1Ptr("sig_mix_seagull_corrected_inclJet_RecoJet_RecoTrack_pTweighted_noCorr", npt, ncent);
		m2j_tg_rr->autoLoad(fdMC_rg);
		m2j_tt_rr->autoLoad(fdMC_rg);
		m2j_in_rr->autoLoad(fdMC_rg);

		m2j_in_rr->scale(1-purity);
		auto m2j_pu_rr = (*m2j_tg_rr)-(*m2j_in_rr);
		m2j_pu_rr->scale(1.0/purity);
		auto m2j_pu_rr_dr = jtc_utility::getDr("purified_rr", m2j_pu_rr, ndrbin,  drbins);
		auto m2j_tg_rr_dr = jtc_utility::getDr("tg_rr", m2j_tg_rr, ndrbin,  drbins);
		auto m2j_tt_rr_dr = jtc_utility::getDr("tt_rr", m2j_tt_rr, ndrbin,  drbins);
		gQA->addm2TH1Pair(m2j_tg_rr_dr, m2j_tt_rr_dr);
//		gQA->addm2TH1Pair(m2j_tg_rr, m2j_tt_rr);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair( "tagged b-jet", "t&t b-jet", 0);
		gQA->setLowerYrange(0.8 , 1.2);
		gAnaIO.saveCanvas(gQA->overlayR(), "validation_before_purification");
		gQA->addm2TH1Pair(m2j_pu_rr_dr, m2j_tt_rr_dr);
//		gQA->addm2TH1Pair(m2j_pu_rr, m2j_tt_rr);
		gQA->bookLegend(.5, .6, .9, .8);
		gQA->addLegendPair( "purified b-jet", "t&t b-jet", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "validation_purification");
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
