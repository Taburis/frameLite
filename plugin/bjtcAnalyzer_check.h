
#include "bjtcAnalyzer_base.h"

int ndetabin = 74;
Double_t  detabin[] = {-1.520000,-1.480000,-1.440000,-1.400000,-1.360000,-1.320000,-1.280000,-1.240000,-1.200000,-1.160000,-1.120000,-1.080000,-1.040000,-1.000000,-0.960000,-0.920000,-0.880000,-0.840000,-0.800000,-0.760000,-0.720000,-0.680000,-0.640000,-0.600000,-0.560000,-0.520000,-0.480000,-0.440000,-0.400000,-0.360000,-0.320000,-0.280000,-0.240000,-0.200000,-0.160000,-0.120000,-0.080000,-0.040000,0.040000,0.080000,0.120000,0.160000,0.200000,0.240000,0.280000,0.320000,0.360000,0.400000,0.440000,0.480000,0.520000,0.560000,0.600000,0.640000,0.680000,0.720000,0.760000,0.800000,0.840000,0.880000,0.920000,0.960000,1.000000,1.040000,1.080000,1.120000,1.160000,1.200000,1.240000,1.280000,1.320000,1.360000,1.400000,1.440000,1.480000, 1.52};

class bjtcAnalyzer_check : public bjtcAnalyzer_base{
		public : bjtcAnalyzer_check(ParaSet & pset) : bjtcAnalyzer_base(pset){};
				 void data_produce_refMC();
				 void side_band_mixing();
				 void ratio_ckpt(TString name, jtcTH1Player *, jtcTH1Player *, TString  , TString);
				 void overlay();
				 void plots();
				 void syst();
				 void an_cont_corr();
				 void make_cont_corr(jtcTH1Player*);
				 void contCorr_check();
				 void js_check();
				 void yield_check();
				 void trackingCorr();
				 void residualCorr();
				 void newAndOldPurity();
				 void jesRawCheck();
				 //		public : 
};

void bjtcAnalyzer_check::ratio_ckpt(TString name, jtcTH1Player *j1, jtcTH1Player *j2, TString leg1, TString leg2){
		gQA->addm2TH1Pair(j1, j2);
		gQA->addRLine(1);
		gQA->setLowerYrange(0.5, 1.5);
		gQA->bookLegend(0.6, 0.6, 0.95, 0.8);
		gQA->addLegendPair(leg1, leg2, 0);
		gAnaIO.saveCanvas(gQA->overlayR(), name);
}

void bjtcAnalyzer_check::make_cont_corr(jtcTH1Player* j2){
		for(int j=0; j<j2->Ncol(); ++j){
				for(int i=0; i<j2->Nrow(); i++){
						auto h = j2->at(i,j);
						for(auto k=1; k<h->GetNbinsX()+1; k++){
								if(h->GetXaxis()->GetBinCenter(k)>0.45) h->SetBinContent(k,1);
						}
				}
		}
}

void bjtcAnalyzer_check::an_cont_corr(){
		auto fcorr = TFile::Open("/Users/tabris/frameLite/output/step3/bjtc_correction.root");
		auto cocorr = new jtcTH1Player("bjtc_jetShape_contCorr");
		auto * m2j_in_gg = new jtcTH1Player("dr_raw_inclJet_GenJet_GenTrack_pTweighted_noCorr");
		auto * m2j_co_gg = new jtcTH1Player("dr_raw_contJet_GenJet_GenTrack_pTweighted_noCorr");
		//auto * m2j_in_gg = new jtcTH1Player("dr_mix_seagull_corrected_inclJet_GenJet_GenTrack_pTweighted_noCorr");
		//auto * m2j_co_gg = new jtcTH1Player("dr_mix_seagull_corrected_contJet_GenJet_GenTrack_pTweighted_noCorr");
		m2j_in_gg->autoLoad(fdMC_gg);
		m2j_co_gg->autoLoad(fdMC_gg);
		auto ratio = (*m2j_co_gg)%(*m2j_in_gg);
		//		ratio->setAxisRange(.06, .99, "X");
		auto ratio2 = (jtcTH1Player*) ratio->clone("ratio2");
		ratio2->at(3, 0 ) ->SetAxisRange(0.051, .99, "X");
		ratio2->smooth(1, "R");

		ratio2->doChi2Test((jtcTH1Player*) ratio, "CHI2/NDF");
		make_cont_corr(ratio2);
		cocorr->autoLoad(fcorr);
		//		gQA->addm2TH1(cocorr);
		gQA->addm2TH1(ratio2, "hist");
		gQA->addm2TH1(ratio);

		gQA->bookLegend(0.5, 0.55, 0.9, 0.75);
		gQA->addLegendEntry("cont. corr.", 0, "l");
		gQA->addLegendEntry("cont/incl", 1);
		gQA->addhLine(1);
		gQA->setXrange(0, .99);
		gQA->setYrange(0.5, 1.5);
		//gQA->setXrange(0, 2.499);
		//gQA->overlay();	
		gAnaIO.saveCanvas(gQA->overlay(), "bjtc_contCorr");
}

void bjtcAnalyzer_check::data_produce_refMC(){
		auto mc_in_raw = new jtcTH1Player("dr_raw_inclJet_GenJet_GenTrack_pTweighted_noCorr");
		auto da_in_raw = new jtcTH1Player("dr_raw_inclJet_Data_pTweighted_noCorr");
		mc_in_raw->autoLoad(fdMC_gg);
		da_in_raw->autoLoad(fdata);
		ratio_ckpt("incl_ckpt_raw", da_in_raw, mc_in_raw, "data", "MC");
		gQA->setXrange(0, 2.499);

		auto mc_in_me = new jtcTH1Player("dr_mix_seagull_corrected_inclJet_GenJet_GenTrack_pTweighted_noCorr");
		auto da_in_me = new jtcTH1Player("dr_mix_seagull_corrected_inclJet_Data_pTweighted_noCorr");
		mc_in_me->autoLoad(fdMC_gg);
		da_in_me->autoLoad(fdata);
		ratio_ckpt("incl_ckpt_me", da_in_me, mc_in_me, "data", "MC");

		auto m2da_in_raw = new jtcTH1Player("sig_mix_seagull_corrected_inclJet_Data_pTweighted_noCorr");
		m2da_in_raw->autoLoad(fdata);
		auto fcorr = TFile::Open("/Users/tabris/frameLite/output/step3/inclJet_correction.root");
		auto trk_corr = new jtcTH1Player("dijtc_jetShape_trkCorr");
		trk_corr->autoLoad(fcorr);
		auto da_in_trk = (jtcTH1Player*) m2da_in_raw->clone("incl_trkcorr");
		da_in_trk->ring_corr(trk_corr , 2.5);
		auto dr_in_trk = da_in_trk->drIntegral("dr_in_trk", ndrbin, drbins);
		ratio_ckpt("incl_ckpt_trk", dr_in_trk, mc_in_me, "data", "MC");


		auto jff = new jtcTH1Player("dijtc_jetShape_jff");
		jff->autoLoad(fcorr);
		auto da_in_jff = (jtcTH1Player*) da_in_trk->clone("incl_jff");
		da_in_jff->ring_corr(jff , 2.5);
		auto dr_in_jff = da_in_jff->drIntegral("dr_in_jff", ndrbin, drbins);
		ratio_ckpt("incl_ckpt_jff", dr_in_jff, mc_in_me, "data", "MC");

		//		gQA->setLowerYrange(0., 3);
		auto mc_in_sig = new jtcTH1Player("dr_signal_inclJet_GenJet_GenTrack_pTweighted_noCorr");
		mc_in_sig->autoLoad(fdMC_gg);
		auto da_in_final = da_in_jff->bkgSub("incl_final", 1.5, 2.5);
		auto dr_in_final = da_in_final->drIntegral("dr_in_final", ndrbin, drbins);
		ratio_ckpt("incl_ckpt_final", dr_in_final, mc_in_sig, "data", "MC");
		//ratio_ckpt("incl_ckpt_final", mc_in_sig, dr_in_final, "MC", "data");
}

void bjtcAnalyzer_check::side_band_mixing(){
		auto fmc_gg = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/pp_28Jan19/bJTC_PYTHIA6_GenGen_5TeV_dijetMC_WTAaxis_csvV2p9_28Jan19.root");
		auto fmc_rg = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/pp_28Jan19/bJTC_PYTHIA6_RecGen_5Tev_dijetMC_WTAaxis_csvV2p9_28Jan19.root");
		auto fmc_rr = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/pp_28Jan19/bJTC_PYTHIA6_RecRec_5TeV_dijetMC_WTAaxis_csvV2p9_28Jan19.root");
		auto frd = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/pp_28Jan19/bJTC_pp_data_5TeV_wtaAxis_28Jan19.root");

		auto mc_gg = new jtcTH1Player("inclJet_GenJet_GenTrack_pTweighted_noCorr");
		auto mc_rg = new jtcTH1Player("inclJet_RecoJet_GenTrack_pTweighted_noCorr");
		auto mc_rr = new jtcTH1Player("inclJet_RecoJet_RecoTrack_pTweighted_noCorr");
		auto rdata = new jtcTH1Player("inclJet_Data_pTweighted_noCorr");
		mc_gg->autoLoad(fmc_gg);
		mc_rg->autoLoad(fmc_rg);
		mc_rr->autoLoad(fmc_rr);
		rdata->autoLoad(frd);
		auto jtpt_gg = (TH1D*)fmc_gg->Get("inclJet_corrpt_0");
		float ss = jtpt_gg->Integral();
		mc_gg->scale(1.0/ss);
		mc_gg->invariant();

		auto jtpt_rg = (TH1D*)fmc_rg->Get("inclJet_corrpt_0");
		ss = jtpt_rg->Integral();
		mc_rg->scale(1.0/ss);
		mc_rg->invariant();

		auto jtpt_rr = (TH1D*)fmc_rr->Get("inclJet_corrpt_0");
		ss = jtpt_rr->Integral();
		mc_rr->scale(1.0/ss);
		mc_rr->invariant();

		auto jtpt_data = (TH1D*)frd->Get("inclJet_corrpt_0");
		ss = jtpt_data->Integral();
		rdata->scale(1.0/ss);
		rdata->invariant();

		auto data_sig = rdata->getSignal_ME_based("inclJet_data_signal", 1.4, 1.8, 1);
		auto dr_sig = data_sig->drIntegral("inclJet_dr_data_signal", ndrbin, drbins);

		//correction
		auto mc_rr_step1 = mc_rr->getSignal_ME_based("inclJet_rr_step1", 1.4, 1.8, 0);
		auto dr_rr_step1 = mc_rr_step1->drIntegral("dr_inclJet_rr_step1", ndrbin, drbins);
		auto mc_rg_step1 = mc_rg->getSignal_ME_based("inclJet_rg_step1", 1.4, 1.8, 0);
		auto dr_rg_step1 = mc_rg_step1->drIntegral("dr_inclJet_rg_step1", ndrbin, drbins);
		auto mc_gg_step1 = mc_gg->getSignal_ME_based("inclJet_gg_step1", 1.4, 1.8, 0);
		auto dr_gg_step1 = mc_gg_step1->drIntegral("dr_inclJet_gg_step1", ndrbin, drbins);
		auto trk = (jtcTH1Player*)((*dr_rr_step1)/(*dr_rg_step1));
		trk->smooth();
		auto fcorr = TFile::Open("/Users/tabris/frameLite/output/step3/inclJet_correction.root");
		auto trkcorr = new jtcTH1Player("dijtc_jetShape_trkCorr");
		trkcorr->autoLoad(fcorr);
		ratio_ckpt("incl_ckpt_trkcorr", trk, trkcorr, "side based", "nominal");

		auto jff = (jtcTH1Player*)((*dr_rg_step1)/(*dr_gg_step1));
		jff->smooth();
		auto jffcorr = new jtcTH1Player("dijtc_jetShape_jff");
		jffcorr->autoLoad(fcorr);
		ratio_ckpt("incl_ckpt_jffcorr", jff, jffcorr, "side based", "nominal");

		auto data_step1 = rdata->getSignal_ME_based("inclJet_data_signal", 1.4, 1.8, 0);
		data_step1->ring_corr(trk, 2.5);	
		data_step1->ring_corr(jff, 2.5);	
		auto data_sig_m2 = data_step1->bkgSub("inclJet_data_final_signal", 1.5, 2.5);
		auto dr_data_sig = data_sig_m2->drIntegral("dr_inclJet_data_signal", ndrbin, drbins);

		auto mc_gg_sig2 = mc_gg_step1->bkgSub("inclJet_gg_signal", 1.5, 2.5);
		auto dr_gg_sig2 = mc_gg_sig2->drIntegral("inclJet_dr_gg_signal", ndrbin, drbins);

		ratio_ckpt("incl_ckpt_final", dr_data_sig, dr_gg_sig2, "side based", "mc");

		auto data_nominal = new jtcTH1Player("dr_signal_dijtc_jetShape_step3");
		auto ffinal = TFile::Open("/Users/tabris/frameLite/output/step3/inclJet_data_final.root");
		data_nominal->autoLoad(ffinal);
		ratio_ckpt("incl_ckpt_data", dr_data_sig, data_nominal, "data(side)", "data(ME)");

		auto data_signal = new jtcTH1Player("dr_signal_inclJet_Data_pTweighted_noCorr");
		data_signal->autoLoad(fdata);

		auto mc_gg_prox = mc_gg->projX("prox", 1.4, 1.8, "e", 1);
		//		auto mc_me_prox = mc_gg_me->projX("me_prox", 1.4, 1.8, "e", 1);
		//		gQA->addm2TH1(mc_gg_prox);
		//		gQA->addm2TH1(mc_me_prox);
		//		gQA->bookLegend();
		//		gQA->setXrange(-2.5, 2.499);
		//		gQA->addLegendEntry("true band", 0);
		//		gQA->addLegendEntry("ME", 1);
		//		gQA->overlay();

		gAnaIO.saveCanvas(gQA->jtc_check001(*mc_gg_sig2), "incl_gg_check001");
		auto dr_gg_sig1 = mc_gg_step1->drIntegral("inclJet_dr_gg_step1", ndrbin, drbins);
		auto mc_in_sig = new jtcTH1Player("dr_signal_inclJet_GenJet_GenTrack_pTweighted_noCorr");
		auto mc_in_me = new jtcTH1Player("dr_mix_seagull_corrected_inclJet_GenJet_GenTrack_pTweighted_noCorr");
		mc_in_sig->autoLoad(fdMC_gg);
		mc_in_me->autoLoad(fdMC_gg);
		ratio_ckpt("incl_ckpt_side_corrected_step1", dr_gg_sig1, mc_in_me, "side", "MC");
		ratio_ckpt("incl_ckpt_side_corrected_step2", dr_gg_sig2, mc_in_sig, "MC(side)", "MC(ME)");

		ratio_ckpt("incl_ckpt_ratio", dr_data_sig, dr_gg_sig2, "data(side)", "MC(side)");
		auto ratio1 = (*dr_data_sig)/(*dr_gg_sig2);
		auto ratio2 = (*data_nominal)/(*mc_in_sig);
		gQA->addm2TH1Pair(ratio1, ratio2);
		gQA->addRLine(1);
		gQA->setLowerYrange(0.5, 1.5);
		gQA->bookLegend(0.6, 0.6, 0.95, 0.8);
		gQA->addLegendPair("side", "ME", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "double_ratio");


		auto dr_data_all = data_nominal->contractX("dr_data_all");
		auto dr_mc_all  = dr_gg_sig2->contractX("dr_data_mc");
		gQA->addm2TH1Pair(dr_data_all, dr_mc_all);
		gQA->bookLegend(0.6, 0.6, 0.95, 0.8);
		gQA->addLegendPair("incl.(WTA)", "MC", 0);
		gQA->ncol = 1;
		gQA->nrow = 1;
		gQA->setXrange(.0, .99);
		gQA->setLowerYrange(.5, 2.);
		bjtc_pp_config::ptString[0]= "p^{track}_{T} > 1 GeV";
		gAnaIO.saveCanvas(gQA->overlayR(), "integratedJS");

}

void bjtcAnalyzer_check::plots(){
		auto fbcorr = TFile::Open("/Users/tabris/frameLite/output/step3/bjtc_correction.root");
		auto fdcorr = TFile::Open("/Users/tabris/frameLite/output/step3/inclJet_correction.root");
		auto j2dtrk = new jtcTH1Player("dijtc_jetShape_trkCorr");
		auto j2btrk = new jtcTH1Player("bjtc_jetShape_trkCorr");
		auto j2djff = new jtcTH1Player("dijtc_jetShape_jff");
		auto j2bjff = new jtcTH1Player("bjtc_jetShape_jff");
		j2dtrk->autoLoad(fdcorr);
		j2btrk->autoLoad(fbcorr);
		j2djff->autoLoad(fdcorr);
		j2bjff->autoLoad(fbcorr);
		gQA->setXrange(.0, .99);
		gQA->setYrange(.5, 1.1);
		gQA->addhLine(1);
		gQA->addm2TH1(j2dtrk);
		gQA->addm2TH1(j2btrk);
		gQA->bookLegend(0.5, 0.3, 0.85, 0.5);
		gQA->addLegendEntry("incl. trk corr.",0);
		gQA->addLegendEntry("b-jet trk corr.",1);
		gAnaIO.saveCanvas(gQA->overlay(), "trkCorr_comparison");

		gQA->addm2TH1(j2djff);
		gQA->addm2TH1(j2bjff);
		gQA->setYrange(.5, 1.5);
		gQA->bookLegend(0.2, 0.55, 0.85, 0.85);
		gQA->addLegendEntry("incl. residual jff corr.",0);
		gQA->addLegendEntry("b-jet residual jff+bias corr.",1);
		gAnaIO.saveCanvas(gQA->overlay(), "jffCorr_comparison");
}

void bjtcAnalyzer_check::contCorr_check(){
		TString folder = ps->getPara<TString>("step3output_folder");
		auto f0 = TFile::Open(folder+"bjtc_correction.root");
		auto j2_in_gg = new jtcTH1Player("sig_mix_seagull_corrected_inclJet_GenJet_GenTrack_pTweighted_noCorr");
		auto j2_co_gg = new jtcTH1Player("sig_mix_seagull_corrected_contJet_GenJet_GenTrack_pTweighted_noCorr");
		j2_in_gg->autoLoad(fdMC_gg);
		j2_co_gg->autoLoad(fdMC_gg);
		auto px_in_gg = j2_in_gg->projX("projEta_incl", -1, 1, "e", 0);
		auto px_co_gg = j2_co_gg->projX("projEta_cont", -1, 1, "e", 0);
		gQA->addm2TH1Pair(px_co_gg, px_in_gg);
		gQA->bookLegend(0.6, 0.6, 0.95, 0.8);
		gQA->setXrange(-.5, .499);
		gQA->addRLine(1);
		gQA->setLowerYrange(0.5, 1.5);
		gQA->addLegendPair("cont.", "incl.", 0);
		gAnaIO.saveCanvas(gQA->overlayR(), "cont_projEta");

		jtcTH1Player* m2j_in_data = new jtcTH1Player("sig_mix_corrected_inclJet_Data_pTweighted_noCorr", npt, ncent);
		m2j_in_data->autoLoad(fdata);
		auto m2j_co_data = (jtcTH1Player*) m2j_in_data->clone("m2j_co_data");
		auto m2j_cont = new jtcTH1Player("bjtc_jetShape_contCorr", npt, ncent);
		m2j_cont->autoLoad(f0);
		m2j_co_data->ring_corr(m2j_cont, 0.4, 1);
		auto deta_co_data = m2j_co_data->projX("deta_contamination", -1, 1, "e", 0);
		auto deta_in_data = m2j_in_data->projX("deta_incl", -1, 1, "e", 0);
		deta_co_data->setAxisRange(-1, .99, "X");
		deta_in_data->setAxisRange(-1, .99, "X");
		gQA->addm2TH1(deta_co_data);
		gQA->addm2TH1(deta_in_data);
		gQA->bookLegend(0.6, 0.6, 0.95, 0.8);
		gQA->addLegendEntry("cont.", 0);
		gQA->addLegendEntry("incl.", 1);
		gQA->setXrange(-1, .99);
		//		gQA->addhLine(0);
		gAnaIO.saveCanvas(gQA->overlay(), "contTemp_projEta");

		auto j2_in_dr = new jtcTH1Player("dr_signal_inclJet_GenJet_GenTrack_pTweighted_noCorr");
		auto j2_co_dr = new jtcTH1Player("dr_signal_contJet_GenJet_GenTrack_pTweighted_noCorr");
		//auto j2_in_dr = new jtcTH1Player("dr_mix_seagull_corrected_inclJet_GenJet_GenTrack_pTweighted_noCorr");
		//auto j2_co_dr = new jtcTH1Player("dr_mix_seagull_corrected_contJet_GenJet_GenTrack_pTweighted_noCorr");
		j2_in_dr->autoLoad(fdMC_gg);
		j2_co_dr->autoLoad(fdMC_gg);
		//		j2_in_dr->smooth();
		//		j2_co_dr->smooth();
		gQA->addm2TH1Pair(j2_in_dr, j2_co_dr);
		gQA->setXrange(0, .99);
		gAnaIO.saveCanvas(gQA->overlayR(), "contTemp_dr2");
}

void bjtcAnalyzer_check::yield_check(){
		auto dfinal = TFile::Open("/Users/tabris/frameLite/output/step3/inclJet_data_final.root");
		auto bfinal = TFile::Open("/Users/tabris/frameLite/output/step3/bjtc_data_final.root");
		auto j2_tg = new jtcTH1Player("signal_bjtc_yield_step3");
		auto j2_in = new jtcTH1Player("dijtc_yield_step3");
		j2_tg->autoLoad(bfinal);
		j2_in->autoLoad(dfinal);
		auto deta_tg0 = j2_tg->projX("deta_bjtc_yield0", -1, 1, "e", 0);
		auto deta_in0 = j2_in->projX("deta_incl_yield0", -1, 1, "e", 0);
		auto deta_tg = deta_tg0->rebin("deta_bjtc_yield", ndetabin, detabin);
		auto deta_in = deta_in0->rebin("deta_incl_yield", ndetabin, detabin);
		deta_tg->setAxisRange(-1, 1, "X");
		deta_in->setAxisRange(-1, 1, "X");
		auto data_in = new jtcTH1Player("dr_signal_dijtc_yield_step3");
		data_in->autoLoad(dfinal);
		auto data_b = new jtcTH1Player("dr_signal_bjtc_yield_step3");
		data_b->autoLoad(bfinal);
		auto diff = (*data_b)-(*data_in);
		gQA->addm2TH1(diff);
		gQA->setXrange(0, .99);
		gQA->setYrange(-5, 5);
		gQA->addhLine(0);
		gAnaIO.saveCanvas(gQA->overlay(), "yield_diff");

		gQA->addm2TH1(deta_tg);
		gQA->addm2TH1(deta_in);
		gQA->bookLegend(0.65, 0.6, 0.925, 0.8);
		gQA->addLegendEntry("b-jet", 0);
		gQA->addLegendEntry("incl. jet", 1);
		gQA->setXrange(-1, .99);
		gQA->autoYrange = 1;
		gQA->fixYrange = 0;
		gQA->addhLine(0);
		gAnaIO.saveCanvas(gQA->overlay(), "yield_deta");
}

void bjtcAnalyzer_check::syst(){
		auto bfmc1 = TFile::Open("/Users/tabris/frameLite/output/step2/Signal_PYTHIA6_bjetSample_JER_RecoJet_GenTrack.root");
		auto bfmc2 = TFile::Open("/Users/tabris/frameLite/output/step2/Signal_PYTHIA6_bjetSample_allJets_RecoJet_GenTrack.root");
		auto j2_tg1 = new jtcTH1Player("dr_signal_taggedBJet_RecoJet_GenTrack_pTweighted_noCorr");
		auto j2_tg2 = new jtcTH1Player("dr_signal_taggedBJet_RecoJet_GenTrack_pTweighted_noCorr");
		j2_tg1->autoLoad(bfmc1);
		j2_tg2->autoLoad(bfmc2);
		gQA->addm2TH1Pair(j2_tg1, j2_tg2);
		gQA->bookLegend(0.65, 0.6, 0.925, 0.8);
		gQA->addLegendPair("Smeared jet", "Reco. jet", 0);
		gQA->setXrange(0, .99);
		gQA->setLowerYrange(0.9, 1.1);
		gQA->autoYrange = 1;
		gQA->fixYrange = 0;
		gQA->addhLine(0);
		gAnaIO.saveCanvas(gQA->overlayR(), "JER");

}

void bjtcAnalyzer_check::js_check(){
		auto dfinal = TFile::Open("/Users/tabris/frameLite/output/step3/inclJet_data_final.root");
		auto bfinal = TFile::Open("/Users/tabris/frameLite/output/step3/bjtc_data_final.root");
		auto j2_tg = new jtcTH1Player("signal_bjtc_jetShape_step3");
		auto j2_in = new jtcTH1Player("dijtc_jetShape_step3");
		j2_tg->autoLoad(bfinal);
		j2_in->autoLoad(dfinal);
		auto deta_tg0 = (jtcTH1Player*) j2_tg->projX("deta_bjtc_js0", -1, 1, "e", 0);
		auto deta_in0 = (jtcTH1Player*) j2_in->projX("deta_incl_js0", -1, 1, "e", 0);
		auto deta_tg = deta_tg0->rebin("deta_bjtc_js", ndetabin, detabin);
		auto deta_in = deta_in0->rebin("deta_incl_js", ndetabin, detabin);
		deta_tg->setAxisRange(-1, 1, "X");
		deta_in->setAxisRange(-1, 1, "X");

		gQA->addm2TH1(deta_tg);
		gQA->addm2TH1(deta_in);
		gQA->bookLegend(0.65, 0.6, 0.925, 0.8);
		gQA->addLegendEntry("b-jet", 0);
		gQA->addLegendEntry("incl. jet", 1);
		gQA->setXrange(-1, .99);
		gQA->autoYrange = 1;
		gQA->fixYrange = 0;
		gQA->addhLine(0);
		gAnaIO.saveCanvas(gQA->overlay(), "js_deta");
}

void bjtcAnalyzer_check::trackingCorr(){
		auto bfinal = TFile::Open("/Users/tabris/frameLite/output/step3/bjtc_correction.root");
		auto j2_trk = new jtcTH1Player("bjtc_yield_trkCorr_*_*");
		j2_trk->autoLoad(bfinal);
		auto j2_err = (jtcTH1Player*) j2_trk->clone("j2_err");
		j2_err->absError(.04);
		gQA->addm2TH1(j2_trk);
		gQA->addm2TH1Error(j2_err);
		gQA->addhLine(1);
		gQA->setXrange(0, .99);
		gQA->Ytitle = "tracking efficiency";
//		gQA->overlay();
		gAnaIO.saveCanvas(gQA->overlay(), "trkCorr");
}

void bjtcAnalyzer_check::residualCorr(){
		auto bfinal = TFile::Open("/Users/tabris/frameLite/output/step3/bjtc_correction.root");
		auto j2_corr = new jtcTH1Player("bjtc_yield_jff");
		j2_corr->autoLoad(bfinal);
		auto j2_err = (jtcTH1Player*) j2_corr->clone("j2_err");
		j2_err->absError(.05);
		gQA->addm2TH1(j2_corr);
		gQA->addm2TH1Error(j2_err);
		gQA->addhLine(1);
		gQA->setXrange(0, .99);
		gQA->setYrange(0.6, 1.4);
		gAnaIO.saveCanvas(gQA->overlay(), "jffCorr");
}

void bjtcAnalyzer_check::newAndOldPurity(){
		TString folder0 = "/Users/tabris/frameLite/output/step3/";
		TString folder1 = "/Users/tabris/frameLite/output/step3_p66/";
		auto f0 = TFile::Open(folder0+"bjtc_data_final.root");
		auto f1 = TFile::Open(folder1+"bjtc_data_final.root");
		auto data0 = new jtcTH1Player("dr_signal_bjtc_jetShape_step3_*_*");
		auto data1 = new jtcTH1Player("dr_signal_bjtc_jetShape_step3_*_*");
		auto err0  = new jtcTH1Player("dr_bjtc_jetShape_systError_*_*", npt, ncent);
		auto err1  = new jtcTH1Player("dr_bjtc_jetShape_systError_*_*", npt, ncent);
		err0->autoLoad(f0);
		err1->autoLoad(f1);
		data0->autoLoad(f0);
		data1->autoLoad(f1);
		err0->add_frac_error(0.05);
		err1->add_frac_error(0.05);
		auto da0 = data0->contractX("dr_data0");
		auto da1 = data1->contractX("dr_data1");
		auto er0 = err0->contractX("dr_error0");
		auto er1 = err1->contractX("dr_error1");

//		gQA->addm2TH1(da0);
//		gQA->addm2TH1Error(er0);
//		gQA->addm2TH1(da1);
//		gQA->addm2TH1Error(er1);
		gQA->addm2TH1Pair(da0, da1);
		gQA->addm2TH1ErrorPair(er0, er1);
		 gQA->bookLegend(0.55, 0.55, 0.9, 0.83);
		 gQA->setLowerYrange(0.85, 1.15);
        gQA->addLegendPair("misid.(MC)", "misid.(data)", 0 );
		gQA->addRLine(1);
//        gQA->addLegendEntry("P=0.7", 0);
//        gQA->addLegendEntry("P=0.65", 1);
        gQA->ncol = 1;
        gQA->nrow = 1;
        gQA->setXrange(.0, .99);
        bjtc_pp_config::ptString[0]= "p^{track}_{T} > 1 GeV";
        gQA->Ytitle = "P(#Delta r)";
        gAnaIO.saveCanvas(gQA->overlayR(), "purity_overlay");

}

void bjtcAnalyzer_check::jesRawCheck(){
		TString folder0 = "/Users/tabris/frameLite/output/step2/";
		TFile * fjes = TFile::Open(folder0+"Signal_PYTHIA6_dijetSample_allJets_JES117_RecoJet_GenTrack.root");
		TFile * fori = TFile::Open(folder0+"Signal_PYTHIA6_dijetSample_allJets_RecoJet_GenTrack.root");
		auto jes117 = new jtcTH1Player("dr_signal_inclJet_RecoJet_GenTrack_pTweighted_noCorr_*_*");
		auto jes120 = new jtcTH1Player("dr_signal_inclJet_RecoJet_GenTrack_pTweighted_noCorr_*_*");
		jes117->autoLoad(fjes);
		jes120->autoLoad(fori);
//		auto dr_jes117 = jes117->drIntegral("dr_jes117", ndrbin, drbins);
//		auto dr_jes120 = jes120->drIntegral("dr_jes120", ndrbin, drbins);
		auto dr_jes117all = jes117->contractX("dr_jes117all");
		auto dr_jes120all = jes120->contractX("dr_jes120all");
		gQA->addm2TH1Pair(dr_jes117all, dr_jes120all);
//		gQA->addm2TH1(jes120);
		gQA->bookLegend(0.6, 0.6, 0.95, 0.8);
		gQA->addLegendPair("shifted", "origin", 0);
//		gQA->addLegendEntry("117", 0);
//		gQA->addLegendEntry("120", 1);
		gQA->fixYrange=0;
		gQA->ncol = 1;
		gQA->nrow = 1;
		gQA->setXrange(.0, .99);
//		gQA->setLowerYrange(.5, 2.);
		bjtc_pp_config::ptString[0]= "p^{track}_{T} > 1 GeV";
		gAnaIO.saveCanvas(gQA->overlayR(), "JESbjetDiff");
}

