
#include "bjtcAnalyzer_base.h"

class bjtcAnalyzer_check : public bjtcAnalyzer_base{
		public : bjtcAnalyzer_check(ParaSet & pset) : bjtcAnalyzer_base(pset){};
				 void data_produce_refMC();
				 void side_band_mixing();
				 void ratio_ckpt(TString name, jtcTH1Player *, jtcTH1Player *, TString  , TString);
};

void bjtcAnalyzer_check::ratio_ckpt(TString name, jtcTH1Player *j1, jtcTH1Player *j2, TString leg1, TString leg2){
		gQA->addm2TH1Pair(j1, j2);
		gQA->addRLine(1);
		gQA->setLowerYrange(0.0, 2);
		gQA->bookLegend(0.6, 0.65, 0.9, 0.8);
		gQA->addLegendPair(leg1, leg2, 0);
		gAnaIO.saveCanvas(gQA->overlayR(), name);
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

		auto jtpt_data = (TH1D*)frd->Get("inclJet_corrpt_0");
		ss = jtpt_data->Integral();
		rdata->scale(1.0/ss);
		rdata->invariant();
		/*
		auto data_me = mc_gg->getSideMixTable("inclJet_data_sideMix", 1.4, 1.8);
		jtcTH1Player* data_me_step1 = (jtcTH1Player*) ((*rdata)/(*data_me));
		auto data_sig2 = data_me_step1->bkgSub("inclJet_data_signal", 1.5, 2.5);
		*/
		auto data_sig2 = rdata->getSignal_ME_based("inclJet_data_signal", 1.5, 2.5, 1);
		auto dr_sig2 = data_sig2->drIntegral("inclJet_dr_data_signal", ndrbin, drbins);

		auto data_signal = new jtcTH1Player("dr_signal_inclJet_Data_pTweighted_noCorr");
		data_signal->autoLoad(fdata);
		ratio_ckpt("incl_ckpt_side_corrected_data", dr_sig2, data_signal, "side corrected", "data");

		auto mc_gg_me = mc_gg->getSideMixTable("inclJet_gg_sideMix", 1.4, 1.8);
		jtcTH1Player* mc_gg_step1 = (jtcTH1Player*) ((*mc_gg)/(*mc_gg_me));
		auto mc_gg_sig2 = mc_gg_step1->bkgSub("inclJet_gg_signal", 1.5, 2.5);
		auto dr_gg_sig2 = mc_gg_sig2->drIntegral("inclJet_dr_gg_signal", ndrbin, drbins);

		auto mc_gg_prox = mc_gg->projX("prox", 1.4, 1.8, "e", 1);
		auto mc_me_prox = mc_gg_me->projX("me_prox", 1.4, 1.8, "e", 1);
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
		ratio_ckpt("incl_ckpt_side_corrected_step2", dr_gg_sig2, mc_in_sig, "side", "MC");
}
