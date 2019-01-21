
#include "analyzerIOServer.h"
#include "jtcQaMonitor.h"
#include "edmJtcAnalyzer.h"
#include "bjtcAnalyzer_Step2.h"
#include "../run/ppjtc_cfg.h"
#include "jtcTH1Player.h"

void jetQA(){
		//plotting the jet pt eta phi comparision between the EScheme jet axis and WTA jet axis for inclusive jets.
		init();
		auto fesm = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/myCorrelation/CSVv2/bJTC_GenGen_5TeV_25mix_csvV2_drCorr_12June18.root");
		auto fwta = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/pp_18Dec18/bJTC_PYTHIA6_GenGen_5TeV_dijetMC_wtaAxis_11Jan19.root");
		auto hwtaeta = (TH1D*) fwta->Get("inclJet_corrpt_0");
		auto m2esm= new matrixTH1Ptr("esmjet_property",  1, 3);
		auto m2wta= new matrixTH1Ptr("wtajet_property",  1, 3);
		m2esm->add((TH1D*) fesm->Get("inclJet_corrpt_0"), 0, 0);
		m2wta->add((TH1D*) fwta->Get("inclJet_corrpt_0"), 0, 0);
		m2esm->add((TH1D*) fesm->Get("inclJet_eta_0"), 0, 1);
		m2wta->add((TH1D*) fwta->Get("inclJet_eta_0"), 0, 1);
		m2esm->add((TH1D*) fesm->Get("inclJet_phi_0"), 0, 2);
		m2wta->add((TH1D*) fwta->Get("inclJet_phi_0"), 0, 2);
		m2esm->at(0,1)->Rebin(4);
		m2wta->at(0,1)->Rebin(4);
		m2esm->at(0,2)->Rebin(4);
		m2wta->at(0,2)->Rebin(4);
		m2esm->at(0,0)->GetXaxis()->SetTitle("p_{T}");
		m2esm->at(0,1)->GetXaxis()->SetTitle("#eta");
		m2esm->at(0,2)->GetXaxis()->SetTitle("#phi");
		m2wta->at(0,0)->GetXaxis()->SetTitle("p_{T}");
		m2wta->at(0,1)->GetXaxis()->SetTitle("#eta");
		m2wta->at(0,2)->GetXaxis()->SetTitle("#phi");
		gQA->nrow = 1;
		gQA->ncol = 3;
	//	gQAmonitor->addm2TH1(m2wta);
	//	gQAmonitor->addm2TH1(m2esm);
//		auto raito = (*m2wta)/ (*m2esm);
		//gQAmonitor->overlayR("" );
		m2wta->normalized_by_area();
		m2esm->normalized_by_area();
		gQA->addm2TH1Pair(m2wta, m2esm);
		gQA->addRLine(1.0);
		gQA->setLowerYrange(0.9, 1.1);
		gAnaIO.saveCanvas(gQA->overlayR(), "jetQA_esmAxis_vs_wtaAxis");
}
