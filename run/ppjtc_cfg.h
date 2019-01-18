
#include "edmJtcUtility.h"
#include "ParaSet.h"
#include "TString.h"
namespace bjtc_pp_config{
		TString ptString[] = {"1 < p_{T} < 2 GeV", "2 < p_{T} < 3 GeV", "3 < p_{T} < 4 GeV", "4 < p_{T} < 8 GeV", "8 < p_{T} < 12 GeV", "p_{T} > 12 GeV"};
		jtc_utility::index2d arrange_pp_bjtc(int i, int j){
				jtc_utility::index2d x;
				x.i1 = i/3;
				x.i2 = i%3;
				return x;
		}
		TString pTtitle(int i, int j){
				return ptString[i];
		}
}

ParaSet makePSet_edmJtcDefault(){
		// config for step2 and step3
		ParaSet st("edmJtcAnalzyer_pset");
		double jetptbin[] = {110, 120, 136, 152, 168, 184, 200, 216, 232, 248, 264, 280, 296, 312, 328, 344, 360, 380, 400, 432, 500};
		//		cout<<"size is "<<sizeof(jetptbin)/sizeof(jetptbin[0])<<endl;
		//float drbin[] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.};
		//st.setPara<int>("ndr", 14);
		float drbin[] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.45,0.6,0.8,1. };
		st.setPara<int>("ndr", 11);
		st.setParaVector<float>("drbins", sizeof(drbin)/sizeof(drbin[0]), drbin);
		st.setParaVector<double>("jet_pt_bin", sizeof(jetptbin)/sizeof(jetptbin[0]), jetptbin);
		st.setPara<int>("npt", 6);
		st.setPara<int>("ncent", 1);
		return st;
}

ParaSet makePSet_bjtc_pp_step23(){
		TString inputfolder = "/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/pp_18Dec18/";
		TString fbjetMC_rg = "bJTC_PYTHIA6_RecGen_5TeV_bjetMC_wtaAxis_7Jan19.root";
		TString fbjetMC_gg = "bJTC_PYTHIA6_GenGen_5TeV_bjetMC_wtaAxis_7Jan19.root";
		TString fbjetMC_rr = "bJTC_PYTHIA6_RecRec_5TeV_bjetMC_wtaAxis_7Jan19.root";
		TString fdijetMC_rg = "bJTC_PYTHIA6_RecGen_5TeV_dijetMC_wtaAxis_11Jan19.root";
		TString fdijetMC_gg = "bJTC_PYTHIA6_GenGen_5TeV_dijetMC_wtaAxis_11Jan19.root";
		TString fdijetMC_rr = "bJTC_PYTHIA6_RecRec_5TeV_dijetMC_wtaAxis_11Jan19.root";
		TString fdata = "bJTC_pp_data_5TeV_wtaAxis_7Jan19.root";
		ParaSet st("bjtc_pp_allstep_pset");
		//st.setPara<TString>("step2output_folder", "/Users/tabris/frameLite/output/step2_noSeagull/");
		st.setPara<TString>("step2output_folder", "/Users/tabris/frameLite/output/step2/");
		st.setPara<TString>("bjetMC_step2output_name", "Signal_PYTHIA6_bjetSample_allJets");
		st.setPara<TString>("bjetMC_step2input_rg_file", inputfolder+fbjetMC_rg);
		st.setPara<TString>("bjetMC_step2input_gg_file", inputfolder+fbjetMC_gg);
		st.setPara<TString>("bjetMC_step2input_rr_file", inputfolder+fbjetMC_rr);
		st.setPara<TString>("dijetMC_step2input_rg_file", inputfolder+fdijetMC_rg);
		st.setPara<TString>("dijetMC_step2input_gg_file", inputfolder+fdijetMC_gg);
		st.setPara<TString>("dijetMC_step2input_rr_file", inputfolder+fdijetMC_rr);
		st.setPara<TString>("pp5TeVData_step2input_file", inputfolder+fdata);
		st.setPara<TString>("dijetMC_step2output_name", "Signal_PYTHIA6_dijetSample_allJets");
		st.setPara<TString>("pp5TeVData_step2output_name", "Signal_pp5TeVData_allJets");
		st.setPara<bool>("doQA", 1);
		st.setPara<bool>("doSeagullCorr", 0);

		st.setPara<mapper_func>("pad_map", bjtc_pp_config::arrange_pp_bjtc);
		st.setPara<TString (*)(int, int)>("pad_title", bjtc_pp_config::pTtitle);
		st.setPara<int>("pad_nrow", 2);
		st.setPara<int>("pad_ncol", 3);
		return st;
}



