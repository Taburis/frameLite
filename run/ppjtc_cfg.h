
#include "ParaSet.h"
#include "edmJtcUtility.h"
#include "TString.h"

ParaSet makePSet_edmJtcDefault(){
// config for step2 and step3
		ParaSet st("edmJtcAnalzyer_pset");
		double jetptbin[] = {110, 120, 136, 152, 168, 184, 200, 216, 232, 248, 264, 280, 296, 312, 328, 344, 360, 380, 400, 432, 500};
		float drbin[] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.};
//		cout<<"size is "<<sizeof(jetptbin)/sizeof(jetptbin[0])<<endl;
		st.setPara<int>("ndr", 14);
		st.setParaVector<float>("drbins", sizeof(drbin)/sizeof(drbin[0]), drbin);
		st.setParaVector<double>("jet_pt_bin", sizeof(jetptbin)/sizeof(jetptbin[0]), jetptbin);
		st.setPara<int>("npt", 6);
		st.setPara<int>("ncent", 1);
		return st;
}

ParaSet makePSet_bjtc_pp_step23(){
		TString inputfolder = "/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/pp_18Dec18/";
		TString fbjetMC_rg = "bJTC_PYTHIA6_RecGen_5TeV_bjetMC_wtaAxis_19Dec18.root";
		TString fbjetMC_gg = "bJTC_PYTHIA6_RecGen_5TeV_bjetMC_wtaAxis_4Jan19.root";
		ParaSet st("bjtc_pp_allstep_pset");
		st.setPara<TString>("bjetMC_step2output_folder", "/Users/tabris/frameLite/run/");
		st.setPara<TString>("bjetMC_step2output_name", "Signal_PYTHIA6_bjetSample_allJets");
		st.setPara<TString>("bjetMC_step2input_rg_file", inputfolder+fbjetMC_rg);
		st.setPara<TString>("bjetMC_step2input_gg_file", inputfolder+fbjetMC_gg);
		return st;
}

jtc_utility::index2d arrange_pp_bjtc(int i, int j){
		jtc_utility::index2d x;
		x.i1 = i/3;
		x.i2 = i%3;
		return x;
}
