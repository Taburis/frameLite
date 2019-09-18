
#include "ParaSet.h"

namespace bjtc_Pythia8_cfg{
	const Float_t jetptbin[] = {110, 120, 136, 152, 168, 184, 200, 216, 232, 248, 264, 280, 296, 312, 328, 344, 360, 380, 400, 425, 450, 475, 500};
	int njtpTbin = 22;
	const Float_t pThatInterval[] = {80, 100, 120, 170, 220, 500};
	int npThat = 5;
	const Double_t xsec[] = {4.919E-04, 1.73E-04, 7.096E-05, 1.265E-05, 3.031E-06};

	bool evtCut(float vz){
		if(fabs(vz)> 15) return 1;
		return 0;
	}
}
//DBCode:
// 0. Pythia8 dijet
// 1. Pythia8 bjete
// 2. Pythia6 bjete
// 3. Pythia6 bjete
// 4. Herwig  dijet 
// 5. Data

ParaSet init_bjetPythia8_cfg(){
	ParaSet st("bjtc_cfg");	
	st.setPara<int>("DBCode", 1);
	st.setPara<int>("njtpTbin", bjtc_Pythia8_cfg::njtpTbin);
	st.setPara<int>("npThat", bjtc_Pythia8_cfg::npThat);
	st.setPara<const Float_t*>("jtpTbin", bjtc_Pythia8_cfg::jetptbin);
	st.setPara<const Float_t*>("pThatInterval", bjtc_Pythia8_cfg::pThatInterval);
	st.setPara<const Double_t*>("xsec", bjtc_Pythia8_cfg::xsec);
	st.setPara<bool (*)(float)>("evtCut", &(bjtc_Pythia8_cfg::evtCut));
	return st;
}

ParaSet init_Herwig_cfg(){
	ParaSet st("bjtc_cfg");	
	st.setPara<int>("DBCode", 4);
	st.setPara<int>("njtpTbin", bjtc_Pythia8_cfg::njtpTbin);
	st.setPara<int>("npThat", bjtc_Pythia8_cfg::npThat);
	st.setPara<const Float_t*>("jtpTbin", bjtc_Pythia8_cfg::jetptbin);
	st.setPara<bool (*)(float)>("evtCut", &(bjtc_Pythia8_cfg::evtCut));
	return st;
}
