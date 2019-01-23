
#include "TROOT.h"
#include "tdrStyle.h"
#include "ppjtc_cfg.h"
#include "bjtcAnalyzer_Step2.h"
#include "bjtcAnalyzer_Step3.h"

int main(){
		init();
		auto style = set_tdrStyle();
		auto ps = makePSet_edmJtcDefault();
		auto ps_step23 = makePSet_bjtc_pp_step23();
		auto ps_all = ps+ps_step23;
		/*
		bjtcAnalyzer_Step2 step2(ps_all);
		step2.dijetMCwf101();
		step2.bjetMCwf101();
		step2.Datawf101();
		step2.evaluate();
		*/
		bjtcAnalyzer_Step3 step3(ps_all);
		step3.recoJetCheck();
		step3.recoTrackCheck();
		std::cout<<"Done"<<std::endl;
		return 0;
}
