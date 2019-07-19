
#include "TROOT.h"
#include "tdrStyle.h"
//#include "ppjtc_Esch_cfg.h"
#include "ppjtc_cfg.h"
#include "bjtcAnalyzer_Step2.h"
#include "bjtcAnalyzer_Step3.h"

int main(){
		init();
		auto style = set_tdrStyle();
		auto ps = makePSet_edmJtcDefault();
		auto ps_step23 = makePSet_bjtc_pp_step23();
		auto ps_all = ps+ps_step23;
		bjtcAnalyzer_Step2 step2(ps_all);
		//step2.inclCalo_Wf001();
		//step2.dijetMCwf101();
		//step2.Datawf101();
		//step2.bjetMCwf101();
		//step2.jerMCwf101();
		//step2.PYTHIAwf101();
		//step2.Syswf101();
		step2.jesMCwf101();
		step2.evaluate();
		/*
		bjtcAnalyzer_Step3 step3(ps_all);
//		step3.fordijet_recoTrackCheck();
//		step3.get_dijtc_correction();
//		step3.apply_dijtc_correction();
//		step3.get_bjtc_correction();
		step3.produce_bjtc();
		step3.overlay_injtc_vs_bjtc();
		std::cout<<"Done"<<std::endl;
		*/
		return 0;
}
