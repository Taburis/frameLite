
#include "tdrStyle.h"
#include "ppjtc_cfg.h"
#include "bjtcAnalyzer_Step2.h"

int main(){
		auto style = set_tdrStyle();
		auto ps = makePSet_edmJtcDefault();
		auto ps_step23 = makePSet_bjtc_pp_step23();
		auto ps_all = ps+ps_step23;
		bjtcAnalyzer_Step2 step2(ps_all);
		step2.dijetMCwf101();
		step2.bjetMCwf101();
		step2.Datawf101();
		step2.evaluate();
		std::cout<<"Done"<<std::endl;
		return 0;
}
