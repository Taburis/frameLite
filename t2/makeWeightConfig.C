
#include "qaAnalyzer.h"
#include "cfg.h"

void makeWeightConfig(){
	auto qa = new qaAnalyzer();
	ParaSet ps = init_dijetPythia8_cfg();
	//ParaSet ps = init_bjetPythia8_cfg();
	qa->loadConfig(ps);
	cout<<ps.getPara<int>("npThat")<<endl;
	qa->loadFiles("Pythia8_dijet_mixingTreeScan.root", "Data2015pp_mixingTreeScan.root");
	qa->makePythiaConfig("config_dijetMc.root");
	//qa->loadFiles("Pythia8_bjet_mixingTreeScan.root", "Data2015pp_mixingTreeScan.root");
	//qa->makePythiaConfig("config_bjetMc.root");
}
