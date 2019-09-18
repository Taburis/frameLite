
#include "qaAnalyzer.h"
#include "cfg.h"

void makeWeightConfig(){
	auto qa = new qaAnalyzer();
	ParaSet ps = init_bjetPythia8_cfg();
	qa->loadConfig(ps);
	qa->loadFiles("Pythia8_bjet_mixingTreeScan.root", "Data2015pp_mixingTreeScan.root");
	qa->makePythiaConfig("config_bjetMc.root");
}
