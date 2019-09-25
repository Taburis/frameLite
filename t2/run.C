
#include "ParaSet.h"
#include "preScan_mixingTree.h"
#include "cfg.h"
#include "/home/wang4187/frameLite/crab/lib/configReader.h"

//void run(TString inf = "/mnt/hadoop/store/user/wangx/bJTC_5TeV_pp/Herwig_5TeV_withCSVv2Skim/CRAB_UserFiles/crab_Herwig_5TeV_dijet_forestSkim_finalState_18Sep19/190918_184257/0002/pp_MC_2999.root", TString ofstr = "output.root"){
void run(TString inf = "/mnt/hadoop/store/user/wangx/bJTC_5TeV_pp/PYTHIA8_5TeV_withCSVv2Skim/CRAB_UserFiles/crab_PYTHIA8_5TeV_bjet_forestSkim_finalState_11Sep18/190911_210747/0000/pp_MC_10.root", TString ofstr = "output.root"){
	//ParaSet ps = init_Herwig_cfg();
	ParaSet ps = init_dijetPythia8_cfg();
	//ParaSet ps = init_bjetPythia8_cfg();
	auto frm = new rootEDMFrame();
	auto pm = new preScan_mixingTree();
//	auto cr = new configReader("config_dijetMc.root");
	//auto cr = new configReader("config_bjetMc.root");
//	cr->loadConfig();
//	pm->cr = cr;
	pm->testWeight = 0;
	pm->loadCfg(ps);
	pm->outfstr = ofstr;
	frm->EventMax = -1;
	pm->isdata = 0;
	frm->open(inf);
	frm->addProducer(pm);
	frm->evaluate();
}
