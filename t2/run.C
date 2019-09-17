
#include "ParaSet.h"
#include "preScan_mixingTree.h"
#include "cfg.h"

void run(TString inf = "/mnt/hadoop/store/user/wangx/bJTC_5TeV_pp/PYTHIA8_5TeV_withCSVv2Skim/CRAB_UserFiles/crab_PYTHIA8_5TeV_bjet_forestSkim_finalState_11Sep18/190911_210747/0000/pp_MC_10.root", TString ofstr = "output.root"){
	ParaSet ps = init_bjetPythia8_cfg();
	auto frm = new rootEDMFrame();
	auto pm = new preScan_mixingTree();
	pm->loadCfg(ps);
	pm->outfstr = ofstr;
	frm->EventMax = -1;
	frm->open(inf);
	frm->addProducer(pm);
	frm->evaluate();
}
