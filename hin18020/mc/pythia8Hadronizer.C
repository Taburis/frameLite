
#include "Pythia8/Pythia.h"
#include <iostream>
#include "genParicle.h"
#include "evtInfo.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "reClusJetSet.h"

using namespace Pythia8;

int main(int argc, char * argv[]){

		std::string input;
		TString	output;
		//		std::cout<<argc<<std::endl; 
		if(argc != 3 ) {
				std::cout << "Usage: hadronizer.o input.lhe output.root" <<std::endl;
				return 0;
		}else {
				input = argv[1];
				output = argv[2];	
		}

		std::string istr = "Beams:LHEF = "+input;
		int nAbort = 10;
		Pythia pythia;
		// Initialize Les Houches Event File run.
		pythia.readString(istr);	
		pythia.readString("Beams:frameType = 4");	

		//pythia.readString("HardQCD:all = on");
		//pythia.readString("Random:setSeed = on");
		//pythia.readString("Random:seed = 0");

		pythia.readString("Tune:preferLHAPDF = 2");
		pythia.readString("Main:timesAllowErrors = 10000");
		pythia.readString("Check:epTolErr = 0.01");
		pythia.readString("Beams:setProductionScalesFromLHEF = off");
		pythia.readString("SLHA:keepSM = on");
		pythia.readString("SLHA:minMassSM = 1000.");
		pythia.readString("ParticleDecays:limitTau0 = on");
		pythia.readString("ParticleDecays:tau0Max = 10");
		pythia.readString("ParticleDecays:allowPhotonRadiation = on");

		pythia.readString("Tune:pp 14");
		pythia.readString("Tune:ee 7");
		pythia.readString("MultipartonInteractions:pT0Ref=2.4024");
		pythia.readString("MultipartonInteractions:ecmPow=0.25208");
		pythia.readString("MultipartonInteractions:expPow=1.6");

		pythia.readString("POWHEG:veto = 1");
		pythia.readString("POWHEG:pTdef = 1");
		pythia.readString("POWHEG:emitted = 0");
		pythia.readString("POWHEG:pTemt = 0");
		pythia.readString("POWHEG:pThard = 0");
		pythia.readString("POWHEG:vetoCount = 100");
		pythia.readString("SpaceShower:pTmaxMatch = 2");
		pythia.readString("TimeShower:pTmaxMatch = 2");

		pythia.init();
		//std::cout<<std::ios_base::failbit;
		//std::cout.clear();

		int iAbort = 0;
		long nevt = 0;

		auto wf = TFile::Open(output, "recreate");
		auto t = new TTree("mixing_tree", "");
		genParticle gp;
		gp.keepParton = 1;
		evtInfo ei;
		gp.init(t);

		std::vector<fastjet::PseudoJet> pfcands;
		reClusterJetSet* rjset = new reClusterJetSet("ak4ES_WTAJet", fastjet::antikt_algorithm, 0.4, fastjet::E_scheme, pfcands);
		rjset->reclusterSetup(fastjet::antikt_algorithm, 1, fastjet::WTA_pt_scheme);
		rjset->makeTree(0, 80, 1.6);




		for (int iEvent = 0; ; ++iEvent) {
				if( !pythia.next()) {
						if(pythia.info.atEndOfFile()) break;
						if (++iAbort < nAbort) continue;
						break;
				}	
				nevt++;
				//if(nevt > 1000) break;

				gp.addEvtInfo(pythia.info.weight(), pythia.info.pTHat());
				for (int i = 0; i < pythia.event.size(); ++i){
						if(!pythia.event[i].isFinal()){
								if(fabs(pythia.event[i].id()) < 6 || fabs(pythia.event[i].id()) == 21){
										gp.addParton(pythia.event[i].pT(), pythia.event[i].eta(), pythia.event[i].phi(), pythia.event[i].id(), pythia.event[i].status());
								}
						//		continue;
						} else {

								if(pythia.event[i].pT() < 1) continue;
								if(fabs(pythia.event[i].eta()) > 2.4) continue;
								//if(fabs(pythia.event[i].charge()) == 0) continue;

								if(fabs(pythia.event[i].id()) == 12) continue;
								if(fabs(pythia.event[i].id()) == 14) continue;
								if(fabs(pythia.event[i].id()) == 16) continue;
								gp.addParticle(pythia.event[i].pT(), pythia.event[i].eta(), pythia.event[i].phi(), pythia.event[i].charge());
								pfcands.push_back(fastjet::PseudoJet(pythia.event[i].px(), pythia.event[i].py(),pythia.event[i].pz(),pythia.event[i].e()));
						}
				}
				t->Fill();
				rjset->fillTree();

				pfcands.clear();
				gp.clear();
		}

		wf->cd();
		t->Write();
		rjset->jtTree->Write();
		wf->Close();

		return 0;
}
