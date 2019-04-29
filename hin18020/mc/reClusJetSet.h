
#include "jetSet.h"

class reClusterJetSet : public jetSet{
		public : reClusterJetSet(const char *name_, const fastjet::JetAlgorithm &jtAlgo, float jtRVal, const fastjet::RecombinationScheme &scheme, std::vector<fastjet::PseudoJet> &p_): 
				 jetSet(name_, jtAlgo, jtRVal, scheme, p_){};
				 void reclusterSetup(const fastjet::JetAlgorithm &jtAlgo, float jtRVal, const fastjet::RecombinationScheme &scheme){
						 jtDef2 = new fastjet::JetDefinition(jtAlgo, jtRVal, scheme);
				 }
				 virtual void makeTree(bool keepCandidate_, float minPtCut_, float maxEtaCut_){
						 keepCandidate = keepCandidate_; minPtCut=minPtCut_; maxEtaCut=maxEtaCut_;
						 jtTree = new TTree(TString(name), "");
						 jtTree->Branch("jt_pt" , &jet_pt );
						 jtTree->Branch("jt_eta", &jet_eta);
						 jtTree->Branch("jt_phi", &jet_phi);
						 jtTree->Branch("rejt_pt" , &jet2_pt );
						 jtTree->Branch("rejt_eta", &jet2_eta);
						 jtTree->Branch("rejt_phi", &jet2_phi);
						 if(keepCandidate){
								 jtTree->Branch("pfcand_indx", &pfcand_indx);
								 jtTree->Branch("pfcand_pt",  &pfcand_pt);
								 jtTree->Branch("pfcand_eta", &pfcand_eta);
								 jtTree->Branch("pfcand_phi", &pfcand_phi);
						 }
				 }
				 virtual void flushSet(){
                         counter = 0;
                         jtTree->Fill();
                         jet_pt ->clear();
                         jet_phi->clear();
                         jet_eta->clear();
                         jet2_pt .clear();
                         jet2_phi.clear();
                         jet2_eta.clear();
                         if(keepCandidate){
                                 pfcand_indx->clear();
                                 pfcand_pt  ->clear();
                                 pfcand_eta ->clear();
                                 pfcand_phi ->clear();
                         }
                 }
				 virtual void fillTree(){
						 cs = new fastjet::ClusterSequence(*particle, *jtDef);
						 std::vector<fastjet::PseudoJet> jets0 = fastjet::sorted_by_pt(cs->inclusive_jets(minPtCut));
						 for(auto &it : jets0) {
								 if(fabs(it.eta())>maxEtaCut) continue;
								 jet_pt->push_back(it.pt());
								 jet_eta->push_back(it.eta());
								 jet_phi->push_back(it.phi_std());
								 std::vector<fastjet::PseudoJet> cands = it.constituents();
								 if(keepCandidate){
										 pfcand_indx->push_back(counter);
										 for(auto & ii: cands){
												 pfcand_pt->push_back (ii.perp());
												 pfcand_eta->push_back(ii.pseudorapidity());
												 pfcand_phi->push_back(ii.phi_std());
												 counter++;
										 }
								 }
								 auto cs2 = new fastjet::ClusterSequence(cands, *jtDef2);
								 std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs2->inclusive_jets(0));
								 if(jets.size()!= 1 ) std::cout<<"Warning: recluster causes jet splitting!!"<<std::endl;
								 //		 std::cout<<jets.size()<<std::endl;
								 for(auto &tt : jets){
										 if(fabs(tt.eta())>maxEtaCut) continue;
										 jet2_pt.push_back(tt.pt());
										 jet2_eta.push_back(tt.eta());
										 jet2_phi.push_back(tt.phi_std());
								 }
								 jets.clear();
								 delete cs2;
						 }
						 flushSet();
						 jets0.clear();
						 delete cs;
				 }

		public : fastjet::JetDefinition *jtDef2;
				 std::vector<float> jet2_pt, jet2_eta, jet2_phi;
};

