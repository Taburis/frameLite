
#include <string>
#include <vector>
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"

class jetSet {
		public : jetSet(const char *name_, const fastjet::JetAlgorithm &jtAlgo, float jtRVal, const fastjet::RecombinationScheme &scheme, std::vector<fastjet::PseudoJet> &p_): name(name_){
						 jtDef = new fastjet::JetDefinition(jtAlgo, jtRVal, scheme);
						 particle = &p_;
						 //						 cs = new fastjet::ClusterSequence(p_, *jtDef);
				 };
				 ~jetSet(){
						 flushSet();
						 delete jtDef;
						 if(cs !=0 ) delete cs;
				 }
				 virtual void makeTree(bool , float, float);
				 virtual void fillTree();
				 virtual void flushSet(){
						 counter = 0;
						 jtTree->Fill();
						 jet_pt ->clear();
						 jet_phi->clear();
						 jet_eta->clear();
						 if(keepCandidate){
								 pfcand_indx->clear();
								 pfcand_pt  ->clear();
								 pfcand_eta ->clear();
								 pfcand_phi ->clear();
						 }
				 }

				 std::string name;
				 float minPtCut, maxEtaCut;
				 TTree* jtTree = 0;
				 std::vector<fastjet::PseudoJet> * particle;
				 bool keepCandidate=0;
				 fastjet::JetDefinition * jtDef;
				 fastjet::ClusterSequence * cs = 0;
				 unsigned int counter = 0;
				 std::vector<float>* jet_pt = new std::vector<float>;
				 std::vector<float>* jet_eta= new std::vector<float>;
				 std::vector<float>* jet_phi= new std::vector<float>;
				 // the pfcand_indx[i] will label the position of the 1st pfcandidation respect to ith jet
				 std::vector<unsigned int>* pfcand_indx = new std::vector<unsigned int>;
				 std::vector<float>* pfcand_pt = new std::vector<float>;
				 std::vector<float>* pfcand_eta= new std::vector<float>;
				 std::vector<float>* pfcand_phi= new std::vector<float>;
};

void jetSet::makeTree(bool keepCandidate_, float minPtCut_, float maxEtaCut_){
		keepCandidate = keepCandidate_; minPtCut=minPtCut_; maxEtaCut=maxEtaCut_;
		jtTree = new TTree(TString(name), "");
		jtTree->Branch("jt_pt" , &jet_pt );
		jtTree->Branch("jt_eta", &jet_eta);
		jtTree->Branch("jt_phi", &jet_phi);
		if(keepCandidate){
				jtTree->Branch("pfcand_indx", &pfcand_indx);
				jtTree->Branch("pfcand_pt",  &pfcand_pt);
				jtTree->Branch("pfcand_eta", &pfcand_eta);
				jtTree->Branch("pfcand_phi", &pfcand_phi);
		}
}

void jetSet::fillTree(){
		cs = new fastjet::ClusterSequence(*particle, *jtDef);
		std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs->inclusive_jets(minPtCut));

		for(unsigned int j =0; j<jets.size(); ++j){
				if(fabs(jets.at(j).eta())>maxEtaCut) continue;
				//				std::cout<<jets.at(j).pt()<<std::endl;
				jet_pt->push_back(jets.at(j).pt());
				jet_eta->push_back(jets.at(j).eta());
				jet_phi->push_back(jets.at(j).phi_std());
				if(keepCandidate){
						pfcand_indx->push_back(counter);
						std::vector<fastjet::PseudoJet> cand = jets.at(j).constituents();
						for(auto & it: cand){
								pfcand_pt->push_back(it.perp());	
								pfcand_eta->push_back(it.pseudorapidity());	
								pfcand_phi->push_back(it.phi_std());	
								counter++;
						}
				}
		}
		flushSet();
		jets.clear();
		delete cs;
}

