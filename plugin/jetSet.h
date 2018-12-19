

#ifndef POG_JETSET_H
#define POG_JETSET_H

class jetSet {
		public : jetSet(){};
				 ~jetSet(){};
				 jetSet(const char * name, TTree* t,  bool keep){ getJetSet(name, t, keep); };
				 virtual void getJetSet(const char * name_, TTree* t, bool keep){
						 name = name_;
						 keepCandidate = keep;
						 jtTree = t;
						 jtTree->SetBranchAddress("jt_pt" , &jt_pt);
						 jtTree->SetBranchAddress("jt_eta", &jt_eta);
						 jtTree->SetBranchAddress("jt_phi", &jt_phi);
						 if(keepCandidate){
								 jtTree->SetBranchAddress("pfcand_indx" , &pfcand_indx);
								 jtTree->SetBranchAddress("pfcand_pt" , &pfcand_pt);
								 jtTree->SetBranchAddress("pfcand_eta", &pfcand_eta);
								 jtTree->SetBranchAddress("pfcand_phi", &pfcand_phi);
						 }
				 }
				 virtual void for_each_jet(unsigned int ){};
				 void loopJets(){
						 for(unsigned int j=0; j<jt_pt->size(); ++j){
								 for_each_jet(j);
						 } 
				 }

				 TString name;
				 bool keepCandidate = 0;
				 TTree* jtTree = 0;
				 std::vector<float>* jt_pt  =0;
				 std::vector<float>* jt_eta =0;
				 std::vector<float>* jt_phi =0;
				 std::vector<unsigned int>* pfcand_indx =0;
				 std::vector<float>* pfcand_pt =0;
				 std::vector<float>* pfcand_eta =0;
				 std::vector<float>* pfcand_phi =0;
};

#endif
