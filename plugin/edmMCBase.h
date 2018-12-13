
#include <algorithm>
#include <unordered_map>
#include <string>
#include "rootEDM.h"
#include "histManager.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "jetSet.h"


class edmMCBase : public rootEDMProducer{
		public : edmMCBase(){};
				 virtual ~edmMCBase(){
						// std::cout<<"deleting edmMCBase"<<std::endl;
						// for(auto & it : js){
						//		 delete it.second;
						// }
				 }
				 virtual void handleGenParticle(const char * name);
				 jetSet* handleJetSet(const char* name, bool keep = 0){
						 auto tt= handle(name);
						 _js[std::string(name)] = new jetSet(name, tt, keep);
						 return _js[std::string(name)];
				 }
		public :
				 TTree* gen_particle_t;
				 float weight =0;
				 float pthat =0;
				 std::vector<float> *pt=0;
				 std::vector<float> *eta=0;
				 std::vector<float> *phi=0;
				 std::vector<int >  *chg=0;
				 std::vector<int >  *pdg=0;
				 std::unordered_map<std::string, jetSet*> _js;

};

void edmMCBase::handleGenParticle(const char * name){
		gen_particle_t = handle(name);
		gen_particle_t->SetBranchAddress("pthat", &pthat);
		gen_particle_t->SetBranchAddress("weight",&weight);
		gen_particle_t->SetBranchAddress("pt" , &pt);
		gen_particle_t->SetBranchAddress("eta", &eta);
		gen_particle_t->SetBranchAddress("phi", &phi);
		gen_particle_t->SetBranchAddress("chg", &chg);
		gen_particle_t->SetBranchAddress("pdg", &pdg);
}
