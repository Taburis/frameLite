
#include <algorithm>
#include <unordered_map>
#include <string>
#include "rootEDM.h"
#include "histManager.h"
#include "TMath.h"
#include "TFile.h"
#include "TString.h"
#include "jetSet.h"

struct track{
		TTree* tree=0;
		int highPurity;
		std::vector<float> *pt  =0 ;
		std::vector<float> *ptErr=0 ;
		std::vector<float> *eta =0 ;
		std::vector<float> *phi =0 ;
		std::vector<float> *chi2=0 ;
		std::vector<float> *dz=0 ;
		std::vector<float> *dzErr=0 ;
		std::vector<float> *dxy=0 ;
		std::vector<float> *dxyErr=0 ;
		std::vector<int>   *ndof=0 ;
		std::vector<int>  *nlayer=0 ;
		std::vector<int>  *nHit=0 ;
		std::vector<float>  *pfEcal=0 ;
		std::vector<float>  *pfHcal=0 ;
};

class edmMCBase : public rootEDMProducer{
		public : edmMCBase(){};
				 virtual ~edmMCBase(){
						 // std::cout<<"deleting edmMCBase"<<std::endl;
						 // for(auto & it : js){
						 //		 delete it.second;
						 // }
				 }
				 virtual void handleGenParticle(const char * name);
				 virtual void handleRecoTrack (const char * name);
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
				 track trk;
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

void edmMCBase::handleRecoTrack(const char * name){
	trk.tree = handle(name);
	trk.tree->SetBranchAddress("highPurity", &trk.highPurity);
	trk.tree->SetBranchAddress("trkPt", &trk.pt);
	trk.tree->SetBranchAddress("trkPtError", &trk.ptErr);
	trk.tree->SetBranchAddress("trkEta", &trk.eta);
	trk.tree->SetBranchAddress("trkPhi", &trk.phi);
	trk.tree->SetBranchAddress("trkDxy", &trk.dxy);
	trk.tree->SetBranchAddress("trkDz", &trk.dz);
	trk.tree->SetBranchAddress("trkDxyError1", &trk.dxyErr);
	trk.tree->SetBranchAddress("trkDzError1", &trk.dzErr);
	trk.tree->SetBranchAddress("trkNdof", &trk.ndof);
	trk.tree->SetBranchAddress("trkNHit", &trk.nHit);
	trk.tree->SetBranchAddress("trkNlayer", &trk.nlayer);
	trk.tree->SetBranchAddress("trkChi2", &trk.chi2);
	trk.tree->SetBranchAddress("pfEcal", &trk.pfEcal);
	trk.tree->SetBranchAddress("pfHcal", &trk.pfHcal);
}


