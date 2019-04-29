
#ifndef GENPARTICLE_H
#define GENPARTICLE_H
#include "TTree.h"

class genParticle {
		public :
				genParticle(){};
				void init(TTree *t);
				void clear();
				void addParticle(float pt, float eta, float phi, int chg);
				void addEvtInfo(double iweight, float ipthat){weight = iweight; ipthat = pthat;};
				void addParton(float pt, float eta, float phi, int id){
						parton_pt.push_back(pt);
						parton_eta.push_back(eta);
						parton_phi.push_back(phi);
						parton_id.push_back(id);
				}

				TTree *ot;	
				double weight;
			   	float	pthat;
				std::vector<float> ppt, peta, pphi;
				std::vector<int> pchg, parton_id;
				std::vector<float> parton_pt, parton_eta, parton_phi;
				bool keepParton = 0;
};

void genParticle::init(TTree *t){
		ot = t;
		ot->Branch("pthat", &pthat, "pthat/F");
		ot->Branch("weight", &weight, "weight/D");
		ot->Branch("gen_pt" , &ppt);
		ot->Branch("gen_eta", &peta);
		ot->Branch("gen_phi", &pphi);
		ot->Branch("gen_chg", &pchg);
		if(keepParton){
				ot->Branch("parton_pt" ,&parton_pt);
				ot->Branch("parton_eta",&parton_eta);
				ot->Branch("parton_phi",&parton_phi);
				ot->Branch("parton_id", &parton_id);
		}
}

void genParticle::clear(){
		ppt.clear();
		peta.clear();
		pphi.clear();
		pchg.clear();
		if(keepParton){
				parton_pt.clear();
				parton_eta.clear();
				parton_phi.clear();
				parton_id.clear();
		}
}

void genParticle::addParticle(float pt, float eta, float phi, int chg){
		ppt.push_back(pt);
		peta.push_back(eta);
		pphi.push_back(phi);
		pchg.push_back(chg);
}

#endif
