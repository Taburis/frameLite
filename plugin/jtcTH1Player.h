
#ifndef JTCTH1PLAYER_H
#define JTCTH1PLAYER_H
#include "edmJtcUtility.h"
#include "matrixTH1Ptr.h"
#include "TString.h"
#include "jtcError.h"

enum jetType {incl, tgb, ttb, trb};	
enum dataSet {data, bjet_rr, bjet_gg, bjet_rg, bjet_gr, dijet_gg, dijet_rg, dijet_rr};
class jtcTH1Player : public matrixTH1Ptr{
		// to use this class, a configure file needs to be loaded in g_cfg;
		public : jtcTH1Player (): matrixTH1Ptr(){};
				 virtual ~jtcTH1Player(){};
				 jtcTH1Player (const char *_name, int n = 6, int m=1): matrixTH1Ptr(_name, n,m) {}
				 jtcTH1Player (const char *_name, jtcTH1Player &j2): matrixTH1Ptr(_name, j2) {}
				 jtcTH1Player (jtcTH1Player &jtc): matrixTH1Ptr(jtc){}
				 jtcTH1Player * phiSideband(float x, float y, TString opt = "");
				 void add_frac_error(float frac);
				 void ring_corr(matrixTH1Ptr * m2, float drmax = -1, bool check = 0);
				 void ringCorr();
				 jtcTH1Player* bkgSub(const char * name, float side1 = 1.5, float side2 = 2.5);
				 jtcTH1Player* getSideMixTable(const char *name, float sidemin, float sidemax);
				 jtcTH1Player* drIntegral(const char *name, int ndr, float *drbins);
				 jtcTH1Player* projX(const char * name, float a0 , float a1 , TString opt="e", bool dorebin = 1);
				 jtcTH1Player* projY(const char * name, float a0 , float a1 , TString opt="e", bool dorebin = 1);
				 jtcTH1Player* getSignal_ME_based(const char *name, float sidemin, float sidemax, bool );
				 jtcTH1Player* contractX(const char *name);
				 jtcTH1Player* rotate2D(const char* name);
				 jtcTH1Player* smoothMixing(const char* name);
				 void doChi2Test(jtcTH1Player *, Option_t* opt);
				 void doGeoCorr(jtcTH1Player* j2);
				 void invariant();
				 void setDrError();
				 void setDrError(jtcTH1Player *jtc);
				 jtcTH1Player* getBkgError();
				 void drawBkgErrorCheck(int i, int j);
				 float bkgError(int i, int j){return m2error->at(i,j)->bg_err;}
				 float meError(int i, int j){return m2error->at(i,j)->me_err;}
				 float bkgLevel(int i, int j){return m2error->at(i,j)->bkg_level;}
				 bool loadBkgError(jtcTH1Player * );
				 void rebinX(int n){
						 for(int j=0; j<Ncol(); ++j){ for(int i=0; i<Nrow(); i++){
								 at(i,j)->RebinX(n);
						 }}
				 }
				 void set2DConst(float x);
				 void absError(float x);

		public : bool kGotError = 0;
				 matrixPtrHolder<jtcError>* m2error = nullptr;
};

jtcTH1Player* jtcTH1Player::phiSideband(float x, float y, TString opt){
		auto m2 = new jtcTH1Player(Form("phisb_%s",this->getName()), nrow, ncol);
		for(auto i=0; i<ref.size() ; ++i) {
				auto h = jtc_utility::projX(1, (TH2D*)ref[i], x, y, opt);
				(m2->ref).at(i) = (TH1*)h;
		}
		return m2;
}

void jtcTH1Player::ring_corr(matrixTH1Ptr * m2, float drmax, bool errCheck){
		for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
				for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
						auto h = at(i,j);
						auto corr = m2->at(i,j);
						for(int k=1; k<h->GetNbinsX()+1; ++k){
								for(int l=1; l<h->GetNbinsY()+1; ++l){
										float deta = h->GetXaxis()->GetBinCenter(k);
										float dphi = h->GetYaxis()->GetBinCenter(l);
										float dr = pow(deta*deta+dphi*dphi, 0.5);
										if(dr > drmax && drmax > 0) continue;
										int nn = corr->GetXaxis()->FindBin(dr);
										if(errCheck){
												float error = corr->GetBinError(nn);
												if( fabs(corr->GetBinContent(nn)-1)<error) continue; 
										}
										h->SetBinContent(k,l, h->GetBinContent(k,l)/corr->GetBinContent(nn));
								}
						}
				}
		}
		return ;
}

jtcTH1Player * jtcTH1Player::bkgSub(const char * name, float side1, float side2){
		// this function will return the bkg
		jtcTH1Player * bkg2 = new jtcTH1Player(Form("bkg_%s",this->getName()), Nrow(), Ncol());
		jtcTH1Player * res = (jtcTH1Player*) this->clone(name);
		for(int j=0; j<matrixTObjPtr<TH1>::ncol; ++j){
				for(int i=0; i<matrixTObjPtr<TH1>::nrow; i++){
						auto bkg = (TH2D*) jtc_utility::getV2Bkg((TH2D*)this->at(i,j), side1, side2);
						res->at(i,j)->Add(res->at(i,j), bkg, 1, -1);
						bkg2->add(bkg, i, j);
				}
		}
		bkg2->doFree = 1;
		delete bkg2;
		return res;
}

jtcTH1Player * jtcTH1Player::projX(const char * name, float a0 , float a1 , TString opt, bool dorebin){
		jtcTH1Player *m2p = new jtcTH1Player(name, *this);
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						auto h = jtc_utility::projX(dorebin, (TH2D*)this->at(i,j), a0, a1, opt);
						m2p->add(h, i, j);
				}
		}
		m2p->setName(name);
		return m2p;
}

jtcTH1Player * jtcTH1Player::projY(const char * name, float a0 , float a1 , TString opt, bool dorebin){
		jtcTH1Player *m2p = new jtcTH1Player(name, *this);
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						auto h = jtc_utility::projY(dorebin, (TH2D*)this->at(i,j), a0, a1, opt);
						m2p->add(h, i, j);
				}
		}
		m2p->setName(name);
		return m2p;
}

jtcTH1Player* jtcTH1Player::getBkgError(){
		auto j2 = projX(Form("proX_%s", getName()), -1, 0.99);
		m2error = new matrixPtrHolder<jtcError>(j2->Nrow(), j2->Ncol());
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						auto er = new jtcError();
						er->getBkgError((TH1D*)j2->at(i,j));
						m2error->add(er, i,j);
				}
		}
		j2->m2error = m2error;
		kGotError = 1;
		return j2;
}

jtcTH1Player* jtcTH1Player::drIntegral(const char* name, int ndr, float *drbins){
		auto j2dr = new jtcTH1Player(name, *this);
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						auto hdr = jtc_utility::doDrIntegral(at(i,j)->GetName(), (TH2D*) at(i,j), ndr, drbins);
						j2dr->add(hdr, i, j);
				}
		}
		j2dr->doFree = 1;
		j2dr->setName(name);
		return j2dr;
}


void jtcTH1Player::setDrError(){
		if( !kGotError ){ 
				std::cout<<"please load or get the error info. for "<<getName()<<" first!"<<std::endl;
				return;
		}
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						m2error->at(i,j)->makeDrError(this->at(i,j));
				}
		}
		return;
}

void jtcTH1Player::add_frac_error(float frac){
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						jtc_utility::add_frac_error((TH1*)this->at(i,j), frac);
				}
		}
		return;
}

void jtcTH1Player::setDrError(jtcTH1Player *j2){
		loadBkgError(j2);
		setDrError();
}

bool jtcTH1Player::loadBkgError(jtcTH1Player * j2){
		m2error = j2->m2error;
		kGotError = j2->kGotError;
		return kGotError;
}

jtcTH1Player* jtcTH1Player::getSideMixTable(const char *name, float sidemin, float sidemax){
		auto me = new jtcTH1Player(name, *this);
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						auto hme = jtc_utility::sideBandMixingTableMaker((TH2D*)this->at(i,j), sidemin, sidemax);
						hme->SetName(Form("%s_%d_%d", name, i, j));
						me->add(hme, i, j);
				}
		}
		me->doFree = 1;
		return me;
}
void jtcTH1Player::invariant(){
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						jtc_utility::invariant_TH2((TH2*)this->at(i,j));
				}
		}
}

jtcTH1Player* jtcTH1Player::getSignal_ME_based(const char *name, float sidemin, float sidemax, bool doBkgSub){
		auto me = this->getSideMixTable(Form("%s_sideME", name), sidemin, sidemax);
		auto step1 = (jtcTH1Player*)((*this)/(*me));
		jtcTH1Player* res;
		if(doBkgSub){
				auto step2 = step1->bkgSub(name, 1.5, 2.5);
				res = step2;
				delete step1;
				delete me;
		} else {
				res = step1; 
				res->setName(name);
				delete me;
		}
		return res;
}

void jtcTH1Player::doGeoCorr(jtcTH1Player* j2){
		auto corr = jtc_utility::drGeoTest((TH2D*) j2->at(0,0), (TH1D*) this->at(0,0));
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						at(i,j)->Divide(corr);
				}
		}
}

jtcTH1Player* jtcTH1Player::contractX(const char *name){
		jtcTH1Player * jc = new jtcTH1Player(name, 1 , Ncol());
		for(int j=0; j<Ncol(); ++j){
				auto h = (TH1*) at(0,j)->Clone(Form("%s_0_%d", name, j));
				jc->add(h, 0, j);
				for(int i=1; i<Nrow(); i++){
						//cout<<i<<" : "<<at(i,j)->GetName()<<endl;
						h->Add(at(i,j));
				}
		}
		return jc;
}

void jtcTH1Player::doChi2Test(jtcTH1Player *j2, Option_t* opt){
		Double_t *par = 0;
		m2error = new matrixPtrHolder<jtcError>(j2->Nrow(), j2->Ncol());
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						auto er = new jtcError();
						m2error->add(er, i,j);
						m2error->at(i,j)->chi2 = at(i,j)->Chi2Test(j2->at(i,j), opt, par);
						std::cout<<"for ("<<i<<", "<<j<<") p value: "<<m2error->at(i,j)->chi2<<std::endl;
				}
		}
		return ;
}

void jtcTH1Player::set2DConst(float x){
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						jtc_utility::set2DConst(x,(TH1*)at(i,j));
				}
		}
}

void jtcTH1Player::absError(float x){
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						auto h = at(i,j);
						for(int k=1; k<h->GetNbinsX()+1; k++){
								h->SetBinError(k, x);
						}
				}
		}
}

jtcTH1Player* jtcTH1Player::rotate2D(const char* name){
		auto m2new = new jtcTH1Player(name, this->Nrow(), this->Ncol());
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						auto h = jtc_utility::rotate2D(Form("%s_%d_%d", name, i,j), (TH2D*)at(i,j));
						m2new->add(h, i,j);
				}
		}
		return m2new;
}

jtcTH1Player* jtcTH1Player::smoothMixing(const char* name){
		auto m2new = new jtcTH1Player(name, this->Nrow(), this->Ncol());
		for(int j=0; j<Ncol(); ++j){
				for(int i=0; i<Nrow(); i++){
						m2new->add(jtc_utility::mixingTableMaker((TH2D*)at(i,j), 1), i,j);
				}
		}
		return m2new;
}

#endif
