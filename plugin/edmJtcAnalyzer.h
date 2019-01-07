
#ifndef EDMJTCANALYZER_H
#define EDMJTCANALYZER_H
#include "rootEDM.h"
#include "edmJtcSignalProducer.h"
#include "matrixTObjPtr.h"
#include "TString.h"

class edmJtcAnalyzer : public rootEDMAnalyzer{
		/*
		   the object to produce the signal from the correlation histgorams, i.e. correct the acceptance, remove bkg, for one set of jtc
		   */
		public : edmJtcAnalyzer(): rootEDMAnalyzer(){
						 _nodeClassName_ = "edmJtcAanlyzer";
				 }
				 virtual ~edmJtcAnalyzer(){
						 m2producer.cleanAll();
						 delete [] drbins;
				 }
				 edmJtcAnalyzer(ParaSet &ps);
				 void open(const char* path, bool mix_same = 1){
						 sigf = TFile::Open(path);
						 if(mix_same) mixf=sigf;
				 }
				 virtual int analyze();
				 void initialization();
				 void setupProducer(edmJtcSignalProducer *p);
				 void check_phi_sideband(TString name);

				 matrixTObjPtr<TH2D> m2sig;
				 matrixTObjPtr<TH2D> m2mix;
				 matrixPtrHolder<edmJtcSignalProducer> m2producer;
				 int npt, ncent, ndrbin;
				 float *drbins;
				 TString sig_name, mix_name; // the prefix name for signal histograms
				 TFile *sigf=0, *mixf=0;
				 jtc_utility::index2d (*pad_map)(int, int);
				 bool saveFile =0 ; TString output, outputOpt = "update";
};

edmJtcAnalyzer::edmJtcAnalyzer(ParaSet &ps): rootEDMAnalyzer()
{
		npt = ps.getPara<int>("npt");
		ncent = ps.getPara<int>("ncent");
		ndrbin = ps.getPara<int>("ndr");
		drbins = ps.getVectorAsArray<float>("drbins");
}

void edmJtcAnalyzer::initialization(){
		m2sig.setup(sig_name, npt, ncent); 
		m2mix.setup(mix_name, npt, ncent); 
		m2producer.setup(npt, ncent); 
}

void edmJtcAnalyzer::setupProducer(edmJtcSignalProducer *p){
		// setup the producer for getting signal
		p->ndrbin = ndrbin;
		p->drbins = drbins;
		p->setFitFunction();
}

int edmJtcAnalyzer::analyze(){

		std::cout<<"staring "<<_nodeClassName_<<"..."<<std::endl;
		initialization();
		if(sigf==0) {
				std::cout<<"no signal file specified yet, abort!"<<std::endl;
				return 1;
		}
		m2sig.autoLoad(sigf);
		m2mix.autoLoad(sigf);
		for(int i=0; i<npt; i++){
				for(int j=0; j<ncent; j++){
						auto ptr = new edmJtcSignalProducer();
						m2producer.add(ptr, i,j);
						m2producer(i,j)->addSig(m2sig(i,j), m2mix(i,j));
						setupProducer(m2producer(i,j));
						m2producer(i,j)->produce();
				}
		}
		if(saveFile){
				TFile wf(output, outputOpt);
				wf.cd();
				for(int i=0; i<npt; i++){
						for(int j=0; j<ncent; j++){
								m2producer(i,j)->write();
						}
				}
				wf.Close();
		}
		return 0;
}

void edmJtcAnalyzer::check_phi_sideband(TString cname){
		multi_canvas<TH1D> canvas(TString(m2sig.name+"_c1"), TString(m2sig.name), 2, 3);
		//canvas.CD(0,0);
		//m2sig(0,0)->ProjectionX()->Draw();
		for(int i=0; i<m2sig.nrow; ++i){
				for(int j=0; j<m2sig.ncol; ++j){
						auto indx = pad_map(i,j);
						canvas.CD(indx.i1, indx.i2);
						//				std::cout<<indx.i1<<", "<<indx.i2<<std::endl;
						m2producer(i,j)->check_seagull();
						//				m2sig(i,j)->ProjectionX()->Draw();
				}
		}
		canvas.SaveAs(cname);
}
#endif
