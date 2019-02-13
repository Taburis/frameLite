
#ifndef EDMJTCANALYZER_H
#define EDMJTCANALYZER_H
#include "analyzerIOServer.h"
#include "rootEDM.h"
#include "edmJtcSignalProducer.h"
#include "matrixTH1Ptr.h"
//#include "matrixTObjPtr.h"
#include "TString.h"
#include "jtcQaMonitor.h"

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
				 matrixTH1Ptr* m2_ring_corr(TString name, matrixTH1Ptr* m2, matrixTH1Ptr* corr){
						 auto m2return = m2->clone(name);
						 for(int i=0; i<m2->Nrow(); i++){
								 for(int j=0; j<m2->Ncol(); j++){
										 jtc_utility::ring_corr((TH2D*) m2return->at(i,j), (TH1D*) corr->at(i,j));
								 }
						 }
						 return m2return;
				 }
				 matrixTH1Ptr* m2_bkgSub(TString name, matrixTH1Ptr* m2){
						 auto m2sig = new matrixTH1Ptr(name, npt, ncent);
						 for(int i=0; i<m2->Nrow(); i++){
								 for(int j=0; j<m2->Ncol(); j++){
										 auto h = m2producer(i,j)->bkgSub((TH2D*)(m2->at(i,j)), name+Form("_%d_%d", i,j));
										 m2sig->add(h, i, j);
								 }
						 }
						 return m2sig;
				 }
				 void add2Step3(matrixTH1Ptr* m2){
						 for(int i=0; i<m2->Nrow(); i++){
								 for(int j=0; j<m2->Ncol(); j++){
										 m2producer(i,j)->sig_step3 = (TH2D*) m2->at(i,j);
										 m2producer(i,j)->name = m2->at(i,j)->GetName();
								 }
						 }
				 }
				 void runProduce(){
						 for(int i=0; i<npt; i++){
								 for(int j=0; j<ncent; j++){
										 m2producer(i,j)->produce();
								 }
						 }
				 }
				 void write(){
						 for(int i=0; i<npt; i++){
								 for(int j=0; j<ncent; j++){
										 m2producer(i,j)->write();
								 }
						 }
				 }

				 virtual int analyze();
				 void initialization();
				 void setupProducer(edmJtcSignalProducer *p);
				 void check_phi_sideband(TString name);
				 void signal_deta_qa();

				 matrixTObjPtr<TH2D> m2sig;
				 matrixTObjPtr<TH2D> m2mix;
				 matrixPtrHolder<edmJtcSignalProducer> m2producer;
				 int npt, ncent, ndrbin;
				 float *drbins;
				 float norm; // the normalization factor comes from integral of jet pt
				 TString sig_name, mix_name, norm_name, dataset_name; // the prefix name for signal histograms
				 TFile *sigf=0, *mixf=0;
				 jtc_utility::index2d (*pad_map)(int, int);
				 bool saveFile =0 ; TString output, outputOpt = "update";
				 bool doSeagullCorr = 0;
				 bool doQA = 0;
				 ParaSet * ps; 
};

edmJtcAnalyzer::edmJtcAnalyzer(ParaSet &ps_): rootEDMAnalyzer()
{
		npt = ps_.getPara<int>("npt");
		ncent = ps_.getPara<int>("ncent");
		ndrbin = ps_.getPara<int>("ndr");
		drbins = ps_.getVectorAsArray<float>("drbins");
		ps = &ps_;
}

void edmJtcAnalyzer::initialization(){
		m2producer.setup(npt, ncent); 
		for(int i=0; i<npt; i++){
				for(int j=0; j<ncent; j++){
						auto ptr = new edmJtcSignalProducer();
						setupProducer(ptr);
						m2producer.add(ptr, i,j);
				}
		}
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
		m2sig.setup(sig_name, npt, ncent); 
		m2mix.setup(mix_name, npt, ncent); 
		m2sig.autoLoad(sigf);
		m2mix.autoLoad(sigf);
		multi_canvas<TH1> *mm=nullptr;
		int nrow = ps->getPara<int>("pad_nrow");
		int ncol = ps->getPara<int>("pad_ncol");
		mapper_func fmap;
		if(doSeagullCorr){
				mm = new multi_canvas<TH1>("sg_"+output, "" , nrow, ncol);
				fmap  = ps->getPara<mapper_func>("pad_map");
		}
		for(int i=0; i<npt; i++){
				for(int j=0; j<ncent; j++){
						TString name = norm_name+Form("_%d", j);
						auto hjtpt = (TH1D*)sigf->Get(name);
						m2sig(i,j)->Scale(1.0/hjtpt->Integral());
						m2producer(i,j)->addSig(m2sig(i,j), m2mix(i,j));
						setupProducer(m2producer(i,j));
						m2producer(i,j)->doSeagullCorr = doSeagullCorr;
						if(doSeagullCorr) mm->CD(fmap(i,j).i1, fmap(i,j).i2);
						m2producer(i,j)->produce();
				}
		}
		if(doSeagullCorr) gAnaIO.saveCanvas((TCanvas*)mm, dataset_name+sig_name+"_s1_check_sg");
		if(doQA){
				signal_deta_qa();
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

void edmJtcAnalyzer::signal_deta_qa(){
		gQA->flash();
		TString label = doSeagullCorr ? "fixSg_" : "";
		matrixTH1Ptr m2sig_deta("deta_"+sig_name, m2sig.Nrow(), m2sig.Ncol() ); 		
		matrixTH1Ptr m2sig_side_deta("deta_side_"+sig_name, m2sig.Nrow(), m2sig.Ncol() ); 		
		gQA->pad_map = ps->getPara<mapper_func>("pad_map");
		gQA->pad_title = ps->getPara<TString (*)(int, int)>("pad_title");
		gQA->makeTitle = 1;
		for(int i=0; i<m2sig.Nrow(); ++i){
				for(int j=0; j<m2sig.Ncol(); ++j){
						auto h = jtc_utility::projX(1, m2producer(i,j)->sig, -1, .99, "");
						h->Scale(0.5);
						h->GetXaxis()->SetTitleSize(0.06);
						h->GetXaxis()->SetLabelSize(0.06);
						m2sig_deta.add(h, i, j);
						//h = jtc_utility::projX(1, m2producer(i,j)->sig_step2, 1.4, 1.8, "e");
						//h->Scale(1./0.4);
						//m2sig_side_deta.add(h, i, j);
						h = jtc_utility::projX(1, m2producer(i,j)->sig, 1.4, 1.8, "e");
						h->Scale(1./0.4);
						h->GetXaxis()->SetTitleSize(0.06);
						h->GetXaxis()->SetLabelSize(0.06);
						gQA->errorDrivenRange(h, -2.5, 2.5);
						m2sig_side_deta.add(h, i, j);
				}
		}
		gQA->autoYrange = 1;
		gQA->addhLine(0);
		gQA->addm2TH1(&m2sig_deta);
		gQA->x1=-3; gQA->x2 = 2.99;
		gQA->doSave=0;
		gAnaIO.saveCanvas(gQA->overlay(dataset_name+sig_name+"qa_Signal_deta_"+label), dataset_name+sig_name+"_s3_QASignal_deta_"+label);
		gQA->autoYrange = 0;
		gQA->flash();
		gQA->addm2TH1(&m2sig_side_deta);
		gQA->addm2TH1(&m2sig_deta);
		gAnaIO.saveCanvas(gQA->overlay(dataset_name+sig_name+"qa_seagull_"+label), dataset_name+sig_name+"_s2_QASeagull_deta_"+label);
		m2sig_deta.cleanAll();
		m2sig_side_deta.cleanAll();
		gQA->flash();
}
#endif
