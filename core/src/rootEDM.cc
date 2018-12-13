
#include "TFile.h"
#include "TTree.h"
#include "rootEDM.h"
rootEDMFrame * _edmFrame = 0;

TFile* rootEDMFrame::open(const char * f){
		_infile = TFile::Open(f);
		return _infile;
}

TTree* rootEDMFrame::regTree(const char* tn){
//		do not to use the add friend in case two branchs in different trees have the same name
		if(_forest.find(tn) == _forest.end()){
		  	_forest[tn]= (TTree*) _infile->Get(tn);
		}
		return _forest[tn];
}

TTree* rootEDMProducer::handle(const char* t){
		return _edmFrame->regTree(t);
}

void rootEDMFrame::addProducer(rootEDMProducer * ana){
		_analyzer.push_back(ana);
}
void rootEDMFrame::addEventFilter(rootEDMEventFilter * filter){
		_evtfilter.push_back(filter);
}

void rootEDMFrame::runBeginSection(){
		if( _analyzer.size() != 0 )
			   for(auto _ip =  _analyzer.begin(); _ip!=_analyzer.end(); ++_ip){
					   (*_ip)->beginJob();
			   }
		if( _evtfilter.size()!= 0 )
			   for(auto _ip =  _evtfilter.begin(); _ip!=_evtfilter.end(); ++_ip) (*_ip)->beginJob();
		std::cout<<"the beginJobs complete."<<std::endl;
}

void rootEDMFrame::runEndSection(){
		if( _analyzer.size() != 0 )
			   for(auto _ip =  _analyzer.begin(); _ip!=_analyzer.end(); ++_ip) (*_ip)->endJob();
		if( _evtfilter.size()!= 0 )
			   for(auto _ip =  _evtfilter.begin(); _ip!=_evtfilter.end(); ++_ip) (*_ip)->endJob();
		std::cout<<"endJobs complete. "<<std::endl;
}

void rootEDMFrame::runEventSection(){
		// filter has been treated as AND path, only the event pass all filter will be selected 
		if( _evtfilter.size()!= 0) 
			   for(auto _ip =  _evtfilter.begin(); _ip!=_evtfilter.end(); ++_ip) 
					   if(!(*_ip)->filter()) return;
		if( _analyzer.size() != 0 )
				for(auto _ip =  _analyzer.begin(); _ip!=_analyzer.end(); ++_ip)
						(*_ip)->analyzer();
}

void rootEDMFrame::eventLoop(){
		_edmFrame = this;
		_infile->cd();
		runBeginSection();
		if( _forest.size() == 0){
				std::cout<<"no tree has been registered, skipped the event loop!"<<std::endl;;
			   	return;
		}
		Long64_t nevt = (*_forest.begin()).second->GetEntries();
		if( EventMax == -1 ) EventMax = nevt;
	//	TTree* t = (TTree*)_infile->Get("genParticleTree");
		std::cout<<"starting the event loop.."<<std::endl;
		for(Long64_t ievt=0; ievt< EventMax ; ++ievt){
				if(ievt % 1000 == 0) std::cout<<"processing "<<ievt<<"th event.."<<std::endl;
				for(auto _ip= _forest.begin(); _ip!= _forest.end(); ++_ip) {
						((*_ip).second)->GetEntry(ievt); 
				}
				runEventSection();
		}

		runEndSection();
		_done();
//		_infile->Close();
}



