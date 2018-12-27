
//#include <mach/mach.h>
#ifndef rootEDM_H
#define rootEDM_H
#include <vector>
#include <unordered_map>
#include <iostream>
#include "TTree.h"
#include "node.h"
#include "ParaSet.h"


class rootEDMPlugin { 
		// root event data model base class, used for any class related to root analysis
		public : rootEDMPlugin(){}; ~rootEDMPlugin(){};
				 //					 virtual void begin(const PSet &) = 0;
};


class rootEDMProducer : public rootEDMPlugin {
		public : rootEDMProducer (){};
				 virtual ~rootEDMProducer () {};
				 TTree* handle(const char * t);
				 virtual void beginJob() = 0; 
				 virtual void analyze() = 0; 
				 virtual void endJob() = 0;
};

class rootEDMAnalyzer : public node{
		public : rootEDMAnalyzer(){_nodeClassName_ = "rootEDMAnalyzer";};
				 virtual ~rootEDMAnalyzer(){};
				 virtual int analyze()=0;
				 virtual int evaluate() { return analyze();}
				 void handle_cfg(ParaSet & ps ) { cfg = &ps;}
				 ParaSet * cfg;
};

class rootEDMEventFilter : public rootEDMPlugin {
		public : rootEDMEventFilter(){};
				 virtual ~rootEDMEventFilter(){};
				 virtual void beginJob() = 0;
				 virtual void endJob() = 0;
				 // convention for the filter:
				 // return 1 means keep the event
				 virtual bool filter() = 0;
};

class rootEDMFrame : public node{
		// every edm has to be run with a frame 
		public : rootEDMFrame(){ _nodeClassName_ = "rootEDMFrame";};
				 ~rootEDMFrame(){
						 std::cout<<"deleting rootEDMFrame"<<std::endl;
						 for(auto & it : _analyzer) delete it;
						 for(auto & it : _evtfilter) delete it;
						 for(auto & it : _forest) delete it.second;
				 };
				 void eventLoop();
				 TFile* open(const char*);
				 TTree* regTree(const char* t);
				 void addProducer(rootEDMProducer * );
				 void addEventFilter(rootEDMEventFilter * );
				 void runBeginSection();
				 void runEndSection();
				 void runEventSection();
				 virtual int evaluate(){
						 eventLoop();
						 _done();
						 return 0;
				 };
		public :
				 TFile* _infile;
				 Long64_t EventMax = 1;
				 std::vector<rootEDMProducer * > _analyzer;
				 std::vector<rootEDMEventFilter * > _evtfilter;
				 std::vector<void (*)()> begin_seq;
				 std::unordered_map<std::string, TTree*> _forest;
};


#endif
