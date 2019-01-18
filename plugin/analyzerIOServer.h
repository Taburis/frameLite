
#ifndef ANALYZERIOSERVER_H
#define ANALYZERIOSERVER_H

#include "jtcQaMonitor.h"  
#include "ParaSet.h"  
#include "TCanvas.h"  
#include "TFile.h"  
#include "TString.h"  

ParaSet* gCfg;

class AnalyzerIOServer {
		public : AnalyzerIOServer(){};
				 ~AnalyzerIOServer(){};
				 void setGCfg(ParaSet &ps){ cfg = &ps; gCfg = &ps;}
				 void saveCanvas(TCanvas *c, TString name){
						 c->SaveAs(output_plot_path+name+plotFormat);
						 delete c;
						 qam.flash();
				 }
				 TFile* bookRootFile(TString name, TString opt){
						auto f = TFile::Open(output_root_path+name, opt);
						return f;
				 }
				 TFile *getRootFile(TString name){
						 return TFile::Open(input_root_path+name);
				 }

				 ParaSet *cfg = nullptr;
				 TString output_plot_path = "./";
				 TString output_root_path = "./";
				 TString input_root_path  = "./";
				 TString plotFormat=".eps";
				 jtcQaMonitor qam;
};

AnalyzerIOServer gAnaIO;
jtcQaMonitor* gQAmonitor = &(gAnaIO.qam);

#endif
