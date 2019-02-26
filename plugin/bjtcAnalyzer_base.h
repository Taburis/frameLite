
#include "edmJtcAnalyzer.h"
#include "jtcTH1Player.h"

class bjtcAnalyzer_base : public rootEDMAnalyzer{
		public : bjtcAnalyzer_base(ParaSet & pset) : rootEDMAnalyzer(){
                         _nodeClassName_ = "bjtcAnalyzer_base";
                         ps = &pset;
                         loadStep2OutputFiles();
                         gAnaIO.output_plot_path = "./figs_check/";
                         ndrbin = ps->getPara<int>("ndr");
						 drbins = ps->getVectorAsArray<float>("drbins");
                 }
				 int analyze(){return 0;}
				 void loadStep2OutputFiles();
				 void loadStep3OutputFiles();
				 void slice_check();
		public : 
				 int npt= 6, ncent=1;
                 int ndrbin;
                 float *drbins;
                 ParaSet *ps;
                 TString inputPath;
				 jtc_utility::index2d (*pad_map) (int, int);
				 TFile *fbMC_gg, *fbMC_rg, *fbMC_rr, *fdMC_gg, *fdMC_rg, *fdMC_rr, *fdata, *fbdata_s3, *fidata_s3;
};

void bjtcAnalyzer_base::loadStep3OutputFiles(){
		TString folder = ps->getPara<TString>("step3output_folder");
		fbdata_s3 = TFile::Open(folder+"bjtc_data_final.root");
		fidata_s3 = TFile::Open(folder+"inclJet_data_final.root");
}

void bjtcAnalyzer_base::loadStep2OutputFiles(){
        TString path = ps->getPara<TString>("step2output_folder");
        TString cap = ps->getPara<TString>("bjetMC_step2output_name");
        TString file = jtc_utility::histNameScheme("_", 0, 1)+".root";
        fbMC_rg = TFile::Open(path+cap+file);
        file = jtc_utility::histNameScheme("_", 1, 1)+".root";
        fbMC_gg = TFile::Open(path+cap+file);
        file = jtc_utility::histNameScheme("_", 0, 0)+".root";
        fbMC_rr = TFile::Open(path+cap+file);
        //dijet MC loading section
        cap = ps->getPara<TString>("dijetMC_step2output_name");
        file = jtc_utility::histNameScheme("_", 0, 1)+".root";
        fdMC_rg = TFile::Open(path+cap+file);
        file = jtc_utility::histNameScheme("_", 1, 1)+".root";
        fdMC_gg = TFile::Open(path+cap+file);
        file = jtc_utility::histNameScheme("_", 0, 0)+".root";
        fdMC_rr = TFile::Open(path+cap+file);
        file = ps->getPara<TString>("pp5TeVData_step2output_name");
        fdata = TFile::Open(path+file+".root");
}



