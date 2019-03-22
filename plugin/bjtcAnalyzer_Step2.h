
#include "edmJtcAnalyzer.h"

//name convention for bjtc_pp_an


class bjtcAnalyzer_Step2 : public rootEDMAnalyzer {
		public : 
				bjtcAnalyzer_Step2(): rootEDMAnalyzer(){
						_nodeClassName_ = "bjtcAnalyzer_Step2";
				}
				bjtcAnalyzer_Step2(ParaSet & p);
				~bjtcAnalyzer_Step2(){// for(auto &it : analyzers) delete it;
				}
				edmJtcAnalyzer * initAnalyzer(const char * name);
				// 101 workflow to pull the signal from a standard bjtc correlation files which include: 1. inclusive jet 2. tagged b jet 3. tagged & true bjet 4. true b jet
				void wf001(const char * dataset_name, const char * input, const char *output_caption, bool isGenJet, bool isGenTrack, bool doCont = 0);
				void datawf001(const char *input, const char *output_caption);
				void bjetMCwf101();
				void dijetMCwf101();
				void PYTHIAwf101();
				void Datawf101();
				void Syswf101();
				void testWf001();
				void inclCalo_Wf001();
				virtual int analyze();
				void jerMCwf101();

		public : std::vector<edmJtcAnalyzer*> analyzers;
				 ParaSet *ps;
				 // config variables
				 TString sig_name, mix_name, norm_name, outputname,output_folder, openOpt, infile;
				 bool doSave = 1;
				 bool doQA = 1;
				 bool doSeagullCorr = 0;
};

bjtcAnalyzer_Step2::bjtcAnalyzer_Step2(ParaSet &ps_): rootEDMAnalyzer(){
		ps = &ps_;
		_nodeClassName_ = "bjtcAnalyzer_Step2";
		doQA = ps_.getPara<bool>("doQA");
		output_folder = ps->getPara<TString>("step2output_folder");
		gAnaIO.output_plot_path =output_folder; 
		doSeagullCorr = ps_.getPara<bool>("doSeagullCorr");
}


edmJtcAnalyzer * bjtcAnalyzer_Step2::initAnalyzer(const char* name){
		auto an = new edmJtcAnalyzer(*ps);
		an->open(infile);
		an->doQA = doQA;
		an->saveFile=doSave;
		an->output =output_folder+outputname;
		an->dataset_name =name;
		an->sig_name = sig_name+"_noCorr";
		an->mix_name = mix_name+"_noCorr";
		an->norm_name = norm_name;
		an->outputOpt = openOpt;
		an->doSeagullCorr = doSeagullCorr;
		analyzers.push_back(an);
		return an;	
}

void bjtcAnalyzer_Step2::datawf001(const char * input, const char *output_caption){
		infile = input;
		outputname = TString(output_caption) + ".root";
		sig_name = "inclJet_Data";
		mix_name = "inclJet_Data_mixing";
		norm_name = "inclJet_corrpt";
		openOpt = "recreate";
		initAnalyzer("pp5TeVjet80_");
		sig_name = "inclJet_Data_pTweighted";
		openOpt = "update";
		initAnalyzer("pp5TeVjet80_");

		outputname = TString(output_caption) + ".root";
		sig_name = "taggedBJet_Data";
		mix_name = "taggedBJet_Data_mixing";
		norm_name = "taggedBJet_corrpt";
		openOpt = "update";
		initAnalyzer("pp5TeVjet80_");
		sig_name = "taggedBJet_Data_pTweighted";
		openOpt = "update";
		initAnalyzer("pp5TeVjet80_");
}

void bjtcAnalyzer_Step2::inclCalo_Wf001(){
		auto output_caption = ps->getPara<TString>("inclCaloJet_MC_step2output_name");

		if( ps->exists("inclCaloJet_step2input_gg_file")){
				infile = ps->getPara<TString> ("inclCaloJet_step2input_gg_file");
				outputname = output_caption + jtc_utility::histNameScheme("_", 1, 1)+".root";
				sig_name = jtc_utility::histNameScheme("inclJet_", 1, 1);
				mix_name = jtc_utility::histNameScheme("inclJet_", 1, 1, 0, 1);
				norm_name = "inclJet_corrpt";
				openOpt = "recreate";
				initAnalyzer("caloJet_");
				sig_name = jtc_utility::histNameScheme("inclJet_", 1, 1, 1);
				openOpt = "update";
				initAnalyzer("caloJet_");	}

		if( ps->exists("inclCaloJet_step2input_rg_file")){
				infile = ps->getPara<TString> ("inclCaloJet_step2input_rg_file");
				outputname = output_caption + jtc_utility::histNameScheme("_", 0, 1)+".root";
				sig_name = jtc_utility::histNameScheme("inclJet_", 0, 1);
				mix_name = jtc_utility::histNameScheme("inclJet_", 0, 1, 0, 1);
				norm_name = "inclCaloJet_corrpt";
				openOpt = "recreate";
				initAnalyzer("caloJet_");
				sig_name = jtc_utility::histNameScheme("inclJet_", 0, 1, 1);
				openOpt = "update";
				initAnalyzer("caloJet_");	}

		if( ps->exists("inclCaloJet_step2input_rr_file")){
				infile = ps->getPara<TString> ("inclCaloJet_step2input_rr_file");
				outputname = output_caption + jtc_utility::histNameScheme("_", 0, 1)+".root";
				sig_name = jtc_utility::histNameScheme("inclJet_", 0, 0);
				mix_name = jtc_utility::histNameScheme("inclJet_", 0, 0, 0, 1);
				norm_name = "inclCaloJet_corrpt";
				openOpt = "recreate";
				initAnalyzer("caloJet_");
				sig_name = jtc_utility::histNameScheme("inclJet_", 0, 0, 1);
				openOpt = "update";
				initAnalyzer("caloJet_");	}
}

void bjtcAnalyzer_Step2::wf001(const char* dataset_name, const char * input, const char *output_caption, bool isGenJet, bool isGenTrack, bool doCont){
		//for inclusive jet
		infile = input;
		outputname = output_caption + jtc_utility::histNameScheme("_", isGenJet, isGenTrack)+".root";
		sig_name = jtc_utility::histNameScheme("inclJet_", isGenJet, isGenTrack);
		mix_name = jtc_utility::histNameScheme("inclJet_", isGenJet, isGenTrack, 0, 1);
		norm_name = "inclJet_corrpt";
		openOpt = "recreate";
		initAnalyzer(dataset_name);
		sig_name = jtc_utility::histNameScheme("inclJet_", isGenJet, isGenTrack, 1);
		openOpt = "update";
		initAnalyzer(dataset_name);

		//for tagged bjet
		outputname = output_caption+jtc_utility::histNameScheme("_", isGenJet, isGenTrack)+".root";
		sig_name = jtc_utility::histNameScheme("taggedBJet_", isGenJet, isGenTrack);
		mix_name = jtc_utility::histNameScheme("taggedBJet_", isGenJet, isGenTrack, 0, 1);
		norm_name = "taggedBJet_corrpt";
		openOpt = "update";
		initAnalyzer(dataset_name);
		sig_name = jtc_utility::histNameScheme("taggedBJet_", isGenJet, isGenTrack, 1);
		openOpt = "update";
		initAnalyzer(dataset_name);

		//for tagged bjet
		outputname = output_caption+jtc_utility::histNameScheme("_", isGenJet, isGenTrack)+".root";
		sig_name = jtc_utility::histNameScheme("trueBJet_", isGenJet, isGenTrack);
		mix_name = jtc_utility::histNameScheme("trueBJet_", isGenJet, isGenTrack, 0, 1);
		norm_name = "trueBJet_corrpt";
		openOpt = "update";
		initAnalyzer(dataset_name);
		sig_name = jtc_utility::histNameScheme("trueBJet_", isGenJet, isGenTrack, 1);
		openOpt = "update";
		initAnalyzer(dataset_name);

		//for tagged and true bjet
		outputname = output_caption+jtc_utility::histNameScheme("_", isGenJet, isGenTrack)+".root";
		sig_name = jtc_utility::histNameScheme("taggedTrueBJet_", isGenJet, isGenTrack);
		mix_name = jtc_utility::histNameScheme("taggedTrueBJet_", isGenJet, isGenTrack, 0, 1);
		norm_name = "taggedTrueBJet_corrpt";
		openOpt = "update";
		initAnalyzer(dataset_name);
		sig_name = jtc_utility::histNameScheme("taggedTrueBJet_", isGenJet, isGenTrack, 1);
		openOpt = "update";
		initAnalyzer(dataset_name);

		if(doCont){
				outputname = output_caption+jtc_utility::histNameScheme("_", isGenJet, isGenTrack)+".root";
				sig_name = jtc_utility::histNameScheme("contJet_", isGenJet, isGenTrack);
				mix_name = jtc_utility::histNameScheme("contJet_", isGenJet, isGenTrack, 0, 1);
				norm_name = "contJet_corrpt";
				openOpt = "update";
				initAnalyzer(dataset_name);
				sig_name = jtc_utility::histNameScheme("contJet_", isGenJet, isGenTrack, 1);
				openOpt = "update";
				initAnalyzer(dataset_name);
		}
}

void bjtcAnalyzer_Step2::testWf001(){
		output_folder = "./";
		TString output_name = ps->getPara<TString>("bjetMC_step2output_name");
		//Reco-Gen
		TString input_file = ps->getPara<TString>("bjetMC_step2input_rg_file");
		wf001("testSet", input_file, output_name, 0, 1);
}

void bjtcAnalyzer_Step2::Datawf101(){
		auto input_file = ps->getPara<TString>("pp5TeVData_step2input_file");
		auto output_name = ps->getPara<TString>("pp5TeVData_step2output_name");
		datawf001(input_file, output_name);
}
void bjtcAnalyzer_Step2::PYTHIAwf101(){
		TString output_name = "Signal_PYTHIA6_bjetSample_noGPS_GenJet_GenTrack.root";
		TString input_file;
		input_file = ps->getPara<TString>("bjetMC_step2input_gg_origin_file");
		wf001("bjetMCOrigin_",input_file, output_name, 1, 1);
}
void bjtcAnalyzer_Step2::jerMCwf101(){
		auto input_file = ps->getPara<TString>("dijetMC_step2input_rg_gA_file");
		auto output_name = ps->getPara<TString>("dijetMC_step2output_name");
		output_name = output_name+"_genAxis";
		wf001("dijetMC_GA_",input_file, output_name, 0, 1);

		input_file = ps->getPara<TString>("bjetMC_step2input_rg_gA_file");
		output_name = ps->getPara<TString>("bjetMC_step2output_name");
		output_name = output_name+"_genAxis";
	//	wf001("bjetMC_GA_",input_file, output_name, 0, 1);
}
void bjtcAnalyzer_Step2::bjetMCwf101(){
		// when all the bjet samples and data are ready and start to pull the signal from all of them:
		TString output_name = ps->getPara<TString>("bjetMC_step2output_name");
		TString input_file;
		input_file = ps->getPara<TString>("bjetMC_step2input_rg_file");
		wf001("bjetMC_",input_file, output_name, 0, 1);
		input_file = ps->getPara<TString>("bjetMC_step2input_gg_file");
		wf001("bjetMC_",input_file, output_name, 1, 1);
		input_file = ps->getPara<TString>("bjetMC_step2input_rr_file");
		wf001("bjetMC_",input_file, output_name, 0, 0);
}

void bjtcAnalyzer_Step2::Syswf101(){
		TString output_name ="Signal_PYTHIA6_bjetSample_JER";
		TString input_file;
		input_file = "/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/pp_28Jan19/bJTC_PYTHIA6_RecoGen_5TeV_bjetMC_GSPweight_WTAaxis_JER15_csvV2p9_14Mar19.root";
		wf001("bjetJER_",input_file, output_name, 0, 1);
}

void bjtcAnalyzer_Step2::dijetMCwf101(){
		TString output_name = ps->getPara<TString>("dijetMC_step2output_name");
		TString input_file;
		input_file = ps->getPara<TString>("dijetMC_step2input_rg_file");
		wf001("dijetMC_", input_file, output_name, 0, 1, 1);
		input_file = ps->getPara<TString>("dijetMC_step2input_gg_file");
		wf001("dijetMC_", input_file, output_name, 1, 1, 1);
		input_file = ps->getPara<TString>("dijetMC_step2input_rr_file");
		wf001("dijetMC_", input_file, output_name, 0, 0, 1);
}

int bjtcAnalyzer_Step2::analyze(){
		printClassName();
		for(auto &it : analyzers){
				it->evaluate();
				delete it;
		}
		return 0;
}
