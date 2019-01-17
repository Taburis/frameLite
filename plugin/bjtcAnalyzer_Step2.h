
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
				edmJtcAnalyzer * initAnalyzer();
				// 101 workflow to pull the signal from a standard bjtc correlation files which include: 1. inclusive jet 2. tagged b jet 3. tagged & true bjet 4. true b jet
				void wf001(const char * input, const char *output_caption, bool isGenJet, bool isGenTrack);
				void datawf001(const char *input, const char *output_caption);
				void bjetMCwf101();
				void dijetMCwf101();
				void Datawf101();
				void testWf001();
				virtual int analyze();

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


edmJtcAnalyzer * bjtcAnalyzer_Step2::initAnalyzer(){
		auto an = new edmJtcAnalyzer(*ps);
		an->open(infile);
		an->doQA = doQA;
		an->saveFile=doSave;
		an->output =output_folder+outputname;
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
		initAnalyzer();
		sig_name = "inclJet_Data_pTweighted";
		openOpt = "update";
		initAnalyzer();

		outputname = TString(output_caption) + ".root";
		sig_name = "taggedBJet_Data";
		mix_name = "taggedBJet_Data_mixing";
		norm_name = "inclJet_corrpt";
		openOpt = "update";
		initAnalyzer();
		sig_name = "taggedBJet_Data_pTweighted";
		openOpt = "update";
		initAnalyzer();
}

void bjtcAnalyzer_Step2::wf001(const char * input, const char *output_caption, bool isGenJet, bool isGenTrack){
		//for inclusive jet
		infile = input;
		outputname = output_caption + jtc_utility::histNameScheme("_", isGenJet, isGenTrack)+".root";
		sig_name = jtc_utility::histNameScheme("inclJet_", isGenJet, isGenTrack);
		mix_name = jtc_utility::histNameScheme("inclJet_", isGenJet, isGenTrack, 0, 1);
		norm_name = "inclJet_corrpt";
		openOpt = "recreate";
		initAnalyzer();
		sig_name = jtc_utility::histNameScheme("inclJet_", isGenJet, isGenTrack, 1);
		openOpt = "update";
		initAnalyzer();

		//for tagged bjet
		outputname = output_caption+jtc_utility::histNameScheme("_", isGenJet, isGenTrack)+".root";
		sig_name = jtc_utility::histNameScheme("taggedBJet_", isGenJet, isGenTrack);
		mix_name = jtc_utility::histNameScheme("taggedBJet_", isGenJet, isGenTrack, 0, 1);
		norm_name = "taggedBJet_corrpt";
		openOpt = "update";
		initAnalyzer();
		sig_name = jtc_utility::histNameScheme("taggedBJet_", isGenJet, isGenTrack, 1);
		openOpt = "update";
		initAnalyzer();

		//for tagged bjet
		outputname = output_caption+jtc_utility::histNameScheme("_", isGenJet, isGenTrack)+".root";
		sig_name = jtc_utility::histNameScheme("trueBJet_", isGenJet, isGenTrack);
		mix_name = jtc_utility::histNameScheme("trueBJet_", isGenJet, isGenTrack, 0, 1);
		norm_name = "trueBJet_corrpt";
		openOpt = "update";
		initAnalyzer();
		sig_name = jtc_utility::histNameScheme("trueBJet_", isGenJet, isGenTrack, 1);
		openOpt = "update";
		initAnalyzer();

		//for tagged and true bjet
		outputname = output_caption+jtc_utility::histNameScheme("_", isGenJet, isGenTrack)+".root";
		sig_name = jtc_utility::histNameScheme("taggedTrueBJet_", isGenJet, isGenTrack);
		mix_name = jtc_utility::histNameScheme("taggedTrueBJet_", isGenJet, isGenTrack, 0, 1);
		norm_name = "taggedTrueBJet_corrpt";
		openOpt = "update";
		initAnalyzer();
		sig_name = jtc_utility::histNameScheme("taggedTrueBJet_", isGenJet, isGenTrack, 1);
		openOpt = "update";
		initAnalyzer();
}

void bjtcAnalyzer_Step2::testWf001(){
		output_folder = "./";
		TString output_name = ps->getPara<TString>("bjetMC_step2output_name");
		//Reco-Gen
		TString input_file = ps->getPara<TString>("bjetMC_step2input_rg_file");
		wf001(input_file, output_name, 0, 1);
}

void bjtcAnalyzer_Step2::Datawf101(){
		auto input_file = ps->getPara<TString>("pp5TeVData_step2input_file");
		auto output_name = ps->getPara<TString>("pp5TeVData_step2output_name");
		datawf001(input_file, output_name);
}

void bjtcAnalyzer_Step2::bjetMCwf101(){
		// when all the bjet samples and data are ready and start to pull the signal from all of them:
		TString output_name = ps->getPara<TString>("bjetMC_step2output_name");
		TString input_file;
		input_file = ps->getPara<TString>("bjetMC_step2input_rg_file");
		wf001(input_file, output_name, 0, 1);
		input_file = ps->getPara<TString>("bjetMC_step2input_gg_file");
		wf001(input_file, output_name, 1, 1);
		input_file = ps->getPara<TString>("bjetMC_step2input_rr_file");
		wf001(input_file, output_name, 0, 0);
}

void bjtcAnalyzer_Step2::dijetMCwf101(){
		TString output_name = ps->getPara<TString>("dijetMC_step2output_name");
		TString input_file;
		input_file = ps->getPara<TString>("dijetMC_step2input_rg_file");
		wf001(input_file, output_name, 0, 1);
		input_file = ps->getPara<TString>("dijetMC_step2input_gg_file");
		wf001(input_file, output_name, 1, 1);
		input_file = ps->getPara<TString>("dijetMC_step2input_rr_file");
		wf001(input_file, output_name, 0, 0);
}

int bjtcAnalyzer_Step2::analyze(){
		printClassName();
		for(auto &it : analyzers){
			   	it->evaluate();
				delete it;
		}
		return 0;
}
