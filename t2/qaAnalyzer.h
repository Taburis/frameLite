
#include "ParaSet.h"

class qaAnalyzer{
	public : 
		qaAnalyzer(){}
		void loadConfig(ParaSet &ps0){ps = &ps0;}
		void getVzWeight();
		void getPThatWeight();
		void loadFiles(TString f0str, TString f1str){
	//usually f0 is mc and f1 is data
			f0 = TFile::Open(f0str);
			f1 = TFile::Open(f1str);
		}
		void check(TString , TString);
		void style(TH1* h, TString ,TString, Color_t);
		void normalize(TH1* );
		void makePythiaConfig(TString);
	public :
		TFile * f0, *f1;
		ParaSet *ps;
		TF1 *vzmcf, *vzdataf;
		TH1F* hpthatRatio, *hvzRatio;
		std::vector<float> pThatInterval, pThatEntries;
		std::vector<Double_t>xsec;
		TTree* cfgTree;
		int DBCode;
};

void qaAnalyzer::style(TH1* h, TString xtitle, TString ytitle, Color_t col){
	h->GetXaxis()->CenterTitle();
	h->GetYaxis()->CenterTitle();
	h->GetXaxis()->SetTitle(xtitle);
	h->GetYaxis()->SetTitle(ytitle);
	h->SetMarkerStyle(20);
	h->SetMarkerSize(0.5);
	h->SetMarkerColor(col);
}

void qaAnalyzer::normalize(TH1* h){
	h->Scale(1.0/h->Integral());
}

void qaAnalyzer::getVzWeight(){
	TH1F* hvzmc = (TH1F*)f0->Get("vz");
	TH1F* hvzdata = (TH1F*)f1->Get("vz");
//	vzmcf = new TF1("vzMC_fit", "gaus(0)", -15, 15);
//	vzdataf = new TF1("vzData_fit", "gaus(0)", -15, 15);
//	normalize(hvzmc);
//	normalize(hvzdata);
//	int bin0= hvzmc->FindBin(0.0);
//	float cont = hvzmc->GetBinContent(bin0);
//	vzmcf->SetParameters(cont, 0, 1);
//	vzdataf->SetParameters(cont, 0, 1);
//	auto c = new TCanvas("c", "vz distribution", 1000, 450);
//	c->Divide(2,1);
//c->cd(1);
//	style(hvzmc, "Vz", "dN/dvz/N", kBlue);
//	hvzmc->Fit(vzmcf);
//c->cd(2);
//	style(hvzdata, "Vz", "dN/dvz/N", kBlue);
//	hvzdata->Fit(vzdataf);
	hvzRatio = (TH1F*) hvzdata ->Clone("vzWeight");
	hvzRatio->Divide(hvzmc);
}

void qaAnalyzer::getPThatWeight(){
	auto hpthatmc = (TH1F*) f0->Get("pthatStat");
	cfgTree = new TTree("cfgTree", "");
	cfgTree->Branch("pThatInterval", &pThatInterval);
	cfgTree->Branch("pThatEntries", &pThatEntries);
	cfgTree->Branch("xsec", &xsec);
	cfgTree->Branch("DBCode", &DBCode);
	int npthatbin = ps->getPara<int>("npThat");
	DBCode = ps->getPara<int>("DBCode");
	const Double_t * xsec0 = ps->getPara<const Double_t*>("xsec");
	const Float_t * pthatEdge = ps->getPara<const Float_t*>("pThatInterval");
	pThatEntries.clear();
	for(int i = 0; i<npthatbin; ++i){
		xsec.push_back(xsec0[i]);				
		pThatInterval.push_back(float(pthatEdge[i]));
		float content = hpthatmc->GetBinContent(i+1);
		if(DBCode == 1 && i == npthatbin-1) content = content / 1.15;
		pThatEntries.push_back(content);
	}
	xsec.push_back(0.0);
	pThatInterval.push_back(99999);
	cfgTree->Fill();
}
void qaAnalyzer::makePythiaConfig(TString ofstr){
	auto wf = TFile::Open(ofstr, "recreate");
	wf->cd();	
	getVzWeight();
	getPThatWeight();

	//hvzRatio->Write();
//	vzdataf->Write();
//	vzmcf->Write();
	wf->Write();
	wf->Close();
}

void qaAnalyzer::check(TString var, TString ytitle){
	auto h0 = (TH1F*) f0->Get(var);
	auto h1 = (TH1F*) f1->Get(var);
	style(h0, var, ytitle, kBlue);
	style(h1, var, ytitle, kRed);
	normalize(h0);
	normalize(h1);
	auto c = new TCanvas("c_"+var, "", 500, 450);
	h0->Draw();
	h1->Draw("same");
}

