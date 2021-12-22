// ================================================================================
//
//		2021.08.24 Li YI
//		caculate and plot Rcp for JETSCAPE & JETSCAPE-noFT
//
// ================================================================================
// 	
// 		2021.12.20 Li YI
// 		Rcp	
//
// ================================================================================
#include <sstream>	// stringstream
#include <fstream>	// ifstream
#include <iostream>	
#include <string>	// string
#include <dirent.h>

#include "TH1.h"
#include "TCanvas.h"
#include "TText.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TFile.h"
#include "TLine.h"
#include "TFrame.h"


#define DRAW		// if define draw plots
#define SAVEROOT	// if define save root file

using std::cout; using std::cin;
using std::endl; using std::vector;

int verbose = 0;

// ------------------ Functions ----------------

// Read Each input files (jet line only) & fill histogram
bool ReadJet(char *filename, TH1D *h, double &Nevent) {
	ifstream in(filename);
	if(in==NULL && !in.is_open()) {cout<<"ERROR:: "<< __PRETTY_FUNCTION__ << " NULL input file " << filename <<endl; return false;}
	if(h==NULL) {cout<<"ERROR:: "<< __PRETTY_FUNCTION__ << " NULL histogram. "<<endl; return false;}
	if(verbose) cout<<"Read "<<filename<<endl;
	string line;
	Nevent = 0;
	while (!in.eof() && getline(in, line)){
		stringstream ss(line);
		string tmp;
		vector<double> v;

		std::istringstream iss(line); 
		for(std::string s; iss >> s; ) 	{
			v.push_back(stod(s)); 
			if(verbose == 10) cout<<stod(s)<<endl;
		}
		if(v.size()==3) {	// read in an event
			Nevent++;
			if(verbose==10) cout<<Nevent<<endl;
		}
		if(v.size()==5) {	// read in a jet information line
			if(verbose==10) cout<<v.at(3)<<'\t'<<v.at(4)<<endl;
			if(v.at(3)>=0) {
				h->Fill(v.at(4));
			}
			else {
				h->Fill(v.at(4),-1);	// if negative (for those from initial medium particles)
			}
		}
	}
	in.close();
	return true;
}

// List all files in directory
bool ListDir(const char *path, const char *tag, vector<string>& fileList, const char *tag2="") {
	struct dirent *entry;
	DIR *dir = opendir(path);

	fileList.clear();	// clear set vector size 0

	if (dir == NULL) {
		cout<<"ERROR:: "<< __PRETTY_FUNCTION__ << " NULL directory " << path <<endl; 
		return false;
	}
	while ((entry = readdir(dir)) != NULL) {
		if(verbose==10) cout << entry->d_name << endl;
		string s = entry->d_name;
                if(s.find(tag)!=string::npos && s.find(tag2)!=string::npos) {	// if found filename with tag
			if(verbose) cout<<"Push back "<<s<<endl;
			fileList.push_back(s);

		}		
	}
	closedir(dir);

	return true;
}

bool ReadXsec(const char *dir, const int Npt, const int *ptmin, const int *ptmax, double *xsec) {
	char name[500];
	ifstream input;

	for(int i = 0; i<Npt; i++) {
		sprintf(name,"%sXsec_pTHat%03d_to%03d.dat",dir, ptmin[i], ptmax[i]);
		input.open(name,std::ifstream::in);
		if(input && input.is_open()) {
			double nevt = 0, ievt = 0;
			double sumxsec = 0, ixsec = 0;
			double tmp;
			input >> ievt >> ixsec >> tmp >> tmp;
			while(input.good()) {
				//cout<<ievt<<'\t'<<ixsec<<endl;
				nevt+=ievt;
				sumxsec+=ievt*ixsec;
				input >> ievt >> ixsec >> tmp >> tmp;
			}
			input.close();
			//cout<<sumxsec<<'\t'<<nevt<<endl;
			if(nevt>0) sumxsec/=nevt;
			xsec[i] = sumxsec;
			//cout<<xsec[i]<<endl;
		}
		else {
			if(verbose) cout<<"ERROR: "<<__PRETTY_FUNCTION__<<" cannot open "<<name<<endl;
			return false;
		}
	}
	return true;
}

// ------------------ Main ---------------------
void compareRcp() {

	enum Type {kPeri, kLBT, kPeriNoFT, kLBTNoFT, kPeri4S, kLBTSimple, kTotal};
	const char *name[kTotal] = {"Peri","LBT","PeriNoFT","LBTNoFT","Peri4S","LBTSimple"};
	const char *tag[kTotal] = {"LBT-jet", "LBT-jet", "LBT-jet", "LBT-jet", "LBT-jet", "LBT-jet"};
	const char *LegendName[kTotal] = {"peripheral","LBT w/ f.t.","peri default","LBT default","peripheral4S","LBT w/ simple f.t."};

	// input 
	const char * dir[kTotal] = {"../JETSCAPE-6080/results/","../JETSCAPE/results/", "../JETSCAPE-noFT-6080/results/", "../JETSCAPE-noFT/results/", "../JETSCAPE-6080/results-SimpleIF/", "../JETSCAPE/results-SimpleIF/"};

	ifstream input[kTotal]; 

	const int Npt = 14+1;		// +1 for the sum
	int arrayMin[Npt]= {0 , 2 , 4 , 6  , 8 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60, 0};		// last bin is sum
	int arrayMax[Npt]={2 , 4 , 6 , 8 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60  ,100, 100};		// last bin is sum
	TH1D *hist[kTotal][Npt];	

	double Xsec[Npt-1] = {28.4852,9.61467,0.472347,0.0635027,0.0134543,0.00527283,0.000456384,6.42e-05,1.16125e-05,2.39501e-06,5.32087e-07,1.49413e-07,7.60858e-09,3.33557e-10};
	double Xsec4S[Npt-1] = {28.4852,9.61467,0.472347,0.0635027,0.0134543,0.00527283,0.000456384,6.42e-05,1.16125e-05,2.39501e-06,5.32087e-07,1.49413e-07,7.60858e-09,3.33557e-10};

	double tmpsec[Npt-1] = {0};
	if(ReadXsec("../JETSCAPE/pythiaInput/",Npt-1,arrayMin,arrayMax,tmpsec)) {	// read xsec file if exists
		for(int i = 0; i<Npt-1; i++) Xsec[i] = tmpsec[i];
	}
	if(ReadXsec("../JETSCAPE/pythiaInput-SimpleIF/",Npt-1,arrayMin,arrayMax,tmpsec)) {	// read xsec file if exists
		for(int i = 0; i<Npt-1; i++) Xsec4S[i] = tmpsec[i];
	}

	//double SumXsec=0;
	//for(int i = 0; i<Npt-1; i++) SumXsec+=Xsec[i];

	if(1) {
		cout<<"ptHat\tCross section"<<endl;
		for(int i = 0; i<Npt-1; i++){
			cout<<arrayMin[i]<<'-'<<arrayMax[i]<<'\t'<<Xsec[i]<<'\t'<<Xsec4S[i]<<endl;
		}
		cout<<endl;
	}



	int Nbins = 100;
	double xmin = 0, xmax = 100;
	for(int t = 0; t<kTotal; t++) {
		for(int i = 0; i<Npt; i++) {
			hist[t][i] = new TH1D(Form("h%s_pt%d",name[t],i),Form("%s jet for %d<ptHat<%d",name[t],arrayMin[i],arrayMax[i]),Nbins,xmin,xmax);
			hist[t][i]->Sumw2();
			hist[t][i]->GetXaxis()->SetTitle("jet_pT_sub");
		}
	}


	// Loop over pt bins
	for(int i = 0; i<Npt-1; i++) {	// the last bin for sum over all pt bins
		vector<string> list[kTotal];
		for(int t = 0; t<kTotal; t++) {
			char ctmp[500] = "";
			sprintf(ctmp,"pTHat%03d_to%03d",arrayMin[i],arrayMax[i]);
			if(verbose) cout<<"Finding "<<ctmp<<endl;
			ListDir(dir[t], tag[t], list[t], ctmp);

			double Nevt = 0;	// total number of events for type t & pt i
			for(auto it : list[t]) {	// for each files in list
				char ifilename[500];
				double iNevt = 0; 
				sprintf(ifilename,"%s%s",dir[t],it.data()); 
				ReadJet(ifilename, hist[t][i], iNevt);
				Nevt=+iNevt;
			}	// End of Loop over list of files for type t & pt i
			if(Nevt>0) hist[t][i]->Scale(1./Nevt);	// Normalized by number of events for type k & pt i
			if(t==kPeri4S || t==kLBTSimple) hist[t][i]->Scale(Xsec4S[i]);
			else hist[t][i]->Scale(Xsec[i]);	// Normalized by cross section
		}	// End of Loop over type t
	}	// End of Loop over pt i

	for(int t = 0; t<kTotal; t++) {
		for(int i = 0; i<Npt-1; i++) {
			hist[t][Npt-1]->Add(hist[t][i]);	// the sum of all pt bins
			//cout<<"Add "<<hist[t][i]->GetName()<<" to "<<hist[t][Npt-1]<<endl;
		}
	}


	// Rcp
	TH1D *hRcp = (TH1D*)hist[kLBT][Npt-1]->Clone("hRcp");hRcp->SetTitle("R_{CP} w/ timing");
	TH1D *hRcpNoFT = (TH1D*)hist[kLBTNoFT][Npt-1]->Clone("hRcpNoFT");hRcpNoFT->SetTitle("R_{CP}");
	TH1D *hRcpSimpleFT = (TH1D*)hist[kLBTSimple][Npt-1]->Clone("hRcpSimpleFT");hRcpSimpleFT->SetTitle("R_{CP} w/ simple timing");
	hRcp->Divide(hist[kPeri][Npt-1]);
	hRcpNoFT->Divide(hist[kPeriNoFT][Npt-1]);
	hRcpSimpleFT->Divide(hist[kPeri4S][Npt-1]);
	hRcp->GetYaxis()->SetTitle("R_{CP}");
	hRcpNoFT->GetYaxis()->SetTitle("R_{CP}");
	hRcpSimpleFT->GetYaxis()->SetTitle("R_{CP}");
	hRcp->GetXaxis()->SetTitle("jet_pT_sub");
	hRcpNoFT->GetXaxis()->SetTitle("jet_pT_sub");
	hRcpSimpleFT->GetXaxis()->SetTitle("jet_pT_sub");




	// Drawing
#ifdef DRAW
	unsigned int color[kTotal] = {kMagenta+4, kBlue+4, kViolet+4, kGreen+4, kMagenta+4, kRed+4};
	for(int t = 0; t<kTotal; t++) {
		for(int i = 0; i<Npt; i++) {
			hist[t][i]->SetLineColor(color[t]-i);
			hist[t][i]->SetMarkerColor(color[t]-i);
		}
		hist[t][Npt-1]->SetMarkerStyle(8);
	}
	hist[kPeri4S][Npt-1]->SetMarkerStyle(4);

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TCanvas *c[kTotal+2];
	for(int t = 0; t<kTotal+2; t++) {
		c[t] = new TCanvas();
	}
	TText *text = new TText();
	for(int t = 0; t<kTotal; t++) {		// for each type draw all pt bins together
		c[t]->cd();
		c[t]->SetLogy();
		hist[t][Npt-1]->Draw("HISTp");
		for(int i = 0; i<Npt-1; i++) {
			hist[t][i]->Draw("HISTsame");
		}
		text->DrawTextNDC(0.5,0.95,name[t]);
	}

	// draw pt sum of each type together
	c[kTotal]->cd();
	c[kTotal]->SetLogy();
	hist[kPeri][Npt-1]->Draw("c");
	for(int t = kLBT; t<kTotal; t++) { 
		hist[t][Npt-1]->Draw("csame");
	}


	TLegend *l = new TLegend(0.65,0.7,0.9,0.9);
	for(int t = 0; t<kTotal ; t++) {
		l->AddEntry(hist[t][Npt-1],LegendName[t],"pl");
	}
	//l->AddEntry(hist[kPeri][Npt-1],"Pythia","pl");
	//l->AddEntry(hist[kPeri4S][Npt-1],"Pythia4S","pl");
	//l->AddEntry(hist[kLBT][Npt-1],"LBT w/ f.t.","pl");
	//l->AddEntry(hist[kLBTSimple][Npt-1],"LBT w/ simple f.t.","pl");
	//l->AddEntry(hist[kLBTNoFT][Npt-1],"LBT default","pl");
	l->Draw();


	c[kTotal+1]->cd();
	c[kTotal+1]->SetFrameLineWidth(3);
	c[kTotal+1]->SetLeftMargin(0.15);
	c[kTotal+1]->SetBottomMargin(0.2);

	hRcp->SetLineColor(hist[kLBT][Npt-1]->GetLineColor());
	hRcp->SetMarkerColor(hist[kLBT][Npt-1]->GetMarkerColor());
	hRcp->SetMarkerStyle(hist[kLBT][Npt-1]->GetMarkerStyle());
	hRcp->SetMarkerSize(1.5);
	hRcpNoFT->SetLineColor(hist[kLBTNoFT][Npt-1]->GetLineColor());
	hRcpNoFT->SetMarkerColor(hist[kLBTNoFT][Npt-1]->GetMarkerColor());
	hRcpNoFT->SetMarkerStyle(hist[kLBTNoFT][Npt-1]->GetMarkerStyle());
	hRcpNoFT->SetMarkerSize(1.5);
	hRcpSimpleFT->SetLineColor(hist[kLBTSimple][Npt-1]->GetLineColor());
	hRcpSimpleFT->SetMarkerColor(hist[kLBTSimple][Npt-1]->GetMarkerColor());
	hRcpSimpleFT->SetMarkerStyle(hist[kLBTSimple][Npt-1]->GetMarkerStyle());
	hRcpSimpleFT->SetMarkerSize(1.5);

	hRcp->GetXaxis()->SetNdivisions(505);
	hRcp->GetYaxis()->SetNdivisions(505);
	hRcp->GetYaxis()->SetRangeUser(0,1.5);
	hRcp->GetXaxis()->SetLabelSize(0.06);
	hRcp->GetXaxis()->SetTitleSize(0.06);
	hRcp->GetYaxis()->SetLabelSize(0.06);
	hRcp->GetYaxis()->SetTitleSize(0.06);


	hRcp->Draw();
	hRcpNoFT->Draw("same");
	hRcpSimpleFT->Draw("same");

	TLegend *l2 = new TLegend(0.65,0.7,0.898,0.898);
	l2->AddEntry(hRcp,LegendName[kLBT],"pl");//"LBT w/ f.t.","pl");
	l2->AddEntry(hRcpSimpleFT,LegendName[kLBTSimple],"pl");//"LBT w/ simple f.t.","pl");
	l2->AddEntry(hRcpNoFT,LegendName[kLBTNoFT],"pl");//"LBT no f.t.","pl");
	l2->Draw();

	TLine *line = new TLine(0,1,100,1);
	line->SetLineStyle(2);
	line->Draw("same");
#endif

#ifdef SAVEROOT
		TFile *fout = new TFile("Rcp_FTvsSimplevsNo.root","RECREATE");
		for(int i = 0; i<kTotal; i++) for(int j = 0; j<Npt; j++) hist[i][j]->Write();
		hRcp->Write(); 
		hRcpSimpleFT->Write(); 
		hRcpNoFT->Write();
#ifdef DRAW
		for(int i = 0; i<kTotal+2; i++) c[i]->Write();
#endif
		fout->Close();
#endif

}

