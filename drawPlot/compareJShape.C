// ================================================================================
//
//		2021.10.08 Li YI
//		Compare jet shapes (pt rho, ntrk yield) vs dr
//		for LBT with different formation time assumption
//
//		modified from compareRaa.C
//		merge histograms in root file for same pt bins
//		add pt bins with their cross section
//		take ratio between LBT/pythia
//
// ================================================================================
//
//		2021.08.24 Li YI
//		caculate and plot Raa for JETSCAPE & JETSCAPE-noFT
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


#define DRAW		// if define draw plots
#define SAVEROOT	// if define save root file

using std::cout; using std::cin;
using std::endl; using std::vector;

int verbose = 0;

// ------------------ Functions ----------------

// Read Each input files & obtain histogram
// nHist for number of histrograms to add
// histName is histogram name
// h is the histogram
bool ReadJetShape(char *filename, const int nHist, const string histName[], TH1D* h[], TDirectory *dir) {
	 Bool_t status = TH1::AddDirectoryStatus();
	 TH1::AddDirectory(kFALSE);

	if(verbose) cout<<"Read "<<filename<<endl;
	TFile *f = new TFile(filename);
	if(!f->IsOpen()) { cout<<"Warning: Cannot open "<<filename <<endl; return false;}

	for(int i = 0 ;i<nHist ; i++ ) {
		TH1D *htmp = (TH1D*)f->Get(histName[i].c_str());
		if(!htmp) { cout<<"Error: Cannot find "<<histName[i]<<endl; return false;}
		if(!h[i])  {				// if this is the first file, read and set the histogram to current directory
			if(verbose) cout<<"Read "<<histName[i]<<endl;
			h[i] = (TH1D*)htmp->Clone(htmp->GetName());
			h[i]->SetDirectory(dir);	
		}
		else {					// this is not the first file, add it to the histogram in first file
			if(verbose) cout<<"Add "<<htmp->GetName()<<" "<<h[i]->GetName()<<endl;
			h[i]->Add(htmp);
		}
	}
	f->Close();

	TH1::AddDirectory(status);

	return true;
}

// List all files in directory
bool ListDir(const char *path, vector<string>& fileList, const char *tag, const char *tag2="", const char *tag3="") {
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
                if((s.find(tag)!=string::npos && s.find(tag2)!=string::npos) && s.find(tag3)!=string::npos) {	// if found filename with tag
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
void compareJShape() {

	enum Type {kLBT, kLBTNoFT, kLBTSimple, kPythia, kPythia4S, kTotal};
	const char *name[kTotal] = {"LBT","LBTNoFT","LBTSimple","Pythia","Pythia4S"};

	enum ShapeHist {krho, kYield, kSpectra, kSumPar, kShapes};
	const string shapeName[kShapes] = {"rho","Yield","pt","Sums"};
	const string histName[kShapes] = {"hPtSumVsdr","hNtrkVsdr","hpt","hSums"};

	// input 
	ifstream input[kTotal]; // kLBT, kLBTNoFT, kLBTSimple, kPythia, kPythia4S;

	const char * dir[kTotal] = {"../JETSCAPE/results/", "../JETSCAPE-noFT/results/", "../JETSCAPE/results-SimpleIF/","../JETSCAPE/results/", "../JETSCAPE/results-SimpleIF/"};

	const int Npt = 14+1;		// +1 for the sum
	int arrayMin[Npt]= {0 , 2 , 4 , 6  , 8 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60, 0};		// last bin is sum
	int arrayMax[Npt]={2 , 4 , 6 , 8 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60  ,100, 100};		// last bin is sum
	TH1D *hist[kTotal][Npt][kShapes]= {{{NULL}}};	

	double Xsec[Npt-1] = {28.4852,9.61573,0.47458,0.0632963,0.0135362,0.00525799,0.00046218,6.42463e-05,1.1666e-05,2.39516e-06,5.34032e-07,1.4737e-07,7.64804e-09,3.3109e-10};
	double Xsec4S[Npt-1] = {28.4852,9.61573,0.47458,0.0632963,0.0135362,0.00525799,0.00046218,6.42463e-05,1.1666e-05,2.39516e-06,5.34032e-07,1.4737e-07,7.64804e-09,3.3109e-10};

	double tmpsec[Npt-1] = {0};
	if(ReadXsec("../JETSCAPE/pythiaInput/",Npt-1,arrayMin,arrayMax,tmpsec)) {	// read xsec file if exists
		for(int i = 0; i<Npt-1; i++) Xsec[i] = tmpsec[i];
	}
	if(ReadXsec("../JETSCAPE/pythiaInput-SimpleIF/",Npt-1,arrayMin,arrayMax,tmpsec)) {	// read xsec file if exists
		for(int i = 0; i<Npt-1; i++) Xsec4S[i] = tmpsec[i];
	}

	double SumXsec=0;
	for(int i = 0; i<Npt-1; i++) SumXsec+=Xsec[i];
	double SumXsec4S=0;
	for(int i = 0; i<Npt-1; i++) SumXsec4S+=Xsec4S[i];

	if(1) {
		cout<<"ptHat\tCross section"<<endl;
		for(int i = 0; i<Npt-1; i++){
			cout<<arrayMin[i]<<'-'<<arrayMax[i]<<'\t'<<Xsec[i]<<'\t'<<Xsec4S[i]<<endl;
		}
		cout<<"sum = "<<SumXsec<<'\t'<<SumXsec4S<<endl;
		cout<<endl;
	}

	const char *tag[kTotal] = {"LBTjet", "LBTjet", "LBTjet", "Pythiajet", "Pythiajet"};


	double sXsecNjets[kTotal] = {0};	// cross section * nJets (from hSums histogram) for all pT
	double sXsecPtSum[kTotal] = {0};	// cross section * PtSums (from hSums histogram) for all pT
	double sXsecNjets4S[kTotal] = {0};	// cross section * nJets (from hSums histogram) for all pT. simple-IF case.
	double sXsecPtSum4S[kTotal] = {0};	// cross section * PtSums (from hSums histogram) for all pT. simple-IF case.
	double sXsecNjetsPerEvt[kTotal] = {0};	// cross section * nJets / Nevents (from hSums histogram) for all pT
	double sXsecPtSumPerEvt[kTotal] = {0};	// cross section * PtSums / Nevents (from hSums histogram) for all pT
	double sXsecNjets4SPerEvt[kTotal] = {0};	// cross section * nJets / Nevents (from hSums histogram) for all pT. simple-IF case.
	double sXsecPtSum4SPerEvt[kTotal] = {0};	// cross section * PtSums / Nevents (from hSums histogram) for all pT. simple-IF case.
	// Loop over pt bins
	for(int i = 0; i<Npt-1; i++) {	// the last bin for sum over all pt bins
		vector<string> list[kTotal];
		for(int t = 0; t<kTotal; t++) {
			char ctmp[500] = "";
			sprintf(ctmp,"pTHat%03d_to%03d",arrayMin[i],arrayMax[i]);
			if(verbose) cout<<"Finding "<<ctmp<<endl;
			ListDir(dir[t], list[t], tag[t], ctmp, ".root");	// root file only

			if(list[t].size()==0) {
				cout << "WARNING: Find no file for " << ctmp << endl;
				continue;
			}
			for(auto it : list[t]) {	// for each files in list
				char ifilename[500];
				sprintf(ifilename,"%s%s",dir[t],it.data()); 
				TDirectory * current_directory = gDirectory->CurrentDirectory();
				ReadJetShape(ifilename, kShapes, histName, hist[t][i], current_directory);
			}	// End of Loop over list of files for type t & pt i
			for(int s = 0; s<kShapes; s++) {
				hist[t][i][s]->SetName(Form("h%s%s_pt%d",shapeName[s].c_str(),name[t],i));
				hist[t][i][s]->SetTitle(Form("%s %s %s",shapeName[s].c_str(),name[t],ctmp));
			}

			// Scale by jets & ptsum
			double i_nevent = hist[t][i][kSumPar]->GetBinContent(1);	// number of events 
			double i_njets = hist[t][i][kSumPar]->GetBinContent(2);	// number of jets 
			double i_ptsums = hist[t][i][kSumPar]->GetBinContent(5);	// sum of all particle pt in jet cone
			if(1) cout<<"type="<<t<<" pt="<<i<<" njet = "<<i_njets<<" ptsum = "<<i_ptsums;
			if(t==kPythia4S || t==kLBTSimple) {
				sXsecNjets4S[t]+=i_njets*Xsec4S[i];
				sXsecPtSum4S[t]+=i_ptsums*Xsec4S[i];
				if(i_nevent>0) {
					sXsecNjets4SPerEvt[t]+=i_njets*Xsec4S[i]/i_nevent;
					sXsecPtSum4SPerEvt[t]+=i_ptsums*Xsec4S[i]/i_nevent;
				}
			}
			else {
				sXsecNjets[t]+=i_njets*Xsec[i];
				sXsecPtSum[t]+=i_ptsums*Xsec[i];
				if(i_nevent>0) {
					sXsecNjetsPerEvt[t]+=i_njets*Xsec[i]/i_nevent;
					sXsecPtSumPerEvt[t]+=i_ptsums*Xsec[i]/i_nevent;
				}
			}
			//if(i_njets>0) {
			//	hist[t][i][kSpectra]->Scale(1./i_njets);

			//	hist[t][i][kYield]->Scale(1./i_njets);
			//	hist[t][i][kYield]->Scale(1,"width");

			//	hist[t][i][kYield]->GetYaxis()->SetTitle(Form("%s#frac{1}{N_{jets}#delta#Deltar}",hist[t][i][kYield]->GetYaxis()->GetTitle()));
			//}
			//if(i_ptsums>0) {
			//	hist[t][i][krho]->Scale(1./i_ptsums);
			//	hist[t][i][krho]->Scale(1,"width");

			//	hist[t][i][krho]->GetYaxis()->SetTitle(Form("%s#frac{1}{#Sigma_{jets}#Sigma_{trk}p_{T}^{trk}#delta#Deltar}",hist[t][i][krho]->GetYaxis()->GetTitle()));
			//}

			// Scale by cross section
			if(1) cout<<" xsec4s="<<Xsec4S[i]<<" per_xsec4s="<<Xsec4S[i]/SumXsec4S<<" xsec="<<Xsec[i]<<" per_xsec="<<Xsec[i]/SumXsec<<endl;
			if(t==kPythia4S || t==kLBTSimple) {
				//hist[t][i][kSpectra]->Scale(Xsec4S[i]);
				//hist[t][i][kYield]->Scale(Xsec4S[i]/SumXsec4S);
				//hist[t][i][krho]->Scale(Xsec4S[i]/SumXsec4S);
				if(i_nevent>0) {
					hist[t][i][kSpectra]->Scale(Xsec4S[i]/i_nevent);
					hist[t][i][kYield]->Scale(Xsec4S[i]/i_nevent);
					hist[t][i][krho]->Scale(Xsec4S[i]/i_nevent);
				}
			}
			else {
				//hist[t][i][kYield]->Scale(Xsec[i]/SumXsec);
				//hist[t][i][krho]->Scale(Xsec[i]/SumXsec);
				if(i_nevent>0) {
					hist[t][i][kSpectra]->Scale(Xsec[i]/i_nevent);
					hist[t][i][kYield]->Scale(Xsec[i]/i_nevent);
					hist[t][i][krho]->Scale(Xsec[i]/i_nevent);
				}
			}
		}	// End of Loop over type t
	}	// End of Loop over pt i

	if(verbose) cout<<"Sum over pt bins"<<endl;
	for(int t = 0; t<kTotal; t++) {
		for(int s = 0; s<kShapes-1; s++) {	// not add kSumPar
			for(int i = 0; i<Npt-1; i++) {
				if(hist[t][i][s] && hist[t][i][s]->GetEntries()>10)   {	// not too few entries to avoid large fluctuation
					if(!hist[t][Npt-1][s]) {
						hist[t][Npt-1][s] = (TH1D*) hist[t][i][s]->Clone(Form("h%s%s_pt%d",shapeName[s].c_str(),name[t],Npt-1));
						hist[t][Npt-1][s]->SetTitle(Form("%s %s",shapeName[s].c_str(),name[t]));
					}
					else hist[t][Npt-1][s]->Add(hist[t][i][s]);	// the sum of all pt bins
					if(verbose) cout<<" add "<<hist[t][i][s]->GetName()<<" type="<<t<<" shape="<<s<<" pt="<<i<<endl;
				}  // End of if hist exists
			}  // End of pt for loop
		}  // End of kShapes for loop
		if(t==kPythia4S || t==kLBTSimple) {
			if(sXsecNjets4SPerEvt[t]>0) {
				//hist[t][Npt-1][kYield]->Scale(1./sXsecNjets4S[t]);
				hist[t][Npt-1][kYield]->Scale(1./sXsecNjets4SPerEvt[t]);
				hist[t][Npt-1][kYield]->Scale(1,"width");
				//hist[t][Npt-1][kYield]->GetYaxis()->SetTitle(Form("%s#frac{1}{N_{jets}#delta#Deltar}",hist[t][Npt-1][kYield]->GetYaxis()->GetTitle()));
			}

			if(sXsecPtSum4SPerEvt[t]>0) {
				//hist[t][Npt-1][krho]->Scale(1./sXsecPtSum4S[t]);
				hist[t][Npt-1][krho]->Scale(1./sXsecPtSum4SPerEvt[t]);
				hist[t][Npt-1][krho]->Scale(1,"width");
				//hist[t][Npt-1][krho]->GetYaxis()->SetTitle(Form("%s#frac{1}{#Sigma_{jets}#Sigma_{trk}p_{T}^{trk}#delta#Deltar}",hist[t][Npt-1][krho]->GetYaxis()->GetTitle()));
			}
		}	// End of if simple-IF
		else {
			if(sXsecNjetsPerEvt[t]>0) {
				hist[t][Npt-1][kYield]->Scale(1./sXsecNjetsPerEvt[t]);
				hist[t][Npt-1][kYield]->Scale(1,"width");
				//hist[t][Npt-1][kYield]->GetYaxis()->SetTitle(Form("%s#frac{1}{N_{jets}#delta#Deltar}",hist[t][Npt-1][kYield]->GetYaxis()->GetTitle()));
			}

			if(sXsecPtSumPerEvt[t]>0) {
				hist[t][Npt-1][krho]->Scale(1./sXsecPtSumPerEvt[t]);
				hist[t][Npt-1][krho]->Scale(1,"width");
				//hist[t][Npt-1][krho]->GetYaxis()->SetTitle(Form("%s#frac{1}{#Sigma_{jets}#Sigma_{trk}p_{T}^{trk}#delta#Deltar}",hist[t][Npt-1][krho]->GetYaxis()->GetTitle()));
			}
		}	// End if

		hist[t][Npt-1][kSpectra]->GetYaxis()->SetTitle("Cross section");
		hist[t][Npt-1][kYield]->GetYaxis()->SetTitle(Form("%s#frac{1}{N_{jets}#delta#Deltar}",hist[t][Npt-1][kYield]->GetYaxis()->GetTitle()));
		hist[t][Npt-1][krho]->GetYaxis()->SetTitle(Form("%s#frac{1}{#Sigma_{jets}#Sigma_{trk}p_{T}^{trk}#delta#Deltar}",hist[t][Npt-1][krho]->GetYaxis()->GetTitle()));
	}  // End of kTotal for loop



	// Raa & JetShapes
	TH1D *hratio[kTotal-2][kShapes-1]={{NULL}};	// krho, kYield, kSpectra; kLBT, kLBTNoFT, kLBTSimple 
	for(int s = 0; s<kShapes-1; s++) {
		for(int t = 0; t<kTotal-2; t++) {
			hratio[t][s] = (TH1D*)hist[t][Npt-1][s]->Clone(Form("hRatio%s%s",shapeName[s].c_str(),name[t]));
			hratio[t][s]->SetTitle(Form("%s_{AA/pp} %s",shapeName[s].c_str(), name[t]));
			if(verbose) cout<<"Create ratio for "<<hratio[t][s]->GetName()<<endl;
			if( t==kLBT || t==kLBTNoFT ) {
				hratio[t][s]->Divide(hist[kPythia][Npt-1][s]);
			}
			else if(t==kLBTSimple) {
				hratio[t][s]->Divide(hist[kPythia4S][Npt-1][s]);
			}
			hratio[t][s]->GetYaxis()->SetTitle(Form("%s_{AA/pp}",shapeName[s].c_str()));
			if(s==kSpectra) {
				hratio[t][s]->GetXaxis()->SetTitle("jet_pT_sub");
			}
			else {
				hratio[t][s]->GetXaxis()->SetTitle("#Deltar");
			}
		}
	}



	// Drawing
#ifdef DRAW
	if(verbose) cout<<"Drawing"<<endl;
	unsigned int color[kTotal] = {kBlue+4, kGreen+4, kRed+4, kMagenta+4, kMagenta+4};
	for(int t = 0; t<kTotal; t++) {
		for(int s = 0; s<kShapes-1; s++) {
			for(int i = 0; i<Npt; i++) {
				hist[t][i][s]->SetLineColor(color[t]-i);
				hist[t][i][s]->SetMarkerColor(color[t]-i);
			}
			hist[t][Npt-1][s]->SetMarkerStyle(8);
		}
	}
	for(int s = 0; s<kShapes-1; s++) {
		hist[kPythia4S][Npt-1][s]->SetMarkerStyle(4);
	}

	gStyle->SetOptStat(0);
	gStyle->SetOptTitle(0);

	TCanvas *c[kTotal+2][kShapes-1];
	for(int t = 0; t<kTotal+2; t++) {
		for(int s = 0; s<kShapes-1; s++) {
			c[t][s] = new TCanvas();
		}
	}
	TText *text = new TText();
	for(int t = 0; t<kTotal; t++) {		// for each type draw all pt bins together
		for(int s = 0; s<kShapes-1; s++) {
			c[t][s]->cd();
			if(s==kSpectra) c[t][s]->SetLogy();
			hist[t][Npt-1][s]->Draw("HISTp");
			for(int i = 0; i<Npt-1; i++) {
				hist[t][i][s]->Draw("HISTsame");
			}
			text->DrawTextNDC(0.5,0.95,Form("%s %s", shapeName[s].c_str(), name[t]));
		}
	}

	// draw pt bin sum of each type together for each kShape
	TLegend *leg[kShapes-1];
	TLegend *leg2[kShapes-1];
	for(int s = 0; s<kShapes-1; s++) {
		c[kTotal][s]->cd();
		c[kTotal][s]->SetLogy();
		hist[kPythia][Npt-1][s]->Draw("c");
		for(int t = kLBT; t<kTotal; t++) { 
			hist[t][Npt-1][s]->Draw("csame");
		}


		leg[s] = new TLegend(0.65,0.7,0.9,0.9);
		leg[s]->AddEntry(hist[kPythia][Npt-1][s],"Pythia","pl");
		leg[s]->AddEntry(hist[kPythia4S][Npt-1][s],"Pythia4S","pl");
		leg[s]->AddEntry(hist[kLBT][Npt-1][s],"LBT w/ f.t.","pl");
		leg[s]->AddEntry(hist[kLBTSimple][Npt-1][s],"LBT w/ simple f.t.","pl");
		leg[s]->AddEntry(hist[kLBTNoFT][Npt-1][s],"LBT no f.t.","pl");
		leg[s]->Draw();


		c[kTotal+1][s]->cd();
		for(int t=0; t<kTotal-2; t++) {
			hratio[t][s]->SetLineColor(hist[t][Npt-1][s]->GetLineColor());
			hratio[t][s]->SetMarkerColor(hist[t][Npt-1][s]->GetMarkerColor());
			hratio[t][s]->SetMarkerStyle(hist[t][Npt-1][s]->GetMarkerStyle());
		}
		hratio[kLBT][s]->Draw();
		for(int t=kLBTNoFT; t<kTotal-2; t++) {
			hratio[t][s]->Draw("same");
		}

		leg2[s] = new TLegend(0.65,0.7,0.9,0.9);
		leg2[s]->AddEntry(hratio[kLBT][s],"LBT w/ f.t.","pl");
		leg2[s]->AddEntry(hratio[kLBTSimple][s],"LBT w/ simple f.t.","pl");
		leg2[s]->AddEntry(hratio[kLBTNoFT][s],"LBT no f.t.","pl");
		leg2[s]->Draw();
	}
#endif

#ifdef SAVEROOT
	if(verbose) cout<<"Saving root file."<<endl;
	TFile *fout = new TFile("out.root","RECREATE");
	for(int i = 0; i<kTotal; i++) for(int j = 0; j<Npt; j++) for(int s=0; s<kShapes; s++) {if(hist[i][j][s]) hist[i][j][s]->Write();}
	for(int i = 0; i<kTotal; i++) {
		for(int s=0; s<kShapes-1; s++) {
			if(hratio[i][s]) hratio[i][s]->Write();
		}
	}
#ifdef DRAW
	for(int i = 0; i<kTotal+2; i++) for(int s=0; s<kShapes-1; s++) c[i][s]->Write();
#endif
	fout->Close();
#endif

	cout<<".. End .."<<endl;
}

