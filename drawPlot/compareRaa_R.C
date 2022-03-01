// ================================================================================
//
//		2022.02.26 Li YI
//		caculate and plot Raa for different jet R
//
//		to run:
//		root -l
//		.L compareRaa_R.C+
//		compareRaa()		// or with args:
//        				// LBTmode = 0: LBT with FT;    1: LBT SimpleIF;        2: LBT noFT
//					// alphasFlag = 0: 0.15;        1: 000; 2:030;  3:050
//		.q
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
#define PDF		// if define as draw plots as pdf
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
		if(line.empty()) { cout<<"WARNING!! skip Empty line in FILE "<<filename<<endl; continue;}
		for(std::string s; iss >> s; ) 	{
			try{
				if(verbose == 10) cout<<s<<" -> ";
				v.push_back(stod(s)); 
				if(verbose == 10) cout<<stod(s)<<endl;
				} catch (std::invalid_argument const& ex) {
					cout<<s<<" in FILE "<<filename<<" at LINE "<<line<<" after reading "<<Nevent<<" events"<<endl;
				        std::cout << "ERROR!!!!!!! " << ex.what() << '\n';
			}
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
	if(verbose) cout<<"Finish read "<<filename<<endl;
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
void compareRaa(int LBTmode=0, int alphasFlag=0) {
	// LBTmode = 0: LBT with FT;	1: LBT SimpleIF;	2: LBT noFT
	// alphasFlag = 0: 0.15;	1: 000;	2:030;	3:050
	if(LBTmode>2||alphasFlag>3||LBTmode<0||alphasFlag<0) {
		cout<<"ERROR: WRONG input arguments!!!"<<endl;
		cout<<"LBTmode = 0: LBT with FT;    1: LBT SimpleIF;        2: LBT noFT"<<endl;
		cout<<"alphasFlag = 0: 0.15;        1: 000; 2:030;  3:050"<<endl;
		return;
	}

	cout<<"LBTmode = "<<LBTmode<<" alphasFlag = "<<alphasFlag<<endl;
	
	enum Type {kPythiaR2, kPythiaR3, kPythiaR4, kLBTR2, kLBTR3, kLBTR4, kTotal};
	const char *name[kTotal] = {"PythiaR2","PythiaR3","PythiaR4","LBTR2","LBTR3","LBTR4"};
	const char *tag[kTotal] = {"Pythia-jet-g.dat", "Pythia-jet-g-R03","Pythia-jet-g-R04","LBT-jet-g.dat", "LBT-jet-g-R03", "LBT-jet-g-R04"};
	const char *LegendName[kTotal] = {"Pythia R=0.2","Pythia R=0.3","Pythia R=0.4","LBT R=0.2","LBT R=0.3","LBT R=0.4"};

	// ---#LY to be continue ----
	// input 
	//const char * dirFT[2] = {"../JETSCAPE/results/","../JETSCAPE/results-alphaS000/","../JETSCAPE/results/","../JETSCAPE/results-alphaS030/", "../JETSCAPE/results-alphaS050/"};
	//const char * dirSimpleIF[2] = {"../JETSCAPE/results-SimpleIF/","../JETSCAPE/results-SimpleIF-alphaS000/","../JETSCAPE/results-SimpleIF/","../JETSCAPE/results-SimpleIF-alphaS030/","../JETSCAPE/results-SimpleIF-alphaS050/"};	// Simple
	//const char * dirNoFT[2] = {"../JETSCAPE/results/","../JETSCAPE-noFT/results-alphaS000/","../JETSCAPE-noFT/results/","../JETSCAPE-noFT/results-alphaS030/", "../JETSCAPE-noFT/results-alphaS050/"};	// noFT
	
	const char *tagLBTMode[3] = {"","-SimpleIF","-noFT"};
	const char *tagAlphaS[4] = {"","-alphaS000","-alphaS030","-alphaS050"};
	char dirPythia[500] = "../JETSCAPE/results/";
	char dirLBT[500] = "../JETSCAPE/results/";
        if(LBTmode==0) {
                sprintf(dirPythia,"../JETSCAPE/results/");
                sprintf(dirLBT,"../JETSCAPE/results%s/",tagAlphaS[alphasFlag]);
        }
        else if(LBTmode==1) {
                sprintf(dirPythia,"../JETSCAPE/results%s/",tagLBTMode[LBTmode]);
                sprintf(dirLBT,"../JETSCAPE/results%s%s/",tagLBTMode[LBTmode],tagAlphaS[alphasFlag]);
        }
        else if(LBTmode==2) {
                sprintf(dirPythia,"../JETSCAPE/results/");
                sprintf(dirLBT,"../JETSCAPE%s/results%s/",tagLBTMode[LBTmode],tagAlphaS[alphasFlag]);       
        }
	cout<<dirPythia<<endl;
	cout<<dirLBT<<endl;

	ifstream input[kTotal]; 

	const int Npt = 14+1;		// +1 for the sum
	int arrayMin[Npt]= {0 , 2 , 4 , 6  , 8 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60, 0};		// last bin is sum
	int arrayMax[Npt]={2 , 4 , 6 , 8 , 10 , 15 , 20 , 25 , 30 , 35 , 40 , 50 , 60  ,100, 100};		// last bin is sum
	TH1D *hist[kTotal][Npt];	

	double Xsec[Npt-1] = {28.4852,9.61467,0.472347,0.0635027,0.0134543,0.00527283,0.000456384,6.42e-05,1.16125e-05,2.39501e-06,5.32087e-07,1.49413e-07,7.60858e-09,3.33557e-10};

	if(0) {
		double tmpsec[Npt-1] = {0};
		if(ReadXsec("../JETSCAPE/pythiaInput/",Npt-1,arrayMin,arrayMax,tmpsec)) {	// read xsec file if exists
			for(int i = 0; i<Npt-1; i++) Xsec[i] = tmpsec[i];
		}
	}

	//double SumXsec=0;
	//for(int i = 0; i<Npt-1; i++) SumXsec+=Xsec[i];


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
			char *dir;
			if(t==kPythiaR2||t==kPythiaR3||t==kPythiaR4) {
				dir = dirPythia;
			}
			else {
				dir = dirLBT;
			}
			//ListDir(dir[t], tag[t], list[t], ctmp);
			ListDir(dir, tag[t], list[t], ctmp);

			double Nevt = 0;	// total number of events for type t & pt i
			for(auto it : list[t]) {	// for each files in list
				char ifilename[500];
				double iNevt = 0; 
				//sprintf(ifilename,"%s%s",dir[t],it.data()); 
				sprintf(ifilename,"%s%s",dir,it.data()); 
				ReadJet(ifilename, hist[t][i], iNevt);
				Nevt=+iNevt;
			}	// End of Loop over list of files for type t & pt i
			if(Nevt>0) hist[t][i]->Scale(Xsec[i]/Nevt);	// Normalized by number of events for type k & pt i & cross section
		}	// End of Loop over type t
	}	// End of Loop over pt i

	for(int t = 0; t<kTotal; t++) {
		for(int i = 0; i<Npt-1; i++) {
			hist[t][Npt-1]->Add(hist[t][i]);	// the sum of all pt bins
			//cout<<"Add "<<hist[t][i]->GetName()<<" to "<<hist[t][Npt-1]<<endl;
		}
	}


	// Raa
	TH1D *hRaa[kLBTR4-kPythiaR4];
	for(int t=0; t<kLBTR4-kPythiaR4; t++) {
		hRaa[t] = (TH1D*)hist[kLBTR2+t][Npt-1]->Clone(Form("hRaa%s",name[t]));
		hRaa[t]->SetTitle(Form("R_{AA} %s",name[t]));
		hRaa[t]->Divide(hist[t][Npt-1]);
		hRaa[t]->GetYaxis()->SetTitle("R_{AA}");
		hRaa[t]->GetXaxis()->SetTitle("jet_pT_sub");
	}




	// Drawing
#ifdef DRAW
	unsigned int color[6] = {kMagenta+4, kBlue+4, kOrange+10, kPink+7, kGreen+4, kCyan};
	for(int t = 0; t<kTotal; t++) {
		for(int i = 0; i<Npt; i++) {
			hist[t][i]->SetLineColor(color[t]-i);
			hist[t][i]->SetMarkerColor(color[t]-i);
		}
		hist[t][Npt-1]->SetMarkerStyle(8);
	}

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
	hist[kPythiaR2][Npt-1]->Draw("c");
	for(int t = 1; t<kTotal; t++) { 
		hist[t][Npt-1]->Draw("csame");
	}


	TLegend *l = new TLegend(0.65,0.7,0.9,0.9);
	for(int t = 0; t<kTotal ; t++) {
		l->AddEntry(hist[t][Npt-1],LegendName[t],"pl");
	}
	l->Draw();


	c[kTotal+1]->cd();
	c[kTotal+1]->SetFrameLineWidth(3);
	c[kTotal+1]->SetLeftMargin(0.15);
	c[kTotal+1]->SetBottomMargin(0.2);

	for(int t = 0; t<kLBTR4-kPythiaR4 ; t++) {
		hRaa[t]->SetLineColor(hist[t+kPythiaR2][Npt-1]->GetLineColor());
		hRaa[t]->SetMarkerColor(hist[t+kPythiaR2][Npt-1]->GetMarkerColor());
		hRaa[t]->SetMarkerStyle(hist[t+kPythiaR2][Npt-1]->GetMarkerStyle());
		hRaa[t]->SetMarkerSize(1.5);

		hRaa[t]->GetXaxis()->SetNdivisions(505);
		hRaa[t]->GetYaxis()->SetNdivisions(505);
		hRaa[t]->GetYaxis()->SetRangeUser(0,1.5);
		hRaa[t]->GetXaxis()->SetLabelSize(0.06);
		hRaa[t]->GetXaxis()->SetTitleSize(0.06);
		hRaa[t]->GetYaxis()->SetLabelSize(0.06);
		hRaa[t]->GetYaxis()->SetTitleSize(0.06);
	}

	hRaa[0]->Draw();
	for(int t = 1; t<kLBTR4-kPythiaR4 ; t++) {
		hRaa[t]->Draw("same");
	}

	TLegend *l2 = new TLegend(0.65,0.7,0.898,0.898);
	for(int t = 0; t<kLBTR4-kPythiaR4 ; t++) {
		l2->AddEntry(hRaa[t],LegendName[t+kLBTR2],"pl");
	}
	l2->Draw();

	TLine *line = new TLine(0,1,100,1);
	line->SetLineStyle(2);
	line->Draw("same");
#endif
#ifdef SAVEROOT
	TFile *fout = new TFile(Form("Raa%s%s_Rdep.root",tagLBTMode[LBTmode],tagAlphaS[alphasFlag]),"RECREATE");
	for(int i = 0; i<kTotal; i++) for(int j = 0; j<Npt; j++) hist[i][j]->Write();
	for(int t = 0; t<kLBTR4-kPythiaR4 ; t++)  hRaa[t]->Write();
#ifdef DRAW
	for(int i = 0; i<kTotal+2; i++) c[i]->Write();

#ifdef PDF
	Int_t ci = TColor::GetFreeColorIndex();
	TColor *mycolor0 = new TColor(ci,0.8,0.24,0.3);
	TColor *mycolor1 = new TColor(ci+1,0.24,0.53,0.15);
	TColor *mycolor2 = new TColor(ci+2,0.20,0.53,0.97);

	for(int t = 0; t<kLBTR4-kPythiaR4 ; t++)  {
		hRaa[t]->SetLineColor(ci+t);
		hRaa[t]->SetMarkerColor(ci+t);
	}
	c[kTotal+1]->Update();
	c[kTotal+1]->SaveAs(Form("Raa%s%s_Rdep.pdf",tagLBTMode[LBTmode],tagAlphaS[alphasFlag]));

#endif

#endif
	fout->Close();
#endif

}

