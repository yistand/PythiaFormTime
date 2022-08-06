//================================================================================ 
//
//	2022.08.06 Li YI
//	just put 3 different Rcp simulations together & compare to data
//
//================================================================================ 

void draw3Rcp(int Rmode=0) {
        const int NR=3;
        const char *Rtag[NR] = {"","-R03","-R04"};
        const char *RLegName[NR] = {"R=0.2","R=0.3","R=0.4"};


	const int N = 3;
	const char *fileNameTag[N] = {"-zeroFT", "-SimpleIF",""};

	const char *alphaS[N] = {"030","030","050"};
	const char *legName[N] = {"Setup 1 #alpha_{S}=0.3",
				"Setup 2 #alpha_{S}=0.3",
				"Setup 3 #alpha_{S}=0.5"};

	// Read Rcp LBT
	TFile *f[N];
	TH1D *h[N];
	for(int i = 0; i<N; i++) {
		f[i] = new TFile(Form("Rcp%s%s_alphaS.root",fileNameTag[i],Rtag[Rmode]));
		h[i] = (TH1D*)f[i]->Get(Form("hRcpPeri%s",alphaS[i]));
	}

	// Get Rcp Data
	TGraphAsymmErrors *gr;
	TCanvas *c_tmp = (TCanvas*)f[N-1]->Get("c1_n10");
	gr = (TGraphAsymmErrors*)c_tmp->FindObject("Graph");


	// Set histogram style
	unsigned int style[N] = {20, 21, 47};
        Int_t ci = TColor::GetFreeColorIndex();
        TColor *mycolor0 = new TColor(ci,0.46,0.37,0.73);
        TColor *mycolor1 = new TColor(ci+1,0.9,0.78,0.96);
        TColor *mycolor2 = new TColor(ci+2,0.35,0.77,0.66);
        TColor *mycolor3 = new TColor(ci+3,0.82,0.97,0.8);

	for(int i = 0; i<N; i++) {
		h[i]->SetMarkerStyle(style[i]);
		h[i]->SetMarkerColor(ci+i);
		h[i]->SetLineColor(ci+i);
	}

	// Draw
	gStyle->SetOptStat(0);
        gStyle->SetOptTitle(0);

	TCanvas *c = new TCanvas();
	c->SetFrameLineWidth(3);
        c->SetLeftMargin(0.15);
        c->SetBottomMargin(0.2);

	h[0]->GetXaxis()->SetTitle("Jet p_{T}");
	h[0]->Draw("p");
	for(int i = 1; i<N; i++) {
		h[i]->Draw("psame");
	}
	gr->Draw("psame");

	TLegend *l = new TLegend(0.65,0.7,0.898,0.898);
	for(int i = 0; i<N; i++) {
		l->AddEntry(h[i],legName[i],"pl");
	}
	l->SetHeader(RLegName[Rmode]);
	l->Draw("same");


        TLine *line = new TLine(0,1,100,1);
        line->SetLineStyle(2);
        line->Draw("same");
}
