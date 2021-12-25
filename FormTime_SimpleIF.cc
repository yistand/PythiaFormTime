//===========================================================
//
// 2021.09.17
// Li Yi
//
// Modified from FormTime.cc.
// Update with a simpler time calculation:
// from OLD cumulate all splittings 2*E_mother*x*(1-x)/kT^2 
// to NEW initial & final parton 2*E_daughter/kT^2
//
//===========================================================
//
// 2021.08.13
// Li Yi
// This program is to calculate/set parton formation time
// 
// based on PYTHIA example main91.cc
//
//===========================================================
// Original message:
// main91.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It studies the charged multiplicity distribution at the LHC.
// Modified by Rene Brun, Axel Naumann and Bernhard Meirose
// to use ROOT for histogramming.

// Stdlib header file for input and output.
#include <iostream>

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"

// ROOT, for histogramming.
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"

#define MINIMALPT 0.1		// if pTHat_Min is too small, do softQCD & cut off on pTHat_Max. just use a value which is small.

using namespace Pythia8;

int verbose = 0;		// Print More information (1), less information (0)
int warning = 0;		// Print warning msg (1), no warning msg (0)


double CalFormTime_SimpleIF(Vec4 mother, Vec4 daughter) {
	double x = dot3(mother,daughter)/mother.pAbs2();
	double kt = cross3(mother,daughter).pAbs()/mother.pAbs();
	double time = 0; 
	if(kt>0 || kt<0) {
		time = 2.*daughter.e()/pow(kt,2);
	}
	if(verbose) cout<<"mother e="<<mother.e()<<" daughter e="<<daughter.e()<<" x="<<x<<" kt="<<kt<<" time="<<time;
	return time;
}



// ===================== Particle species accepted by LBT (u,d,s,c,g) ==========
bool Accept4LBT(int id) {
	const int Nok = 9;
	int acceptedId[Nok] = {1, -1, 2, -2, 3, -3, 4, -4, 21};
	for(int i = 0; i<Nok; i++) {
		if(id == acceptedId[i]) {
			return true;
		}
	}
	return false;
}


//======================  Main Program  ==========================
int main(int argc, char* argv[]) {

	// Create the ROOT application environment.
	//TApplication theApp("hist", &argc, argv);

	int Nevents = 1;
	int pTHat_Min = 15;
	int pTHat_Max = 20;

	if(argc==4) {
		Nevents = atoi(argv[1]);
		pTHat_Min = atoi(argv[2]);
		pTHat_Max = atoi(argv[3]);
	}

	int JobId = -1;		// If not -1, will update output .dat name for muliple jobs & use this as random seed

	if(argc==5) {		// use this only for job on cluster 
		Nevents = atoi(argv[1]);
		pTHat_Min = atoi(argv[2]);
		pTHat_Max = atoi(argv[3]);
		JobId = atoi(argv[4]);
	}

	cout<<Nevents<<" events with pT: "<<pTHat_Min<<" - "<<pTHat_Max<<endl<<endl;
	if(argc==5) {		// if on cluster
		cout<<"JobId as random seed = "<<JobId<<endl;
	}
	

	// Create Pythia instance and set it up to generate hard QCD processes
	// pTHat = 20 GeV - 25 GeV for pp collisions at 200 GeV.
	// Stop the generation after the hard process and parton-level activity has been generated, but before the hadron-level steps.
	Pythia pythia;

	pythia.readString("Beams:eCM = 200.");
	pythia.readString("Beams:idA = 2212");
	pythia.readString("Beams:idB = 2212");

	char StringName[400];
	if(pTHat_Min<MINIMALPT) {		// do softqcd & cut off at pTHat_Max
		pythia.readString("SoftQCD:nonDiffractive = on");
		pythia.readString("HardQCD:all = off");
	}
	else {
		pythia.readString("SoftQCD:nonDiffractive = off");
		pythia.readString("HardQCD:all = on");
		sprintf(StringName,"PhaseSpace:pTHatMin = %i",pTHat_Min);
		pythia.readString(StringName);
		sprintf(StringName,"PhaseSpace:pTHatMax = %i",pTHat_Max);
		pythia.readString(StringName);
		//pythia.readString("PhaseSpace:bias2Selection = on");	// Switch on a biased phase space sampling, with compensatingly weighted events, for 2 → 2 processes.
		//pythia.readString("PhaseSpace:bias2SelectionPow = 6");	// A 2 → 2 process at a scale pTHat will be oversampled in phase space by an amount (pTHat/pTRef)^pow, where you set the power pow here. Events are assigned a compensating weight the inverse of this, i.e. Info::weight() will return (pTRef/pTHat)^pow. This weight should then be used in the histogramming of event properties. The final overall normalization also involves the Info::weightSum() value.
		pythia.readString("PromptPhoton:all = on");	// this one is a hard process, should not mix with softqcd
	}

	pythia.readString("PartonLevel:ISR = on");
	pythia.readString("PartonLevel:MPI = on");
	pythia.readString("PartonLevel:FSR = on");

	pythia.readString("HadronLevel:all = off");

	pythia.readString("WeakSingleBoson:all=off");
	pythia.readString("WeakDoubleBoson:all=off");

	pythia.readString("ParticleDecays:limitTau0=on");
	pythia.readString("ParticleDecays:tau0Max = 10");



	// Less screen print out
	if(!verbose) pythia.readString("Print:quiet = on");

	// Set Random seed for replicable result
	if(JobId!=-1) {
		pythia.readString("Random:setSeed = on");
		sprintf(StringName,"Random:seed = %i",JobId);
		pythia.readString(StringName);
	}

	pythia.init();

	// Create file on which histogram(s) can be saved && dat files
	const char *dir = "../JETSCAPE/pythiaInput-SimpleIF/";
	if(JobId==-1) {
		sprintf(StringName,"hist_pTHat%03d_to%03d.root",pTHat_Min,pTHat_Max);
	}
	else 	{
		sprintf(StringName,"%shist_pTHat%03d_to%03d_Job%d.root",dir,pTHat_Min,pTHat_Max,JobId);
	}
	TFile* outFile = new TFile(StringName, "RECREATE");
	ofstream fParton;
	if(JobId==-1) {
		sprintf(StringName,"parton_pTHat%03d_to%03d.dat",pTHat_Min,pTHat_Max);
	}
	else 	{
		sprintf(StringName,"%sparton_pTHat%03d_to%03d_Job%d.dat",dir,pTHat_Min,pTHat_Max,JobId);
	}
	fParton.open(StringName,ios::out);
	cout<<"Output: "<<StringName<<endl;


	// Book histogram.
	TH1F *mult = new TH1F("mult","charged multiplicity", 100, -0.5, 799.5);
	TH1D *hTime = new TH1D("hTime","Formation Time of Final Particles (fm) from Hardest Subprocess",1000,-100,100);
	TH2D *hTimeVspT = new TH2D("hTimeVspT","Formation Time of Final Particles (fm) from Hardest Subprocess vs pT (GeV)",100,0,100,100,0,100);
	TProfile *pTimeVspT = new TProfile("pTimeVspT","Formation Time of Final Particles (fm) from Hardest Subprocess vs pT (GeV)",100,0,100);
	TH2D *hTimeVsJetPt = new TH2D("hTimeVsJetPt","Formation Time of Final Particles from Hardest Subprocess vs its Ancestor's pT",100,0,100,100,0,100);
	TProfile *pTimeVsJetPt = new TProfile("pTimeVsJetPt","Formation Time of Final Particles from Hardest Subprocess vs its Ancestor's pT",100,0,100);

	// Begin event loop. Generate event; skip if generation aborted.
	//for (int iEvent = 0; iEvent < Nevents; ++iEvent) {
	int iEvent = 0;
	for (     ; iEvent < Nevents;     ) {
		if (!pythia.next()) continue;
		if ( pTHat_Min<MINIMALPT && pythia.info.pTHat()>pTHat_Max) continue;	// if softQcd, need to cut off at pTHat_Max


		if(iEvent%1000==0) cout<<"event "<<iEvent<<endl;



		// Number of final partons 
		int Nfparton = 0;
		int nCharged = 0;
		for (int k = 0; k < pythia.event.size(); ++k) {
			if(pythia.event[k].isFinal() && Accept4LBT(pythia.event[k].id())) Nfparton++;		// Final particle && can be accepted by LBT program
			if (pythia.event[k].isFinal() && pythia.event[k].isCharged()) nCharged++;
		}

		// Write dat information
		fParton << iEvent+1 << "\t" << Nfparton << "\t" <<pythia.info.pTHat()<<endl;
	
		// Fill number of charged particle
		mult->Fill( nCharged );



		// =============== Section for Formation Time Calculation ===========
		for (int k = 1; k < pythia.event.size(); ++k) {	// skip 0 for system itself
			if(pythia.event[k].isFinal() && Accept4LBT(pythia.event[k].id())) {	// final partons only, should be accepted by LBT program later

				double time = 0;		// if not from hardest parton, use time = 0
				int ancestor = 0;
				if ( pythia.event[k].isAncestor(5) && pythia.event[5].status()==-23  ) {		// From one of hardest subprocess outgoing parton 5 or 6
					ancestor = 5;
					time  = 0.1973*CalFormTime_SimpleIF(pythia.event[5].p(),pythia.event[k].p());
				}	// End of If originated from Hardest subprocess	
				else if ( pythia.event[k].isAncestor(6) && pythia.event[6].status()==-23 ) {
					ancestor = 6;
					time  = 0.1973*CalFormTime_SimpleIF(pythia.event[6].p(),pythia.event[k].p());

				}	// End of If originated from Hardest subprocess	

				hTime->Fill(time);
				if(ancestor) {
					// Fill histogram
					hTimeVspT->Fill(pythia.event[k].pT(),time);
					pTimeVspT->Fill(pythia.event[k].pT(),time);
					hTimeVsJetPt->Fill(pythia.event[ancestor].pT(),time);
					pTimeVsJetPt->Fill(pythia.event[ancestor].pT(),time);
				}	// End of If flaghard

				// Write dat file
				// if not from hardest, time is 0, still record this parton
				fParton<< k+1 <<"\t"<<pythia.event[k].id()<<"\t"<<pythia.event[k].px()<<"\t"<<pythia.event[k].py()<<"\t"<<pythia.event[k].pz()<<"\t"<<pythia.event[k].e()<<"\t"<<pythia.event[k].m()<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<time<<endl;		
				if(verbose) {
					if(ancestor) cout<<endl<<"Final particle "<<k<<" ("<<pythia.event[k].id()<<") orignated from "<<ancestor<<" formation time = "<<time<<endl;
					else cout<<"Final particle "<<k<<" ("<<pythia.event[k].id()<<") formation time = "<<time<<endl;
				}


			}	// End of If final 
		}	// End of pythia particle loop

		// =============== End of Section for Formation Time Calculation ===========
		++iEvent;

	}	// end for loop events


	// Show histogram. Possibility to close it.
	//mult->Draw();
	//hTime->Draw();
	//std::cout << "\nDouble click on the histogram window to quit.\n";
	//gPad->WaitPrimitive();

	// Save histogram on file and close file.
	mult->Write();
	hTime->Write();
	hTimeVspT->Write();
	pTimeVspT->Write();
        hTimeVsJetPt->Write();
        pTimeVsJetPt->Write();
	delete outFile;
	
	fParton.close();

	// Statistics on event generation.
	//pythia.stat();
	ofstream fXsec;
	if(JobId==-1) {
		sprintf(StringName,"Xsec_pTHat%03d_to%03d.dat",pTHat_Min,pTHat_Max);
	}
	else {
		sprintf(StringName,"%sXsec_pTHat%03d_to%03d.dat",dir,pTHat_Min,pTHat_Max);
	}
	fXsec.open(StringName,std::ofstream::out | std::ofstream::app);
	//fXsec << pythia.info.nAccepted() << '\t' << pythia.info.sigmaGen() << '\t' << pythia.info.sigmaErr() << '\t' << JobId << endl;
	fXsec << iEvent << '\t' << pythia.info.sigmaGen() << '\t' << pythia.info.sigmaErr() << '\t' << JobId << endl;
	fXsec.close();
	
	// Done.
	cout<<"====== Done ======"<<endl;
	return 0;
}
