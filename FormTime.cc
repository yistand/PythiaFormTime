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

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"

using namespace Pythia8;

int verbose = 0;		// Print More information (1), less information (0)
int warning = 0;		// Print warning msg (1), no warning msg (0)

//Initializing Structure: myParticle (no, formation time, integrated formation time) 
struct myParticle {
	int no;
	double time;
	double totalTime;
};

void PrintMyParticle(vector<myParticle>::iterator mp) {
	cout<<"no="<<mp->no<<" time="<<mp->time<<" totalTime="<<mp->totalTime<<endl;
}

double CalFormTimeFrom4Da(Vec4 daughter1, Vec4 daughter2) {	// calculate formation time from two daughters. mother are reconstructed by two daughters
	Vec4 mother = daughter1+daughter2;
	//double x = dot3(mother,daughter1)/mother.pAbs2();
	double x = 1.*daughter1.e()/mother.e();
	double kt = 1.*cross3(mother,daughter1).pAbs()/mother.pAbs();
	double time = 0; 
	if(kt>0 || kt<0) {
		time = 2.*mother.e()*x*(1.-x)/pow(kt,2);
	}
	if(warning) {
		if(x>1 || (time > 1000 && fabs(mother.eta())<2)) {
			if(x>1) cout<<"WARNING!!! Mother energy less than daughter"<<endl;
			if(time>1000) cout<<"WARNING!!! Long time at mid-rapidity"<<endl;
			cout<<" mother: "<<mother<<" daughter1: "<<daughter1<<" daugther2: "<<daughter2;
			cout<<"mother x daughter1 = "<<cross3(mother,daughter1).pAbs()<<" mother.pAbs = "<<mother.pAbs()<<endl;
			cout<<"x="<<x<<" kt="<<kt<<" time="<<time<<endl<<endl;
		}
	}
	return time;
}


double CalFormTime(Vec4 mother, Vec4 daughter) {
	double x = dot3(mother,daughter)/mother.pAbs2();
	double kt = cross3(mother,daughter).pAbs()/mother.pAbs();
	double time = 0; 
	if(kt>0 || kt<0) {
		time = 2.*mother.e()*x*(1.-x)/pow(kt,2);
	}
	if(verbose) cout<<"mother e="<<mother.e()<<" daughter e="<<daughter.e()<<" x="<<x<<" kt="<<kt<<" time="<<time;
	return time;
}

// Tracing shower develop using mother-daughter relationship. make use of trace option int Particle::iTopCopy()   && int Particle::iBotCopy()   provided by pythia
bool TraceShower(Event &event, int ip, double currentTime, vector<myParticle> &tc) {
	//pythia event, current particle location in event, current particle cumulated formation time, Shower chain
	//vector<myParticle> &tc only keep one copy if carbon copies of particle exist
	
	int botCopy = event[ip].iBotCopy(); // are used to trace carbon copies of the particle down to its bottom daughter. If there are no such carbon copies, the index of the particle itself will be returned. A carbon copy is when the "same" particle appears several times in the event record, but with changed momentum owing to recoil effects.

	Particle &p = event[botCopy];

	// as we already trace down to its last carbon copy, it will not have one and only one daughter
	if(p.daughterList().size()==0) {	// there are no daughters. The particle is the final one. End this trace
		if(verbose) cout<<"Final "<<botCopy<<" time="<<currentTime<<endl;

		if(botCopy!=ip) {	// is carbon copy. add into TimeChain. Otherwise, it should already in the record in the last function call as its mother's daughter.
			struct myParticle fp;
			fp.no = botCopy;
			fp.time = 0;
			fp.totalTime = fp.time+currentTime;

			tc.push_back(fp);
		}

		return true;			// END of this chain
	}
	else if(p.daughterList().size()==2) {	// 2 daughters

		if(event[p.daughter1()].motherList().size()==1) {	// 1->2, will caculate formation time for this process

			if(verbose) cout<<"1->2: "<<botCopy<<" - ";

			struct myParticle da1;
			da1.no = p.daughter1();

			struct myParticle da2;
			da2.no = p.daughter2();

			double deltaTime = CalFormTimeFrom4Da(event[p.daughter1()].p(),event[p.daughter2()].p());
			if(verbose) cout<<"Delta time = "<<deltaTime<<" -> ";

			da1.time = deltaTime;
			da2.time = deltaTime;

			da1.totalTime = da1.time + currentTime;
			tc.push_back(da1);
			if(verbose) cout<<" ("<<da1.no<<", "<<da1.time<<" ,"<<da1.totalTime<<") ";

			da2.totalTime = da2.time + currentTime;
			tc.push_back(da2);
			if(verbose) cout<<" ("<<da2.no<<", "<<da2.time<<" ,"<<da2.totalTime<<") "<<endl;


			TraceShower(event, p.daughter1(), da1.totalTime, tc);
			TraceShower(event, p.daughter2(), da2.totalTime, tc);
		}
		else if(event[p.daughter1()].motherList().size()==2) {	// 2->2. Set formation time as 0 for this process. 
		// Daughters will the longer time of two mothers. The program will record the time of daugther the first time encounter, stop futher trace until this daughter is encountered again from its another mother, then the program will decide what is the correct time for this daughter and its sister, and continue from there.
			if(verbose) {
				cout<<"2->2: ("<<event[p.daughter1()].motherList()[0]<<" + " <<event[p.daughter1()].motherList()[1]<<") -> ("<<p.daughter1()<<" + "<<p.daughter2()<<")"<< endl;
			}
			// Need to check whether daughters have been looped already
			bool Checked = false;
			vector<myParticle>::iterator idl = tc.begin();	// location of the daughter if in the TimeChain
			for(  ; idl!=tc.end(); idl++) {
				if(idl->no==p.daughter1()) {
					Checked = true;
					break;
				}
			}
			if(Checked) {	// if found, compare time record and proceed with the larger total time
				if(currentTime < idl->totalTime) {	
					currentTime = idl->totalTime;	// use the time in TimeChain

				}	// end if
				else {
					idl->totalTime = currentTime;	// update the record in TimeChain
				}
				struct myParticle da2;
				da2.no = p.daughter2();
				da2.time = 0;
				da2.totalTime = da2.time + currentTime;

				tc.push_back(da2);	// add daughter2, too

				if(verbose) cout<<"Was here before. currentTime = "<<currentTime<<endl;

				// Continue with two daugthers' tracing
				TraceShower(event, p.daughter1(), idl->totalTime, tc);
				TraceShower(event, p.daughter2(), da2.totalTime, tc);



			}	// end if 
			else {	// Checked = false. the first time to meet this daugther, push back the daughter info into TimeChain for later comparison of 2->2's two mother times
				struct myParticle da;
				da.no = p.daughter1();
				da.time = 0; 
				da.totalTime = da.time + currentTime;

				tc.push_back(da);

				if(verbose) cout<<"first time here. currentTime = "<<currentTime<<endl;

				return true;	// Hold the trace for now. Until we find the other mother particle for comparison
			}
		}	// end 2->2
		else {	// 2 daughters, more than 2 mothers
			if(warning) cout<<"WARNING!!! "<<botCopy<<" invovled in "<< event[p.daughter1()].motherList().size()<<" -> "<< p.daughterList().size()<<" process"<<endl;
		}

	}
	else {
		if(warning) cout<<"WARNNING!! Number of daughters: "<<p.daughterList().size()<<" Number of mothers: "<<event[p.daughter1()].motherList().size()<<endl;
	}
	return false;
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

	pythia.readString("HardQCD:all = on");
	char StringName[400];
	sprintf(StringName,"PhaseSpace:pTHatMin = %i",pTHat_Min);
	pythia.readString(StringName);
	sprintf(StringName,"PhaseSpace:pTHatMax = %i",pTHat_Max);
	pythia.readString(StringName);

	pythia.readString("PartonLevel:ISR = on");
	pythia.readString("PartonLevel:MPI = on");
	pythia.readString("PartonLevel:FSR = on");

	pythia.readString("HadronLevel:all = off");

	pythia.readString("WeakSingleBoson:all=off");
	pythia.readString("WeakDoubleBoson:all=off");

	pythia.readString("ParticleDecays:limitTau0=on");
	pythia.readString("ParticleDecays:tau0Max = 10");

	pythia.readString("PromptPhoton:all = on");


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
	//TFile* outFile = new TFile("hist.root", "RECREATE");
	ofstream fParton;
	if(JobId==-1) {
		sprintf(StringName,"../JETSCAPE/pythiaInput/parton_pTHat%03d_to%03d.dat",pTHat_Min,pTHat_Max);
	}
	else 	{
		sprintf(StringName,"../JETSCAPE/pythiaInput/parton_pTHat%03d_to%03d_Job%d.dat",pTHat_Min,pTHat_Max,JobId);
	}
	fParton.open(StringName,ios::out);
	cout<<"Output: "<<StringName<<endl;


	// Book histogram.
	TH1F *mult = new TH1F("mult","charged multiplicity", 100, -0.5, 799.5);
	TH1D *hTime = new TH1D("hTime","Formation Time of Final Particles (fm)",1000,-100,100);
	TH2D *hTimeVspT = new TH2D("hTimeVspT","Formation Time of Final Particles (fm) vs pT (GeV)",100,0,50,200,0,200);

	// Begin event loop. Generate event; skip if generation aborted.
	for (int iEvent = 0; iEvent < Nevents; ++iEvent) {
		if (!pythia.next()) continue;


		if(iEvent%1000==0) cout<<"event "<<iEvent<<endl;


		// =============== Section for Formation Time Calculation ===========

		// Time ordering chain 
		vector<myParticle> TimeChain;


		// Reset myParticle
		struct myParticle iTC;
		iTC.no = -1;
		iTC.time = 0;
		iTC.totalTime = 0;

		for (int ib = 1; ib <= 2; ++ib) {	// 1, 2 for incoming beam particle (proton)
			Particle& p = pythia.event[ib];	// particle ib

			if(verbose) {
				cout<<endl<<"===== Tracing START HERE for particle"<<ib<<" "<<p.name()<<" ====="<<endl;
				cout<<" Daughter parton: "<<endl;
				for (unsigned i=0; i<p.daughterList().size(); ++i)
					std::cout << ' ' << p.daughterList()[i];
				std::cout << '\n';
			}


			// Loop beam daughters' shower chains
			vector<int> beamDaughter = p.daughterList();
			for(vector<int>::iterator ibd = beamDaughter.begin(); ibd!=beamDaughter.end(); ibd++){
				if(verbose) cout<<endl<<"Trace beam parton "<<*ibd<<" "<<pythia.event[*ibd].name()<<" in beam particle"<<endl;
				iTC.no = *ibd;
				iTC.time = 0;
				iTC.totalTime = 0; 	// beam's daugther formation time sets to 0
				TimeChain.push_back(iTC);

				// Start to follow this shower
				TraceShower(pythia.event, iTC.no, iTC.totalTime, TimeChain);

			}	// end for loop beam daughter 

		}	// end for loop beam particle

		// =============== End of Section for Formation Time Calculation ===========




		// Find number of all final charged particles.
		int nCharged = 0;
		for (int i = 1; i < pythia.event.size(); ++i) {	// skip 0 for system itself

			if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
				++nCharged;	// count charged particle


			// Print mother daughter detailed list
			if(verbose) {
				Particle& p = pythia.event[i];	// particle i
				Vec4 mom = p.p();
				vector<int> motherlist = p.motherList();
				vector<int> daughterlist = p.daughterList();

				cout<<endl<<"Particle "<<i<<": ";
				cout<<" id = "<<p.id();
				cout<<" name = "<<p.name();
				cout<<" status = "<<p.status();
				//cout<<" tau = "<<p.tau();
				//cout<<" tau0 = "<<p.tau0();
				//cout<<" scale = "<<p.scale();
				//cout<<" mother1 = "<<p.mother1();
				//cout<<" mother2 = "<<p.mother1();
				//cout<<" daughter1 = "<<p.daughter1();
				//cout<<" daughter2 = "<<p.daughter1();
				cout<<endl;

				cout<<" motherlist: "<<endl;
				for (unsigned i=0; i<motherlist.size(); ++i)
					std::cout << ' ' << motherlist[i];
				std::cout << '\n';

				cout<<" daughterlist: "<<endl;
				for (unsigned i=0; i<daughterlist.size(); ++i)
					std::cout << ' ' << daughterlist[i];
				std::cout << '\n';
			}

		}	// end for loop particles


		// Fill charged multiplicity in histogram. End event loop.
		mult->Fill( nCharged );

		// Final list summary
		int Nfparton = 0;
		if(verbose)  cout<<"Shower Particles"<<endl;
		for(vector<myParticle>::iterator ipf = TimeChain.begin(); ipf!=TimeChain.end(); ipf++) {
			if(verbose) PrintMyParticle(ipf);
			if(pythia.event[ipf->no].isFinal() && Accept4LBT(pythia.event[ipf->no].id())) Nfparton++;		// Final particle && can be accepted by LBT program
		}

		// Write dat information
		fParton << iEvent+1 << "\t" << Nfparton << "\t" <<pythia.info.pTHat()<<endl;
		
		// Check whether all final particles are looped
		for (int k = 0; k < pythia.event.size(); ++k) {
			if(pythia.event[k].isFinal() && Accept4LBT(pythia.event[k].id())) {
				bool ifound = false;
				for(vector<myParticle>::iterator ipf = TimeChain.begin(); ipf!=TimeChain.end(); ipf++) {
					if(ipf->no==k){		// found a match
						ifound = true;

						// Fill histogram
						hTime->Fill(ipf->totalTime*0.1973);
						hTimeVspT->Fill(pythia.event[k].pT(), ipf->totalTime*0.1973);
						// Write dat file
						fParton<< k+1 <<"\t"<<pythia.event[k].id()<<"\t"<<pythia.event[k].px()<<"\t"<<pythia.event[k].py()<<"\t"<<pythia.event[k].pz()<<"\t"<<pythia.event[k].e()<<"\t"<<pythia.event[k].m()<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<"0"<<"\t"<<ipf->totalTime*0.1973<<endl;		// convert into fm (was 1/GeV)

						if(verbose) {
							cout<<"Final particle "<<k<<" ("<<pythia.event[k].id()<<") formation time = "<<ipf->totalTime<<endl;
						}
						break;
					}
				}
				if(!ifound && warning) cout<<"WARNING!! Final particle "<<k<<" not in TimeChain"<<endl;

			}	// end if final
		}	// end for loop particles

	}	// end for loop events

	// Statistics on event generation.
	//pythia.stat();
	
	// Show histogram. Possibility to close it.
	//mult->Draw();
	//hTime->Draw();
	//std::cout << "\nDouble click on the histogram window to quit.\n";
	//gPad->WaitPrimitive();

	// Save histogram on file and close file.
	//mult->Write();
	//hTime->Write();
	//hTimeVspT->Write();
	//delete outFile;
	
	fParton.close();

	// Done.
	cout<<"====== Done ======"<<endl;
	return 0;
}
