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

// ROOT, for interactive graphics.
#include "TVirtualPad.h"
#include "TApplication.h"

// ROOT, for saving file.
#include "TFile.h"

using namespace Pythia8;

int verbose = 0;		// Print More information (1), less information (0)
int warning = 1;		// Print warning msg (1), no warning msg (0)

//Initializing Structure: myParticle (no, formation time, integrated formation time) 
struct myParticle {
	int no;
	double time;
	double totalTime;
};

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

// Tracing shower develop using mother-daughter relationship
bool TraceShower(Event &event, int ip, double currentTime, int i_currentMaxE, vector<myParticle> &tc, vector<double> &FinalList) {
	//pythia event, current particle location in event, current particle cumulated formation time, location in event for max energy in all carbon copies of current particle, Shower chain, list final particle location in event
	
	Particle &p = event[ip];

	if(p.daughterList().size()==0) {	// there are no daughters
		if(verbose) cout<<"Final "<<ip<<" time="<<currentTime<<endl;
		FinalList.push_back(currentTime);
		return true;			// END of this chain
	}
	else if(p.daughterList().size()==1) {	// the particle has a "carbon copy" as its sole daughter, but with changed momentum as a "recoil" effect, e.g. in a shower;
		if(event[p.daughter1()].motherList().size()>1) {	// more than 1 mom
			if(warning) cout<<"WARNING!!! "<<ip<<" invovled in "<< event[p.daughter1()].motherList().size()<<" -> 1 process"<<endl;
		}

		//Use formation time = 0 for 1->1
		struct myParticle da;
		da.no = p.daughter1();
		da.time = 0;
		da.totalTime = da.time + currentTime;

		tc.push_back(da);
		if(verbose) cout<<"1->1:"<<ip<<"-> ("<<da.no<<", "<<da.time<<" ,"<<da.totalTime<<")"<<endl;
		if(event[da.no].e()>event[i_currentMaxE].e()) {	// found a carbon copy with higher energy, set it as the new current max energy copy of this particle
			i_currentMaxE = da.no;
		}

		if(verbose) cout<<"Copy/Mother Energy = "<<event[da.no].e()<<"/"<<p.e()<<"="<<1.*event[da.no].e()/p.e()<<" i_currentMaxE="<<i_currentMaxE<<endl;

		TraceShower(event, p.daughter1(), da.totalTime, i_currentMaxE, tc, FinalList);
	}
	else if(p.daughterList().size()==2) {	// 2 daughters

		if(event[p.daughter1()].motherList().size()==1) {	// 1->2

			if(verbose) cout<<"1->2: "<<ip<<"(currentMaxE="<<i_currentMaxE<<") - ";

			struct myParticle da1;
			da1.no = p.daughter1();
			//da1.time = CalFormTime(event[i_currentMaxE].p(), event[p.daughter1()].p());
			//da1.time = CalFormTime(p.p(), event[p.daughter1()].p());
			//da1.time = CalFormTime(event[p.iTopCopy()].p(), event[p.daughter1()].p());
			//da1.totalTime = da1.time + currentTime;
			//tc.push_back(da1);
			//if(verbose) cout<<" ("<<da1.no<<", "<<da1.time<<" ,"<<da1.totalTime<<") ";

			struct myParticle da2;
			da2.no = p.daughter2();
			//da2.time = CalFormTime(event[i_currentMaxE].p(), event[p.daughter2()].p());
			//da2.time = CalFormTime(p.p(), event[p.daughter2()].p());
			//da2.time = CalFormTime(event[p.iTopCopy()].p(), event[p.daughter2()].p());

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

			//cout<<"Sum Vec4 Direct"<<p.p()-(event[p.daughter1()].p()+event[p.daughter2()].p())<<endl;
			//cout<<"Top: "<<p.iTopCopy()<<endl;
			//cout<<"Sum Vec4 Top"<<event[p.iTopCopy()].p()-(event[p.daughter1()].p()+event[p.daughter2()].p())<<endl;

			
			//if(verbose) {
			//	if(p.e()<event[da1.no].e()) {
			//		cout<<"WARNING!!! MOTHER energy less than daughter. Mother: "<<ip<<" i_currentMaxE "<<i_currentMaxE<<" "<<event[i_currentMaxE].name()<<" "<<event[i_currentMaxE].p()<<" Daughter1: "<<da1.no<<" "<<event[da1.no].name()<<" "<<event[da1.no].p()<<" Daughter2: "<<da2.no<<" "<<event[da2.no].name()<<" "<<event[da2.no].p()<<endl;
			//		//event.list();
			//	}
			//}




			// i_currentMaxE shall be reset for each daughter to find their max energy copies
			TraceShower(event, p.daughter1(), da1.totalTime, p.daughter1(), tc, FinalList);
			TraceShower(event, p.daughter2(), da2.totalTime, p.daughter2(), tc, FinalList);
		}
		else if(event[p.daughter1()].motherList().size()==2) {	// 2->2
			if(verbose) {
				cout<<"2->2: ("<<event[p.daughter1()].motherList()[0]<<" + " <<event[p.daughter1()].motherList()[1]<<") -> ("<<p.daughter1()<<" + "<<p.daughter2()<<")"<< endl;
			}
			// Need to check whether daughters have been looped already
			bool Checked = false;
			vector<myParticle>::iterator idl = tc.begin();
			for(  ; idl!=tc.end(); idl++) {
				if(idl->no==p.daughter1()) {
					Checked = true;
					break;
				}
			}
			if(Checked) {	// if found, compare time record and proceed with the larger total time
				if(currentTime < idl->totalTime) {
					currentTime = idl->totalTime;

				}	// end if
				else {
					idl->totalTime = currentTime;
				}
				struct myParticle da2;
				da2.no = p.daughter2();
				da2.time = 0;
				da2.totalTime = da2.time + currentTime;

				tc.push_back(da2);	// add daughter2, too

				if(verbose) cout<<"Was here before. currentTime = "<<currentTime<<endl;

				// Continue with two daugthers' tracing
				TraceShower(event, p.daughter1(), idl->totalTime, p.daughter1(), tc, FinalList);
				TraceShower(event, p.daughter2(), da2.totalTime, p.daughter2(), tc, FinalList);



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
		else {
			if(warning) cout<<"WARNING!!! "<<ip<<" invovled in "<< event[p.daughter1()].motherList().size()<<" -> "<< p.daughterList().size()<<" process"<<endl;
		}

	}
	else {
		if(warning) cout<<"WARNNING!! Number of daughters: "<<p.daughterList().size()<<endl;
	}
	return false;
}

int main(int argc, char* argv[]) {

	// Create the ROOT application environment.
	TApplication theApp("hist", &argc, argv);

	// Create Pythia instance and set it up to generate hard QCD processes
	// pTHat = 20 GeV - 25 GeV for pp collisions at 200 GeV.
	// Stop the generation after the hard process and parton-level activity has been generated, but before the hadron-level steps.
	Pythia pythia;
	pythia.readString("Beams:eCM = 200.");
	pythia.readString("HardQCD:all = on");
	pythia.readString("PhaseSpace:pTHatMin = 20.");
	pythia.readString("PhaseSpace:pTHatMax = 25.");
	pythia.readString("HadronLevel:all = off");

	// Less screen print out
	if(!verbose) pythia.readString("Print:quiet = on");

	// Set Random seed for replicable result
	//pythia.readString("Random:setSeed = on");
	//pythia.readString("Random:seed = 2021");
	

	pythia.init();

	// Create file on which histogram(s) can be saved.
	TFile* outFile = new TFile("hist.root", "RECREATE");

	// Book histogram.
	TH1F *mult = new TH1F("mult","charged multiplicity", 100, -0.5, 799.5);
	TH1D *hTime = new TH1D("hTime","Formation Time of Final Particles",1000,-100,100);

	// Begin event loop. Generate event; skip if generation aborted.
	for (int iEvent = 0; iEvent < 1000; ++iEvent) {
		if (!pythia.next()) continue;

		if(iEvent%1000==0) cout<<"event "<<iEvent<<endl;

		// Time ordering chain 
		vector<myParticle> TimeChain;

		// Shower Final Particle cumulated Formation Time
		vector<double> FinalList;



		// Reset myParticle
		struct myParticle iTC;
		iTC.no = -1;
		iTC.time = 0;
		iTC.totalTime = 0;

		// Find number of all final charged particles.
		int nCharged = 0;
		for (int i = 1; i < pythia.event.size(); ++i) {	// skip 0 for system itself
			if (pythia.event[i].isFinal() && pythia.event[i].isCharged())
				++nCharged;
			Particle& p = pythia.event[i];	// particle i
			Vec4 mom = p.p();
			vector<int> motherlist = p.motherList();
			vector<int> daughterlist = p.daughterList();

			if(verbose) {
				cout<<i<<" --> ";
				cout<<" id = "<<p.id();
				cout<<" name = "<<p.name();
				cout<<" status = "<<p.status();
				cout<<" tau = "<<p.tau();
				cout<<" tau0 = "<<p.tau0();
				cout<<" scale = "<<p.scale();
				cout<<" mother1 = "<<p.mother1();
				cout<<" mother2 = "<<p.mother1();
				cout<<" daughter1 = "<<p.daughter1();
				cout<<" daughter2 = "<<p.daughter1();
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

			//if(p.status()==-23) {		// outgoing particles of the hardest subprocess
			if(p.status()==-12) {		// incoming beam 
			if(verbose) cout<<endl<<"===== Tracing START HERE for particle"<<i<<" "<<p.name()<<" ====="<<endl;
				iTC.no = i;
				iTC.time = 0;
				iTC.totalTime = 0;

				TimeChain.push_back(iTC);	

				// Loop beam daughters' shower chains
				vector<int> beamDaughter = p.daughterList();
				for(vector<int>::iterator ibd = beamDaughter.begin(); ibd!=beamDaughter.end(); ibd++){
					if(verbose) cout<<endl<<"Trace parton "<<*ibd<<" "<<pythia.event[*ibd].name()<<" in beam particle"<<endl;
					iTC.no = *ibd;
					iTC.time = 0;
					iTC.totalTime = 0; 	// beam's daugther formation time sets to 0
					TimeChain.push_back(iTC);

					// Start to follow this shower
					TraceShower(pythia.event, iTC.no, iTC.totalTime, iTC.no, TimeChain, FinalList);
					// this current max energy copy is itself

				}	// end for loop beam daughter 

			}	// end for loop beam

		}	// end for loop pythia


		// Fill charged multiplicity in histogram. End event loop.
		mult->Fill( nCharged );

		// Final list summary
		if(verbose) cout<<"Final Shower Particles"<<endl;
		for(vector<double>::iterator it = FinalList.begin(); it!=FinalList.end(); it++) {
			if(verbose) cout<<*it<<endl;
			hTime->Fill(*it);
		}
		
		// Check whether all final particles are looped
		for (int k = 0; k < pythia.event.size(); ++k) {
			if(pythia.event[k].isFinal()) {
				bool ifound = false;
				for(vector<myParticle>::iterator ipf = TimeChain.begin(); ipf!=TimeChain.end(); ipf++) {
					if(ipf->no==k){		// found a match
						ifound = true;
						if(verbose) {
							cout<<"Final particle "<<k<<" formation time = "<<ipf->totalTime<<endl;
						}
						break;
					}
				}
				if(!ifound && warning) cout<<"WARNING!! Final particle "<<k<<" not in TimeChain"<<endl;

			}
		}

	}	// end for loop events

	// Statistics on event generation.
	//pythia.stat();
	
	// Show histogram. Possibility to close it.
	//mult->Draw();
	//hTime->Draw();
	//std::cout << "\nDouble click on the histogram window to quit.\n";
	//gPad->WaitPrimitive();

	// Save histogram on file and close file.
	mult->Write();
	hTime->Write();
	delete outFile;
	
	// Done.
	cout<<"====== Done ======"<<endl;
	return 0;
}
