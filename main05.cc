// main05.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// It studies jet production at the LHC, using SlowJet and CellJet.
// Note: the two finders are intended to construct approximately the same
// jet properties, but provides output in slightly different format,
// and have here not been optimized to show maximum possible agreement.

#include <iostream>
#include "Pythia8/Pythia.h"
#include "Pythia8/Basics.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"

using namespace Pythia8;

int PID = 1;

double Emin = 50;
double Emax = 50; //500.0;	//#LY

double TimeCutOff = 100.;	//#LY	(fm)

int verbose = 0;
int warning = 1;

	int nEvent    = 100000;	//#LY

//Initializing Structure: myParticle (no, formation time, integrated formation time) 
struct myParticle {
	int no;
	double time;
	double totalTime;
};

void PrintMyParticle(vector<myParticle>::iterator mp) {
	cout<<"no="<<mp->no<<" time="<<mp->time<<" totalTime="<<mp->totalTime<<endl;
}
//---------------------- 'w/ f.t.' Calculation ----------------------
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
	if(time*0.1973>TimeCutOff) time = TimeCutOff/0.1973;	// cut off
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
	if(time*0.1973>TimeCutOff) time = TimeCutOff/0.1973;	// cut off
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


//---------------------- 'w/ simple f.t.' Calculation ------------------------
double CalFormTime_SimpleIF(Vec4 mother, Vec4 daughter, double &x, double &kt) {
	x = dot3(mother,daughter)/mother.pAbs2();
	kt = cross3(mother,daughter).pAbs()/mother.pAbs();
	double time = 0; 
	if(kt>0 || kt<0) {
		//time = 2.*daughter.e()/pow(kt,2);
		//time = 2.*mother.e()/pow(kt,2);
		time = 2.*mother.e()*x*(1.-x)/pow(kt,2);
	}
	if(time*0.1973>TimeCutOff) time = TimeCutOff/0.1973;	// cut off
	if(verbose) cout<<"mother e="<<mother.e()<<" daughter e="<<daughter.e()<<" x="<<x<<" kt="<<kt<<" time="<<time;
	return time;
}

double CalFormTime_SimpleIF(Vec4 mother, Vec4 daughter) {
	double x = 0, kt = 0;
	return CalFormTime_SimpleIF(mother,  daughter, x, kt);
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

int main() {

	// output file
	ofstream fout;
	stringstream filename;
	filename << "jet_" << PID << "_R04.dat";
	fout.open(filename.str().c_str(), ios::out);

	TFile* outFile = new TFile("hist05.root", "RECREATE");
	int Nbins = 100;
	double xmin = 0, xmax = 35;
	// w/ f.t.
	TH1D *hTime = new TH1D("hTime","Formation Time of Final Particles (fm)",Nbins,xmin,xmax);
	hTime->GetXaxis()->SetTitle("Formation Time (fm)");
	hTime->GetYaxis()->SetTitle("Number of Partons");
	TH2D *hTimeVspT = new TH2D("hTimeVspT","Formation Time of Final Particles (fm) vs pT (GeV)",100,0,Emax,200,0,100);
	TProfile *pTimeVspT = new TProfile("pTimeVspT","Formation Time of Final Particles (fm) vs pT (GeV)",100,0,Emax);
	// w/ simple f.t.
	TH1D *hTime_SimpleIF = new TH1D("hTime_SimpleIF","Simple Formation Time of Final Particles (fm)",Nbins,xmin,xmax);
	hTime_SimpleIF->GetXaxis()->SetTitle("Formation Time (fm)");
	hTime_SimpleIF->GetYaxis()->SetTitle("Number of Partons");
	TH2D *hTimeVspT_SimpleIF = new TH2D("hTimeVspT_SimpleIF","Simple Formation Time of Final Particles (fm) vs pT (GeV)",100,0,Emax,200,0,100);
	TProfile *pTimeVspT_SimpleIF = new TProfile("pTimeVspT_SimpleIF","Simple Formation Time of Final Particles (fm) vs pT (GeV)",100,0,Emax);
	TProfile *pkTVspT = new TProfile("pkTVspT","k_{T} vs p_{T}",100,0,Emax);
	TProfile *pEx1xVspT = new TProfile("pEx1xVspT","E*x(1-x) vs p_{T}",100,0,Emax);
	// default LBT time
	TH1D *hTime_default = new TH1D("hTime_default","LBT default Formation Time of Final Particles (fm)",Nbins,xmin,xmax);
	TH2D *hTimeVspT_default = new TH2D("hTimeVspT_default","LBT default Formation Time of Final Particles (fm) vs pT (GeV)",100,0,Emax,200,0,100);
	hTime_default->GetXaxis()->SetTitle("Formation Time (fm)");
	hTime_default->GetYaxis()->SetTitle("Number of Partons");
	TProfile *pTimeVspT_default = new TProfile("pTimeVspT_default","LBT default Formation Time of Final Particles (fm) vs pT (GeV)",100,0,Emax);


	// Number of events, generated and listed ones.
	//#LY int nEvent    = 1;
	//int nListJets = 5;

	// Generator. LHC process and output selection. Initialization.
	Pythia pythia;

	pythia.readString("Beams:eCM = 200.");	// #Ly

	pythia.readString("ProcessLevel:all=off");
	pythia.readString("PartonLevel:FSR=on");
	pythia.readString("HadronLevel:Hadronize=off");

	pythia.readString("Random:setSeed = on");
	pythia.readString("Random:seed = 0");

	pythia.init();

	Rndm pythiaRandom;
	pythiaRandom.init(0);

	// Common parameters for the two jet finders.
	double etaMax   = 4.;
	double radius   = 0.4;
	double pTjetMin = 10.;
	// Exclude neutrinos (and other invisible) from study.
	int    nSel     = 2;

	// Set up SlowJet jet finder, with anti-kT clustering
	// and pion mass assumed for non-photons..
	//SlowJet slowJet( -1, radius, pTjetMin, etaMax, nSel, 1);

	// Begin event loop. Generate event. Skip if error.
	for ( int iEvent = 0; iEvent < nEvent; ++iEvent ) {

		pythia.event.reset();//reset the event

		double Einit = Emax; 	//#LY 
		if(Emax!=Emin ) Einit = Emin + pythiaRandom.flat() * (Emax - Emin);	// #LY
		//#LY double Einit = Emin + pythiaRandom.flat() * (Emax - Emin);	

		if(verbose) cout << "Einit: " << Einit << endl;

		if(PID==21) { //append 2 back to back gluons
			pythia.event.append(21,23,101,102,Einit,0.0,0.0,Einit);
			pythia.event.append(21,23,102,101,-Einit,0.0,0.0,Einit);
		} else if(PID<=3 && PID>0) { // append 2 light quarks
			int lightID;
			double rr = pythiaRandom.flat();
			if ( rr < 1.0/3.0 ) lightID = 1;
			else if ( rr < 2.0/3.0) lightID = 2;
			else lightID = 3;
			pythia.event.append(lightID,23,101,0,Einit,0.0,0.0,Einit);
			pythia.event.append(-lightID,23,0,101,-Einit,0.0,0.0,Einit);
		} else if(PID==4 || PID==5) { // append 2 light quarks
			pythia.event.append(PID,23,101,0,Einit,0.0,0.0,Einit);
			pythia.event.append(-PID,23,0,101,-Einit,0.0,0.0,Einit);
		} else {
			cout << "Only accept 1, 2, 3, 4, 5, 21 for PID." << endl;
		}

		pythia.event[1].scale(Einit);
		pythia.event[2].scale(Einit);

		pythia.forceTimeShower(1, 2, Einit);

		if (!pythia.next()) continue;      

		if(verbose) pythia.event.list();

		vector<myParticle> TimeChain;

		for (int ib = 1; ib <= 2; ++ib) {	// 1, 2 for incoming beam particle (proton)
			Particle& p = pythia.event[ib];	// particle ib


			struct myParticle iTC;
			iTC.no = ib;
			iTC.time = 0;
			iTC.totalTime = 0;

			TimeChain.push_back(iTC);

			TraceShower(pythia.event, iTC.no, iTC.totalTime, TimeChain);
		}

		// Check whether all final particles are looped
		for (int k = 0; k < pythia.event.size(); ++k) {
			if(pythia.event[k].isFinal() && Accept4LBT(pythia.event[k].id())) {
				// --------- w/ f.t.
				bool ifound = false;
				for(vector<myParticle>::iterator ipf = TimeChain.begin(); ipf!=TimeChain.end(); ipf++) {
					if(ipf->no==k){		// found a match
						ifound = true;

						// Fill histogram
						hTime->Fill(ipf->totalTime*0.1973);
						hTimeVspT->Fill(pythia.event[k].pT(), ipf->totalTime*0.1973);
						pTimeVspT->Fill(pythia.event[k].pT(), ipf->totalTime*0.1973);
			
						if(verbose) {
							cout<<"Final particle "<<k<<" ("<<pythia.event[k].id()<<") formation time = "<<ipf->totalTime<<endl;
						}
						break;
					}
				}
				if(!ifound && warning) cout<<"WARNING!! Final particle "<<k<<" not in TimeChain"<<endl;

				// End of w/ f.t.

				// -------  w/ simple f.t.
				double timeSIF = 0;		// if not from hardest parton, use time = 0
				int ancestor = 0;
				double tmpx = 0, tmpkt = 0;
				if ( pythia.event[k].isAncestor(1) ) {		// Back to embedded parton 1 or 2
					ancestor = 1;
					timeSIF  = 0.1973*CalFormTime_SimpleIF(pythia.event[ancestor].p(),pythia.event[k].p(), tmpx, tmpkt);
				}	// End of If 	
				else if ( pythia.event[k].isAncestor(2)  ) {
					ancestor = 2;
					timeSIF  = 0.1973*CalFormTime_SimpleIF(pythia.event[ancestor].p(),pythia.event[k].p(), tmpx, tmpkt);

				}	// End of If	

				if(ancestor) {
					hTime_SimpleIF->Fill(timeSIF);
					pTimeVspT_SimpleIF->Fill(pythia.event[k].pT(),timeSIF);
					hTimeVspT_SimpleIF->Fill(pythia.event[k].pT(),timeSIF);
					pkTVspT->Fill(pythia.event[k].pT(),tmpkt);
					pEx1xVspT->Fill(pythia.event[k].pT(),Einit*tmpx*(1.-tmpx));
				}
				// End w/ simple IF


				// ------ default LBT time
				double timedef = 0;
				if(pythia.event[k].pT()>0) timedef= 0.1973*2*pythia.event[k].e()/pythia.event[k].pT();
				//cout<<"pt = "<<pythia.event[k].pT()<<" time = "<<timedef<<endl;
				hTime_default->Fill(timedef);
				pTimeVspT_default->Fill(pythia.event[k].pT(),timedef);
				hTimeVspT_default->Fill(pythia.event[k].pT(),timedef);
			}	// end if final
		}	// end for check loop particles


		/*
		// Analyze Slowet jet properties. List first few.
		slowJet. analyze( pythia.event );

		//cout << "size: " << slowJet.sizeJet() << endl;

		if ( slowJet.sizeJet() == 0 ) continue;
		fout << "# Event " << iEvent << "  " << slowJet.sizeJet() << endl;

		//for (int iJet = 0; iJet < min(2, slowJet.sizeJet()); iJet++) {
		for (int iJet = 0; iJet < slowJet.sizeJet(); iJet++) {

		double tot_px = 0.0;
		double tot_py = 0.0;
		double tot_pz = 0.0;
		double tot_p0 = 0.0;
		double tot_m  = 0.0;

		fout << "# Jet " << iJet << "  " << slowJet.multiplicity(iJet) << "  " << slowJet.p(iJet).px() << "  " << slowJet.p(iJet).py() << "  " << slowJet.p(iJet).pz() << "  " << slowJet.p(iJet).e() << "  " << slowJet.m(iJet) << endl;

		vector<int> parIndices = slowJet.constituents(iJet);
//        cout << parIndices.size() << endl;

for (int iPar = 0; iPar < slowJet.multiplicity(iJet); iPar++) {
int parIndex = parIndices[iPar];
fout << iPar << "  " << pythia.event[parIndex].id() << "  " << pythia.event[parIndex].px() << "  " << pythia.event[parIndex].py() << "  " << pythia.event[parIndex].pz() << "  " << pythia.event[parIndex].e() << "  " << pythia.event[parIndex].m() << endl;
tot_px += pythia.event[parIndex].px();
tot_py += pythia.event[parIndex].py();
tot_pz += pythia.event[parIndex].pz();
tot_p0 += pythia.event[parIndex].e();
}

tot_m   = sqrt(tot_p0*tot_p0 - tot_px*tot_px - tot_py*tot_py - tot_pz*tot_pz);
//cout << "check 4-p:  " << tot_px << "  " << tot_py << "  " << tot_pz << "  " << tot_p0 << "  " << tot_m << endl;

}
*/

//    cout << slowJet.constituents(1).size() << endl;
//    slowJet.list();

	} // end event loop

	if(1) { 	// Normalize
		hTime->Scale(1./nEvent);
		hTime_SimpleIF->Scale(1./nEvent);
		hTime_default->Scale(1./nEvent);
	}


	TCanvas *c = new TCanvas();
        c->SetFrameLineWidth(3);
        c->SetLeftMargin(0.15);
        c->SetBottomMargin(0.2);
	c->SetLogy();

	hTime->SetLineWidth(3);
	hTime_SimpleIF->SetLineWidth(3);
	hTime_default->SetLineWidth(3);
	hTime_SimpleIF->SetLineColor(2);
	hTime_default->SetLineColor(3);
	hTime_SimpleIF->SetMarkerColor(2);
	hTime_default->SetMarkerColor(3);
        hTime_default->GetXaxis()->SetNdivisions(505);
        hTime_default->GetYaxis()->SetNdivisions(505);
        hTime_default->GetXaxis()->SetLabelSize(0.06);
        hTime_default->GetXaxis()->SetTitleSize(0.06);
        hTime_default->GetYaxis()->SetLabelSize(0.06);
        hTime_default->GetYaxis()->SetTitleSize(0.06);
        hTime->GetXaxis()->SetNdivisions(505);
        hTime->GetYaxis()->SetNdivisions(505);
        hTime->GetXaxis()->SetLabelSize(0.06);
        hTime->GetXaxis()->SetTitleSize(0.06);
        hTime->GetYaxis()->SetLabelSize(0.06);
        hTime->GetYaxis()->SetTitleSize(0.06);
        hTime->GetYaxis()->SetTitleOffset(1.2);
        hTime->GetXaxis()->SetRangeUser(0,25);
        hTime->GetYaxis()->SetRangeUser(2e-4,9);
	pTimeVspT->SetLineColor(hTime->GetLineColor());
	pTimeVspT_SimpleIF->SetLineColor(hTime_SimpleIF->GetLineColor());
	pTimeVspT_default->SetLineColor(hTime_default->GetLineColor());
	hTimeVspT->SetLineColor(hTime->GetLineColor());
	hTimeVspT_SimpleIF->SetLineColor(hTime_SimpleIF->GetLineColor());
	hTimeVspT_default->SetLineColor(hTime_default->GetLineColor());

	//hTime_default->Draw();
	hTime->Draw("");
	hTime_SimpleIF->Draw("same");

	TLegend *leg = new TLegend(0.7,0.7,0.88,0.88);
	if(Emax>Emin) leg->SetHeader(Form("%g<E_{init}<%gGeV/#it{c}",Emin, Emax));
	else leg->SetHeader(Form("E_{init} = %gGeV/#it{c}",Emax));
	leg->AddEntry(hTime_SimpleIF,"Setup 2","l");//w/ simple f.t.","l");
	leg->AddEntry(hTime,"Setup 3","l");//w/ f.t.","l");
	//leg->AddEntry(hTime_default,"default","l");
	leg->Draw("same");

	TCanvas *c2 = new TCanvas();
        c2->SetFrameLineWidth(3);
        c2->SetLeftMargin(0.15);
        c2->SetBottomMargin(0.2);
	c2->SetLogy();
	pTimeVspT->GetYaxis()->SetRangeUser(0.3,500);
        pTimeVspT->GetXaxis()->SetTitle("Final parton p_{T} (GeV/#it{c})");
        pTimeVspT->GetYaxis()->SetTitle("Formation Time (fm)");
        pTimeVspT->GetXaxis()->SetNdivisions(505);
        pTimeVspT->GetYaxis()->SetNdivisions(505);
        pTimeVspT->GetXaxis()->SetLabelSize(0.06);
        pTimeVspT->GetXaxis()->SetTitleSize(0.06);
        pTimeVspT->GetYaxis()->SetLabelSize(0.06);
        pTimeVspT->GetYaxis()->SetTitleSize(0.06);
	pTimeVspT->Draw();
	pTimeVspT_default->Draw("same");
	pTimeVspT_SimpleIF->Draw("same");
	leg->Draw("same");

	c->Write();
	c2->Write();

	hTime->Write();
	hTimeVspT->Write();
	pTimeVspT->Write();

	hTime_SimpleIF->Write();
	pTimeVspT_SimpleIF->Write();
	hTimeVspT_SimpleIF->Write();
	pkTVspT->Write();
	pEx1xVspT->Write();

	hTime_default->Write();
	pTimeVspT_default->Write();
	hTimeVspT_default->Write();

	outFile->Close();
	//delete outFile;
	
	// Done.
	return 0;
}

