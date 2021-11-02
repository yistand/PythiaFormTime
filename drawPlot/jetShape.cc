#include <iostream>    //....C++ headers
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

#include <cstdio>    //....C headers
#include <cstdlib>
#include <cmath>
#include <ctime>

//#include "fastjet/ClusterSequence.hh"    //....FASTJET headers
#include "fastjet/ClusterSequenceArea.hh"  // use this instead of the "usual" ClusterSequence to get area support
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "TH1.h"	//....ROOT headers
#include "TFile.h"
#include "TMath.h"

using namespace std;
using namespace fastjet;


//....Classes definition
//....FASTJET class dealing with negative particles
typedef fastjet::JetDefinition::Recombiner Recombiner;
/// Recombiner class that propagates the user index and arranges the
/// recombination accordingly
class NegativeEnergyRecombiner : public  Recombiner {
	public:
		NegativeEnergyRecombiner(const int ui) : _ui(ui) {}

		virtual std::string description() const {return "E-scheme Recombiner that checks a flag for a 'negative momentum' particle, and subtracts the 4-momentum in recombinations.";}

		/// recombine pa and pb and put result into pab
		virtual void recombine(const fastjet::PseudoJet & pa, 
				const fastjet::PseudoJet & pb, 
				fastjet::PseudoJet & pab) const {

			int ai=1,bi=1;

			// If a particle is flagged, restore its real negative energy. 
			// The -1 will flip the full 4-momentum, reversing the convention for 
			// negative energy particles.
			if (pa.user_index() < 0) { ai = -1;}
			if (pb.user_index() < 0) { bi = -1;}

			// recombine particles
			pab = ai*pa+bi*pb;

			// if the combination has negative energy, flip the whole 4-momentum and flag it, 
			// so that we have the same convention as for input particles
			if(pab.E() < 0) { 
				pab.set_user_index(_ui); 
				pab.reset_momentum(-pab.px(),
						-pab.py(),
						-pab.pz(),
						-pab.E());
			} else { pab.set_user_index(0);}

		}

	private:
		const int _ui;  
};


double GetPtSum(const PseudoJet &jet, const vector<fastjet::PseudoJet> &input_particles, TH1D *hdr, double ptmin = 0.2, double ptmax = 30){	// all particles in the event
	double PtSum = 0;
	int nCon = input_particles.size();
	for(int i_ncon=0; i_ncon<nCon; i_ncon++){
		double i_pt = input_particles[i_ncon].perp();	// scalar transverse momentum 
		if( i_pt<ptmin || i_pt>=ptmax) continue;
		if(input_particles[i_ncon].user_index()<0) i_pt=-i_pt;	// if hole (negative particle), subtract its pt
		PtSum+=i_pt;
		double dr = input_particles[i_ncon].delta_R(jet);	// return the cylinder (rap-phi) distance between this jet and another
		hdr->Fill(dr,i_pt);
	}
	return PtSum;
}

double GetPtSum(const PseudoJet &jet, TH1D *hdr, double ptmin = 0.2, double ptmax = 30){	// jet particles only
	double PtSum = 0;
	vector<PseudoJet> jconstituents = jet.constituents();
	int nCon = jconstituents.size();
	for(int i_ncon=0; i_ncon<nCon; i_ncon++){
		double i_pt = jconstituents[i_ncon].perp();	// scalar transverse momentum 
		if( i_pt<ptmin || i_pt>=ptmax) continue;
		if(jconstituents[i_ncon].user_index()<0) i_pt=-i_pt;	// if hole (negative particle), subtract its pt
		PtSum+=i_pt;
		double dr = jconstituents[i_ncon].delta_R(jet);	// return the cylinder (rap-phi) distance between this jet and another
		hdr->Fill(dr,i_pt);
	}
	return PtSum;
}


double GetPtInConeSum(const PseudoJet &jet, double R, const vector<fastjet::PseudoJet> &input_particles, double ptmin = 0.2, double ptmax = 30){	// all particles inside jet cone 
	double PtSum = 0;
	int nCon = input_particles.size();
	for(int i_ncon=0; i_ncon<nCon; i_ncon++){
		double i_pt = input_particles[i_ncon].perp();	// scalar transverse momentum 
		if( i_pt<ptmin || i_pt>=ptmax) continue;
		if(input_particles[i_ncon].user_index()<0) i_pt=-i_pt;	// if hole (negative particle), subtract its pt
		double dr = input_particles[i_ncon].delta_R(jet);	// return the cylinder (rap-phi) distance between this jet and another
		if(dr<=R) PtSum+=i_pt;
	}
	return PtSum;
}
bool Fill4dNtrk(const PseudoJet &jet, const vector<fastjet::PseudoJet> &input_particles, TH1D *hdr, double ptmin = 0.2, double ptmax = 30){
	int nCon = input_particles.size();
	for(int i_ncon=0; i_ncon<nCon; i_ncon++){
		double i_pt = input_particles[i_ncon].perp();	// scalar transverse momentum 
		if( i_pt<ptmin || i_pt>=ptmax) continue;
		double dr = input_particles[i_ncon].delta_R(jet);	// return the cylinder (rap-phi) distance between this jet and another
		if(input_particles[i_ncon].user_index()<0) {
			hdr->Fill(dr,-1);	// subtract hole (negative particle)
		}
		else {
			hdr->Fill(dr,1);
		}
	}
	return true;
}

bool Fill4dNtrk(const PseudoJet &jet, TH1D *hdr, double ptmin = 0.2, double ptmax = 30){
	vector<PseudoJet> jconstituents = jet.constituents();
	int nCon = jconstituents.size();
	for(int i_ncon=0; i_ncon<nCon; i_ncon++){
		double i_pt = jconstituents[i_ncon].perp();	// scalar transverse momentum 
		if( i_pt<ptmin || i_pt>=ptmax) continue;
		double dr = jconstituents[i_ncon].delta_R(jet);	// return the cylinder (rap-phi) distance between this jet and another
		if(jconstituents[i_ncon].user_index()<0) {
			hdr->Fill(dr,-1);	// subtract hole (negative particle)
		}
		else {
			hdr->Fill(dr,1);
		}
	}
	return true;
}

int main(int argc, char* argv[]) {		//  outputFile .root file
	//double pt_cut=0.01;
	double pt_cut_min=0.2;
	//#LY 2021.09.15 double pt_cut_max=30.;
	//    double R=0.4;
	double R=0.4;
	double AreaCut = 0.4; 	
	// PRC 102, 054913 (2020)
	if(R==0.2) AreaCut=0.07;
	if(R==0.3) AreaCut=0.2;
	if(R==0.4) AreaCut=0.4;
	double pTsubCut = 20;//5.0;
	cout<<"R = "<<R<<" AreaCut = "<<AreaCut<<" jetpt>"<<pTsubCut<<" pt_cut_min="<<pt_cut_min<<endl;

	//double y[5];
	//y[0]=0.0, y[1]=0.3, y[2]=0.8, y[3]=1.2, y[4]=2.1;

	//....time counting begins
	struct tm *local_start;
	time_t time_start;
	time_start=time(NULL);
	local_start=localtime(&time_start);

	char buf1[80];
	strftime(buf1,80,"Current Time: %Y-%m-%d %H:%M:%S",local_start);
	cout << "the program starts at:" <<endl;
	cout << buf1 << endl;

	ifstream parton1,parton2;
	string inputFile_p,inputFile_n,outputFile;
	string line;

	//....Generator. Initialization. 

	if (argc == 3) {
		inputFile_p = argv[1];
		outputFile = argv[2];
	} else if (argc == 4) {
		inputFile_p = argv[1];
		inputFile_n = argv[2];
		outputFile = argv[3];
	} else {
		cout << "./exe input_file_p (input_file_n) output_file" << endl;
		exit(EXIT_FAILURE);
	}

	parton1.open(inputFile_p.c_str());
	if (!parton1) {
		cout << "parton1.dat not open!!!" << endl;
		exit(1);
	}

	if (argc == 4) {
		parton2.open(inputFile_n.c_str());
		if (!parton2) {
			cout << "parton2.dat not open!!!" << endl;
			exit(1);
		}
	}


	unsigned nEvent=0, iEvent1=0, iEvent2=0, np1=0, np2=0;
	unsigned njets=0;// , nrap[6] = {0};
	//bool eventbool = false;
	//unsigned ntriggerjet = 0;
	int dummyInt;
	double dummyF;
	double initPT; // ,initRange;

	double PtSum=0;		// input particles in event. sum taken per jet
	double PtInConeSum=0;	// input particles in jet cone
	double PtJetConSum=0;	// jet consituents
	double PtJetSum=0;	// sum of jet pt
	double JetsSum=0;
	TH1D *hpt = new TH1D("hpt","hpt;jet_pT_sub;nJets",500,0,100);
	// hPtSumVsdr & hNtrkVsdr need to be scale by parameters in hSums (see later) and bin/pt widths for rho & Yield. Use these two after merging files with hadd
	TH1D *hPtSumVsdr = new TH1D("hPtSumVsdr","hPtSumVsdr;#Deltar;#Sigma_{jets}#Sigma_{trk#in#delta#Deltar}p_{T}^{trk}",100,0,TMath::Pi());
	TH1D *hNtrkVsdr = new TH1D("hNtrkVsdr","hNtrkVsdr;#Deltar;#Sigma_{jets}N_{trk#in#delta#Deltar}",100,0,TMath::Pi());
	TH1D *hPtJetSumVsdr = new TH1D("hPtJetSumVsdr","hPtJetSumVsdr;#Deltar;#Sigma_{jets}#Sigma_{trk#in#delta#Deltar}p_{T}^{jet trk}",100,0,TMath::Pi());
	TH1D *hNtrkJetVsdr = new TH1D("hNtrkJetVsdr","hNtrkJetVsdr;#Deltar;#Sigma_{jets}N_{jet trk#in#delta#Deltar}",100,0,TMath::Pi());

	//    while (!parton1.eof()&&!parton2.eof()) {
	while (!parton1.eof()) {
		//        double tot_ener=0.0;
		//        double tot_px=0.0;
		//        double tot_py=0.0;
		//        double tot_pz=0.0; 
		//        parton1 >> iEvent1 >> np1 >> initPT >> initRange;
		getline(parton1, line);
		istringstream iss(line);
		iss >> iEvent1 >> np1 >> initPT;
		if (argc == 4) {
			//            parton2 >> iEvent2 >> np2 >> dummyF >> dummyF;
			getline(parton2, line);
			istringstream iss2(line);
			iss2 >> iEvent2 >> np2 >> dummyF;
		}
		//        parton2 >> iEvent2 >> np2;
		//        np2=0;
		if(parton1.eof()) break;
		if (argc == 4) {
			if( parton2.eof() ) break;
		        if(iEvent1!=iEvent2) {
		           cout << "ERROR!!!! Event ID's don't match between positive and negative ... in jetShape code" << endl;
			   break;
		        }
		}
		nEvent += 1; 
		//        cout << iEvent1 << "  " << iEvent2 << "  " << np1 << "  " << np2 << "  " << dummyF << endl;

		vector<double> e, px, py, pz, eta;
		vector<int> id;

		double ie=0.0 , ipx=0.0 , ipy=0.0 , ipz=0.0;
		int iid=0;// , icat=0 ;

		for (unsigned i = 0; i < np1; ++i) {
			//            parton1 >> dummyInt >> iid >> ipx >> ipy >> ipz >> ie >> dummyF >> dummyF >> dummyF >> dummyF >> dummyF;
			getline(parton1, line);
			istringstream iss1(line);
			iss1 >> dummyInt >> iid >> ipx >> ipy >> ipz >> ie >> dummyF >> dummyF >> dummyF >> dummyF >> dummyF;

			id.push_back(iid);
			e.push_back(ie);
			px.push_back(ipx);
			py.push_back(ipy);
			pz.push_back(ipz);
		}

		if (argc == 4) {
			for (unsigned i = 0; i < np2; ++i) {
				//                parton2 >> dummyInt >> iid >> ipx >> ipy >> ipz >> ie >> dummyF >> dummyF >> dummyF >> dummyF >> dummyF;
				getline(parton2, line);
				istringstream iss2(line);
				iss2 >> dummyInt >> iid >> ipx >> ipy >> ipz >> ie >> dummyF >> dummyF >> dummyF >> dummyF >> dummyF;

				id.push_back(iid);
				e.push_back(-ie);  // modified by Shanshan, put - sign
				px.push_back(-ipx);
				py.push_back(-ipy);
				pz.push_back(-ipz);
			}
		}

		vector<fastjet::PseudoJet> input_particles;
		for (unsigned k0 = 0; k0 < np1 + np2; ++k0) {
			fastjet::PseudoJet pj;
			//            cout << px[k0] << "  " << py[k0] << endl;
			//            tot_ener=tot_ener+e[k0];
			//            tot_px=tot_px+px[k0];
			//            tot_py=tot_py+py[k0];
			//            tot_pz=tot_pz+pz[k0];
			if (e[k0] > 0) {
				pj.reset_momentum(px[k0], py[k0], pz[k0], e[k0]);
				pj.set_user_index(k0);
			}
			else  {
				pj.reset_momentum(-px[k0], -py[k0], -pz[k0], -e[k0]);
				pj.set_user_index(-k0);
			}

			double _pt= sqrt(px[k0]*px[k0] + py[k0]*py[k0]);

			//if (px[k0]*px[k0] + py[k0]*py[k0] > pt_cut*pt_cut) {
			//#LY 2021.09.15 if(_pt>pt_cut_min && _pt<pt_cut_max){
			if(_pt>pt_cut_min){//#LY 2021.09.15 
				double _eta=0.5*log((e[k0]+pz[k0])/(e[k0]-pz[k0]));
				if(fabs(_eta)<1.0){
					input_particles.push_back(pj); 
				}
			}

		}

		//        cout << "size of input particles: " << input_particles.size() << endl;

		//....Create a jet definition: a jet algorithm with a given radius parameter
		fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

		//....Create an instance of the negative energy recombiner, with a given flag ui
		NegativeEnergyRecombiner uir(-123456);
		//....Tell jet_def to use this new recombiner
		jet_def.set_recombiner(&uir);


		//      //....Run the jet clustering with the above jet definition
		//      fastjet::ClusterSequence clust_seq(input_particles, jet_def);

		double maxrap = 1.0;
		unsigned int n_repeat = 1; // default is 1
		double ghost_area = 0.01; // this is the default
		fastjet::GhostedAreaSpec area_spec(maxrap, n_repeat, ghost_area);

		//fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);
		fastjet::AreaDefinition area_def(active_area_explicit_ghosts, area_spec);
		// run the jet clustering with the above jet and area definitions
		//
		// The only change is the usage of a ClusterSequenceArea rather than
		//a ClusterSequence
		//----------------------------------------------------------
		fastjet::ClusterSequenceArea clust_seq(input_particles, jet_def, area_def);

		//....Get the resulting jets ordered in pt
		vector<PseudoJet> jets;
		double pt_Min = 0.2;
		vector<PseudoJet> jets_cs = sorted_by_pt(clust_seq.inclusive_jets(pt_Min));
		//select jet with eta<1-R
		Selector Fiducial_cut_selector = SelectorAbsEtaMax(1.0 - R);
		jets = Fiducial_cut_selector(jets_cs);

		//background substract_____________________________________________
		JetDefinition jet_def_bkgd(kt_algorithm, R);
		NegativeEnergyRecombiner uir_bkg(-123456);
		jet_def_bkgd.set_recombiner(&uir_bkg);
		AreaDefinition area_def_bkgd(active_area_explicit_ghosts,GhostedAreaSpec(maxrap, n_repeat, ghost_area));
		//select eta<1,not the hardest
		int Remove_N_hardest=1;
		Selector selector = SelectorAbsEtaMax(1.0) * (!SelectorNHardest(Remove_N_hardest));
		//Background estimation
		JetMedianBackgroundEstimator bkgd_estimator(selector, jet_def_bkgd, area_def_bkgd);
		bkgd_estimator.set_particles(input_particles);

		//setting subtractor
		Subtractor subtractor(&bkgd_estimator);
		double jets_rho = bkgd_estimator.rho();         //median of pt/area


		//        cout << "size of jets: " << jets.size() << endl;
		//        if (jets.size() == 0 ) continue;
		njets += jets.size();

		// print the jets
		//        cout << "evt: " << iEvent1 << "  pt y phi" << endl;
		//        cout << "Initial Energy: " << tot_ener << "  " << tot_px << "  " << tot_py << "  " << tot_pz << endl;
		//        tot_ener=0.0;
		//        tot_px=0.0;
		//        tot_py=0.0;
		//        tot_pz=0.0;

		int nJet=0;
		for (unsigned i = 0; i < jets.size(); i++) {

			double jet_Area = jets[i].area();
			double jet_pT_sub= jets[i].perp() - jet_Area*jets_rho;

			if(jet_Area<AreaCut) continue;

			if(jets[i].user_index()<0) hpt->Fill(jet_pT_sub,-1);
			else hpt->Fill(jet_pT_sub,1);

			if(jet_pT_sub>pTsubCut) nJet++;
			else continue;


			PtSum+=GetPtSum(jets[i],input_particles,hPtSumVsdr);
			Fill4dNtrk(jets[i],input_particles,hNtrkVsdr);

			PtJetConSum+=GetPtSum(jets[i],hPtJetSumVsdr);
			Fill4dNtrk(jets[i],hNtrkJetVsdr);

			PtInConeSum+=GetPtInConeSum(jets[i],R,input_particles);
			PtJetSum+=jet_pT_sub;

			//vector<PseudoJet> _constituents = jets[i].constituents();
			//int nCon = _constituents.size();
			//int Nlead=0;
			//for(int i_ncon=0; i_ncon<nCon; i_ncon++){
			//if(_constituents[i_ncon].perp()>5.0) Nlead++;
			//}
			//if(Nlead<1) continue;
		}
		JetsSum+=nJet;

		//        cout << "done" << endl;

	}
	
	TH1D *hSums = new TH1D("hSums","hSums",6,0,6);
	hSums->GetXaxis()->SetBinLabel(1,"N_{evt}");
	hSums->GetXaxis()->SetBinLabel(2,"N_{jets}");
	hSums->GetXaxis()->SetBinLabel(3,"#Sigma_{jets}#Sigma_{trk}p_{T}^{trk}");
	hSums->GetXaxis()->SetBinLabel(4,"#Sigma_{jets}#Sigma_{trk}p_{T}^{jet trk}");
	hSums->GetXaxis()->SetBinLabel(5,Form("#Sigma_{jets}#Sigma_{trk#in[0,R=%.1f]}p_{T}^{trk}",R));
	hSums->GetXaxis()->SetBinLabel(6,"#Sigma_{jets}p_{T}^{jet}");

	hSums->SetBinContent(1,nEvent);
	hSums->SetBinContent(2,JetsSum);
	hSums->SetBinContent(3,PtSum);
	hSums->SetBinContent(4,PtJetConSum);
	hSums->SetBinContent(5,PtInConeSum);
	hSums->SetBinContent(6,PtJetSum);

	// Scale histogram for quick single file check
	// hrho & hYield for single root file only, will become meanless after merging files with hadd
	TH1D *hrho = (TH1D*)hPtSumVsdr->Clone("hrho");
	hrho->SetTitle(Form("hrho;#Deltar;#frac{1}{#scale[1.2]{#Sigma_{jets}#Sigma_{trk#in[0,R=%.1f]}}p_{T}^{trk}} #frac{#scale[1.2]{#Sigma_{jets}}dp_{T}^{trk}}{d#Deltar}",R));
	TH1D *hYield = (TH1D*)hNtrkVsdr->Clone("hYield");
	hYield->SetTitle("hYield;#Deltar;#frac{1}{N_{jets}} #frac{d^{2}N_{trk}}{d#Deltardp_{T}^{trk}}");
	TH1D *hjrho = (TH1D*)hPtJetSumVsdr->Clone("hjrho");
	hjrho->SetTitle("hjrho;#Deltar;#frac{1}{#scale[1.2]{#Sigma_{jets}#Sigma_{trk}}p_{T}^{jet trk}} #frac{#scale[1.2]{#Sigma_{jets}}dp_{T}^{jet trk}}{d#Deltar}");
	TH1D *hjYield = (TH1D*)hNtrkJetVsdr->Clone("hjYield");
	hjYield->SetTitle("hjYield;#Deltar;#frac{1}{N_{jets}} #frac{d^{2}N_{jet trk}}{d#Deltardp_{T}^{trk}}");
	double DeltaPt = 30-0.2;	// pt range for yield & rho calculation. the default value. Need to check what is used
	hYield->Scale(JetsSum>0 ? 1./(JetsSum*DeltaPt) : 1);
	hYield->Scale(1,"width");	// option "width": bin contents and errors divided by the bin width
	hrho->Scale(PtInConeSum>0 ? 1./PtInConeSum : 1);
	hrho->Scale(1,"width");
	hjYield->Scale(JetsSum>0 ? 1./(JetsSum*DeltaPt) : 1);
	hjYield->Scale(1,"width");	// option "width": bin contents and errors divided by the bin width
	hjrho->Scale(PtJetConSum>0 ? 1./PtJetConSum : 1);
	hjrho->Scale(1,"width");


	//....save to ROOT file
	string outputRootFile = outputFile;
	size_t findroot = outputRootFile.find(".root");
	ostringstream oscut;
	oscut << "_jetpT";
	oscut.width(2);
	oscut.fill('0');
	oscut << pTsubCut;
	oscut.width(1);
	oscut << "R" << R*10;
	string scut = oscut.str();
	if(findroot!=string::npos) outputRootFile.insert(findroot,scut);
	if(findroot==string::npos) {outputRootFile+=scut; outputRootFile+=".root";}
	TFile *fout = new TFile(outputRootFile.c_str(),"RECREATE");
	hpt->Write();
	hNtrkVsdr->Write();
	hPtSumVsdr->Write();
	hNtrkJetVsdr->Write();
	hPtJetSumVsdr->Write();
	hSums->Write();
	hYield->Write();
	hjrho->Write();
	hjYield->Write();
	hrho->Write();
	fout->Close();

	//....time counting ends
	struct tm *local_end;
	time_t time_end;
	time_end=time(NULL);
	local_end=localtime(&time_end);

	char buf2[80];
	strftime(buf2,80,"Current Time: %Y-%m-%d %H:%M:%S",local_end);
	cout << "the program ends at:" << endl;
	cout << buf2 << endl;

	unsigned cost,nh,nm,ns;
	cost=difftime(time_end,time_start);

	nh=cost/3600;
	nm=(cost%3600)/60;
	ns=(cost%3600)%60;

	cout << "the program costs:" << endl;
	cout << cost << "s:" << " " << nh << "h" << " " << nm << "m" << " " << ns << "s" << endl;

	return 0;
}
