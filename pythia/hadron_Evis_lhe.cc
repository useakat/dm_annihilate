#include <iostream>
#include "Pythia8/Pythia.h"

using namespace Pythia8;

//int main(int argc, char *argv[]) {
int main(int argc, char *argv[]) {
  Pythia pythia;
  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;

  char *ends;
  double mDM = strtod(argv[1], &ends);
  //  int inevents = atoi(argv[2]);

  std::ofstream ofsmass("mass.dat");
  ofsmass << "mass " << argv[1] << " " << mDM << endl;
  // Initialize Les Houches Event File run. List initialization information.
  //pythia.readString("ProcessLevel:all = off");
  //pythia.readString("PartonLevel:all = on");
  //  pythia.readString("HadronLevel:Decay = off");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0"); //pick new random seed for each run, based on clock

/////////////// LHE input mode ////////////////
  pythia.readString("Beams:frameType = 4");
  pythia.readString("Beams:LHEF = unweighted_events.lhe");
  //  pythia.readString("SLHA:useDecayTable = true");

/////////////// Generic resonance mode ////////////////
  // // A class to generate the fictitious resonance initial state.
  // SigmaProcess* sigma1GenRes = new Sigma1GenRes();
  // // Hand pointer to Pythia.
  // pythia.setSigmaPtr( sigma1GenRes);
  // // Read in the rest of the settings and data from a separate file.  
  // pythia.settings.parm("Beams:eCM", 2*mDM);
  // pythia.readFile("generic_resonance.cmnd");

//////////////// Evis data ////////////////////// 

//// Life time < 10^-10 s ////////////
  pdt.mayDecay(15,true);     // tau
  pdt.mayDecay(-15,true);
  pdt.mayDecay(310,true);    // KS
  pdt.mayDecay(111,true);    // pi0
  pdt.mayDecay(311,true);    // K0
  pdt.mayDecay(-311,true);   // K0_bar
  pdt.mayDecay(3122,true);   // Lambda
//// Life time > 10^-10 s ////////////
  pdt.mayDecay(13,true);     // muon
  pdt.mayDecay(-13,true);
  pdt.mayDecay(211,true);   // pi+
  pdt.mayDecay(-211,true);  // pi-
  pdt.mayDecay(130,true);    // KL
  pdt.mayDecay(321,true);    // K+
  pdt.mayDecay(-321,true);   // K-
  pdt.mayDecay(2112,false);  // n
  pdt.mayDecay(-2112,false); // n_bar

  pythia.init();

  // Allow for possibility of a few faulty events.
  // Extract settings to be used in the main program.
  //  int nAbort  = pythia.mode("Main:timesAllowErrors");
  int nAbort  = 10;
  int iAbort = 0;

  //  double nevents = 90000.;
  //  int inevents = pythia.mode("Main:numberOfEvents");;
  double nevents = 0;

  double me = 0.000511;
  double mp = 0.938272;
  double mn = 0.939565;
  //  double mpi = 0.13957018;
  //double mkpm = 0.493677;
  //double mkL = 0.497614; // K0 mass

  double Ekin = 0;
  double Ekin_neu = 0;
  double Ekin_other = 0;
  // Begin event loop; generate until none left in input file.
  //  for (int iEvent = 0; iEvent < inevents; ++iEvent) {
  for (int iEvent = 0; ; ++iEvent) {
    // Generate events, and check whether generation failed.
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;
      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }

    for (int i = 0; i < event.size(); ++i) {
      if (event[i].isFinal()) {
	if (abs(event[i].id()) == 12) {
	  Ekin_neu = Ekin_neu +event[i].e();
	}
	else if (abs(event[i].id()) == 14) {
	  Ekin_neu = Ekin_neu +event[i].e();
	}
	else if (abs(event[i].id()) == 16) {
	  Ekin_neu = Ekin_neu +event[i].e();
	}
	else if (abs(event[i].id()) == 11) {
	  Ekin = Ekin +event[i].e() -me;
	  //Ekin = Ekin +event[i].e();
	}
	else if (abs(event[i].id()) == 22) {
	  Ekin = Ekin +event[i].e();
	}
	else if (abs(event[i].id()) == 2112) {
	  Ekin = Ekin +event[i].e() -mn;
	  //Ekin = Ekin +event[i].e();
	}
	else if (abs(event[i].id()) == 2212) {
	  Ekin = Ekin +event[i].e() -mp;
	  //Ekin = Ekin +event[i].e();
	}
	else {
	  Ekin_other = Ekin_other +event[i].e();
	}
      }
    }
    nevents = nevents +1;
  // End of event loop.
  }

  // Give statistics. Print histogram.
  pythia.stat();

  std::ofstream ofsEvis("Evis.dat");
  ofsEvis << fixed << 2*mDM << " " << Ekin/nevents << " "  << Ekin/(2*mDM)/nevents << endl;

  cout << 2*mDM << " " << Ekin/nevents << " "  << Ekin/(2*mDM)/nevents << " " << Ekin_other/nevents << " " << nevents << endl;

  return 0;
}
