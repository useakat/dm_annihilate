#include <iostream>
#include <iomanip>
#include "Pythia8/Pythia.h"

using namespace Pythia8;

class Sigma1GenRes : public Sigma1Process {

public:

  // Constructor.
  Sigma1GenRes() {}

  // Evaluate sigmaHat(sHat): dummy unit cross section.
  virtual double sigmaHat() {return 1.;}

  // Select flavour. No colour or anticolour.
  virtual void setIdColAcol() {setId( -11, 11, 999999);
    setColAcol( 0, 0, 0, 0, 0, 0);}

  // Info on the subprocess.
  virtual string name()    const {return "GenericResonance";}
  virtual int    code()    const {return 9001;}
  virtual string inFlux()  const {return "ffbarSame";}

};

int main(int argc, char *argv[]) {
  Pythia pythia;
  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;
  char *ends;
  int irun_mode = atoi(argv[1]);
  double mDM = strtod(argv[2], &ends);
  int inevents = atoi(argv[3]);
  int idecay = atoi(argv[4]);

  std::ofstream ofsmass("mass.dat");
  ofsmass << mDM << endl;

  double ECM = 0;
  if (irun_mode == 1){
    ECM = mDM; // DM decay process
  } else if (irun_mode == 2){
    ECM = 2*mDM; // DM annihilation process
  } else {
    cout << "ERROR: run_mode is invalid value. 1 or 2 is allowed." << endl;
  } 

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = 0"); //pick new random seed for each run, based on clock
  pythia.readString("Next:numberShowEvent = 10"); // print event record n times
  pythia.readString("Next:numberCount = 0"); // Print a line telling how many events have been generated so far, once every numberCount events. If set zero then no lines are ever printed. 
  //  pythia.readString("Tune:ee = 2");

  //pythia.readString("ProcessLevel:all = off");
  //pythia.readString("ProcessLevel:resonanceDecays = off");
  //pythia.readString("PartonLevel:all = off");
  //  pythia.readString("PartonLevel:MPI = off");
  //  pythia.readString("PartonLevel:FSRinResonances = off");
  //  pythia.readString("PartonLevel:ISR = off");
  //pythia.readString("PartonLevel:FSR = off");
  //pythia.readString("PartonLevel:FSRinProcess = off");
  //  pythia.readString("PartonLevel:Remnants = off");
  //  pythia.readString("HadronLevel:Decay = off");
  //  pythia.readString("PDFinProcess:nQuarkIn = 1");
  pythia.readString("PhaseSpace:useBreitWigners = off");
  pythia.readString("PhaseSpace:mHatMin = 2.");

/////////////// LHE input mode ////////////////
//  pythia.readString("Beams:frameType = 4");
  // pythia.readString("Beams:LHEF = unweighted_events_mDM10000.lhe");
  //pythia.readString("SLHA:useDecayTable = false");

  //pythia.readString("Beams:frameType = 1");

  pythia.readString("Beams:frameType = 2");
  pythia.readString("PDF:lepton = off"); 
  pythia.readString("Beams:idA = 11");
  pythia.readString("Beams:idB = -11");
  pythia.settings.parm("Beams:eA", ECM/2.);
  pythia.settings.parm("Beams:eB", ECM/2.);

  //  pythia.settings.parm("Beams:eCM", ECM);
  //  pythia.readString("11:m0 = 5000.");
  //pythia.readString("-11:m0 = 5000.");

/////////////// Generic resonance mode ////////////////
  if (idecay == 0) {
    // test mode with Z resonance
    pythia.readString("WeakSingleBoson:ffbar2ffbar(s:gmZ) = on"); 
    pythia.readString("WeakZ0:gmZmode = 2");
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 1 2 3 4 5");
  } else {   
    // generic resonance decay mode
    SigmaProcess* sigma1GenRes = new Sigma1GenRes();
    pythia.setSigmaPtr(sigma1GenRes);
    pythia.readFile("generic_resonance.cmnd");
  }

/////////////// Pythia internal process mode ////////////////
  // if (idecay == 24) {
  //   pythia.readString("WeakDoubleBoson:ffbar2WW = on");
  // } else if (idecay == 23) {
  //   pythia.readString("WeakDoubleBoson:ffbar2gmZgmZ = on");
  //   pythia.readString("WeakZ0:gmZmode = 2");
  // } else {
  //   pythia.readString("WeakSingleBoson:ffbar2ffbar(s:gmZ) = on"); 
  //   pythia.readString("WeakZ0:gmZmode = 2");
  //   pythia.readString("23:onMode = off");
  //   if (idecay == 15) {
  //     pythia.readString("23:onIfAny = 15");
  //   } else if (idecay == 1) {
  //     pythia.readString("23:onIfAny = 1");
  //   } else if (idecay == 2) {
  //     pythia.readString("23:onIfAny = 2");
  //   } else if (idecay == 3) {
  //     pythia.readString("23:onIfAny = 3");
  //   } else if (idecay == 4) {
  //     pythia.readString("23:onIfAny = 4");
  //   } else if (idecay == 5) {
  //     pythia.readString("23:onIfAny = 5");
  //   } else if (idecay == 0) {
  //     pythia.readString("23:onIfAny = 1 2 3 4 5");
  //   }   
  // } 

  if (idecay == 0) {
    pdt.mayDecay(13,false);     // muon
    pdt.mayDecay(-13,false);
    pdt.mayDecay(15,true);     // tau
    pdt.mayDecay(-15,true);
    pdt.mayDecay(130,false);    // KL
    pdt.mayDecay(310,true);    // KS
    pdt.mayDecay(111,true);    // pi0
    pdt.mayDecay(211,false);   // pi+
    pdt.mayDecay(-211,false);  // pi-
    pdt.mayDecay(311,false);    // K0
    pdt.mayDecay(-311,false);   // K0_bar
    pdt.mayDecay(321,false);    // K+
    pdt.mayDecay(-321,false);   // K-
    pdt.mayDecay(2112,false);  // n
    pdt.mayDecay(-2112,false); // n_bar
    pdt.mayDecay(3122,true);   // Lambda
  } else {
    //// Life time < 10^-10 s ////////////
    pdt.mayDecay(15,true);     // tau
    pdt.mayDecay(-15,true);
    pdt.mayDecay(310,true);    // KS
    pdt.mayDecay(111,true);    // pi0
    pdt.mayDecay(3122,true);   // Lambda
    //// Life time > 10^-10 s ////////////
    pdt.mayDecay(13,true);     // muon
    pdt.mayDecay(-13,true);
    pdt.mayDecay(211,true);   // pi+
    pdt.mayDecay(-211,true);  // pi-
    pdt.mayDecay(130,true);    // KL
    pdt.mayDecay(311,true);    // K0
    pdt.mayDecay(-311,true);   // K0_bar
    pdt.mayDecay(321,true);    // K+
    pdt.mayDecay(-321,true);   // K-
    pdt.mayDecay(2112,false);  // n
    pdt.mayDecay(-2112,false); // n_bar
  }
  double me = 0.000511;
  double mp = 0.938272;
  double mn = 0.939565;
  double mpi = 0.13957018;
  double mkpm = 0.493677;
  double mkL = 0.497614; // K0 mass

  pythia.init();

  // Allow for possibility of a few faulty events.
  // Extract settings to be used in the main program.
  //  int nAbort  = pythia.mode("Main:timesAllowErrors");
  int nAbort = 10;
  int iAbort = 0;

  double nevents = 0;
  double Ekin = 0;
  double Ekin_neu = 0;
  double Ekin_other = 0;
  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; iEvent < inevents; ++iEvent) {
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
	}
	else if (abs(event[i].id()) == 22) {
	  Ekin = Ekin +event[i].e();
	}
	else if (abs(event[i].id()) == 2112) {
	  Ekin = Ekin +event[i].e() -mn;
	}
	else if (abs(event[i].id()) == 2212) {
	  Ekin = Ekin +event[i].e() -mp;
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
  ofsEvis << setiosflags(ios::scientific) << ECM << " " << Ekin/nevents << " "  << Ekin/ECM/nevents << endl;

  cout << setiosflags(ios::scientific) << ECM << " " << Ekin/nevents << " "  << Ekin/ECM/nevents << " " << Ekin_other/nevents << " " << nevents << endl;

  std::ofstream ofsGamE("Evis_tot.dat");
  ofsGamE << Ekin/nevents << endl;

  return 0;
}
