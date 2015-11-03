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

//int main(int argc, char *argv[]) {
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
    pdt.mayDecay(13,false);     // muon
    pdt.mayDecay(-13,false);
    pdt.mayDecay(211,false);   // pi+
    pdt.mayDecay(-211,false);  // pi-
    pdt.mayDecay(130,false);    // KL
    pdt.mayDecay(311,false);    // K0
    pdt.mayDecay(-311,false);   // K0_bar
    pdt.mayDecay(321,false);    // K+
    pdt.mayDecay(-321,false);   // K-
    pdt.mayDecay(311,false);    // K0
    pdt.mayDecay(-311,false);   // K0_bar
    pdt.mayDecay(2112,false);  // n
    pdt.mayDecay(-2112,false); // n_bar
  }
  pythia.init();

  // Allow for possibility of a few faulty events.
  // Extract settings to be used in the main program.
  //  int nAbort  = pythia.mode("Main:timesAllowErrors");
  int nAbort = 10;
  int iAbort = 0;

  //  double nevents = 90000.;
  //  int inevents = pythia.mode("Main:numberOfEvents");;
  double nevents = 0;

  int nbins;
  double xmin;
  double xmax;
  if (idecay == 0) {
    xmin = 0.01;
    xmax = 100.;
    nbins = 72;
  } else {
    xmin = 0.01;
    xmax = 1000000.;
    nbins = 800;
  }
  double *x = new double[nbins];
  for (int i = 0; i < nbins+1; ++i) {
    //    x[i] = xmin +(xmax -xmin)/nbins*i;
    x[i] = pow(10,log10(xmin) +(log10(xmax) -log10(xmin))/float(nbins)*i);
  }

  double me = 0.000511;
  double mp = 0.938272;
  double mn = 0.939565;
  double mpi = 0.13957018;
  double mkpm = 0.493677;
  double mkL = 0.497614; // K0 mass

  int *nnGam = new int[nbins];
  int *nnEle = new int[nbins];
  int *nnNeue = new int[nbins];
  int *nnNeumu = new int[nbins];
  int *nnNeutau = new int[nbins];
  int *nnP = new int[nbins];
  int *nnaP = new int[nbins];
  int *nnN = new int[nbins];
  int *nnaN = new int[nbins];
  int *nnpip = new int[nbins];
  int *nnpim = new int[nbins];
  int *nnKp = new int[nbins];
  int *nnKm = new int[nbins];
  int *nnKL = new int[nbins];
  int *nnK0 = new int[nbins];
  for (int i = 0; i < nbins+1; ++i) {
    nnGam[i] = 0;
    nnEle[i] = 0;
    nnNeue[i] = 0;
    nnNeumu[i] = 0;
    nnNeutau[i] = 0;
    nnP[i] = 0;
    nnaP[i] = 0;
    nnN[i] = 0;
    nnaN[i] = 0;
    nnpip[i] = 0;
    nnpim[i] = 0;
    nnKp[i] = 0;
    nnKm[i] = 0;
    nnKL[i] = 0;
    nnK0[i] = 0;
  }
  int totGam = 0;
  int totEle = 0;
  int totNeue = 0;
  int totNeumu = 0;
  int totNeutau = 0;
  int totP = 0;
  int totaP = 0;
  int totN = 0;
  int totaN = 0;
  int totpip = 0;
  int totpim = 0;
  int totKp = 0;
  int totKm = 0;
  int totKL = 0;
  int totK0 = 0;
  int totOther = 0;

  // Begin event loop; generate until none left in input file.
  for (int iEvent = 0; iEvent < inevents; ++iEvent) {
  //for (int iEvent = 0; ; ++iEvent) {
    // Set up parton-level configuration.
    //    fillPartons( 91.188, event, pdt, pythia.rndm);

    // Generate events, and check whether generation failed.
    if (!pythia.next()) {
      // If failure because reached end of file then exit event loop.
      if (pythia.info.atEndOfFile()) break;
      // First few failures write off as "acceptable" errors, then quit.
      if (++iAbort < nAbort) continue;
      break;
    }

    //    pythia.event.list();

    double Ekin;
    // Fill in histogram.
    for (int i = 0; i < event.size(); ++i) {
      if (event[i].isFinal()) {
	if (event[i].id() == 22) {
	  Ekin = event[i].e();
	  totGam = totGam +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnGam[ii] = nnGam[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (abs(event[i].id()) == 11) {
	  Ekin = event[i].e() -me;
	  totEle = totEle +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnEle[ii] = nnEle[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (abs(event[i].id()) == 12) {
	  Ekin = event[i].e();
	  totNeue = totNeue +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnNeue[ii] = nnNeue[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (abs(event[i].id()) == 14) {
	  Ekin = event[i].e();
	  totNeumu = totNeumu +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnNeumu[ii] = nnNeumu[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (abs(event[i].id()) == 16) {
	  Ekin = event[i].e();
	  totNeutau = totNeutau +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnNeutau[ii] = nnNeutau[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (event[i].id() == 2212) {
	  Ekin = event[i].e() -mp;
	  totP = totP +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnP[ii] = nnP[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (event[i].id() == -2212) {
	  Ekin = event[i].e() -mp;
	  totaP = totaP +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnaP[ii] = nnaP[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (event[i].id() == 2112) {
	  Ekin = event[i].e() -mn;
	  totN = totN +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnN[ii] = nnN[ii] +1;
		//		cout << "neutron is generated!!!" << "" << nnN[ii] << endl;
		//		pythia.event.list();
		break;
	      }
	    }
	  }
	}
	else if (event[i].id() == -2112) {
	  Ekin = event[i].e() -mn;
	  totaN = totaN +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnaN[ii] = nnaN[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (event[i].id() == 211) {
	  Ekin = event[i].e() -mpi;
	  totpip = totpip +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnpip[ii] = nnpip[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (event[i].id() == -211) {
	  Ekin = event[i].e() -mpi;
	  totpim = totpim +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnpim[ii] = nnpim[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (event[i].id() == 321) {
	  Ekin = event[i].e() -mkpm;
	  totKp = totKp +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnKp[ii] = nnKp[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (event[i].id() == -321) {
	  Ekin = event[i].e() -mkpm;
	  totKm = totKm +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnKm[ii] = nnKm[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (event[i].id() == 130) {
	  Ekin = event[i].e() -mkL;
	  totKL = totKL +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnKL[ii] = nnKL[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (event[i].id() == 311) {
	  Ekin = event[i].e() -mkL;
	  totK0 = totK0 +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnK0[ii] = nnK0[ii] +1;
		break;
	      }
	    }
	  }
	}
	else if (event[i].id() == -311) {
	  Ekin = event[i].e() -mkL;
	  totK0 = totK0 +1;
	  for (int ii = 0; ii < nbins; ++ii) {
	    if (x[ii] < Ekin) {
	      if (Ekin < x[ii+1]) {
		nnK0[ii] = nnK0[ii] +1;
		break;
	      }
	    }
	  }
	}
	else {
	  totOther = totOther +1;
	}
      }
    }
    nevents = nevents +1;
  // End of event loop.
  }
  // Give statistics. Print histogram.
  pythia.stat();
  
  std::ofstream ofsnp("np_sptrm.dat");
  for (int i = 0; i < nbins+1; ++i) {
    ofsnp << setiosflags(ios::scientific) << ECM << " " << x[i] << " " << nnN[i]/nevents << " " << nnP[i]/nevents << " " 
	    << nnpip[i]/nevents << " " << nnpim[i]/nevents << " " << nnKp[i]/nevents << " "
	  << nnKm[i]/nevents << " " << (nnKL[i]+nnK0[i]/2.)/nevents << " " << nnaN[i]/nevents << " " 
            << nnaP[i]/nevents << endl;
  }

  std::ofstream ofsnpnorm("np_sptrm_norm.dat");
  for (int i = 0; i < nbins+1; ++i) {
    double dx = x[i+1] -x[i];
    ofsnpnorm << setiosflags(ios::scientific) << ECM << " " << x[i] << " " << nnN[i]/nevents/dx << " " << nnP[i]/nevents/dx << " " 
	    << nnpip[i]/nevents/dx << " " << nnpim[i]/nevents/dx << " " << nnKp[i]/nevents/dx << " "
	  << nnKm[i]/nevents/dx << " " << (nnKL[i]+nnK0[i]/2.)/nevents/dx << " " << nnaN[i]/nevents/dx << " " 
            << nnaP[i]/nevents/dx << endl;
  }

  std::ofstream ofsEdist("Edist.dat");
  ofsEdist << setiosflags(ios::scientific) << (totpip+totpim)/nevents << " " << (totN+totaN)/nevents << " " << (totP+totaP)/nevents << endl;
  for (int i = 0; i < nbins+1; ++i) {
    ofsEdist << setiosflags(ios::scientific) << x[i] << " " << (nnpip[i]+nnpim[i])/nevents << " " << (nnN[i]+nnaN[i])/nevents << " " << (nnP[i]+nnaP[i])/nevents << endl;
  }

  std::ofstream ofsnini("nini.dat");
  ofsnini << setiosflags(ios::scientific) << ECM << " " << totN/nevents << " " << totP/nevents << " " << totpip/nevents << " " 
	  << totpim/nevents << " " << totKp/nevents << " " << totKm/nevents << " " << (totKL+totK0/2.)/nevents << " " 
	  << totaN/nevents << " " << totaP/nevents << endl;

  cout << setiosflags(ios::scientific) << totGam/nevents << " " << totEle/nevents << " " << totN/nevents << " " << totP/nevents << " " 
       << totpip/nevents << " " << totpim/nevents << " " << totKp/nevents << " " << totKm/nevents << " " 
       << (totKL+totK0/2.)/nevents << " " << totaN/nevents << " " << totaP/nevents << " " << totOther/nevents << endl;


////////////////           ////////////////////// 
//////////////// Evis data ////////////////////// 
  if (idecay != 0) {
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
    pdt.mayDecay(311,true);    // K0
    pdt.mayDecay(-311,true);   // K0_bar
    pdt.mayDecay(2112,false);  // n
    pdt.mayDecay(-2112,false); // n_bar
  }
  pythia.init();

  nevents = 0;
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
