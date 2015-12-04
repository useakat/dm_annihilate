*HOW TO USE

1.prepare param_card_temp.dat:

<Using SUSYHIT>
  cd susyhit
  modify susyhit.in files
  ./run.sh

<Using SPheno>
  cd SUSYspc
  SPheno LesHouches.in

<Using SOFTSUSY>
  cd SUSYspc
  softpoint.x leshouches < [softsusy].in > [softsusy].spc

<modify output SLHA file>
  add 1000049 MMGRV # Mgrv to MASS COMMON BLOCK
  add the following lines to the end of the file
    DECAY 1000049 0.000000e+00 # Wgrv
    #===========================================================
    # QUANTUM NUMBERS OF NEW STATE(S) (NON SM PDG CODE)
    #===========================================================
    Block QNUMBERS 1000049  # grv
            1 0  # 3 times electric charge
            2 4  # number of spin states (2S+1)
            3 1  # colour rep (1: singlet, 3: triplet, 8: octet)
            4 0  # Particle/Antiparticle distinction (0=own anti)

<Comments>
Please check your parameter card so that MODSEL and MINPAR is set 
appropreately. Inappropreate settings would cause many abortions in
Pythia run.

In DECAY blocks, width should not be too small, otherwise pythia 
doesn't work correctly (width should be larger at least 10^-12 GeV)

It is highly recommended to check allprocess.log in par_n folder and 
ensure that there is not many abortion or erros in Pythia run. 


2. Modify relevant files
   * run_grv_decay_parallel.sh
      chi0mass: change to the LSP mass
      nevents (L35 & L38)

   * run_grv_decay.sh
      set param_card_XXX.dat appropreately (~L21)

   * run_mass-parallel.sh
      modify parameters part
   
3. Run
   ./bin/run_mass-parallel.sh
   * Do not generate more than one billion events in general DM mode, otherwise pythia8 couldn't terminate event generation...

3. Notes

