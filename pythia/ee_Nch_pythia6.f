      PROGRAM MAIN81
C...Exle program to illustrate how to read in LHEF files, e.g., from 
C...a matrix-element generator, into Pythia 6 for further processing 
C...(resonance decays, parton showering, underlying event, and hadronization)
C...     For LHEF, see: hep-ph/0609017
C...     For LHEF+SLHA/BSM, see also: arXiv:0712.3311

C--------------- PREAMBLE: COMMON BLOCK DECLARATIONS ETC -------------
C...All real arithmetic done in double precision.
      IMPLICIT DOUBLE PRECISION(A-H, O-Z)
C...The PYTHIA event record:
      COMMON/PYJETS/N,NPAD,K(4000,5),P(4000,5),V(4000,5)
C...Pythia parameters
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
C...Parameters.
      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
C...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
C...Decay information.
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
C...Selection of hard scattering subprocesses.
      COMMON/PYSUBS/MSEL,MSELPD,MSUB(500),KFIN(2,-40:40),CKIN(200)
C...EXTERNAL statement links PYDATA on most machines.
      CHARACTER*3 chlun
      EXTERNAL PYDATA

      integer i,j
      integer nbins,PYCOMP
      parameter (nbins=60)
      integer nnGam(nbins),nnEle(nbins),nnNeue(nbins),nnNeumu(nbins)
     &     ,nnNeutau(nbins),nnP(nbins),nnN(nbins)
      integer totGam,totEle,totNeue,totNeumu,totNeutau,totP,totN
      integer totOther
      real*8 xmin,xmax,x(nbins),dx
      real*8 me,mp,mn
      real*8 nevents
      real*8 ECM,Nch
      integer NEV,IDC
      character*10 cECM,cNEV

      integer MSTU,MSTJ
      double precision PARU,PARJ
      integer MDCY,MDME,KFDP
      double precision BRAT
C-------------------------- PYTHIA SETUP -----------------------------

C...Main parameters of run: c.m. energy and number of events.
      call getarg(1,cECM)
      call getarg(2,cNEV)
      read (cECM,*) ECM
      read (cNEV,*) NEV

c$$$C...Select gamma*/Z0 production process.
      MSEL=0     ! process are selected by process number
      MSUB(1)=1  ! gam/Z0 production
c$$$      MSTP(42) = 1 ! 0:resonance is on the mass shell 1:Breigt Wigner
c$$$      MSTP(43) = 3 ! 1:only gam*  2:only Z0*
c$$$
C...Only allow Z0 decay to quarks (i.e. no leptonic final states).
      DO 100 IDC = MDCY(23,2),MDCY(23,2)+MDCY(23,3)-1
        IF(IABS(KFDP(IDC,1)).gt.5) MDME(IDC,1)=MIN(0,MDME(IDC,1))
  100 CONTINUE

c$$$C...1) Open LHEF file on unit LUN, and tell Pythia where to find it.
c$$$      LUN=88
c$$$c      OPEN(LUN,FILE='unweighted_events_mod.lhe')
c$$$      OPEN(LUN,FILE='unweighted_events.lhe')
c$$$      WRITE(CHLUN,'(I3)') LUN
c$$$      CALL PYGIVE('MSTP(161)='//CHLUN)
c$$$      CALL PYGIVE('MSTP(162)='//CHLUN)

C...2) Initialize Pythia for user process  
      MDCY(PYCOMP(13),1) = 0
      MDCY(PYCOMP(15),1) = 1
      MDCY(PYCOMP(130),1) = 0
      MDCY(PYCOMP(310),1) = 1
      MDCY(PYCOMP(111),1) = 1
      MDCY(PYCOMP(211),1) = 0
      MDCY(PYCOMP(311),1) = 0
      MDCY(PYCOMP(321),1) = 0
      MDCY(PYCOMP(2112),1) = 0
c      MDCY(PYCOMP(221),1)=1
c      MDCY(PYCOMP(421),1)=1
c      MDCY(PYCOMP(511),1)=1
c      MDCY(PYCOMP(521),1)=1

c      MSTJ(104)=6 ! allow tt~ production
c      MSTJ(24)=0  ! discrete mass value
c      MSTP(122)=0 ! 
      MSTU(21)=1 ! check on possible errors during program execution
      MSTU(22)=0 ! maximum number of errors that are printed
      MSTU(25)=0 ! printing of warning messages
      CALL PYINRE
      CALL PYGIVE("MSTJ(41)=12")
      CALL PYINIT('CMS','e+','e-',ECM)
c      CALL PYINIT('USER',' ',' ',0D0)

C...Check that Z0 decay channels set correctly.
      CALL PYSTAT(2)
      
C------------------------- GENERATE EVENTS ---------------------------

C...Initial values for number of events and cumulative charged multiplicity
      Nch = 0
      nevents = 0

C...Get next event from file and process it
      DO 200 IEV=1,NEV
         CALL PYEVNT

C...If event generation failed, quit loop
         IF(MSTI(51).EQ.1) THEN
            GOTO 999
         ENDIF
         
C...  Print first few events, both LHEF input and Pythia output.
         IF(IEV.LE.3) THEN
            CALL PYLIST(7)
            CALL PYLIST(2)
         ENDIF
         
c     CALL PYEDIT(1) ! decayed particles are removed
         CALL PYEDIT(1)         ! neutral particles are removed
         do i = 1,N
            if (abs(K(i,2)).eq.211) then ! pi+-
               Nch = Nch +1
c            elseif (abs(K(i,2)).eq.2112) then ! neutron
c               Nch = Nch +1
            elseif (abs(K(i,2)).eq.2212) then ! proton
               Nch = Nch +1
c            elseif (abs(K(i,2)).eq.130) then ! K0L
c               Nch = Nch +1
c            elseif (abs(K(i,2)).eq.311) then ! K0
c               Nch = Nch +1
            elseif (abs(K(i,2)).eq.321) then ! K+-
               Nch = Nch +1
            endif
         enddo
         nevents = nevents +1
         
C...  Loop back to look for next event
 200  CONTINUE
         
C...  Print final statistics.
 999  CALL PYSTAT(1)
         
      open(1,file="Nch.dat",status="replace")
      write(1,*) ECM,Nch/nevents
      close(1)

      write(*,*) "counted evnets:", IEV
      write(*,*) "N charged particles:", Nch/nevents

      END
