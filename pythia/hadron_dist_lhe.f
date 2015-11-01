      PROGRAM MAIN81
C...Example program to illustrate how to read in LHEF files, e.g., from 
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
C...EXTERNAL statement links PYDATA on most machines.
      CHARACTER*3 chlun
      EXTERNAL PYDATA

      integer i,j
      integer nbins,PYCOMP
      parameter (nbins=800)
      integer nnGam(nbins),nnEle(nbins),nnNeue(nbins),nnNeumu(nbins)
     &     ,nnNeutau(nbins),nnP(nbins),nnaP(nbins),nnN(nbins)
     &     ,nnaN(nbins),nnpip(nbins),nnpim(nbins),nnKp(nbins)
     &     ,nnKm(nbins),nnKL(nbins)
      integer totGam,totEle,totNeue,totNeumu,totNeutau,totP,totaP
      integer totN,totaN,totpip,totpim,totKp,totKm,totKL,totOther
      real*8 xmin,xmax,x(nbins),dx
      real*8 me,mp,mn,mpi,mkpm,mkL
      real*8 nnevents,mDM

      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer MDCY,MDME,KFDP
      double precision BRAT

      mDM = 398d0
C-------------------------- PYTHIA SETUP -----------------------------

C...1) Open LHEF file on unit LUN, and tell Pythia where to find it.
      LUN=88
c      OPEN(LUN,FILE='unweighted_events_mod.lhe')
      OPEN(LUN,FILE='unweighted_events.lhe')
      WRITE(CHLUN,'(I3)') LUN
      CALL PYGIVE('MSTP(161)='//CHLUN)
      CALL PYGIVE('MSTP(162)='//CHLUN)
C...2) Initialize Pythia for user process  
ccc Life time < 10^-10 s
      MDCY(PYCOMP(15),1) = 1
      MDCY(PYCOMP(310),1) = 1
      MDCY(PYCOMP(111),1) = 1
      MDCY(PYCOMP(311),1) = 1
      MDCY(PYCOMP(3122),1) = 1
ccc Life time > 10^-10 s
      MDCY(PYCOMP(13),1) = 0
      MDCY(PYCOMP(211),1) = 0
      MDCY(PYCOMP(130),1) = 0
      MDCY(PYCOMP(321),1) = 0
      MDCY(PYCOMP(2112),1) = 0
c      MDCY(PYCOMP(221),1)=1
c      MDCY(PYCOMP(421),1)=1
c      MDCY(PYCOMP(511),1)=1
c      MDCY(PYCOMP(521),1)=1
c      MSTJ(104)=6 ! include top-pair production
c      MSTJ(24)=0 ! 
c      MSTP(122)=0 ! suppress some output
c      MSTU(21)=1 
c      MSTU(22)=0
c      MSTU(25)=0 
      CALL PYINRE
      CALL PYINIT('USER',' ',' ',0D0)
      
C------------------------- Initialize Parameters ------------------------

      me = 0.000511
      mp = 0.938272
      mn = 0.939565
      mpi = 0.13957018
      mkpm = 0.493677
      mkL = 0.497614

C------------------------- Initialize Histogram ------------------------
      xmin = 0.01
      xmax = 1000000
      do i = 1, nbins
         x(i) = 10**( dlog10(xmin) +( dlog10(xmax) -dlog10(xmin) )/nbins
     &        *i)
      enddo

      do i = 1,nbins
         nnGam(i) = 0
         nnEle(i) = 0
         nnNeue(i) = 0
         nnNeumu(i) = 0
         nnNeutau(i) = 0
         nnP(i) = 0
         nnaP(i) = 0
         nnN(i) = 0
         nnaN(i) = 0
         nnpip(i) = 0
         nnpim(i) = 0
         nnKp(i) = 0
         nnKm(i) = 0
         nnKL(i) = 0
      enddo
      totGam = 0
      totEle = 0
      totNeue = 0
      totNeumu = 0
      totNeutau = 0
      totP = 0
      totaP = 0
      totN = 0
      totaN = 0
      totpip = 0
      totpim = 0
      totKp = 0
      totKm = 0
      totKL = 0
      totOther = 0

C------------------------- GENERATE EVENTS ---------------------------

C...Initial values for number of events and cumulative charged multiplicity
      IEV=0
c      DNSUM=0D0
c      DN2SUM=0D0

C...Get next event from file and process it
 100  CALL PYEVNT

C...If event generation failed, quit loop
      IF(MSTI(51).EQ.1) THEN
        GOTO 999
      ENDIF

C...Else count up number of generated events
      IEV=IEV+1

C...Print first few events, both LHEF input and Pythia output.
      IF(IEV.LE.3) THEN
        CALL PYLIST(7)
        CALL PYLIST(2)
      ENDIF

C.../PYJETS/ now contains a fully generated event.
C...Insert user analysis here (or save event to output) 
C...(example: count charged multiplicity)
      CALL PYEDIT(1) ! decayed particles are removed
c      write(*,*) IEV,N
      do i = 1,N
         if (K(i,2).eq.22) then
            Ekin = P(i,4)
            totGam = totGam +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnGam(j) = nnGam(j) +1
                  exit
               endif
            enddo
         elseif (abs(K(i,2)).eq.11) then
            Ekin = P(i,4) -me
            totEle = totEle +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnEle(j) = nnEle(j) +1
                  exit
               endif
            enddo
         elseif (abs(K(i,2)).eq.12) then
            Ekin = P(i,4)
            totNeue = totNeue +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnNeue(j) = nnNeue(j) +1
                  exit
               endif
            enddo
         elseif (abs(K(i,2)).eq.14) then
            Ekin = P(i,4)
            totNeumu = totNeumu +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnNeumu(j) = nnNeumu(j) +1
                  exit
               endif
            enddo
         elseif (abs(K(i,2)).eq.16) then
            Ekin = P(i,4)
            totNeutau = totNeutau +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnNeutau(j) = nnNeutau(j) +1
                  exit
               endif
            enddo
         elseif (K(i,2).eq.2212) then
            Ekin = P(i,4) -mp
            totP = totP +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnP(j) = nnP(j) +1
                  exit
               endif
            enddo
         elseif (K(i,2).eq.-2212) then
            Ekin = P(i,4) -mp
            totaP = totaP +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnaP(j) = nnaP(j) +1
                  exit
               endif
            enddo
         elseif (K(i,2).eq.2112) then
            Ekin = P(i,4) -mn
            totN = totN +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnN(j) = nnN(j) +1
                  exit
               endif
            enddo
         elseif (K(i,2).eq.-2112) then
            Ekin = P(i,4) -mn
            totaN = totaN +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnaN(j) = nnaN(j) +1
                  exit
               endif
            enddo
         elseif (K(i,2).eq.211) then
            Ekin = P(i,4) -mpi
            totpip = totpip +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnpip(j) = nnpip(j) +1
                  exit
               endif
            enddo
         elseif (K(i,2).eq.-211) then
            Ekin = P(i,4) -mpi
            totpim = totpim +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnpim(j) = nnpim(j) +1
                  exit
               endif
            enddo
         elseif (K(i,2).eq.321) then
            Ekin = P(i,4) -mkpm
            totKp = totKp +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnKp(j) = nnKp(j) +1
                  exit
               endif
            enddo
         elseif (K(i,2).eq.-321) then
            Ekin = P(i,4) -mkpm
            totKm = totKm +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnKm(j) = nnKm(j) +1
                  exit
               endif
            enddo
         elseif (K(i,2).eq.130) then
            Ekin = P(i,4) -mkL
            totKL = totKL +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnKL(j) = nnKL(j) +1
                  exit
               endif
            enddo
         else
            totOther = totOther +1
         endif
      enddo
                              
C...Loop back to look for next event
      GOTO 100

C...Jump point when end-of-file reached (or other problem encountered)
C...Print final statistics.
 999  CALL PYSTAT(1)

c      DNAVG=DNSUM/IEV
c      DNRMS=SQRT(DN2SUM/IEV)
c      SIGMA2=MAX(0D0,DNRMS**2-DNAVG**2)

      open(1,file="np_sptrm_pythia6.dat",status="replace")
      do i = 1,nbins
         write(1,*) 2*mDM,x(i), nnN(i)/dble(IEV),nnP(i)/dble(IEV)
     &        ,nnpip(i)/dble(IEV),nnpim(i)/dble(IEV),nnKp(i)/dble(IEV)
     &        ,nnKm(i)/dble(IEV),nnKL(i)/dble(IEV),nnaN(i)/dble(IEV)
     &        ,nnaP(i)/dble(IEV)
      enddo
      close(1)

      open(1,file="nini_pythia6.dat",status="replace")
      write(1,*) 2*mDM,totN/dble(IEV),totP/dble(IEV)
     &     ,totpip/dble(IEV),totpim/dble(IEV),totKp/dble(IEV)
     &     ,totKm/dble(IEV),totKL/dble(IEV),totaN/dble(IEV)
     &     ,totaP/dble(IEV)      
      close(1)

      write(*,*) 2*mDM,totN/dble(IEV),totP/dble(IEV)
     &     ,totpip/dble(IEV),totpim/dble(IEV),totKp/dble(IEV)
     &     ,totKm/dble(IEV),totKL/dble(IEV),totaN/dble(IEV)
     &     ,totaP/dble(IEV)

      END
