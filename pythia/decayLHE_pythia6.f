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
      parameter (nbins=60)
      integer nnGam(nbins),nnEle(nbins),nnNeue(nbins),nnNeumu(nbins)
     &     ,nnNeutau(nbins),nnP(nbins),nnN(nbins)
      integer totGam,totEle,totNeue,totNeumu,totNeutau,totP,totN
      integer totOther
      real*8 xmin,xmax,x(nbins),dx
      real*8 me,mp,mn
      real*8 nnevents

      COMMON/PYDAT1/MSTU(200),PARU(200),MSTJ(200),PARJ(200)
      integer MSTU,MSTJ
      double precision PARU,PARJ
      COMMON/PYDAT3/MDCY(500,3),MDME(8000,2),BRAT(8000),KFDP(8000,5)
      integer MDCY,MDME,KFDP
      double precision BRAT
C-------------------------- PYTHIA SETUP -----------------------------

C...1) Open LHEF file on unit LUN, and tell Pythia where to find it.
      LUN=88
c      OPEN(LUN,FILE='unweighted_events_mod.lhe')
      OPEN(LUN,FILE='unweighted_events_mod_900k.lhe')
      WRITE(CHLUN,'(I3)') LUN
      CALL PYGIVE('MSTP(161)='//CHLUN)
      CALL PYGIVE('MSTP(162)='//CHLUN)
C...2) Initialize Pythia for user process  
      MDCY(PYCOMP(13),1) = 1
      MDCY(PYCOMP(15),1) = 1
      MDCY(PYCOMP(130),1) = 1
      MDCY(PYCOMP(310),1) = 1
      MDCY(PYCOMP(111),1) = 1
      MDCY(PYCOMP(211),1) = 1
      MDCY(PYCOMP(311),1) = 1
      MDCY(PYCOMP(321),1) = 1
c      MDCY(PYCOMP(2112),1) = 0
      MDCY(PYCOMP(221),1)=1
      MDCY(PYCOMP(421),1)=1
      MDCY(PYCOMP(511),1)=1
      MDCY(PYCOMP(521),1)=1
      MSTJ(104)=6 
      MSTJ(24)=0
      MSTP(122)=0 
      MSTU(21)=1 
      MSTU(22)=0
      MSTU(25)=0 
      CALL PYINRE
      CALL PYINIT('USER',' ',' ',0D0)
      
C------------------------- Initialize Parameters ------------------------

      me = 0.000511
      mp = 0.938272
      mn = 0.939565

C------------------------- Initialize Histogram ------------------------
      xmin = 0.0001
      xmax = 100
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
         nnN(i) = 0
      enddo
      totGam = 0
      totEle = 0
      totNeue = 0
      totNeumu = 0
      totNeutau = 0
      totP = 0
      totN = 0

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
         elseif (abs(K(i,2)).eq.2212) then
            Ekin = P(i,4) -mp
            totP = totP +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnP(j) = nnP(j) +1
                  exit
               endif
            enddo
         elseif (abs(K(i,2)).eq.2112) then
            Ekin = P(i,4) -mn
            totN = totN +1
            do j = 1,nbins
               if ((x(j).gt.Ekin).and.(Ekin.lt.x(j+1))) then
                  nnN(j) = nnN(j) +1
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

      open(1,file="dist_takaesu_pythia6.dat",status="replace")
      do i = 1,nbins
         write(1,*) x(i), nnEle(i)/dble(IEV),nnGam(i)/dble(IEV)
     &        ,nnP(i)/dble(IEV),nnN(i)/dble(IEV),nnNeue(i)/dble(IEV)
     &        ,nnNeumu(i)/dble(IEV),nnNeutau(i)/dble(IEV)
      enddo
      close(1)

      write(*,*) "counted evnets:", IEV
      write(*,*) "N particles:", totGam,totEle,totP,totN,totOther

      END
