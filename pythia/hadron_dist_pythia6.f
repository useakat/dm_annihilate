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
      COMMON/PYJETS/N,NPAD,K(30000,5),P(30000,5),V(30000,5)
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

      integer i,j,i1
      integer nbins,maxnbins,PYCOMP,nerrors
      parameter (maxnbins=1000)
      integer nnGam(maxnbins),nnEle(maxnbins),nnNeue(maxnbins)
     &     ,nnNeumu(maxnbins),nnNeutau(maxnbins),nnP(maxnbins)
      integer nnKp(maxnbins),nnKm(maxnbins),nnKL(maxnbins)
      integer nnaN(maxnbins),nnN(maxnbins),nnaP(maxnbins)
      integer nnpip(maxnbins),nnpim(maxnbins),nnK0(maxnbins)
      integer totGam,totEle,totNeue,totNeumu,totNeutau,totP,totN
      integer totpip,totpim,totKp,totKm,totKL,totaN,totaP,totK0
      integer totOther
      integer NEV,IDC,idecay,ipyevnw
      real*8 xmin,xmax,x(maxnbins),dx
      real*8 me,mp,mn,mpi,mkpm,mkL,nevents,Nch,mDM,ECM
      character*10 cNEV,cmDM,cdecay
      real*8 Ekin,Ekin_neu,Ekin_other

      integer MSTU,MSTJ
      double precision PARU,PARJ
      integer MDCY,MDME,KFDP
      double precision BRAT
C-------------------------- PYTHIA SETUP -----------------------------

C...Main parameters of run: c.m. energy and number of events.
      call getarg(1,cmDM)
      call getarg(2,cNEV)
      call getarg(3,cdecay)
      read (cdecay,*) idecay
      if (idecay.eq.0) then
         read (cmDM,*) ECM
      else
         read (cmDM,*) mDM
         ECM = 2*mDM
         open(1,file='mass.dat',status='replace')
         write(1,'(e11.4)') mDM
         close(1)
      endif
      read (cNEV,*) NEV
      
      ipyevnw = 1
      MSTU(4) = 30000 ! change the array size of K,P,V
      MSTU(5) = 30000 ! change the array size of K,P,V
c      MSTU(21)=1 ! check on possible errors during program execution
      MSTU(22) = 10000 ! The maxium number of errors to be accepted
c      MSTU(25)=0 ! printing of warning messages
      PARP(2) = 2d0 ! lowest c.m. energy for the event as a whole that the program will accept to simulate
                    ! e+ e- collision can be simulated for PARP(2) > 4 GeV with some confidence (P.217 in Pythia6 manual). 

      if (ipyevnw.eq.1) then ! PYEVNW tune according to https://indico.cern.ch/event/11784/session/8/contribution/97/material/slides/0.ppt
         MSTP(5) = 350 ! Perugia 11
      endif
c$$$         MSTP(70) = 2           ! smooth ISR
c$$$         MSTP(72) = 0           ! high FSR
c$$$         MSTP(81) = 20          ! no MPI, PYEVNW
c$$$         MSTP(82) = 4           ! structure of MPI
c$$$         MSTP(84) = 1
c$$$         MSTP(85) = 1
c$$$         MSTP(86) = 2
c$$$         MSTP(87) = 4
c$$$         MSTP(88) = 0
c$$$         MSTP(89) = 1
c$$$         MSTP(90) = 1
c$$$         MSTP(95) = 1           ! reconnect
c$$$         
c$$$c      PARP(63) = 1.8d0
c$$$         PARP(78) = 0.2d0
c$$$c     PARP(79) = 2d0
c$$$         PARP(80) = 0.01d0
c$$$         PARP(82) = 1.9d0       ! pT0
c$$$         PARP(83) = 0.3d0
c$$$         PARP(84) = 0.5d0
c$$$         PARP(89) = 1800d0
c$$$         PARP(90) = 0.22d0
c$$$         
c$$$         PARJ(81) = 0.14d0
c$$$      endif

      MSEL=0     ! process are selected by process number
      MSTP(42) = 0              ! 0:resonance is on the mass shell 1:Breigt Wigner
      if (idecay.eq.0) then ! for checking 
         MSUB(1) = 1            ! gam/Z0 production
         MSTP(43) = 2           ! 1:only gam*  2:only Z0* 3:both
         DO IDC = MDCY(23,2),MDCY(23,2)+MDCY(23,3)-1
            IF(IABS(KFDP(IDC,1)).gt.5) MDME(IDC,1)=MIN(0,MDME(IDC,1))
         enddo
c         DO IDC = MDCY(22,2),MDCY(22,2)+MDCY(22,3)-1
c            IF(IABS(KFDP(IDC,1)).gt.5) MDME(IDC,1)=MIN(0,MDME(IDC,1))
c         enddo
      elseif (idecay.eq.24) then
         MSUB(25) = 1           ! W+W-
      elseif (idecay.eq.23) then
         MSUB(22) = 1           ! ZZ
         MSTP(43) = 2           ! 1:only gam*  2:only Z0*
      elseif (idecay.eq.25) then
         MSUB(27) = 1           ! ff > hh
c$$$         DO IDC = MDCY(35,2),MDCY(35,2)+MDCY(35,3)-1
c$$$            IF(IABS(KFDP(IDC,1)).ne.idecay) 
c$$$     &           MDME(IDC,1)=MIN(0,MDME(IDC,1))
c$$$         enddo
      else
         MSUB(1) = 1            ! gam/Z0 production
         MSTP(43) = 2           ! 1:only gam*  2:only Z0*
         DO IDC = MDCY(23,2),MDCY(23,2)+MDCY(23,3)-1
            IF(IABS(KFDP(IDC,1)).ne.idecay) 
     &           MDME(IDC,1)=MIN(0,MDME(IDC,1))
         enddo
      endif
      
c$$$
c$$$C...1) Open LHEF file on unit LUN, and tell Pythia where to find it.
c$$$      LUN=88
c$$$      OPEN(LUN,FILE='unweighted_events_mod.lhe')
c$$$      OPEN(LUN,FILE='unweighted_events.lhe')
c$$$      WRITE(CHLUN,'(I3)') LUN
c$$$      CALL PYGIVE('MSTP(161)='//CHLUN)
c$$$      CALL PYGIVE('MSTP(162)='//CHLUN)

C...2) Initialize Pythia for user process  
      if (idecay.eq.0) then
         MDCY(PYCOMP(13),1) = 0  ! mu
         MDCY(PYCOMP(15),1) = 1  ! tau
         MDCY(PYCOMP(130),1) = 0 ! KL
         MDCY(PYCOMP(310),1) = 1 ! Ks
         MDCY(PYCOMP(111),1) = 1 ! pi0
         MDCY(PYCOMP(211),1) = 0 ! pi+-
         MDCY(PYCOMP(311),1) = 0 ! K0
         MDCY(PYCOMP(321),1) = 0 ! K+-
         MDCY(PYCOMP(2112),1) = 0 ! n
         MDCY(PYCOMP(3122),1) = 1 ! Lambda0
      else      
ccc   Life time < 10^-10 s
         MDCY(PYCOMP(15),1) = 1
         MDCY(PYCOMP(310),1) = 1
         MDCY(PYCOMP(111),1) = 1
         MDCY(PYCOMP(311),1) = 1
         MDCY(PYCOMP(3122),1) = 1
ccc   Life time > 10^-10 s
         MDCY(PYCOMP(13),1) = 0
         MDCY(PYCOMP(211),1) = 0
         MDCY(PYCOMP(130),1) = 0
         MDCY(PYCOMP(321),1) = 0
         MDCY(PYCOMP(2112),1) = 0
      endif
c      do i = 1,500
c         MDCY(i,1) = 0
c      enddo

c      MSTJ(104)=6 ! allow tt~ production
c      MSTJ(24)=0  ! discrete mass value
      MSTP(11)=0  ! no ISR photon radiation
      CALL PYINRE
      CALL PYINIT('CMS','e+','e-',ECM)
c      CALL PYINIT('USER',' ',' ',0D0)

C...Check that Z0 decay channels set correctly.
      CALL PYSTAT(2)
      
C------------------------- Initialize Parameters ------------------------
      me = 0.000511
      mp = 0.938272
      mn = 0.939565
      mpi = 0.13957018;
      mkpm = 0.493677;
      mkL = 0.497614; ! K0 mass

C------------------------- Initialize Histogram ------------------------
      if (idecay.eq.0) then
         xmin = 0.01
         xmax = 100
         nbins = 72
      else
         xmin = 0.01
         xmax = 1000000
         nbins = 800
      endif
      do i = 1, nbins+1
         x(i) = 10**( dlog10(xmin) +( dlog10(xmax) -dlog10(xmin) )/nbins
     &        *(i-1))
      enddo

      do i = 1,nbins+1
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
         nnK0(i) = 0
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
      totK0 = 0
      totOther = 0

C------------------------- GENERATE EVENTS ---------------------------

C...Initial values for number of events and cumulative charged multiplicity
      nevents = 0
      nerrors = 0
C...  Get next event from file and process it
      DO IEV = 1,NEV
         if (ipyevnw.eq.1) then
            CALL PYEVNW
         else
            CALL PYEVNT
         endif
C...  If event generation failed, quit loop
         IF(MSTI(51).EQ.1) THEN
            GOTO 999
         ENDIF
         IF(MSTU(23).gt.nerrors) then
            nerrors = MSTU(23)
c     GOTO 999
            cycle
c            continue
         endif
         
C...  Print first few events, both LHEF input and Pythia output.
         IF(IEV.LE.10) THEN
            CALL PYLIST(7)
            CALL PYLIST(2)
         ENDIF
         
c     CALL PYEDIT(1) ! decayed particles are removed
         CALL PYEDIT(1)         ! neutral particles are removed
         do i = 1,N
            if (K(i,2).eq.22) then
               Ekin = P(i,4)
               totGam = totGam +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnGam(j) = nnGam(j) +1
                     exit
                  endif
               enddo
            elseif (abs(K(i,2)).eq.11) then
               Ekin = P(i,4) -me
               totEle = totEle +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnEle(j) = nnEle(j) +1
                     exit
                  endif
               enddo
            elseif (abs(K(i,2)).eq.12) then
               Ekin = P(i,4)
               totNeue = totNeue +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then                 
                     nnNeue(j) = nnNeue(j) +1
                     exit
                  endif
               enddo
            elseif (abs(K(i,2)).eq.14) then
               Ekin = P(i,4)
               totNeumu = totNeumu +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnNeumu(j) = nnNeumu(j) +1
                     exit
                  endif
               enddo
            elseif (abs(K(i,2)).eq.16) then
               Ekin = P(i,4)
               totNeutau = totNeutau +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnNeutau(j) = nnNeutau(j) +1
                     exit
                  endif
               enddo
            elseif (K(i,2).eq.2212) then
               Ekin = P(i,4) -mp
               totP = totP +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnP(j) = nnP(j) +1
                     exit
                  endif
               enddo
            elseif (K(i,2).eq.-2212) then
               Ekin = P(i,4) -mp
               totaP = totaP +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnaP(j) = nnaP(j) +1
                     exit
                  endif
               enddo
            elseif (K(i,2).eq.2112) then
               Ekin = P(i,4) -mn
               totN = totN +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnN(j) = nnN(j) +1
                     exit
                  endif
               enddo
            elseif (K(i,2).eq.-2112) then
               Ekin = P(i,4) -mn
               totaN = totaN +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnaN(j) = nnaN(j) +1
                     exit
                  endif
               enddo
            elseif (K(i,2).eq.211) then
               Ekin = P(i,4) -mpi
               totpip = totpip +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnpip(j) = nnpip(j) +1
                     exit
                  endif
               enddo
            elseif (K(i,2).eq.-211) then
               Ekin = P(i,4) -mpi
               totpim = totpim +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnpim(j) = nnpim(j) +1
                     exit
                  endif
               enddo
            elseif (K(i,2).eq.321) then
               Ekin = P(i,4) -mkpm
               totKp = totKp +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnKp(j) = nnKp(j) +1
                     exit
                  endif
               enddo
            elseif (K(i,2).eq.-321) then
               Ekin = P(i,4) -mkpm
               totKm = totKm +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnKm(j) = nnKm(j) +1
                     exit
                  endif
               enddo
            elseif (K(i,2).eq.130) then
               Ekin = P(i,4) -mkL
               totKL = totKL +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnKL(j) = nnKL(j) +1
                     exit
                  endif
               enddo
            elseif (abs(K(i,2)).eq.311) then
               Ekin = P(i,4) -mkL
               totK0 = totK0 +1
               do j = 1,nbins
                  if ((x(j).lt.Ekin).and.(Ekin.lt.x(j+1))) then
                     nnK0(j) = nnK0(j) +1
                     exit
                  endif
               enddo
            else
               totOther = totOther +1
            endif
         enddo
         nevents = nevents +1         
C...  Loop back to look for next event
      enddo
         
C...  Print final statistics.
 999  CALL PYSTAT(1)
         
      open(1,file="np_sptrm.dat",status="replace")
      do i = 1,nbins
         write(1,'(11e11.4)') ECM,x(i),nnN(i)/dble(nevents)
     &        ,nnP(i)/dble(nevents)
     &        ,nnpip(i)/dble(nevents),nnpim(i)/dble(nevents)
     &        ,nnKp(i)/dble(nevents),nnKm(i)/dble(nevents)
     &        ,(nnKL(i)+nnK0(i)/2d0)/dble(nevents),nnaN(i)/dble(nevents)
     &        ,nnaP(i)/dble(nevents)
      enddo
      close(1)

      open(1,file="np_sptrm_norm.dat",status="replace")
      do i = 1,nbins+1
         dx = x(i+1) -x(i)
         write(1,'(11e11.4)') ECM,x(i),nnN(i)/dble(nevents)/dx
     &        ,nnP(i)/dble(nevents)/dx
     &        ,nnpip(i)/dble(nevents)/dx,nnpim(i)/dble(nevents)/dx
     &        ,nnKp(i)/dble(nevents)/dx,nnKm(i)/dble(nevents)/dx
     &        ,(nnKL(i)+nnK0(i)/2d0)/dble(nevents)/dx
     &        ,nnaN(i)/dble(nevents)/dx,nnaP(i)/dble(nevents)/dx
      enddo
      close(1)

      open(1,file="Edist.dat",status="replace")
      write(1,'(3e11.4)') (totpip+totpim)/dble(nevents)
     &     ,(totN+totaN)/dble(nevents),(totP+totaP)/dble(nevents)
      do i = 1,nbins+1
         write(1,'(4e11.4)') x(i),(nnpip(i)+nnpim(i))/dble(nevents)
     &        ,(nnN(i)+nnaN(i))/dble(nevents)
     &        ,(nnP(i)+nnaP(i))/dble(nevents)
      enddo
      close(1)

      open(1,file="nini.dat",status="replace")
      write(1,'(10e11.4)') ECM,totN/dble(nevents),totP/dble(nevents)
     &     ,totpip/dble(nevents),totpim/dble(nevents)
     &     ,totKp/dble(nevents),totKm/dble(nevents)
     &     ,(totKL+totK0/2d0)/dble(nevents)
     &     ,totaN/dble(nevents),totaP/dble(nevents)
      close(1)

      write(*,*) "counted evnets:", nevents
      write(*,*) "N particles:", totGam,totEle,totP,totN,totOther
      write(*,*) "avg particles:", totGam/dble(nevents)
     &     ,totEle/dble(nevents),totP/dble(nevents),totN/dble(nevents)
     &     ,totOther/dble(nevents)
      write(*,*)

cccccccccccccccccc            ccccccccccccccccccccccccc 
cccccccccccccccccc  Evis data ccccccccccccccccccccccccc
      if (idecay.ne.0) then
ccc   Life time > 10^-10 s
         MDCY(PYCOMP(13),1) = 1
         MDCY(PYCOMP(211),1) = 1
         MDCY(PYCOMP(130),1) = 1
         MDCY(PYCOMP(311),1) = 1
         MDCY(PYCOMP(321),1) = 1
         MDCY(PYCOMP(2112),1) = 0
      endif         

      CALL PYINIT('CMS','e+','e-',ECM)

C...Initial values for number of events and cumulative charged multiplicity
      nevents = 0
      Ekin = 0
      Ekin_neu = 0
      Ekin_other = 0
      nerrors = MSTU(23)
C...Get next event from file and process it
      DO IEV = 1,NEV
         if (ipyevnw.eq.1) then
            CALL PYEVNW
         else
            CALL PYEVNT
         endif
C...If event generation failed, quit loop
         IF(MSTI(51).EQ.1) THEN
            GOTO 1999
         ENDIF
         IF(MSTU(23).gt.nerrors) then
            nerrors = MSTU(23)
c     GOTO 1999
            cycle
c            continue
         endif
         
C...  Print first few events, both LHEF input and Pythia output.
         IF(IEV.LE.10) THEN
            CALL PYLIST(7)
            CALL PYLIST(2)
         ENDIF
         
c     CALL PYEDIT(1) ! decayed particles are removed
         CALL PYEDIT(1)         ! decayed particles are removed
         do i = 1,N
            if (abs(K(i,2)).eq.12) then
               Ekin_neu = Ekin_neu +P(i,4)
            elseif (abs(K(i,2)).eq.14) then
               Ekin_neu = Ekin_neu +P(i,4)
            elseif (abs(K(i,2)).eq.16) then
               Ekin_neu = Ekin_neu +P(i,4)
            elseif (abs(K(i,2)).eq.11) then
               Ekin = Ekin +P(i,4) -me
            elseif (K(i,2).eq.22) then
               Ekin = Ekin +P(i,4)
            elseif (abs(K(i,2)).eq.2112) then
               Ekin = Ekin +P(i,4) -mn
            elseif (abs(K(i,2)).eq.2212) then
               Ekin = Ekin +P(i,4) -mp
            else
               Ekin_other = Ekin_Other +P(i,4)
            endif
         enddo
         nevents = nevents +1           
C...  Loop back to look for next event
      enddo

 1999 CALL PYSTAT(1)

      open(1,file="Evis.dat",status="replace")
      write(1,'(3e11.4)') ECM,Ekin/dble(nevents)
     &     ,Ekin/ECM/dble(nevents)
      close(1)
      write(*,'(a5,e11.4,a15,e11.4,a19,e11.4)') "ECM ="
     &     ,ECM," Ekin/nevents ="
     &     ,Ekin/dble(nevents)," Ekin/ECM/nevents ="
     &     ,Ekin/ECM/dble(nevents)
      write(*,*)

      open(1,file="Evis_tot.dat",status="replace")
      write(1,'(e11.4)') Ekin/dble(nevents)
      close(1)

      END
