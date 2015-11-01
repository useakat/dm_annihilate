c the following routines were written by ramona on 19/6/2013
c ===================================================================== c
c               light stop 4-body decays				c
c ===================================================================== c
      subroutine SD_lightstop4bod(totwidth, widthbtau,widthbbjet,
     .widthbjets,widthbmu,widthbelec,widthjettau,widthjetb,
     .widthjetjet,widthjetmu,widthjetelec)
      implicit none
      integer ifavvio, i, j, k, s, iiter, kk, jj, ifer, ifer2, ifer3
      double precision sdgf, amz, amw, pi, g2, mf
      double precision amsupq(6), amsdownq(6), amslepton(6)
      double precision sw, cw, alp_mssm, tanbeta, amt, amb, amtau
      double precision uu(2,2), vv(2,2), zz(4,4), zp(4,4)
      double precision vckm(3,3), msq2(3,3), msu2(3,3), 
     .msd2(3,3), td(3,3), tu(3,3),
     . usqmix(6,6), dsqmix(6,6)
      double precision sinbeta,alsew,g2ew,g1ew, time
      double precision sigma, DcharL(3,6,2), EcharR(3,6,2), 
     . cchaneutL(4,2), cchaneutR(4,2) 
      double precision sigbtau
      double precision mfin(100),pfout(4,100), etot,  wt, ampwchabtau
      double precision amneut(4),xmneut(4),amchar(2),xmchar(2)
      double precision p1p2, p1p3, p1p4, p2p3, p2p4, p3p4,pI2, 
     . p1pi,p2pi, p3pi, p4pi, m1, m3, m4, g1L1, g1L2, g1R1, g1r2,g3L, 
     . g2l1, g2r1, g2l2, g2r2, p1(4), p2(4), p3(4), p4(4), pInt(4),
     .SD_scal, pfer(4), p1pfer, p2pfer, p3pfer, p4pfer, pfer2, pipfer
      double precision ampwchabmu, sigbmu, ampwchar, zwi, 
     . ampbsquarkjets, ampjetssquarkb, ampjetsquarkjets
      double precision widthbtau, widthbmu, widthbjets,
     . sigbjets, ampwchabjets, widthbbjet, widthjetb, 
     . widthjetjet, widthjetmu, widthjettau, widthjetelec
       double precision ampSMfer, ampcharsup
      double precision sigbbjet, sigjetb, sigjetjet,
     . sigjetmu, sigjettau, sigjetelec 
      double precision ampwchabbjet, ampwchajetb, ampwchajetjet, 
     . ampwchajetmu, ampwchajettau, totwidth, ampcharsupsdownint
      double precision mfd(3), mfu(3), selmix(6,6), mfe(3)
      double precision cchaneutslep(3,2,6), dchalepsneut(3,2), 
     . echalepsneut(3,2), ampslepbtau, gneutneut(3,4), 
     . gneutelecl(3,4,6), gneutelecr(3,4,6), g3l1, g3l2, g3r1, g3r2
      double precision ccharl(3,6,2), fcharr(3,6,2), gneutul(3,4,6), 
     . gneutur(3,4,6), gneutdl(3,4,6), gneutdr(3,4,6), QchncR(4,2), 
     . QchncL(4,2)
      double precision ampcharslepton, ampintslepsneut, ampintwslep
      double precision gXb2R, gXc2L, gXc2R,gbd3R, gcd3L, gxb2l,
     . gcd3R, amsneutrino(3), ampinterferenzwslep, ampintwsneut,
     . ampcharsneutrino, ampinterferenzwselec, ampinterferenzwsmu,
     . ampslepbelec, ampslepbmu, sigbelec, widthbelec,
     . ampinterferenzwjetslep, ampslepjettau, ampinterferenzwjetelec,
     . ampinterferenzwjetmu, ampslepjetelec, ampslepjetmu, wwidth,
     . ampbsquarkb, ampSMfermion, ampintWSMfersneut, ampintWSMferslep,
     . ampintWSMferWchar, ampSMfermionjet, ampSMfermionel,
     . ampsquarkpropW, ampvecscalar,ampchahi,amphiggssquark, 
     . ampinthiggssquarkWsquark,ampSMfermionmu, ampintcharsneutWsquark,
     . ampintcharslepWsquark, ampintWcharWsquark, ampintWferWsquark, 
     . ampsquarkpropWe, ampsquarkpropWmu, ampsquarkpropWjet,
     . ampintchaHisquarkWchar, ampintchaHisquarkWSMfer, ampchaHiSMfer,
     . ampintchaHisquarkslepton, ampintchaHisquarksneutrino, 
     . ampchaHichar, ampintchaHichargslepton, ampintchahichargsneutrino,
     . ampintchaHiWSMfer, ampintchaHiWchar, ampintchaHichargchaHisquark,
     . ampintchaHichargWsquark, ampintchaHiSmferslepton, 
     . ampintchaHiSmfersneutrino, ampintchaHiSmferWcharg, 
     . ampintchaHiSmferWsquark, ampintchaHiSmferWSMfer, ampchahib, 
     . ampintchaHiSmferchaHicharg, ampintchaHiSmferchaHisquark,
     . ampintchaHisquarksup
      double precision sdatop,sdabot,sdatau,sdmu
      double precision gsqsqchaHi(6),zwi2,gsqsqW(6)
      double precision ama,aml,amh,amch,amar
      double precision test, test2


      COMMON/SD_param/sdgf,amz,amw,pi,g2
      COMMON/flavviolation/vckm, msq2, msd2, msu2, td, 
     .tu, usqmix, ifavvio, dsqmix
      COMMON/msfermion/ amsupq, amsdownq, amslepton, amsneutrino
      COMMON/SD_weinberg/sw,cw
      COMMON/SD_mixang/alp_mssm,tanbeta
      COMMON/SD_fermion/amt,amb,amtau
      COMMON/SD_mixmat/uu,vv,zz,zp
      COMMON/SD_coupewsb/alsew,g2ew,g1ew
      COMMON/SD_massgino/amneut,xmneut,amchar,xmchar
      common/SD_selectron/selmix
      COMMON/SD_break/sdatop,sdabot,sdatau,sdmu
      COMMON/SD_hmass/ama,aml,amh,amch,amar

	PI = 4*DATAN(1D0)
c------ initialize random number generator
      call SD_rstart(12,34,56,78)


c---- set fermion masses
      mfd(1)=0d0
      mfd(2)=0d0
      mfd(3)=amb
      mfu(1)=0d0
      mfu(2)=0d0
      mfu(3)=amt
      mfe(1)=0d0
      mfe(2)=0d0
      mfe(3)=amtau
c---- W boson width      
      wwidth=2.085d0



c---- couplings, notation see SPARTICLES book
c---- chargino-squark quark coupling
      do i=1,3
      do k=1,2
      do s=1,6
! Attention in sparticles book on which one is with left handed projector and which one is with right handed
      DcharL(i,s,k)=g2ew*vv(k,2)/Sqrt(2d0)/amw/dsin(datan(tanbeta))
     .		    *mfu(3)*Vckm(3,i)*USQMix(s,6)-g2ew*vv(k,1)*
     .		    (Usqmix(s,1)*Vckm(1,i)+Usqmix(s,2)*Vckm(2,i)
     .		    +Usqmix(s,3)*Vckm(3,i))
      EcharR(i,s,k)=g2ew*uu(k,2)/dsqrt(2d0)/amw
     .		    /dcos(datan(tanbeta))*mfd(i)*usqmix(s,3)
     .		    *Vckm(3,i)
      FcharR(i,s,k)=g2ew*vv(k,2)/dsqrt(2d0)/amw
     .		    /dsin(datan(tanbeta))*mfu(i)*dsqmix(s,3)
     .		    *Vckm(i,3)
       CcharL(i,s,k)=g2ew*uu(k,2)/Sqrt(2d0)/amw/dcos(datan(tanbeta))
     .		    *mfd(3)*Vckm(i,3)*DSQMix(s,6)-g2ew*uu(k,1)*
     .		    (dsqmix(s,1)*Vckm(i,1)+dsqmix(s,2)*Vckm(i,2)
     .		    +dsqmix(s,3)*Vckm(i,3))
      
      end do
      end do
      end do
      
      
      do i=1,3
      do k=1,4
      do s=1,6
      Gneutul(i,k,s)=-Sqrt(2d0)*(1d0/2d0*ZZ(k,2)*g2ew
     .               +1d0/6d0*g1ew*ZZ(k,1))*Usqmix(s,i)-g2ew*mfu(i)
     .	             /(dsqrt(2d0)*amw*Sin(datan(tanbeta)))
     .	             *ZZ(k,4)*Usqmix(s,i+3)
      Gneutur(i,k,s)=2d0*dsqrt(2d0)/3d0*g1ew*ZZ(k,1)
     .		     *Usqmix(s,i+3)-g2ew*mfu(i)
     . 	      	     /(dsqrt(2d0)*amw*Sin(datan(tanbeta)))
     .		     *ZZ(k,4)*Usqmix(s,i)
      
      Gneutdl(i,k,s)=Sqrt(2d0)*(1d0/2d0*ZZ(k,2)*g2ew
     .      	      -1d0/6d0*g1ew*ZZ(k,1))*dsqmix(s,i)-g2ew*mfd(i)
     .	              /(dsqrt(2d0)*amw*Cos(datan(tanbeta)))
     .		      *ZZ(k,3)*dsqmix(s,i+3)
      Gneutdr(i,k,s)=-dsqrt(2d0)/3d0*g1ew*ZZ(k,1)
     .		     *dsqmix(s,i+3)-g2ew*mfd(i)
     . 	      	      /(dsqrt(2d0)*amw*Cos(datan(tanbeta)))
     .		      *ZZ(k,3)*dsqmix(s,i)
      end do
      end do
      end do


c------ chargino -lepton-slepton couplings
      do i=1,3
      do k=1,2
      do s=1,6
      cchaneutslep(i,k,s)=-g2ew*selmix(s,i)*UU(k,1)
     .			  +g2ew/dSqrt(2d0)/amw/dcos(datan(tanbeta))
     .			  *mfe(i)*selmix(s, i+3)*UU(k,2)


      end do
      dchalepsneut(i,k)=-g2ew*VV(k,1)
      echalepsneut(i,k)=g2ew*mfe(i)/dSqrt(2d0)
     .			  /amw/dcos(datan(tanbeta))*UU(k,2)
      end do
      end do

      do i=1,3
      do k=1,4
      do s=1,6
      
      gneutelecl(i,k,s)=g2ew/dsqrt(2d0)*ZZ(k,2)*selmix(s,i)
     .                  +g1ew/dsqrt(2d0)*ZZ(k,1)*selmix(s,i)
     .			-g2ew*mfe(i)/dsqrt(2d0)/amw/dcos(datan(tanbeta))
     .			*selmix(s, i+3)*ZZ(k,3)
      gneutelecr(i,k,s)=-dsqrt(2d0)*g1ew*ZZ(k,1)*selmix(s,i+3)
     .			  -g2ew*mfe(i)/dSqrt(2d0)
     .			  /amw/dcos(datan(tanbeta))
     .			 *ZZ(k,3)*selmix(s,i)
      end do
      gneutneut(i,k)=-g2ew/dsqrt(2d0)*ZZ(k,2)+g1ew/dsqrt(2d0)*ZZ(k,1)
      end do
      end do


c--- neutralino-chargino-W-coupling
      do k=1,2
      do i=1,4
      CchaneutL(i,k)=g2ew*(ZZ(i,2)*VV(k,1)
     .		     -1d0/dsqrt(2d0)*ZZ(i,4)*VV(k,2))
      CChaneutR(i,k)=g2ew*(ZZ(i,2)*UU(k,1)
     .		     +1d0/dsqrt(2d0)*ZZ(i,3)*UU(k,2))
      end do
      end do


c---- W-squark-squark coupling
      do k=1,6
      gsqsqW(k)=-g2ew/dsqrt(2d0)*((USQmix(1,1)*VCKM(1,1)
     . +USQMIX(1,2)*VCKM(2,1)+USQMIX(1,3)*VCKM(3,1))*Dsqmix(k,1)
     . +(USQmix(1,1)*VCKM(1,2)
     . +USQMIX(1,2)*VCKM(2,2)+USQMIX(1,3)*VCKM(3,2))*Dsqmix(k,2)
     . +(USQmix(1,1)*VCKM(1,3)
     . +USQMIX(1,2)*VCKM(2,3)+USQMIX(1,3)*VCKM(3,3))*Dsqmix(k,3))
      enddo


c---- squark squark charged Higgs coupling 
      
      do k=1,6
      zwi2=0d0
      zwi=0d0
      do i=1,3
      do j=1,3
      zwi2=0d0
      do s=1,3
      zwi2=zwi2-g2ew/dsqrt(2d0)/amw*
     .(tanbeta*td(s,j)*Vckm(i, s)*dsqmix(k,j+3)*usqmix(1,i)
     . +1d0/tanbeta*tu(s,j)*VCKM(s,i)*dsqmix(k,i)*usqmix(1, j+3))


      enddo
       zwi=zwi+zwi2-g2ew/dsqrt(2d0)/amw*
     .(sdmu*USQMIX(1,i)*DSQMIX(k,j+3)*VCKM(i,j)*mfd(j)
     .+sdmu*USQMIX(1,j+3)*DSQMIX(k,i)*VCKM(j,i)*mfu(j)
     . +USQMIX(1,i)*mfu(i)**2*DSQMIX(k,j)*VCKM(i,j)*1d0/tanbeta
     . +USQMIX(1,j)*mfd(i)**2*DSQMIX(k,i)*VCKM(j,i)*tanbeta
     . +(tanbeta+1d0/tanbeta)*USQMIX(1,i+3)
     . *DSQMIX(k,j+3)*mfu(i)*mfd(j)*VCKM(i,j))

c---- compare with SLHA if mass needs to be defined or not


      zwi=zwi+g2ew/dsqrt(2d0)/amw*USQMIX(1,i)*DSQMIX(k,j)
     .*amw**2*sin(2d0*DATAN(tanbeta))*VCKM(i,j)
      enddo
      enddo
      gsqsqchaHi(k)=zwi
      
      enddo

c----- charged Higgs neutralino chargino coupling
      do i=1,4
      do k=1,2
      QchncL(i,k)=dcos(Datan(tanbeta))*(ZZ(i,4)*VV(k,1)+1d0/dsqrt(2d0)
     .		  *VV(k,2)*(ZZ(i,2)+g1ew/g2ew*ZZ(i,1)))
      QchncR(i,k)=dsin(Datan(tanbeta))*(ZZ(i,3)*UU(k,1)-1d0/dsqrt(2d0)
     .		  *UU(k,2)*(ZZ(i,2)+g1ew/g2ew*ZZ(i,1)))
      enddo
      enddo

c----- W exchange
c----- W with chargino
      sigbtau=0d0
      sigbjets=0d0
      sigbmu=0d0
      sigjetb=0d0
      sigjetjet=0d0
      sigjetmu=0d0
      sigjettau=0d0
      sigbbjet=0d0
      sigjetelec=0d0
      sigbelec=0d0



      if(amsupq(1).le.amneut(1))then
      totwidth=0d0
      print*,"Four-body decay kinematically forbidden"
      return
      endif 


      if(amsupq(1)-amneut(1).gt.amb+amw)then
      totwidth=0d0
      print*,"mstop>mw+mb+mneut, three body decay not implemented"
      stop
      endif 

c------ ramona changed 21/8/14    
c------ check whether stop is NLSP
	if(amsupq(1).gt.amsdownq(1)) then
	print*, "Stop not NLSP"
        stop
        endif
        if(amsupq(1).gt.amslepton(1))then
        print*, "Stop not NLSP"
        stop
        endif
       if(amsupq(1).gt.amsneutrino(1))then
        print*, "Stop not NLSP"
        stop
        endif
        if(amsupq(1).gt.amchar(1))then
        print*, "Stop not NLSP"
        stop
        endif
	if(amsupq(1).gt.amneut(2))then
        print*, "Stop not NLSP"
        stop
        endif
c------ end ramona changed
c-------------------------------------------
c-----FIRST CASE INITIAL: B_QUARK FINAL: TAU
c-------------------------------------------
      Etot=amsupq(1)
      Mfin(1)=amb
      Mfin(2)=0d0
      Mfin(3)=amtau
      MFin(4)=amneut(1)
      
     

      
c---- simple Monte Carlo
      iiter=100000
      if(etot.gt.(mfin(1)+mfin(2)+mfin(3)+mfin(4)))then
!       Go To 14
      do i=1, iiter
      call SD_rambo(4,etot,mfin,pfout,wt)

c---- take the momenta of Rambo
      do j=1,4
      p1(j)=pfout(j,1)	!quark
      p2(j)=pfout(j,2)	!fermion
      p3(j)=pfout(j,3)	!antifermion
      p4(j)=pfout(j,4)	!neutralino
   
      pInt(j)=p2(j)+p3(j)+p4(j)	!intermediate particle (chargino)
      pfer(j)=p1(j)+p2(j)+p3(j) !intemdeiate particle(fermion such as top)
      end do
      
! define scalar products in amplitude
      p1p2=SD_scal(p1,p2)
      p1p3=SD_scal(p1,p3)
      p1p4=SD_scal(p1,p4)
      p2p3=SD_scal(p2,p3)
      p2p4=SD_scal(p2,p4)
      p3p4=SD_scal(p3,p4)
      p1pI=SD_scal(p1,pInt)
      p2pI=SD_scal(p2,pInt)
      p3pI=SD_scal(p3,pInt)
      p4pI=SD_scal(p4,pInt)
      pI2=SD_scal(pInt,pInt)
      
      m1=Mfin(1)
      m3=amtau
      m4=amneut(1)
! W amplitude with chargino as intermediate particle and final state tau and bottom
      ampwchabtau=0d0

c--- sum up different charginos as intermediate particles
      do j=1,2
      do k=1,2
      g1R1=DcharL(3,1,j)
      g1L1=EcharR(3,1,j)
      g1R2=DcharL(3,1,k)
      g1L2=EcharR(3,1,k)
      g3L=-g2ew/dSqrt(2d0)
      g2L1=CchaneutL(1,j)
      g2L2=CchaneutL(1,k)
      g2R1=CchaneutR(1,j)
      g2R2=CchaneutR(1,k)
	
       ampwchabtau=ampwchabtau+ampwchar(g1L1, g1L2, g2l1, g2l2, g1r1, 
     . g1r2, g2r1, g2r2, g3l, m1, m3, m4, 
     . amchar(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(pI2-amchar(j)**2)
     - /(pI2-amchar(k)**2)/
     - Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
c----- multiply with propagators

      end do
      enddo
      


c---- slepton contributions
      ampinterferenzwslep=0d0
      ampslepbtau=0d0

      do j=1,2
      do k=1,2
      do jj=1,6
      do kk=1,6
      g1R1=DcharL(3,1,j)
      g1L1=EcharR(3,1,j)
      g1R2=DcharL(3,1,k)
      g1L2=EcharR(3,1,k)
      g2r1=cchaneutslep(3,k,kk)
      g2r2=cchaneutslep(3,j,jj)
      g3l1=gneutelecl(3,1,kk)
      g3l2=gneutelecl(3,1,jj)
      g3r1=gneutelecr(3,1,kk)
      g3r2=gneutelecr(3,1,jj)
      g2l1=0d0
      g2l2=0d0
c------- selectron,... contribtuions
      ampslepbtau=ampslepbtau+ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)

      
      end do
      
c---- interference slepton and sneutrino diagrams
       gXb2R=cchaneutslep(3,j,jj)
       gXc2L=-dchalepsneut(3,k)
       gXc2R=-echalepsneut(3,k)
       gbd3R=gneutneut(3,1)
       gcd3L=gneutelecl(3,1,jj)
       gcd3R=gneutelecr(3,1,jj)
       gXb2l=0d0

      ampslepbtau=ampslepbtau+ampintslepsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXb2l,gXb2R, gXc2L, gXc2R,gbd3R, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2)


c----- interference w slepton diagrams
      g2l1=CchaneutL(1,k)
      g2r1=CchaneutR(1,k)
      gXb2R=cchaneutslep(3,j,jj)
      gcd3L=gneutelecl(3,1,jj)
      gcd3R=gneutelecr(3,1,jj)
      g3L=-g2ew/dSqrt(2d0)
      gxb2l=0d0 
      
      ampinterferenzwslep=ampinterferenzwslep+
     . 2d0*RealPart(ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2, gxb2l,gXb2R, g2l1, g2r1, g3L, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      
      end do
      
c----- sneutrino diagram
      g3l1=0d0
      g2l1=-dchalepsneut(3,j)
      g2l2=-dchalepsneut(3,k)
      g2r1=-echalepsneut(3,j)
      g2r2=-echalepsneut(3,k)
      g3l2=0d0
      g3r1=gneutneut(3,1)
      g3r2=gneutneut(3,1)
 


      ampslepbtau=ampslepbtau+ampcharsneutrino(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2, g3l1,g3l2, g3r1, g3r2, m1,
     . m3, m4, amchar(j),
     .  amchar(k), p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(Pi2-amchar(j)**2)
     ./(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsneutrino(3)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2)


      




c---- interference with W-contribution
      gXc2L=dchalepsneut(3,j)
      gXC2r=echalepsneut(3,j)
      g2l1=CchaneutL(1,k)
      g2r1=CchaneutR(1,k)
      g3L=-g2ew/dSqrt(2d0)
      gbd3R=gneutneut(3,1)
     
      ampinterferenzwslep=ampinterferenzwslep+
     . 2d0*RealPart(ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, g2l1, g2r1, gXc2R, g3L,
     . gbd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(3)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      end do
      end do
      
c----- contribtutions from diagram invlolving a SM fermion
      ampSMfermion=0d0
c----- SM fermion matrix elemnt squared
       p1pfer=SD_scal(p1,pfer)
      p2pfer=SD_scal(p2,pfer)
      p3pfer=SD_scal(p3,pfer)
      p4pfer=SD_scal(p4,pfer)
      pfer2=SD_scal(pfer,pfer)
      pIpfer=SD_scal(pInt, pfer)

     
      
      do j=1,3
      do k=1,3

      g1L1=-g2ew/dsqrt(2d0)*Vckm(j,3)
      g1L2=-g2ew/dsqrt(2d0)*Vckm(k,3)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2l2=gneutur(k,1,1)
      g2r2=gneutul(k,1,1)
      g3l=-g2ew/dsqrt(2d0)
      


c------ diagram involovin a W boson and a SM fermion
      ampSMfermion=ampSMfermion+ampSMfer(g1L1, g1L2, g2l1, g2l2, 
     . g2r1, g2r2,g3l,
     . m1, m3, m4, mfu(j), mfu(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2)/(pfer2-mfu(j)**2)
     . /(pfer2-mfu(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2


      end do
      enddo
c---- interference  W chargino - W SM fermion diagrams
      

      do k=1,2
      do j=1,3
      g1r2=DcharL(3,1,k)
      g1l2=EcharR(3,1,k)
      g2L2=CchaneutL(1,k)
      g2R2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)
      g1L1=-g2ew/dsqrt(2d0)*Vckm(j,3)
!       g1L2=-g2ew/dsqrt(2d0)*Vckm(3,k)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)



      ampSMfermion=ampSMfermion+2d0*ampintWSMferWchar(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l,  m1, m3, m4, 
     .mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2 

! c---- interfernce slepton  and W SM diagrams
! 
! 
      do kk=1,6
      g2l2=0d0
      g2r2=cchaneutslep(3,k,kk)
      g3l1=-g2ew/dsqrt(2d0)
      g3r2=gneutelecr(3,1,kk)
      g3l2=gneutelecl(3,1,kk)
      ampSMfermion=ampSMfermion+2d0*RealPart(
     . ampintWSMferslep(g1L1, g1L2, g1r2,
     .  g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,m1,m3,
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2))


      enddo
c---- interfernce sneutrino  and W SM diagrams
      g2l2=-dchalepsneut(3,k)
      g2r2=-echalepsneut(3,k)
      g3l2=0d0
      g3r2=gneutneut(3,1)
      ampSMfermion=ampSMfermion+2d0*RealPart(
     . ampintWSMfersneut(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2, m1, m3, 
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2))
 
       enddo
! 
        end do
! 	
c------- the diagram with W and squark propagator

      ampsquarkpropW=0d0
      ampchahi=0d0
      do j=1,6
      do k=1,6
      g1l1=gsqsqW(j)
      g1l2=gsqsqW(k)
      g2l1=gneutdr(3,1,j)
      g2l2=gneutdr(3,1,k)
      g2r1=gneutdl(3,1,j)
      g2r2=gneutdl(3,1,k)
      g3l=-g2ew/dsqrt(2d0)
      

      ampsquarkpropW=ampsquarkpropW+ampvecscalar(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,g3l,0d0,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, amw)/
     . Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)

      g1l1=gsqsqchaHi(j)
      g1l2=gsqsqchaHi(k)
      g2r1=gneutdl(3,1,j)
      g2r2=gneutdl(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=gneutdr(3,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta

      
c----- diagram charged Higgs squark propagator
      ampchahi=ampchahi+amphiggssquark(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,0d0,g3r1,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4)/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/
     . (2d0*p2p3+m3**2-amch**2)**2


	
       

! 
! c---- interference Higgs-squark with W squark diagram
      g1l1=gsqsqW(j)
      g1l2=gsqsqchaHi(k)
      g2r1=gneutdl(3,1,j)
      g2r2=gneutdl(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=gneutdr(3,1,k)
      g3l1=-g2ew/dsqrt(2d0)
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta

      ampchahi=ampchahi+2d0*RealPart(
     .ampinthiggssquarkWsquark(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,g3l1,0d0,g3r1,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4)/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/
     . (2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
      enddo
    
! 
c----- interference charged Higgs squark- W chargino diagram
      do j=1,6
      do k=1,2
      g1l1=gsqsqchaHi(j)
      g1l2=EcharR(3,1,k)
      g1r2=DcharL(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=CchaneutL(1,k)
      g2r1=gneutdl(3,1,j)
      g2r2=CchaneutR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/dsqrt(2d0)

      ampchahi=ampchahi+2d0*RealPart(ampintchaHisquarkWchar(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
      enddo
! 
! c------ interference charged Higgs squark -W SM fermion diagram
      do j=1,6
      do k=1,3
      g1l1=gsqsqchaHi(j)
      g1l2=gneutur(k,1,1)
      g1r2=gneutul(k,1,1)
      g2l1=gneutdr(3,1,j)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(k,3)
      g2r1=gneutdl(3,1,j)
      g2r2=0d0
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/dsqrt(2d0)

      ampchahi=ampchahi+2d0*RealPart(ampintchaHisquarkWSMfer(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(k))/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pfer2-mfu(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
      enddo

c---- interference charged Higgs squark- slepton chargino diagram
      do j=1,6
      do k=1,2
      do kk=1,6
      g1l1=gsqsqchaHi(j)
      g1l2=EcharR(3,1,k)
      g1r2=DcharL(3,1,k)
      g2r1=gneutdr(3,1,j)
      g2l1=gneutdl(3,1,j)
      g2r2=cchaneutslep(3,k,kk)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=gneutelecl(3,1,kk)
      g3r2=gneutelecr(3,1,kk)


      ampchahi=ampchahi
     . +2d0*RealPart(ampintchaHisquarkslepton(g1l1, g1l2, 
     . g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/(pi2-amchar(k)**2)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)
     ./(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      g2l2=dchalepsneut(3,k)
      g2r2=echalepsneut(3,k)
      g3r2=gneutneut(3,1)


            ampchahi=ampchahi
     . -2d0*RealPart(ampintchaHisquarksneutrino(g1l1, g1l2, 
     . g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/(pi2-amchar(k)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2)
     ./(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo


! 
! c----- interference W squark- W chargino diagrams
      do j=1,6
      do k=1,2
      g1l1=gsqsqW(j)
      g1l2=EcharR(3,1,j)
      g1R2=DcharL(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=CchaneutL(1,k)
      g2r1=gneutdl(3,1,j)
      g2r2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)

      ampsquarkpropW=ampsquarkpropW+2d0*RealPart(
     . ampintWcharWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pi2-amchar(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2)
      enddo
      enddo
! c------ interfernce W squark - W SM Fermion
      do j=1,6
      do k=1,3
      g1l1=gsqsqW(j)
      g1l2=gneutur(k,1,1)
      g1R2=gneutul(k,1,1)
      g2l1=gneutdr(3,1,j)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(k,3)
      g2r1=gneutdl(3,1,j)
      g3l=-g2ew/dsqrt(2d0)
      ampsquarkpropW=ampsquarkpropW+2d0*RealPart(
     . ampintWferWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pfer2-mfu(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2) 
      enddo
      enddo
! c--------interfernce W squark Slepton diagram
      do j=1,6
      do k=1,2
      do kk=1,6
      g1l1=gsqsqW(j)
      g1l2=EcharR(3,1,k)
      g1r2=DcharL(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=0d0
      g2r1=gneutdl(3,1,j)
      g2r2=cchaneutslep(3,k,kk)
      g3l1=-g2ew/dsqrt(2d0)
      g3l2=gneutelecl(3,1,kk)
      g3r2=gneutelecr(3,1,kk)

          
   
      ampsquarkpropW=ampsquarkpropW+2d0*RealPart(
     .ampintcharslepWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     ./(pi2-amchar(k)**2)/(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth) )

      enddo
!       
c---- interference W squark- chargino sneutrino diagram      
      g2l2=dchalepsneut(3,k)
      g2r2=echalepsneut(3,k)
      g3r2=gneutneut(3,1)
      g3l2=0d0

      ampsquarkpropW=ampsquarkpropW+2d0*
     . RealPart(ampintcharsneutWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     . /(pi2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(3)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
! 
! c----- charged Higgs diagram: SM fermion (only top diagrams because no FV on vertex with charged Higgs)
      
      g2l1=gneutur(3,1,1)
      g2r1=gneutul(3,1,1)
      g2r2=gneutul(3,1,1)
      g2l2=gneutur(3,1,1)
      g1r1=g2ew/dsqrt(2d0)/amw*mfu(3)*1d0/tanbeta*VCKM(3,3)
      g1r2=g2ew/dsqrt(2d0)/amw*mfu(3)*1d0/tanbeta*VCKM(3,3)
      g1l1=g2ew/dsqrt(2d0)/amw*mfd(3)*tanbeta*VCKM(3,3)
      g1l2=g2ew/dsqrt(2d0)/amw*mfd(3)*tanbeta*VCKM(3,3)
      g3r1= g2ew/dsqrt(2d0)/amw*mfe(3)*tanbeta
! 
      ampchahi=ampchahi+ampchaHiSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(3), mfu(3))
     . /(pfer2-mfu(3)**2)/(pfer2-mfu(3)**2)/(2d0*p2p3+m3**2-amch**2)**2
! 
! 
! c----- charged Higgs diagram: chargino

      do j=1,2
      do k=1,2
      g1l1=EcharR(3,1,j)
      g1r1=DcharL(3,1,j)
      g1l2=EcharR(3,1,k)
      g1r2=DcharL(3,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=-g2ew*QchncL(1,k)
      g2r2=-g2ew*QchncR(1,k)
      g3r1= g2ew/dsqrt(2d0)/amw*mfe(3)*tanbeta

      ampchahi=ampchahi+ampchaHichar(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(j), amchar(k))/
     . (pi2-amchar(j)**2)/(pI2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)**2
      enddo
      enddo
! ! c----- interference charged Higgs chargino and charged Higgs squark diagrams
      do j=1,2
      do k=1,6
      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=gsqsqchaHi(k)
      g2r1=-g2ew*QchncR(1,j)
      g2l1=-g2ew*QchncL(1,j)
      g2l2=gneutdr(3,1,k)
      g2r2=gneutdl(3,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta

      ampchahi=ampchahi+2d0*ampintchaHichargchaHisquark(g1l2,
     .  g1l1,  g1r1,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j))/(pI2-amchar(j)**2)
     . /(2d0*p2p3+m3**2-amch**2)**2
     . /(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)  
      enddo
      enddo
! c------ interference charged Higgs chargino - W squark diagrams
       
       do j=1,2
       do k=1,6

      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=gsqsqW(k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=gneutdr(3,1,k)
      g2r2=gneutdl(3,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))
      ampchahi=ampchahi+2d0*RealPart(ampintchaHichargWsquark(g1l2,
     .  g1l1,  g1r1,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi,p2pI, p3pI, p4pI, pI2,amchar(j))/(pI2-amchar(j)**2)
     ./(2d0*p2p3+m3**2-amch**2)/(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)/
     . (2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
     
       enddo
       enddo

! ! 
! c----- interfernce charged Higgs chargino-W chargino
       do j=1,2
      do k=1,2
      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=CchaneutL(1,k)
      g2r2=CchaneutR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))

      ampchahi=ampchahi+2d0*RealPart(ampintchaHiWchar(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))
     . /(pi2-amchar(j)**2)/(pi2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
      enddo
! ! ! 
! c----- interfernce charged Higgs chargino-W Sm fermion diagrams
       do j=1,2
      do k=1,3
      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=gneutur(k,1,1)
      g1r2=gneutul(k,1,1)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=-g2ew/dsqrt(2d0)*Vckm(k,3)
      g2r2=-g2ew/dsqrt(2d0)*Vckm(k,3)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))

      ampchahi=ampchahi-2d0*RealPart(ampintchaHiWSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, pIpFer, mfu(k), amchar(j))
     . /(pfer2-mfu(k)**2)/(pi2-amchar(j)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
! 
! c----- interference charged Higgs chargino+ slepton diagram
      do kk=1,6
      do j=1,2
      do k=1,2
      
      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
 
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2r2=cchaneutslep(3,k,kk)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=gneutelecl(3,1,kk)
      g3r2=gneutelecr(3,1,kk)

      
      ampchahi=ampchahi+2d0*ampintchaHichargslepton(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))
     ./(pi2-amchar(j)**2)
     . /(pi2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)


      enddo
      enddo
      enddo


! 
! c---- interference charged Higgs chargino+sneutrino diagrams
      do j=1,2
      do k=1,2
      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=-dchalepsneut(3,k)
      g2r2=-echalepsneut(3,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3r2=gneutneut(3,1)




      ampchahi=ampchahi+2d0*ampintchaHichargsneutrino(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     . m1,m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))/
     .(pi2-amchar(j)**2)/(pi2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p4+m4**2-amsneutrino(3)**2)

 
      enddo
      enddo

! 
! c---- interference charged Higgs Sm fermion with slepton diagram
      do k=1,2
      do kk=1,6
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r2=cchaneutslep(3,k,kk)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3r2=gneutelecr(3,1,kk)
      g3l2=gneutelecl(3,1,kk)

      ampchahi=ampchahi+2d0*ampintchaHiSmferslepton(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfu(3), amchar(k))
     ./(pi2-amchar(k)**2)/(pfer2-mfu(3)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2) 
      enddo
      enddo

c---- interference charged Higgs Sm fermion with sneutrino diagram
      do k=1,2
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)
      g2l2=dchalepsneut(3,k)
      g2r2=echalepsneut(3,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3r2=gneutneut(3,1)


      ampchahi=ampchahi+2d0*ampintchaHiSmfersneutrino(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfu(3), amchar(k))
     ./(pi2-amchar(k)**2)/(pfer2-mfu(3)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p4+m4**2-amsneutrino(3)**2) 
      enddo
c----- interference charged Higgs chargino and charged Higgs SM fermion
      do k=1,2
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)
      g2l2=-g2ew*QchncL(1,k)
      g2r2=-g2ew*QchncR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta


      ampchahi=ampchahi-2d0*ampintchaHiSmferchaHicharg(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfu(3), amchar(k))
     ./(pi2-amchar(k)**2)/(pfer2-mfu(3)**2)
     . /(2d0*p2p3+m3**2-amch**2)**2
      enddo
! 
c----- interference charged Higgs squark and charged Higgs SM fermion diagrams
      do k=1,6
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1l2=gsqsqchaHi(k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)
      g2l2=gneutdr(3,1,k)
      g2r2=gneutdl(3,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta


      ampchahi=ampchahi+2d0*ampintchaHiSmferchaHisquark(g1l2,
     . g1r1, g1l1,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer,p4pfer, pfer2, mfu(3))
     ./(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)/(pfer2-mfu(3)**2)
     . /(2d0*p2p3+m3**2-amch**2)**2
      enddo
! ! 
! c----- interference charged Higgs SM fermion and W chargino diagrams
      do k=1,2
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1r2=EcharR(3,1,k)
      g1l2=DcharL(3,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2l2=CchaneutL(1,k)
      g2r2=CchaneutR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))

      ampchahi=ampchahi-2d0*RealPart(ampintchaHiSmferWcharg(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pipfer, amchar(k), mfu(3))
     . /(pfer2-mfu(3)**2)/(pi2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
! ! 
! c---- interference charged Higgs SM fermion and W squark diagrams
      do k=1,6
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1l2=gsqsqW(k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2l2=gneutdr(3,1,k)
      g2r2=gneutdl(3,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))


	ampchahi=ampchahi-2d0*RealPart(ampintchaHiSmferWsquark(g1l1,
     . g1r1, g1l2, g2l1, g2r1, g2l2,g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, mfu(3))  
     . /(pfer2-mfu(3)**2)/(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
c---- interference SM fermion W and SM fermion charged Higgs diagrams
      do k=1,3
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1l2=gneutur(k,1,1)
      g1r2=gneutur(k,1,1)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)
      g2l2=-g2ew/dsqrt(2d0)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))


      ampchahi=ampchahi+2d0*Realpart(ampintchaHiSmferWSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, mfu(3), mfu(k))    
     . /(pfer2-mfu(3)**2)/(pfer2-mfu(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo 


      sigbtau=sigbtau+wt*ampwchabtau/iiter
     .+wt*ampslepbtau/iiter+wt*ampinterferenzwslep/iiter
     .+wt*ampSMfermion/iiter
     .+wt*ampsquarkpropW/iiter
     .+wt*ampchahi/iiter


       enddo
!         GoTo 13
      else
      sigbtau=0d0
      print*, 
     ."decay stop->neutralino tau b neutrino kinematically forbidden"
      endif

     

c-------------------------------------------------
c------ SECOND CASE: INITIAL B_QUARK, FINAL B_QUARK
c-------------------------------------------------
      Etot=amsupq(1)
      Mfin(1)=amb
      Mfin(2)=0d0
      Mfin(3)=amb
      MFin(4)=amneut(1)

      if(etot.gt.(amb+amb+amneut(1)))then
      do i=1, iiter
      call SD_rambo(4,etot,mfin,pfout,wt)

c---- take the momenta of Rambo 
      do j=1,4
      p1(j)=pfout(j,1)	!quark
      p2(j)=pfout(j,2)	!fermion
      p3(j)=pfout(j,3)	!antifermion
      p4(j)=pfout(j,4)	!neutralino
      pInt(j)=p2(j)+p3(j)+p4(j)	!intermediate particle (chargino)
      pfer(j)=p2(j)+p3(j)+p1(j) !intermediate partcile (SM fermion e.g. top)
      end do

! define scalar products in amplitude
      p1p2=SD_scal(p1,p2)
      p1p3=SD_scal(p1,p3)
      p1p4=SD_scal(p1,p4)
      p2p3=SD_scal(p2,p3)
      p2p4=SD_scal(p2,p4)
      p3p4=SD_scal(p3,p4)
      p1pI=SD_scal(p1,pInt)
      p2pI=SD_scal(p2,pInt)
      p3pI=SD_scal(p3,pInt)
      p4pI=SD_scal(p4,pInt)
      pI2=SD_scal(pInt,pInt)
      
      m1=Mfin(1)
      m3=amb
      m4=amneut(1)

      ampwchabbjet=0d0
c--- sum up different charginos as intermediate particles
      do j=1,2
      do k=1,2
      g1R1=DcharL(3,1,j)
      g1L1=EcharR(3,1,j)
      g1R2=DcharL(3,1,k)
      g1L2=EcharR(3,1,k)
      g3L=-g2ew/dSqrt(2d0)
      g2L1=CchaneutL(1,j)
      g2L2=CchaneutL(1,k)
      g2R1=CchaneutR(1,j)
      g2R2=CchaneutR(1,k)
     



c--- propagators, colour factor and summing up bc and bu contribution
       ampwchabbjet=ampwchabbjet+3d0*(Vckm(2,3)**2+VCKM(1,3)**2)*
     . ampwchar(g1L1, g1L2, g2l1, g2l2, g1r1, 
     . g1r2, g2r1, g2r2, g3l, m1, m3, m4, 
     . amchar(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     . /(pI2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2

      end do
      end do

c------ squark contributions
      ampbsquarkb=0d0
      do j=1,2
      do k=1,2
      do jj=1,6
      do ifer=1,2
      do kk=1,6
      
      g1R1=DcharL(3,1,j)
      g1L1=EcharR(3,1,j)
      g1R2=DcharL(3,1,k)
      g1L2=EcharR(3,1,k)
      g2l1=-DcharL(3,jj,j)
      g2l2=-DcharL(3,kk,k)
      g2r1=-EcharR(3,jj,j)
      g2r2=-EcharR(3,kk,k)
      g3l1=gneutur(ifer,1,jj)
      g3l2=gneutur(ifer,1,kk)
      g3r1=gneutul(ifer,1,jj)
      g3r2=gneutul(ifer,1,kk)

      
      ampbsquarkb=ampbsquarkb+3d0*ampcharsup(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2,g3l1, g3l2, g3r1, g3r2, 
     . m1, m3, m4, amchar(j),
     .  amchar(k), p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     ./(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p2p4+m4**2-amsupq(kk)**2)

      g2l1=fcharr(ifer, jj, j)
      g2l2=fcharr(ifer, kk, k)
      g2r1=ccharl(ifer, jj, j)
      g2r2=ccharl(ifer, kk, k)
      g3l1=gneutdL(3,1,jj)
      g3l2=gneutdL(3,1,kk)
      g3r1=gneutdr(3,1,jj)
      g3r2=gneutdr(3,1,kk)

      ampbsquarkb=ampbsquarkb+3d0*ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)


      g2l1=-DcharL(3,jj,j)
      g2r1=-EcharR(3,jj,j)
      g3r1=gneutuL(ifer,1,jj)
      g3l1=gneutur(ifer,1,jj)
      g2l2=fcharr(ifer, kk, k)
      g2r2=ccharl(ifer, kk, k)
      g3l2=gneutdL(3,1,kk)
      g3r2=gneutdr(3,1,kk)

c----- squark exchange diagrams, same amplitude than slepton
c----- and sneutrino diagrams but couplings need to be changed
      ampbsquarkb=ampbsquarkb
     .+6d0*ampcharsupsdownint(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)    
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)


c---- interfernz squark-W-Terme
      enddo
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      g3L=-g2ew/dSqrt(2d0)*Vckm(ifer,3)
      gXb2r=ccharl(ifer, jj, j)
      gXb2l=fcharr(ifer, jj, j)
      gcd3l=gneutdl(3,1,jj)
      gcd3r=gneutdr(3,1,jj)

      ampbsquarkb=ampbsquarkb+
     . 6d0*RealPart(ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2, gxb2l,gXb2R, g2l1, g2r1, g3L, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
 

      gXc2R=dcharl(3, jj, j)
      gXC2L=echarr(3, jj, j)
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
 
      g3L=-g2ew/dSqrt(2d0)*VCKM(ifer,3)
      gbd3R=gneutul(ifer, 1,jj)


      ampbsquarkb=ampbsquarkb+
     .6d0*RealPart(ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, g2l1, g2r1, gXc2R, g3L,
     . gbd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      
      enddo
      enddo
      enddo
      enddo


c----- contribtutions from diagram invlolving a SM fermion
      ampSMfermion=0d0
c----- SM fermion matrix elemnt squared
       p1pfer=SD_scal(p1,pfer)
      p2pfer=SD_scal(p2,pfer)
      p3pfer=SD_scal(p3,pfer)
      p4pfer=SD_scal(p4,pfer)
      pfer2=SD_scal(pfer,pfer)
      pIpfer=SD_scal(pInt, pfer)



      do j=1,3
      do k=1,3
      
      g1L1=-g2ew/dsqrt(2d0)*Vckm(j,3)
      g1L2=-g2ew/dsqrt(2d0)*Vckm(k,3)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2l2=gneutur(k,1,1)
      g2r2=gneutul(k,1,1)
      g3l=-g2ew/dsqrt(2d0)


      ampSMfermion=ampSMfermion+3d0*(Vckm(1,3)**2+VCKM(2,3)**2)
     .*ampSMfer(g1L1,g1L2, g2l1, g2l2, 
     . g2r1, g2r2,g3l,
     . m1, m3, m4, mfu(j), mfu(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2)/(pfer2-mfu(j)**2)
     . /(pfer2-mfu(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
      
      end do
      
c---- interference term with W-chargino diagram


      do k=1,2
      g1r2=DcharL(3,1,k)
      g1l2=EcharR(3,1,k)
      g2L2=CchaneutL(1,k)
      g2R2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)
     
      ampSMfermion=ampSMfermion
     .+6d0*(Vckm(1,3)**2+VCKM(2,3)**2)*ampintWSMferWchar(g1L1,g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l,  m1, m3, m4, 
     .mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2 

c---- interfernce with sdown contribtuions


      do kk=1,6
      do ifer=1,2
      g2l2=fcharr(ifer, kk, k)
      g2r2=ccharl(ifer, kk, k)
      g3l1=-g2ew/dsqrt(2d0)*Vckm(ifer,3)
      g3r2=gneutdR(3,1,kk)
      g3l2=gneutdL(3,1,kk)
      ampSMfermion=ampSMfermion+6d0*RealPart(
     . ampintWSMferslep(g1L1, g1L2, g1r2,
     .  g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,m1, m3,
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2))

      
      g2l2=-Dcharl(3,kk,k)
      g2r2=-EcharR(3,kk,k)
      g3l2=gneutur(ifer,1,kk)
      g3r2=gneutul(ifer,1,kk)
      ampSMfermion=ampSMfermion+6d0*RealPart(
     . ampintWSMfersneut(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2, m1, m3, 
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p2p4+m4**2-amsupq(kk)**2))
      enddo
      end do
      enddo
      end do


c---- diagram with W squark propagator
      ampsquarkpropW=0d0

      do j=1,6
      do k=1,6
      g1l1=gsqsqW(j)
      g1l2=gsqsqW(k)
      g2l1=gneutdr(3,1,j)
      g2l2=gneutdr(3,1,k)
      g2r1=gneutdl(3,1,j)
      g2r2=gneutdl(3,1,k)
      g3l=-g2ew/dsqrt(2d0)
      

     
      ampsquarkpropW=ampsquarkpropW+3d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*ampvecscalar(g1l1,g1l2, 
     . g2l1, g2l2, g2r1, g2r2,g3l,0d0,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, amw)/
     . Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
      enddo
      enddo
c---- diagrams with charged Higgs not necessary since no FV interactions

c----- interference W squark- W chargino diagrams
      do j=1,6
      do k=1,2
      g1l1=gsqsqW(j)
      g1l2=EcharR(3,1,j)
      g1R2=DcharL(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=CchaneutL(1,k)
      g2r1=gneutdl(3,1,j)
      g2r2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)

      ampsquarkpropW=ampsquarkpropW+
     . 6d0*(VCKM(1,3)**2+VCKM(2,3)**2)*RealPart(
     . ampintWcharWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pi2-amchar(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2)
      enddo
      enddo
c------ interfernce W squark - W SM Fermion
      do j=1,6
      do k=1,3
      do ifer=1,2
      g1l1=gsqsqW(j)
      g1l2=gneutur(k,1,1)
      g1R2=gneutul(k,1,1)
      g2l1=gneutdr(3,1,j)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(k,3)
      g2r1=gneutdl(3,1,j)
      g3l=-g2ew/dsqrt(2d0)*(VCKM(ifer,3))
      ampsquarkpropW=ampsquarkpropW+6d0*RealPart(
     . ampintWferWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pfer2-mfu(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2) 
      enddo
      enddo
      enddo
c--------interfernce W squark Sdown diagram
      do j=1,6
      do k=1,2
      do kk=1,6
      do ifer=1,2
      g1l1=gsqsqW(j)
      g1l2=EcharR(3,1,k)
      g1R2=DcharL(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=fcharr(ifer, kk, k)
      g2r1=gneutdl(3,1,j)
      g2r2=ccharl(ifer, kk, k)
      g3l1=-g2ew/dsqrt(2d0)*VCKM(ifer,3)
      g3l2=gneutdl(3,1,kk)
      g3r2=gneutdr(3,1,kk)
          
   
      ampsquarkpropW=ampsquarkpropW+6d0*RealPart(
     .ampintcharslepWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     ./(pi2-amchar(k)**2)/(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth) )

      
      
c---- interference W squark- chargino sup diagram      

      g1R2=DcharL(3,1,k)
      g1L2=EcharR(3,1,k)
      g2l2=-DcharL(3,kk,k)
      g2r2=-EcharR(3,kk,k)
      g3l2=gneutuL(ifer,1,kk)
      g3r2=gneutur(ifer,1,kk)

      ampsquarkpropW=ampsquarkpropW+6d0*
     . RealPart(ampintcharsneutWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     . /(pi2-amchar(k)**2)/(2d0*p2p4+m4**2-amsupq(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo
      enddo
       ampchahib=0d0
      do k=1,6
      do j=1,6
      
 
      g1l1=gsqsqchaHi(j)
      g1l2=gsqsqchaHi(k)
      g2r1=gneutdl(3,1,j)
      g2r2=gneutdl(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=gneutdr(3,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta

      
c----- diagram charged Higgs squark propagator
      ampchahib=ampchahib+(VCKM(1,3)**2+VCKM(2,3)**2)*
     . 3d0*amphiggssquark(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,0d0,g3r1,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4)/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/
     . (2d0*p2p3+m3**2-amch**2)**2


	
       

! 
! ! c---- interference Higgs-squark with W squark diagram
      g1l1=gsqsqW(j)
      g1l2=gsqsqchaHi(k)
      g2r1=gneutdl(3,1,j)
      g2r2=gneutdl(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=gneutdr(3,1,k)
      g3l1=-g2ew/dsqrt(2d0)
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta

      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)*RealPart(
     .ampinthiggssquarkWsquark(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,g3l1,0d0,g3r1,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4)/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/
     . (2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
! 
      enddo
      enddo
    
! 
! c----- interference charged Higgs squark- W chargino diagram
      do j=1,6
      do k=1,2
      g1l1=gsqsqchaHi(j)
      g1l2=EcharR(3,1,k)
      g1r2=DcharL(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=CchaneutL(1,k)
      g2r1=gneutdl(3,1,j)
      g2r2=CchaneutR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/dsqrt(2d0)

     
      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)*2)*RealPart(
     . ampintchaHisquarkWchar ( g1l1 , g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
      enddo
! 
! c------ interference charged Higgs squark -W SM fermion diagram
      do j=1,6
      do k=1,3
      g1l1=gsqsqchaHi(j)
      g1l2=gneutur(k,1,1)
      g1r2=gneutul(k,1,1)
      g2l1=gneutdr(3,1,j)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(3,k)
      g2r1=gneutdl(3,1,j)
      g2r2=0d0
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/dsqrt(2d0)

     
      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)*RealPart(
     . ampintchaHisquarkWSMfer( g1l1 , g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(k))/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pfer2-mfu(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
      enddo

c---- interference charged Higgs squark- slepton chargino diagram
      do j=1,6
      do k=1,2
      do kk=1,6
      do ifer=1,2
      g1l1=gsqsqchaHi(j)
      g1l2=EcharR(3,1,k)
      g1r2=DcharL(3,1,k)
      g2r1=gneutdr(3,1,j)
      g2l1=gneutdl(3,1,j)
      g2r2=ccharl(ifer,kk,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta*VCKM(ifer,3)
      g3l2=gneutdl(3,1,kk)
      g3r2=gneutdr(3,1,kk)


      ampchahib=ampchahib
     . +6d0*RealPart(ampintchaHisquarkslepton(g1l1, g1l2, 
     . g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/(pi2-amsdownq(k)**2)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)
     ./(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      
      g2l2=-DcharL(3,kk,j)
      g2r2=-EcharR(3,kk,j)
      g3r2=gneutul(ifer,1,kk)
      g3l2=gneutur(ifer,1,kk)


            ampchahib=ampchahib
     . -6d0*RealPart(ampintchaHisquarksup(g1l1, g1l2, 
     . g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/(pi2-amchar(k)**2)
     . /(2d0*p2p4+m4**2-amsupq(kk)**2)
     ./(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo
      enddo


      g2l1=gneutur(3,1,1)
      g2r1=gneutul(3,1,1)
      g2r2=gneutul(3,1,1)
      g2l2=gneutur(3,1,1)
      g1r1=g2ew/dsqrt(2d0)/amw*mfu(3)*1d0/tanbeta*VCKM(3,3)
      g1r2=g2ew/dsqrt(2d0)/amw*mfu(3)*1d0/tanbeta*VCKM(3,3)
      g1l1=g2ew/dsqrt(2d0)/amw*mfd(3)*tanbeta*VCKM(3,3)
      g1l2=g2ew/dsqrt(2d0)/amw*mfd(3)*tanbeta*VCKM(3,3)
      g3r1= g2ew/dsqrt(2d0)/amw*mfd(3)*tanbeta
! 
      ampchahib=ampchahib+3d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     . *ampchaHiSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(3), mfu(3))
     . /(pfer2-mfu(3)**2)/(pfer2-mfu(3)**2)/(2d0*p2p3+m3**2-amch**2)**2

! 
! 
! c----- charged Higgs diagram: chargino

      do j=1,2
      do k=1,2
      g1l1=EcharR(3,1,j)
      g1r1=DcharL(3,1,j)
      g1l2=EcharR(3,1,k)
      g1r2=DcharL(3,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=-g2ew*QchncL(1,k)
      g2r2=-g2ew*QchncR(1,k)
      g3r1= g2ew/dsqrt(2d0)/amw*mfd(3)*tanbeta

      ampchahib=ampchahib+3d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     . *ampchaHichar(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(j), amchar(k))/
     . (pi2-amchar(j)**2)/(pI2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)**2
      enddo
      enddo

c----- interference charged Higgs chargino and charged Higgs squark diagrams
      do j=1,2
      do k=1,6
      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=gsqsqchaHi(k)
      g2r1=-g2ew*QchncR(1,j)
      g2l1=-g2ew*QchncL(1,j)
      g2l2=gneutdr(3,1,k)
      g2r2=gneutdl(3,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta

    
 
 
      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*ampintchaHichargchaHisquark( g1l2 ,
     .  g1l1,  g1r1,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j))/(pI2-amchar(j)**2)
     . /(2d0*p2p3+m3**2-amch**2)**2
     . /(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)  
      enddo
      enddo
! ! c------ interference charged Higgs chargino - W squark diagrams
       
       do j=1,2
       do k=1,6

      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=gsqsqW(k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=gneutdr(3,1,k)
      g2r2=gneutdl(3,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))
     
      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*RealPart(ampintchaHichargWsquark(g1l2,
     .  g1l1,  g1r1,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi,p2pI, p3pI, p4pI, pI2,amchar(j))/(pI2-amchar(j)**2)
     ./(2d0*p2p3+m3**2-amch**2)/(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)/
     . (2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
     
       enddo
       enddo

! ! 
! c----- interfernce charged Higgs chargino-W chargino
       do j=1,2
      do k=1,2
      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=CchaneutL(1,k)
      g2r2=CchaneutR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))

     
      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*RealPart(ampintchaHiWchar(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))
     . /(pi2-amchar(j)**2)/(pi2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
      enddo
! ! ! ! 
! ! c----- interfernce charged Higgs chargino-W Sm fermion diagrams
       do j=1,2
      do k=1,3
      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=gneutur(k,1,1)
      g1r2=gneutul(k,1,1)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=-g2ew/dsqrt(2d0)*Vckm(k,3)
      g2r2=0d0
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))

     
 

 
      ampchahib=ampchahib-6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*RealPart(ampintchaHiWSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, pIpFer, mfu(k), amchar(j))
     . /(pfer2-mfu(k)**2)/(pi2-amchar(j)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
! ! 
! ! c----- interference charged Higgs chargino+ slepton diagram
      do kk=1,6
      do j=1,2
      do k=1,2
      do ifer=1,2
      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
 
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2r2=ccharl(ifer,kk,j)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta*Vckm(ifer,3)
      g3l2=gneutdl(3,1,kk)
      g3r2=gneutdl(3,1,kk)


      ampchahib=ampchahib+6d0*ampintchaHichargslepton(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))
     ./(pi2-amchar(j)**2)
     . /(pi2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)


      enddo
      enddo
      enddo
      enddo
! 
! ! 
! c---- interference charged Higgs chargino+sneutrino diagrams
c--- sneutrino diagram can be used, because there is no dependence on g3l2
      do j=1,2
      do k=1,2
      do ifer=1,2
      do kk=1,6
      g1l1=EcharR(3,1,j)
      g1r1=Dcharl(3,1,j)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=-dcharl(3,kk,k)
      g2r2=-echarr(3,kk,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta*VCKM(ifer,3)
      g3r2=gneutul(ifer,1,kk)


      

      ampchahib=ampchahib+6d0*ampintchaHichargsneutrino(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     . m1,m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))/
     .(pi2-amchar(j)**2)/(pi2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p4+m4**2-amsupq(kk)**2)

 
      enddo
      enddo
      enddo
      enddo
! 
! c---- interference charged Higgs Sm fermion with slepton diagram
      do k=1,2
      do kk=1,6
      do ifer=1,2
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r2=ccharl(3,kk,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta*VCKM(ifer,3)
      g3r2=gneutdr(3,1,kk)
      g3l2=gneutdl(3,1,kk)
      
      ampchahib=ampchahib+6d0*ampintchaHiSmferslepton(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfu(3), amchar(k))
     ./(pi2-amchar(k)**2)/(pfer2-mfu(3)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2) 
      enddo
      enddo
      enddo
c---- interference charged Higgs Sm fermion with sneutrino diagram
      do k=1,2
      do kk=1,6
      do ifer=1,2
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2l2=-dcharl(3,kk,k)
      g2r2=-echarr(3,kk,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta*VCKM(ifer,3)
      g3r2=gneutul(ifer,kk,1)


      ampchahib=ampchahib+6d0*ampintchaHiSmfersneutrino(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfu(3), amchar(k))
     ./(pi2-amchar(k)**2)/(pfer2-mfu(3)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p4+m4**2-amsupq(kk)**2) 
      enddo
      enddo
      enddo
! c----- interference charged Higgs chargino and charged Higgs SM fermion
      do k=1,2
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1l2=EcharR(3,1,k)
      g1r2=Dcharl(3,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2l2=-g2ew*QchncL(1,k)
      g2r2=-g2ew*QchncR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta


     
       ampchahib=ampchahib-6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*ampintchaHiSmferchaHicharg(g1l1 ,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfu(3), amchar(k))
     ./(pi2-amchar(k)**2)/(pfer2-mfu(3)**2)
     . /(2d0*p2p3+m3**2-amch**2)**2
      enddo
! ! 
! ! c----- interference charged Higgs squark and charged Higgs SM fermion diagrams
      do k=1,6
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1l2=gsqsqchaHi(k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2l2=gneutdr(3,1,k)
      g2r2=gneutdl(3,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta


     
       ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*ampintchaHiSmferchaHisquark(g1l2,
     . g1r1, g1l1,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer,p4pfer, pfer2, mfu(3))
     ./(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)/(pfer2-mfu(3)**2)
     . /(2d0*p2p3+m3**2-amch**2)**2
      enddo

! ! 
! c----- interference charged Higgs SM fermion and W chargino diagrams
      do k=1,2
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1r2=EcharR(3,1,k)
      g1l2=DcharL(3,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2l2=CchaneutL(1,k)
      g2r2=CchaneutR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))

     
       ampchahib=ampchahib-6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*RealPart(ampintchaHiSmferWcharg(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pipfer, amchar(k), mfu(3))
     . /(pfer2-mfu(3)**2)/(pi2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
! ! 
! c---- interference charged Higgs SM fermion and W squark diagrams
      do k=1,6
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1l2=gsqsqW(k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2l2=gneutdr(3,1,k)
      g2r2=gneutdl(3,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))


	ampchahib=ampchahib-6d0*(VCKM(1,3)**2+VCKM(2,3)**2)*RealPart(
     .ampintchaHiSmferWsquark(g1l1,
     . g1r1, g1l2, g2l1, g2r1, g2l2,g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, mfu(3))  
     . /(pfer2-mfu(3)**2)/(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
c---- interference SM fermion W and SM fermion charged Higgs diagrams
      do k=1,3
      g1l1=gneutur(3,1,1)
      g1r1=gneutul(3,1,1)
      g1r2=EcharR(3,1,k)
      g1l2=DcharL(3,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2r1=mfu(3)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(3,3)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(3,3)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))


     
      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*Realpart(ampintchaHiSmferWSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, mfu(3), mfu(k))    
     . /(pfer2-mfu(3)**2)/(pfer2-mfu(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo 





      sigbbjet=sigbbjet+wt*ampwchabbjet/iiter+wt*ampbsquarkb/iiter
     . +wt*ampSMfermion/iiter+wt*ampsquarkpropW/iiter+wt*ampchahib/iiter

      enddo
      else
      sigbbjet=0d0
      print*,
     ."Decay stop-> neutralino b b jet kinematically forbidden"
      endif
!       Go To 13
c--------------------------------------------------------------
c----- THIRD CASE: INITIAL: B_QUARK< FINAL: MASSLESS PARTICLES
c-------------------------------------------------------------


c---- the same but with another mass zero (for mu or e-, light jets)

      Mfin(1)=amb
      Mfin(2)=0d0
      Mfin(3)=0d0
      MFin(4)=amneut(1)
      call SD_rstart(12,34,56,78)
      if(Etot.gt.(Mfin(1)+MFIN(2)+MFIN(3)+MFIN(4)))then
      do i=1, iiter
      call SD_rambo(4,etot,mfin,pfout,wt)

c---- take the momenta of Rambo 
      do j=1,4
      p1(j)=pfout(j,1)	!quark
      p2(j)=pfout(j,2)	!fermion
      p3(j)=pfout(j,3)	!antifermion
      p4(j)=pfout(j,4)	!neutralino
      pInt(j)=p2(j)+p3(j)+p4(j)	!intermediate particle (chargino)
      pfer(j)=p1(j)+p2(j)+p3(j)
      end do

! define scalar products in amplitude
      p1p2=SD_scal(p1,p2)
      p1p3=SD_scal(p1,p3)
      p1p4=SD_scal(p1,p4)
      p2p3=SD_scal(p2,p3)
      p2p4=SD_scal(p2,p4)
      p3p4=SD_scal(p3,p4)
      p1pI=SD_scal(p1,pInt)
      p2pI=SD_scal(p2,pInt)
      p3pI=SD_scal(p3,pInt)
      p4pI=SD_scal(p4,pInt)
      pI2=SD_scal(pInt,pInt)
      
      m1=Mfin(1)
      m3=0d0
      m4=amneut(1)
! W amplitude with chargino as intermediate particle and final state tau and bottom
      ampwchabmu=0d0
      ampwchabjets=0d0
c--- sum up different charginos as intermediate particles
      do j=1,2
      do k=1,2
      g1R1=DcharL(3,1,j)
      g1L1=EcharR(3,1,j)
      g1R2=DcharL(3,1,k)
      g1L2=EcharR(3,1,k)
      g3L=-g2ew/dSqrt(2d0)
      g2L1=CchaneutL(1,j)
      g2L2=CchaneutL(1,k)
      g2R1=CchaneutR(1,j)
      g2R2=CchaneutR(1,k)
    
      zwi=ampwchar(g1L1, g1L2, g2l1, g2l2, g1r1, 
     . g1r2, g2r1, g2r2, g3l, m1, m3, m4, 
     . amchar(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
      
    
c---- light jet results can be obtained by multiplying with colourfactor and ckmmatrixelements	
	ampwchabjets=ampwchabjets+zwi*(Vckm(1,1)**2+
     . Vckm(2,2)**2+vckm(1,2)**2+vckm(2,1)**2)*3d0/(pI2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/Abs(2d0*p2p3-amw**2+(0d0,1d0)*amw*wwidth)**2
c----- multiply with propagators
      ampwchabmu=ampwchabmu+zwi/(pI2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
      end do
      end do





c---- slepton contributions
      ampinterferenzwsmu=0d0
      ampinterferenzwselec=0d0
      ampslepbmu=0d0
      ampslepbelec=0d0
      do j=1,2
      do k=1,2
      do jj=1,6
      do kk=1,6
      g1R1=DcharL(3,1,j)
      g1L1=EcharR(3,1,j)
      g1R2=DcharL(3,1,k)
      g1L2=EcharR(3,1,k)
      g2r1=cchaneutslep(2,k,kk)
      g2r2=cchaneutslep(2,j,jj)
      g3l1=gneutelecl(2,1,kk)
      g3l2=gneutelecl(2,1,jj)
      g3r1=gneutelecr(2,1,kk)
      g3r2=gneutelecr(2,1,jj)
      g2l1=0d0
      g2l2=0d0

c------- selectron,... contribtuions
      ampslepbmu=ampslepbmu+ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)

      g2r1=cchaneutslep(1,k,kk)
      g2r2=cchaneutslep(1,j,jj)
      g3l1=gneutelecl(1,1,kk)
      g3l2=gneutelecl(1,1,jj)
      g3r1=gneutelecr(1,1,kk)
      g3r2=gneutelecr(1,1,jj)
       ampslepbelec=ampslepbelec+ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)

      end do
c---- interference slepton and sneutrino diagrams
       gXb2R=cchaneutslep(2,j,jj)
       gXc2L=-dchalepsneut(2,k)
       gXc2R=-echalepsneut(2,k)
       gbd3R=gneutneut(2,1)
       gcd3L=gneutelecl(2,1,jj)
       gcd3R=gneutelecr(2,1,jj)
       gXb2l=0d0

      ampslepbmu=ampslepbmu+ampintslepsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXb2l, gXb2R, gXc2L, gXc2R,gbd3R, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(2)**2)

       gXb2R=cchaneutslep(1,j,jj)
       gXc2L=-dchalepsneut(1,k)
       gXc2R=-echalepsneut(1,k)
       gbd3R=gneutneut(1,1)
       gcd3L=gneutelecl(1,1,jj)
       gcd3R=gneutelecr(1,1,jj)
       gXb2l=0d0

      ampslepbelec=ampslepbelec+ampintslepsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXb2l, gXb2R, gXc2L, gXc2R,gbd3R, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(1)**2)

c----- interference w slepton diagrams
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      gXb2R=cchaneutslep(2,j,jj)
      gcd3L=gneutelecl(2,1,jj)
      gcd3R=gneutelecr(2,1,jj)
      g3L=-g2ew/dSqrt(2d0)
      gxb2l=0d0 
      
      ampinterferenzwsmu=ampinterferenzwsmu+
     . 2d0*RealPart(ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2, gxb2l,gXb2R, g2l1, g2r1, g3L, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
c----- interference w slepton diagrams
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      gXb2R=cchaneutslep(1,j,jj)
      gcd3L=gneutelecl(1,1,jj)
      gcd3R=gneutelecr(1,1,jj)
      g3L=-g2ew/dSqrt(2d0)
      gxb2l=0d0 
      
      ampinterferenzwselec=ampinterferenzwselec+
     . 2d0*RealPart(ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2, gXb2l, gXb2R, g2l1, g2r1, g3L, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))      
      end do

c----- sneutrino diagram
      g3r1=gneutneut(2,1)
      g3r2=gneutneut(2,1)
      g3l1=0d0
      g3l2=0d0
      g2l1=-dchalepsneut(2,j)
      g2l2=-dchalepsneut(2,k)
      g2r1=-echalepsneut(2,j)
      g2r2=-echalepsneut(2,k)

      ampslepbmu=ampslepbmu+ampcharsneutrino(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2,g3l1,g3l2, g3r1, g3r2, 
     .  m1, m3, m4, amchar(j),
     .  amchar(k), p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     ./(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsneutrino(2)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(2)**2)

c----- sneutrino diagram
      g3r1=gneutneut(1,1)
      g3r2=gneutneut(1,1)
      g2l1=-dchalepsneut(1,j)
      g2l2=-dchalepsneut(1,k)
      g2r1=-echalepsneut(1,j)
      g2r2=-echalepsneut(1,k)

      ampslepbelec=ampslepbelec+ampcharsneutrino(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2,  g3l1, g3l2,g3r1, g3r2,
     .  m1, m3, m4, amchar(j),
     .  amchar(k), p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     ./(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsneutrino(1)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(1)**2)

c---- interference with W-contribution
      gXc2L=dchalepsneut(2,j)
      gXC2r=echalepsneut(2,j)
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      g3L=-g2ew/dSqrt(2d0)
      gbd3R=gneutneut(2,1)
     
      ampinterferenzwsmu=ampinterferenzwsmu+
     . 2d0*RealPart(ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, g2l1, g2r1, gXc2R, g3L,
     . gbd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(2)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
c---- interference with W-contribution
      gXc2L=dchalepsneut(1,j)
      gXC2r=echalepsneut(1,j)
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      g3L=-g2ew/dSqrt(2d0)
      gbd3R=gneutneut(1,1)
     
      ampinterferenzwselec=ampinterferenzwselec+
     . 2d0*RealPart(ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, g2l1, g2r1, gXc2R, g3L,
     . gbd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(3)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      end do
      end do


c------ squark contributions
      ampbsquarkjets=0d0
      do j=1,2
      do k=1,2
      do jj=1,6     
      do ifer=1,2
      do ifer2=1,2
      do kk=1,6
      
      g1R1=DcharL(3,1,j)
      g1L1=EcharR(3,1,j)
      g1R2=DcharL(3,1,k)
      g1L2=EcharR(3,1,k)
      g2l1=-DcharL(ifer,jj,j)
      g2l2=-DcharL(ifer,kk,k)
      g2r1=-EcharR(ifer,jj,j)
      g2r2=-EcharR(ifer,kk,k)
      g3r1=gneutuL(ifer2,1,jj)
      g3r2=gneutuL(ifer2,1,kk)
      g3l1=gneutur(ifer2,1,jj)
      g3l2=gneutur(ifer2,1,kk)



      ampbsquarkjets=ampbsquarkjets
     . +3d0*ampcharsup(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2,g3l1, g3l2, g3r1, g3r2, 
     . m1, 0d0, m4, amchar(j),
     . amchar(k), p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     ./(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p2p4+m4**2-amsupq(kk)**2)


      g2l1=fcharr(ifer, jj, j)
      g2l2=fcharr(ifer, kk, k)
      g2r1=ccharl(ifer, jj, j)
      g2r2=ccharl(ifer, kk, k)
      g3l1=gneutdL(ifer2,1,jj)
      g3l2=gneutdL(ifer2,1,kk)
      g3r1=gneutdr(ifer2,1,jj)
      g3r2=gneutdr(ifer2,1,kk)
	


      

      ampbsquarkjets=ampbsquarkjets
     . +3d0*ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)

      
    
      g2l1=-DcharL(ifer,jj,j)
      g2r1=-EcharR(ifer,jj,j)
      g3r1=gneutuL(ifer2,1,jj)
      g3l1=gneutur(ifer2,1,jj)
      g2l2=fcharr(ifer, kk, k)
      g2r2=ccharl(ifer, kk, k)
      g3l2=gneutdL(ifer2,1,kk)
      g3r2=gneutdr(ifer2,1,kk)


c----- squark exchange diagrams, same amplitude than slepton
c----- and sneutrino diagrams but couplings need to be changed
      ampbsquarkjets=ampbsquarkjets
     .+6d0*ampcharsupsdownint(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)     
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)
	

c---- interfernz squark-W-Terme
      enddo
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      g3L=-g2ew/dSqrt(2d0)*Vckm(ifer,ifer2)
      gXb2l=ccharl(ifer, jj, j)
      gXb2R=fcharr(ifer, jj, j)
      gcd3l=gneutdl(ifer2,1,jj)
      gcd3r=gneutdr(ifer2,1,jj)

      ampbsquarkjets=ampbsquarkjets+
     .6d0 *RealPart(ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2, gxb2l,gXb2R, g2l1, g2r1, g3L, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
 

      gXc2R=dcharl(ifer2, jj, j)
      gXC2L=echarr(ifer2, jj, j)
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
 
      g3L=-g2ew/dSqrt(2d0)*Vckm(ifer,ifer2)
      gbd3R=gneutul(ifer, 1,jj)


      ampbsquarkjets=ampbsquarkjets+
     .6d0*RealPart(ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, g2l1, g2r1, gXc2R, g3L,
     . gbd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      
      enddo
      enddo
      enddo
      enddo
      enddo



c----- contribtutions from diagram invlolving a SM fermion
      ampSMfermionmu=0d0
      ampSMfermionel=0d0
c----- SM fermion matrix elemnt squared
       p1pfer=SD_scal(p1,pfer)
      p2pfer=SD_scal(p2,pfer)
      p3pfer=SD_scal(p3,pfer)
      p4pfer=SD_scal(p4,pfer)
      pfer2=SD_scal(pfer,pfer)
      pIpfer=SD_scal(pInt, pfer)


      zwi=0d0
      zwi2=0d0
      do j=1,3
      do k=1,3
      g1L1=-g2ew/dsqrt(2d0)*Vckm(j,3)
      g1L2=-g2ew/dsqrt(2d0)*Vckm(k,3)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2l2=gneutur(k,1,1)
      g2r2=gneutul(k,1,1)
      g3l=-g2ew/dsqrt(2d0)

      zwi=zwi+ampSMfer(g1L1, g1L2, g2l1, g2l2, 
     . g2r1, g2r2,g3l,
     . m1, m3, m4, mfu(j), mfu(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2)/(pfer2-mfu(j)**2)
     . /(pfer2-mfu(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
      
      
      enddo
      end do
c---- interference term with W-chargino diagram
      ampSMfermionel=zwi
      ampSMfermionmu=zwi

      do k=1,2
      do j=1,3
      g1r2=DcharL(3,1,k)
      g1l2=EcharR(3,1,k)
      g2L2=CchaneutL(1,k)
      g2R2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)
      g1L1=-g2ew/dsqrt(2d0)*Vckm(j,3)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      ampSMfermionel=ampSMfermionel+
     . 2d0*ampintWSMferWchar(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l,  m1, m3, m4, 
     .mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2 


	ampSMfermionmu=ampSmfermionel
c---- interfernce with slpeton contribtuions

      do kk=1,6
      g2l2=0d0
      g2r2=cchaneutslep(1,k,kk)
      g3l1=-g2ew/dsqrt(2d0)
      g3r2=gneutelecr(1,1,kk)
      g3l2=gneutelecl(1,1,kk)
      ampSMfermionel=ampSMfermionel+2d0*RealPart(
     . ampintWSMferslep(g1L1, g1L2, g1r2,
     .  g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,m1, m3,
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2))
      
      g2r2=cchaneutslep(2,k,kk)
      g3l1=-g2ew/dsqrt(2d0)
      g3r2=gneutelecr(2,1,kk)
      g3l2=gneutelecl(2,1,kk)
      ampSMfermionmu=ampSMfermionmu+2d0*RealPart(
     . ampintWSMferslep(g1L1, g1L2, g1r2,
     .  g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,m1, m3,
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2))
      enddo


      g2l2=-dchalepsneut(1,k)
      g2r2=-echalepsneut(1,k)
      g3r2=0d0
      g3l2=gneutneut(1,1)
      ampSMfermionel=ampSMfermionel+2d0*RealPart(
     . ampintWSMfersneut(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2, m1, m3, 
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2))
      

      g2l2=-dchalepsneut(2,k)
      g2r2=-echalepsneut(2,k)
      g3r2=0d0
      g3l2=gneutneut(2,1)
      ampSMfermionmu=ampSMfermionmu+2d0*RealPart(
     . ampintWSMfersneut(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2, m1, m3, 
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2))
      enddo
      enddo
      


c----- contribtutions from diagram invlolving a SM fermion
      ampSMfermionjet=0d0
c----- SM fermion matrix elemnt squared

      do j=1,3
      do k=1,3
      
      g1L1=-g2ew/dsqrt(2d0)*(Vckm(j,3))
      g1L2=-g2ew/dsqrt(2d0)*Vckm(k,3)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2l2=gneutur(k,1,1)
      g2r2=gneutul(k,1,1)
      g3l=-g2ew/dsqrt(2d0)


      ampSMfermionjet=ampSMfermionjet+3d0*ampSMfer(g1L1, g1L2, 
     . g2l1, g2l2, g2r1, g2r2,g3l,
     . m1, m3, m4, mfu(j), mfu(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2)/(pfer2-mfu(j)**2)
     . /(pfer2-mfu(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
      
      end do
      enddo
c---- interference term with W-chargino diagram
      
      do j=1,3
      do k=1,2
       g1L1=-g2ew/dsqrt(2d0)*(Vckm(j,3))
      g1r2=DcharL(3,1,k)
      g1l2=EcharR(3,1,k)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2L2=CchaneutL(1,k)
      g2R2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)
      ampSMfermionjet=ampSMfermionjet+6d0*ampintWSMferWchar(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l,  m1, m3, m4, 
     .mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2 
	
      enddo
      enddo
      ampSMfermionjet=ampSMfermionjet*(Vckm(1,1)**2+VCKM(2,1)**2
     . +VCKM(1,2)**2+VCKM(2,2)**2)
c---- interfernce with sdown contribtuions

      do j=1,3
      do k=1,2
      do kk=1,6
      do ifer=1,2
      do ifer2=1,2
      g1L1=-g2ew/dsqrt(2d0)*(Vckm(j,3))
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2l2=ccharl(ifer, kk, k)
      g2r2=ccharl(ifer, kk, k)
      g3l1=-g2ew/dsqrt(2d0)*(Vckm(1,1)+VCKM(2,1)
     . +VCKM(1,2)+VCKM(2,2))
      g3r2=gneutdR(ifer2,1,kk)
      g3l2=gneutdL(ifer2,1,kk)
      ampSMfermionjet=ampSMfermionjet+6d0*RealPart(
     . ampintWSMferslep(g1L1, g1L2, g1r2,
     .  g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,m1, m3,
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2))

      
      g2l2=Dcharl(ifer2,kk,k)
      g2r2=EcharR(ifer2,kk,k)
      g3r2=gneutur(ifer,1,kk)
      g3l2=gneutul(ifer,1,kk)
      ampSMfermionjet=ampSMfermionjet+6d0*RealPart(
     . ampintWSMfersneut(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2, m1, m3, 
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p2p4+m4**2-amsupq(kk)**2))
      enddo
      end do
      enddo
      end do
      enddo
c---- diagram with W squark propagator
      ampsquarkpropW=0d0

      do j=1,6
      do k=1,6
      g1l1=gsqsqW(j)
      g1l2=gsqsqW(k)
      g2l1=gneutdr(3,1,j)
      g2l2=gneutdr(3,1,k)
      g2r1=gneutdl(3,1,j)
      g2r2=gneutdl(3,1,k)
      g3l=-g2ew/dsqrt(2d0)
      

      ampsquarkpropW=ampsquarkpropW+ampvecscalar(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,g3l,0d0,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, amw)/
     . Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
      enddo
      enddo


c----- interference W squark- W chargino diagrams
      do j=1,6
      do k=1,2
      g1l1=gsqsqW(j)
      g1l2=EcharR(3,1,j)
      g1R2=DcharL(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=CchaneutL(1,k)
      g2r1=gneutdl(3,1,j)
      g2r2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)

      ampsquarkpropW=ampsquarkpropW+2d0*RealPart(
     . ampintWcharWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pi2-amchar(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2)
      enddo
      enddo
c------ interfernce W squark - W SM Fermion
      do j=1,6
      do k=1,3
      g1l1=gsqsqW(j)
      g1l2=gneutur(k,1,1)
      g1R2=gneutul(k,1,1)
      g2l1=gneutdr(3,1,j)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(k,3)
      g2r1=gneutdl(3,1,j)
      g3l=-g2ew/dsqrt(2d0)
      ampsquarkpropW=ampsquarkpropW+2d0*RealPart(
     . ampintWferWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pfer2-mfu(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2) 
      enddo
      enddo

      ampsquarkpropWjet=0d0
c--------interfernce W squark Sdown diagram
      do j=1,6
      do k=1,2
      do kk=1,6
      do ifer=1,2
      do ifer2=1,2
      g1l1=gsqsqW(j)
      g1l2=EcharR(3,1,k)
      g1R2=DcharL(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=fcharr(ifer, kk, k)
      g2r1=gneutdl(3,1,j)
      g2r2=ccharl(ifer, kk, k)
      g3l1=-g2ew/dsqrt(2d0)*(VCKM(ifer,ifer2))
      g3l2=gneutdl(ifer2,1,kk)
      g3r2=gneutdr(ifer2,1,kk)
          
   
      ampsquarkpropWjet=ampsquarkpropWjet+6d0*RealPart(
     .ampintcharslepWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     ./(pi2-amchar(k)**2)/(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth) )

      
      
c---- interference W squark- chargino sup diagram      

      g1R2=DcharL(3,1,k)
      g1L2=EcharR(3,1,k)
      g2l2=DcharL(ifer2,kk,k)
      g2r2=EcharR(ifer2,kk,k)
      g3l2=gneutuL(ifer,1,kk)
      g3r2=gneutur(ifer,1,kk)

      ampsquarkpropWjet=ampsquarkpropWjet+6d0*
     . RealPart(ampintcharsneutWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     . /(pi2-amchar(k)**2)/(2d0*p2p4+m4**2-amsupq(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo
      enddo
      enddo


c--------interfernce W squark Slepton diagram

      ampsquarkpropWmu=0d0
      ampsquarkpropWe=0d0
      do j=1,6
      do k=1,2
      do kk=1,6
      g1l1=gsqsqW(j)
      g1l2=EcharR(3,1,k)
      g1R2=DcharL(3,1,k)
      g2l1=gneutdr(3,1,j)
      g2l2=0d0
      g2r1=gneutdl(2,1,j)
      g2r2=cchaneutslep(3,k,kk)
      g3l1=-g2ew/dsqrt(2d0)
      g3l2=gneutelecl(2,1,kk)
      g3r2=gneutelecr(2,1,kk)

          
   
      ampsquarkpropWmu=ampsquarkpropWmu+2d0*RealPart(
     .ampintcharslepWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     ./(pi2-amchar(k)**2)/(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth) )

      g2r1=gneutdl(1,1,j)
      g3l2=gneutelecl(1,1,kk)
      g3r2=gneutelecr(1,1,kk)
      ampsquarkpropWe=ampsquarkpropWe+2d0*RealPart(
     .ampintcharslepWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     ./(pi2-amchar(k)**2)/(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth) )
      enddo
      
c---- interference W squark- chargino sneutrino diagram      
      g2l2=dchalepsneut(2,k)
      g2r2=echalepsneut(2,k)
      g3l2=gneutneut(2,1)
      g3r2=0d0



      ampsquarkpropWmu=ampsquarkpropWmu+2d0*
     . RealPart(ampintcharsneutWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     . /(pi2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(3)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      g2l2=dchalepsneut(1,k)
      g2r2=echalepsneut(1,k)
      g3l2=gneutneut(1,1)
      g3r2=0d0



      ampsquarkpropWe=ampsquarkpropWe+2d0*
     . RealPart(ampintcharsneutWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     . /(pi2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(3)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo

c---- pay attention here, in sigbmu also e- part of ampSmfermion inside...
      sigbmu=sigbmu+wt*ampwchabmu/iiter
     .+wt*ampinterferenzwsmu/iiter+wt*ampslepbmu/iiter
     . +wt*ampSMfermionmu/iiter+wt*ampsquarkpropW/iiter
     .+wt*ampsquarkpropWmu/iiter


      sigbelec=sigbelec+wt*ampwchabmu/iiter
     . +wt*ampinterferenzwselec/iiter+wt*ampslepbelec/iiter
     . +wt*ampSMfermionel/iiter+wt*ampsquarkpropW/iiter
     . +wt*ampsquarkpropWe/iiter


      sigbjets=sigbjets+wt*ampwchabjets/iiter+wt*ampbsquarkjets/iiter
     . +wt*ampSMfermionjet/iiter+wt*ampsquarkpropW*
     .(VCKM(1,2)**2+VCKM(1,1)**2+VCKM(2,1)**2+VCKM(2,2)**2)*
     . 3d0/iiter+wt*ampsquarkpropWjet/iiter

      enddo
      else
      print*, "decay stop---> 
     . neutralino b + mass less fermions forbidden"
      endif
!        Go To 13
C-------------------------------------------------------
c------- FOURTH CASE: INITIAL MASSLESS JET, FINAL TAU
C--------------------------------------------------------
14      Etot=amsupq(1)
      Mfin(1)=0d0
      Mfin(2)=0d0
      Mfin(3)=amtau
      MFin(4)=amneut(1)
      if(Etot.gt.(amtau+amneut(1)))then
      do i=1, iiter
      call SD_rambo(4,etot,mfin,pfout,wt)

c---- take the momenta of Rambo 
      do j=1,4
      p1(j)=pfout(j,1)	!quark
      p2(j)=pfout(j,2)	!fermion
      p3(j)=pfout(j,3)	!antifermion
      p4(j)=pfout(j,4)	!neutralino
      pInt(j)=p2(j)+p3(j)+p4(j)	!intermediate particle (chargino)
      pfer(j)=p1(j)+p2(j)+p3(j) !intemdeiate particle(fermion such as top)
      end do

! define scalar products in amplitude
      p1p2=SD_scal(p1,p2)
      p1p3=SD_scal(p1,p3)
      p1p4=SD_scal(p1,p4)
      p2p3=SD_scal(p2,p3)
      p2p4=SD_scal(p2,p4)
      p3p4=SD_scal(p3,p4)
      p1pI=SD_scal(p1,pInt)
      p2pI=SD_scal(p2,pInt)
      p3pI=SD_scal(p3,pInt)
      p4pI=SD_scal(p4,pInt)
      pI2=SD_scal(pInt,pInt)
      
      m1=Mfin(1)
      m3=MFIN(3)
      m4=amneut(1)
! W amplitude with chargino as intermediate particle and final state tau and bottom
      ampwchajettau=0d0
c--- sum up different charginos as intermediate particles
      do j=1,2
      do k=1,2
      g1R1=DcharL(2,1,j)
      g1L1=EcharR(2,1,j)
      g1R2=DcharL(2,1,k)
      g1L2=EcharR(2,1,k)
      g3L=-g2ew/dSqrt(2d0)
      g2L1=CchaneutL(1,j)
      g2L2=CchaneutL(1,k)
      g2R1=CchaneutR(1,j)
      g2R2=CchaneutR(1,k) 
     
       ampwchajettau=ampwchajettau
     . +ampwchar(g1L1, g1L2, g2l1, g2l2, g1r1, 
     . g1r2, g2r1, g2r2, g3l, m1, m3, m4, 
     . amchar(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(pI2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/
     . Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
c---- also sum the up-contribtuion
      g1R1=DcharL(1,1,j)
      g1L1=EcharR(1,1,j)
      g1R2=DcharL(1,1,k)
      g1L2=EcharR(1,1,k)
	ampwchajettau=ampwchajettau+ampwchar(g1L1, g1L2, g2l1, g2l2, g1r1, 
     . g1r2, g2r1, g2r2, g3l, m1, m3, m4, 
     . amchar(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(pI2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/
     . Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
      end do
      end do



c---- slepton contributions
      ampinterferenzwjetslep=0d0
      ampslepjettau=0d0
      do j=1,2
      do k=1,2
      do jj=1,6
      do kk=1,6
      g1R1=DcharL(2,1,j)
      g1L1=EcharR(2,1,j)
      g1R2=DcharL(2,1,k)
      g1L2=EcharR(2,1,k)
      g2r1=cchaneutslep(3,k,kk)
      g2r2=cchaneutslep(3,j,jj)
      g3l1=gneutelecl(3,1,kk)
      g3l2=gneutelecl(3,1,jj)
      g3r1=gneutelecr(3,1,kk)
      g3r2=gneutelecr(3,1,jj)
      g2l1=0d0 
      g2l2=0d0

c------- selectron,... contribtuions
      ampslepjettau=ampslepjettau+ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2,g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)
      g1R1=DcharL(1,1,j)
      g1L1=EcharR(1,1,j)
      g1R2=DcharL(1,1,k)
      g1L2=EcharR(1,1,k)
       ampslepjettau=ampslepjettau+ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)
      
      end do
c---- interference slepton and sneutrino diagrams
       g1R1=DcharL(1,1,j)
       g1L1=EcharR(1,1,j)
       g1R2=DcharL(1,1,k)
       g1L2=EcharR(1,1,k)
       gXb2R=cchaneutslep(3,j,jj)
       gXc2L=-dchalepsneut(3,k)
       gXc2R=-echalepsneut(3,k)
       gbd3R=gneutneut(3,1)
       gcd3L=gneutelecl(3,1,jj)
       gcd3R=gneutelecr(3,1,jj)
       gXb2l=0d0

      ampslepjettau=ampslepjettau+ampintslepsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXb2l,gXb2R, gXc2L, gXc2R,gbd3R, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2)

      g1R1=DcharL(2,1,j)
      g1L1=EcharR(2,1,j)
      g1R2=DcharL(2,1,k)
      g1L2=EcharR(2,1,k)
       ampslepjettau=ampslepjettau+ampintslepsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXb2l,gXb2R, gXc2L, gXc2R,gbd3R, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2)

c----- interference w slepton diagrams
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      gXb2R=cchaneutslep(3,j,jj)
      gcd3L=gneutelecl(3,1,jj)
      gcd3R=gneutelecr(3,1,jj)
      g3L=-g2ew/dSqrt(2d0)
       gxb2l=0d0
      
      ampinterferenzwjetslep=ampinterferenzwjetslep+
     . 2d0*RealPart(ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2, gXb2l,gXb2R, g2l1, g2r1, g3L, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      
      g1R1=DcharL(1,1,j)
      g1L1=EcharR(1,1,j)
      g1R2=DcharL(1,1,k)
      g1L2=EcharR(1,1,k)
      ampinterferenzwjetslep=ampinterferenzwjetslep+
     . 2d0*RealPart(ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2, gXb2l,gXb2R, g2l1, g2r1, g3L, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      
      end do
c----- sneutrino diagram
      g1R1=DcharL(1,1,j)
      g1L1=EcharR(1,1,j)
      g1R2=DcharL(1,1,k)
      g1L2=EcharR(1,1,k)
      g3r1=gneutneut(3,1)
      g3r2=gneutneut(3,1)
      g3l1=0d0
      g3l2=0d0
      g2l1=-dchalepsneut(3,j)
      g2l2=-dchalepsneut(3,k)
      g2r1=-echalepsneut(3,j)
      g2r2=-echalepsneut(3,k)

      ampslepjettau=ampslepjettau+ampcharsneutrino(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2, g3l1,
     . g3l2,g3r1, g3r2,  m1, m3, m4, amchar(j),
     .  amchar(k), p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     ./(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsneutrino(3)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2)
      g1R1=DcharL(2,1,j)
      g1L1=EcharR(2,1,j)
      g1R2=DcharL(2,1,k)
      g1L2=EcharR(2,1,k)
       ampslepjettau=ampslepjettau+ampcharsneutrino(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2, g3l1, g3l2,g3r1, g3r2, 
     . m1, m3, m4, amchar(j),
     .  amchar(k), p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     ./(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsneutrino(3)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2)

c---- interference with W-contribution
      gXc2L=dchalepsneut(3,j)
      gXC2r=echalepsneut(3,j)
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      g3L=-g2ew/dSqrt(2d0)
      gbd3R=gneutneut(3,1)
     
      ampinterferenzwjetslep=ampinterferenzwjetslep+
     . 2d0*RealPart(ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, g2l1, g2r1, gXc2R, g3L,
     . gbd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(3)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

            
      g1R1=DcharL(1,1,j)
      g1L1=EcharR(1,1,j)
      g1R2=DcharL(1,1,k)
      g1L2=EcharR(1,1,k)

      ampinterferenzwjetslep=ampinterferenzwjetslep+
     . 2d0*RealPart(ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, g2l1, g2r1, gXc2R, g3L,
     . gbd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(3)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      end do
      end do

c----- contribtutions from diagram invlolving a SM fermion
      ampSMfermion=0d0
c----- SM fermion matrix elemnt squared
       p1pfer=SD_scal(p1,pfer)
      p2pfer=SD_scal(p2,pfer)
      p3pfer=SD_scal(p3,pfer)
      p4pfer=SD_scal(p4,pfer)
      pfer2=SD_scal(pfer,pfer)
      pIpfer=SD_scal(pInt, pfer)



      do j=1,3
      do k=1,3
      do jj=1,2
      g1L1=-g2ew/dsqrt(2d0)*Vckm(j,jj)
      g1L2=-g2ew/dsqrt(2d0)*Vckm(k,jj)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2l2=gneutur(k,1,1)
      g2r2=gneutul(k,1,1)
      g3l=-g2ew/dsqrt(2d0)


      ampSMfermion=ampSMfermion+ampSMfer(g1L1, g1L2, g2l1, g2l2, 
     . g2r1, g2r2,g3l,
     . m1, m3, m4, mfu(j), mfu(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2)/(pfer2-mfu(j)**2)
     . /(pfer2-mfu(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
      
      enddo
      end do
      end do
c---- interference term with W-chargino diagram

      do ifer=1,2
      do k=1,2
      do j=1,3
      g1r2=DcharL(ifer,1,k)
      g1l2=EcharR(ifer,1,k)
      g2L2=CchaneutL(1,k)
      g2R2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)
      g1L1=-g2ew/dsqrt(2d0)*Vckm(j,ifer)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      ampSMfermion=ampSMfermion+2d0*ampintWSMferWchar(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l,  m1, m3, m4, 
     .mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2 

c---- interfernce with slpeton contribtuions


      do kk=1,6
      g2l2=0d0
      g2r2=cchaneutslep(3,k,kk)
      g3l1=-g2ew/dsqrt(2d0)
      g3r2=gneutelecr(3,1,kk)
      g3l2=gneutelecl(3,1,kk)
      ampSMfermion=ampSMfermion+2d0*RealPart(
     . ampintWSMferslep(g1L1, g1L2, g1r2,
     .  g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,m1, m3,
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2))


      enddo
      g2l2=-dchalepsneut(3,k)
      g2r2=-echalepsneut(3,k)
      g3r2=0d0
      g3l2=gneutneut(3,1)
      ampSMfermion=ampSMfermion+2d0*RealPart(
     . ampintWSMfersneut(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2, m1, m3, 
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2))
      enddo
      enddo
      end do


c------- the diagram with W and squark propagator
      ampsquarkpropW=0d0
      ampchahi=0d0
      do ifer=1,2
      do j=1,6
      do k=1,6
      g1l1=gsqsqW(j)
      g1l2=gsqsqW(k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=gneutdr(ifer,1,k)
      g2r1=gneutdl(ifer,1,j)
      g2r2=gneutdl(ifer,1,k)
      g3l=-g2ew/dsqrt(2d0)
      

      ampsquarkpropW=ampsquarkpropW+ampvecscalar(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,g3l,0d0,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, amw)/
     . Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)

      g1l1=gsqsqchaHi(j)
      g1l2=gsqsqchaHi(k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      
      
c----- digram charged Higgs squark propagator
      ampchahi=ampchahi+amphiggssquark(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,0d0,g3r1,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4)/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/
     . (2d0*p2p3+m3**2-amch**2)**2

c---- interference Higgs-squark with W squark diagram
      g1l1=gsqsqW(j)
      g1l2=gsqsqchaHi(k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=gneutdr(ifer,1,k)
      g2r1=gneutdl(ifer,1,j)
      g2r2=gneutdl(ifer,1,k)
      g3l1=-g2ew/dsqrt(2d0)
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta

      ampchahi=ampchahi+2d0*RealPart(
     .ampinthiggssquarkWsquark(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,g3l1,0d0, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4)/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/
     . (2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo

c----- interference charged Higgs squark- W chargino diagram
      do j=1,6
      do k=1,2
      do ifer=1,2
      g1l1=gsqsqchaHi(j)
      g1l2=EcharR(ifer,1,j)
      g1r2=DcharL(ifer,1,j)
      g2l1=gneutdr(3,1,j)
      g2l2=CchaneutL(1,k)
      g2r1=gneutdl(3,1,j)
      g2r2=CchaneutR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/dsqrt(2d0)

      ampchahi=ampchahi+2d0*RealPart(ampintchaHisquarkWchar(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo

c------ interference charged Higgs squark -W SM fermion diagram
      do j=1,6
      do k=1,3
      do ifer=1,2
      g1l1=gsqsqchaHi(j)
      g1l2=gneutur(k,1,1)
      g1r2=gneutul(k,1,1)
      g2l1=gneutdr(ifer,1,j)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(k,ifer)
      g2r1=gneutdl(ifer,1,j)
      g2r2=0d0
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/dsqrt(2d0)

      ampchahi=ampchahi+2d0*RealPart(ampintchaHisquarkWSMfer(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(k))/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pfer2-mfu(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
      enddo
      enddo
c---- interference charged Higgs squark- slepton chargino diagram
      do ifer=1,2
      do j=1,6
      do k=1,2
      do kk=1,6
      g1l1=gsqsqchaHi(j)
      g1l2=EcharR(ifer,1,k)
      g1r2=DcharL(ifer,1,k)
      g2l1=gneutdr(ifer,1,j)
      g2r1=gneutdl(ifer,1,j)
      g2r2=cchaneutslep(3,k,kk)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=gneutelecl(3,1,kk)
      g3r2=gneutelecr(3,1,kk)


      ampchahi=ampchahi
     . +2d0*RealPart(ampintchaHisquarkslepton(g1l1, g1l2, 
     . g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/(pi2-amchar(k)**2)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)
     ./(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      g2l2=dchalepsneut(3,k)
      g2r2=echalepsneut(3,k)
      g3r2=gneutneut(3,1)


            ampchahi=ampchahi
     . +2d0*RealPart(ampintchaHisquarksneutrino(g1l1, g1l2, 
     . g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/(pi2-amchar(k)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(3)**2)
     ./(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo



c----- interference W squark- W chargino diagrams
      do j=1,6
      do k=1,2
      do ifer=1,2
      g1l1=gsqsqW(j)
      g1l2=EcharR(ifer,1,j)
      g1R2=DcharL(ifer,1,k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=CchaneutL(1,k)
      g2r1=gneutdl(ifer,1,j)
      g2r2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)

      ampsquarkpropW=ampsquarkpropW+2d0*RealPart(
     . ampintWcharWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pi2-amchar(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2)
      enddo
      enddo
      enddo
c------ interfernce W squark - W SM Fermion
      do j=1,6
      do k=1,3
      do ifer=1,2
      g1l1=gsqsqW(j)
      g1l2=gneutur(k,1,1)
      g1R2=gneutul(k,1,1)
      g2l1=gneutdr(ifer,1,j)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(k,ifer)
      g2r1=gneutdl(ifer,1,j)
      g3l=-g2ew/dsqrt(2d0)
      ampsquarkpropW=ampsquarkpropW+2d0*RealPart(
     . ampintWferWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pfer2-mfu(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2) 

      enddo
      enddo
      enddo
c--------interfernce W squark Slepton diagram
      do j=1,6
      do k=1,2
      do ifer=1,2
      do kk=1,6
      
      g1l1=gsqsqW(j)
      g1l2=EcharR(ifer,1,k)
      g1R2=DcharL(ifer,1,k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=0d0
      g2r1=gneutdl(ifer,1,j)
      g2r2=cchaneutslep(3,k,kk)
      g3l1=-g2ew/dsqrt(2d0)
      g3l2=gneutelecl(3,1,kk)
      g3r2=gneutelecr(3,1,kk)

          
   
      ampsquarkpropW=ampsquarkpropW+2d0*RealPart(
     .ampintcharslepWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     ./(pi2-amchar(k)**2)/(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth) )

      enddo
      
c---- interference W squark- chargino sneutrino diagram      
      g2l2=dchalepsneut(3,k)
      g2r2=echalepsneut(3,k)
      g3l2=gneutneut(3,1)
      g3r2=0d0

      ampsquarkpropW=ampsquarkpropW+2d0*
     . RealPart(ampintcharsneutWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     . /(pi2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(3)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo



c----- charged Higgs diagram: chargino
      do ifer=1,2 
      do j=1,2
      do k=1,2
      g1l1=EcharR(ifer,1,j)
      g1r1=DcharL(ifer,1,j)
      g1l2=EcharR(ifer,1,k)
      g1r2=DcharL(ifer,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=-g2ew*QchncL(1,k)
      g2r2=-g2ew*QchncR(1,k)
      g3r1= g2ew/dsqrt(2d0)/amw*mfe(3)*tanbeta

      ampchahi=ampchahi+ampchaHichar(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(j), amchar(k))/
     . (pi2-amchar(j)**2)/(pI2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)**2
      enddo
      enddo
      enddo


c----- interference charged Higgs chargino and charged Higgs squark diagrams
      do ifer=1,2
      do j=1,2
      do k=1,6
      g1l1=EcharR(ifer,1,j)
      g1r1=Dcharl(ifer,1,j)
      g1l2=gsqsqchaHi(k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=gneutdr(ifer,1,k)
      g2r2=gneutdl(ifer,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta

      ampchahi=ampchahi+2d0*ampintchaHichargchaHisquark(g1l2,
     .  g1l1,  g1r1,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j))/(pI2-amchar(j)**2)
     . /(2d0*p2p3+m3**2-amch**2)**2
     . /(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)  
      enddo
      enddo
      enddo
c------ interference charged Higgs chargino - W squark diagrams
       do ifer=1,2
       do j=1,2
       do k=1,6

      g1l1=EcharR(ifer,1,j)
      g1r1=Dcharl(ifer,1,j)
      g1l2=gsqsqW(k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=gneutdr(ifer,1,k)
      g2r2=gneutdl(ifer,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))
      ampchahi=ampchahi+2d0*RealPart(ampintchaHichargWsquark(g1l2,
     .  g1l1,  g1r1,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j))/(pI2-amchar(j)**2)
     ./(2d0*p2p3+m3**2-amch**2)/(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)/
     . (2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
       enddo
       enddo
      enddo
c----- interfernce charged Higgs chargino-W chargino
      do ifer=1,2
      do j=1,2
      do k=1,2
      g1l1=EcharR(ifer,1,j)
      g1r1=Dcharl(ifer,1,j)
      g1l2=EcharR(ifer,1,k)
      g1r2=Dcharl(ifer,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=CchaneutL(1,k)
      g2r2=CchaneutR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))

      ampchahi=ampchahi+2d0*RealPart(ampintchaHiWchar(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))
     . /(pi2-amchar(j)**2)/(pi2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
      enddo
      enddo
c----- interfernce charged Higgs chargino-W Sm fermion diagrams
      do ifer=1,2
      do j=1,2
      do k=1,3
      g1l1=EcharR(ifer,1,j)
      g1r1=Dcharl(ifer,1,j)
      g1l2=gneutur(k,1,1)
      g1r2=gneutul(k,1,1)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=-g2ew/dsqrt(2d0)*Vckm(ifer,k)
      g2r2=-g2ew/dsqrt(2d0)*Vckm(ifer,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))

      ampchahi=ampchahi-2d0*RealPart(ampintchaHiWSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, pIpFer, mfu(k), amchar(j))
     . /(pfer2-mfu(k)**2)/(pi2-amchar(j)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo
c----- interference charged Higgs chargino+ slepton diagram
      do ifer=1,2
      do j=1,2
      do k=1,2
      do kk=1,6
      g1l1=EcharR(ifer,1,j)
      g1r1=Dcharl(ifer,1,j)
      g1l2=EcharR(ifer,1,k)
      g1r2=Dcharl(ifer,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2r2=cchaneutslep(3,k,kk)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=gneutelecl(3,1,kk)
      g3r2=gneutelecr(3,1,kk)


      ampchahi=ampchahi+2d0*ampintchaHichargslepton(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))/
     .(pi2-amchar(j)**2)/(pi2-amchar(k)**2)/(2d0*p2p3-amch**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)


      enddo
      enddo
      enddo
      enddo
c---- interference charged Higgs chargino+sneutrino diagrams
      do ifer=1,2
      do j=1,2
      do k=1,2
      g1l1=EcharR(ifer,1,j)
      g1r1=Dcharl(ifer,1,j)
      g1l2=EcharR(ifer,1,k)
      g1r2=Dcharl(ifer,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=dchalepsneut(3,k)
      g2r2=echalepsneut(3,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3r2=gneutneut(3,1)

      ampchahi=ampchahi+2d0*ampintchaHichargsneutrino(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))/
     .(pi2-amchar(j)**2)/(pi2-amchar(k)**2)/(2d0*p2p3-amch**2)
     ./(2d0*p2p4+m4**3-amsneutrino(3)**2)
      enddo
      enddo
      enddo


      sigjettau=sigjettau+wt*ampwchajettau/iiter
     .+wt*ampinterferenzwjetslep/iiter+ampslepjettau*wt/iiter
     .+wt*ampSMfermion/iiter+wt*ampchahi/iiter+wt*ampsquarkpropW/iiter
     


      end do
      else
      sigjettau=0d0
      print*, "decay stop---> neutralino
     . tau + massless fermions forbidden"
      endif

c-------------------------------------------------------------
c----- FIFTH CASE: INITIAL: MASSLESS QUARK, FINAL: B_QUARK
C---------------------------------------------------------
      Etot=amsupq(1)
      Mfin(1)=0d0
      Mfin(2)=0d0
      Mfin(3)=amb
      MFin(4)=amneut(1)

      if(Etot.gt.(amb+amneut(1)))then
      do i=1, iiter
      call SD_rambo(4,etot,mfin,pfout,wt)

c---- take the momenta of Rambo 
      do j=1,4
      p1(j)=pfout(j,1)	!quark
      p2(j)=pfout(j,2)	!fermion
      p3(j)=pfout(j,3)	!antifermion
      p4(j)=pfout(j,4)	!neutralino
      pInt(j)=p2(j)+p3(j)+p4(j)	!intermediate particle (chargino)
      pfer(j)=p1(j)+p2(j)+p3(j)
      end do

! define scalar products in amplitude
      p1p2=SD_scal(p1,p2)
      p1p3=SD_scal(p1,p3)
      p1p4=SD_scal(p1,p4)
      p2p3=SD_scal(p2,p3)
      p2p4=SD_scal(p2,p4)
      p3p4=SD_scal(p3,p4)
      p1pI=SD_scal(p1,pInt)
      p2pI=SD_scal(p2,pInt)
      p3pI=SD_scal(p3,pInt)
      p4pI=SD_scal(p4,pInt)
      pI2=SD_scal(pInt,pInt)
      
      m1=Mfin(1)
      m3=amb
      m4=amneut(1)
! W amplitude with chargino as intermediate particle and final state tau and bottom
      ampwchajetb=0d0
c--- sum up different charginos as intermediate particles
      do j=1,2
      do k=1,2
      g1R1=DcharL(2,1,j)
      g1L1=EcharR(2,1,j)
      g1R2=DcharL(2,1,k)
      g1L2=EcharR(2,1,k)
      g3L=-g2ew/dSqrt(2d0)
      g2L1=CchaneutL(1,j)
      g2L2=CchaneutL(1,k)
      g2R1=CchaneutR(1,j)
      g2R2=CchaneutR(1,k)      

      zwi=ampwchar(g1L1, g1L2, g2l1, g2l2, g1r1, 
     . g1r2, g2r1, g2r2, g3l, m1, m3, m4, 
     . amchar(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)

c--- propagators, colour factor and summing up bc and bu contribution
       ampwchajetb=ampwchajetb+zwi*3d0*(Vckm(1,3)**2+Vckm(2,3)**2)
     . /(pI2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
c----- the same but up as initail state

      g1R1=DcharL(1,1,j)
      g1L1=EcharR(1,1,j)
      g1R2=DcharL(1,1,k)
      g1L2=EcharR(1,1,k)
      zwi=ampwchar(g1L1, g1L2, g2l1, g2l2, g1r1, 
     . g1r2, g2r1, g2r2, g3l, m1, m3, m4, 
     . amchar(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)

c--- propagators, colour factor and summing up bc and bu contribution
       ampwchajetb=ampwchajetb+zwi*3d0*(Vckm(1,3)**2+Vckm(2,3)**2)
     . /(pI2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
      end do
      end do
c------ squark contributions
      ampjetssquarkb=0d0
      do ifer2=1,2
      do j=1,2
      do k=1,2
      do jj=1,6
      do ifer=1,2
      do kk=1,6
      
      g1R1=DcharL(ifer2,1,j)
      g1L1=EcharR(ifer2,1,j)
      g1R2=DcharL(ifer2,1,k)
      g1L2=EcharR(ifer2,1,k)
      g2l1=-DcharL(3,jj,j)
      g2l2=-DcharL(3,kk,k)
      g2r1=-EcharR(3,jj,j)
      g2r2=-EcharR(3,kk,k)
      g3l1=gneutur(ifer,1,jj)
      g3l2=gneutur(ifer,1,kk)
      g3r1=gneutul(ifer,1,jj)
      g3r2=gneutul(ifer,1,kk)
      
      ampjetssquarkb=ampjetssquarkb+
     .3d0*ampcharsup(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2,g3l1, g3l2, g3r1, g3r2, 
     . m1, m3, m4, amchar(j),
     .  amchar(k), p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     ./(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p2p4+m4**2-amsupq(kk)**2)

      g2l1=fcharr(ifer, jj, j)
      g2l2=fcharr(ifer, kk, k)
      g2r1=ccharl(ifer, jj, j)
      g2r2=ccharl(ifer, kk, k)
      g3l1=gneutdL(3,1,jj)
      g3l2=gneutdL(3,1, kk)
      g3r1=gneutdr(3,1, jj)
      g3r2=gneutdr(3,1, kk)

      ampjetssquarkb=ampjetssquarkb+3d0*ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)


      g2l1=-DcharL(3,jj,j)
      g2r1=-EcharR(3,jj,j)
      g3r1=gneutuL(ifer,1,jj)
      g3l1=gneutur(ifer,1,jj)
      g2l2=fcharr(ifer, kk, k)
      g2r2=ccharl(ifer, kk, k)
      g3l2=gneutdL(3,1,kk)
      g3r2=gneutdr(3,1,kk)

c----- squark exchange diagrams, same amplitude than slepton
c----- and sneutrino diagrams but couplings need to be changed
      ampjetssquarkb=ampjetssquarkb
     .+6d0*ampcharsupsdownint(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)    
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)


c---- interfernz squark-W-Terme
      enddo
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      g3L=-g2ew/dSqrt(2d0)*(Vckm(ifer,3))
      gXb2r=ccharl(ifer, jj, j)
      gXb2l=fcharr(ifer, jj, j)
      gcd3l=gneutdl(3,1, jj)
      gcd3r=gneutdr(3,1, jj)

      ampjetssquarkb=ampjetssquarkb+
     . 6d0*RealPart(ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2, gxb2l,gXb2R, g2l1, g2r1, g3L, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
 

      gXc2l=dcharl(3, jj, j)
      gXC2r=echarr(3, jj, j)
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
 
      g3L=-g2ew/dSqrt(2d0)*Vckm(ifer,3)
      gbd3R=gneutul(ifer,1, jj)


      ampjetssquarkb=ampjetssquarkb+
     . 6d0*RealPart(ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, g2l1, g2r1, gXc2R, g3L,
     . gbd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      
      enddo
      enddo
      enddo
      enddo
      enddo

c----- contribtutions from diagram invlolving a SM fermion
      ampSMfermion=0d0
c----- SM fermion matrix elemnt squared
       p1pfer=SD_scal(p1,pfer)
      p2pfer=SD_scal(p2,pfer)
      p3pfer=SD_scal(p3,pfer)
      p4pfer=SD_scal(p4,pfer)
      pfer2=SD_scal(pfer,pfer)
      pIpfer=SD_scal(pInt, pfer)



      do j=1,3
      do k=1,3
      do ifer=1,2
      g1L1=-g2ew/dsqrt(2d0)*Vckm(j,ifer)
      g1L2=-g2ew/dsqrt(2d0)*Vckm(k,ifer)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2l2=gneutur(k,1,1)
      g2r2=gneutul(k,1,1)
      g3l=-g2ew/dsqrt(2d0)


      ampSMfermion=ampSMfermion+3d0*(Vckm(1,3)**2+VCKM(2,3)**2)
     .*ampSMfer(g1L1, g1L2,g2l1, g2l2, 
     . g2r1, g2r2,g3l,
     . m1, m3, m4, mfu(j), mfu(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2)/(pfer2-mfu(j)**2)
     . /(pfer2-mfu(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
      
      end do
      enddo
c---- interference term with W-chargino diagram

      do ifer2=1,2
      do k=1,2
      g1R2=DcharL(ifer2,1,k)
      g1L2=EcharR(ifer2,1,k)
      g1L1=-g2ew/dsqrt(2d0)*Vckm(j,ifer2)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2L2=CchaneutL(1,k)
      g2R2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)
      ampSMfermion=ampSMfermion+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*ampintWSMferWchar(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l,  m1, m3, m4, 
     .mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2 

c---- interfernce with sdown contribtuions


      do kk=1,6
      do ifer=1,2
      g2l2=fcharr(ifer, kk, k)
      g2r2=ccharl(ifer, kk, k)
      g3l1=-g2ew/dsqrt(2d0)*Vckm(ifer,3)
      g3r2=gneutdR(3,1, kk)
      g3l2=gneutdL(3,1, kk)
      ampSMfermion=ampSMfermion+6d0*RealPart(
     . ampintWSMferslep(g1L1, g1L2, g1r2,
     .  g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,m1, m3,
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2))

      
      g2l2=-Dcharl(3,kk,k)
      g2r2=-EcharR(3,kk,k)
      g3l2=gneutur(ifer,1, kk)
      g3r2=gneutul(ifer,1, kk)
      ampSMfermion=ampSMfermion+6d0*RealPart(
     . ampintWSMfersneut(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2, m1, m3, 
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p2p4+m4**2-amsupq(kk)**2))
      enddo
      end do
      enddo
      end do
      enddo
c---- diagram with W squark propagator
      ampsquarkpropW=0d0
      do ifer=1,2
      do j=1,6
      do k=1,6
      g1l1=gsqsqW(j)
      g1l2=gsqsqW(k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=gneutdr(ifer,1,k)
      g2r1=gneutdl(ifer,1,j)
      g2r2=gneutdl(ifer,1,k)
      g3l=-g2ew/dsqrt(2d0)
      

     
      ampsquarkpropW=ampsquarkpropW+3d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*ampvecscalar(g1l1,g1l2, 
     . g2l1, g2l2, g2r1, g2r2,g3l,0d0,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, amw)/
     . Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
      enddo
      enddo
      enddo

c----- interference W squark- W chargino diagrams
      do ifer=1,2
      do j=1,6
      do k=1,2
      g1l1=gsqsqW(j)
      g1l2=EcharR(ifer,1,j)
      g1R2=DcharL(ifer,1,k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=CchaneutL(1,k)
      g2r1=gneutdl(ifer,1,j)
      g2r2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)

      ampsquarkpropW=ampsquarkpropW+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*RealPart(
     . ampintWcharWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pi2-amchar(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2)
      enddo
      enddo
      enddo
c------ interfernce W squark - W SM Fermion
      do ifer=1,2
      do j=1,6
      do k=1,3
      
      g1l1=gsqsqW(j)
      g1l2=gneutur(k,1,1)
      g1R2=gneutul(k,1,1)
      g2l1=gneutdr(ifer,1,j)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(k,ifer)
      g2r1=gneutdl(ifer,1,j)
      g3l=-g2ew/dsqrt(2d0)
      ampsquarkpropW=ampsquarkpropW+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     . *RealPart(
     . ampintWferWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pfer2-mfu(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2) 
      enddo
      enddo
      enddo
      
c--------interfernce W squark Sdown diagram
      do ifer2=1,2
      do j=1,6
      do k=1,2
      do kk=1,6
      do ifer=1,2
      g1l1=gsqsqW(j)
      g1l2=EcharR(ifer2,1,k)
      g1R2=DcharL(ifer2,1,k)
      g2l1=gneutdr(ifer2,1,j)
      g2l2=fcharr(ifer, kk, k)
      g2r1=gneutdl(ifer2,1,j)
      g2r2=ccharl(ifer, kk, k)
      g3l1=-g2ew/dsqrt(2d0)*VCKM(ifer,3)
      g3l2=gneutdl(3,1, kk)
      g3r2=gneutdr(3,1, kk)
          
   
      ampsquarkpropW=ampsquarkpropW+6d0*RealPart(
     .ampintcharslepWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     ./(pi2-amchar(k)**2)/(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth) )

      
      
c---- interference W squark- chargino sup diagram      

      g1R2=DcharL(ifer2,1,k)
      g1L2=EcharR(ifer2,1,k)
      g2l2=DcharL(3,kk,k)
      g2r2=EcharR(3,kk,k)
      g3l2=gneutuL(ifer,1, kk)
      g3r2=gneutur(ifer,1, kk)

      ampsquarkpropW=ampsquarkpropW+6d0*
     . RealPart(ampintcharsneutWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     . /(pi2-amchar(k)**2)/(2d0*p2p4+m4**2-amsupq(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo
      enddo
      enddo
        ampchahib=0d0
      do k=1,6
      do j=1,6
      do ifer=1,2
 
      g1l1=gsqsqchaHi(j)
      g1l2=gsqsqchaHi(k)
      g2r1=gneutdl(ifer,1,j)
      g2r2=gneutdl(ifer,1,k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=gneutdr(ifer,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta

      
c----- diagram charged Higgs squark propagator
      ampchahib=ampchahib+(VCKM(1,3)**2+VCKM(2,3)**2)*
     . 3d0*amphiggssquark(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,0d0,g3r1,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4)/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/
     . (2d0*p2p3+m3**2-amch**2)**2


	
       

! 
! ! c---- interference Higgs-squark with W squark diagram
      g1l1=gsqsqW(j)
      g1l2=gsqsqchaHi(k)
      g2r1=gneutdl(ifer,1,j)
      g2r2=gneutdl(ifer,1,k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=gneutdr(ifer,1,k)
      g3l1=-g2ew/dsqrt(2d0)
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta

      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)*RealPart(
     .ampinthiggssquarkWsquark(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,g3l1,0d0,g3r1,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4)/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/
     . (2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
! 
      enddo
      enddo
       enddo
! 
! c----- interference charged Higgs squark- W chargino diagram
      do j=1,6
      do k=1,2
      do ifer=1,2
      g1l1=gsqsqchaHi(j)
      g1l2=EcharR(ifer,1,k)
      g1r2=DcharL(ifer,1,k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=CchaneutL(1,k)
      g2r1=gneutdl(ifer,1,j)
      g2r2=CchaneutR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/dsqrt(2d0)

     
      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)*2)*RealPart(
     . ampintchaHisquarkWchar ( g1l1 , g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo
! 
! c------ interference charged Higgs squark -W SM fermion diagram
      do j=1,6
      do k=1,3
      do ifer=1,2
      g1l1=gsqsqchaHi(j)
      g1l2=gneutur(k,1,1)
      g1r2=gneutul(k,1,1)
      g2l1=gneutdr(ifer,1,j)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(3,k)
      g2r1=gneutdl(ifer,1,j)
      g2r2=0d0
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfe(3)*tanbeta
      g3l2=-g2ew/dsqrt(2d0)

     
      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)*RealPart(
     . ampintchaHisquarkWSMfer( g1l1 , g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(k))/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pfer2-mfu(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo

c---- interference charged Higgs squark- slepton chargino diagram
      do j=1,6
      do k=1,2
      do kk=1,6
      do ifer=1,2
      do ifer2=1,2
      g1l1=gsqsqchaHi(j)
      g1l2=EcharR(ifer2,1,k)
      g1r2=DcharL(ifer2,1,k)
      g2r1=gneutdr(ifer2,1,j)
      g2l1=gneutdl(ifer2,1,j)
      g2r2=ccharl(ifer,kk,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta*VCKM(ifer,3)
      g3l2=gneutdl(3,1,kk)
      g3r2=gneutdr(3,1,kk)


      ampchahib=ampchahib
     . +6d0*RealPart(ampintchaHisquarkslepton(g1l1, g1l2, 
     . g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/(pi2-amsdownq(k)**2)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)
     ./(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      
      g2l2=-DcharL(3,kk,j)
      g2r2=-EcharR(3,kk,j)
      g3r2=gneutul(ifer,1,kk)
      g3l2=gneutur(ifer,1,kk)


            ampchahib=ampchahib
     . -6d0*RealPart(ampintchaHisquarksup(g1l1, g1l2, 
     . g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/(pi2-amchar(k)**2)
     . /(2d0*p2p4+m4**2-amsupq(kk)**2)
     ./(2d0*p2p3+m3**2-amch**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo
      enddo
      enddo

      do ifer=1,2
      
      g2l1=gneutur(3,1,1)
      g2r1=gneutul(3,1,1)
      g2r2=gneutul(3,1,1)
      g2l2=gneutur(3,1,1)
      g1r1=g2ew/dsqrt(2d0)/amw*mfu(3)*1d0/tanbeta*VCKM(3,ifer)
      g1r2=g2ew/dsqrt(2d0)/amw*mfu(3)*1d0/tanbeta*VCKM(3,ifer)
      g1l1=g2ew/dsqrt(2d0)/amw*mfd(ifer)*tanbeta*VCKM(3,ifer)
      g1l2=g2ew/dsqrt(2d0)/amw*mfd(ifer)*tanbeta*VCKM(3,ifer)
      g3r1= g2ew/dsqrt(2d0)/amw*mfd(3)*tanbeta

      ampchahib=ampchahib+3d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     . *ampchaHiSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(3), mfu(3))
     . /(pfer2-mfu(3)**2)/(pfer2-mfu(3)**2)/(2d0*p2p3+m3**2-amch**2)**2
      
      enddo

! ! 
! ! 
! ! c----- charged Higgs diagram: chargino
! 
 
      do j=1,2
      do k=1,2
      do ifer=1,2
      g1l1=EcharR(ifer,1,j)
      g1r1=DcharL(ifer,1,j)
      g1l2=EcharR(ifer,1,k)
      g1r2=DcharL(ifer,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=-g2ew*QchncL(1,k)
      g2r2=-g2ew*QchncR(1,k)
      g3r1= g2ew/dsqrt(2d0)/amw*mfd(3)*tanbeta

      ampchahib=ampchahib+3d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     . *ampchaHichar(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(j), amchar(k))/
     . (pi2-amchar(j)**2)/(pI2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)**2
      enddo
      enddo
       enddo
! ! 
! c----- interference charged Higgs chargino and charged Higgs squark diagrams
      do j=1,2
      do k=1,6
      do ifer=1,2
      g1l1=EcharR(ifer,1,j)
      g1r1=Dcharl(ifer,1,j)
      g1l2=gsqsqchaHi(k)
      g2r1=-g2ew*QchncR(1,j)
      g2l1=-g2ew*QchncL(1,j)
      g2l2=gneutdr(ifer,1,k)
      g2r2=gneutdl(ifer,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta

    
 
 
      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*ampintchaHichargchaHisquark( g1l2 ,
     .  g1l1,  g1r1,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j))/(pI2-amchar(j)**2)
     . /(2d0*p2p3+m3**2-amch**2)**2
     . /(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)  
      enddo
      enddo
      enddo
! ! c------ interference charged Higgs chargino - W squark diagrams
       
       do j=1,2
       do k=1,6
      do ifer=1,2
      g1l1=EcharR(ifer,1,j)
      g1r1=Dcharl(ifer,1,j)
      g1l2=gsqsqW(k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=gneutdr(ifer,1,k)
      g2r2=gneutdl(ifer,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))
     
      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*RealPart(ampintchaHichargWsquark(g1l2,
     .  g1l1,  g1r1,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi,p2pI, p3pI, p4pI, pI2,amchar(j))/(pI2-amchar(j)**2)
     ./(2d0*p2p3+m3**2-amch**2)/(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)/
     . (2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
       enddo
       enddo
       enddo

! ! 
! c----- interfernce charged Higgs chargino-W chargino
       do j=1,2
      do k=1,2
      do ifer=1,2
      g1l1=EcharR(ifer,1,j)
      g1r1=Dcharl(ifer,1,j)
      g1l2=EcharR(ifer,1,k)
      g1r2=Dcharl(ifer,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=CchaneutL(1,k)
      g2r2=CchaneutR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))

     
      ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*RealPart(ampintchaHiWchar(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))
     . /(pi2-amchar(j)**2)/(pi2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo
! ! ! ! 
! ! c----- interfernce charged Higgs chargino-W Sm fermion diagrams
       do j=1,2
      do k=1,3
      do ifer=1,2
      g1l1=EcharR(ifer,1,j)
      g1r1=Dcharl(ifer,1,j)
      g1l2=gneutur(k,1,1)
      g1r2=gneutul(k,1,1)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=-g2ew/dsqrt(2d0)*Vckm(k,ifer)
      g2r2=0d0
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))

     
 

 
      ampchahib=ampchahib-6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*RealPart(ampintchaHiWSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, pIpFer, mfu(k), amchar(j))
     . /(pfer2-mfu(k)**2)/(pi2-amchar(j)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo
! ! 
! ! c----- interference charged Higgs chargino+ slepton diagram
      do kk=1,6
      do j=1,2
      do k=1,2
      do ifer=1,2
      do ifer2=1,2
      g1l1=EcharR(ifer2,1,j)
      g1r1=Dcharl(ifer2,1,j)
      g1l2=EcharR(ifer2,1,k)
      g1r2=Dcharl(ifer2,1,k)
 
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2r2=ccharl(ifer,kk,j)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta*Vckm(ifer,3)
      g3l2=gneutdl(3,1,kk)
      g3r2=gneutdl(3,1,kk)


      ampchahib=ampchahib+6d0*ampintchaHichargslepton(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))
     ./(pi2-amchar(j)**2)
     . /(pi2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)

      enddo
      enddo
      enddo
      enddo
      enddo
! 
! ! 
! c---- interference charged Higgs chargino+sneutrino diagrams
c--- sneutrino diagram can be used, because there is no dependence on g3l2
      do j=1,2
      do k=1,2
      do ifer=1,2
      do kk=1,6
      do ifer2=1,2
      g1l1=EcharR(ifer2,1,j)
      g1r1=Dcharl(ifer2,1,j)
      g1l2=EcharR(ifer2,1,k)
      g1r2=Dcharl(ifer2,1,k)
      g2l1=-g2ew*QchncL(1,j)
      g2r1=-g2ew*QchncR(1,j)
      g2l2=-dcharl(3,kk,k)
      g2r2=-echarr(3,kk,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta*VCKM(ifer,3)
      g3r2=gneutul(ifer,1,kk)


      

      ampchahib=ampchahib+6d0*ampintchaHichargsneutrino(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     . m1,m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amchar(j), amchar(k))/
     .(pi2-amchar(j)**2)/(pi2-amchar(k)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p4+m4**2-amsupq(kk)**2)

      enddo
      enddo
      enddo
      enddo
      enddo

c---- interference charged Higgs Sm fermion with slepton diagram
      do k=1,2
      do kk=1,6
      do ifer=1,2
      do j=1,3
      do ifer2=1,2
      g1l1=gneutur(j,1,1)
      g1r1=gneutul(j,1,1)
      g1l2=EcharR(ifer2,1,k)
      g1r2=Dcharl(ifer2,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j, ifer2)
      g2r1=mfu(j)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j,ifer2)
      g2r2=ccharl(ifer,kk,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta*VCKM(ifer,3)
      g3r2=gneutdr(3,1,kk)
      g3l2=gneutdl(3,1,kk)
      
      ampchahib=ampchahib+6d0*ampintchaHiSmferslepton(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfu(3), amchar(k))
     ./(pi2-amchar(k)**2)/(pfer2-mfu(3)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2) 
      enddo
      enddo
      enddo
      enddo
      enddo
c---- interference charged Higgs Sm fermion with sneutrino diagram
      do k=1,2
      do kk=1,6
      do ifer=1,2
      do j=1,3
      do ifer2=1,2
      g1l1=gneutur(j,1,1)
      g1r1=gneutul(j,1,1)
      g1l2=EcharR(ifer2,1,k)
      g1r2=Dcharl(ifer2,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j,ifer2)
      g2r1=mfu(j)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j, ifer2)
      g2l2=-dcharl(3,kk,k)
      g2r2=-echarr(3,kk,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta*VCKM(ifer,3)
      g3r2=gneutul(ifer,kk,1)


      ampchahib=ampchahib+6d0*ampintchaHiSmfersneutrino(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfu(3), amchar(k))
     ./(pi2-amchar(k)**2)/(pfer2-mfu(3)**2)/(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p4+m4**2-amsupq(kk)**2) 
      enddo
      enddo
      enddo
      enddo
      enddo
! ! c----- interference charged Higgs chargino and charged Higgs SM fermion
      do k=1,2
      do j=1,3
      do ifer=1,2
      g1l1=gneutur(j,1,1)
      g1r1=gneutul(j,1,1)
      g1l2=EcharR(ifer,1,k)
      g1r2=Dcharl(ifer,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j, ifer)
      g2r1=mfu(j)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j, ifer)
      g2l2=-g2ew*QchncL(1,k)
      g2r2=-g2ew*QchncR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta


     
       ampchahib=ampchahib-6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*ampintchaHiSmferchaHicharg(g1l1 ,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfu(3), amchar(k))
     ./(pi2-amchar(k)**2)/(pfer2-mfu(3)**2)
     . /(2d0*p2p3+m3**2-amch**2)**2
      enddo
      enddo
      enddo
! ! ! 
! ! ! c----- interference charged Higgs squark and charged Higgs SM fermion
! !diagrams
      do k=1,6
      do j=1,3
      do ifer=1,2
      g1l1=gneutur(j,1,1)
      g1r1=gneutul(j,1,1)
      g1l2=gsqsqchaHi(k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j,ifer)
      g2r1=mfu(j)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j,ifer)
      g2l2=gneutdr(ifer,1,k)
      g2r2=gneutdl(ifer,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3r2=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta


     
       ampchahib=ampchahib+6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*ampintchaHiSmferchaHisquark(g1l2,
     . g1r1, g1l1,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer,p4pfer, pfer2, mfu(3))
     ./(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)/(pfer2-mfu(3)**2)
     . /(2d0*p2p3+m3**2-amch**2)**2
      enddo
      enddo
      enddo
! ! 
! c----- interference charged Higgs SM fermion and W chargino diagrams
      do k=1,2
      do j=1,2
      do ifer=1,2

      g1l1=gneutur(j,1,1)
      g1r1=gneutul(j,1,1)
      g1r2=EcharR(ifer,1,k)
      g1l2=DcharL(ifer,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j, ifer)
      g2r1=mfu(j)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j, ifer)
      g2l2=CchaneutL(1,k)
      g2r2=CchaneutR(1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))

     
       ampchahib=ampchahib-6d0*(VCKM(1,3)**2+VCKM(2,3)**2)
     .*RealPart(ampintchaHiSmferWcharg(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pipfer, amchar(k), mfu(3))
     . /(pfer2-mfu(3)**2)/(pi2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo

! ! 
! c---- interference charged Higgs SM fermion and W squark diagrams
      do k=1,6
      do ifer=1,2
      do j=1,3
      g1l1=gneutur(j,1,1)
      g1r1=gneutul(j,1,1)
      g1l2=gsqsqW(k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j,ifer)
      g2r1=mfu(j)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j,ifer)
      g2l2=gneutdr(ifer,1,k)
      g2r2=gneutdl(ifer,1,k)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta
      g3l2=-g2ew/(dsqrt(2d0))


	ampchahib=ampchahib-6d0*(VCKM(1,3)**2+VCKM(2,3)**2)*RealPart(
     .ampintchaHiSmferWsquark(g1l1,
     . g1r1, g1l2, g2l1, g2r1, g2l2,g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, mfu(3))  
     . /(pfer2-mfu(3)**2)/(2d0*p1p4+m1**2+m4**2-amsdownq(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo

c---- interference SM fermion W and SM fermion charged Higgs diagrams
      do k=1,3
      do ifer=1,2
      do ifer2=1,2
      do j=1,3
      g1l1=gneutur(j,1,1)
      g1r1=gneutul(j,1,1)
      g1l2=EcharR(ifer,1,k)
      g1r2=Dcharl(ifer,1,k)
      g2l1=mfd(3)*tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j,ifer)
      g2r1=mfu(j)*1d0/tanbeta*g2ew/(dsqrt(2d0)*amw)*VCKM(j,ifer)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(k,3)
      g3r1=g2ew/(dsqrt(2d0)*amw)*mfd(3)*tanbeta*VCKM(ifer2,3)
      g3l2=-g2ew/(dsqrt(2d0))*VCKM(ifer2,3)


     
      ampchahib=ampchahib+6d0
     .*Realpart(ampintchaHiSmferWSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, mfu(3), mfu(k))    
     . /(pfer2-mfu(3)**2)/(pfer2-mfu(k)**2)
     . /(2d0*p2p3+m3**2-amch**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo 
      enddo
      enddo
      enddo
      sigjetb=sigjetb+wt*ampwchajetb/iiter+wt*ampjetssquarkb/iiter
     . +wt*ampSMfermion/iiter+wt*ampsquarkpropW/iiter
     .+wt*ampchahib/iiter




      enddo
      else
      sigjetb=0d0
      print*,"decay stop---> b neutralino + massless fermions forbidden"
      endif

c-------------------------------------------------------------
c----- SIXTH CASE: INITIAL: MASSLESS QUARK, FINAL: MASSLESS
C---------------------------------------------------------
      Etot=amsupq(1)
      Mfin(1)=0d0
      Mfin(2)=0d0
      Mfin(3)=0d0
      MFin(4)=amneut(1)


      do i=1, iiter
      call SD_rambo(4,etot,mfin,pfout,wt)

c---- take the momenta of Rambo 
      do j=1,4
      p1(j)=pfout(j,1)	!quark
      p2(j)=pfout(j,2)	!fermion
      p3(j)=pfout(j,3)	!antifermion
      p4(j)=pfout(j,4)	!neutralino
      pInt(j)=p2(j)+p3(j)+p4(j)	!intermediate particle (chargino)
      pfer(j)=p1(j)+p2(j)+p3(j)
      end do

! define scalar products in amplitude
      p1p2=SD_scal(p1,p2)
      p1p3=SD_scal(p1,p3)
      p1p4=SD_scal(p1,p4)
      p2p3=SD_scal(p2,p3)
      p2p4=SD_scal(p2,p4)
      p3p4=SD_scal(p3,p4)
      p1pI=SD_scal(p1,pInt)
      p2pI=SD_scal(p2,pInt)
      p3pI=SD_scal(p3,pInt)
      p4pI=SD_scal(p4,pInt)
      pI2=SD_scal(pInt,pInt)
      
      m1=Mfin(1)
      m3=0d0
      m4=amneut(1)

      ampwchajetjet=0d0
      ampwchajetmu=0d0
c--- sum up different charginos as intermediate particles
      do j=1,2
      do k=1,2
      g1R1=DcharL(2,1,j)
      g1L1=EcharR(2,1,j)
      g1R2=DcharL(2,1,k)
      g1L2=EcharR(2,1,k)
      g3L=-g2ew/dSqrt(2d0)
      g2L1=CchaneutL(1,j)
      g2L2=CchaneutL(1,k)
      g2R1=CchaneutR(1,j)
      g2R2=CchaneutR(1,k)      
	
      zwi=ampwchar(g1L1, g1L2, g2l1, g2l2, g1r1, 
     . g1r2, g2r1, g2r2, g3l, m1, m3, m4, 
     . amchar(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)

c--- propagators, colour factor and summing up bc and bu contribution
       ampwchajetmu=ampwchajetmu+zwi/(pI2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/Abs(2d0*p2p3-amw**2+(0d0,1d0)*amw*wwidth)**2
       ampwchajetjet=ampwchajetjet+zwi*3d0*(Vckm(1,1)**2+Vckm(2,2)**2+
     . VCKM(1,2)**2+VCKM(2,1)**2)
     . /(pI2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p3-amw**2)**2
c----- the same but up as initail state

      g1R1=DcharL(1,1,j)
      g1L1=EcharR(1,1,j)
      g1R2=DcharL(1,1,k)
      g1L2=EcharR(1,1,k)
      
      zwi=ampwchar(g1L1, g1L2, g2l1, g2l2, g1r1, 
     . g1r2, g2r1, g2r2, g3l, m1, m3, m4, 
     . amchar(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
      
c--- propagators, colour factor and summing up bc and bu contribution
	ampwchajetmu=ampwchajetmu+zwi/(pI2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p3-amw**2)**2
       ampwchajetjet=ampwchajetjet+zwi*3d0*(Vckm(1,1)**2+Vckm(2,2)**2
     . +VCKM(1,2)**2+VCKM(2,1)**2)
     . /(pI2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/Abs(2d0*p2p3-amw**2+(0d0,1d0)*amw*wwidth)**2
      end do
      end do




c---- slepton contributions
      ampinterferenzwjetmu=0d0
      ampinterferenzwjetelec=0d0
      ampslepjetmu=0d0
      ampslepjetelec=0d0
      do j=1,2
      do k=1,2
      do jj=1,6
      do ifer=1,2
      do kk=1,6
      
      g1R1=DcharL(ifer,1,j)
      g1L1=EcharR(ifer,1,j)
      g1R2=DcharL(ifer,1,k)
      g1L2=EcharR(ifer,1,k)
      g2r1=cchaneutslep(2,k,kk)
      g2r2=cchaneutslep(2,j,jj)
      g3l1=gneutelecl(2,1,kk)
      g3l2=gneutelecl(2,1,jj)
      g3r1=gneutelecr(2,1,kk)
      g3r2=gneutelecr(2,1,jj)
      g2l1=0d0
      g2l2=0d0
c------- selectron,... contribtuions
      ampslepjetmu=ampslepjetmu+ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     . m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)

      g2r1=cchaneutslep(1,k,kk)
      g2r2=cchaneutslep(1,j,jj)
      g3l1=gneutelecl(1,1,kk)
      g3l2=gneutelecl(1,1,jj)
      g3r1=gneutelecr(1,1,kk)
      g3r2=gneutelecr(1,1,jj)
       ampslepjetelec=ampslepjetelec+ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)

      enddo
      enddo


c---- interference slepton and sneutrino diagrams
       do ifer=1,2
       gXb2R=cchaneutslep(2,j,jj)
       gXc2L=-dchalepsneut(2,k)
       gXc2R=-echalepsneut(2,k)
       gbd3R=gneutneut(2,1)
       gcd3L=gneutelecl(2,1,jj)
       gcd3R=gneutelecr(2,1,jj)
       gXb2l=0d0
      g1R1=DcharL(ifer,1,j)
      g1L1=EcharR(ifer,1,j)
      g1R2=DcharL(ifer,1,k)
      g1L2=EcharR(ifer,1,k)
      ampslepjetmu=ampslepjetmu+ampintslepsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXb2l,gXb2R, gXc2L, gXc2R,gbd3R, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(2)**2)

       gXb2R=cchaneutslep(1,j,jj)
       gXc2L=-dchalepsneut(1,k)
       gXc2R=-echalepsneut(1,k)
       gbd3R=gneutneut(1,1)
       gcd3L=gneutelecl(1,1,jj)
       gcd3R=gneutelecr(1,1,jj)
       gxb2l=0d0

      ampslepjetelec=ampslepjetelec+ampintslepsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXb2l, gXb2R, gXc2L, gXc2R,gbd3R, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(1)**2)

c----- interference w slepton diagrams
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      gXb2R=cchaneutslep(2,j,jj)
      gcd3L=gneutelecl(2,1,jj)
      gcd3R=gneutelecr(2,1,jj)
      g3L=-g2ew/dSqrt(2d0)
       gXb2l=0d0
      
      ampinterferenzwjetmu=ampinterferenzwjetmu+
     . 2d0*RealPart(ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2, gXb2l,gXb2R, g2l1, g2r1, g3L, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

c----- interference w slepton diagrams
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      gXb2R=cchaneutslep(1,j,jj)
      gcd3L=gneutelecl(1,1,jj)
      gcd3R=gneutelecr(1,1,jj)
      g3L=-g2ew/dSqrt(2d0)
       gXb2l=0d0
      
      ampinterferenzwjetelec=ampinterferenzwjetelec+
     . 2d0*RealPart(ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2,gXb2L, gXb2R, g2l1, g2r1, g3L, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amslepton(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      enddo
      enddo
c----- sneutrino diagram
      do ifer=1,2
      g3r1=gneutneut(2,1)
      g3r2=gneutneut(2,1)
      g3l1=0d0
      g3l2=0d0
      g2l1=-dchalepsneut(2,j)
      g2l2=-dchalepsneut(2,k)
      g2r1=-echalepsneut(2,j)
      g2r2=-echalepsneut(2,k)
      g1R1=DcharL(ifer,1,j)
      g1L1=EcharR(ifer,1,j)
      g1R2=DcharL(ifer,1,k)
      g1L2=EcharR(ifer,1,k)
      zwi=ampcharsneutrino(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2, g3l1, g3l2, g3r1, g3r2,
     . m1, m3, m4, amchar(j),
     .  amchar(k), p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)

      ampslepjetmu=ampslepjetmu+zwi/(Pi2-amchar(j)**2)
     ./(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsneutrino(2)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(2)**2)
      
c----- sneutrino diagram
c---- electron and muon will give the same result for matrixelement since
c---- gneuneut and dchslepsneut and echaslepsneut do not change 

      ampslepjetelec=ampslepjetelec+zwi/(Pi2-amchar(j)**2)
     ./(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsneutrino(1)**2)
     . /(2d0*p2p4+m4**2-amsneutrino(1)**2)


c---- interference with W-contribution
      gXc2L=dchalepsneut(2,j)
      gXC2r=echalepsneut(2,j)
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      g3L=-g2ew/dSqrt(2d0)
      gbd3R=gneutneut(2,1)
     
      ampinterferenzwjetmu=ampinterferenzwjetmu+
     . 2d0*RealPart(ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, g2l1, g2r1, gXc2R, g3L,
     . gbd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(2)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
c---- interference with W-contribution
      gXc2L=dchalepsneut(1,j)
      gXC2r=echalepsneut(1,j)
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      g3L=-g2ew/dSqrt(2d0)
      gbd3R=gneutneut(1,1)
     
      ampinterferenzwjetelec=ampinterferenzwjetelec+
     . 2d0*RealPart(ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, g2l1, g2r1, gXc2R, g3L,
     . gbd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(1)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo

c------ squark contributions
      ampjetsquarkjets=0d0
      do j=1,2
      do k=1,2
      do jj=1,6
      do ifer=1,2
      do ifer2=1,2
      do ifer3=1,2
      do kk=1,6
      
      g1R1=DcharL(ifer3,1,j)
      g1L1=EcharR(ifer3,1,j)
      g1R2=DcharL(ifer3,1,k)
      g1L2=EcharR(ifer3,1,k)
      g2l1=-DcharL(ifer,jj,j)
      g2l2=-DcharL(ifer,kk,k)
      g2r1=-EcharR(ifer,jj,j)
      g2r2=-EcharR(ifer,kk,k)
      g3r1=gneutuL(ifer2,1,jj)
      g3r2=gneutuL(ifer2,1, kk)
      g3l1=gneutur(ifer2,1, jj)
      g3l2=gneutur(ifer2,1, kk)
      
      ampjetsquarkjets=ampjetsquarkjets
     .+3d0*ampcharsup(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2,g3l1, g3l2, g3r1, g3r2, 
     . m1, m3, m4, amchar(j),
     .  amchar(k), p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     ./(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p2p4+m4**2-amsupq(kk)**2)

      g2l1=fcharr(ifer, jj, j)
      g2l2=fcharr(ifer, kk, k)
      g2r1=ccharl(ifer, jj, j)
      g2r2=ccharl(ifer, kk, k)
      g3l1=gneutdL(ifer2,1, jj)
      g3l2=gneutdL(ifer2,1, kk)
      g3r1=gneutdr(ifer2,1, jj)
      g3r2=gneutdr(ifer2,1, kk)

      ampjetsquarkjets=ampjetsquarkjets
     . +3d0*ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)


      g2l1=-DcharL(ifer2,jj,j)
      g2r1=-EcharR(ifer2,jj,j)
      g3r1=gneutuL(ifer,1,jj)
      g3l1=gneutur(ifer,1,jj)
      g2l2=fcharr(ifer, kk, k)
      g2r2=ccharl(ifer, kk, k)
      g3l2=gneutdL(ifer2,1,kk)
      g3r2=gneutdr(ifer2,1,kk)

c----- squark exchange diagrams, same amplitude than slepton
c----- and sneutrino diagrams but couplings need to be changed
      ampjetsquarkjets=ampjetsquarkjets
     .+6d0*ampcharsupsdownint(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)    
     ./(Pi2-amchar(j)**2)/(pI2-amchar(k)**2)
     ./(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)


c---- interfernz squark-W-Terme
      enddo
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
      g3L=-g2ew/dSqrt(2d0)*VCKM(ifer,ifer2)
      gXb2l=ccharl(ifer, jj, j)
      gXb2R=fcharr(ifer, jj, j)
      gcd3l=gneutdl(ifer2,jj,1)
      gcd3r=gneutdr(ifer2,jj,1)

      ampjetsquarkjets=ampjetsquarkjets+
     . 6d0*RealPart(ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2, gxb2l,gXb2R, g2l1, g2r1, g3L, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)
     ./(2d0*p3p4+m3**2+m4**2-amsdownq(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
 

      gXc2R=dcharl(ifer2, jj, j)
      gXC2L=echarr(ifer2, jj, j)
      g2l1=CchaneutL(1,k)
      g2R1=CchaneutR(1,k)
 
      g3L=-g2ew/dSqrt(2d0)*VCKM(ifer,ifer2)
      gbd3R=gneutul(ifer, jj, 1)


      ampjetsquarkjets=ampjetsquarkjets+
     . 6d0*(Vckm(1,ifer2)**2+Vckm(2,ifer2)**2)
     .*RealPart(ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, g2l1, g2r1, gXc2R, g3L,
     . gbd3R, m1, m3,
     .  m4, amchar(j), amchar(k), p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)/(Pi2-amchar(j)**2)
     . /(pI2-amchar(k)**2)/(2d0*p2p4+m4**2-amsupq(jj)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo


c----- contribtutions from diagram invlolving a SM fermion
      ampSMfermionel=0d0
      ampSmfermionmu=0d0
c----- SM fermion matrix elemnt squared
       p1pfer=SD_scal(p1,pfer)
      p2pfer=SD_scal(p2,pfer)
      p3pfer=SD_scal(p3,pfer)
      p4pfer=SD_scal(p4,pfer)
      pfer2=SD_scal(pfer,pfer)
      pIpfer=SD_scal(pInt, pfer)

      zwi=0d0
      zwi2=0d0

      do j=1,3
      do k=1,3
      do ifer=1,2
      g1L1=-g2ew/dsqrt(2d0)*VCKM(j,ifer)
      g1L2=-g2ew/dsqrt(2d0)*VCKM(k,ifer)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2l2=gneutur(k,1,1)
      g2r2=gneutul(k,1,1)
      g3l=-g2ew/dsqrt(2d0)


       zwi=zwi+ampSMfer(g1L1, g1L2, g2l1, g2l2, 
     . g2r1, g2r2,g3l,
     . m1, m3, m4, mfu(j), mfu(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2)/(pfer2-mfu(j)**2)
     . /(pfer2-mfu(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
      
      enddo
      end do
      enddo
c---- interference term with W-chargino diagram
      ampSMfermionel=zwi
      ampSMfermionmu=zwi

      do ifer2=1,2
      do k=1,2
      do j=1,3
      g1L1=-g2ew/dsqrt(2d0)*VCKM(j,ifer2)
      g1R2=DcharL(ifer2,1,k)
      g1L2=EcharR(ifer2,1,k)
      g2L2=CchaneutL(1,k)
      g2R2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)

      ampSMfermionmu=ampSMfermionmu+2d0*ampintWSMferWchar(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l,  m1, m3, m4, 
     .mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2 
  
      enddo
      enddo
      end do
      ampSMfermionel=ampSMfermionmu
c---- interfernce with slpeton contribtuions

      do ifer2=1,2
      do k=1,2
      do j=1,3
      do kk=1,6
      g1L1=-g2ew/dsqrt(2d0)*VCKM(j,ifer2)
      g1r2=DcharL(ifer2,1,k)
      g1l2=EcharR(ifer2,1,k)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2l2=0d0
      g2r2=cchaneutslep(1,k,kk)
      g3l1=-g2ew/dsqrt(2d0)
      g3r2=gneutelecr(1,1,kk)
      g3l2=gneutelecl(1,1,kk)
      ampSMfermionel=ampSMfermionel+2d0*RealPart(
     . ampintWSMferslep(g1L1, g1L2, g1r2,
     .  g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,m1, m3,
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2))

      g2r2=cchaneutslep(2,k,kk)
      g3l1=-g2ew/dsqrt(2d0)
      g3r2=gneutelecr(2,1,kk)
      g3l2=gneutelecl(2,1,kk)
      ampSMfermionmu=ampSMfermionmu+2d0*RealPart(
     . ampintWSMferslep(g1L1, g1L2, g1r2,
     .  g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,m1, m3,
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2))

      enddo
      g2l2=-dchalepsneut(1,k)
      g2r2=-echalepsneut(1,k)
      g3l2=0d0
      g3r2=gneutneut(1,1)
      ampSMfermionel=ampSMfermionel+2d0*RealPart(
     . ampintWSMfersneut(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2, m1, m3, 
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p2p4+m4**2-amsneutrino(1)**2))

      g2l2=-dchalepsneut(2,k)
      g2r2=-echalepsneut(2,k)
      g3l2=0d0
      g3r2=gneutneut(2,1)
      ampSMfermionmu=ampSMfermionmu+2d0*RealPart(
     . ampintWSMfersneut(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2, m1, m3, 
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p2p4+m4**2-amsneutrino(2)**2))

      enddo
      enddo
      end do
      

c----- contribtutions from diagram invlolving a SM fermion
      ampSMfermionjet=0d0
c----- SM fermion matrix elemnt squared
       p1pfer=SD_scal(p1,pfer)
      p2pfer=SD_scal(p2,pfer)
      p3pfer=SD_scal(p3,pfer)
      p4pfer=SD_scal(p4,pfer)
      pfer2=SD_scal(pfer,pfer)
      pIpfer=SD_scal(pInt, pfer)



      do j=1,3
      do k=1,3
      do ifer=1,2
      g1L1=-g2ew/dsqrt(2d0)*(VCKM(j,ifer))
      g1L2=-g2ew/dsqrt(2d0)*VCKM(k,ifer)
      g2L1=gneutur(j,1,1)
      g2r1=gneutul(j,1,1)
      g2l2=gneutur(k,1,1)
      g2r2=gneutul(k,1,1)
      g3l=-g2ew/dsqrt(2d0)


      ampSMfermionjet=ampSMfermionjet+
     .3d0*(Vckm(1,1)**2+VCKM(2,1)**2+VCKM(1,2)**2+VCKM(2,2)**2)*
     .ampSMfer(g1L1, g1L2, g2l1, g2l2, 
     . g2r1, g2r2,g3l,
     . m1, m3, m4, mfu(j), mfu(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2)/(pfer2-mfu(j)**2)
     . /(pfer2-mfu(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2
      enddo
      end do
      
c---- interference term with W-chargino diagram

      do jj=1,2
      do k=1,2
      g1L1=-g2ew/dsqrt(2d0)*(VCKM(j,jj))
      g1L2=-g2ew/dsqrt(2d0)*VCKM(k,jj)
      g1R2=DcharL(jj,1,k)
      g1L2=EcharR(jj,1,k)
      g2L2=CchaneutL(1,k)
      g2R2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)
      ampSMfermionjet=ampSMfermionjet+
     .6d0*(Vckm(1,1)**2+VCKM(2,1)**2+VCKM(1,2)**2+VCKM(2,2)**2)*
     .ampintWSMferWchar(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l,  m1, m3, m4, 
     .mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2 

c---- interfernce with sdown contribtuions


      do kk=1,6
      do ifer=1,2
      do ifer2=1,2
      g2l2=ccharl(ifer, kk, k)
      g2r2=ccharl(ifer, kk, k)
      g3l1=-g2ew/dsqrt(2d0)*VCKM(ifer,ifer2)
      g3r2=gneutdR(ifer2,1, kk)
      g3l2=gneutdL(ifer2,1, kk)
      ampSMfermionjet=ampSMfermionjet+6d0*RealPart(
     . ampintWSMferslep(g1L1, g1L2, g1r2,
     .  g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,m1, m3,
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2))

      
      g2L2=Dcharl(ifer2,kk,k)
      g2r2=EcharR(ifer2,kk,k)
      g3r2=gneutur(ifer,1, kk)
      g3l2=gneutul(ifer,1, kk)
      ampSMfermionjet=ampSMfermionjet+6d0*RealPart(
     . ampintWSMfersneut(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2, m1, m3, 
     . m4, mfu(j), amchar(k), p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)/(pfer2-mfu(j)**2)/(pI2-amchar(k)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)
     . /(2d0*p2p4+m4**2-amsupq(kk)**2))
      enddo
      end do
      enddo
      end do
      enddo
      enddo

c---- diagram with W squark propagator
      ampsquarkpropW=0d0
      do ifer=1,2
      do j=1,6
      do k=1,6
      g1l1=gsqsqW(j)
      g1l2=gsqsqW(k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=gneutdr(ifer,1,k)
      g2r1=gneutdl(ifer,1,j)
      g2r2=gneutdl(ifer,1,k)
      g3l=-g2ew/dsqrt(2d0)
      

      ampsquarkpropW=ampsquarkpropW+ampvecscalar(g1l1, g1l2, 
     . g2l1, g2l2, g2r1, g2r2,g3l,0d0,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, amw)/
     . Abs(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(k)**2)
     . /(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
      enddo
      enddo
      end do
c---- no charged Higgs diagrams


c----- interference W squark- W chargino diagrams
      do ifer=1,2
      do j=1,6
      do k=1,2
      g1l1=gsqsqW(j)
      g1l2=EcharR(ifer,1,j)
      g1R2=DcharL(ifer,1,k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=CchaneutL(1,k)
      g2r1=gneutdl(ifer,1,j)
      g2r2=CchaneutR(1,k)
      g3l=-g2ew/dsqrt(2d0)

      ampsquarkpropW=ampsquarkpropW+2d0*RealPart(
     . ampintWcharWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pi2-amchar(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2)
      enddo
      enddo
      enddo
c------ interfernce W squark - W SM Fermion
      do ifer=1,2
      do j=1,6
      do k=1,3
      g1l1=gsqsqW(j)
      g1l2=gneutur(k,1,1)
      g1R2=gneutul(k,1,1)
      g2l1=gneutdr(ifer,1,j)
      g2l2=-g2ew/dsqrt(2d0)*VCKM(k,ifer)
      g2r1=gneutdl(ifer,1,j)
      g3l=-g2ew/dsqrt(2d0)
      ampsquarkpropW=ampsquarkpropW+2d0*RealPart(
     . ampintWferWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfu(k))/
     .(2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)/(pfer2-mfu(k)**2)
     ./(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth)**2) 
      enddo
      enddo
      enddo

      ampsquarkpropWjet=0d0
c--------interfernce W squark Sdown diagram
      do jj=1,2
      do j=1,6
      do k=1,2
      do kk=1,6
      do ifer=1,2
      do ifer2=1,2
      g1l1=gsqsqW(j)
      g1l2=EcharR(jj,1,k)
      g1R2=DcharL(jj,1,k)
      g2l1=gneutdr(jj,1,j)
      g2l2=fcharr(ifer, kk, k)
      g2r1=gneutdl(jj,1,j)
      g2r2=ccharl(ifer, kk, k)
      g3l1=-g2ew/dsqrt(2d0)*(VCKM(ifer,ifer2))
      g3l2=gneutdl(ifer2,kk,1)
      g3r2=gneutdr(ifer2,kk,1)
          
   
      ampsquarkpropWjet=ampsquarkpropWjet+6d0*RealPart(
     .ampintcharslepWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     ./(pi2-amchar(k)**2)/(2d0*p3p4+m3**2+m4**2-amsdownq(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth) )

      
      
c---- interference W squark- chargino sup diagram      

      g1R2=DcharL(jj,1,k)
      g1L2=EcharR(jj,1,k)
      g2l2=-DcharL(ifer2,kk,k)
      g2r2=-EcharR(ifer2,kk,k)
      g3l2=gneutuL(ifer,kk,1)
      g3r2=gneutur(ifer,kk,1)

      ampsquarkpropWjet=ampsquarkpropWjet+6d0*
     . RealPart(ampintcharsneutWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     . /(pi2-amchar(k)**2)/(2d0*p2p4+m4**2-amsupq(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo

c--------interfernce W squark Slepton diagram

      ampsquarkpropWmu=0d0
      ampsquarkpropWe=0d0
      
      do ifer=1,2
      do j=1,6
      do k=1,2
      do kk=1,6
      g1l1=gsqsqW(j)
      g1l2=EcharR(ifer,1,k)
      g1R2=DcharL(ifer,1,k)
      g2l1=gneutdr(ifer,1,j)
      g2l2=0d0
      g2r1=gneutdl(ifer,1,j)
      g2r2=cchaneutslep(2,k,kk)
      g3l1=-g2ew/dsqrt(2d0)
      g3l2=gneutelecl(2,1,kk)
      g3r2=gneutelecr(2,1,kk)

          
   
      ampsquarkpropWmu=ampsquarkpropWmu+2d0*RealPart(
     .ampintcharslepWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     ./(pi2-amchar(k)**2)/(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth) )

      g2r2=cchaneutslep(1,k,kk)
      g3l2=gneutelecl(1,1,kk)
      g3r2=gneutelecr(1,1,kk)
      ampsquarkpropWe=ampsquarkpropWe+2d0*RealPart(
     .ampintcharslepWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     ./(pi2-amchar(k)**2)/(2d0*p3p4+m3**2+m4**2-amslepton(kk)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth) )
      
      enddo
c---- interference W squark- chargino sneutrino diagram      
      g2l2=dchalepsneut(2,k)
      g2r2=echalepsneut(2,k)
      g3r2=gneutneut(2,1)
      g3l2=0d0



      ampsquarkpropWmu=ampsquarkpropWmu+2d0*
     . RealPart(ampintcharsneutWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     . /(pi2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(2)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))

      g2l2=dchalepsneut(1,k)
      g2r2=echalepsneut(1,k)
      g3r2=gneutneut(1,1)
      g3l2=0d0



      ampsquarkpropWe=ampsquarkpropWe+2d0*
     . RealPart(ampintcharsneutWsquark(g1l1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amchar(k))
     . /(pi2-amchar(k)**2)/(2d0*p2p4+m4**2-amsneutrino(1)**2)/
     . (2d0*p1p4+m4**2+m1**2-amsdownq(j)**2)
     . /(2d0*p2p3+m3**2-amw**2+(0d0,1d0)*amw*wwidth))
      
      enddo
      enddo
      enddo
      sigjetelec=sigjetelec+wt*ampwchajetmu/iiter
     .+wt*ampslepjetelec/iiter
     .+wt*ampinterferenzwjetelec/iiter+wt*ampSMfermionel/iiter
     .+wt*ampsquarkpropW/iiter+wt*ampsquarkpropWe/iiter

      sigjetjet=sigjetjet+wt*ampwchajetjet/iiter+
     . wt*ampjetsquarkjets/iiter+wt*ampSMfermionjet/iiter
     .+wt*ampsquarkpropW*3d0
     .*(VCKM(1,1)**2+VCKM(1,2)**2+VCKM(2,1)**2+VCKM(2,2)**2)/iiter
     . +wt*ampsquarkpropWjet/iiter




      sigjetmu=sigjetmu+wt*ampwchajetmu/iiter
     .+wt*ampslepjetmu/iiter
     .+wt*ampinterferenzwjetmu/iiter+wt*ampSMfermionmu/iiter
     .+wt*ampsquarkpropW/iiter+wt*ampsquarkpropWmu/iiter

      enddo

c---------------------------------------------------------
c---- DECAY WIDTHS
C---------------------------------------------------------


13      widthbtau=1d0/(2d0*pi)**8*sigbtau
     .*1d0/(2d0*amsupq(1))
      widthbmu=1d0/(2d0*pi)**8*sigbmu
     .*1d0/(2d0*amsupq(1)) 
      widthbjets=1d0/(2d0*pi)**8*sigbjets
     .*1d0/(2d0*amsupq(1)) 
      widthbbjet=1d0/(2d0*pi)**8*sigbbjet
     .*1d0/(2d0*amsupq(1))
      widthjettau=1d0/(2d0*pi)**8*sigjettau
     .*1d0/(2d0*amsupq(1))
      widthjetb=1d0/(2d0*pi)**8*sigjetb
     .*1d0/(2d0*amsupq(1))
       widthjetjet=1d0/(2d0*pi)**8*sigjetjet
     .*1d0/(2d0*amsupq(1))
       widthjetmu=1d0/(2d0*pi)**8*sigjetmu
     .*1d0/(2d0*amsupq(1))
       widthbelec=1d0/(2d0*pi)**8*sigbelec
     .*1d0/(2d0*amsupq(1)) 
       widthjetelec=1d0/(2d0*pi)**8*sigjetelec
     .*1d0/(2d0*amsupq(1)) 
       totwidth=widthbtau+widthbbjet+widthbjets
     .+widthbmu+widthbelec+widthjettau+widthjetb+
     .widthjetjet+widthjetmu+widthjetelec

!       print*, "result", totwidth, widthbtau,
!      .widthbbjet,widthbjets,
!      ."mu", widthbmu,"elec",widthbelec,widthjettau,widthjetb,
!      .widthjetjet,"MU",widthjetmu, "ELEC",widthjetelec
!       call CPU_TIME(time)
!       print*, "total runtime", time
      end
c-----------------------------------------------------------------
c------- Interferenz W slepton Beitraege
c------- Attention: No factor 2 inside yet
c----------------------------------------------------------------
      double precision function ampintwslep(g1L1, g1L2, g1r1, 
     . g1r2, gbX2L, gbX2R, gXdW2R, gXdW2L, gwcbL, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amcha1, amcha2, p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
      implicit none
      double precision g1L1, g1L2,g3l,gbX2L, gbX2R, gXdW2L, gXdW2R,
     . gcd3R, g1r1, gwcbl, gcd3L, 
     .  m1, m3, m4, amcha1, amcha2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4
      double precision sdgf,amz,amw,pi,g2
      COMMON/SD_param/sdgf,amz,amw,pi,g2
      ampintwslep=
     -  -8*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2L*m3*m4*amcha1*amcha2*
     -   p1p2 + (2*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2L*m3**3*
     -     m4*amcha1*amcha2*p1p2)/amw**2 + 
     -  4*gbX2R*gcd3L*gWcbL*g1l2*g1r1*gXdW2R*m1*m4*amcha1*amcha2*
     -   p2p3 + (2*gbX2R*gcd3L*gWcbL*g1l2*g1r1*gXdW2R*m1*
     -     m3**2*m4*amcha1*amcha2*p2p3)/amw**2 + 
     -  (4*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2L*m3*m4*amcha1*amcha2*
     -     p1p2*p2p3)/amw**2 + 
     -  (2*gbX2R*gcd3L*gWcbL*g1r1*g1r2*gXdW2L*m3**2*amcha1*amcha2*
     -     p1p4*p2p3)/amw**2 - 
     -  4*gbX2R*gcd3L*gWcbL*g1r1*g1r2*gXdW2R*m4*amcha1*p1pI*
     -   p2p3 - 4*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2R*m4*amcha2*
     -   p1pI*p2p3 - (2*gbX2R*gcd3L*gWcbL*g1r1*g1r2*gXdW2R*
     -     m3**2*m4*amcha1*p1pI*p2p3)/amw**2 - 
     -  (2*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2R*m3**2*m4*amcha2*
     -     p1pI*p2p3)/amw**2 - 
     -  4*gbX2R*gcd3R*gWcbL*g1l2*g1r1*gXdW2R*m1*m3*amcha1*amcha2*
     -   p2p4 + (2*gbX2R*gcd3R*gWcbL*g1l2*g1r1*gXdW2R*m1*
     -     m3**3*amcha1*amcha2*p2p4)/amw**2 - 
     -  (4*gbX2R*gcd3L*gWcbL*g1r1*g1r2*gXdW2L*m3**2*amcha1*amcha2*
     -     p1p2*p2p4)/amw**2 - 
     -  (2*gbX2R*gcd3L*gWcbL*g1r1*g1r2*gXdW2L*m3**2*amcha1*amcha2*
     -     p1p3*p2p4)/amw**2 + 
     -  4*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*m3*amcha1*p1pI*
     -   p2p4 + 4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3*amcha2*
     -   p1pI*p2p4 - (2*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*
     -     m3**3*amcha1*p1pI*p2p4)/amw**2 - 
     -  (2*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3**3*amcha2*
     -     p1pI*p2p4)/amw**2 + 
     -  8*gbX2R*gcd3R*gWcbL*g1l2*g1r1*gXdW2L*m1*m3*m4*amcha1*
     -   p2pI + 8*gbX2R*gcd3R*gWcbL*g1l1*g1r2*gXdW2L*m1*m3*
     -   m4*amcha2*p2pI - (2*gbX2R*gcd3R*gWcbL*g1l2*g1r1*
     -     gXdW2L*m1*m3**3*m4*amcha1*p2pI)/amw**2 - 
     -  (2*gbX2R*gcd3R*gWcbL*g1l1*g1r2*gXdW2L*m1*m3**3*m4*
     -     amcha2*p2pI)/amw**2 + 
     -  4*gbX2R*gcd3L*gWcbL*g1r1*g1r2*gXdW2R*m4*amcha1*p1p3*
     -   p2pI - 4*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2R*m4*amcha2*
     -   p1p3*p2pI + (2*gbX2R*gcd3L*gWcbL*g1r1*g1r2*gXdW2R*
     -     m3**2*m4*amcha1*p1p3*p2pI)/amw**2 - 
     -  (2*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2R*m3**2*m4*amcha2*
     -     p1p3*p2pI)/amw**2 - 
     -  4*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*m3*amcha1*p1p4*
     -   p2pI + 4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3*amcha2*
     -   p1p4*p2pI + (2*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*
     -     m3**3*amcha1*p1p4*p2pI)/amw**2 - 
     -  (2*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3**3*amcha2*
     -     p1p4*p2pI)/amw**2 - 
     -  16*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2L*m3*m4*p1pI*
     -   p2pI + (4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2L*m3**3*
     -     m4*p1pI*p2pI)/amw**2 - 
     -  (4*gbX2R*gcd3R*gWcbL*g1l2*g1r1*gXdW2L*m1*m3*m4*amcha1*
     -     p2p3*p2pI)/amw**2 - 
     -  (4*gbX2R*gcd3R*gWcbL*g1l1*g1r2*gXdW2L*m1*m3*m4*amcha2*
     -     p2p3*p2pI)/amw**2 + 
     -  (4*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*m3*amcha1*p1p4*
     -     p2p3*p2pI)/amw**2 - 
     -  (4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3*amcha2*p1p4*
     -     p2p3*p2pI)/amw**2 + 
     -  (8*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2L*m3*m4*p1pI*
     -     p2p3*p2pI)/amw**2 + 
     -  (4*gbX2R*gcd3L*gWcbL*g1l2*g1r1*gXdW2L*m1*m3**2*amcha1*
     -     p2p4*p2pI)/amw**2 + 
     -  (4*gbX2R*gcd3L*gWcbL*g1l1*g1r2*gXdW2L*m1*m3**2*amcha2*
     -     p2p4*p2pI)/amw**2 - 
     -  (4*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*m3*amcha1*p1p3*
     -     p2p4*p2pI)/amw**2 + 
     -  (4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3*amcha2*p1p3*
     -     p2p4*p2pI)/amw**2 - 
     -  (8*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2L*m3**2*p1pI*
     -     p2p4*p2pI)/amw**2 + 
     -  8*gbX2R*gcd3L*gWcbL*g1r1*g1r2*gXdW2L*amcha1*amcha2*p1p2*
     -   p3p4 - (2*gbX2R*gcd3L*gWcbL*g1r1*g1r2*gXdW2L*m3**2*
     -     amcha1*amcha2*p1p2*p3p4)/amw**2 - 
     -  (4*gbX2R*gcd3R*gWcbL*g1l2*g1r1*gXdW2R*m1*m3*amcha1*amcha2*
     -     p2p3*p3p4)/amw**2 + 
     -  (4*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*m3*amcha1*p1pI*
     -     p2p3*p3p4)/amw**2 + 
     -  (4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3*amcha2*p1pI*
     -     p2p3*p3p4)/amw**2 - 
     -  8*gbX2R*gcd3L*gWcbL*g1l2*g1r1*gXdW2L*m1*amcha1*p2pI*
     -   p3p4 - 8*gbX2R*gcd3L*gWcbL*g1l1*g1r2*gXdW2L*m1*amcha2*
     -   p2pI*p3p4 + (2*gbX2R*gcd3L*gWcbL*g1l2*g1r1*gXdW2L*
     -     m1*m3**2*amcha1*p2pI*p3p4)/amw**2 + 
     -  (2*gbX2R*gcd3L*gWcbL*g1l1*g1r2*gXdW2L*m1*m3**2*amcha2*
     -     p2pI*p3p4)/amw**2 - 
     -  (4*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*m3*amcha1*p1p3*
     -     p2pI*p3p4)/amw**2 + 
     -  (4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3*amcha2*p1p3*
     -     p2pI*p3p4)/amw**2 + 
     -  16*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2L*p1pI*p2pI*
     -   p3p4 - (4*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2L*m3**2*
     -     p1pI*p2pI*p3p4)/amw**2 - 
     -  4*gbX2R*gcd3L*gWcbL*g1r1*g1r2*gXdW2R*m4*amcha1*p1p2*
     -   p3pI + 4*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2R*m4*amcha2*
     -   p1p2*p3pI - (2*gbX2R*gcd3L*gWcbL*g1r1*g1r2*gXdW2R*
     -     m3**2*m4*amcha1*p1p2*p3pI)/amw**2 + 
     -  (2*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2R*m3**2*m4*amcha2*
     -     p1p2*p3pI)/amw**2 + 
     -  (2*gbX2R*gcd3L*gWcbL*g1l2*g1r1*gXdW2L*m1*m3**2*amcha1*
     -     p2p4*p3pI)/amw**2 + 
     -  (2*gbX2R*gcd3L*gWcbL*g1l1*g1r2*gXdW2L*m1*m3**2*amcha2*
     -     p2p4*p3pI)/amw**2 + 
     -  (4*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*m3*amcha1*p1p2*
     -     p2p4*p3pI)/amw**2 - 
     -  (4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3*amcha2*p1p2*
     -     p2p4*p3pI)/amw**2 - 
     -  (4*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2L*m3**2*p1pI*
     -     p2p4*p3pI)/amw**2 + 
     -  (4*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*m3*amcha1*p1p2*
     -     p3p4*p3pI)/amw**2 - 
     -  (4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3*amcha2*p1p2*
     -     p3p4*p3pI)/amw**2 + 
     -  4*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*m3*amcha1*p1p2*
     -   p4pI - 4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3*amcha2*
     -   p1p2*p4pI - (2*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*
     -     m3**3*amcha1*p1p2*p4pI)/amw**2 + 
     -  (2*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3**3*amcha2*
     -     p1p2*p4pI)/amw**2 - 
     -  (2*gbX2R*gcd3L*gWcbL*g1l2*g1r1*gXdW2L*m1*m3**2*amcha1*
     -     p2p3*p4pI)/amw**2 - 
     -  (2*gbX2R*gcd3L*gWcbL*g1l1*g1r2*gXdW2L*m1*m3**2*amcha2*
     -     p2p3*p4pI)/amw**2 - 
     -  (4*gbX2R*gcd3R*gWcbL*g1r1*g1r2*gXdW2R*m3*amcha1*p1p2*
     -     p2p3*p4pI)/amw**2 + 
     -  (4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2R*m3*amcha2*p1p2*
     -     p2p3*p4pI)/amw**2 + 
     -  (4*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2L*m3**2*p1pI*
     -     p2p3*p4pI)/amw**2 + 
     -  8*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2L*m3*m4*p1p2*
     -   pI2 - (2*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2L*m3**3*
     -     m4*p1p2*pI2)/amw**2 + 
     -  4*gbX2R*gcd3L*gWcbL*g1l1*g1r2*gXdW2R*m1*m4*p2p3*
     -   pI2 + (2*gbX2R*gcd3L*gWcbL*g1l1*g1r2*gXdW2R*m1*
     -     m3**2*m4*p2p3*pI2)/amw**2 - 
     -  (4*gbX2R*gcd3R*gWcbL*g1l1*g1l2*gXdW2L*m3*m4*p1p2*
     -     p2p3*pI2)/amw**2 - 
     -  (2*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2L*m3**2*p1p4*
     -     p2p3*pI2)/amw**2 - 
     -  4*gbX2R*gcd3R*gWcbL*g1l1*g1r2*gXdW2R*m1*m3*p2p4*
     -   pI2 + (2*gbX2R*gcd3R*gWcbL*g1l1*g1r2*gXdW2R*m1*
     -     m3**3*p2p4*pI2)/amw**2 + 
     -  (4*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2L*m3**2*p1p2*
     -     p2p4*pI2)/amw**2 + 
     -  (2*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2L*m3**2*p1p3*
     -     p2p4*pI2)/amw**2 - 
     -  8*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2L*p1p2*p3p4*
     -   pI2 + (2*gbX2R*gcd3L*gWcbL*g1l1*g1l2*gXdW2L*m3**2*
     -     p1p2*p3p4*pI2)/amw**2 - 
     -  (4*gbX2R*gcd3R*gWcbL*g1l1*g1r2*gXdW2R*m1*m3*p2p3*
     -     p3p4*pI2)/amw**2
      end
c------------------------------------------------------------------
c------- Interferenz W-sneutrino Beitraege
c----- ATTENTION: NO FACTOR 2 INSIDE YET!!!
c------------------------------------------------------------------
	double precision function ampintwsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXc2L, gXdW2R, gXdW2L, gXc2R, gwcbL, 
     . gbd3R, m1, m3,
     .  m4, amcha1, amcha2, p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
      implicit none
      double precision g1L1, g1L2,g3l,gXc2L, gXc2R, gXdW2L, gXdW2R, 
     . gcd3R, g1r1, gbd3R, gbd3L, gcd3L, gwcbl,
     .  m1, m3, m4, amcha1, amcha2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4
      double precision sdgf,amz,amw,pi,g2
      COMMON/SD_param/sdgf,amz,amw,pi,g2
c---- g1-coupling of first vertex (stop-quark-chargino, 
c--- distinguing between different charginos necessary), 
c-----gXb2-coupling of second vertex (chargino-slepton-neutrino)
c---- gXc2-coupling of second vertex (chargino-sneutrino-lepton)
c---- gcd3 coupling of thrid vertex (slepton-lepton-neutralino)
c----- gbd3 coupling of third vertex (sneutrino-neutrino-neutralino)
c--- m1, p1-> quark
c--- m2, p2-> fermion
c--- m3, p3-> antifermion
c--- m4, p4 -> neutralino
      ampintwsneut=
     -       4*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2L*m3*m4*amcha1*amcha2*
     -   p1p2 - (2*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2L*m3**3*
     -     m4*amcha1*amcha2*p1p2)/amw**2 - 
     -  4*gbd3R*gWcbL*g1L1*g1R2*gXc2L*gXdW2L*m1*m4*amcha1*amcha2*
     -   p2p3 - (2*gbd3R*gWcbL*g1L1*g1R2*gXc2L*gXdW2L*m1*
     -     m3**2*m4*amcha1*amcha2*p2p3)/amw**2 + 
     -  (4*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2L*m3*m4*amcha1*amcha2*
     -     p1p3*p2p3)/amw**2 - 
     -  (2*gbd3R*gWcbL*g1L1*g1L2*gXc2L*gXdW2R*m3**2*amcha1*amcha2*
     -     p1p4*p2p3)/amw**2 + 
     -  4*gbd3R*gWcbL*g1L1*g1L2*gXc2L*gXdW2L*m4*amcha1*p1pI*
     -   p2p3 + 4*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2L*m4*amcha2*
     -   p1pI*p2p3 + (2*gbd3R*gWcbL*g1L1*g1L2*gXc2L*gXdW2L*
     -     m3**2*m4*amcha1*p1pI*p2p3)/amw**2 + 
     -  (2*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2L*m3**2*m4*amcha2*
     -     p1pI*p2p3)/amw**2 + 
     -  8*gbd3R*gWcbL*g1L2*g1R1*gXc2R*gXdW2R*m1*m3*amcha1*amcha2*
     -   p2p4 - (2*gbd3R*gWcbL*g1L2*g1R1*gXc2R*gXdW2R*m1*
     -     m3**3*amcha1*amcha2*p2p4)/amw**2 + 
     -  (4*gbd3R*gWcbL*g1L1*g1L2*gXc2L*gXdW2R*m3**2*amcha1*amcha2*
     -     p1p2*p2p4)/amw**2 - 
     -  8*gbd3R*gWcbL*g1L1*g1L2*gXc2L*gXdW2R*amcha1*amcha2*p1p3*
     -   p2p4 + (2*gbd3R*gWcbL*g1L1*g1L2*gXc2L*gXdW2R*m3**2*
     -     amcha1*amcha2*p1p3*p2p4)/amw**2 - 
     -  8*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2R*m3*amcha1*p1pI*
     -   p2p4 - 8*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2R*m3*amcha2*
     -   p1pI*p2p4 + (2*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2R*
     -     m3**3*amcha1*p1pI*p2p4)/amw**2 + 
     -  (2*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2R*m3**3*amcha2*
     -     p1pI*p2p4)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1L2*g1R1*gXc2R*gXdW2R*m1*m3*amcha1*amcha2*
     -     p2p3*p2p4)/amw**2 + 
     -  (4*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2R*m3*amcha1*p1pI*
     -     p2p3*p2p4)/amw**2 + 
     -  (4*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2R*m3*amcha2*p1pI*
     -     p2p3*p2p4)/amw**2 - 
     -  4*gbd3R*gWcbL*g1L2*g1R1*gXc2R*gXdW2L*m1*m3*m4*amcha1*
     -   p2pI - 4*gbd3R*gWcbL*g1L1*g1R2*gXc2R*gXdW2L*m1*m3*
     -   m4*amcha2*p2pI + (2*gbd3R*gWcbL*g1L2*g1R1*gXc2R*
     -     gXdW2L*m1*m3**3*m4*amcha1*p2pI)/amw**2 + 
     -  (2*gbd3R*gWcbL*g1L1*g1R2*gXc2R*gXdW2L*m1*m3**3*m4*
     -     amcha2*p2pI)/amw**2 + 
     -  4*gbd3R*gWcbL*g1L1*g1L2*gXc2L*gXdW2L*m4*amcha1*p1p3*
     -   p2pI - 4*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2L*m4*amcha2*
     -   p1p3*p2pI - (2*gbd3R*gWcbL*g1L1*g1L2*gXc2L*gXdW2L*
     -     m3**2*m4*amcha1*p1p3*p2pI)/amw**2 + 
     -  (2*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2L*m3**2*m4*amcha2*
     -     p1p3*p2pI)/amw**2 - 
     -  (2*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2R*m3**3*amcha1*
     -     p1p4*p2pI)/amw**2 + 
     -  (2*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2R*m3**3*amcha2*
     -     p1p4*p2pI)/amw**2 + 
     -  8*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2L*m3*m4*p1pI*
     -   p2pI - (4*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2L*m3**3*
     -     m4*p1pI*p2pI)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1L1*g1R2*gXc2L*gXdW2R*m1*m3**2*amcha1*
     -     p2p4*p2pI)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1L2*g1R1*gXc2L*gXdW2R*m1*m3**2*amcha2*
     -     p2p4*p2pI)/amw**2 + 
     -  (4*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2R*m3*amcha1*p1p3*
     -     p2p4*p2pI)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2R*m3*amcha2*p1p3*
     -     p2p4*p2pI)/amw**2 + 
     -  (8*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2R*m3**2*p1pI*
     -     p2p4*p2pI)/amw**2 + 
     -  (2*gbd3R*gWcbL*g1L1*g1L2*gXc2L*gXdW2R*m3**2*amcha1*amcha2*
     -     p1p2*p3p4)/amw**2 - 
     -  (2*gbd3R*gWcbL*g1L1*g1R2*gXc2L*gXdW2R*m1*m3**2*amcha1*
     -     p2pI*p3p4)/amw**2 - 
     -  (2*gbd3R*gWcbL*g1L2*g1R1*gXc2L*gXdW2R*m1*m3**2*amcha2*
     -     p2pI*p3p4)/amw**2 + 
     -  (4*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2R*m3*amcha1*p1p3*
     -     p2pI*p3p4)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2R*m3*amcha2*p1p3*
     -     p2pI*p3p4)/amw**2 + 
     -  (4*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2R*m3**2*p1pI*
     -     p2pI*p3p4)/amw**2 - 
     -  4*gbd3R*gWcbL*g1L1*g1L2*gXc2L*gXdW2L*m4*amcha1*p1p2*
     -   p3pI + 4*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2L*m4*amcha2*
     -   p1p2*p3pI + (2*gbd3R*gWcbL*g1L1*g1L2*gXc2L*gXdW2L*
     -     m3**2*m4*amcha1*p1p2*p3pI)/amw**2 - 
     -  (2*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2L*m3**2*m4*amcha2*
     -     p1p2*p3pI)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1L2*g1R1*gXc2R*gXdW2L*m1*m3*m4*amcha1*
     -     p2p3*p3pI)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1L1*g1R2*gXc2R*gXdW2L*m1*m3*m4*amcha2*
     -     p2p3*p3pI)/amw**2 + 
     -  (4*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2R*m3*amcha1*p1p4*
     -     p2p3*p3pI)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2R*m3*amcha2*p1p4*
     -     p2p3*p3pI)/amw**2 + 
     -  (8*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2L*m3*m4*p1pI*
     -     p2p3*p3pI)/amw**2 + 
     -  8*gbd3R*gWcbL*g1L1*g1R2*gXc2L*gXdW2R*m1*amcha1*p2p4*
     -   p3pI + 8*gbd3R*gWcbL*g1L2*g1R1*gXc2L*gXdW2R*m1*amcha2*
     -   p2p4*p3pI - (2*gbd3R*gWcbL*g1L1*g1R2*gXc2L*gXdW2R*
     -     m1*m3**2*amcha1*p2p4*p3pI)/amw**2 - 
     -  (2*gbd3R*gWcbL*g1L2*g1R1*gXc2L*gXdW2R*m1*m3**2*amcha2*
     -     p2p4*p3pI)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2R*m3*amcha1*p1p2*
     -     p2p4*p3pI)/amw**2 + 
     -  (4*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2R*m3*amcha2*p1p2*
     -     p2p4*p3pI)/amw**2 - 
     -  16*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2R*p1pI*p2p4*
     -   p3pI + (4*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2R*m3**2*
     -     p1pI*p2p4*p3pI)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2R*m3*amcha1*p1p2*
     -     p3p4*p3pI)/amw**2 + 
     -  (4*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2R*m3*amcha2*p1p2*
     -     p3p4*p3pI)/amw**2 + 
     -  (2*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2R*m3**3*amcha1*
     -     p1p2*p4pI)/amw**2 - 
     -  (2*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2R*m3**3*amcha2*
     -     p1p2*p4pI)/amw**2 + 
     -  (2*gbd3R*gWcbL*g1L1*g1R2*gXc2L*gXdW2R*m1*m3**2*amcha1*
     -     p2p3*p4pI)/amw**2 + 
     -  (2*gbd3R*gWcbL*g1L2*g1R1*gXc2L*gXdW2R*m1*m3**2*amcha2*
     -     p2p3*p4pI)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1R1*g1R2*gXc2R*gXdW2R*m3*amcha1*p1p3*
     -     p2p3*p4pI)/amw**2 + 
     -  (4*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2R*m3*amcha2*p1p3*
     -     p2p3*p4pI)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2R*m3**2*p1pI*
     -     p2p3*p4pI)/amw**2 - 
     -  4*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2L*m3*m4*p1p2*
     -   pI2 + (2*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2L*m3**3*
     -     m4*p1p2*pI2)/amw**2 - 
     -  4*gbd3R*gWcbL*g1L2*g1R1*gXc2L*gXdW2L*m1*m4*p2p3*
     -   pI2 - (2*gbd3R*gWcbL*g1L2*g1R1*gXc2L*gXdW2L*m1*
     -     m3**2*m4*p2p3*pI2)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1L1*g1L2*gXc2R*gXdW2L*m3*m4*p1p3*
     -     p2p3*pI2)/amw**2 + 
     -  (2*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2R*m3**2*p1p4*
     -     p2p3*pI2)/amw**2 + 
     -  8*gbd3R*gWcbL*g1L1*g1R2*gXc2R*gXdW2R*m1*m3*p2p4*
     -   pI2 - (2*gbd3R*gWcbL*g1L1*g1R2*gXc2R*gXdW2R*m1*
     -     m3**3*p2p4*pI2)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2R*m3**2*p1p2*
     -     p2p4*pI2)/amw**2 + 
     -  8*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2R*p1p3*p2p4*
     -   pI2 - (2*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2R*m3**2*
     -     p1p3*p2p4*pI2)/amw**2 - 
     -  (4*gbd3R*gWcbL*g1L1*g1R2*gXc2R*gXdW2R*m1*m3*p2p3*
     -     p2p4*pI2)/amw**2 - 
     -  (2*gbd3R*gWcbL*g1R1*g1R2*gXc2L*gXdW2R*m3**2*p1p2*
     -     p3p4*pI2)/amw**2
       end
c-------------------------------------------------------------------
c-------- Interferenz sneutrino-slepton beitraege
c------- ATTENTION: This is 2*Re()!!!
c-------------------------------------------------------------------
	double precision function ampintslepsneut(g1L1, g1L2, g1r1, 
     . g1r2, gXb2L, gXb2R, gXc2L, gXc2R,gbd3R, gcd3L, 
     . gcd3R, m1, m3,
     .  m4, amcha1, amcha2, p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
      implicit none
      double precision g1L1, g1L2,g3l,gXb2L, gXb2R, gXc2L, gXc2R,
     . gcd3R, g1r1, gbd3R, gbd3L, gcd3L, 
     .  m1, m3, m4, amcha1, amcha2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4
      double precision sdgf,amz,amw,pi,g2
c---- g1-coupling of first vertex (stop-quark-chargino, 
c--- distinguing between different charginos necessary), 
c-----gXb2-coupling of second vertex (chargino-slepton-neutrino)
c---- gXc2-coupling of second vertex (chargino-sneutrino-lepton)
c---- gcd3 coupling of thrid vertex (slepton-lepton-neutralino)
c----- gbd3 coupling of third vertex (sneutrino-neutrino-neutralino)
c--- m1, p1-> quark
c--- m2, p2-> fermion
c--- m3, p3-> antifermion
c--- m4, p4 -> neutralino

      ampintslepsneut=
     -  -(4*gbd3L*gcd3R*g1L1*g1L2*gXb2L*gXc2L*m3*m4*amcha1*amcha2*
     -   p1p2 + 4*gbd3R*gcd3L*g1R1*g1R2*gXb2R*gXc2R*m3*m4*
     -   amcha1*amcha2*p1p2 - 4*gbd3R*gcd3L*g1L2*g1R1*gXb2R*gXc2L*
     -   m1*m4*amcha1*amcha2*p2p3 - 
     -  4*gbd3L*gcd3R*g1L1*g1R2*gXb2L*gXc2R*m1*m4*amcha1*amcha2*
     -   p2p3 + 4*gbd3L*gcd3L*g1L1*g1L2*gXb2L*gXc2L*amcha1*amcha2*
     -   p1p4*p2p3 + 4*gbd3R*gcd3R*g1R1*g1R2*gXb2R*gXc2R*
     -   amcha1*amcha2*p1p4*p2p3 + 
     -  4*gbd3R*gcd3L*g1R1*g1R2*gXb2R*gXc2L*m4*amcha1*p1pI*
     -   p2p3 + 4*gbd3L*gcd3R*g1L1*g1L2*gXb2L*gXc2R*m4*amcha1*
     -   p1pI*p2p3 + 4*gbd3R*gcd3L*g1L1*g1L2*gXb2R*gXc2L*m4*
     -   amcha2*p1pI*p2p3 + 4*gbd3L*gcd3R*g1R1*g1R2*gXb2L*
     -   gXc2R*m4*amcha2*p1pI*p2p3 + 
     -  4*gbd3R*gcd3R*g1L2*g1R1*gXb2R*gXc2L*m1*m3*amcha1*amcha2*
     -   p2p4 + 4*gbd3L*gcd3L*g1L1*g1R2*gXb2L*gXc2R*m1*m3*
     -   amcha1*amcha2*p2p4 - 4*gbd3L*gcd3L*g1L1*g1L2*gXb2L*gXc2L*
     -   amcha1*amcha2*p1p3*p2p4 - 
     -  4*gbd3R*gcd3R*g1R1*g1R2*gXb2R*gXc2R*amcha1*amcha2*p1p3*
     -   p2p4 - 4*gbd3R*gcd3R*g1R1*g1R2*gXb2R*gXc2L*m3*amcha1*
     -   p1pI*p2p4 - 4*gbd3L*gcd3L*g1L1*g1L2*gXb2L*gXc2R*m3*
     -   amcha1*p1pI*p2p4 - 4*gbd3R*gcd3R*g1L1*g1L2*gXb2R*
     -   gXc2L*m3*amcha2*p1pI*p2p4 - 
     -  4*gbd3L*gcd3L*g1R1*g1R2*gXb2L*gXc2R*m3*amcha2*p1pI*
     -   p2p4 - 4*gbd3L*gcd3R*g1L1*g1R2*gXb2L*gXc2L*m1*m3*
     -   m4*amcha1*p2pI - 4*gbd3R*gcd3L*g1L2*g1R1*gXb2R*gXc2R*
     -   m1*m3*m4*amcha1*p2pI - 
     -  4*gbd3L*gcd3R*g1L2*g1R1*gXb2L*gXc2L*m1*m3*m4*amcha2*
     -   p2pI - 4*gbd3R*gcd3L*g1L1*g1R2*gXb2R*gXc2R*m1*m3*
     -   m4*amcha2*p2pI - 4*gbd3R*gcd3L*g1R1*g1R2*gXb2R*gXc2L*
     -   m4*amcha1*p1p3*p2pI - 
     -  4*gbd3L*gcd3R*g1L1*g1L2*gXb2L*gXc2R*m4*amcha1*p1p3*
     -   p2pI + 4*gbd3R*gcd3L*g1L1*g1L2*gXb2R*gXc2L*m4*amcha2*
     -   p1p3*p2pI + 4*gbd3L*gcd3R*g1R1*g1R2*gXb2L*gXc2R*m4*
     -   amcha2*p1p3*p2pI + 4*gbd3R*gcd3R*g1R1*g1R2*gXb2R*
     -   gXc2L*m3*amcha1*p1p4*p2pI + 
     -  4*gbd3L*gcd3L*g1L1*g1L2*gXb2L*gXc2R*m3*amcha1*p1p4*
     -   p2pI - 4*gbd3R*gcd3R*g1L1*g1L2*gXb2R*gXc2L*m3*amcha2*
     -   p1p4*p2pI - 4*gbd3L*gcd3L*g1R1*g1R2*gXb2L*gXc2R*m3*
     -   amcha2*p1p4*p2pI + 8*gbd3L*gcd3R*g1R1*g1R2*gXb2L*
     -   gXc2L*m3*m4*p1pI*p2pI + 
     -  8*gbd3R*gcd3L*g1L1*g1L2*gXb2R*gXc2R*m3*m4*p1pI*
     -   p2pI - 4*gbd3L*gcd3L*g1L1*g1L2*gXb2L*gXc2L*amcha1*amcha2*
     -   p1p2*p3p4 - 4*gbd3R*gcd3R*g1R1*g1R2*gXb2R*gXc2R*
     -   amcha1*amcha2*p1p2*p3p4 + 
     -  4*gbd3L*gcd3L*g1L1*g1R2*gXb2L*gXc2L*m1*amcha1*p2pI*
     -   p3p4 + 4*gbd3R*gcd3R*g1L2*g1R1*gXb2R*gXc2R*m1*amcha1*
     -   p2pI*p3p4 + 4*gbd3L*gcd3L*g1L2*g1R1*gXb2L*gXc2L*m1*
     -   amcha2*p2pI*p3p4 + 4*gbd3R*gcd3R*g1L1*g1R2*gXb2R*
     -   gXc2R*m1*amcha2*p2pI*p3p4 - 
     -  8*gbd3L*gcd3L*g1R1*g1R2*gXb2L*gXc2L*p1pI*p2pI*
     -   p3p4 - 8*gbd3R*gcd3R*g1L1*g1L2*gXb2R*gXc2R*p1pI*
     -   p2pI*p3p4 + 4*gbd3R*gcd3L*g1R1*g1R2*gXb2R*gXc2L*m4*
     -   amcha1*p1p2*p3pI + 4*gbd3L*gcd3R*g1L1*g1L2*gXb2L*
     -   gXc2R*m4*amcha1*p1p2*p3pI - 
     -  4*gbd3R*gcd3L*g1L1*g1L2*gXb2R*gXc2L*m4*amcha2*p1p2*
     -   p3pI - 4*gbd3L*gcd3R*g1R1*g1R2*gXb2L*gXc2R*m4*amcha2*
     -   p1p2*p3pI + 4*gbd3L*gcd3L*g1L1*g1R2*gXb2L*gXc2L*m1*
     -   amcha1*p2p4*p3pI + 4*gbd3R*gcd3R*g1L2*g1R1*gXb2R*
     -   gXc2R*m1*amcha1*p2p4*p3pI + 
     -  4*gbd3L*gcd3L*g1L2*g1R1*gXb2L*gXc2L*m1*amcha2*p2p4*
     -   p3pI + 4*gbd3R*gcd3R*g1L1*g1R2*gXb2R*gXc2R*m1*amcha2*
     -   p2p4*p3pI - 8*gbd3L*gcd3L*g1R1*g1R2*gXb2L*gXc2L*
     -   p1pI*p2p4*p3pI - 
     -  8*gbd3R*gcd3R*g1L1*g1L2*gXb2R*gXc2R*p1pI*p2p4*
     -   p3pI - 4*gbd3R*gcd3R*g1R1*g1R2*gXb2R*gXc2L*m3*amcha1*
     -   p1p2*p4pI - 4*gbd3L*gcd3L*g1L1*g1L2*gXb2L*gXc2R*m3*
     -   amcha1*p1p2*p4pI + 4*gbd3R*gcd3R*g1L1*g1L2*gXb2R*
     -   gXc2L*m3*amcha2*p1p2*p4pI + 
     -  4*gbd3L*gcd3L*g1R1*g1R2*gXb2L*gXc2R*m3*amcha2*p1p2*
     -   p4pI - 4*gbd3L*gcd3L*g1L1*g1R2*gXb2L*gXc2L*m1*amcha1*
     -   p2p3*p4pI - 4*gbd3R*gcd3R*g1L2*g1R1*gXb2R*gXc2R*m1*
     -   amcha1*p2p3*p4pI - 4*gbd3L*gcd3L*g1L2*g1R1*gXb2L*
     -   gXc2L*m1*amcha2*p2p3*p4pI - 
     -  4*gbd3R*gcd3R*g1L1*g1R2*gXb2R*gXc2R*m1*amcha2*p2p3*
     -   p4pI + 8*gbd3L*gcd3L*g1R1*g1R2*gXb2L*gXc2L*p1pI*
     -   p2p3*p4pI + 8*gbd3R*gcd3R*g1L1*g1L2*gXb2R*gXc2R*
     -   p1pI*p2p3*p4pI - 
     -  4*gbd3L*gcd3R*g1R1*g1R2*gXb2L*gXc2L*m3*m4*p1p2*
     -   pI2 - 4*gbd3R*gcd3L*g1L1*g1L2*gXb2R*gXc2R*m3*m4*
     -   p1p2*pI2 - 4*gbd3R*gcd3L*g1L1*g1R2*gXb2R*gXc2L*m1*
     -   m4*p2p3*pI2 - 4*gbd3L*gcd3R*g1L2*g1R1*gXb2L*gXc2R*
     -   m1*m4*p2p3*pI2 - 
     -  4*gbd3L*gcd3L*g1R1*g1R2*gXb2L*gXc2L*p1p4*p2p3*pI2 - 
     -  4*gbd3R*gcd3R*g1L1*g1L2*gXb2R*gXc2R*p1p4*p2p3*pI2 + 
     -  4*gbd3R*gcd3R*g1L1*g1R2*gXb2R*gXc2L*m1*m3*p2p4*
     -   pI2 + 4*gbd3L*gcd3L*g1L2*g1R1*gXb2L*gXc2R*m1*m3*
     -   p2p4*pI2 + 4*gbd3L*gcd3L*g1R1*g1R2*gXb2L*gXc2L*
     -   p1p3*p2p4*pI2 + 4*gbd3R*gcd3R*g1L1*g1L2*gXb2R*
     -   gXc2R*p1p3*p2p4*pI2 + 
     -  4*gbd3L*gcd3L*g1R1*g1R2*gXb2L*gXc2L*p1p2*p3p4*pI2 + 
     -  4*gbd3R*gcd3R*g1L1*g1L2*gXb2R*gXc2R*p1p2*p3p4*pI2)
        return
        end

c---------------------------------------------------------------------
c--------- slepton amplitude
c--------------------------------------------------------------------
      double precision function ampcharslepton(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amcha1, amcha2, p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
      implicit none
      double precision g1L1, g1L2, g2l1, g2l2, g3l1,g3l2, g3r1, g3r2,
     . g2r2, g2r1,g1r1, 
     .  m1, m3, m4, amcha1, amcha2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4
      double precision sdgf,amz,amw,pi,g2 
c---- g1-coupling of first vertex (stop-quark-chargino, 
c--- distinguing between different charginos necessary), 
c-----g2-coupling of second vertex (chargino-slepton-neutrino)
c---- g3 coupling of thrid vertex (slepton-lepton-neutralino)
c--- m1, p1-> quark
c--- m2, p2-> fermion
c--- m3, p3-> antifermion
c--- m4, p4 -> neutralino
      ampcharslepton=
     -   -4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*m4*amcha1*amcha2*p1p2 - 
     -  4*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*amcha1*amcha2*p1p2 - 
     -  4*g1L1*g1L2*g2L1*g2L2*g3L1*g3R2*m3*m4*amcha1*amcha2*p1p2 - 
     -  4*g1R1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*m4*amcha1*amcha2*p1p2 + 
     -  4*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*amcha1*p2pI + 
     -  4*g1L2*g1R1*g2R1*g2R2*g3L2*g3R1*m1*m3*m4*amcha1*p2pI + 
     -  4*g1L1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m3*m4*amcha1*p2pI + 
     -  4*g1L2*g1R1*g2R1*g2R2*g3L1*g3R2*m1*m3*m4*amcha1*p2pI + 
     -  4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*amcha2*p2pI + 
     -  4*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m1*m3*m4*amcha2*p2pI + 
     -  4*g1L2*g1R1*g2L1*g2L2*g3L1*g3R2*m1*m3*m4*amcha2*p2pI + 
     -  4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m1*m3*m4*amcha2*p2pI - 
     -  8*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pI*p2pI - 
     -  8*g1L1*g1L2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1pI*p2pI - 
     -  8*g1R1*g1R2*g2L1*g2L2*g3L1*g3R2*m3*m4*p1pI*p2pI - 
     -  8*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m3*m4*p1pI*p2pI + 
     -  4*g1L1*g1L2*g2L1*g2L2*g3L1*g3L2*amcha1*amcha2*p1p2*p3p4 + 
     -  4*g1R1*g1R2*g2R1*g2R2*g3L1*g3L2*amcha1*amcha2*p1p2*p3p4 + 
     -  4*g1L1*g1L2*g2L1*g2L2*g3R1*g3R2*amcha1*amcha2*p1p2*p3p4 + 
     -  4*g1R1*g1R2*g2R1*g2R2*g3R1*g3R2*amcha1*amcha2*p1p2*p3p4 - 
     -  4*g1L1*g1R2*g2L1*g2L2*g3L1*g3L2*m1*amcha1*p2pI*p3p4 - 
     -  4*g1L2*g1R1*g2R1*g2R2*g3L1*g3L2*m1*amcha1*p2pI*p3p4 - 
     -  4*g1L1*g1R2*g2L1*g2L2*g3R1*g3R2*m1*amcha1*p2pI*p3p4 - 
     -  4*g1L2*g1R1*g2R1*g2R2*g3R1*g3R2*m1*amcha1*p2pI*p3p4 - 
     -  4*g1L2*g1R1*g2L1*g2L2*g3L1*g3L2*m1*amcha2*p2pI*p3p4 - 
     -  4*g1L1*g1R2*g2R1*g2R2*g3L1*g3L2*m1*amcha2*p2pI*p3p4 - 
     -  4*g1L2*g1R1*g2L1*g2L2*g3R1*g3R2*m1*amcha2*p2pI*p3p4 - 
     -  4*g1L1*g1R2*g2R1*g2R2*g3R1*g3R2*m1*amcha2*p2pI*p3p4 + 
     -  8*g1R1*g1R2*g2L1*g2L2*g3L1*g3L2*p1pI*p2pI*p3p4 + 
     -  8*g1L1*g1L2*g2R1*g2R2*g3L1*g3L2*p1pI*p2pI*p3p4 + 
     -  8*g1R1*g1R2*g2L1*g2L2*g3R1*g3R2*p1pI*p2pI*p3p4 + 
     -  8*g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p1pI*p2pI*p3p4 + 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1p2*pI2 + 
     -  4*g1L1*g1L2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1p2*pI2 + 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L1*g3R2*m3*m4*p1p2*pI2 + 
     -  4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m3*m4*p1p2*pI2 - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L1*g3L2*p1p2*p3p4*pI2 - 
     -  4*g1L1*g1L2*g2R1*g2R2*g3L1*g3L2*p1p2*p3p4*pI2 - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3R1*g3R2*p1p2*p3p4*pI2 - 
     -  4*g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p1p2*p3p4*pI2
        return
	end
c---------------------------------------------------------------------
c--------- sup sdown interference
c--------------------------------------------------------------------
      double precision function ampcharsupsdownint(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2l2, g2r1, g2r2, g3l1, g3l2, g3r1, g3r2, m1, m3,
     .  m4, amcha1, amcha2, p1p2, 
     . p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
      implicit none
      double precision g1L1, g1L2, g2l1, g2l2, g3l1,g3l2, g3r1, g3r2,
     . g2r2, g2r1,g1r1, 
     .  m1, m3, m4, amcha1, amcha2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, mI1, mI2
      double precision sdgf,amz,amw,pi,g2 
    
      mI1=amcha1
      mI2=amcha2
c---- g1-coupling of first vertex (stop-quark-chargino, 
c--- distinguing between different charginos necessary), 
c-----g2-coupling of second vertex (chargino-squark-quark)
c---- g3 coupling of thrid vertex (squark quark-neutralino)
c--- m1, p1-> quark
c--- m2, p2-> fermion
c--- m3, p3-> antifermion
c--- m4, p4 -> neutralino

      ampcharsupsdownint=
     -  -2*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*mI1*mI2*p1p2 - 
     -  2*g1L1*g1L2*g2L1*g2L2*g3L1*g3R2*m3*m4*mI1*mI2*p1p2 + 
     -  2*g1L1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m4*mI1*mI2*p2p3 + 
     -  2*g1L2*g1R1*g2L2*g2R1*g3L1*g3R2*m1*m4*mI1*mI2*p2p3 - 
     -  2*g1L1*g1L2*g2L1*g2L2*g3L1*g3L2*mI1*mI2*p1p4*p2p3 - 
     -  2*g1R1*g1R2*g2R1*g2R2*g3R1*g3R2*mI1*mI2*p1p4*p2p3 - 
     -  2*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m4*mI1*p1pI*p2p3 - 
     -  2*g1R1*g1R2*g2L2*g2R1*g3L1*g3R2*m4*mI1*p1pI*p2p3 - 
     -  2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m4*mI2*p1pI*p2p3 - 
     -  2*g1L1*g1L2*g2L2*g2R1*g3L1*g3R2*m4*mI2*p1pI*p2p3 - 
     -  2*g1L2*g1R1*g2L2*g2R1*g3L1*g3L2*m1*m3*mI1*mI2*p2p4 - 
     -  2*g1L1*g1R2*g2L1*g2R2*g3R1*g3R2*m1*m3*mI1*mI2*p2p4 + 
     -  2*g1L1*g1L2*g2L1*g2L2*g3L1*g3L2*mI1*mI2*p1p3*p2p4 + 
     -  2*g1R1*g1R2*g2R1*g2R2*g3R1*g3R2*mI1*mI2*p1p3*p2p4 + 
     -  2*g1R1*g1R2*g2L2*g2R1*g3L1*g3L2*m3*mI1*p1pI*p2p4 + 
     -  2*g1L1*g1L2*g2L1*g2R2*g3R1*g3R2*m3*mI1*p1pI*p2p4 + 
     -  2*g1L1*g1L2*g2L2*g2R1*g3L1*g3L2*m3*mI2*p1pI*p2p4 + 
     -  2*g1R1*g1R2*g2L1*g2R2*g3R1*g3R2*m3*mI2*p1pI*p2p4 + 
     -  2*g1L2*g1R1*g2R1*g2R2*g3L2*g3R1*m1*m3*m4*mI1*p2pI + 
     -  2*g1L1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m3*m4*mI1*p2pI + 
     -  2*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m1*m3*m4*mI2*p2pI + 
     -  2*g1L2*g1R1*g2L1*g2L2*g3L1*g3R2*m1*m3*m4*mI2*p2pI - 
     -  2*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m4*mI1*p1p3*p2pI - 
     -  2*g1R1*g1R2*g2L2*g2R1*g3L1*g3R2*m4*mI1*p1p3*p2pI + 
     -  2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m4*mI2*p1p3*p2pI + 
     -  2*g1L1*g1L2*g2L2*g2R1*g3L1*g3R2*m4*mI2*p1p3*p2pI + 
     -  2*g1R1*g1R2*g2L2*g2R1*g3L1*g3L2*m3*mI1*p1p4*p2pI + 
     -  2*g1L1*g1L2*g2L1*g2R2*g3R1*g3R2*m3*mI1*p1p4*p2pI - 
     -  2*g1L1*g1L2*g2L2*g2R1*g3L1*g3L2*m3*mI2*p1p4*p2pI - 
     -  2*g1R1*g1R2*g2L1*g2R2*g3R1*g3R2*m3*mI2*p1p4*p2pI - 
     -  4*g1L1*g1L2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1pI*p2pI - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L1*g3R2*m3*m4*p1pI*p2pI + 
     -  2*g1L1*g1L2*g2L1*g2L2*g3L1*g3L2*mI1*mI2*p1p2*p3p4 + 
     -  2*g1R1*g1R2*g2R1*g2R2*g3R1*g3R2*mI1*mI2*p1p2*p3p4 - 
     -  2*g1L1*g1R2*g2L1*g2L2*g3L1*g3L2*m1*mI1*p2pI*p3p4 - 
     -  2*g1L2*g1R1*g2R1*g2R2*g3R1*g3R2*m1*mI1*p2pI*p3p4 - 
     -  2*g1L2*g1R1*g2L1*g2L2*g3L1*g3L2*m1*mI2*p2pI*p3p4 - 
     -  2*g1L1*g1R2*g2R1*g2R2*g3R1*g3R2*m1*mI2*p2pI*p3p4 + 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L1*g3L2*p1pI*p2pI*p3p4 + 
     -  4*g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p1pI*p2pI*p3p4 + 
     -  2*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m4*mI1*p1p2*p3pI + 
     -  2*g1R1*g1R2*g2L2*g2R1*g3L1*g3R2*m4*mI1*p1p2*p3pI - 
     -  2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m4*mI2*p1p2*p3pI - 
     -  2*g1L1*g1L2*g2L2*g2R1*g3L1*g3R2*m4*mI2*p1p2*p3pI - 
     -  2*g1L1*g1R2*g2L1*g2L2*g3L1*g3L2*m1*mI1*p2p4*p3pI - 
     -  2*g1L2*g1R1*g2R1*g2R2*g3R1*g3R2*m1*mI1*p2p4*p3pI - 
     -  2*g1L2*g1R1*g2L1*g2L2*g3L1*g3L2*m1*mI2*p2p4*p3pI - 
     -  2*g1L1*g1R2*g2R1*g2R2*g3R1*g3R2*m1*mI2*p2p4*p3pI + 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L1*g3L2*p1pI*p2p4*p3pI + 
     -  4*g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p1pI*p2p4*p3pI - 
     -  2*g1R1*g1R2*g2L2*g2R1*g3L1*g3L2*m3*mI1*p1p2*p4pI - 
     -  2*g1L1*g1L2*g2L1*g2R2*g3R1*g3R2*m3*mI1*p1p2*p4pI + 
     -  2*g1L1*g1L2*g2L2*g2R1*g3L1*g3L2*m3*mI2*p1p2*p4pI + 
     -  2*g1R1*g1R2*g2L1*g2R2*g3R1*g3R2*m3*mI2*p1p2*p4pI + 
     -  2*g1L1*g1R2*g2L1*g2L2*g3L1*g3L2*m1*mI1*p2p3*p4pI + 
     -  2*g1L2*g1R1*g2R1*g2R2*g3R1*g3R2*m1*mI1*p2p3*p4pI + 
     -  2*g1L2*g1R1*g2L1*g2L2*g3L1*g3L2*m1*mI2*p2p3*p4pI + 
     -  2*g1L1*g1R2*g2R1*g2R2*g3R1*g3R2*m1*mI2*p2p3*p4pI - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L1*g3L2*p1pI*p2p3*p4pI - 
     -  4*g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p1pI*p2p3*p4pI + 
     -  2*g1L1*g1L2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1p2*pI2 + 
     -  2*g1R1*g1R2*g2L1*g2L2*g3L1*g3R2*m3*m4*p1p2*pI2 + 
     -  2*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m1*m4*p2p3*pI2 + 
     -  2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m1*m4*p2p3*pI2 + 
     -  2*g1R1*g1R2*g2L1*g2L2*g3L1*g3L2*p1p4*p2p3*pI2 + 
     -  2*g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p1p4*p2p3*pI2 - 
     -  2*g1L1*g1R2*g2L2*g2R1*g3L1*g3L2*m1*m3*p2p4*pI2 - 
     -  2*g1L2*g1R1*g2L1*g2R2*g3R1*g3R2*m1*m3*p2p4*pI2 - 
     -  2*g1R1*g1R2*g2L1*g2L2*g3L1*g3L2*p1p3*p2p4*pI2 - 
     -  2*g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p1p3*p2p4*pI2 - 
     -  2*g1R1*g1R2*g2L1*g2L2*g3L1*g3L2*p1p2*p3p4*pI2 - 
     -  2*g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p1p2*p3p4*pI2
        
        return
	end



c------------------------------------------------------------------
c------ sneutrino interchange amplitude
c--------------------------------------------------------------------
      double precision function ampcharsneutrino(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2, g3l1, g3l2, g3r1,g3r2, 
     .	m1, m3, m4, amcha1,
     .  amcha2, p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
      implicit none
      double precision g1L1, g1L2, g2l1, g2l2, g3r1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, amcha1, amcha2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, g3r2, g3l1, g3l2
      double precision sdgf,amz,amw,pi,g2, mI1, mI2 
      mI1=amcha1
      mI2=amcha2
c---- g1-coupling of first vertex (stop-quark-chargino, 
c--- distinguing between different charginos necessary), 
c-----g2-coupling of second vertex (chargino-sneutrino-lepton)
c---- g3 coupling of thrid vertex (sneutrino-neutrino-neutralino)
c--- m1, p1-> quark
c--- m2, p2-> fermion
c--- m3, p3-> antifermion
c--- m4, p4 -> neutralino
      ampcharsneutrino=

     -  -4*g1L2*g1R1*g2L2*g2R1*g3R1**2*m1*m3*mI1*mI2*p2p4 - 
     -  4*g1L1*g1R2*g2L1*g2R2*g3R1**2*m1*m3*mI1*mI2*p2p4 + 
     -  4*g1L1*g1L2*g2L1*g2L2*g3R1**2*mI1*mI2*p1p3*p2p4 + 
     -  4*g1R1*g1R2*g2R1*g2R2*g3R1**2*mI1*mI2*p1p3*p2p4 + 
     -  4*g1R1*g1R2*g2L2*g2R1*g3R1**2*m3*mI1*p1pI*p2p4 + 
     -  4*g1L1*g1L2*g2L1*g2R2*g3R1**2*m3*mI1*p1pI*p2p4 + 
     -  4*g1L1*g1L2*g2L2*g2R1*g3R1**2*m3*mI2*p1pI*p2p4 + 
     -  4*g1R1*g1R2*g2L1*g2R2*g3R1**2*m3*mI2*p1pI*p2p4 - 
     -  4*g1L1*g1R2*g2L1*g2L2*g3R1**2*m1*mI1*p2p4*p3pI - 
     -  4*g1L2*g1R1*g2R1*g2R2*g3R1**2*m1*mI1*p2p4*p3pI - 
     -  4*g1L2*g1R1*g2L1*g2L2*g3R1**2*m1*mI2*p2p4*p3pI - 
     -  4*g1L1*g1R2*g2R1*g2R2*g3R1**2*m1*mI2*p2p4*p3pI + 
     -  8*g1R1*g1R2*g2L1*g2L2*g3R1**2*p1pI*p2p4*p3pI + 
     -  8*g1L1*g1L2*g2R1*g2R2*g3R1**2*p1pI*p2p4*p3pI - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3R1**2*m1*m3*p2p4*pI2 - 
     -  4*g1L2*g1R1*g2L1*g2R2*g3R1**2*m1*m3*p2p4*pI2 - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3R1**2*p1p3*p2p4*pI2 - 
     -  4*g1L1*g1L2*g2R1*g2R2*g3R1**2*p1p3*p2p4*pI2


       
       return
        end   
c------------------------------------------------------------------
c------ sneutrino interchange amplitude
c--------------------------------------------------------------------
      double precision function ampcharsup(g1L1, g1L2, g1r1, 
     . g1r2, g2l1, g2r1, g2l2, g2r2, g3l1, g3l2, g3r1,g3r2, 
     .	m1, m3, m4, amcha1,
     .  amcha2, p1p2,p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
      implicit none
      double precision g1L1, g1L2, g2l1, g2l2, g3r1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, amcha1, amcha2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, g3r2, g3l1, g3l2
      double precision sdgf,amz,amw,pi,g2, mI1, mI2 
      mI1=amcha1
      mI2=amcha2

      ampcharsup=
     - 16d0*((g3L1*g3L2*p2p4)/2. + (g3R1*g3R2*p2p4)/2.)*
     -  (-(g1L2*g1R1*g2L2*g2R1*m1*m3*mI1*mI2)/2. - 
     -    (g1L1*g1R2*g2L1*g2R2*m1*m3*mI1*mI2)/2. + 
     -    (g1L1*g1L2*g2L1*g2L2*mI1*mI2*p1p3)/2. + 
     -    (g1R1*g1R2*g2R1*g2R2*mI1*mI2*p1p3)/2. + 
     -    (g1R1*g1R2*g2L2*g2R1*m3*mI1*p1pI)/2. + 
     -    (g1L1*g1L2*g2L1*g2R2*m3*mI1*p1pI)/2. + 
     -    (g1L1*g1L2*g2L2*g2R1*m3*mI2*p1pI)/2. + 
     -    (g1R1*g1R2*g2L1*g2R2*m3*mI2*p1pI)/2. - 
     -    (g1L1*g1R2*g2L1*g2L2*m1*mI1*p3pI)/2. - 
     -    (g1L2*g1R1*g2R1*g2R2*m1*mI1*p3pI)/2. - 
     -    (g1L2*g1R1*g2L1*g2L2*m1*mI2*p3pI)/2. - 
     -    (g1L1*g1R2*g2R1*g2R2*m1*mI2*p3pI)/2. + 
     -    g1R1*g1R2*g2L1*g2L2*p1pI*p3pI +
     -     g1L1*g1L2*g2R1*g2R2*p1pI*p3pI - 
     -    (g1L1*g1R2*g2L2*g2R1*m1*m3*pI2)/2. - 
     -    (g1L2*g1R1*g2L1*g2R2*m1*m3*pI2)/2. - 
     -    (g1R1*g1R2*g2L1*g2L2*p1p3*pI2)/2. - 
     -    (g1L1*g1L2*g2R1*g2R2*p1p3*pI2)/2.)

     
       
       return
        end   
c---------------------------------------------------
c------ W-amplitude one mass zero
c---------------------------------------------------
      double precision function ampwchar(g1L1, g1L2, g2l1, g2l2, 
     . g1r1, g1r2, g2r1, g2r2,g3l,
     .  m1, m3, m4, amcha1, amcha2, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
      implicit none
      double precision g1L1, g1L2, g2l1, g2l2, g3l,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, amcha1, amcha2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4
      double precision sdgf,amz,amw,pi,g2, mI1, mI2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2
c---- g1-coupling of first vertex (stop-quark-chargino, 
c--- distinguing between different charginos necessary), 
c-----g2-coupling of second vertex (chargino-neutralino-W)
c---- g3 coupling of thrid vertex (W-fermion -antifermion)
c--- m1, p1-> quark
c--- m2, p2-> fermion
c--- m3, p3-> antifermion
c--- m4, p4 -> neutralino
      mI1=amcha1
      mI2=amcha2
      
      
      ampwchar=   8*g1L2*g1R1*g2L2*g2R1*g3L**2*m1*m4*mI1*mI2*p2p3 + 
     -  8*g1L1*g1R2*g2L1*g2R2*g3L**2*m1*m4*mI1*mI2*p2p3 + 
     -  (8*g1L2*g1R1*g2L2*g2R1*g3L**2*m1*m3**2*m4*mI1*mI2*p2p3)/
     -   amw**2 + (8*g1L1*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*m4*mI1*
     -     mI2*p2p3)/amw**2 - 
     -  (4*g1L2*g1R1*g2L2*g2R1*g3L**2*m1*m3**4*m4*mI1*mI2*p2p3)/
     -   amw**4 - (4*g1L1*g1R2*g2L1*g2R2*g3L**2*m1*m3**4*m4*mI1*
     -     mI2*p2p3)/amw**4 + 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mI1*mI2*p1p4*p2p3)/
     -   amw**2 + (8*g1R1*g1R2*g2R1*g2R2*g3L**2*m3**2*mI1*mI2*
     -     p1p4*p2p3)/amw**2 - 
     -  (4*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**4*mI1*mI2*p1p4*p2p3)/
     -   amw**4 - (4*g1R1*g1R2*g2R1*g2R2*g3L**2*m3**4*mI1*mI2*
     -     p1p4*p2p3)/amw**4 - 
     -  8*g1R1*g1R2*g2L2*g2R1*g3L**2*m4*mI1*p1pI*p2p3 - 
     -  8*g1L1*g1L2*g2L1*g2R2*g3L**2*m4*mI1*p1pI*p2p3 - 
     -  (8*g1R1*g1R2*g2L2*g2R1*g3L**2*m3**2*m4*mI1*p1pI*p2p3)/
     -   amw**2 - (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**2*m4*mI1*
     -     p1pI*p2p3)/amw**2 + 
     -  (4*g1R1*g1R2*g2L2*g2R1*g3L**2*m3**4*m4*mI1*p1pI*p2p3)/
     -   amw**4 + (4*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**4*m4*mI1*
     -     p1pI*p2p3)/amw**4 - 
     -  8*g1L1*g1L2*g2L2*g2R1*g3L**2*m4*mI2*p1pI*p2p3 - 
     -  8*g1R1*g1R2*g2L1*g2R2*g3L**2*m4*mI2*p1pI*p2p3 - 
     -  (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mI2*p1pI*p2p3)/
     -   amw**2 - (8*g1R1*g1R2*g2L1*g2R2*g3L**2*m3**2*m4*mI2*
     -     p1pI*p2p3)/amw**2 + 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**4*m4*mI2*p1pI*p2p3)/
     -   amw**4 + (4*g1R1*g1R2*g2L1*g2R2*g3L**2*m3**4*m4*mI2*
     -     p1pI*p2p3)/amw**4 - 
     -  (8*g1L2*g1R1*g2L2*g2R1*g3L**2*m1*m3**2*m4*mI1*mI2*
     -     p2p3**2)/amw**4 - 
     -  (8*g1L1*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*m4*mI1*mI2*
     -     p2p3**2)/amw**4 - 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mI1*mI2*p1p4*
     -     p2p3**2)/amw**4 - 
     -  (8*g1R1*g1R2*g2R1*g2R2*g3L**2*m3**2*mI1*mI2*p1p4*
     -     p2p3**2)/amw**4 + 
     -  (8*g1R1*g1R2*g2L2*g2R1*g3L**2*m3**2*m4*mI1*p1pI*
     -     p2p3**2)/amw**4 + 
     -  (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**2*m4*mI1*p1pI*
     -     p2p3**2)/amw**4 + 
     -  (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mI2*p1pI*
     -     p2p3**2)/amw**4 + 
     -  (8*g1R1*g1R2*g2L1*g2R2*g3L**2*m3**2*m4*mI2*p1pI*
     -     p2p3**2)/amw**4 - 
     -  (16*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mI1*mI2*p1p2*p2p4)/
     -   amw**2 - (16*g1R1*g1R2*g2R1*g2R2*g3L**2*m3**2*mI1*mI2*
     -     p1p2*p2p4)/amw**2 + 
     -  8*g1L1*g1L2*g2L1*g2L2*g3L**2*mI1*mI2*p1p3*p2p4 + 
     -  8*g1R1*g1R2*g2R1*g2R2*g3L**2*mI1*mI2*p1p3*p2p4 - 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mI1*mI2*p1p3*p2p4)/
     -   amw**2 - (8*g1R1*g1R2*g2R1*g2R2*g3L**2*m3**2*mI1*mI2*
     -     p1p3*p2p4)/amw**2 + 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mI1*mI2*p1p2*p2p3*
     -     p2p4)/amw**4 + 
     -  (8*g1R1*g1R2*g2R1*g2R2*g3L**2*m3**2*mI1*mI2*p1p2*p2p3*
     -     p2p4)/amw**4 + 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mI1*mI2*p1p3*p2p3*
     -     p2p4)/amw**4 + 
     -  (8*g1R1*g1R2*g2R1*g2R2*g3L**2*m3**2*mI1*mI2*p1p3*p2p3*
     -     p2p4)/amw**4 + 
     -  (16*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*mI1*p2p4*p2pI)/
     -   amw**2 + (16*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*m3**2*mI1*
     -     p2p4*p2pI)/amw**2 + 
     -  (16*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*m3**2*mI2*p2p4*p2pI)/
     -   amw**2 + (16*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*m3**2*mI2*
     -     p2p4*p2pI)/amw**2 - 
     -  (32*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1pI*p2p4*p2pI)/
     -   amw**2 - (32*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1pI*
     -     p2p4*p2pI)/amw**2 - 
     -  (8*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*mI1*p2p3*p2p4*
     -     p2pI)/amw**4 - 
     -  (8*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*m3**2*mI1*p2p3*p2p4*
     -     p2pI)/amw**4 - 
     -  (8*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*m3**2*mI2*p2p3*p2p4*
     -     p2pI)/amw**4 - 
     -  (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*m3**2*mI2*p2p3*p2p4*
     -     p2pI)/amw**4 + 
     -  (16*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1pI*p2p3*p2p4*
     -     p2pI)/amw**4 + 
     -  (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1pI*p2p3*p2p4*
     -     p2pI)/amw**4 + 
     -  8*g1L1*g1L2*g2L1*g2L2*g3L**2*mI1*mI2*p1p2*p3p4 + 
     -  8*g1R1*g1R2*g2R1*g2R2*g3L**2*mI1*mI2*p1p2*p3p4 - 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mI1*mI2*p1p2*p3p4)/
     -   amw**2 - (8*g1R1*g1R2*g2R1*g2R2*g3L**2*m3**2*mI1*mI2*
     -     p1p2*p3p4)/amw**2 + 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mI1*mI2*p1p2*p2p3*
     -     p3p4)/amw**4 + 
     -  (8*g1R1*g1R2*g2R1*g2R2*g3L**2*m3**2*mI1*mI2*p1p2*p2p3*
     -     p3p4)/amw**4 + 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mI1*mI2*p1p3*p2p3*
     -     p3p4)/amw**4 + 
     -  (8*g1R1*g1R2*g2R1*g2R2*g3L**2*m3**2*mI1*mI2*p1p3*p2p3*
     -     p3p4)/amw**4 - 
     -  8*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*mI1*p2pI*p3p4 - 
     -  8*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*mI1*p2pI*p3p4 + 
     -  (8*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*mI1*p2pI*p3p4)/
     -   amw**2 + (8*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*m3**2*mI1*
     -     p2pI*p3p4)/amw**2 - 
     -  8*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*mI2*p2pI*p3p4 - 
     -  8*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*mI2*p2pI*p3p4 + 
     -  (8*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*m3**2*mI2*p2pI*p3p4)/
     -   amw**2 + (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*m3**2*mI2*
     -     p2pI*p3p4)/amw**2 + 
     -  16*g1R1*g1R2*g2L1*g2L2*g3L**2*p1pI*p2pI*p3p4 + 
     -  16*g1L1*g1L2*g2R1*g2R2*g3L**2*p1pI*p2pI*p3p4 - 
     -  (16*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1pI*p2pI*p3p4)/
     -   amw**2 - (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1pI*
     -     p2pI*p3p4)/amw**2 - 
     -  (8*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*mI1*p2p3*p2pI*
     -     p3p4)/amw**4 - 
     -  (8*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*m3**2*mI1*p2p3*p2pI*
     -     p3p4)/amw**4 - 
     -  (8*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*m3**2*mI2*p2p3*p2pI*
     -     p3p4)/amw**4 - 
     -  (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*m3**2*mI2*p2p3*p2pI*
     -     p3p4)/amw**4 + 
     -  (16*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1pI*p2p3*p2pI*
     -     p3p4)/amw**4 + 
     -  (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1pI*p2p3*p2pI*
     -     p3p4)/amw**4 - 
     -  8*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*mI1*p2p4*p3pI - 
     -  8*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*mI1*p2p4*p3pI + 
     -  (8*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*mI1*p2p4*p3pI)/
     -   amw**2 + (8*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*m3**2*mI1*
     -     p2p4*p3pI)/amw**2 - 
     -  8*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*mI2*p2p4*p3pI - 
     -  8*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*mI2*p2p4*p3pI + 
     -  (8*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*m3**2*mI2*p2p4*p3pI)/
     -   amw**2 + (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*m3**2*mI2*
     -     p2p4*p3pI)/amw**2 + 
     -  16*g1R1*g1R2*g2L1*g2L2*g3L**2*p1pI*p2p4*p3pI + 
     -  16*g1L1*g1L2*g2R1*g2R2*g3L**2*p1pI*p2p4*p3pI - 
     -  (16*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1pI*p2p4*p3pI)/
     -   amw**2 - (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1pI*
     -     p2p4*p3pI)/amw**2 - 
     -  (8*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*mI1*p2p3*p2p4*
     -     p3pI)/amw**4 - 
     -  (8*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*m3**2*mI1*p2p3*p2p4*
     -     p3pI)/amw**4 - 
     -  (8*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*m3**2*mI2*p2p3*p2p4*
     -     p3pI)/amw**4 - 
     -  (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*m3**2*mI2*p2p3*p2p4*
     -     p3pI)/amw**4 + 
     -  (16*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1pI*p2p3*p2p4*
     -     p3pI)/amw**4 + 
     -  (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1pI*p2p3*p2p4*
     -     p3pI)/amw**4 - 
     -  (8*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*mI1*p2p3*p3p4*
     -     p3pI)/amw**4 - 
     -  (8*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*m3**2*mI1*p2p3*p3p4*
     -     p3pI)/amw**4 - 
     -  (8*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*m3**2*mI2*p2p3*p3p4*
     -     p3pI)/amw**4 - 
     -  (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*m3**2*mI2*p2p3*p3p4*
     -     p3pI)/amw**4 + 
     -  (16*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1pI*p2p3*p3p4*
     -     p3pI)/amw**4 + 
     -  (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1pI*p2p3*p3p4*
     -     p3pI)/amw**4 - 
     -  (8*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*mI1*p2p3*p4pI)/
     -   amw**2 - (8*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*m3**2*mI1*
     -     p2p3*p4pI)/amw**2 + 
     -  (4*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*m3**4*mI1*p2p3*p4pI)/
     -   amw**4 + (4*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*m3**4*mI1*
     -     p2p3*p4pI)/amw**4 - 
     -  (8*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*m3**2*mI2*p2p3*p4pI)/
     -   amw**2 - (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*m3**2*mI2*
     -     p2p3*p4pI)/amw**2 + 
     -  (4*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*m3**4*mI2*p2p3*p4pI)/
     -   amw**4 + (4*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*m3**4*mI2*
     -     p2p3*p4pI)/amw**4 + 
     -  (16*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1pI*p2p3*p4pI)/
     -   amw**2 + (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1pI*
     -     p2p3*p4pI)/amw**2 - 
     -  (8*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**4*p1pI*p2p3*p4pI)/
     -   amw**4 - (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**4*p1pI*p2p3*
     -     p4pI)/amw**4 + 
     -  (8*g1L1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*mI1*p2p3**2*
     -     p4pI)/amw**4 + 
     -  (8*g1L2*g1R1*g2R1*g2R2*g3L**2*m1*m3**2*mI1*p2p3**2*
     -     p4pI)/amw**4 + 
     -  (8*g1L2*g1R1*g2L1*g2L2*g3L**2*m1*m3**2*mI2*p2p3**2*
     -     p4pI)/amw**4 + 
     -  (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m1*m3**2*mI2*p2p3**2*
     -     p4pI)/amw**4 - 
     -  (16*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1pI*p2p3**2*p4pI)/
     -   amw**4 - (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1pI*
     -     p2p3**2*p4pI)/amw**4 + 
     -  8*g1L1*g1R2*g2L2*g2R1*g3L**2*m1*m4*p2p3*pI2 + 
     -  8*g1L2*g1R1*g2L1*g2R2*g3L**2*m1*m4*p2p3*pI2 + 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p3*pI2)/
     -   amw**2 + (8*g1L2*g1R1*g2L1*g2R2*g3L**2*m1*m3**2*m4*
     -     p2p3*pI2)/amw**2 - 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L**2*m1*m3**4*m4*p2p3*pI2)/
     -   amw**4 - (4*g1L2*g1R1*g2L1*g2R2*g3L**2*m1*m3**4*m4*
     -     p2p3*pI2)/amw**4 - 
     -  (8*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p3*pI2)/
     -   amw**2 - (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3*
     -     pI2)/amw**2 + (4*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**4*
     -     p1p4*p2p3*pI2)/amw**4 + 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**4*p1p4*p2p3*pI2)/
     -   amw**4 - (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m1*m3**2*m4*
     -     p2p3**2*pI2)/amw**4 - 
     -  (8*g1L2*g1R1*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p3**2*pI2)/
     -   amw**4 + (8*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1p4*
     -     p2p3**2*pI2)/amw**4 + 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3**2*pI2)/
     -   amw**4 + (16*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1p2*
     -     p2p4*pI2)/amw**2 + 
     -  (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*p2p4*pI2)/
     -   amw**2 - 8*g1R1*g1R2*g2L1*g2L2*g3L**2*p1p3*p2p4*pI2 - 
     -  8*g1L1*g1L2*g2R1*g2R2*g3L**2*p1p3*p2p4*pI2 + 
     -  (8*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1p3*p2p4*pI2)/
     -   amw**2 + (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p3*p2p4*
     -     pI2)/amw**2 - (8*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*
     -     p1p2*p2p3*p2p4*pI2)/amw**4 - 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*p2p3*p2p4*pI2)/
     -   amw**4 - (8*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1p3*p2p3*
     -     p2p4*pI2)/amw**4 - 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p3*p2p3*p2p4*pI2)/
     -   amw**4 - 8*g1R1*g1R2*g2L1*g2L2*g3L**2*p1p2*p3p4*pI2 - 
     -  8*g1L1*g1L2*g2R1*g2R2*g3L**2*p1p2*p3p4*pI2 + 
     -  (8*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1p2*p3p4*pI2)/
     -   amw**2 + (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*p3p4*
     -     pI2)/amw**2 - (8*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*
     -     p1p2*p2p3*p3p4*pI2)/amw**4 - 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*p2p3*p3p4*pI2)/
     -   amw**4 - (8*g1R1*g1R2*g2L1*g2L2*g3L**2*m3**2*p1p3*p2p3*
     -     p3p4*pI2)/amw**4 - 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p3*p2p3*p3p4*pI2)/
     -   amw**4

       
	return
        end

c--------------------------------------------------------------
c---- amplitude squared to SM fermion exchange
c---- g1 is the coupling bottom-top-W
c---- g3 is the coupling fermion-fermion W

      double precision function ampSMfer(g1L1, g1L2, g2l1, g2l2, 
     . g2r1, g2r2,g3l,
     .  m1, m3, m4, mf1, mf2, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2)
      implicit none
      double precision g1L1, g1L2, g2l1, g2l2, g3l,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mf1, mf2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

     

        ampSMfer=
     - (-4*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**4*mf1*mf2*p1p4*p2p3)/
     -   amw**4 + (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mf1*mf2*
     -     p1p4*p2p3)/amw**2 + 
     -  (4*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**4*m4*mf1*p1pI*p2p3)/
     -   amw**4 + (4*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**4*m4*mf2*
     -     p1pI*p2p3)/amw**4 - 
     -  (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**2*m4*mf1*p1pI*p2p3)/
     -   amw**2 - (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mf2*
     -     p1pI*p2p3)/amw**2 - 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mf1*mf2*p1p4*
     -     p2p3**2)/amw**4 + 
     -  (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**2*m4*mf1*p1pI*
     -     p2p3**2)/amw**4 + 
     -  (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mf2*p1pI*
     -     p2p3**2)/amw**4 - 
     -  (16*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mf1*mf2*p1p2*p2p4)/
     -   amw**2 + 8*g1L1*g1L2*g2L1*g2L2*g3L**2*mf1*mf2*p1p3*
     -   p2p4 - (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mf1*mf2*
     -     p1p3*p2p4)/amw**2 + 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mf1*mf2*p1p2*p2p3*
     -     p2p4)/amw**4 + (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*
     -     mf1*mf2*p1p3*p2p3*p2p4)/amw**4 + 
     -  (16*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**2*m4*mf1*p1p2*p2pI)/
     -   amw**2 + (16*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mf2*
     -     p1p2*p2pI)/amw**2 - 
     -  8*g1L1*g1L2*g2L1*g2R2*g3L**2*m4*mf1*p1p3*p2pI - 
     -  8*g1L1*g1L2*g2L2*g2R1*g3L**2*m4*mf2*p1p3*p2pI + 
     -  (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**2*m4*mf1*p1p3*p2pI)/
     -   amw**2 + (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mf2*
     -     p1p3*p2pI)/amw**2 - 
     -  (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**2*m4*mf1*p1p2*p2p3*
     -     p2pI)/amw**4 - (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*
     -     mf2*p1p2*p2p3*p2pI)/amw**4 - 
     -  (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**2*m4*mf1*p1p3*p2p3*
     -     p2pI)/amw**4 - (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*
     -     mf2*p1p3*p2p3*p2pI)/amw**4 + 
     -  8*g1L1*g1L2*g2L1*g2L2*g3L**2*mf1*mf2*p1p2*p3p4 - 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mf1*mf2*p1p2*p3p4)/
     -   amw**2 + (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mf1*mf2*
     -     p1p2*p2p3*p3p4)/amw**4 + 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m3**2*mf1*mf2*p1p3*p2p3*
     -     p3p4)/amw**4 - 8*g1L1*g1L2*g2L1*g2R2*g3L**2*m4*mf1*
     -   p1p2*p3pI - 8*g1L1*g1L2*g2L2*g2R1*g3L**2*m4*mf2*p1p2*
     -   p3pI + (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**2*m4*mf1*p1p2*
     -     p3pI)/amw**2 + (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*
     -     mf2*p1p2*p3pI)/amw**2 - 
     -  (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**2*m4*mf1*p1p2*p2p3*
     -     p3pI)/amw**4 - (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*
     -     mf2*p1p2*p2p3*p3pI)/amw**4 - 
     -  (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m3**2*m4*mf1*p1p3*p2p3*
     -     p3pI)/amw**4 - (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*
     -     mf2*p1p3*p2p3*p3pI)/amw**4 - 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**4*p1pI*p2p3*p4pI)/
     -   amw**4 + (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1pI*p2p3*
     -     p4pI)/amw**2 - (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*
     -     p1pI*p2p3**2*p4pI)/amw**4 - 
     -  (32*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*p2pI*p4pI)/
     -   amw**2 + 16*g1L1*g1L2*g2R1*g2R2*g3L**2*p1p3*p2pI*p4pI - 
     -  (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p3*p2pI*p4pI)/
     -   amw**2 + (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*p2p3*
     -     p2pI*p4pI)/amw**4 + 
     -  (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p3*p2p3*p2pI*
     -     p4pI)/amw**4 + 16*g1L1*g1L2*g2R1*g2R2*g3L**2*p1p2*
     -   p3pI*p4pI - (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*
     -     p3pI*p4pI)/amw**2 + 
     -  (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*p2p3*p3pI*
     -     p4pI)/amw**4 + (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*
     -     p1p3*p2p3*p3pI*p4pI)/amw**4 + 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**4*p1p4*p2p3*pI2)/
     -   amw**4 - (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3*
     -     pI2)/amw**2 + (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*
     -     p1p4*p2p3**2*pI2)/amw**4 + 
     -  (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*p2p4*pI2)/
     -   amw**2 - 8*g1L1*g1L2*g2R1*g2R2*g3L**2*p1p3*p2p4*pI2 + 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p3*p2p4*pI2)/
     -   amw**2 - (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*p2p3*
     -     p2p4*pI2)/amw**4 - 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p3*p2p3*p2p4*pI2)/
     -   amw**4 - 8*g1L1*g1L2*g2R1*g2R2*g3L**2*p1p2*p3p4*pI2 + 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*p3p4*pI2)/
     -   amw**2 - (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p2*p2p3*
     -     p3p4*pI2)/amw**4 - 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m3**2*p1p3*p2p3*p3p4*pI2)/
     -   amw**4
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c---- Interfernce between diagrams with W-Prop chargino Prop and top Prop
c---- g1L1 =W-b-top kopplung
c---- g1R2, g1L2 stop top neutralino
c---  g2l1, g2r1 neutralino stop top Kopplung
c---- g2l2, g2r2 chargino bottom stop Kopplung
c---- g3  W Kopplung fermionen
c---- pfer= ptop
c---- mI1=top mass mI2=chargino mass


       double precision function ampintWSMferWchar(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l,
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)
      implicit none
      double precision g1L1, g1L2, g2l1, g2l2, g3l,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4,pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintWSMferWchar=
     - (-8*g1L1*g1L2*g2L1*g2L2*g3L**2*m1*m4*mI1*mI2*p2p3 + 
     -  (4*g1L1*g1L2*g2L1*g2L2*g3L**2*m1*m3**4*m4*mI1*mI2*p2p3)/
     -   amw**4 - (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*m4*mI1*
     -     mI2*p2p3)/amw**2 + 
     -  (4*g1L1*g1R2*g2L1*g2R2*g3L**2*m3**4*mI1*mI2*p1p4*p2p3)/
     -   amw**4 - (8*g1L1*g1R2*g2L1*g2R2*g3L**2*m3**2*mI1*mI2*
     -     p1p4*p2p3)/amw**2 - 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L**2*m3**4*m4*mI2*p1pfer*p2p3)/
     -   amw**4 + (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m3**2*m4*mI2*
     -     p1pfer*p2p3)/amw**2 + 
     -  8*g1L1*g1R2*g2L1*g2L2*g3L**2*m4*mI1*p1pI*p2p3 - 
     -  (4*g1L1*g1R2*g2L1*g2L2*g3L**2*m3**4*m4*mI1*p1pI*p2p3)/
     -   amw**4 + (8*g1L1*g1R2*g2L1*g2L2*g3L**2*m3**2*m4*mI1*
     -     p1pI*p2p3)/amw**2 + 
     -  (8*g1L1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*m4*mI1*mI2*
     -     p2p3**2)/amw**4 + 
     -  (8*g1L1*g1R2*g2L1*g2R2*g3L**2*m3**2*mI1*mI2*p1p4*
     -     p2p3**2)/amw**4 - 
     -  (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m3**2*m4*mI2*p1pfer*
     -     p2p3**2)/amw**4 - 
     -  (8*g1L1*g1R2*g2L1*g2L2*g3L**2*m3**2*m4*mI1*p1pI*
     -     p2p3**2)/amw**4 + 
     -  (16*g1L1*g1R2*g2L1*g2R2*g3L**2*m3**2*mI1*mI2*p1p2*p2p4)/
     -   amw**2 - 8*g1L1*g1R2*g2L1*g2R2*g3L**2*mI1*mI2*p1p3*
     -   p2p4 + (8*g1L1*g1R2*g2L1*g2R2*g3L**2*m3**2*mI1*mI2*
     -     p1p3*p2p4)/amw**2 - 
     -  (8*g1L1*g1R2*g2L1*g2R2*g3L**2*m3**2*mI1*mI2*p1p2*p2p3*
     -     p2p4)/amw**4 - (8*g1L1*g1R2*g2L1*g2R2*g3L**2*m3**2*
     -     mI1*mI2*p1p3*p2p3*p2p4)/amw**4 - 
     -  (16*g1L1*g1R2*g2R1*g2R2*g3L**2*m3**2*m4*mI2*p1p2*
     -     p2pfer)/amw**2 + 
     -  8*g1L1*g1R2*g2R1*g2R2*g3L**2*m4*mI2*p1p3*p2pfer - 
     -  (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m3**2*m4*mI2*p1p3*p2pfer)/
     -   amw**2 + (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m3**2*m4*mI2*
     -     p1p2*p2p3*p2pfer)/amw**4 + 
     -  (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m3**2*m4*mI2*p1p3*p2p3*
     -     p2pfer)/amw**4 - 
     -  (16*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*m3**2*mI1*p2p4*p2pI)/
     -   amw**2 + (16*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*
     -     p2p4*p2pI)/amw**2 + 
     -  (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*m3**2*mI1*p2p3*p2p4*
     -     p2pI)/amw**4 - (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*
     -     p1pfer*p2p3*p2p4*p2pI)/amw**4 + 
     -  (16*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m3**2*m4*p2pfer*p2pI)/
     -   amw**2 - (16*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*
     -     p2pfer*p2pI)/amw**2 - 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m3**2*m4*p2p3*p2pfer*
     -     p2pI)/amw**4 + (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*
     -     p1p4*p2p3*p2pfer*p2pI)/amw**4 - 
     -  8*g1L1*g1R2*g2L1*g2R2*g3L**2*mI1*mI2*p1p2*p3p4 + 
     -  (8*g1L1*g1R2*g2L1*g2R2*g3L**2*m3**2*mI1*mI2*p1p2*p3p4)/
     -   amw**2 - (8*g1L1*g1R2*g2L1*g2R2*g3L**2*m3**2*mI1*mI2*
     -     p1p2*p2p3*p3p4)/amw**4 - 
     -  (8*g1L1*g1R2*g2L1*g2R2*g3L**2*m3**2*mI1*mI2*p1p3*p2p3*
     -     p3p4)/amw**4 + 8*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*mI1*
     -   p2pI*p3p4 - (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*m3**2*mI1*
     -     p2pI*p3p4)/amw**2 - 
     -  8*g1L1*g1R2*g2L2*g2R1*g3L**2*p1pfer*p2pI*p3p4 + 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*p2pI*p3p4)/
     -   amw**2 + (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*m3**2*mI1*
     -     p2p3*p2pI*p3p4)/amw**4 - 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*p2p3*p2pI*
     -     p3p4)/amw**4 + 8*g1L1*g1R2*g2R1*g2R2*g3L**2*m4*mI2*
     -   p1p2*p3pfer - (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m3**2*m4*
     -     mI2*p1p2*p3pfer)/amw**2 + 
     -  (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m3**2*m4*mI2*p1p2*p2p3*
     -     p3pfer)/amw**4 + 
     -  (8*g1L1*g1R2*g2R1*g2R2*g3L**2*m3**2*m4*mI2*p1p3*p2p3*
     -     p3pfer)/amw**4 - 
     -  8*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m4*p2pI*p3pfer + 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m3**2*m4*p2pI*p3pfer)/
     -   amw**2 + 8*g1L1*g1R2*g2L2*g2R1*g3L**2*p1p4*p2pI*
     -   p3pfer - (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2pI*
     -     p3pfer)/amw**2 - 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m3**2*m4*p2p3*p2pI*
     -     p3pfer)/amw**4 + 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p3*p2pI*
     -     p3pfer)/amw**4 + 
     -  8*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*mI1*p2p4*p3pI - 
     -  (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*m3**2*mI1*p2p4*p3pI)/
     -   amw**2 - 8*g1L1*g1R2*g2L2*g2R1*g3L**2*p1pfer*p2p4*
     -   p3pI + (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*p2p4*
     -     p3pI)/amw**2 + (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*m3**2*
     -     mI1*p2p3*p2p4*p3pI)/amw**4 - 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*p2p3*p2p4*
     -     p3pI)/amw**4 - 8*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m4*
     -   p2pfer*p3pI + (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m3**2*
     -     m4*p2pfer*p3pI)/amw**2 + 
     -  8*g1L1*g1R2*g2L2*g2R1*g3L**2*p1p4*p2pfer*p3pI - 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2pfer*p3pI)/
     -   amw**2 - (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m3**2*m4*p2p3*
     -     p2pfer*p3pI)/amw**4 + 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p3*p2pfer*
     -     p3pI)/amw**4 + (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*m3**2*
     -     mI1*p2p3*p3p4*p3pI)/amw**4 - 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*p2p3*p3p4*
     -     p3pI)/amw**4 - (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m3**2*
     -     m4*p2p3*p3pfer*p3pI)/amw**4 + 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p3*p3pfer*
     -     p3pI)/amw**4 + 8*g1L1*g1L2*g2L2*g2R1*g3L**2*m1*mI2*
     -   p2p3*p4pfer - (4*g1L1*g1L2*g2L2*g2R1*g3L**2*m1*m3**4*
     -     mI2*p2p3*p4pfer)/amw**4 + 
     -  (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m1*m3**2*mI2*p2p3*p4pfer)/
     -   amw**2 - 8*g1L1*g1R2*g2L2*g2R1*g3L**2*p1pI*p2p3*
     -   p4pfer + (4*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1pI*p2p3*
     -     p4pfer)/amw**4 - 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pI*p2p3*p4pfer)/
     -   amw**2 - (8*g1L1*g1L2*g2L2*g2R1*g3L**2*m1*m3**2*mI2*
     -     p2p3**2*p4pfer)/amw**4 + 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pI*p2p3**2*
     -     p4pfer)/amw**4 - 
     -  (4*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*m3**4*mI1*p2p3*p4pI)/
     -   amw**4 + (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*m3**2*mI1*
     -     p2p3*p4pI)/amw**2 + 
     -  8*g1L1*g1R2*g2L2*g2R1*g3L**2*p1pfer*p2p3*p4pI + 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1pfer*p2p3*p4pI)/
     -   amw**4 - (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*
     -     p2p3*p4pI)/amw**2 - 
     -  (8*g1L1*g1L2*g2L1*g2R2*g3L**2*m1*m3**2*mI1*p2p3**2*
     -     p4pI)/amw**4 + (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*
     -     p1pfer*p2p3**2*p4pI)/amw**4 + 
     -  (16*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p2pfer*p4pI)/
     -   amw**2 - 8*g1L1*g1R2*g2L2*g2R1*g3L**2*p1p3*p2pfer*
     -   p4pI + (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p2pfer*
     -     p4pI)/amw**2 - (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*
     -     p1p2*p2p3*p2pfer*p4pI)/amw**4 - 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p2p3*p2pfer*
     -     p4pI)/amw**4 - 8*g1L1*g1R2*g2L2*g2R1*g3L**2*p1p2*
     -   p3pfer*p4pI + (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*
     -     p3pfer*p4pI)/amw**2 - 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p2p3*p3pfer*
     -     p4pI)/amw**4 - (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*
     -     p1p3*p2p3*p3pfer*p4pI)/amw**4 + 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m3**4*m4*p2p3*pIpfer)/
     -   amw**4 - (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m3**2*m4*p2p3*
     -     pIpfer)/amw**2 - 
     -  8*g1L1*g1R2*g2L2*g2R1*g3L**2*p1p4*p2p3*pIpfer - 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1p4*p2p3*pIpfer)/
     -   amw**4 + (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p3*
     -     pIpfer)/amw**2 + 
     -  (8*g1L1*g1L2*g2R1*g2R2*g3L**2*m1*m3**2*m4*p2p3**2*
     -     pIpfer)/amw**4 - 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p3**2*
     -     pIpfer)/amw**4 - 
     -  (16*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p2p4*pIpfer)/
     -   amw**2 + 8*g1L1*g1R2*g2L2*g2R1*g3L**2*p1p3*p2p4*
     -   pIpfer - (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p2p4*
     -     pIpfer)/amw**2 + 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p2p3*p2p4*
     -     pIpfer)/amw**4 + 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p2p3*p2p4*
     -     pIpfer)/amw**4 + 
     -  8*g1L1*g1R2*g2L2*g2R1*g3L**2*p1p2*p3p4*pIpfer - 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p3p4*pIpfer)/
     -   amw**2 + (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p2p3*
     -     p3p4*pIpfer)/amw**4 + 
     -  (8*g1L1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p2p3*p3p4*
     -     pIpfer)/amw**4)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c---- Interfernce between diagrams with slepton and top Prop

       double precision function ampintWSMferslep(g1L1, g1L2, g1r2,
     .  g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, mfer, amcha, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)
      implicit none
      double precision g1L1, g1L2, g2l1, g2l2, g3r1, g3r2,g2r2, g2r1, 
     .g1r1,m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha,  pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer, g3l1, g3l2, mfer
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2
c     typo in matrixelement...
      m3=-m3

      ampintWSMferslep= 
     -   -(     8*g1L1*g1R2*g2L1*g2R2*g3L1*g3R2*m3*m4*mfer*amcha*p1p2 - 
     -  (2*g1L1*g1R2*g2L1*g2R2*g3L1*g3R2*m3**3*m4*mfer*amcha*
     -     p1p2)/amw**2 - 
     -  (4*g1L1*g1R2*g2L1*g2R2*g3L1*g3R2*m3*m4*mfer*amcha*p1p2*
     -     p2p3)/amw**2 + 
     -  (2*g1L1*g1R2*g2L1*g2R2*g3L1*g3L2*m3**2*mfer*amcha*p1p4*
     -     p2p3)/amw**2 - 
     -  (2*g1L1*g1R2*g2R1*g2R2*g3L1*g3L2*m3**2*m4*amcha*p1pfer*
     -     p2p3)/amw**2 - 
     -  (4*g1L1*g1R2*g2L1*g2R2*g3L1*g3L2*m3**2*mfer*amcha*p1p2*
     -     p2p4)/amw**2 - 
     -  (2*g1L1*g1R2*g2L1*g2R2*g3L1*g3L2*m3**2*mfer*amcha*p1p3*
     -     p2p4)/amw**2 + 
     -  (2*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*amcha*p1pfer*p2p4)/
     -   amw**2 + (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3L2*m3**2*m4*amcha*
     -     p1p2*p2pfer)/amw**2 + 
     -  (2*g1L1*g1R2*g2R1*g2R2*g3L1*g3L2*m3**2*m4*amcha*p1p3*
     -     p2pfer)/amw**2 - 
     -  (2*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*amcha*p1p4*p2pfer)/
     -   amw**2 - 8*g1L1*g1L2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mfer*
     -   p2pI + (2*g1L1*g1L2*g2L1*g2R2*g3L1*g3R2*m1*m3**3*m4*
     -     mfer*p2pI)/amw**2 + 
     -  (4*g1L1*g1L2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mfer*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*g1L1*g1L2*g2L1*g2R2*g3L1*g3L2*m1*m3**2*mfer*p2p4*
     -     p2pI)/amw**2 - 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3**2*m4*p2pfer*
     -     p2pI)/amw**2 + 
     -  8*g1L1*g1R2*g2L1*g2R2*g3L1*g3L2*mfer*amcha*p1p2*p3p4 - 
     -  (2*g1L1*g1R2*g2L1*g2R2*g3L1*g3L2*m3**2*mfer*amcha*p1p2*
     -     p3p4)/amw**2 - 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1pfer*p2p3*
     -     p3p4)/amw**2 + 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2*p2pfer*
     -     p3p4)/amw**2 + 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p3*p2pfer*
     -     p3p4)/amw**2 - 
     -  8*g1L1*g1L2*g2L1*g2R2*g3L1*g3L2*m1*mfer*p2pI*p3p4 + 
     -  (2*g1L1*g1L2*g2L1*g2R2*g3L1*g3L2*m1*m3**2*mfer*p2pI*
     -     p3p4)/amw**2 - 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2pfer*p2pI*
     -     p3p4)/amw**2 - 
     -  8*g1L1*g1R2*g2R1*g2R2*g3L1*g3L2*m4*amcha*p1p2*p3pfer + 
     -  (2*g1L1*g1R2*g2R1*g2R2*g3L1*g3L2*m3**2*m4*amcha*p1p2*
     -     p3pfer)/amw**2 + 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p4*p2p3*
     -     p3pfer)/amw**2 - 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2*p2p4*
     -     p3pfer)/amw**2 - 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p3*p2p4*
     -     p3pfer)/amw**2 + 
     -  8*g1L1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m4*p2pI*p3pfer - 
     -  (2*g1L1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3**2*m4*p2pI*
     -     p3pfer)/amw**2 + 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*p2pI*
     -     p3pfer)/amw**2 + 
     -  (2*g1L1*g1L2*g2L1*g2R2*g3L1*g3L2*m1*m3**2*mfer*p2p4*
     -     p3pI)/amw**2 - 
     -  (2*g1L1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3**2*m4*p2pfer*
     -     p3pI)/amw**2 - 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2pfer*p3p4*
     -     p3pI)/amw**2 + 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*p3pfer*
     -     p3pI)/amw**2 - 
     -  8*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2*p4pfer + 
     -  (2*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*amcha*p1p2*p4pfer)/
     -   amw**2 + (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2*
     -     p2p3*p4pfer)/amw**2 + 
     -  8*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2pI*p4pfer - 
     -  (2*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2pI*p4pfer)/
     -   amw**2 - (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*
     -     p2pI*p4pfer)/amw**2 - 
     -  (2*g1L1*g1L2*g2L1*g2R2*g3L1*g3L2*m1*m3**2*mfer*p2p3*
     -     p4pI)/amw**2 + 
     -  (2*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2pfer*p4pI)/
     -   amw**2 - (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*
     -     p3pfer*p4pI)/amw**2 + 
     -  (2*g1L1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3**2*m4*p2p3*
     -     pIpfer)/amw**2 - 
     -  (2*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2p4*pIpfer)/
     -   amw**2 + (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*
     -     p3p4*pIpfer)/amw**2)
       return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c---- Interfernce between diagrams with W-Prop chargino Prop and top Prop

       double precision function ampintWSMfersneut(g1L1, g1L2, 
     . g1r2, g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, mfer, amcha, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer)
      implicit none
      double precision g1L1, g1L2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer, g3l2, g3r1, g3r2, mfer
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2


      ampintWSMfersneut= 
     - 4*g1L1*g1R2*g2L1*g2R2*g3L1*g3R2*m3*m4*mI1*mI2*p1p2 - 
     -  (2*g1L1*g1R2*g2L1*g2R2*g3L1*g3R2*m3**3*m4*mI1*mI2*p1p2)/
     -   amw**2 - 4*g1L1*g1L2*g2L1*g2L2*g3L1*g3R2*m1*m4*mI1*mI2*
     -   p2p3 - (2*g1L1*g1L2*g2L1*g2L2*g3L1*g3R2*m1*m3**2*m4*
     -     mI1*mI2*p2p3)/amw**2 + 
     -  (4*g1L1*g1R2*g2L1*g2R2*g3L1*g3R2*m3*m4*mI1*mI2*p1p3*
     -     p2p3)/amw**2 + 
     -  4*g1L1*g1R2*g2L1*g2L2*g3L1*g3R2*m4*mI1*p1pI*p2p3 + 
     -  (2*g1L1*g1R2*g2L1*g2L2*g3L1*g3R2*m3**2*m4*mI1*p1pI*
     -     p2p3)/amw**2 - 
     -  4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI2*p1pfer*p2p4 + 
     -  (2*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*mI2*p1pfer*p2p4)/
     -   amw**2 + (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI2*
     -     p1pfer*p2p3*p2p4)/amw**2 + 
     -  4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI2*p1p4*p2pfer - 
     -  (2*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*mI2*p1p4*p2pfer)/
     -   amw**2 - (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI2*p1p4*
     -     p2p3*p2pfer)/amw**2 - 
     -  4*g1L1*g1L2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mI1*p2pI + 
     -  (2*g1L1*g1L2*g2L1*g2R2*g3L1*g3R2*m1*m3**3*m4*mI1*p2pI)/
     -   amw**2 - 4*g1L1*g1R2*g2L1*g2L2*g3L1*g3R2*m4*mI1*p1p3*
     -   p2pI + (2*g1L1*g1R2*g2L1*g2L2*g3L1*g3R2*m3**2*m4*mI1*
     -     p1p3*p2pI)/amw**2 + 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1pfer*p2p4*
     -     p2pI)/amw**2 - 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p4*p2pfer*
     -     p2pI)/amw**2 - 
     -  4*g1L1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*mI2*p2pfer*p3p4 - 
     -  (2*g1L1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*m3**2*mI2*p2pfer*
     -     p3p4)/amw**2 + 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI2*p1p2*p2pfer*
     -     p3p4)/amw**2 + 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI2*p1p3*p2pfer*
     -     p3p4)/amw**2 + 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1pI*p2pfer*p3p4 + 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1pI*p2pfer*
     -     p3p4)/amw**2 - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1pfer*p2pI*p3p4 + 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1pfer*p2pI*
     -     p3p4)/amw**2 - 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2pfer*p2pI*
     -     p3p4)/amw**2 + 
     -  4*g1L1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*mI2*p2p4*p3pfer + 
     -  (2*g1L1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*m3**2*mI2*p2p4*
     -     p3pfer)/amw**2 - 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI2*p1p2*p2p4*
     -     p3pfer)/amw**2 - 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI2*p1p3*p2p4*
     -     p3pfer)/amw**2 - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1pI*p2p4*p3pfer - 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1pI*p2p4*
     -     p3pfer)/amw**2 + 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p4*p2pI*p3pfer - 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p4*p2pI*
     -     p3pfer)/amw**2 + 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*p2pI*
     -     p3pfer)/amw**2 + 
     -  4*g1L1*g1R2*g2L1*g2L2*g3L1*g3R2*m4*mI1*p1p2*p3pI - 
     -  (2*g1L1*g1R2*g2L1*g2L2*g3L1*g3R2*m3**2*m4*mI1*p1p2*
     -     p3pI)/amw**2 - 
     -  (4*g1L1*g1L2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mI1*p2p3*
     -     p3pI)/amw**2 - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1pfer*p2p4*p3pI + 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1pfer*p2p4*
     -     p3pI)/amw**2 + 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p4*p2pfer*p3pI - 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p4*p2pfer*
     -     p3pI)/amw**2 - 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2pfer*p3p4*
     -     p3pI)/amw**2 + 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*p3pfer*
     -     p3pI)/amw**2 - 
     -  4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI2*p1p2*p4pfer + 
     -  (2*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*mI2*p1p2*p4pfer)/
     -   amw**2 + 4*g1L1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*mI2*p2p3*
     -   p4pfer + (2*g1L1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*m3**2*mI2*
     -     p2p3*p4pfer)/amw**2 - 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI2*p1p3*p2p3*
     -     p4pfer)/amw**2 - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1pI*p2p3*p4pfer - 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1pI*p2p3*
     -     p4pfer)/amw**2 + 
     -  4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2pI*p4pfer - 
     -  (2*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2pI*p4pfer)/
     -   amw**2 + 4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p3*p2pI*
     -   p4pfer - (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p3*
     -     p2pI*p4pfer)/amw**2 - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p2*p3pI*p4pfer + 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*p3pI*
     -     p4pfer)/amw**2 + 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*p3pI*
     -     p4pfer)/amw**2 + 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1pfer*p2p3*p4pI - 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1pfer*p2p3*
     -     p4pI)/amw**2 - 
     -  4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2pfer*p4pI + 
     -  (2*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2pfer*p4pI)/
     -   amw**2 + (4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*
     -     p2pfer*p4pI)/amw**2 - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p3*p2pfer*p4pI + 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p3*p2pfer*
     -     p4pI)/amw**2 + 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*p2pfer*
     -     p4pI)/amw**2 - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p2*p3pfer*p4pI + 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*p3pfer*
     -     p4pI)/amw**2 - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p4*p2p3*pIpfer + 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p4*p2p3*
     -     pIpfer)/amw**2 + 
     -  4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*pIpfer - 
     -  (2*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2p4*pIpfer)/
     -   amw**2 - (4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*
     -     p2p4*pIpfer)/amw**2 + 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p3*p2p4*pIpfer - 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p3*p2p4*
     -     pIpfer)/amw**2 - 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*p2p4*
     -     pIpfer)/amw**2 + 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p2*p3p4*pIpfer - 
     -  (2*g1L1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*p3p4*
     -     pIpfer)/amw**2
      
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c---- Diagram squared with W boson squark propagators
c---- g1 kopplung squark squark W
c---  g2l1, g2r1 neutralino squark quark Kopplung
c---- g2l2, g2r2 neutralino squark quark Kopplung
c---- g3  Vector boson Kopplung fermionen
c---- MV= mass of vector boson
       double precision function ampvecscalar(g11, g12, 
     . g2l1, g2l2, g2r1, g2r2,g3l,g3r,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, mv)
      implicit none
      double precision g11, g12, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer, g3l, g3r, g3r2, mfer, mv




      ampvecscalar=
     - g11*g12* ( -16*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p2 - 
     -  16*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p2 - 
     -  16*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p2 - 
     -  16*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p2 + 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**4*m4*p1p2)/MV**2 + 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**4*m4*p1p2)/MV**2 + 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**4*m4*p1p2)/MV**2 + 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**4*m4*p1p2)/MV**2 + 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p2**2)/MV**2 + 
     -  (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p2**2)/MV**2 + 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p2**2)/MV**2 + 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p2**2)/MV**2 - 
     -  32*g2L2*g2R1*g3L**2*m1*m4*p1p2*p1p3 - 
     -  32*g2L1*g2R2*g3L**2*m1*m4*p1p2*p1p3 - 
     -  32*g2L2*g2R1*g3R**2*m1*m4*p1p2*p1p3 - 
     -  32*g2L1*g2R2*g3R**2*m1*m4*p1p2*p1p3 + 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p1p3)/MV**2 + 
     -  (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p2*p1p3)/MV**2 + 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p2*p1p3)/MV**2 + 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p2*p1p3)/MV**2 + 
     -  16*g2L1*g2L2*g3L**2*m3**2*p1p2*p1p4 + 
     -  16*g2R1*g2R2*g3L**2*m3**2*p1p2*p1p4 + 
     -  16*g2L1*g2L2*g3R**2*m3**2*p1p2*p1p4 + 
     -  16*g2R1*g2R2*g3R**2*m3**2*p1p2*p1p4 - 
     -  (16*g2L1*g2L2*g3L**2*m3**4*p1p2*p1p4)/MV**2 - 
     -  (16*g2R1*g2R2*g3L**2*m3**4*p1p2*p1p4)/MV**2 - 
     -  (16*g2L1*g2L2*g3R**2*m3**4*p1p2*p1p4)/MV**2 - 
     -  (16*g2R1*g2R2*g3R**2*m3**4*p1p2*p1p4)/MV**2 - 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p2**2*p1p4)/MV**2 - 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p2**2*p1p4)/MV**2 - 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p2**2*p1p4)/MV**2 - 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p2**2*p1p4)/MV**2 + 
     -  32*g2L1*g2L2*g3L**2*p1p2*p1p3*p1p4 + 
     -  32*g2R1*g2R2*g3L**2*p1p2*p1p3*p1p4 + 
     -  32*g2L1*g2L2*g3R**2*p1p2*p1p3*p1p4 + 
     -  32*g2R1*g2R2*g3R**2*p1p2*p1p3*p1p4 - 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p2*p1p3*p1p4)/MV**2 - 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p2*p1p3*p1p4)/MV**2 - 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p2*p1p3*p1p4)/MV**2 - 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p2*p1p3*p1p4)/MV**2 + 
     -  16*g2L2*g2R1*g3L**2*m1**3*m4*p2p3 + 
     -  16*g2L1*g2R2*g3L**2*m1**3*m4*p2p3 + 
     -  16*g2L2*g2R1*g3R**2*m1**3*m4*p2p3 + 
     -  16*g2L1*g2R2*g3R**2*m1**3*m4*p2p3 - 
     -  4*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p3 - 
     -  4*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p3 - 
     -  4*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p3 - 
     -  4*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p3 + 
     -  16*g2L2*g2R1*g3L**2*m1*m4**3*p2p3 + 
     -  16*g2L1*g2R2*g3L**2*m1*m4**3*p2p3 + 
     -  16*g2L2*g2R1*g3R**2*m1*m4**3*p2p3 + 
     -  16*g2L1*g2R2*g3R**2*m1*m4**3*p2p3 - 
     -  (4*g2L2*g2R1*g3L**2*m1*m3**6*m4*p2p3)/MV**4 - 
     -  (4*g2L1*g2R2*g3L**2*m1*m3**6*m4*p2p3)/MV**4 - 
     -  (4*g2L2*g2R1*g3R**2*m1*m3**6*m4*p2p3)/MV**4 - 
     -  (4*g2L1*g2R2*g3R**2*m1*m3**6*m4*p2p3)/MV**4 + 
     -  (8*g2L2*g2R1*g3L**2*m1*m3**4*m4*p2p3)/MV**2 + 
     -  (8*g2L1*g2R2*g3L**2*m1*m3**4*m4*p2p3)/MV**2 + 
     -  (8*g2L2*g2R1*g3R**2*m1*m3**4*m4*p2p3)/MV**2 + 
     -  (8*g2L1*g2R2*g3R**2*m1*m3**4*m4*p2p3)/MV**2 - 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**4*m4*p1p2*p2p3)/MV**4 - 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**4*m4*p1p2*p2p3)/MV**4 - 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**4*m4*p1p2*p2p3)/MV**4 - 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**4*m4*p1p2*p2p3)/MV**4 + 
     -  (48*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p2p3)/MV**2 + 
     -  (48*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p2*p2p3)/MV**2 + 
     -  (48*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p2*p2p3)/MV**2 + 
     -  (48*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p2*p2p3)/MV**2 - 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p2**2*p2p3)/MV**4 - 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p2**2*p2p3)/MV**4 - 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p2**2*p2p3)/MV**4 - 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p2**2*p2p3)/MV**4 - 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**4*m4*p1p3*p2p3)/MV**4 - 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**4*m4*p1p3*p2p3)/MV**4 - 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**4*m4*p1p3*p2p3)/MV**4 - 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**4*m4*p1p3*p2p3)/MV**4 + 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p3*p2p3)/MV**2 + 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p3*p2p3)/MV**2 + 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p3*p2p3)/MV**2 + 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p3*p2p3)/MV**2 - 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p1p3*p2p3)/
     -   MV**4 - (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p2*p1p3*
     -     p2p3)/MV**4 - (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p2*
     -     p1p3*p2p3)/MV**4 - 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p2*p1p3*p2p3)/
     -   MV**4 - (16*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p3**2*p2p3)/
     -   MV**4 - (16*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p3**2*p2p3)/
     -   MV**4 - (16*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p3**2*p2p3)/
     -   MV**4 - (16*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p3**2*p2p3)/
     -   MV**4 - 16*g2L1*g2L2*g3L**2*m1**2*p1p4*p2p3 - 
     -  16*g2R1*g2R2*g3L**2*m1**2*p1p4*p2p3 - 
     -  16*g2L1*g2L2*g3R**2*m1**2*p1p4*p2p3 - 
     -  16*g2R1*g2R2*g3R**2*m1**2*p1p4*p2p3 + 
     -  4*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p3 + 
     -  4*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3 + 
     -  4*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p3 + 
     -  4*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p3 + 
     -  32*g2L2*g2R1*g3L**2*m1*m4*p1p4*p2p3 + 
     -  32*g2L1*g2R2*g3L**2*m1*m4*p1p4*p2p3 + 
     -  32*g2L2*g2R1*g3R**2*m1*m4*p1p4*p2p3 + 
     -  32*g2L1*g2R2*g3R**2*m1*m4*p1p4*p2p3 - 
     -  16*g2L1*g2L2*g3L**2*m4**2*p1p4*p2p3 - 
     -  16*g2R1*g2R2*g3L**2*m4**2*p1p4*p2p3 - 
     -  16*g2L1*g2L2*g3R**2*m4**2*p1p4*p2p3 - 
     -  16*g2R1*g2R2*g3R**2*m4**2*p1p4*p2p3 + 
     -  (4*g2L1*g2L2*g3L**2*m3**6*p1p4*p2p3)/MV**4 + 
     -  (4*g2R1*g2R2*g3L**2*m3**6*p1p4*p2p3)/MV**4 + 
     -  (4*g2L1*g2L2*g3R**2*m3**6*p1p4*p2p3)/MV**4 + 
     -  (4*g2R1*g2R2*g3R**2*m3**6*p1p4*p2p3)/MV**4 - 
     -  (8*g2L1*g2L2*g3L**2*m3**4*p1p4*p2p3)/MV**2 - 
     -  (8*g2R1*g2R2*g3L**2*m3**4*p1p4*p2p3)/MV**2 - 
     -  (8*g2L1*g2L2*g3R**2*m3**4*p1p4*p2p3)/MV**2 - 
     -  (8*g2R1*g2R2*g3R**2*m3**4*p1p4*p2p3)/MV**2 + 
     -  (16*g2L1*g2L2*g3L**2*m3**4*p1p2*p1p4*p2p3)/MV**4 + 
     -  (16*g2R1*g2R2*g3L**2*m3**4*p1p2*p1p4*p2p3)/MV**4 + 
     -  (16*g2L1*g2L2*g3R**2*m3**4*p1p2*p1p4*p2p3)/MV**4 + 
     -  (16*g2R1*g2R2*g3R**2*m3**4*p1p2*p1p4*p2p3)/MV**4 - 
     -  (48*g2L1*g2L2*g3L**2*m3**2*p1p2*p1p4*p2p3)/MV**2 - 
     -  (48*g2R1*g2R2*g3L**2*m3**2*p1p2*p1p4*p2p3)/MV**2 - 
     -  (48*g2L1*g2L2*g3R**2*m3**2*p1p2*p1p4*p2p3)/MV**2 - 
     -  (48*g2R1*g2R2*g3R**2*m3**2*p1p2*p1p4*p2p3)/MV**2 + 
     -  (16*g2L1*g2L2*g3L**2*m3**2*p1p2**2*p1p4*p2p3)/MV**4 + 
     -  (16*g2R1*g2R2*g3L**2*m3**2*p1p2**2*p1p4*p2p3)/MV**4 + 
     -  (16*g2L1*g2L2*g3R**2*m3**2*p1p2**2*p1p4*p2p3)/MV**4 + 
     -  (16*g2R1*g2R2*g3R**2*m3**2*p1p2**2*p1p4*p2p3)/MV**4 + 
     -  (16*g2L1*g2L2*g3L**2*m3**4*p1p3*p1p4*p2p3)/MV**4 + 
     -  (16*g2R1*g2R2*g3L**2*m3**4*p1p3*p1p4*p2p3)/MV**4 + 
     -  (16*g2L1*g2L2*g3R**2*m3**4*p1p3*p1p4*p2p3)/MV**4 + 
     -  (16*g2R1*g2R2*g3R**2*m3**4*p1p3*p1p4*p2p3)/MV**4 - 
     -  (16*g2L1*g2L2*g3L**2*m3**2*p1p3*p1p4*p2p3)/MV**2 - 
     -  (16*g2R1*g2R2*g3L**2*m3**2*p1p3*p1p4*p2p3)/MV**2 - 
     -  (16*g2L1*g2L2*g3R**2*m3**2*p1p3*p1p4*p2p3)/MV**2 - 
     -  (16*g2R1*g2R2*g3R**2*m3**2*p1p3*p1p4*p2p3)/MV**2 + 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p2*p1p3*p1p4*p2p3)/MV**4 + 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p2*p1p3*p1p4*p2p3)/MV**4 + 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p2*p1p3*p1p4*p2p3)/MV**4 + 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p2*p1p3*p1p4*p2p3)/MV**4 + 
     -  (16*g2L1*g2L2*g3L**2*m3**2*p1p3**2*p1p4*p2p3)/MV**4 + 
     -  (16*g2R1*g2R2*g3L**2*m3**2*p1p3**2*p1p4*p2p3)/MV**4 + 
     -  (16*g2L1*g2L2*g3R**2*m3**2*p1p3**2*p1p4*p2p3)/MV**4 + 
     -  (16*g2R1*g2R2*g3R**2*m3**2*p1p3**2*p1p4*p2p3)/MV**4 - 
     -  32*g2L1*g2L2*g3L**2*p1p4**2*p2p3 - 
     -  32*g2R1*g2R2*g3L**2*p1p4**2*p2p3 - 
     -  32*g2L1*g2L2*g3R**2*p1p4**2*p2p3 - 
     -  32*g2R1*g2R2*g3R**2*p1p4**2*p2p3 - 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**4*m4*p2p3**2)/MV**4 - 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**4*m4*p2p3**2)/MV**4 - 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**4*m4*p2p3**2)/MV**4 - 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**4*m4*p2p3**2)/MV**4 + 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p3**2)/MV**2 + 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p3**2)/MV**2 + 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p3**2)/MV**2 + 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p3**2)/MV**2 - 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p2p3**2)/MV**4 - 
     -  (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p2*p2p3**2)/MV**4 - 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p2*p2p3**2)/MV**4 - 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p2*p2p3**2)/MV**4 - 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p3*p2p3**2)/MV**4 - 
     -  (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p3*p2p3**2)/MV**4 - 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p3*p2p3**2)/MV**4 - 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p3*p2p3**2)/MV**4 + 
     -  (16*g2L1*g2L2*g3L**2*m3**4*p1p4*p2p3**2)/MV**4 + 
     -  (16*g2R1*g2R2*g3L**2*m3**4*p1p4*p2p3**2)/MV**4 + 
     -  (16*g2L1*g2L2*g3R**2*m3**4*p1p4*p2p3**2)/MV**4 + 
     -  (16*g2R1*g2R2*g3R**2*m3**4*p1p4*p2p3**2)/MV**4 - 
     -  (16*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p3**2)/MV**2 - 
     -  (16*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3**2)/MV**2 - 
     -  (16*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p3**2)/MV**2 - 
     -  (16*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p3**2)/MV**2 + 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p2*p1p4*p2p3**2)/MV**4 + 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p2*p1p4*p2p3**2)/MV**4 + 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p2*p1p4*p2p3**2)/MV**4 + 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p2*p1p4*p2p3**2)/MV**4 + 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p3*p1p4*p2p3**2)/MV**4 + 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p3*p1p4*p2p3**2)/MV**4 + 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p3*p1p4*p2p3**2)/MV**4 + 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p3*p1p4*p2p3**2)/MV**4 - 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p3**3)/MV**4 - 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p3**3)/MV**4 - 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p3**3)/MV**4 - 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p3**3)/MV**4 + 
     -  (16*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p3**3)/MV**4 + 
     -  (16*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3**3)/MV**4 + 
     -  (16*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p3**3)/MV**4 + 
     -  (16*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p3**3)/MV**4 - 
     -  16*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p4 - 
     -  16*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p4 - 
     -  16*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p4 - 
     -  16*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p4 + 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**4*m4*p2p4)/MV**2 + 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**4*m4*p2p4)/MV**2 + 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**4*m4*p2p4)/MV**2 + 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**4*m4*p2p4)/MV**2 + 
     -  (64*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p2p4)/MV**2 + 
     -  (64*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p2*p2p4)/MV**2 + 
     -  (64*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p2*p2p4)/MV**2 + 
     -  (64*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p2*p2p4)/MV**2 - 
     -  32*g2L2*g2R1*g3L**2*m1*m4*p1p3*p2p4 - 
     -  32*g2L1*g2R2*g3L**2*m1*m4*p1p3*p2p4 - 
     -  32*g2L2*g2R1*g3R**2*m1*m4*p1p3*p2p4 - 
     -  32*g2L1*g2R2*g3R**2*m1*m4*p1p3*p2p4 + 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p3*p2p4)/MV**2 + 
     -  (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p3*p2p4)/MV**2 + 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p3*p2p4)/MV**2 + 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p3*p2p4)/MV**2 + 
     -  16*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p4 + 
     -  16*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p4 + 
     -  16*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p4 + 
     -  16*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p4 - 
     -  (16*g2L1*g2L2*g3L**2*m3**4*p1p4*p2p4)/MV**2 - 
     -  (16*g2R1*g2R2*g3L**2*m3**4*p1p4*p2p4)/MV**2 - 
     -  (16*g2L1*g2L2*g3R**2*m3**4*p1p4*p2p4)/MV**2 - 
     -  (16*g2R1*g2R2*g3R**2*m3**4*p1p4*p2p4)/MV**2 - 
     -  (64*g2L1*g2L2*g3L**2*m3**2*p1p2*p1p4*p2p4)/MV**2 - 
     -  (64*g2R1*g2R2*g3L**2*m3**2*p1p2*p1p4*p2p4)/MV**2 - 
     -  (64*g2L1*g2L2*g3R**2*m3**2*p1p2*p1p4*p2p4)/MV**2 - 
     -  (64*g2R1*g2R2*g3R**2*m3**2*p1p2*p1p4*p2p4)/MV**2 + 
     -  32*g2L1*g2L2*g3L**2*p1p3*p1p4*p2p4 + 
     -  32*g2R1*g2R2*g3L**2*p1p3*p1p4*p2p4 + 
     -  32*g2L1*g2L2*g3R**2*p1p3*p1p4*p2p4 + 
     -  32*g2R1*g2R2*g3R**2*p1p3*p1p4*p2p4 - 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p3*p1p4*p2p4)/MV**2 - 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p3*p1p4*p2p4)/MV**2 - 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p3*p1p4*p2p4)/MV**2 - 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p3*p1p4*p2p4)/MV**2 - 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**4*m4*p2p3*p2p4)/MV**4 - 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**4*m4*p2p3*p2p4)/MV**4 - 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**4*m4*p2p3*p2p4)/MV**4 - 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**4*m4*p2p3*p2p4)/MV**4 + 
     -  (48*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p3*p2p4)/MV**2 + 
     -  (48*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p3*p2p4)/MV**2 + 
     -  (48*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p3*p2p4)/MV**2 + 
     -  (48*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p3*p2p4)/MV**2 - 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p2p3*p2p4)/
     -   MV**4 - (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p2*p2p3*
     -     p2p4)/MV**4 - (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p2*
     -     p2p3*p2p4)/MV**4 - 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p2*p2p3*p2p4)/
     -   MV**4 - (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p3*p2p3*
     -     p2p4)/MV**4 - (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p3*
     -     p2p3*p2p4)/MV**4 - 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p3*p2p3*p2p4)/
     -   MV**4 - (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p3*p2p3*
     -     p2p4)/MV**4 + (16*g2L1*g2L2*g3L**2*m3**4*p1p4*p2p3*
     -     p2p4)/MV**4 + (16*g2R1*g2R2*g3L**2*m3**4*p1p4*p2p3*
     -     p2p4)/MV**4 + (16*g2L1*g2L2*g3R**2*m3**4*p1p4*p2p3*
     -     p2p4)/MV**4 + (16*g2R1*g2R2*g3R**2*m3**4*p1p4*p2p3*
     -     p2p4)/MV**4 - (48*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p3*
     -     p2p4)/MV**2 - (48*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3*
     -     p2p4)/MV**2 - (48*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p3*
     -     p2p4)/MV**2 - (48*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p3*
     -     p2p4)/MV**2 + (32*g2L1*g2L2*g3L**2*m3**2*p1p2*p1p4*
     -     p2p3*p2p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p2*p1p4*p2p3*p2p4)/MV**4 + 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p2*p1p4*p2p3*p2p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p2*p1p4*p2p3*p2p4)/MV**4 + 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p3*p1p4*p2p3*p2p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p3*p1p4*p2p3*p2p4)/MV**4 + 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p3*p1p4*p2p3*p2p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p3*p1p4*p2p3*p2p4)/MV**4 - 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p3**2*p2p4)/MV**4 - 
     -  (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p3**2*p2p4)/MV**4 - 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p3**2*p2p4)/MV**4 - 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p3**2*p2p4)/MV**4 + 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p3**2*p2p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3**2*p2p4)/MV**4 + 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p3**2*p2p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p3**2*p2p4)/MV**4 + 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p4**2)/MV**2 + 
     -  (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p4**2)/MV**2 + 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p4**2)/MV**2 + 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p4**2)/MV**2 - 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p4**2)/MV**2 - 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p4**2)/MV**2 - 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p4**2)/MV**2 - 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p4**2)/MV**2 - 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p3*p2p4**2)/MV**4 - 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p3*p2p4**2)/MV**4 - 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p3*p2p4**2)/MV**4 - 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p3*p2p4**2)/MV**4 + 
     -  (16*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p3*p2p4**2)/MV**4 + 
     -  (16*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3*p2p4**2)/MV**4 + 
     -  (16*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p3*p2p4**2)/MV**4 + 
     -  (16*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p3*p2p4**2)/MV**4 - 
     -  32*g2L2*g2R1*g3L**2*m1*m4*p1p2*p3p4 - 
     -  32*g2L1*g2R2*g3L**2*m1*m4*p1p2*p3p4 - 
     -  32*g2L2*g2R1*g3R**2*m1*m4*p1p2*p3p4 - 
     -  32*g2L1*g2R2*g3R**2*m1*m4*p1p2*p3p4 + 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p3p4)/MV**2 + 
     -  (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p2*p3p4)/MV**2 + 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p2*p3p4)/MV**2 + 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p2*p3p4)/MV**2 + 
     -  32*g2L1*g2L2*g3L**2*p1p2*p1p4*p3p4 + 
     -  32*g2R1*g2R2*g3L**2*p1p2*p1p4*p3p4 + 
     -  32*g2L1*g2L2*g3R**2*p1p2*p1p4*p3p4 + 
     -  32*g2R1*g2R2*g3R**2*p1p2*p1p4*p3p4 - 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p2*p1p4*p3p4)/MV**2 - 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p2*p1p4*p3p4)/MV**2 - 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p2*p1p4*p3p4)/MV**2 - 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p2*p1p4*p3p4)/MV**2 - 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**4*m4*p2p3*p3p4)/MV**4 - 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**4*m4*p2p3*p3p4)/MV**4 - 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**4*m4*p2p3*p3p4)/MV**4 - 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**4*m4*p2p3*p3p4)/MV**4 + 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p3*p3p4)/MV**2 + 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p3*p3p4)/MV**2 + 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p3*p3p4)/MV**2 + 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p3*p3p4)/MV**2 - 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p2p3*p3p4)/
     -   MV**4 - (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p2*p2p3*
     -     p3p4)/MV**4 - (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p2*
     -     p2p3*p3p4)/MV**4 - 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p2*p2p3*p3p4)/
     -   MV**4 - (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p1p3*p2p3*
     -     p3p4)/MV**4 - (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p1p3*
     -     p2p3*p3p4)/MV**4 - 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p1p3*p2p3*p3p4)/
     -   MV**4 - (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p1p3*p2p3*
     -     p3p4)/MV**4 + (16*g2L1*g2L2*g3L**2*m3**4*p1p4*p2p3*
     -     p3p4)/MV**4 + (16*g2R1*g2R2*g3L**2*m3**4*p1p4*p2p3*
     -     p3p4)/MV**4 + (16*g2L1*g2L2*g3R**2*m3**4*p1p4*p2p3*
     -     p3p4)/MV**4 + (16*g2R1*g2R2*g3R**2*m3**4*p1p4*p2p3*
     -     p3p4)/MV**4 - (16*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p3*
     -     p3p4)/MV**2 - (16*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3*
     -     p3p4)/MV**2 - (16*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p3*
     -     p3p4)/MV**2 - (16*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p3*
     -     p3p4)/MV**2 + (32*g2L1*g2L2*g3L**2*m3**2*p1p2*p1p4*
     -     p2p3*p3p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p2*p1p4*p2p3*p3p4)/MV**4 + 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p2*p1p4*p2p3*p3p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p2*p1p4*p2p3*p3p4)/MV**4 + 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p3*p1p4*p2p3*p3p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p3*p1p4*p2p3*p3p4)/MV**4 + 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p3*p1p4*p2p3*p3p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p3*p1p4*p2p3*p3p4)/MV**4 - 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p3**2*p3p4)/MV**4 - 
     -  (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p3**2*p3p4)/MV**4 - 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p3**2*p3p4)/MV**4 - 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p3**2*p3p4)/MV**4 + 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p3**2*p3p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3**2*p3p4)/MV**4 + 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p3**2*p3p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p3**2*p3p4)/MV**4 - 
     -  32*g2L2*g2R1*g3L**2*m1*m4*p2p4*p3p4 - 
     -  32*g2L1*g2R2*g3L**2*m1*m4*p2p4*p3p4 - 
     -  32*g2L2*g2R1*g3R**2*m1*m4*p2p4*p3p4 - 
     -  32*g2L1*g2R2*g3R**2*m1*m4*p2p4*p3p4 + 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p4*p3p4)/MV**2 + 
     -  (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p4*p3p4)/MV**2 + 
     -  (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p4*p3p4)/MV**2 + 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p4*p3p4)/MV**2 + 
     -  32*g2L1*g2L2*g3L**2*p1p4*p2p4*p3p4 + 
     -  32*g2R1*g2R2*g3L**2*p1p4*p2p4*p3p4 + 
     -  32*g2L1*g2L2*g3R**2*p1p4*p2p4*p3p4 + 
     -  32*g2R1*g2R2*g3R**2*p1p4*p2p4*p3p4 - 
     -  (32*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p4*p3p4)/MV**2 - 
     -  (32*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p4*p3p4)/MV**2 - 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p4*p3p4)/MV**2 - 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p4*p3p4)/MV**2 - 
     -  (32*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p3*p2p4*p3p4)/
     -   MV**4 - (32*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p3*p2p4*
     -     p3p4)/MV**4 - (32*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p3*
     -     p2p4*p3p4)/MV**4 - 
     -  (32*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p3*p2p4*p3p4)/
     -   MV**4 + (32*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p3*p2p4*
     -     p3p4)/MV**4 + (32*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3*
     -     p2p4*p3p4)/MV**4 + 
     -  (32*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p3*p2p4*p3p4)/MV**4 + 
     -  (32*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p3*p2p4*p3p4)/MV**4 - 
     -  (16*g2L2*g2R1*g3L**2*m1*m3**2*m4*p2p3*p3p4**2)/MV**4 - 
     -  (16*g2L1*g2R2*g3L**2*m1*m3**2*m4*p2p3*p3p4**2)/MV**4 - 
     -  (16*g2L2*g2R1*g3R**2*m1*m3**2*m4*p2p3*p3p4**2)/MV**4 - 
     -  (16*g2L1*g2R2*g3R**2*m1*m3**2*m4*p2p3*p3p4**2)/MV**4 + 
     -  (16*g2L1*g2L2*g3L**2*m3**2*p1p4*p2p3*p3p4**2)/MV**4 + 
     -  (16*g2R1*g2R2*g3L**2*m3**2*p1p4*p2p3*p3p4**2)/MV**4 + 
     -  (16*g2L1*g2L2*g3R**2*m3**2*p1p4*p2p3*p3p4**2)/MV**4 + 
     -  (16*g2R1*g2R2*g3R**2*m3**2*p1p4*p2p3*p3p4**2)/MV**4)
      
        
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c---- Diagram squared with W boson squark propagators
c---- g1 kopplung squark squark W
c---  g2l1, g2r1 neutralino squark quark Kopplung
c---- g2l2, g2r2 neutralino squark quark Kopplung
c---- g3  Vector boson Kopplung fermionen
c---- MV= mass of vector boson
       double precision function amphiggssquark(g11, g12, 
     . g2l1, g2l2, g2r1, g2r2,g3l,g3r,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4)
      implicit none
      double precision g11, g12, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha, pfer2, p1pfer, p2pfer,
     .  p3pfer, p4pfer, pIpfer, g3l, g3r, g3r2, mfer, mv
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      amphiggssquark=g11*g12*(-4*g2L2*g2R1*g3L**2*m1*m4*p2p3 - 
     -  4*g2L1*g2R2*g3L**2*m1*m4*p2p3 - 
     -  4*g2L2*g2R1*g3R**2*m1*m4*p2p3 - 
     -  4*g2L1*g2R2*g3R**2*m1*m4*p2p3 + 
     -  4*g2L1*g2L2*g3L**2*p1p4*p2p3 + 
     -  4*g2R1*g2R2*g3L**2*p1p4*p2p3 + 
     -  4*g2L1*g2L2*g3R**2*p1p4*p2p3 + 
     -  4*g2R1*g2R2*g3R**2*p1p4*p2p3)
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c---- Interference term Higgs squark- W squark diagram
c---- g11 kopplung squark squark W
c---- g12 Kopplung squark squark charged Higgs
c---  g2l1, g2r1 neutralino squark quark Kopplung
c---- g2l2, g2r2 neutralino squark quark Kopplung
c---- g3l1  Vector boson Kopplung fermionen
c--- g3l2 g3r2 charged higgs fermion kopplung

       double precision function ampinthiggssquarkWsquark(g11, g12, 
     . g2l1, g2l2, g2r1, g2r2,g3l1,g3l2,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4)
      implicit none
      double precision g11, g12, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha,
     .g3l2, g3r1, g3r2, mfer, mv
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2
      
      
      ampinthiggssquarkWsquark=
     -   -4*g11*g12*(-(g2L2*g2R1*m1*m4) - g2L1*g2R2*m1*m4 + 
     -    g2L1*g2L2*p1p4 + g2R1*g2R2*p1p4)*
     -  (-(g3L1*g3R2*m3*(2*p1p2 + p2p3 + 2*p2p4)) + 
     -    (g3L1*g3R2*m3*p2p3*
     -       (m3**2 + 2*(p1p2 + p1p3 + p2p3 + p2p4 + p3p4)))/
     -     amw**2)
       
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c---- Interference term W chargino- W squark diagram
c---- g1 kopplung squark squark W
c---- g1l2, g1r2 Kopplung squark chargino quark
c---  g2l1, g2r1 neutralino squark quark Kopplung
c---- g2l2, g2r2 neutralino W chargino Kopplung
c---- g3l1  Vector boson Kopplung fermionen
c--- g3l2 g3r2 charged higgs fermion kopplung

       double precision function ampintWcharWsquark(g1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amcha)
      implicit none
      double precision g1, g1l2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha,
     .g3l2, g3r1, g3r2, mfer, mv, g3l
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintWcharWsquark=
     - -g1*(-4*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p2 - 
     -  4*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p2 + 
     -  (4*g1R2*g2R2*g2R1*g3L**2*m3**4*m4*amcha*p1p2)/amw**2 + 
     -  (4*g1L2*g2L1*g2L2*g3L**2*m3**4*m4*amcha*p1p2)/amw**2 + 
     -  (16*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p2**2)/amw**2 + 
     -  (16*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p2**2)/amw**2 - 
     -  16*g1R2*g2R2*g2R1*g3L**2*m4*amcha*p1p2*p1p3 - 
     -  16*g1L2*g2L1*g2L2*g3L**2*m4*amcha*p1p2*p1p3 + 
     -  (16*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p2*p1p3)/
     -   amw**2 + (16*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p2*
     -     p1p3)/amw**2 + 
     -  8*g1R2*g2R2*g2R1*g3L**2*m1**2*m4*amcha*p2p3 + 
     -  8*g1L2*g2L1*g2L2*g3L**2*m1**2*m4*amcha*p2p3 - 
     -  8*g1R2*g2L1*g2R2*g3L**2*m1*m4**2*amcha*p2p3 - 
     -  8*g1L2*g2R1*g2L2*g3L**2*m1*m4**2*amcha*p2p3 + 
     -  (12*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p2*p2p3)/
     -   amw**2 + (12*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p2*
     -     p2p3)/amw**2 - 
     -  (4*g1R2*g2R2*g2R1*g3L**2*m3**4*m4*amcha*p1p2*p2p3)/
     -   amw**4 - (4*g1L2*g2L1*g2L2*g3L**2*m3**4*m4*amcha*p1p2*
     -     p2p3)/amw**4 - 
     -  (8*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p2**2*p2p3)/
     -   amw**4 - (8*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p2**2*
     -     p2p3)/amw**4 + 
     -  (4*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p3*p2p3)/
     -   amw**2 + (4*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p3*
     -     p2p3)/amw**2 - 
     -  (4*g1R2*g2R2*g2R1*g3L**2*m3**4*m4*amcha*p1p3*p2p3)/
     -   amw**4 - (4*g1L2*g2L1*g2L2*g3L**2*m3**4*m4*amcha*p1p3*
     -     p2p3)/amw**4 - 
     -  (16*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p2*p1p3*p2p3)/
     -   amw**4 - (16*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p2*
     -     p1p3*p2p3)/amw**4 - 
     -  (8*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p3**2*p2p3)/
     -   amw**4 - (8*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p3**2*
     -     p2p3)/amw**4 - 
     -  8*g1R2*g2L1*g2R2*g3L**2*m1*amcha*p1p4*p2p3 - 
     -  8*g1L2*g2R1*g2L2*g3L**2*m1*amcha*p1p4*p2p3 + 
     -  8*g1R2*g2R2*g2R1*g3L**2*m4*amcha*p1p4*p2p3 + 
     -  8*g1L2*g2L1*g2L2*g3L**2*m4*amcha*p1p4*p2p3 - 
     -  8*g1L2*g2R2*g2R1*g3L**2*m1*m4*p1pI*p2p3 - 
     -  8*g1R2*g2L1*g2L2*g3L**2*m1*m4*p1pI*p2p3 + 
     -  8*g1L2*g2L1*g2R2*g3L**2*m4**2*p1pI*p2p3 + 
     -  8*g1R2*g2R1*g2L2*g3L**2*m4**2*p1pI*p2p3 + 
     -  16*g1L2*g2L1*g2R2*g3L**2*p1p4*p1pI*p2p3 + 
     -  16*g1R2*g2R1*g2L2*g3L**2*p1p4*p1pI*p2p3 - 
     -  (8*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p2*p2p3**2)/
     -   amw**4 - (8*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p2*
     -     p2p3**2)/amw**4 - 
     -  (8*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p3*p2p3**2)/
     -   amw**4 - (8*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p3*
     -     p2p3**2)/amw**4 + 
     -  4*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p2p4 + 
     -  4*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p2p4 - 
     -  (4*g1R2*g2L1*g2R2*g3L**2*m1*m3**4*amcha*p2p4)/amw**2 - 
     -  (4*g1L2*g2R1*g2L2*g3L**2*m1*m3**4*amcha*p2p4)/amw**2 - 
     -  (16*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p1p2*p2p4)/
     -   amw**2 - (16*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p1p2*
     -     p2p4)/amw**2 + 
     -  (16*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p2*p2p4)/
     -   amw**2 + (16*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p2*
     -     p2p4)/amw**2 + 
     -  8*g1R2*g2L1*g2R2*g3L**2*m1*amcha*p1p3*p2p4 + 
     -  8*g1L2*g2R1*g2L2*g3L**2*m1*amcha*p1p3*p2p4 - 
     -  (8*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p1p3*p2p4)/
     -   amw**2 - (8*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p1p3*
     -     p2p4)/amw**2 - 
     -  8*g1R2*g2R2*g2R1*g3L**2*m4*amcha*p1p3*p2p4 - 
     -  8*g1L2*g2L1*g2L2*g3L**2*m4*amcha*p1p3*p2p4 + 
     -  (8*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p3*p2p4)/
     -   amw**2 + (8*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p3*
     -     p2p4)/amw**2 - 
     -  4*g1L2*g2L1*g2R2*g3L**2*m3**2*p1pI*p2p4 - 
     -  4*g1R2*g2R1*g2L2*g3L**2*m3**2*p1pI*p2p4 + 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**4*p1pI*p2p4)/amw**2 + 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**4*p1pI*p2p4)/amw**2 + 
     -  (16*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p1pI*p2p4)/
     -   amw**2 + (16*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p1pI*
     -     p2p4)/amw**2 - 
     -  8*g1L2*g2L1*g2R2*g3L**2*p1p3*p1pI*p2p4 - 
     -  8*g1R2*g2R1*g2L2*g3L**2*p1p3*p1pI*p2p4 + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3*p1pI*p2p4)/amw**2 + 
     -  (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3*p1pI*p2p4)/amw**2 - 
     -  (12*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p2p3*p2p4)/
     -   amw**2 - (12*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p2p3*
     -     p2p4)/amw**2 + 
     -  (4*g1R2*g2L1*g2R2*g3L**2*m1*m3**4*amcha*p2p3*p2p4)/
     -   amw**4 + (4*g1L2*g2R1*g2L2*g3L**2*m1*m3**4*amcha*p2p3*
     -     p2p4)/amw**4 + 
     -  (8*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p1p2*p2p3*p2p4)/
     -   amw**4 + (8*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p1p2*
     -     p2p3*p2p4)/amw**4 - 
     -  (8*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p2*p2p3*p2p4)/
     -   amw**4 - (8*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p2*
     -     p2p3*p2p4)/amw**4 + 
     -  (8*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p1p3*p2p3*p2p4)/
     -   amw**4 + (8*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p1p3*
     -     p2p3*p2p4)/amw**4 - 
     -  (8*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p3*p2p3*p2p4)/
     -   amw**4 - (8*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p3*
     -     p2p3*p2p4)/amw**4 + 
     -  (12*g1L2*g2L1*g2R2*g3L**2*m3**2*p1pI*p2p3*p2p4)/
     -   amw**2 + (12*g1R2*g2R1*g2L2*g3L**2*m3**2*p1pI*p2p3*
     -     p2p4)/amw**2 - 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**4*p1pI*p2p3*p2p4)/amw**4 - 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**4*p1pI*p2p3*p2p4)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p1pI*p2p3*p2p4)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p1pI*p2p3*
     -     p2p4)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3*p1pI*p2p3*p2p4)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3*p1pI*p2p3*
     -     p2p4)/amw**4 + 
     -  (8*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p2p3**2*p2p4)/
     -   amw**4 + (8*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p2p3**2*
     -     p2p4)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1pI*p2p3**2*p2p4)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1pI*p2p3**2*
     -     p2p4)/amw**4 - 
     -  (16*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p2p4**2)/amw**2 - 
     -  (16*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p2p4**2)/amw**2 + 
     -  (16*g1L2*g2L1*g2R2*g3L**2*m3**2*p1pI*p2p4**2)/amw**2 + 
     -  (16*g1R2*g2R1*g2L2*g3L**2*m3**2*p1pI*p2p4**2)/amw**2 + 
     -  (8*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p2p3*p2p4**2)/
     -   amw**4 + (8*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p2p3*
     -     p2p4**2)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1pI*p2p3*p2p4**2)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1pI*p2p3*
     -     p2p4**2)/amw**4 + 
     -  4*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2pI + 
     -  4*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2pI - 
     -  (4*g1L2*g2R2*g2R1*g3L**2*m1*m3**4*m4*p2pI)/amw**2 - 
     -  (4*g1R2*g2L1*g2L2*g3L**2*m1*m3**4*m4*p2pI)/amw**2 - 
     -  (16*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p2pI)/
     -   amw**2 - (16*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p2*
     -     p2pI)/amw**2 + 
     -  8*g1L2*g2R2*g2R1*g3L**2*m1*m4*p1p3*p2pI + 
     -  8*g1R2*g2L1*g2L2*g3L**2*m1*m4*p1p3*p2pI - 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p1p3*p2pI)/
     -   amw**2 - (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p3*
     -     p2pI)/amw**2 - 
     -  4*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2pI - 
     -  4*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2pI + 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**4*p1p4*p2pI)/amw**2 + 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**4*p1p4*p2pI)/amw**2 + 
     -  (16*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p1p4*p2pI)/
     -   amw**2 + (16*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p1p4*
     -     p2pI)/amw**2 - 
     -  8*g1L2*g2L1*g2R2*g3L**2*p1p3*p1p4*p2pI - 
     -  8*g1R2*g2R1*g2L2*g3L**2*p1p3*p1p4*p2pI + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3*p1p4*p2pI)/amw**2 + 
     -  (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3*p1p4*p2pI)/amw**2 - 
     -  (12*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2p3*p2pI)/
     -   amw**2 - (12*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*g1L2*g2R2*g2R1*g3L**2*m1*m3**4*m4*p2p3*p2pI)/
     -   amw**4 + (4*g1R2*g2L1*g2L2*g3L**2*m1*m3**4*m4*p2p3*
     -     p2pI)/amw**4 + 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p2p3*p2pI)/
     -   amw**4 + (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p2*
     -     p2p3*p2pI)/amw**4 + 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p1p3*p2p3*p2pI)/
     -   amw**4 + (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p3*
     -     p2p3*p2pI)/amw**4 + 
     -  (12*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2p3*p2pI)/
     -   amw**2 + (12*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2p3*
     -     p2pI)/amw**2 - 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**4*p1p4*p2p3*p2pI)/amw**4 - 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**4*p1p4*p2p3*p2pI)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p1p4*p2p3*p2pI)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p1p4*p2p3*
     -     p2pI)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3*p1p4*p2p3*p2pI)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3*p1p4*p2p3*
     -     p2pI)/amw**4 + 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2p3**2*p2pI)/
     -   amw**4 + (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3**2*
     -     p2pI)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2p3**2*p2pI)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2p3**2*
     -     p2pI)/amw**4 - 
     -  (16*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2p4*p2pI)/
     -   amw**2 - (16*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p4*
     -     p2pI)/amw**2 + 
     -  (16*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2p4*p2pI)/
     -   amw**2 + (16*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2p4*
     -     p2pI)/amw**2 + 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2p3*p2p4*p2pI)/
     -   amw**4 + (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*
     -     p2p4*p2pI)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2p3*p2p4*p2pI)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2p3*p2p4*
     -     p2pI)/amw**4 + 
     -  8*g1R2*g2L1*g2R2*g3L**2*m1*amcha*p1p2*p3p4 + 
     -  8*g1L2*g2R1*g2L2*g3L**2*m1*amcha*p1p2*p3p4 - 
     -  (8*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p1p2*p3p4)/
     -   amw**2 - (8*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p1p2*
     -     p3p4)/amw**2 - 
     -  8*g1R2*g2R2*g2R1*g3L**2*m4*amcha*p1p2*p3p4 - 
     -  8*g1L2*g2L1*g2L2*g3L**2*m4*amcha*p1p2*p3p4 + 
     -  (8*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p2*p3p4)/
     -   amw**2 + (8*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p2*
     -     p3p4)/amw**2 - 
     -  8*g1L2*g2L1*g2R2*g3L**2*p1p2*p1pI*p3p4 - 
     -  8*g1R2*g2R1*g2L2*g3L**2*p1p2*p1pI*p3p4 + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p1pI*p3p4)/amw**2 + 
     -  (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p1pI*p3p4)/amw**2 - 
     -  (4*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p2p3*p3p4)/
     -   amw**2 - (4*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p2p3*
     -     p3p4)/amw**2 + 
     -  (4*g1R2*g2L1*g2R2*g3L**2*m1*m3**4*amcha*p2p3*p3p4)/
     -   amw**4 + (4*g1L2*g2R1*g2L2*g3L**2*m1*m3**4*amcha*p2p3*
     -     p3p4)/amw**4 + 
     -  (8*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p1p2*p2p3*p3p4)/
     -   amw**4 + (8*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p1p2*
     -     p2p3*p3p4)/amw**4 - 
     -  (8*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p2*p2p3*p3p4)/
     -   amw**4 - (8*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p2*
     -     p2p3*p3p4)/amw**4 + 
     -  (8*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p1p3*p2p3*p3p4)/
     -   amw**4 + (8*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p1p3*
     -     p2p3*p3p4)/amw**4 - 
     -  (8*g1R2*g2R2*g2R1*g3L**2*m3**2*m4*amcha*p1p3*p2p3*p3p4)/
     -   amw**4 - (8*g1L2*g2L1*g2L2*g3L**2*m3**2*m4*amcha*p1p3*
     -     p2p3*p3p4)/amw**4 + 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**2*p1pI*p2p3*p3p4)/amw**2 + 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**2*p1pI*p2p3*p3p4)/amw**2 - 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**4*p1pI*p2p3*p3p4)/amw**4 - 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**4*p1pI*p2p3*p3p4)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p1pI*p2p3*p3p4)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p1pI*p2p3*
     -     p3p4)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3*p1pI*p2p3*p3p4)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3*p1pI*p2p3*
     -     p3p4)/amw**4 + 
     -  (8*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p2p3**2*p3p4)/
     -   amw**4 + (8*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p2p3**2*
     -     p3p4)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1pI*p2p3**2*p3p4)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1pI*p2p3**2*
     -     p3p4)/amw**4 + 
     -  16*g1R2*g2L1*g2R2*g3L**2*m1*amcha*p2p4*p3p4 + 
     -  16*g1L2*g2R1*g2L2*g3L**2*m1*amcha*p2p4*p3p4 - 
     -  (16*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p2p4*p3p4)/
     -   amw**2 - (16*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p2p4*
     -     p3p4)/amw**2 - 
     -  16*g1L2*g2L1*g2R2*g3L**2*p1pI*p2p4*p3p4 - 
     -  16*g1R2*g2R1*g2L2*g3L**2*p1pI*p2p4*p3p4 + 
     -  (16*g1L2*g2L1*g2R2*g3L**2*m3**2*p1pI*p2p4*p3p4)/
     -   amw**2 + (16*g1R2*g2R1*g2L2*g3L**2*m3**2*p1pI*p2p4*
     -     p3p4)/amw**2 + 
     -  (16*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p2p3*p2p4*p3p4)/
     -   amw**4 + (16*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p2p3*
     -     p2p4*p3p4)/amw**4 - 
     -  (16*g1L2*g2L1*g2R2*g3L**2*m3**2*p1pI*p2p3*p2p4*p3p4)/
     -   amw**4 - (16*g1R2*g2R1*g2L2*g3L**2*m3**2*p1pI*p2p3*
     -     p2p4*p3p4)/amw**4 + 
     -  8*g1L2*g2R2*g2R1*g3L**2*m1*m4*p2pI*p3p4 + 
     -  8*g1R2*g2L1*g2L2*g3L**2*m1*m4*p2pI*p3p4 - 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2pI*p3p4)/
     -   amw**2 - (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2pI*
     -     p3p4)/amw**2 - 
     -  8*g1L2*g2L1*g2R2*g3L**2*p1p4*p2pI*p3p4 - 
     -  8*g1R2*g2R1*g2L2*g3L**2*p1p4*p2pI*p3p4 + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2pI*p3p4)/amw**2 + 
     -  (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2pI*p3p4)/amw**2 + 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2p3*p2pI*p3p4)/
     -   amw**4 + (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*
     -     p2pI*p3p4)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2p3*p2pI*p3p4)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2p3*p2pI*
     -     p3p4)/amw**4 + 
     -  (8*g1R2*g2L1*g2R2*g3L**2*m1*m3**2*amcha*p2p3*p3p4**2)/
     -   amw**4 + (8*g1L2*g2R1*g2L2*g3L**2*m1*m3**2*amcha*p2p3*
     -     p3p4**2)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1pI*p2p3*p3p4**2)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1pI*p2p3*
     -     p3p4**2)/amw**4 + 
     -  8*g1L2*g2R2*g2R1*g3L**2*m1*m4*p1p2*p3pI + 
     -  8*g1R2*g2L1*g2L2*g3L**2*m1*m4*p1p2*p3pI - 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p3pI)/
     -   amw**2 - (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p2*
     -     p3pI)/amw**2 - 
     -  8*g1L2*g2L1*g2R2*g3L**2*p1p2*p1p4*p3pI - 
     -  8*g1R2*g2R1*g2L2*g3L**2*p1p2*p1p4*p3pI + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p1p4*p3pI)/amw**2 + 
     -  (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p1p4*p3pI)/amw**2 - 
     -  (4*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2p3*p3pI)/
     -   amw**2 - (4*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*
     -     p3pI)/amw**2 + 
     -  (4*g1L2*g2R2*g2R1*g3L**2*m1*m3**4*m4*p2p3*p3pI)/
     -   amw**4 + (4*g1R2*g2L1*g2L2*g3L**2*m1*m3**4*m4*p2p3*
     -     p3pI)/amw**4 + 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p1p2*p2p3*p3pI)/
     -   amw**4 + (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p2*
     -     p2p3*p3pI)/amw**4 + 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p1p3*p2p3*p3pI)/
     -   amw**4 + (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p3*
     -     p2p3*p3pI)/amw**4 + 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2p3*p3pI)/amw**2 + 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2p3*p3pI)/amw**2 - 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**4*p1p4*p2p3*p3pI)/amw**4 - 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**4*p1p4*p2p3*p3pI)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p1p4*p2p3*p3pI)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p1p4*p2p3*
     -     p3pI)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3*p1p4*p2p3*p3pI)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3*p1p4*p2p3*
     -     p3pI)/amw**4 + 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2p3**2*p3pI)/
     -   amw**4 + (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3**2*
     -     p3pI)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2p3**2*p3pI)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2p3**2*
     -     p3pI)/amw**4 + 
     -  8*g1L2*g2R2*g2R1*g3L**2*m1*m4*p2p4*p3pI + 
     -  8*g1R2*g2L1*g2L2*g3L**2*m1*m4*p2p4*p3pI - 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2p4*p3pI)/
     -   amw**2 - (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p4*
     -     p3pI)/amw**2 - 
     -  8*g1L2*g2L1*g2R2*g3L**2*p1p4*p2p4*p3pI - 
     -  8*g1R2*g2R1*g2L2*g3L**2*p1p4*p2p4*p3pI + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2p4*p3pI)/amw**2 + 
     -  (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2p4*p3pI)/amw**2 + 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2p3*p2p4*p3pI)/
     -   amw**4 + (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*
     -     p2p4*p3pI)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2p3*p2p4*p3pI)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2p3*p2p4*
     -     p3pI)/amw**4 + 
     -  (8*g1L2*g2R2*g2R1*g3L**2*m1*m3**2*m4*p2p3*p3p4*p3pI)/
     -   amw**4 + (8*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*
     -     p3p4*p3pI)/amw**4 - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p4*p2p3*p3p4*p3pI)/
     -   amw**4 - (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p4*p2p3*p3p4*
     -     p3pI)/amw**4 + 
     -  4*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p4pI + 
     -  4*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p4pI - 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**4*p1p2*p4pI)/amw**2 - 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**4*p1p2*p4pI)/amw**2 - 
     -  (16*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2**2*p4pI)/amw**2 - 
     -  (16*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2**2*p4pI)/amw**2 + 
     -  16*g1L2*g2L1*g2R2*g3L**2*p1p2*p1p3*p4pI + 
     -  16*g1R2*g2R1*g2L2*g3L**2*p1p2*p1p3*p4pI - 
     -  (16*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p1p3*p4pI)/
     -   amw**2 - (16*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p1p3*
     -     p4pI)/amw**2 - 
     -  8*g1L2*g2L1*g2R2*g3L**2*m1**2*p2p3*p4pI - 
     -  8*g1R2*g2R1*g2L2*g3L**2*m1**2*p2p3*p4pI - 
     -  8*g1L2*g2R2*g2R1*g3L**2*m1*m4*p2p3*p4pI - 
     -  8*g1R2*g2L1*g2L2*g3L**2*m1*m4*p2p3*p4pI - 
     -  (12*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p2p3*p4pI)/
     -   amw**2 - (12*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p2p3*
     -     p4pI)/amw**2 + 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**4*p1p2*p2p3*p4pI)/amw**4 + 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**4*p1p2*p2p3*p4pI)/amw**4 + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2**2*p2p3*p4pI)/
     -   amw**4 + (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2**2*p2p3*
     -     p4pI)/amw**4 - 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3*p2p3*p4pI)/amw**2 - 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3*p2p3*p4pI)/amw**2 + 
     -  (4*g1L2*g2L1*g2R2*g3L**2*m3**4*p1p3*p2p3*p4pI)/amw**4 + 
     -  (4*g1R2*g2R1*g2L2*g3L**2*m3**4*p1p3*p2p3*p4pI)/amw**4 + 
     -  (16*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p1p3*p2p3*p4pI)/
     -   amw**4 + (16*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p1p3*
     -     p2p3*p4pI)/amw**4 + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3**2*p2p3*p4pI)/
     -   amw**4 + (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3**2*p2p3*
     -     p4pI)/amw**4 + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p2p3**2*p4pI)/
     -   amw**4 + (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p2p3**2*
     -     p4pI)/amw**4 + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3*p2p3**2*p4pI)/
     -   amw**4 + (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3*p2p3**2*
     -     p4pI)/amw**4 - 
     -  (16*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p2p4*p4pI)/
     -   amw**2 - (16*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p2p4*
     -     p4pI)/amw**2 + 
     -  8*g1L2*g2L1*g2R2*g3L**2*p1p3*p2p4*p4pI + 
     -  8*g1R2*g2R1*g2L2*g3L**2*p1p3*p2p4*p4pI - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3*p2p4*p4pI)/amw**2 - 
     -  (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3*p2p4*p4pI)/amw**2 + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p2p3*p2p4*p4pI)/
     -   amw**4 + (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p2p3*p2p4*
     -     p4pI)/amw**4 + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3*p2p3*p2p4*p4pI)/
     -   amw**4 + (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3*p2p3*p2p4*
     -     p4pI)/amw**4 + 
     -  8*g1L2*g2L1*g2R2*g3L**2*p1p2*p3p4*p4pI + 
     -  8*g1R2*g2R1*g2L2*g3L**2*p1p2*p3p4*p4pI - 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p3p4*p4pI)/amw**2 - 
     -  (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p3p4*p4pI)/amw**2 + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p2*p2p3*p3p4*p4pI)/
     -   amw**4 + (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p2*p2p3*p3p4*
     -     p4pI)/amw**4 + 
     -  (8*g1L2*g2L1*g2R2*g3L**2*m3**2*p1p3*p2p3*p3p4*p4pI)/
     -   amw**4 + (8*g1R2*g2R1*g2L2*g3L**2*m3**2*p1p3*p2p3*p3p4*
     -     p4pI)/amw**4)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c---- Interference term W SM fermion- W squark diagram
c---- g1 kopplung squark squark W
c---- g1l2, g1r2 Kopplung squark neutralino quark
c---  g2l1, g2r1 neutralino squark quark Kopplung
c---- g2l2,  W fermion fermion Kopplung
c---- g3l Vector boson Kopplung fermionen


       double precision function ampintWferWsquark(g1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3l,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfer)
      implicit none
      double precision g1, g1l2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,p3p4, amcha,
     .g3l2, g3r1, g3r2, mfer, mv, g3l
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintWferWsquark=
     - -4*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p2 + 
     -  (4*g1*g1L2*g2L2*g2R1*g3L**2*m3**4*m4*mfer*p1p2)/amw**2 + 
     -  (16*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p2**2)/
     -   amw**2 - 16*g1*g1L2*g2L2*g2R1*g3L**2*m4*mfer*p1p2*p1p3 + 
     -  (16*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p2*p1p3)/
     -   amw**2 + 8*g1*g1L2*g2L2*g2R1*g3L**2*m1**2*m4*mfer*p2p3 - 
     -  8*g1*g1L2*g2L1*g2L2*g3L**2*m1*m4**2*mfer*p2p3 + 
     -  (12*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p2*p2p3)/
     -   amw**2 - (4*g1*g1L2*g2L2*g2R1*g3L**2*m3**4*m4*mfer*p1p2*
     -     p2p3)/amw**4 - 
     -  (8*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p2**2*p2p3)/
     -   amw**4 + (4*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p3*
     -     p2p3)/amw**2 - 
     -  (4*g1*g1L2*g2L2*g2R1*g3L**2*m3**4*m4*mfer*p1p3*p2p3)/
     -   amw**4 - (16*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p2*
     -     p1p3*p2p3)/amw**4 - 
     -  (8*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p3**2*p2p3)/
     -   amw**4 - 8*g1*g1L2*g2L1*g2L2*g3L**2*m1*mfer*p1p4*p2p3 + 
     -  8*g1*g1L2*g2L2*g2R1*g3L**2*m4*mfer*p1p4*p2p3 + 
     -  8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m4*p1pfer*p2p3 + 
     -  8*g1*g1R2*g2L2*g2R1*g3L**2*m4**2*p1pfer*p2p3 - 
     -  (8*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p2*p2p3**2)/
     -   amw**4 - (8*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p3*
     -     p2p3**2)/amw**4 + 
     -  4*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p2p4 - 
     -  (4*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**4*mfer*p2p4)/amw**2 - 
     -  (16*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p1p2*p2p4)/
     -   amw**2 + (16*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p2*
     -     p2p4)/amw**2 + 
     -  8*g1*g1L2*g2L1*g2L2*g3L**2*m1*mfer*p1p3*p2p4 - 
     -  (8*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p1p3*p2p4)/
     -   amw**2 - 8*g1*g1L2*g2L2*g2R1*g3L**2*m4*mfer*p1p3*p2p4 + 
     -  (8*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p3*p2p4)/
     -   amw**2 - 4*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*
     -   p2p4 + (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1pfer*p2p4)/
     -   amw**2 + (16*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*
     -     p1pfer*p2p4)/amw**2 - 
     -  8*g1*g1R2*g2L2*g2R1*g3L**2*p1p3*p1pfer*p2p4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p1pfer*p2p4)/
     -   amw**2 - (12*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p2p3*
     -     p2p4)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**4*mfer*p2p3*p2p4)/
     -   amw**4 + (8*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p1p2*
     -     p2p3*p2p4)/amw**4 - 
     -  (8*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p2*p2p3*p2p4)/
     -   amw**4 + (8*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p1p3*
     -     p2p3*p2p4)/amw**4 - 
     -  (8*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p3*p2p3*p2p4)/
     -   amw**4 + (12*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*
     -     p2p3*p2p4)/amw**2 - 
     -  (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1pfer*p2p3*p2p4)/
     -   amw**4 - (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p1pfer*
     -     p2p3*p2p4)/amw**4 - 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p1pfer*p2p3*
     -     p2p4)/amw**4 + 
     -  (8*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p2p3**2*p2p4)/
     -   amw**4 - (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*
     -     p2p3**2*p2p4)/amw**4 - 
     -  (16*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p2p4**2)/
     -   amw**2 + (16*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*
     -     p2p4**2)/amw**2 + 
     -  (8*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p2p3*p2p4**2)/
     -   amw**4 - (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*p2p3*
     -     p2p4**2)/amw**4 - 
     -  4*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2pfer + 
     -  (4*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**4*m4*p2pfer)/
     -   amw**2 + (16*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p2*
     -     p2pfer)/amw**2 - 
     -  8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m4*p1p3*p2pfer + 
     -  (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p3*p2pfer)/
     -   amw**2 + 4*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*
     -   p2pfer - (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1p4*
     -     p2pfer)/amw**2 - 
     -  (16*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p1p4*p2pfer)/
     -   amw**2 + 8*g1*g1R2*g2L2*g2R1*g3L**2*p1p3*p1p4*p2pfer - 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p1p4*p2pfer)/
     -   amw**2 + (12*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*
     -     p2pfer)/amw**2 - 
     -  (4*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**4*m4*p2p3*p2pfer)/
     -   amw**4 - (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p2*
     -     p2p3*p2pfer)/amw**4 - 
     -  (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p3*p2p3*
     -     p2pfer)/amw**4 - 
     -  (12*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p3*p2pfer)/
     -   amw**2 + (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1p4*p2p3*
     -     p2pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p1p4*p2p3*
     -     p2pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p1p4*p2p3*
     -     p2pfer)/amw**4 - 
     -  (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3**2*p2pfer)/
     -   amw**4 + (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*
     -     p2p3**2*p2pfer)/amw**4 + 
     -  (16*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p4*p2pfer)/
     -   amw**2 - (16*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p4*
     -     p2pfer)/amw**2 - 
     -  (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*p2p4*
     -     p2pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p3*p2p4*
     -     p2pfer)/amw**4 + 
     -  8*g1*g1L2*g2L1*g2L2*g3L**2*m1*mfer*p1p2*p3p4 - 
     -  (8*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p1p2*p3p4)/
     -   amw**2 - 8*g1*g1L2*g2L2*g2R1*g3L**2*m4*mfer*p1p2*p3p4 + 
     -  (8*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p2*p3p4)/
     -   amw**2 - 8*g1*g1R2*g2L2*g2R1*g3L**2*p1p2*p1pfer*p3p4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p1pfer*p3p4)/
     -   amw**2 - (4*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p2p3*
     -     p3p4)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**4*mfer*p2p3*p3p4)/
     -   amw**4 + (8*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p1p2*
     -     p2p3*p3p4)/amw**4 - 
     -  (8*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p2*p2p3*p3p4)/
     -   amw**4 + (8*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p1p3*
     -     p2p3*p3p4)/amw**4 - 
     -  (8*g1*g1L2*g2L2*g2R1*g3L**2*m3**2*m4*mfer*p1p3*p2p3*p3p4)/
     -   amw**4 + (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*p2p3*
     -     p3p4)/amw**2 - 
     -  (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1pfer*p2p3*p3p4)/
     -   amw**4 - (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p1pfer*
     -     p2p3*p3p4)/amw**4 - 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p1pfer*p2p3*
     -     p3p4)/amw**4 + 
     -  (8*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p2p3**2*p3p4)/
     -   amw**4 - (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*
     -     p2p3**2*p3p4)/amw**4 + 
     -  16*g1*g1L2*g2L1*g2L2*g3L**2*m1*mfer*p2p4*p3p4 - 
     -  (16*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p2p4*p3p4)/
     -   amw**2 - 16*g1*g1R2*g2L2*g2R1*g3L**2*p1pfer*p2p4*
     -   p3p4 + (16*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*p2p4*
     -     p3p4)/amw**2 + 
     -  (16*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p2p3*p2p4*
     -     p3p4)/amw**4 - 
     -  (16*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*p2p3*p2p4*
     -     p3p4)/amw**4 - 
     -  8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m4*p2pfer*p3p4 + 
     -  (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2pfer*p3p4)/
     -   amw**2 + 8*g1*g1R2*g2L2*g2R1*g3L**2*p1p4*p2pfer*p3p4 - 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2pfer*p3p4)/
     -   amw**2 - (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*
     -     p2pfer*p3p4)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p3*p2pfer*
     -     p3p4)/amw**4 + 
     -  (8*g1*g1L2*g2L1*g2L2*g3L**2*m1*m3**2*mfer*p2p3*p3p4**2)/
     -   amw**4 - (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1pfer*p2p3*
     -     p3p4**2)/amw**4 - 
     -  8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m4*p1p2*p3pfer + 
     -  (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p2*p3pfer)/
     -   amw**2 + 8*g1*g1R2*g2L2*g2R1*g3L**2*p1p2*p1p4*p3pfer - 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p1p4*p3pfer)/
     -   amw**2 + (4*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*
     -     p3pfer)/amw**2 - 
     -  (4*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**4*m4*p2p3*p3pfer)/
     -   amw**4 - (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p2*
     -     p2p3*p3pfer)/amw**4 - 
     -  (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p1p3*p2p3*
     -     p3pfer)/amw**4 - 
     -  (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p3*p3pfer)/
     -   amw**2 + (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1p4*p2p3*
     -     p3pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p1p4*p2p3*
     -     p3pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p1p4*p2p3*
     -     p3pfer)/amw**4 - 
     -  (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3**2*p3pfer)/
     -   amw**4 + (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*
     -     p2p3**2*p3pfer)/amw**4 - 
     -  8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m4*p2p4*p3pfer + 
     -  (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p4*p3pfer)/
     -   amw**2 + 8*g1*g1R2*g2L2*g2R1*g3L**2*p1p4*p2p4*p3pfer - 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p4*p3pfer)/
     -   amw**2 - (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*
     -     p2p4*p3pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p3*p2p4*
     -     p3pfer)/amw**4 - 
     -  (8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m3**2*m4*p2p3*p3p4*
     -     p3pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p4*p2p3*p3p4*
     -     p3pfer)/amw**4 + 
     -  4*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p4pfer - 
     -  (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1p2*p4pfer)/amw**2 - 
     -  (16*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2**2*p4pfer)/
     -   amw**2 + 16*g1*g1R2*g2L2*g2R1*g3L**2*p1p2*p1p3*
     -   p4pfer - (16*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p1p3*
     -     p4pfer)/amw**2 - 
     -  8*g1*g1R2*g2L2*g2R1*g3L**2*m1**2*p2p3*p4pfer + 
     -  8*g1*g1R2*g2L1*g2L2*g3L**2*m1*m4*p2p3*p4pfer - 
     -  (12*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p2p3*p4pfer)/
     -   amw**2 + (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1p2*p2p3*
     -     p4pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2**2*p2p3*p4pfer)/
     -   amw**4 - (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p2p3*
     -     p4pfer)/amw**2 + 
     -  (4*g1*g1R2*g2L2*g2R1*g3L**2*m3**4*p1p3*p2p3*p4pfer)/
     -   amw**4 + (16*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p1p3*
     -     p2p3*p4pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3**2*p2p3*p4pfer)/
     -   amw**4 - 16*g1*g1R2*g2L2*g2R1*g3L**2*p1p4*p2p3*
     -   p4pfer + (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*
     -     p2p3**2*p4pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p2p3**2*p4pfer)/
     -   amw**4 - (16*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p2p4*
     -     p4pfer)/amw**2 + 
     -  8*g1*g1R2*g2L2*g2R1*g3L**2*p1p3*p2p4*p4pfer - 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p2p4*p4pfer)/
     -   amw**2 + (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p2p3*
     -     p2p4*p4pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p2p3*p2p4*
     -     p4pfer)/amw**4 + 
     -  8*g1*g1R2*g2L2*g2R1*g3L**2*p1p2*p3p4*p4pfer - 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p3p4*p4pfer)/
     -   amw**2 + (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p2*p2p3*
     -     p3p4*p4pfer)/amw**4 + 
     -  (8*g1*g1R2*g2L2*g2R1*g3L**2*m3**2*p1p3*p2p3*p3p4*
     -     p4pfer)/amw**4
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- g1= W squark squark coupling
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= neutralino quark squark
c---- g2l2, g2r2= slepton chargino neutrino
c---- g3l1 W fermion fermion
c---- g3l2, g3r2= slepton lepton neutralino
       double precision function ampintcharslepWsquark(g1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amcha)
      implicit none
      double precision g1, g1l2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha,
     .g3l2, g3r1, g3r2, mfer, mv
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2
      
	mI=amcha

      ampintcharslepWsquark=
     - -(-4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mI*p1p2 - 
     -  2*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m3**2*m4*mI*p1p2 + 
     -  (2*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m3**4*m4*mI*p1p2)/
     -   amw**2 + 4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*m4**2*mI*
     -   p1p2 + (4*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m3**2*m4*mI*
     -     p1p2**2)/amw**2 - 
     -  8*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m4*mI*p1p2*p1p3 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m3**2*m4*mI*p1p2*p1p3)/
     -   amw**2 + 8*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p2*
     -   p1p4 + 4*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m1**2*m4*mI*
     -   p2p3 - 2*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mI*
     -   p2p3 + (2*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3**3*m4*mI*
     -     p2p3)/amw**2 - 
     -  4*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*m4**2*mI*p2p3 + 
     -  (4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mI*p1p2*p2p3)/
     -   amw**2 + (4*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m3**2*m4*mI*
     -     p1p2*p2p3)/amw**2 + 
     -  (4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mI*p1p3*p2p3)/
     -   amw**2 - 4*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*mI*p1p4*
     -   p2p3 + 2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p4*
     -   p2p3 - (2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*mI*p1p4*
     -     p2p3)/amw**2 + 
     -  4*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m4*mI*p1p4*p2p3 - 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p2*p1p4*p2p3)/
     -   amw**2 - (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p3*
     -     p1p4*p2p3)/amw**2 - 
     -  4*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m4*p1pI*p2p3 + 
     -  2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1pI*p2p3 - 
     -  (2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3**3*m4*p1pI*p2p3)/
     -   amw**2 + 4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m4**2*p1pI*
     -   p2p3 - (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*
     -     p1pI*p2p3)/amw**2 - 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3*p1pI*p2p3)/
     -   amw**2 + 8*g1*g1L2*g2L1*g2R2*g3L1*g3L2*p1p4*p1pI*
     -   p2p3 + (4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mI*
     -     p2p3**2)/amw**2 - 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p4*p2p3**2)/
     -   amw**2 - (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1pI*
     -     p2p3**2)/amw**2 - 
     -  4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m1**2*m3*mI*p2p4 + 
     -  2*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*m3**2*mI*p2p4 - 
     -  (2*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*m3**4*mI*p2p4)/
     -   amw**2 - 4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mI*
     -   p2p4 - (4*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*m3**2*mI*
     -     p1p2*p2p4)/amw**2 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m3**2*m4*mI*p1p2*p2p4)/
     -   amw**2 + 4*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*mI*p1p3*
     -   p2p4 - 2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p3*
     -   p2p4 - (4*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*m3**2*mI*
     -     p1p3*p2p4)/amw**2 + 
     -  (2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*mI*p1p3*p2p4)/
     -   amw**2 - 4*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m4*mI*p1p3*
     -   p2p4 + (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p2*
     -     p1p3*p2p4)/amw**2 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p3**2*p2p4)/
     -   amw**2 + 4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1pI*
     -   p2p4 - 2*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1pI*
     -   p2p4 + (2*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**4*p1pI*
     -     p2p4)/amw**2 + 
     -  4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1pI*p2p4 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p2*p1pI*p2p4)/
     -   amw**2 - 8*g1*g1L2*g2L1*g2R2*g3L1*g3L2*p1p3*p1pI*
     -   p2p4 + (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p3*
     -     p1pI*p2p4)/amw**2 - 
     -  (4*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*m3**2*mI*p2p3*p2p4)/
     -   amw**2 + (4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mI*
     -     p2p3*p2p4)/amw**2 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p3*p2p3*p2p4)/
     -   amw**2 - (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p4*
     -     p2p3*p2p4)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1pI*p2p3*p2p4)/
     -   amw**2 - (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1pI*
     -     p2p3*p2p4)/amw**2 - 
     -  (4*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*m3**2*mI*p2p4**2)/
     -   amw**2 + (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p3*
     -     p2p4**2)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1pI*p2p4**2)/
     -   amw**2 + 4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m1**2*m3*m4*
     -   p2pI + 2*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3**2*m4*
     -   p2pI - (2*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3**4*m4*
     -     p2pI)/amw**2 - 
     -  4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*m4**2*p2pI - 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3**2*m4*p1p2*p2pI)/
     -   amw**2 + 4*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m4*p1p3*
     -   p2pI + 2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3*
     -   p2pI - (4*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3**2*m4*
     -     p1p3*p2pI)/amw**2 - 
     -  (2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3**3*m4*p1p3*p2pI)/
     -   amw**2 + 4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m4**2*p1p3*
     -   p2pI - (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*
     -     p1p3*p2pI)/amw**2 - 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3**2*p2pI)/
     -   amw**2 - 4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p4*
     -   p2pI - 2*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p4*
     -   p2pI + (2*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**4*p1p4*
     -     p2pI)/amw**2 + 
     -  4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p4*p2pI + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p2*p1p4*p2pI)/
     -   amw**2 + (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p3*
     -     p1p4*p2pI)/amw**2 - 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3**2*m4*p2p3*p2pI)/
     -   amw**2 - (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3*
     -     p2p3*p2pI)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p4*p2p3*p2pI)/
     -   amw**2 - (4*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3**2*m4*
     -     p2p4*p2pI)/amw**2 - 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3*p2p4*p2pI)/
     -   amw**2 + (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p4*
     -     p2p4*p2pI)/amw**2 + 
     -  4*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*mI*p1p2*p3p4 + 
     -  2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p2*p3p4 - 
     -  (2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*mI*p1p2*p3p4)/
     -   amw**2 - 4*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m4*mI*p1p2*
     -   p3p4 + (4*g1*g1R2*g2R1*g2R2*g3L1*g3L2*m3**2*m4*mI*
     -     p1p2*p3p4)/amw**2 - 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p2**2*p3p4)/
     -   amw**2 - (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p2*
     -     p1p3*p3p4)/amw**2 + 
     -  (4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*mI*p2p3*p3p4)/
     -   amw**2 - (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p2*
     -     p2p3*p3p4)/amw**2 - 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p4*p2p3*p3p4)/
     -   amw**2 - (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1pI*
     -     p2p3*p3p4)/amw**2 + 
     -  8*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*mI*p2p4*p3p4 - 
     -  (4*g1*g1R2*g2L1*g2R2*g3L1*g3L2*m1*m3**2*mI*p2p4*p3p4)/
     -   amw**2 - (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p2*
     -     p2p4*p3p4)/amw**2 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p3*p2p4*p3p4)/
     -   amw**2 - 8*g1*g1L2*g2L1*g2R2*g3L1*g3L2*p1pI*p2p4*
     -   p3p4 + (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1pI*
     -     p2p4*p3p4)/amw**2 - 
     -  4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m1**2*p2pI*p3p4 - 
     -  2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2pI*p3p4 + 
     -  (2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2pI*p3p4)/
     -   amw**2 + 4*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m4*p2pI*
     -   p3p4 - (4*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3**2*m4*
     -     p2pI*p3p4)/amw**2 + 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p2*p2pI*p3p4)/
     -   amw**2 + (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p3*
     -     p2pI*p3p4)/amw**2 - 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3*p2pI*p3p4)/
     -   amw**2 - 8*g1*g1L2*g2L1*g2R2*g3L1*g3L2*p1p4*p2pI*
     -   p3p4 + (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p4*
     -     p2pI*p3p4)/amw**2 + 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*p2pI*p3p4)/
     -   amw**2 + (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*
     -     p2pI*p3p4)/amw**2 - 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*mI*p1p2*p3p4**2)/
     -   amw**2 + (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2pI*
     -     p3p4**2)/amw**2 + 
     -  4*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m4*p1p2*p3pI - 
     -  2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*p3pI + 
     -  (2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3**3*m4*p1p2*p3pI)/
     -   amw**2 - 4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m4**2*p1p2*
     -   p3pI + (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2**2*
     -     p3pI)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*p1p3*p3pI)/
     -   amw**2 - 8*g1*g1L2*g2L1*g2R2*g3L1*g3L2*p1p2*p1p4*
     -   p3pI + (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*
     -     p2p3*p3pI)/amw**2 + 
     -  4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m1**2*p2p4*p3pI + 
     -  2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*p3pI - 
     -  (2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2p4*p3pI)/
     -   amw**2 + 4*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m4*p2p4*
     -   p3pI - (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p2*
     -     p2p4*p3pI)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*p2p4*p3pI)/
     -   amw**2 - (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p3*
     -     p2p4*p3pI)/amw**2 - 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*p2p4*p3pI)/
     -   amw**2 - (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*
     -     p2p4**2*p3pI)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*p3p4*p3pI)/
     -   amw**2 - (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*
     -     p3p4*p3pI)/amw**2 - 
     -  4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p2*p4pI + 
     -  2*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p2*p4pI - 
     -  (2*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**4*p1p2*p4pI)/
     -   amw**2 - 4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*
     -   p4pI - (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p2**2*
     -     p4pI)/amw**2 + 
     -  8*g1*g1L2*g2L1*g2R2*g3L1*g3L2*p1p2*p1p3*p4pI - 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p2*p1p3*p4pI)/
     -   amw**2 - 4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m1**2*p2p3*
     -   p4pI - 2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*
     -   p4pI + (2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2p3*
     -     p4pI)/amw**2 - 
     -  4*g1*g1L2*g2R1*g2R2*g3L1*g3L2*m1*m4*p2p3*p4pI + 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p2*p2p3*p4pI)/
     -   amw**2 - (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p2*
     -     p2p3*p4pI)/amw**2 + 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p3*p2p3*p4pI)/
     -   amw**2 + (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*
     -     p2p3**2*p4pI)/amw**2 - 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p2*p2p4*p4pI)/
     -   amw**2 + (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*
     -     p2p4*p4pI)/amw**2 + 
     -  8*g1*g1L2*g2L1*g2R2*g3L1*g3L2*p1p2*p3p4*p4pI - 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3L2*m3**2*p1p2*p3p4*p4pI)/
     -   amw**2 + (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*
     -     p3p4*p4pI)/amw**2)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- g1= W squark squark coupling
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= neutralino quark squark
c---- g2l2, g2r2= sneutrino chargino lepton
c---- g3l1 W fermion fermion
c---- g3l2, g3r2= sneutrino neutrino neutralino
       double precision function ampintcharsneutWsquark(g1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3l1,g3l2, g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amcha)
      implicit none
      double precision g1, g1l2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha,
     .g3l2, g3r1, g3r2, mfer, mv
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintcharsneutWsquark=
     -   4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*amcha*p1p2 - 
     -  2*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m3**2*m4*amcha*p1p2 + 
     -  (2*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m3**4*m4*amcha*p1p2)/
     -   amw**2 + 4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*m4**2*amcha*
     -   p1p2 + (4*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m3**2*m4*amcha*
     -     p1p2**2)/amw**2 - 
     -  8*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m4*amcha*p1p2*p1p3 + 
     -  (4*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m3**2*m4*amcha*p1p2*p1p3)/
     -   amw**2 - 8*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*
     -   p1pI + 4*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m1**2*m4*amcha*
     -   p2p3 + 2*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*amcha*
     -   p2p3 - (2*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3**3*m4*amcha*
     -     p2p3)/amw**2 - 
     -  4*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*m4**2*amcha*p2p3 - 
     -  (4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*amcha*p1p2*p2p3)/
     -   amw**2 + (4*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m3**2*m4*amcha*
     -     p1p2*p2p3)/amw**2 - 
     -  (4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*amcha*p1p3*p2p3)/
     -   amw**2 - 4*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*amcha*p1p4*
     -   p2p3 - 2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p4*p2p3 + 
     -  (2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*amcha*p1p4*p2p3)/
     -   amw**2 + 4*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m4*amcha*p1p4*
     -   p2p3 + (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2*p1p4*
     -     p2p3)/amw**2 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p3*p1p4*p2p3)/
     -   amw**2 - 4*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m4*p1pI*
     -   p2p3 - 2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1pI*p2p3 + 
     -  (2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3**3*m4*p1pI*p2p3)/
     -   amw**2 + 4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m4**2*p1pI*
     -   p2p3 + (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*p1pI*
     -     p2p3)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3*p1pI*p2p3)/
     -   amw**2 + 8*g1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p4*p1pI*
     -   p2p3 - (4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*amcha*
     -     p2p3**2)/amw**2 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p4*p2p3**2)/
     -   amw**2 + (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1pI*
     -     p2p3**2)/amw**2 - 
     -  4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m1**2*m3*amcha*p2p4 + 
     -  2*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*m3**2*amcha*p2p4 - 
     -  (2*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*m3**4*amcha*p2p4)/
     -   amw**2 + 4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*amcha*
     -   p2p4 - (4*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*m3**2*amcha*p1p2*
     -     p2p4)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m3**2*m4*amcha*p1p2*p2p4)/
     -   amw**2 + 4*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*amcha*p1p3*
     -   p2p4 - 2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p3*p2p4 - 
     -  (4*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*m3**2*amcha*p1p3*p2p4)/
     -   amw**2 + (2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*amcha*p1p3*
     -     p2p4)/amw**2 - 
     -  4*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m4*amcha*p1p3*p2p4 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2*p1p3*p2p4)/
     -   amw**2 + (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p3**2*
     -     p2p4)/amw**2 - 
     -  8*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p4*p2p4 + 
     -  4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1pI*p2p4 - 
     -  2*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1pI*p2p4 + 
     -  (2*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**4*p1pI*p2p4)/
     -   amw**2 - 4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1pI*
     -   p2p4 + (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*p1pI*
     -     p2p4)/amw**2 + 
     -  (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p3*p1pI*p2p4)/
     -   amw**2 - (4*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*m3**2*amcha*
     -     p2p3*p2p4)/amw**2 - 
     -  (4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*amcha*p2p3*p2p4)/
     -   amw**2 + (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p3*
     -     p2p3*p2p4)/amw**2 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p4*p2p3*p2p4)/
     -   amw**2 + (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1pI*
     -     p2p3*p2p4)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1pI*p2p3*p2p4)/
     -   amw**2 - (4*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*m3**2*amcha*
     -     p2p4**2)/amw**2 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p3*p2p4**2)/
     -   amw**2 + (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1pI*
     -     p2p4**2)/amw**2 + 
     -  4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m1**2*m3*m4*p2pI + 
     -  2*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m3**2*m4*p2pI - 
     -  (2*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m3**4*m4*p2pI)/
     -   amw**2 - 4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*m4**2*
     -   p2pI - (4*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m3**2*m4*p1p2*
     -     p2pI)/amw**2 + 
     -  4*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m4*p1p3*p2pI + 
     -  2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3*p2pI - 
     -  (4*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m3**2*m4*p1p3*p2pI)/
     -   amw**2 - (2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3**3*m4*p1p3*
     -     p2pI)/amw**2 - 
     -  4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m4**2*p1p3*p2pI - 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*p1p3*p2pI)/
     -   amw**2 - (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3**2*
     -     p2pI)/amw**2 - 
     -  4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p4*p2pI - 
     -  2*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p4*p2pI + 
     -  (2*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**4*p1p4*p2pI)/
     -   amw**2 + 4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p4*
     -   p2pI + (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*p1p4*
     -     p2pI)/amw**2 - 
     -  8*g1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p3*p1p4*p2pI + 
     -  (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p3*p1p4*p2pI)/
     -   amw**2 - (4*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m3**2*m4*
     -     p2p3*p2pI)/amw**2 - 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3*p2p3*p2pI)/
     -   amw**2 + (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p4*
     -     p2p3*p2pI)/amw**2 - 
     -  (4*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m3**2*m4*p2p4*p2pI)/
     -   amw**2 - (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3*
     -     p2p4*p2pI)/amw**2 + 
     -  (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p4*p2p4*p2pI)/
     -   amw**2 + 4*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*amcha*p1p2*
     -   p3p4 + 2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2*p3p4 - 
     -  (2*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3**3*amcha*p1p2*p3p4)/
     -   amw**2 - 4*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m4*amcha*p1p2*
     -   p3p4 + (4*g1*g1L2*g2L1*g2L2*g3L1*g3R2*m3**2*m4*amcha*p1p2*
     -     p3p4)/amw**2 - 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2**2*p3p4)/
     -   amw**2 - (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2*
     -     p1p3*p3p4)/amw**2 - 
     -  8*g1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p2*p1pI*p3p4 - 
     -  (4*g1*g1R2*g2L1*g2R2*g3L1*g3R2*m1*m3*m4*amcha*p2p3*p3p4)/
     -   amw**2 - (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2*
     -     p2p3*p3p4)/amw**2 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p4*p2p3*p3p4)/
     -   amw**2 + (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1pI*
     -     p2p3*p3p4)/amw**2 + 
     -  8*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*amcha*p2p4*p3p4 - 
     -  (4*g1*g1L2*g2L2*g2R1*g3L1*g3R2*m1*m3**2*amcha*p2p4*p3p4)/
     -   amw**2 - (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2*
     -     p2p4*p3p4)/amw**2 + 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p3*p2p4*p3p4)/
     -   amw**2 - 8*g1*g1R2*g2L2*g2R1*g3L1*g3R2*p1pI*p2p4*
     -   p3p4 + (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1pI*p2p4*
     -     p3p4)/amw**2 + 
     -  4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m1**2*p2pI*p3p4 - 
     -  2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2pI*p3p4 + 
     -  (2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2pI*p3p4)/
     -   amw**2 + 4*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m4*p2pI*
     -   p3p4 - (4*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m3**2*m4*p2pI*
     -     p3p4)/amw**2 + 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p2*p2pI*p3p4)/
     -   amw**2 + (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p3*
     -     p2pI*p3p4)/amw**2 - 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p3*p2pI*p3p4)/
     -   amw**2 + (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p4*
     -     p2pI*p3p4)/amw**2 + 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*p2pI*p3p4)/
     -   amw**2 + (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*
     -     p2pI*p3p4)/amw**2 - 
     -  (4*g1*g1R2*g2R1*g2R2*g3L1*g3R2*m3*amcha*p1p2*p3p4**2)/
     -   amw**2 + (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2pI*
     -     p3p4**2)/amw**2 + 
     -  4*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m4*p1p2*p3pI - 
     -  2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*p3pI + 
     -  (2*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3**3*m4*p1p2*p3pI)/
     -   amw**2 + 4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m4**2*p1p2*
     -   p3pI + (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2**2*
     -     p3pI)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*p1p3*p3pI)/
     -   amw**2 + (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*
     -     p2p3*p3pI)/amw**2 - 
     -  4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m1**2*p2p4*p3pI + 
     -  2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*p3pI - 
     -  (2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2p4*p3pI)/
     -   amw**2 + 4*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m4*p2p4*
     -   p3pI - (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p2*p2p4*
     -     p3pI)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*p2p4*p3pI)/
     -   amw**2 - (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p3*
     -     p2p4*p3pI)/amw**2 - 
     -  8*g1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p4*p2p4*p3pI - 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*p2p4*p3pI)/
     -   amw**2 - (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4**2*
     -     p3pI)/amw**2 + 
     -  (4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*p3p4*p3pI)/
     -   amw**2 - (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*
     -     p3p4*p3pI)/amw**2 + 
     -  4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p2*p4pI + 
     -  2*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*p4pI - 
     -  (2*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**4*p1p2*p4pI)/
     -   amw**2 - 4*g1*g1L2*g2L1*g2R2*g3L1*g3R2*m3*m4*p1p2*
     -   p4pI - (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2**2*
     -     p4pI)/amw**2 + 
     -  8*g1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p2*p1p3*p4pI - 
     -  (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*p1p3*p4pI)/
     -   amw**2 - 4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m1**2*p2p3*
     -   p4pI + 2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*p4pI - 
     -  (2*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3**3*p2p3*p4pI)/
     -   amw**2 - 4*g1*g1R2*g2L1*g2L2*g3L1*g3R2*m1*m4*p2p3*
     -   p4pI - (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p2*p2p3*
     -     p4pI)/amw**2 - 
     -  (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*p2p3*p4pI)/
     -   amw**2 - (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p1p3*
     -     p2p3*p4pI)/amw**2 - 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3**2*p4pI)/
     -   amw**2 + 8*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p4*
     -   p4pI - (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*p2p4*
     -     p4pI)/amw**2 + 
     -  8*g1*g1R2*g2L2*g2R1*g3L1*g3R2*p1p3*p2p4*p4pI - 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*p2p4*p4pI)/
     -   amw**2 - (4*g1*g1R2*g2L2*g2R1*g3L1*g3R2*m3**2*p1p2*
     -     p3p4*p4pI)/amw**2 - 
     -  (4*g1*g1L2*g2R1*g2R2*g3L1*g3R2*m1*m3*p2p3*p3p4*p4pI)/
     -   amw**2
       return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- g1=charged Higgs squark squark coupling
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= neutralino quark squark
c---- g2l2, g2r2= neutralino chargino W
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 W fermion fermion

       double precision function ampintchaHisquarkWchar(g1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amcha)
      implicit none
      double precision g1, g1l2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha,
     .g3l2, g3r1, g3r2, mfer, mv
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2
       

      ampintchaHisquarkWchar=
     - -4*g1R2*g2L2*g2R1*g3L2*g3R1*m3*m4*amcha*p1p2 - 
     -  4*g1L2*g2L1*g2R2*g3L2*g3R1*m3*m4*amcha*p1p2 + 
     -  (4*g1R2*g2L2*g2R1*g3L2*g3R1*m3*m4*amcha*p1p2*p2p3)/
     -   amw**2 + (4*g1L2*g2L1*g2R2*g3L2*g3R1*m3*m4*amcha*p1p2*
     -     p2p3)/amw**2 + 
     -  (4*g1R2*g2L2*g2R1*g3L2*g3R1*m3*m4*amcha*p1p3*p2p3)/
     -   amw**2 + (4*g1L2*g2L1*g2R2*g3L2*g3R1*m3*m4*amcha*p1p3*
     -     p2p3)/amw**2 + 
     -  4*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*amcha*p2p4 + 
     -  4*g1L2*g2R1*g2R2*g3L2*g3R1*m1*m3*amcha*p2p4 - 
     -  4*g1L2*g2L1*g2L2*g3L2*g3R1*m3*p1pI*p2p4 - 
     -  4*g1R2*g2R1*g2R2*g3L2*g3R1*m3*p1pI*p2p4 - 
     -  (4*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*amcha*p2p3*p2p4)/
     -   amw**2 - (4*g1L2*g2R1*g2R2*g3L2*g3R1*m1*m3*amcha*p2p3*
     -     p2p4)/amw**2 + 
     -  (4*g1L2*g2L1*g2L2*g3L2*g3R1*m3*p1pI*p2p3*p2p4)/amw**2 + 
     -  (4*g1R2*g2R1*g2R2*g3L2*g3R1*m3*p1pI*p2p3*p2p4)/amw**2 + 
     -  4*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*p2pI + 
     -  4*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*m4*p2pI - 
     -  4*g1L2*g2L1*g2L2*g3L2*g3R1*m3*p1p4*p2pI - 
     -  4*g1R2*g2R1*g2R2*g3L2*g3R1*m3*p1p4*p2pI - 
     -  (4*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*p2p3*p2pI)/
     -   amw**2 - (4*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*g1L2*g2L1*g2L2*g3L2*g3R1*m3*p1p4*p2p3*p2pI)/amw**2 + 
     -  (4*g1R2*g2R1*g2R2*g3L2*g3R1*m3*p1p4*p2p3*p2pI)/amw**2 - 
     -  (4*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*amcha*p2p3*p3p4)/
     -   amw**2 - (4*g1L2*g2R1*g2R2*g3L2*g3R1*m1*m3*amcha*p2p3*
     -     p3p4)/amw**2 + 
     -  (4*g1L2*g2L1*g2L2*g3L2*g3R1*m3*p1pI*p2p3*p3p4)/amw**2 + 
     -  (4*g1R2*g2R1*g2R2*g3L2*g3R1*m3*p1pI*p2p3*p3p4)/amw**2 - 
     -  (4*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*p2p3*p3pI)/
     -   amw**2 - (4*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p3pI)/amw**2 + 
     -  (4*g1L2*g2L1*g2L2*g3L2*g3R1*m3*p1p4*p2p3*p3pI)/amw**2 + 
     -  (4*g1R2*g2R1*g2R2*g3L2*g3R1*m3*p1p4*p2p3*p3pI)/amw**2 + 
     -  4*g1L2*g2L1*g2L2*g3L2*g3R1*m3*p1p2*p4pI + 
     -  4*g1R2*g2R1*g2R2*g3L2*g3R1*m3*p1p2*p4pI - 
     -  (4*g1L2*g2L1*g2L2*g3L2*g3R1*m3*p1p2*p2p3*p4pI)/amw**2 - 
     -  (4*g1R2*g2R1*g2R2*g3L2*g3R1*m3*p1p2*p2p3*p4pI)/amw**2 - 
     -  (4*g1L2*g2L1*g2L2*g3L2*g3R1*m3*p1p3*p2p3*p4pI)/amw**2 - 
     -  (4*g1R2*g2R1*g2R2*g3L2*g3R1*m3*p1p3*p2p3*p4pI)/amw**2

       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference charged Higgs squark- W Sm fermion diagrams
c----- g1=charged Higgs squark squark coupling
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= neutralino quark squark
c---- g2l2, g2r2= neutralino chargino W
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 W fermion fermion

       double precision function ampintchaHisquarkWSMfer(g1, g1l2, 
     . g1r2,g2l1, g2l2, g2r1,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfer)
      implicit none
      double precision g1, g1l2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,p3p4, amcha,
     .g3l2, g3r1, g3r2, mfer, mv
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHisquarkWSMfer=
     -g1*(   -4*g1L2*g2L2*g2R1*g3L2*g3R1*m3*m4*mfer*p1p2 + 
     -  (4*g1L2*g2L2*g2R1*g3L2*g3R1*m3*m4*mfer*p1p2*p2p3)/
     -   amw**2 + (4*g1L2*g2L2*g2R1*g3L2*g3R1*m3*m4*mfer*p1p3*
     -     p2p3)/amw**2 + 
     -  4*g1L2*g2L1*g2L2*g3L2*g3R1*m1*m3*mfer*p2p4 - 
     -  4*g1R2*g2L1*g2L2*g3L2*g3R1*m3*p1pfer*p2p4 - 
     -  (4*g1L2*g2L1*g2L2*g3L2*g3R1*m1*m3*mfer*p2p3*p2p4)/
     -   amw**2 + (4*g1R2*g2L1*g2L2*g3L2*g3R1*m3*p1pfer*p2p3*
     -     p2p4)/amw**2 - 
     -  4*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*p2pfer + 
     -  4*g1R2*g2L1*g2L2*g3L2*g3R1*m3*p1p4*p2pfer + 
     -  (4*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*p2p3*p2pfer)/
     -   amw**2 - (4*g1R2*g2L1*g2L2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p2pfer)/amw**2 - 
     -  (4*g1L2*g2L1*g2L2*g3L2*g3R1*m1*m3*mfer*p2p3*p3p4)/
     -   amw**2 + (4*g1R2*g2L1*g2L2*g3L2*g3R1*m3*p1pfer*p2p3*
     -     p3p4)/amw**2 + 
     -  (4*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*p2p3*p3pfer)/
     -   amw**2 - (4*g1R2*g2L1*g2L2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p3pfer)/amw**2 + 
     -  4*g1R2*g2L1*g2L2*g3L2*g3R1*m3*p1p2*p4pfer - 
     -  (4*g1R2*g2L1*g2L2*g3L2*g3R1*m3*p1p2*p2p3*p4pfer)/
     -   amw**2 - (4*g1R2*g2L1*g2L2*g3L2*g3R1*m3*p1p3*p2p3*
     -     p4pfer)/amw**2)
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference charged Higgs squark- chargino slepton diagrams
c----- g1=charged Higgs squark squark coupling
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= neutralino quark squark
c---- g2l2, g2r2= slepton neutrino chargino
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 slepton lepton neutralino

       double precision function ampintchaHisquarkslepton(g1, g1l2, 
     . g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amcha)
      implicit none
      double precision g1, g1l2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha,
     .g3l2, g3r1, g3r2, mfer, mv
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHisquarkslepton=
     -   4*g1*(-(g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*amcha*p1p2)/2. - 
     -(g1R2*g2L1*g2R2*g3R1*g3R2*m1*m4*amcha*p2p3)/2. + 
     -(g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*amcha*p2p4)/2. + 
     -(g1L2*g2R1*g2R2*g3L2*g3R1*m1*m3*m4*p2pI)/2. + 
     -(g1R2*g2R1*g2R2*g3R1*g3R2*amcha*(p1p4*p2p3-p1p3*p2p4+p1p2*p3p4))/
     -2. + (g1L2*g2L1*g2R2*g3R1*g3R2*m4*
     - (p1pI*p2p3 + p1p3*p2pI - p1p2*p3pI))/2. - 
     -(g1L2*g2L1*g2R2*g3L2*g3R1*m3*(p1pI*p2p4 + p1p4*p2pI-p1p2*p4pI))/
     -2. - (g1L2*g2R1*g2R2*g3R1*g3R2*m1*
     - (p2pI*p3p4 - p2p4*p3pI + p2p3*p4pI))/2.)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference charged Higgs squark- chargino sneutrino diagrams
c----- g1=charged Higgs squark squark coupling
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= neutralino quark squark
c---- g2l2, g2r2= lepton sneutrino chargino
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 sneutrino neutrino neutralino

       double precision function ampintchaHisquarksneutrino(g1, g1l2, 
     . g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amcha)
      implicit none
      double precision g1, g1l2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha,
     .g3l2, g3r1, g3r2, mfer, mv
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHisquarksneutrino=

     -  g1* (-2*g1L2*g2L1*g2L2*g3R1*g3R2*m3*m4*amcha*p1p2 + 
     -  2*g1R2*g2L1*g2R2*g3R1*g3R2*m1*m4*amcha*p2p3 - 
     -  2*g1R2*g2R1*g2R2*g3R1*g3R2*amcha*p1p4*p2p3 - 
     -  2*g1L2*g2L1*g2R2*g3R1*g3R2*m4*p1pI*p2p3 + 
     -  2*g1L2*g2L2*g2R1*g3R1*g3R2*m1*m3*amcha*p2p4 - 
     -  2*g1R2*g2R1*g2R2*g3R1*g3R2*amcha*p1p3*p2p4 - 
     -  2*g1R2*g2L2*g2R1*g3R1*g3R2*m3*p1pI*p2p4 + 
     -  2*g1R2*g2L1*g2L2*g3R1*g3R2*m1*m3*m4*p2pI + 
     -  2*g1L2*g2L1*g2R2*g3R1*g3R2*m4*p1p3*p2pI - 
     -  2*g1R2*g2L2*g2R1*g3R1*g3R2*m3*p1p4*p2pI + 
     -  2*g1R2*g2R1*g2R2*g3R1*g3R2*amcha*p1p2*p3p4 - 
     -  2*g1L2*g2R1*g2R2*g3R1*g3R2*m1*p2pI*p3p4 - 
     -  2*g1L2*g2L1*g2R2*g3R1*g3R2*m4*p1p2*p3pI + 
     -  2*g1L2*g2R1*g2R2*g3R1*g3R2*m1*p2p4*p3pI + 
     -  2*g1R2*g2L2*g2R1*g3R1*g3R2*m3*p1p2*p4pI + 
     -  2*g1L2*g2R1*g2R2*g3R1*g3R2*m1*p2p3*p4pI)
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference charged Higgs squark- chargino sup diagrams
c----- g1=charged Higgs squark squark coupling
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= neutralino quark squark
c---- g2l2, g2r2= sup down chargino
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 sup up neutralino

       double precision function ampintchaHisquarksup(g1, g1l2, 
     . g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amcha)
      implicit none
      double precision g1, g1l2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha,
     .g3l2, g3r1, g3r2, mfer, mv
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHisquarksup=
     -g1*4d0*((amcha*g1R2*g2R1*g2R2*g3L1*g3L2*m3*m4*p1p2)/2. + 
     -    (amcha*g1L2*g2L1*g2L2*g3R1*g3R2*m3*m4*p1p2)/2. - 
     -    (amcha*g1L2*g2L2*g2R1*g3L1*g3L2*m1*m4*p2p3)/2. - 
     -    (amcha*g1R2*g2L1*g2R2*g3R1*g3R2*m1*m4*p2p3)/2. - 
     -    (amcha*g1R2*g2L1*g2R2*g3L1*g3L2*m1*m3*p2p4)/2. - 
     -    (amcha*g1L2*g2L2*g2R1*g3R1*g3R2*m1*m3*p2p4)/2. - 
     -    (g1L2*g2R1*g2R2*g3L1*g3L2*m1*m3*m4*p2pI)/2. - 
     -    (g1R2*g2L1*g2L2*g3R1*g3R2*m1*m3*m4*p2pI)/2. + 
     -    (amcha*g1L2*g2L1*g2L2*g3L1*g3L2*
     -       (p1p4*p2p3 + p1p3*p2p4 - p1p2*p3p4))/2. + 
     -    (amcha*g1R2*g2R1*g2R2*g3R1*g3R2*
     -       (p1p4*p2p3 + p1p3*p2p4 - p1p2*p3p4))/2. + 
     -(g1R2*g2L2*g2R1*g3L1*g3L2*m4*(p1pI*p2p3 - p1p3*p2pI + p1p2*p3pI))/
     -     2. + (g1L2*g2L1*g2R2*g3R1*g3R2*m4*
     -       (p1pI*p2p3 - p1p3*p2pI + p1p2*p3pI))/2. + 
     -(g1L2*g2L1*g2R2*g3L1*g3L2*m3*(p1pI*p2p4 + p1p4*p2pI - p1p2*p4pI))/
     -     2. + (g1R2*g2L2*g2R1*g3R1*g3R2*m3*
     -       (p1pI*p2p4 + p1p4*p2pI - p1p2*p4pI))/2. - 
     -    (g1R2*g2L1*g2L2*g3L1*g3L2*m1*
     -       (-(p2pI*p3p4) + p2p4*p3pI + p2p3*p4pI))/2. - 
     -    (g1L2*g2R1*g2R2*g3R1*g3R2*m1*
     -       (-(p2pI*p3p4) + p2p4*p3pI + p2p3*p4pI))/2.)


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- squared matrix element of chargino charged Higgs diagram
c----- g1l1 g1r1 =chargino squark quark
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= chargino neutralino charged higgs
c---- g2l2, g2r2= chargino neutralino charged higgs
c---- g3r charged Higgs fermion fermion

       double precision function ampchaHichar(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pi, p3pi, p4pi, pi2,amcha1, amcha2)
      implicit none
      double precision g1l1, g1l2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pi, p2pi, p3pi, p4pi, pi2,p3p4, amcha1,amcha2,
     . g3r
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampchaHichar=
     -     -4*amcha1*amcha2*g1L2*g1R1*g2L2*g2R1*g3R**2*m1*m4*
     -   p2p3 - 4*amcha1*amcha2*g1L1*g1R2*g2L1*g2R2*g3R**2*m1*
     -   m4*p2p3 + 4*amcha1*amcha2*g1L1*g1L2*g2L1*g2L2*g3R**2*
     -   p1p4*p2p3 + 4*amcha1*amcha2*g1R1*g1R2*g2R1*g2R2*g3R**2*
     -   p1p4*p2p3 + 4*amcha2*g1L1*g1L2*g2L2*g2R1*g3R**2*m4*
     -   p1pI*p2p3 + 4*amcha1*g1R1*g1R2*g2L2*g2R1*g3R**2*m4*
     -   p1pI*p2p3 + 4*amcha1*g1L1*g1L2*g2L1*g2R2*g3R**2*m4*
     -   p1pI*p2p3 + 4*amcha2*g1R1*g1R2*g2L1*g2R2*g3R**2*m4*
     -   p1pI*p2p3 - 4*amcha2*g1L2*g1R1*g2L1*g2L2*g3R**2*m1*
     -   p2p3*p4pI - 4*amcha1*g1L1*g1R2*g2L1*g2L2*g3R**2*m1*
     -   p2p3*p4pI - 4*amcha1*g1L2*g1R1*g2R1*g2R2*g3R**2*m1*
     -   p2p3*p4pI - 4*amcha2*g1L1*g1R2*g2R1*g2R2*g3R**2*m1*
     -   p2p3*p4pI + 8*g1R1*g1R2*g2L1*g2L2*g3R**2*p1pI*p2p3*
     -   p4pI + 8*g1L1*g1L2*g2R1*g2R2*g3R**2*p1pI*p2p3*p4pI - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3R**2*m1*m4*p2p3*pI2 - 
     -  4*g1L2*g1R1*g2L1*g2R2*g3R**2*m1*m4*p2p3*pI2 - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3R**2*p1p4*p2p3*pI2 - 
     -  4*g1L1*g1L2*g2R1*g2R2*g3R**2*p1p4*p2p3*pI2
       return 
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- squared matrix element of SMfer charged Higgs diagram
c----- g1l1 g1r1 =neutralino squark quark
c---- g1l2, g1r2= neutralino squark quark
c---- g2l1, g2r1= fermion charged higgs
c---- g2l2, g2r2= fermion charged higgs
c---- g3r charged Higgs fermion fermion

       double precision function ampchaHiSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,mfer1, mfer2)
      implicit none
      double precision g1l1, g1l2, g2l1, g2l2, g3l1,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,p3p4,mfer1,mfer2,
     . g3r
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampchaHiSMfer=

     -   -4*g1L2*g1R1*g2L2*g2R1*g3R**2*m1*m4*mfer1*mfer2*p2p3 - 
     -  4*g1L1*g1R2*g2L1*g2R2*g3R**2*m1*m4*mfer1*mfer2*p2p3 + 
     -  4*g1L1*g1L2*g2L1*g2L2*g3R**2*mfer1*mfer2*p1p4*p2p3 + 
     -  4*g1R1*g1R2*g2R1*g2R2*g3R**2*mfer1*mfer2*p1p4*p2p3 - 
     -  4*g1R1*g1R2*g2L2*g2R1*g3R**2*m4*mfer1*p1pfer*p2p3 - 
     -  4*g1L1*g1L2*g2L1*g2R2*g3R**2*m4*mfer1*p1pfer*p2p3 - 
     -  4*g1L1*g1L2*g2L2*g2R1*g3R**2*m4*mfer2*p1pfer*p2p3 - 
     -  4*g1R1*g1R2*g2L1*g2R2*g3R**2*m4*mfer2*p1pfer*p2p3 + 
     -  4*g1L1*g1R2*g2L1*g2L2*g3R**2*m1*mfer1*p2p3*p4pfer + 
     -  4*g1L2*g1R1*g2R1*g2R2*g3R**2*m1*mfer1*p2p3*p4pfer + 
     -  4*g1L2*g1R1*g2L1*g2L2*g3R**2*m1*mfer2*p2p3*p4pfer + 
     -  4*g1L1*g1R2*g2R1*g2R2*g3R**2*m1*mfer2*p2p3*p4pfer + 
     -  8*g1R1*g1R2*g2L1*g2L2*g3R**2*p1pfer*p2p3*p4pfer + 
     -  8*g1L1*g1L2*g2R1*g2R2*g3R**2*p1pfer*p2p3*p4pfer - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3R**2*m1*m4*p2p3*pfer2 - 
     -  4*g1L2*g1R1*g2L1*g2R2*g3R**2*m1*m4*p2p3*pfer2 - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3R**2*p1p4*p2p3*pfer2 - 
     -  4*g1L1*g1L2*g2R1*g2R2*g3R**2*p1p4*p2p3*pfer2
       return
       end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference W chargino and charged Higgs chargino
c----- g1l1 g1r1 =chargino squark quark
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= chargino neutralino charged higgs
c---- g2l2, g2r2= chargino neutralino W
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 W fermion fermion

       double precision function ampintchaHiWchar(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amcha1, amcha2)
      implicit none
      double precision g1l1, g1l2, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pI, p2pI, p3pI, p4pI, pI2,p3p4,amcha1,amcha2,
     . g3r1
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2
     
      ampintchaHiWchar=

     - (4*amcha1*amcha2*g1R1*g1R2*g2R2*g2R1*g3L2*g3R1*m3*m4*
     -   p1p2 + 4*amcha1*amcha2*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*
     -   m3*m4*p1p2 - (4*amcha1*amcha2*g1R1*g1R2*g2R2*g2R1*g3L2*
     -     g3R1*m3*m4*p1p2*p2p3)/amw**2 - 
     -  (4*amcha1*amcha2*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*m4*
     -     p1p2*p2p3)/amw**2 - 
     -  (4*amcha1*amcha2*g1R1*g1R2*g2R2*g2R1*g3L2*g3R1*m3*m4*
     -     p1p3*p2p3)/amw**2 - 
     -  (4*amcha1*amcha2*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*m4*
     -     p1p3*p2p3)/amw**2 - 
     -  4*amcha1*amcha2*g1L1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*
     -   p2p4 - 4*amcha1*amcha2*g1L2*g1R1*g2R1*g2L2*g3L2*g3R1*
     -   m1*m3*p2p4 + 4*amcha1*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*
     -   p1pI*p2p4 + 4*amcha2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m3*
     -   p1pI*p2p4 + 4*amcha2*g1L1*g1L2*g2R1*g2L2*g3L2*g3R1*m3*
     -   p1pI*p2p4 + 4*amcha1*g1R1*g1R2*g2R1*g2L2*g3L2*g3R1*m3*
     -   p1pI*p2p4 + (4*amcha1*amcha2*g1L1*g1R2*g2L1*g2R2*g3L2*
     -     g3R1*m1*m3*p2p3*p2p4)/amw**2 + 
     -  (4*amcha1*amcha2*g1L2*g1R1*g2R1*g2L2*g3L2*g3R1*m1*m3*
     -     p2p3*p2p4)/amw**2 - 
     -  (4*amcha1*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*p1pI*p2p3*
     -     p2p4)/amw**2 - 
     -  (4*amcha2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m3*p1pI*p2p3*
     -     p2p4)/amw**2 - 
     -  (4*amcha2*g1L1*g1L2*g2R1*g2L2*g3L2*g3R1*m3*p1pI*p2p3*
     -     p2p4)/amw**2 - 
     -  (4*amcha1*g1R1*g1R2*g2R1*g2L2*g3L2*g3R1*m3*p1pI*p2p3*
     -     p2p4)/amw**2 - 
     -  4*amcha1*g1L2*g1R1*g2R2*g2R1*g3L2*g3R1*m1*m3*m4*p2pI - 
     -  4*amcha2*g1L1*g1R2*g2R2*g2R1*g3L2*g3R1*m1*m3*m4*p2pI - 
     -  4*amcha2*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*p2pI - 
     -  4*amcha1*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*p2pI + 
     -  4*amcha1*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*p1p4*p2pI - 
     -  4*amcha2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m3*p1p4*p2pI - 
     -  4*amcha2*g1L1*g1L2*g2R1*g2L2*g3L2*g3R1*m3*p1p4*p2pI + 
     -  4*amcha1*g1R1*g1R2*g2R1*g2L2*g3L2*g3R1*m3*p1p4*p2pI + 
     -  8*g1L1*g1L2*g2R2*g2R1*g3L2*g3R1*m3*m4*p1pI*p2pI + 
     -  8*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pI*p2pI + 
     -  (4*amcha1*g1L2*g1R1*g2R2*g2R1*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*amcha2*g1L1*g1R2*g2R2*g2R1*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*amcha2*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*amcha1*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p2pI)/amw**2 - 
     -  (4*amcha1*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*amcha2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*amcha2*g1L1*g1L2*g2R1*g2L2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p2pI)/amw**2 - 
     -  (4*amcha1*g1R1*g1R2*g2R1*g2L2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p2pI)/amw**2 - 
     -  (8*g1L1*g1L2*g2R2*g2R1*g3L2*g3R1*m3*m4*p1pI*p2p3*p2pI)/
     -   amw**2 - (8*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pI*
     -     p2p3*p2pI)/amw**2 + 
     -  (4*amcha1*amcha2*g1L1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*
     -     p2p3*p3p4)/amw**2 + 
     -  (4*amcha1*amcha2*g1L2*g1R1*g2R1*g2L2*g3L2*g3R1*m1*m3*
     -     p2p3*p3p4)/amw**2 - 
     -  (4*amcha1*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*p1pI*p2p3*
     -     p3p4)/amw**2 - 
     -  (4*amcha2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m3*p1pI*p2p3*
     -     p3p4)/amw**2 - 
     -  (4*amcha2*g1L1*g1L2*g2R1*g2L2*g3L2*g3R1*m3*p1pI*p2p3*
     -     p3p4)/amw**2 - 
     -  (4*amcha1*g1R1*g1R2*g2R1*g2L2*g3L2*g3R1*m3*p1pI*p2p3*
     -     p3p4)/amw**2 + 
     -  (4*amcha1*g1L2*g1R1*g2R2*g2R1*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p3pI)/amw**2 + 
     -  (4*amcha2*g1L1*g1R2*g2R2*g2R1*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p3pI)/amw**2 + 
     -  (4*amcha2*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p3pI)/amw**2 + 
     -  (4*amcha1*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p3pI)/amw**2 - 
     -  (4*amcha1*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p3pI)/amw**2 + 
     -  (4*amcha2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p3pI)/amw**2 + 
     -  (4*amcha2*g1L1*g1L2*g2R1*g2L2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p3pI)/amw**2 - 
     -  (4*amcha1*g1R1*g1R2*g2R1*g2L2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p3pI)/amw**2 - 
     -  (8*g1L1*g1L2*g2R2*g2R1*g3L2*g3R1*m3*m4*p1pI*p2p3*p3pI)/
     -   amw**2 - (8*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pI*
     -     p2p3*p3pI)/amw**2 - 
     -  4*amcha1*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*p1p2*p4pI + 
     -  4*amcha2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m3*p1p2*p4pI + 
     -  4*amcha2*g1L1*g1L2*g2R1*g2L2*g3L2*g3R1*m3*p1p2*p4pI - 
     -  4*amcha1*g1R1*g1R2*g2R1*g2L2*g3L2*g3R1*m3*p1p2*p4pI + 
     -  (4*amcha1*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*p1p2*p2p3*
     -     p4pI)/amw**2 - 
     -  (4*amcha2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m3*p1p2*p2p3*
     -     p4pI)/amw**2 - 
     -  (4*amcha2*g1L1*g1L2*g2R1*g2L2*g3L2*g3R1*m3*p1p2*p2p3*
     -     p4pI)/amw**2 + 
     -  (4*amcha1*g1R1*g1R2*g2R1*g2L2*g3L2*g3R1*m3*p1p2*p2p3*
     -     p4pI)/amw**2 + 
     -  (4*amcha1*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*p1p3*p2p3*
     -     p4pI)/amw**2 - 
     -  (4*amcha2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m3*p1p3*p2p3*
     -     p4pI)/amw**2 - 
     -  (4*amcha2*g1L1*g1L2*g2R1*g2L2*g3L2*g3R1*m3*p1p3*p2p3*
     -     p4pI)/amw**2 + 
     -  (4*amcha1*g1R1*g1R2*g2R1*g2L2*g3L2*g3R1*m3*p1p3*p2p3*
     -     p4pI)/amw**2 - 
     -  4*g1L1*g1L2*g2R2*g2R1*g3L2*g3R1*m3*m4*p1p2*pI2 - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1p2*pI2 + 
     -  (4*g1L1*g1L2*g2R2*g2R1*g3L2*g3R1*m3*m4*p1p2*p2p3*pI2)/
     -   amw**2 + (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1p2*
     -     p2p3*pI2)/amw**2 + 
     -  (4*g1L1*g1L2*g2R2*g2R1*g3L2*g3R1*m3*m4*p1p3*p2p3*pI2)/
     -   amw**2 + (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1p3*
     -     p2p3*pI2)/amw**2 - 
     -  4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p4*pI2 - 
     -  4*g1L1*g1R2*g2R1*g2L2*g3L2*g3R1*m1*m3*p2p4*pI2 + 
     -  (4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p3*p2p4*pI2)/
     -   amw**2 + (4*g1L1*g1R2*g2R1*g2L2*g3L2*g3R1*m1*m3*p2p3*
     -     p2p4*pI2)/amw**2 + 
     -  (4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p3*p3p4*pI2)/
     -   amw**2 + (4*g1L1*g1R2*g2R1*g2L2*g3L2*g3R1*m1*m3*p2p3*
     -     p3p4*pI2)/amw**2)
        
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference W SM fermion and charged Higgs chargino
c----- g1l1 g1r1 =chargino squark quark
c---- g1l2, g1r2= neutralino squark quark
c---- g2l1, g2r1= chargino neutralino charged higgs
c---- g2l2, g2r2= Sm fermion W
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 W fermion fermion

       double precision function ampintchaHiWSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, pIpFer, amfer, amcha)
      implicit none
      double precision g1l1, g1l2, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,p1pi, p2pI, p3pI, 
     .p4pI, pI2,p3p4,amcha,amfer,pipfer,
     . g3r1
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHiWSMfer=
     --(4*amcha*amfer*g1L2*g1R1*g2L2*g2R1*g3L2*g3R1*m3*m4*
     -   p1p2 - (4*amcha*amfer*g1L2*g1R1*g2L2*g2R1*g3L2*g3R1*m3*
     -     m4*p1p2*p2p3)/amw**2 - 
     -  (4*amcha*amfer*g1L2*g1R1*g2L2*g2R1*g3L2*g3R1*m3*m4*p1p3*
     -     p2p3)/amw**2 - 
     -  4*amcha*amfer*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m1*m3*
     -   p2p4 + 4*amcha*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*p1pfer*
     -   p2p4 + 4*amfer*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*p1pI*
     -   p2p4 + (4*amcha*amfer*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m1*
     -     m3*p2p3*p2p4)/amw**2 - 
     -  (4*amcha*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*p1pfer*p2p3*
     -     p2p4)/amw**2 - 
     -  (4*amfer*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*p1pI*p2p3*
     -     p2p4)/amw**2 + 
     -  4*amcha*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*p2pfer - 
     -  4*amcha*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*p1p4*p2pfer - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pI*p2pfer - 
     -  (4*amcha*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p2pfer)/amw**2 + 
     -  (4*amcha*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*p1p4*p2p3*
     -     p2pfer)/amw**2 + 
     -  (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pI*p2p3*
     -     p2pfer)/amw**2 - 
     -  4*amfer*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*p2pI - 
     -  4*amfer*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*p1p4*p2pI + 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pfer*p2pI + 
     -  (4*amfer*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*amfer*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p2pI)/amw**2 - 
     -  (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pfer*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*amcha*amfer*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m1*m3*p2p3*
     -     p3p4)/amw**2 - 
     -  (4*amcha*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*p1pfer*p2p3*
     -     p3p4)/amw**2 - 
     -  (4*amfer*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*p1pI*p2p3*
     -     p3p4)/amw**2 - 
     -  (4*amcha*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p3pfer)/amw**2 + 
     -  (4*amcha*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*p1p4*p2p3*
     -     p3pfer)/amw**2 + 
     -  (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pI*p2p3*
     -     p3pfer)/amw**2 + 
     -  (4*amfer*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*p2p3*
     -     p3pI)/amw**2 + 
     -  (4*amfer*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*p1p4*p2p3*
     -     p3pI)/amw**2 - 
     -  (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pfer*p2p3*
     -     p3pI)/amw**2 - 
     -  4*amcha*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*p1p2*p4pfer + 
     -  (4*amcha*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*p1p2*p2p3*
     -     p4pfer)/amw**2 + 
     -  (4*amcha*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*p1p3*p2p3*
     -     p4pfer)/amw**2 + 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2pI*p4pfer - 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p2pI*
     -     p4pfer)/amw**2 - 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p3pI*
     -     p4pfer)/amw**2 + 
     -  4*amfer*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*p1p2*p4pI - 
     -  (4*amfer*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*p1p2*p2p3*
     -     p4pI)/amw**2 - 
     -  (4*amfer*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*p1p3*p2p3*
     -     p4pI)/amw**2 + 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2pfer*p4pI - 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p2pfer*
     -     p4pI)/amw**2 - 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p3pfer*
     -     p4pI)/amw**2 - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1p2*pIpfer + 
     -  (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1p2*p2p3*
     -     pIpfer)/amw**2 + 
     -  (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m3*m4*p1p3*p2p3*
     -     pIpfer)/amw**2 - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p4*pIpfer + 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p2p4*
     -     pIpfer)/amw**2 + 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p3p4*
     -     pIpfer)/amw**2)
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference slepton chargino and charged Higgs chargino
c----- g1l1 g1r1 =chargino squark quark
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= chargino neutralino charged higgs
c---- g2r2= chargino slepton neutrino
c---- g3r1 charged Higgs fermion fermion
c---- g3l2, g3r2 neutralion slepton lepton

       double precision function ampintchaHichargslepton(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amcha1, amcha2)
      implicit none
      double precision g1l1, g1l2, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pI, p2pI, p3pI, p4pI, pI2,p3p4,amcha1,amcha2,
     . g3r1, g3r2
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

       
      ampintchaHichargslepton=
     --4d0*(-(amcha1*amcha2*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1p2)/2.
     -- (amcha1*amcha2*g1L1*g1R2*g2L1*g2R2*g3R1*g3R2*m1*m4*p2p3)/2. + 
     -    (amcha1*amcha2*g1L1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p4)/2. + 
     -    (amcha1*g1L2*g1R1*g2R1*g2R2*g3L2*g3R1*m1*m3*m4*p2pI)/2. + 
     -    (amcha2*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m1*m3*m4*p2pI)/2. - 
     -    g1L1*g1L2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1pI*p2pI + 
     -    (amcha1*amcha2*g1R1*g1R2*g2R1*g2R2*g3R1*g3R2*
     -       (p1p4*p2p3 - p1p3*p2p4 + p1p2*p3p4))/2. + 
     -    (amcha1*g1L1*g1L2*g2L1*g2R2*g3R1*g3R2*m4*
     -       (p1pI*p2p3 + p1p3*p2pI - p1p2*p3pI))/2. + 
     -    (amcha2*g1R1*g1R2*g2L1*g2R2*g3R1*g3R2*m4*
     -       (p1pI*p2p3 - p1p3*p2pI + p1p2*p3pI))/2. + 
     -g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*(p1pI*p2p3 + p1p3*p2pI- p1p2*p3pI)*
     -     p4pI - (amcha1*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*
     -       (p1pI*p2p4 + p1p4*p2pI - p1p2*p4pI))/2. - 
     -    g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p3pI*
     -     (p1pI*p2p4 + p1p4*p2pI - p1p2*p4pI) - 
     -    (amcha2*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m3*
     -       (p1pI*p2p4 - p1p4*p2pI + p1p2*p4pI))/2. + 
     -    g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p2pI*
     -     (p1pI*p3p4 + p1p4*p3pI - p1p3*p4pI) - 
     -    (amcha1*g1L2*g1R1*g2R1*g2R2*g3R1*g3R2*m1*
     -       (p2pI*p3p4 - p2p4*p3pI + p2p3*p4pI))/2. - 
     -    (amcha2*g1L1*g1R2*g2R1*g2R2*g3R1*g3R2*m1*
     -       (p2pI*p3p4 - p2p4*p3pI + p2p3*p4pI))/2. + 
     -    (g1L1*g1L2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1p2*pI2)/2. - 
     -    (g1L2*g1R1*g2L1*g2R2*g3R1*g3R2*m1*m4*p2p3*pI2)/2. + 
     -    (g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p4*pI2)/2. - 
     -    (g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*
     - (p1p4*p2p3 - p1p3*p2p4 + p1p2*p3p4)*pI2)/2.)


       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference sneutrino chargino and charged Higgs chargino
c----- g1l1 g1r1 =chargino squark quark
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= chargino neutralino charged higgs
c---- g2l2, g2r2= chargino lepton sneutrino
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 neutralion sneutrino neutrino

       double precision function ampintchaHichargsneutrino(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amcha1, amcha2)
      implicit none
      double precision g1l1, g1l2, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pI, p2pI, p3pI, p4pI, pI2,p3p4,amcha1,amcha2,
     . g3r1, g3r2
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHichargsneutrino=
     -4d0*((amcha1*amcha2*g1L1*g1L2*g2L1*g2L2*g3R1*g3R2*m3*m4*p1p2)/2.- 
     -    (amcha1*amcha2*g1L1*g1R2*g2L1*g2R2*g3R1*g3R2*m1*m4*p2p3)/2. - 
     -    (amcha1*amcha2*g1L2*g1R1*g2L2*g2R1*g3R1*g3R2*m1*m3*p2p4)/2. - 
     -    (amcha2*g1L2*g1R1*g2L1*g2L2*g3R1*g3R2*m1*m3*m4*p2pI)/2. - 
     -    (amcha1*g1L1*g1R2*g2L1*g2L2*g3R1*g3R2*m1*m3*m4*p2pI)/2. + 
     -    g1R1*g1R2*g2L1*g2L2*g3R1*g3R2*m3*m4*p1pI*p2pI + 
     -    (amcha1*amcha2*g1R1*g1R2*g2R1*g2R2*g3R1*g3R2*
     -       (p1p4*p2p3 + p1p3*p2p4 - p1p2*p3p4))/2. + 
     -    (amcha2*g1R1*g1R2*g2L1*g2R2*g3R1*g3R2*m4*
     -       (p1pI*p2p3 + p1p3*p2pI - p1p2*p3pI))/2. + 
     -    (amcha1*g1L1*g1L2*g2L1*g2R2*g3R1*g3R2*m4*
     -       (p1pI*p2p3 - p1p3*p2pI + p1p2*p3pI))/2. + 
     -g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*(p1pI*p2p3 - p1p3*p2pI + p1p2*p3pI)*
     -     p4pI + (amcha1*g1R1*g1R2*g2L2*g2R1*g3R1*g3R2*m3*
     -       (p1pI*p2p4 + p1p4*p2pI - p1p2*p4pI))/2. + 
     -    g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p3pI*
     -     (p1pI*p2p4 + p1p4*p2pI - p1p2*p4pI) + 
     -    (amcha2*g1L1*g1L2*g2L2*g2R1*g3R1*g3R2*m3*
     -       (p1pI*p2p4 - p1p4*p2pI + p1p2*p4pI))/2. - 
     -    g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*p2pI*
     -     (p1pI*p3p4 + p1p4*p3pI - p1p3*p4pI) - 
     -    (amcha1*g1L2*g1R1*g2R1*g2R2*g3R1*g3R2*m1*
     -       (-(p2pI*p3p4) + p2p4*p3pI + p2p3*p4pI))/2. - 
     -    (amcha2*g1L1*g1R2*g2R1*g2R2*g3R1*g3R2*m1*
     -       (-(p2pI*p3p4) + p2p4*p3pI + p2p3*p4pI))/2. - 
     -    (g1R1*g1R2*g2L1*g2L2*g3R1*g3R2*m3*m4*p1p2*pI2)/2. - 
     -    (g1L2*g1R1*g2L1*g2R2*g3R1*g3R2*m1*m4*p2p3*pI2)/2. - 
     -    (g1L1*g1R2*g2L2*g2R1*g3R1*g3R2*m1*m3*p2p4*pI2)/2. - 
     -    (g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*
     -       (p1p4*p2p3 + p1p3*p2p4 - p1p2*p3p4)*pI2)/2.)





        return
        end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference W squuark and charged Higgs chargino
c----- g1l1 g1r1 =chargino squark quark
c---- g1= W squark squark
c---- g2l1, g2r1= chargino neutralino charged higgs
c---- g2l2, g2r2= neutralino squark quark
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 W fermion fermion

       double precision function ampintchaHichargWsquark(g1,
     .  g1l1,  g1r1,g2l1, g2r1, g2l2, g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amcha)
      implicit none
      double precision g1, g1l1, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pI, p2pI, p3pI, p4pI, pI2,p3p4,amcha,
     . g3r1
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2
      

      
      ampintchaHichargWsquark=
     - 16*g1*((g3L2*g3R1*m3*(2*p1p2 + p2p3 + 2*p2p4))/2d0 - 
     - (g3L2*g3R1*m3*p2p3*(m3**2 + 2*p1p2 + 2*p1p3 + 2*p2p3 + 2*p2p4 + 
     -         2*p3p4))/(2.*amw**2))*
     -  (-(amcha*g1R1*g2L2*g2R1*m1*m4)/2d0 - 
     -    (amcha*g1L1*g2L1*g2R2*m1*m4)/2d0 + 
     -(amcha*g1L1*g2L1*g2L2*p1p4)/2d0 +(amcha*g1R1*g2R1*g2R2*p1p4)/2d0+ 
     -    (g1L1*g2L2*g2R1*m4*p1pI)/2d0 + (g1R1*g2L1*g2R2*m4*p1pI)/2d0 - 
     -    (g1R1*g2L1*g2L2*m1*p4pI)/2d0 - (g1L1*g2R1*g2R2*m1*p4pI)/2d0)


       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference charged Higgs squuark and charged Higgs chargino
c----- g1l1 g1r1 =chargino squark quark
c---- g12= charged Higgs squark squark
c---- g2l1, g2r1= chargino neutralino charged higgs
c---- g2l2, g2r2= neutralino squark quark
c---- g3r1 charged Higgs fermion fermion
c---- g322 charged Higgs fermion fermion

       double precision function ampintchaHichargchaHisquark(g12,
     .  g1l1,  g1r1,g2l1, g2r1, g2l2, g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pi, p2pI, p3pI, p4pI, pI2,amcha)
      implicit none
      double precision g12, g1l1, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pI, p2pI, p3pI, p4pI, pI2,p3p4,amcha,
     . g3r1, g3r2
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHichargchaHisquark=


     - g12*(-4*g1r1*g2l2*g2R1*g3R1*g3R2*m1*m4*amcha*p2p3 - 
     -  4*g1l1*g2L1*g2r2*g3R1*g3R2*m1*m4*amcha*p2p3 + 
     -  4*g1l1*g2L1*g2l2*g3R1*g3R2*amcha*p1p4*p2p3 + 
     -  4*g1r1*g2R1*g2r2*g3R1*g3R2*amcha*p1p4*p2p3 + 
     -  4*g1l1*g2l2*g2R1*g3R1*g3R2*m4*p1pI*p2p3 + 
     -  4*g1r1*g2L1*g2r2*g3R1*g3R2*m4*p1pI*p2p3 - 
     -  4*g1r1*g2L1*g2l2*g3R1*g3R2*m1*p2p3*p4pI - 
     -  4*g1l1*g2R1*g2r2*g3R1*g3R2*m1*p2p3*p4pI)

       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference W SM fermion and charged Higgs Sm fermion
c----- g1l1 g1r1 =neutralino squark quark
c---- g1l2, g1r2= neutralino squark quark
c---- g2l1, g2r1= quarks charged higgs
c---- g2l2, g2r2= Sm fermion W
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 W fermion fermion

       double precision function ampintchaHiSmferWSMfer(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, mfer1, mfer2)
      implicit none
      double precision g1l1, g1l2, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,p1pi, p2pI, p3pI, 
     .p4pI, pI2,p3p4,mfer1, mfer2, 
     . g3r1
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHiSmferWSMfer=


     - -4*g1L2*g1R1*g2L2*g2R1*g3L2*g3R1*m3*m4*mfer1*mfer2*
     -   p1p2 + (4*g1L2*g1R1*g2L2*g2R1*g3L2*g3R1*m3*m4*mfer1*
     -     mfer2*p1p2*p2p3)/amw**2 + 
     -  (4*g1L2*g1R1*g2L2*g2R1*g3L2*g3R1*m3*m4*mfer1*mfer2*p1p3*
     -     p2p3)/amw**2 + 
     -  4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m1*m3*mfer1*mfer2*
     -   p2p4 - 4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*mfer1*p1pfer*
     -   p2p4 + 4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m3*mfer2*p1pfer*
     -   p2p4 - (4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m1*m3*mfer1*
     -     mfer2*p2p3*p2p4)/amw**2 + 
     -  (4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*mfer1*p1pfer*p2p3*
     -     p2p4)/amw**2 - 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m3*mfer2*p1pfer*p2p3*
     -     p2p4)/amw**2 - 
     -  4*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*mfer1*p2pfer - 
     -  4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*mfer2*p2pfer + 
     -  4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*mfer1*p1p4*p2pfer - 
     -  4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m3*mfer2*p1p4*p2pfer + 
     -  (4*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*mfer1*p2p3*
     -     p2pfer)/amw**2 + 
     -  (4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*mfer2*p2p3*
     -     p2pfer)/amw**2 - 
     -  (4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*mfer1*p1p4*p2p3*
     -     p2pfer)/amw**2 + 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m3*mfer2*p1p4*p2p3*
     -     p2pfer)/amw**2 - 
     -  (4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m1*m3*mfer1*mfer2*p2p3*
     -     p3p4)/amw**2 + 
     -  (4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*mfer1*p1pfer*p2p3*
     -     p3p4)/amw**2 - 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m3*mfer2*p1pfer*p2p3*
     -     p3p4)/amw**2 + 
     -  (4*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*mfer1*p2p3*
     -     p3pfer)/amw**2 + 
     -  (4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*mfer2*p2p3*
     -     p3pfer)/amw**2 - 
     -  (4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*mfer1*p1p4*p2p3*
     -     p3pfer)/amw**2 + 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m3*mfer2*p1p4*p2p3*
     -     p3pfer)/amw**2 + 
     -  4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*mfer1*p1p2*p4pfer + 
     -  4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m3*mfer2*p1p2*p4pfer - 
     -  (4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*mfer1*p1p2*p2p3*
     -     p4pfer)/amw**2 - 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m3*mfer2*p1p2*p2p3*
     -     p4pfer)/amw**2 - 
     -  (4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*mfer1*p1p3*p2p3*
     -     p4pfer)/amw**2 - 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m3*mfer2*p1p3*p2p3*
     -     p4pfer)/amw**2 + 
     -  8*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*p2pfer*p4pfer - 
     -  (8*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*p2p3*p2pfer*
     -     p4pfer)/amw**2 - 
     -  (8*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*p2p3*p3pfer*
     -     p4pfer)/amw**2 - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*m4*p1p2*pfer2 + 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*m4*p1p2*p2p3*pfer2)/
     -   amw**2 + (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*m4*p1p3*
     -     p2p3*pfer2)/amw**2 - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*p2p4*pfer2 + 
     -  (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*p2p3*p2p4*pfer2)/
     -   amw**2 + (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*p2p3*
     -     p3p4*pfer2)/amw**2
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference W squark and charged Higgs Sm fermion
c----- g1l1 g1r1 =neutralino squark quark
c---- g1l2, g1r2= neutralino squark quark
c---- g2l1, g2r1= quarks charged higgs
c---- g2l2, g2r2= Sm fermion W
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 W fermion fermion

       double precision function ampintchaHiSmferWsquark(g1l1,
     . g1r1, g11, g2l1, g2r1, g2l2,g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, 
     . p3pfer, p4pfer, pfer2, mfer)

      implicit none
      double precision g1l1, g11, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,p1pi, p2pI, p3pI, 
     .p4pI, pI2,p3p4,mfer,  
     . g3r1
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHiSmferWsquark=

     - 8*g11*g1r1*g2l2*g2r1*g3l2*g3R1*m1*m3*m4*mfer*p1p2 + 
     -  8*g11*g1l1*g2l1*g2r2*g3l2*g3R1*m1*m3*m4*mfer*p1p2 - 
     -  8*g11*g1l1*g2l1*g2l2*g3l2*g3R1*m3*mfer*p1p2*p1p4 - 
     -  8*g11*g1r1*g2r1*g2r2*g3l2*g3R1*m3*mfer*p1p2*p1p4 + 
     -  8*g11*g1r1*g2l1*g2l2*g3l2*g3R1*m3*m4*p1p2*p1pfer + 
     -  8*g11*g1l1*g2r1*g2r2*g3l2*g3R1*m3*m4*p1p2*p1pfer + 
     -  4*g11*g1r1*g2l2*g2r1*g3l2*g3R1*m1*m3*m4*mfer*p2p3 + 
     -  4*g11*g1l1*g2l1*g2r2*g3l2*g3R1*m1*m3*m4*mfer*p2p3 - 
     -  (4*g11*g1r1*g2l2*g2r1*g3l2*g3R1*m1*m3**3*m4*mfer*p2p3)/amw**2 - 
     -  (4*g11*g1l1*g2l1*g2r2*g3l2*g3R1*m1*m3**3*m4*mfer*p2p3)/amw**2 - 
     -(8*g11*g1r1*g2l2*g2r1*g3l2*g3R1*m1*m3*m4*mfer*p1p2*p2p3)/amw**2 - 
     -(8*g11*g1l1*g2l1*g2r2*g3l2*g3R1*m1*m3*m4*mfer*p1p2*p2p3)/amw**2 - 
     -(8*g11*g1r1*g2l2*g2r1*g3l2*g3R1*m1*m3*m4*mfer*p1p3*p2p3)/amw**2 - 
     -(8*g11*g1l1*g2l1*g2r2*g3l2*g3R1*m1*m3*m4*mfer*p1p3*p2p3)/amw**2 - 
     -  4*g11*g1l1*g2l1*g2l2*g3l2*g3R1*m3*mfer*p1p4*p2p3 - 
     -  4*g11*g1r1*g2r1*g2r2*g3l2*g3R1*m3*mfer*p1p4*p2p3 + 
     -  (4*g11*g1l1*g2l1*g2l2*g3l2*g3R1*m3**3*mfer*p1p4*p2p3)/amw**2 + 
     -  (4*g11*g1r1*g2r1*g2r2*g3l2*g3R1*m3**3*mfer*p1p4*p2p3)/amw**2 + 
     -(8*g11*g1l1*g2l1*g2l2*g3l2*g3R1*m3*mfer*p1p2*p1p4*p2p3)/amw**2 + 
     -(8*g11*g1r1*g2r1*g2r2*g3l2*g3R1*m3*mfer*p1p2*p1p4*p2p3)/amw**2 + 
     -(8*g11*g1l1*g2l1*g2l2*g3l2*g3R1*m3*mfer*p1p3*p1p4*p2p3)/amw**2 + 
     -(8*g11*g1r1*g2r1*g2r2*g3l2*g3R1*m3*mfer*p1p3*p1p4*p2p3)/amw**2 + 
     -  4*g11*g1r1*g2l1*g2l2*g3l2*g3R1*m3*m4*p1pfer*p2p3 + 
     -  4*g11*g1l1*g2r1*g2r2*g3l2*g3R1*m3*m4*p1pfer*p2p3 - 
     -  (4*g11*g1r1*g2l1*g2l2*g3l2*g3R1*m3**3*m4*p1pfer*p2p3)/amw**2 - 
     -  (4*g11*g1l1*g2r1*g2r2*g3l2*g3R1*m3**3*m4*p1pfer*p2p3)/amw**2 - 
     -(8*g11*g1r1*g2l1*g2l2*g3l2*g3R1*m3*m4*p1p2*p1pfer*p2p3)/amw**2 - 
     -(8*g11*g1l1*g2r1*g2r2*g3l2*g3R1*m3*m4*p1p2*p1pfer*p2p3)/amw**2 - 
     -(8*g11*g1r1*g2l1*g2l2*g3l2*g3R1*m3*m4*p1p3*p1pfer*p2p3)/amw**2 - 
     -(8*g11*g1l1*g2r1*g2r2*g3l2*g3R1*m3*m4*p1p3*p1pfer*p2p3)/amw**2 - 
     -(8*g11*g1r1*g2l2*g2r1*g3l2*g3R1*m1*m3*m4*mfer*p2p3**2)/amw**2 - 
     -(8*g11*g1l1*g2l1*g2r2*g3l2*g3R1*m1*m3*m4*mfer*p2p3**2)/amw**2 + 
     -(8*g11*g1l1*g2l1*g2l2*g3l2*g3R1*m3*mfer*p1p4*p2p3**2)/amw**2 + 
     -(8*g11*g1r1*g2r1*g2r2*g3l2*g3R1*m3*mfer*p1p4*p2p3**2)/amw**2 - 
     -(8*g11*g1r1*g2l1*g2l2*g3l2*g3R1*m3*m4*p1pfer*p2p3**2)/amw**2 - 
     -(8*g11*g1l1*g2r1*g2r2*g3l2*g3R1*m3*m4*p1pfer*p2p3**2)/amw**2 + 
     -8*g11*g1r1*g2l2*g2r1*g3l2*g3R1*m1*m3*m4*mfer*p2p4 + 
     -8*g11*g1l1*g2l1*g2r2*g3l2*g3R1*m1*m3*m4*mfer*p2p4 - 
     -  8*g11*g1l1*g2l1*g2l2*g3l2*g3R1*m3*mfer*p1p4*p2p4 - 
     -  8*g11*g1r1*g2r1*g2r2*g3l2*g3R1*m3*mfer*p1p4*p2p4 + 
     -  8*g11*g1r1*g2l1*g2l2*g3l2*g3R1*m3*m4*p1pfer*p2p4 + 
     -  8*g11*g1l1*g2r1*g2r2*g3l2*g3R1*m3*m4*p1pfer*p2p4 - 
     -(8*g11*g1r1*g2l2*g2r1*g3l2*g3R1*m1*m3*m4*mfer*p2p3*p2p4)/amw**2 - 
     -(8*g11*g1l1*g2l1*g2r2*g3l2*g3R1*m1*m3*m4*mfer*p2p3*p2p4)/amw**2 + 
     -(8*g11*g1l1*g2l1*g2l2*g3l2*g3R1*m3*mfer*p1p4*p2p3*p2p4)/amw**2 + 
     -(8*g11*g1r1*g2r1*g2r2*g3l2*g3R1*m3*mfer*p1p4*p2p3*p2p4)/amw**2 - 
     -(8*g11*g1r1*g2l1*g2l2*g3l2*g3R1*m3*m4*p1pfer*p2p3*p2p4)/amw**2 - 
     -(8*g11*g1l1*g2r1*g2r2*g3l2*g3R1*m3*m4*p1pfer*p2p3*p2p4)/amw**2 - 
     -(8*g11*g1r1*g2l2*g2r1*g3l2*g3R1*m1*m3*m4*mfer*p2p3*p3p4)/amw**2 - 
     -(8*g11*g1l1*g2l1*g2r2*g3l2*g3R1*m1*m3*m4*mfer*p2p3*p3p4)/amw**2 + 
     -(8*g11*g1l1*g2l1*g2l2*g3l2*g3R1*m3*mfer*p1p4*p2p3*p3p4)/amw**2 + 
     -(8*g11*g1r1*g2r1*g2r2*g3l2*g3R1*m3*mfer*p1p4*p2p3*p3p4)/amw**2 - 
     -(8*g11*g1r1*g2l1*g2l2*g3l2*g3R1*m3*m4*p1pfer*p2p3*p3p4)/amw**2 - 
     -(8*g11*g1l1*g2r1*g2r2*g3l2*g3R1*m3*m4*p1pfer*p2p3*p3p4)/amw**2 - 
     -8*g11*g1l1*g2l2*g2r1*g3l2*g3R1*m1*m3*p1p2*p4pfer - 
     -8*g11*g1r1*g2l1*g2r2*g3l2*g3R1*m1*m3*p1p2*p4pfer - 
     -4*g11*g1l1*g2l2*g2r1*g3l2*g3R1*m1*m3*p2p3*p4pfer - 
     -4*g11*g1r1*g2l1*g2r2*g3l2*g3R1*m1*m3*p2p3*p4pfer + 
     -(4*g11*g1l1*g2l2*g2r1*g3l2*g3R1*m1*m3**3*p2p3*p4pfer)/amw**2 + 
     -(4*g11*g1r1*g2l1*g2r2*g3l2*g3R1*m1*m3**3*p2p3*p4pfer)/amw**2 + 
     -(8*g11*g1l1*g2l2*g2r1*g3l2*g3R1*m1*m3*p1p2*p2p3*p4pfer)/amw**2 + 
     -(8*g11*g1r1*g2l1*g2r2*g3l2*g3R1*m1*m3*p1p2*p2p3*p4pfer)/amw**2 + 
     -(8*g11*g1l1*g2l2*g2r1*g3l2*g3R1*m1*m3*p1p3*p2p3*p4pfer)/amw**2 + 
     -(8*g11*g1r1*g2l1*g2r2*g3l2*g3R1*m1*m3*p1p3*p2p3*p4pfer)/amw**2 + 
     -(8*g11*g1l1*g2l2*g2r1*g3l2*g3R1*m1*m3*p2p3**2*p4pfer)/amw**2 + 
     - (8*g11*g1r1*g2l1*g2r2*g3l2*g3R1*m1*m3*p2p3**2*p4pfer)/amw**2 - 
     -  8*g11*g1l1*g2l2*g2r1*g3l2*g3R1*m1*m3*p2p4*p4pfer - 
     -  8*g11*g1r1*g2l1*g2r2*g3l2*g3R1*m1*m3*p2p4*p4pfer + 
     -(8*g11*g1l1*g2l2*g2r1*g3l2*g3R1*m1*m3*p2p3*p2p4*p4pfer)/amw**2 + 
     -(8*g11*g1r1*g2l1*g2r2*g3l2*g3R1*m1*m3*p2p3*p2p4*p4pfer)/amw**2 + 
     -(8*g11*g1l1*g2l2*g2r1*g3l2*g3R1*m1*m3*p2p3*p3p4*p4pfer)/amw**2 + 
     -(8*g11*g1r1*g2l1*g2r2*g3l2*g3R1*m1*m3*p2p3*p3p4*p4pfer)/amw**2



      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference W chargino and charged Higgs Sm fermion
c----- g1l1 g1r1 =neutralino squark quark
c---- g1l2, g1r2= chargino squark quark
c---- g2l1, g2r1= quarks charged higgs
c---- g2l2, g2r2= neutralino chargino W
c---- g3r1 charged Higgs fermion fermion
c---- g3l2 W fermion fermion

       double precision function ampintchaHiSmferWcharg(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3l2,
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2, 
     . p1pi, p2pi, p3pi, p4pi, pi2, pipfer, amcha, mfer)
      implicit none
      double precision g1l1, g1l2, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,p1pi, p2pI, p3pI, 
     .p4pI, pI2,p3p4,mfer1, mfer2, pipfer,amcha, mfer,
     . g3r1
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHiSmferWcharg=

     -   -4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*m4*mfer1*amcha*p1p2 - 
     -  4*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*m4*mfer1*amcha*p1p2 + 
     -  (4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*m4*mfer1*amcha*p1p2*
     -     p2p3)/amw**2 + 
     -  (4*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*m4*mfer1*amcha*p1p2*
     -     p2p3)/amw**2 + 
     -  (4*g1R1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*m4*mfer1*amcha*p1p3*
     -     p2p3)/amw**2 + 
     -  (4*g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*m4*mfer1*amcha*p1p3*
     -     p2p3)/amw**2 + 
     -  4*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*mfer1*amcha*p2p4 + 
     -  4*g1L2*g1R1*g2R1*g2R2*g3L2*g3R1*m1*m3*mfer1*amcha*p2p4 + 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*amcha*p1pfer*p2p4 + 
     -  4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m3*amcha*p1pfer*p2p4 - 
     -  4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*mfer1*p1pI*p2p4 - 
     -  4*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*mfer1*p1pI*p2p4 - 
     -  (4*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*mfer1*amcha*p2p3*
     -     p2p4)/amw**2 - 
     -  (4*g1L2*g1R1*g2R1*g2R2*g3L2*g3R1*m1*m3*mfer1*amcha*p2p3*
     -     p2p4)/amw**2 - 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*amcha*p1pfer*p2p3*
     -     p2p4)/amw**2 - 
     -  (4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m3*amcha*p1pfer*p2p3*
     -     p2p4)/amw**2 + 
     -  (4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*mfer1*p1pI*p2p3*
     -     p2p4)/amw**2 + 
     -  (4*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*mfer1*p1pI*p2p3*
     -     p2p4)/amw**2 - 
     -  4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*amcha*p2pfer - 
     -  4*g1L1*g1L2*g2R1*g2R2*g3L2*g3R1*m1*m3*m4*amcha*p2pfer - 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*amcha*p1p4*p2pfer - 
     -  4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m3*amcha*p1p4*p2pfer + 
     -  4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pI*p2pfer + 
     -  4*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1pI*p2pfer + 
     -  (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*amcha*p2p3*
     -     p2pfer)/amw**2 + 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L2*g3R1*m1*m3*m4*amcha*p2p3*
     -     p2pfer)/amw**2 + 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*amcha*p1p4*p2p3*
     -     p2pfer)/amw**2 + 
     -  (4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m3*amcha*p1p4*p2p3*
     -     p2pfer)/amw**2 - 
     -  (4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pI*p2p3*
     -     p2pfer)/amw**2 - 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1pI*p2p3*
     -     p2pfer)/amw**2 + 
     -  4*g1L2*g1R1*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*mfer1*p2pI + 
     -  4*g1L1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*m4*mfer1*p2pI - 
     -  4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*mfer1*p1p4*p2pI - 
     -  4*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*mfer1*p1p4*p2pI + 
     -  4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pfer*p2pI + 
     -  4*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1pfer*p2pI - 
     -  (4*g1L2*g1R1*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*mfer1*p2p3*
     -     p2pI)/amw**2 - 
     -  (4*g1L1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*m4*mfer1*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*mfer1*p1p4*p2p3*
     -     p2pI)/amw**2 + 
     -  (4*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*mfer1*p1p4*p2p3*
     -     p2pI)/amw**2 - 
     -  (4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pfer*p2p3*
     -     p2pI)/amw**2 - 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1pfer*p2p3*
     -     p2pI)/amw**2 - 
     -  (4*g1L1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*mfer1*amcha*p2p3*
     -     p3p4)/amw**2 - 
     -  (4*g1L2*g1R1*g2R1*g2R2*g3L2*g3R1*m1*m3*mfer1*amcha*p2p3*
     -     p3p4)/amw**2 - 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*amcha*p1pfer*p2p3*
     -     p3p4)/amw**2 - 
     -  (4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m3*amcha*p1pfer*p2p3*
     -     p3p4)/amw**2 + 
     -  (4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*mfer1*p1pI*p2p3*
     -     p3p4)/amw**2 + 
     -  (4*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*mfer1*p1pI*p2p3*
     -     p3p4)/amw**2 + 
     -  (4*g1R1*g1R2*g2L1*g2L2*g3L2*g3R1*m1*m3*m4*amcha*p2p3*
     -     p3pfer)/amw**2 + 
     -  (4*g1L1*g1L2*g2R1*g2R2*g3L2*g3R1*m1*m3*m4*amcha*p2p3*
     -     p3pfer)/amw**2 + 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*amcha*p1p4*p2p3*
     -     p3pfer)/amw**2 + 
     -  (4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m3*amcha*p1p4*p2p3*
     -     p3pfer)/amw**2 - 
     -  (4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pI*p2p3*
     -     p3pfer)/amw**2 - 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1pI*p2p3*
     -     p3pfer)/amw**2 - 
     -  (4*g1L2*g1R1*g2L2*g2R1*g3L2*g3R1*m1*m3*m4*mfer1*p2p3*
     -     p3pI)/amw**2 - 
     -  (4*g1L1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*m4*mfer1*p2p3*
     -     p3pI)/amw**2 + 
     -  (4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*mfer1*p1p4*p2p3*
     -     p3pI)/amw**2 + 
     -  (4*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*mfer1*p1p4*p2p3*
     -     p3pI)/amw**2 - 
     -  (4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*m4*p1pfer*p2p3*
     -     p3pI)/amw**2 - 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1pfer*p2p3*
     -     p3pI)/amw**2 + 
     -  4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*amcha*p1p2*p4pfer + 
     -  4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m3*amcha*p1p2*p4pfer - 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*amcha*p1p2*p2p3*
     -     p4pfer)/amw**2 - 
     -  (4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m3*amcha*p1p2*p2p3*
     -     p4pfer)/amw**2 - 
     -  (4*g1L1*g1R2*g2L2*g2R1*g3L2*g3R1*m3*amcha*p1p3*p2p3*
     -     p4pfer)/amw**2 - 
     -  (4*g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m3*amcha*p1p3*p2p3*
     -     p4pfer)/amw**2 - 
     -  4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2pI*p4pfer - 
     -  4*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*p2pI*p4pfer + 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p2pI*
     -     p4pfer)/amw**2 + 
     -  (4*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p3*p2pI*
     -     p4pfer)/amw**2 + 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p3pI*
     -     p4pfer)/amw**2 + 
     -  (4*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p3*p3pI*
     -     p4pfer)/amw**2 + 
     -  4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*mfer1*p1p2*p4pI + 
     -  4*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*mfer1*p1p2*p4pI - 
     -  (4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*mfer1*p1p2*p2p3*
     -     p4pI)/amw**2 - 
     -  (4*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*mfer1*p1p2*p2p3*
     -     p4pI)/amw**2 - 
     -  (4*g1L1*g1L2*g2L1*g2L2*g3L2*g3R1*m3*mfer1*p1p3*p2p3*
     -     p4pI)/amw**2 - 
     -  (4*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*mfer1*p1p3*p2p3*
     -     p4pI)/amw**2 + 
     -  4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2pfer*p4pI + 
     -  4*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*p2pfer*p4pI - 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p2pfer*
     -     p4pI)/amw**2 - 
     -  (4*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p3*p2pfer*
     -     p4pI)/amw**2 - 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p3pfer*
     -     p4pI)/amw**2 - 
     -  (4*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p3*p3pfer*
     -     p4pI)/amw**2 - 
     -  4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*m4*p1p2*pIpfer - 
     -  4*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1p2*pIpfer + 
     -  (4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*m4*p1p2*p2p3*
     -     pIpfer)/amw**2 + 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1p2*p2p3*
     -     pIpfer)/amw**2 + 
     -  (4*g1L2*g1R1*g2L1*g2L2*g3L2*g3R1*m3*m4*p1p3*p2p3*
     -     pIpfer)/amw**2 + 
     -  (4*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*p1p3*p2p3*
     -     pIpfer)/amw**2 - 
     -  4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p4*pIpfer - 
     -  4*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p4*pIpfer + 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p2p4*
     -     pIpfer)/amw**2 + 
     -  (4*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p3*p2p4*
     -     pIpfer)/amw**2 + 
     -  (4*g1L1*g1L2*g2L2*g2R1*g3L2*g3R1*m1*m3*p2p3*p3p4*
     -     pIpfer)/amw**2 + 
     -  (4*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*p2p3*p3p4*
     -     pIpfer)/amw**2
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference charged Higgs squark and charged Higgs Sm fermion
c----- g1l1 g1r1 =neutralino squark quark
c---- g12 chraged Higgs squark squark
c---- g2l1, g2r1= quarks charged higgs
c---- g2l2, g2r2= neutralino squark quark
c---- g3r1 charged Higgs fermion fermion
c---- g3r2 charged Higgs fermion fermion

       double precision function ampintchaHiSmferchaHisquark(g12,
     . g1r1, g1l1,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, 
     . p1pfer,p4pfer, pfer2, mfer)
      implicit none
      double precision g1l1, g12, g2l1, g2l2, g3r2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,p1pi, p2pI, p3pI, 
     .p4pI, pI2,p3p4,mfer,
     . g3r1
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHiSmferchaHisquark=

     - 8*g12*g3R2*g3R1*p2p3*
     -  (-(g1R1*g2L2*g2R1*m1*m4*mfer)/2. - 
     -    (g1L1*g2L1*g2R2*m1*m4*mfer)/2. + 
     -    (g1L1*g2L1*g2L2*mfer*p1p4)/2. + 
     -    (g1R1*g2R1*g2R2*mfer*p1p4)/2. - 
     -    (g1R1*g2L1*g2L2*m4*p1pfer)/2. - 
     -    (g1L1*g2R1*g2R2*m4*p1pfer)/2. + 
     -    (g1L1*g2L2*g2R1*m1*p4pfer)/2. + 
     -    (g1R1*g2L1*g2R2*m1*p4pfer)/2.)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference charged Higgs chargino and charged Higgs Sm fermion
c----- g1l1 g1r1 =neutralino squark quark
c---- g1l2, g1r2 chargino squark squark
c---- g2l1, g2r1= quarks charged higgs
c---- g2l2, g2r2= neutralino chargino charged Higgs
c---- g3r1 charged Higgs fermion fermion
c---- g3r2 charged Higgs fermion fermion

       double precision function ampintchaHiSmferchaHicharg(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfer, amcha)
      implicit none
      double precision g1l1, g1l2, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,p1pi, p2pI, p3pI, 
     .p4pI, pI2,p3p4,mfer, amcha, pipfer,
     . g3r1, g3r2
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHiSmferchaHicharg=

     -  -8*g3R2*g3R1*p2p3*
     -  (-(amcha*g1L2*g1R1*g2L2*g2R1*m1*m4*mfer)/2. - 
     -    (amcha*g1L1*g1R2*g2L1*g2R2*m1*m4*mfer)/2. + 
     -    (amcha*g1L1*g1L2*g2L1*g2L2*mfer*p1p4)/2. + 
     -    (amcha*g1R1*g1R2*g2R1*g2R2*mfer*p1p4)/2. - 
     -    (amcha*g1L2*g1R1*g2L1*g2L2*m4*p1pfer)/2. - 
     -    (amcha*g1L1*g1R2*g2R1*g2R2*m4*p1pfer)/2. + 
     -    (g1R1*g1R2*g2L2*g2R1*m4*mfer*p1pI)/2. + 
     -    (g1L1*g1L2*g2L1*g2R2*m4*mfer*p1pI)/2. + 
     -    (amcha*g1L1*g1L2*g2L2*g2R1*m1*p4pfer)/2. + 
     -    (amcha*g1R1*g1R2*g2L1*g2R2*m1*p4pfer)/2. - 
     -    (g1L1*g1R2*g2L1*g2L2*m1*mfer*p4pI)/2. - 
     -    (g1L2*g1R1*g2R1*g2R2*m1*mfer*p4pI)/2. + 
     -    (g1R1*g1R2*g2L1*g2L2*m1*m4*pIpfer)/2. + 
     -    (g1L1*g1L2*g2R1*g2R2*m1*m4*pIpfer)/2. - 
     -    (g1L1*g1R2*g2L2*g2R1*
     -       (p1pI*p4pfer + p1pfer*p4pI - p1p4*pIpfer))/2. - 
     -    (g1L2*g1R1*g2L1*g2R2*
     -       (p1pI*p4pfer + p1pfer*p4pI - p1p4*pIpfer))/2.)
       return
       end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference slepton chargino and charged Higgs Sm fermion
c----- g1l1 g1r1 =neutralino squark quark
c---- g1l2, g1r2 chargino squark squark
c---- g2l1, g2r1= quarks charged higgs
c---- g2l2, g2r2= neutrino chargino slepton
c---- g3r1 charged Higgs fermion fermion
c---- g3r2 g3l2 lepton slepton neutralino

       double precision function ampintchaHiSmferslepton(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2r2,g3r1,g3l2,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfer, amcha)
      implicit none
      double precision g1l1, g1l2, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,p1pi, p2pI, p3pI, 
     .p4pI, pI2,p3p4,mfer, amcha, pipfer, g3r2,
     . g3r1
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHiSmferslepton=


     -  -4*(-(amcha*g1R1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*m4*mfer*
     -        p1p2)/2. - (amcha*g1L1*g1R2*g2L1*g2R2*g3R1*g3R2*
     -       m1*m4*mfer*p2p3)/2. + 
     -    (amcha*g1L1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*mfer*p2p4)/
     -     2. - (amcha*g1R1*g1R2*g2L1*g2R2*g3L2*g3R1*m1*m3*m4*
     -       p2pfer)/2. + 
     -    (g1L2*g1R1*g2R1*g2R2*g3L2*g3R1*m1*m3*m4*mfer*p2pI)/
     -     2. + (amcha*g1R1*g1R2*g2R1*g2R2*g3R1*g3R2*mfer*
     -       (p1p4*p2p3 - p1p3*p2p4 + p1p2*p3p4))/2. - 
     -    (amcha*g1L1*g1R2*g2R1*g2R2*g3R1*g3R2*m4*
     -       (p1pfer*p2p3 - p1p3*p2pfer + p1p2*p3pfer))/2. + 
     -    (g1L1*g1L2*g2L1*g2R2*g3R1*g3R2*m4*mfer*
     -       (p1pI*p2p3 + p1p3*p2pI - p1p2*p3pI))/2. + 
     -    (amcha*g1L1*g1R2*g2R1*g2R2*g3L2*g3R1*m3*
     -       (p1pfer*p2p4 - p1p4*p2pfer + p1p2*p4pfer))/2. + 
     -    (amcha*g1R1*g1R2*g2L1*g2R2*g3R1*g3R2*m1*
     -       (p2pfer*p3p4 - p2p4*p3pfer + p2p3*p4pfer))/2. - 
     -    (g1L1*g1L2*g2L1*g2R2*g3L2*g3R1*m3*mfer*
     -       (p1pI*p2p4 + p1p4*p2pI - p1p2*p4pI))/2. - 
     -    (g1L2*g1R1*g2R1*g2R2*g3R1*g3R2*m1*mfer*
     -       (p2pI*p3p4 - p2p4*p3pI + p2p3*p4pI))/2. + 
     -    (g1L2*g1R1*g2L1*g2R2*g3L2*g3R1*m3*m4*
     -       (p1pI*p2pfer + p1pfer*p2pI - p1p2*pIpfer))/2. + 
     -    (g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*m1*m4*
     -       (p2pI*p3pfer - p2pfer*p3pI + p2p3*pIpfer))/2. - 
     -    (g1L1*g1L2*g2R1*g2R2*g3L2*g3R1*m1*m3*
     -       (p2pI*p4pfer - p2pfer*p4pI + p2p4*pIpfer))/2. - 
     -    (g1L2*g1R1*g2L1*g2R2*g3R1*g3R2*
     -       (p1pI*p2pfer*p3p4 + p1pfer*p2pI*p3p4 - 
     -         p1pI*p2p4*p3pfer - p1p4*p2pI*p3pfer - 
     -         p1pfer*p2p4*p3pI + p1p4*p2pfer*p3pI + 
     -         p1pI*p2p3*p4pfer + p1p3*p2pI*p4pfer - 
     -         p1p2*p3pI*p4pfer + p1pfer*p2p3*p4pI - 
     -         p1p3*p2pfer*p4pI + p1p2*p3pfer*p4pI - 
     -         p1p4*p2p3*pIpfer + p1p3*p2p4*pIpfer - 
     -         p1p2*p3p4*pIpfer))/2.)
         return
         end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c----- interference sneutrino chargino and charged Higgs Sm fermion
c----- g1l1 g1r1 =neutralino squark quark
c---- g1l2, g1r2 chargino squark squark
c---- g2l1, g2r1= quarks charged higgs
c---- g2l2, g2r2= sneutrino chargino lepton
c---- g3r1 charged Higgs fermion fermion
c---- g3r2 g3l2 neutrino sneutrino neutralino

       double precision function ampintchaHiSmfersneutrino(g1l1,
     . g1r1, g1l2,  g1r2,g2l1, g2r1, g2l2,g2r2,g3r1,g3r2,
     .  m1, m4, p1p2, p1p3, p1p4, p2p3, p2p4, p3p4, p1pfer,
     . p4pfer, pfer2,p1pI, p2pI, p3pI, p4pI, pIpfer, mfer, amcha)
      implicit none
      double precision g1l1, g1l2, g2l1, g2l2, g3l2,g2r2, g2r1,g1r1, 
     .  m1, m3, m4, mI1, mI2, p1p2, p1p3, p1p4, p2p3, p2p4, g1r2, 
     . p1pfer, p2pfer, p3pfer, p4pfer, pfer2,p1pi, p2pI, p3pI, 
     .p4pI, pI2,p3p4,mfer, amcha, pipfer,
     . g3r1, g3r2
      double precision sdgf,amz,amw,pi,g2 
      COMMON/SD_param/sdgf,amz,amw,pi,g2

      ampintchaHiSmfersneutrino=

     - -4*(-(amcha*g1L1*g1R2*g2L1*g2R2*g3R1*g3R2*m1*m4*mfer*
     -        p2p3)/2. + (amcha*g1R1*g1R2*g2R1*g2R2*g3R1*g3R2*
     -       mfer*(p1p4*p2p3 - p1p3*p2p4 + p1p2*p3p4))/2. - 
     -    (amcha*g1L1*g1R2*g2R1*g2R2*g3R1*g3R2*m4*
     -       (p1pfer*p2p3 - p1p3*p2pfer + p1p2*p3pfer))/2. + 
     -    (g1L1*g1L2*g2L1*g2R2*g3R1*g3R2*m4*mfer*
     -       (p1pI*p2p3 + p1p3*p2pI - p1p2*p3pI))/2. + 
     -    (amcha*g1R1*g1R2*g2L1*g2R2*g3R1*g3R2*m1*
     -       (p2pfer*p3p4 - p2p4*p3pfer + p2p3*p4pfer))/2. - 
     -    (g1L2*g1R1*g2R1*g2R2*g3R1*g3R2*m1*mfer*
     -       (p2pI*p3p4 - p2p4*p3pI + p2p3*p4pI))/2. + 
     -    (g1L1*g1L2*g2R1*g2R2*g3R1*g3R2*m1*m4*
     -       (p2pI*p3pfer - p2pfer*p3pI + p2p3*pIpfer))/2. - 
     -    (g1L2*g1R1*g2L1*g2R2*g3R1*g3R2*
     -       (p1pI*p2pfer*p3p4 + p1pfer*p2pI*p3p4 - 
     -         p1pI*p2p4*p3pfer - p1p4*p2pI*p3pfer - 
     -         p1pfer*p2p4*p3pI + p1p4*p2pfer*p3pI + 
     -         p1pI*p2p3*p4pfer + p1p3*p2pI*p4pfer - 
     -         p1p2*p3pI*p4pfer + p1pfer*p2p3*p4pI - 
     -         p1p3*p2pfer*p4pI + p1p2*p3pfer*p4pI - 
     -         p1p4*p2p3*pIpfer + p1p3*p2p4*pIpfer - 
     -         p1p2*p3p4*pIpfer))/2.)
        return
        end
