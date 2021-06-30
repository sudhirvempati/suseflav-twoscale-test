!==========================================================================================
****f* SuSeFLAV/oneloophiggs.f/higgs_analytical 
*  NAME
*    subroutine higgs_analytical
*  SYNOPSIS
*    Computes One loop correction(analytical) to cp even higgs boson and 
*    light higgs. 
*
*  FUNCTION
*     The routine calcultes one loop correction to cp even higgs boson.
*     hep-ph/ 9903404, hep-ph/0002213  - analytical expression for light higgs and 
*     CP- even higgs
*
*  INPUTS
*
*     alph3                          -  g3/(16 * pi^2)
*     mt                             -  running masses of top, bottom and tau
*     yuRG,ydRG,yeRG                 - (3 X 3) Yukawas
*     AURG                           - (3 X 3) Trilinear couplings
*     modmu                          - modulus of the \mu paramter 
*     vev1,vev2                      - vacuum expectation values of the two 
*                                      higgs doublet fields
*     M3t                            - Gaugino mass at msusy
*     tanbeta                        - the ratio of the vevs of the 
*                                      two Higgs doublet fields.
*     SUegg   =  6 eigenvalues of UP-Squark mass matrix.
*
*  RESULT
*     sigphi1,sigphi2,
*     sigphi12,sig2phi2
*     sig2phi2yuk         - elements of One loop corrections.
*  
*  EXAMPLE
*
*     SUBROUTINE higgs_analytical(MT,tanbeta,SUegg,AURG,sgnmu,modmu,
*     $     pizzT,piwwT,alph3,sigphi1,sigphi2,
*     $     sigphi12,sig2phi2,sig2phi2yuk)
*     
*  NOTES
*    1. q, the energy scale at which the corrections are added = msusy.
*    2. Conventions and notations followed are that of BPMZ.
*
*  BUGS
*    ---
*  SEE ALSO
*
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C

!======================================================================================

      SUBROUTINE higgs_analytical(MT,tanbeta,SUegg,AURG,sgnmu,modmu,
     $     alph3,sigphi1,sigphi2,
     $     sigphi12,sig2phi2,sig2phi2yuk)

      IMPLICIT NONE
      
      DOUBLE PRECISION MT, tanbeta, SUegg(6), AURG(3,3), modmu
      DOUBLE PRECISION sigphi1,sigphi2,sigphi12,alph3,beta
      DOUBLE PRECISION Msy, lda, sw2, Mzr, Mtlr,sgnmu
      double precision sig2phi2,sig2phi12,sig2phi1,sig2phi2yuk
      DOUBLE PRECISION sinsqthw_susy,alph,Gf,alphas

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
!      common/mwpole/ MWpole

      common/sinsq_susy/sinsqthw_susy

      common/gauge/alph,Gf,alphas

!----------------------------------------
      include 'stdinputs.h'
!---------------------------------------

      pi = 4.d0 *  datan(1.d0)

      MZ = MZpole

      beta = datan(tanbeta)
      
      Mzr = dsqrt((MZ*MZ))


      sw2 = sinsqthw_susy 

      lda = ((1.d0/8.d0) - sw2/3.d0 + ((4.d0 * sw2**2.d0)/9.d0))

      Msy = (SUegg(1)*SUegg(2) + Mt**2.d0 * (SUegg(1) + SUegg(2))
     $     + Mt**4.d0)**(0.25d0)

      Mtlr = (AURG(3,3) - sgnmu*modmu/tanbeta)
      

      sigphi1 = (Gf * dsqrt(2.d0) * Mzr**4.d0 * lda * 
     $     (dcos(beta))**2.d0 * 2.d0 * dlog(MT/Msy))
     $     / (pi * pi)

      sigphi12 = Gf * dsqrt(2.d0) * Mzr * Mzr * 
     $     ( (- 3.d0 * Mt**2.d0 /8.d0) +  (Mzr * Mzr * lda *
     $     (dsin(datan(tanbeta)))**2.d0)) * 2.d0 * dlog(Mt/Msy) 
     $     / (pi * pi * tanbeta)

      sigphi2 = (Gf * dsqrt(2.d0) * Mt**4.d0/
     $     (pi * pi * 8.d0 * dsin(datan(tanbeta))**2.d0)) *
     $     ((-2.d0*(Mzr/Mt)**2.d0) + (11.d0/10.d0 *(Mzr/Mt)**4.d0) +
     $     ((12.d0 - (6.d0*(Mzr*dsin(datan(tanbeta))/Mt)**2.d0) +
     $     (8.d0 * (Mzr*dsin(datan(tanbeta))/Mt)**4.d0 * lda) ) *
     $     2.d0 * dlog(Mt/Msy)) + 
     $     ((Mtlr/Msy)**2.d0 * (-12.d0 + (2.d0 * Mzr/ Mt)**2.d0 +
     $     6.d0 * (Mt/Msy)**2.d0)) +
     $     ((Mtlr/Msy)**4.d0 * (1.d0 - (2.d0*Mt/Msy)**2.d0 + 
     $     3.d0 * (Mt/Msy)**4.d0)) +
     $     ((Mtlr/Msy)**6.d0 * (3.d0*(Mt/Msy)**2.d0 - 
     $     12.d0*(Mt/Msy)**4.d0 +
     $     10.d0*(Mt/Msy)**6.d0))/5.d0 + 
     $     ((Mtlr/Msy)**8.d0 * (3.d0*(Mt/Msy)**4.d0 - 
     $     12.d0*(Mt/Msy)**6.d0+
     $     10.5d0*(Mt/Msy)**8.d0))/7.d0)

      sig2phi1  = 0.d0

      sig2phi12 = 0.d0

      sig2phi2  = ((Gf * dsqrt(2.d0) * 4.d0 * alph3 *  mt**4.d0)/
     $     (pi**3.d0 * dsin(beta)**2.d0)) * (3.d0 * 
     $     (dlog((mt/Msy)**2.d0))**2.d0 - 12.d0 * dlog(mt/Msy) - 
     $     6.d0 * (MtLR/Msy) - 6.d0 * (MtLR/Msy)**2.d0 * dlog(mt/Msy) +
     $     (3.d0/4.d0) * (MtLR/Msy)**4.d0)

      sig2phi2yuk = - (9.d0/(16.d0 * pi**4.d0)) * 
     $     ((Gf*Gf * mt**6.d0)/dsin(beta)**2.d0) * 
     $     ((2.d0 * dlog(mt/Msy))**2.d0 + 4.d0 * (MtLR/Msy)**2.d0 * !<--- sign before the second term reversed... 
     $     dlog(mt/Msy) + (1.d0/3.d0) * (MtLR/Msy)**4.d0 * dlog(mt/Msy))
      

      RETURN
      END SUBROUTINE 

!===========================================================================

!--------------------------------------------------------------------------------
!
!      ONE LOOP SELF-ENERGY TERMS FOR THE TWO CP-EVEN HIGGS BOSONS
!
!--------------------------------------------------------------------------------
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  1. couplings checked.
C  2. sign typos fixed.
C  3. all expressions checked - 20th May, 2010 
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
!-------------------------------------------------------------------------------
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  1.  check sign of \mu. ? 
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

!===============================================================================

      subroutine pis1s1(p,q,g,gp,mt,mb,mtau,tanbeta,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,ON,OCL,
     $     OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,mA0sq,vev1,vev2,
     $     pis1s1ans)

      implicit none
      integer i,j
      double precision p,q,g,gp,pis1s1ans,stops,sbots,staus
      double precision sneutrinos
      double precision sleps,sups,sdowns,fermions,higgs,charginos
      double precision neutralinos

      double precision mT, mB, mTau, tanbeta,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      double precision rot2dt(2,2),trot2dt(2,2),rot2db(2,2),trot2db(2,2)
      double precision rot2dtau(2,2),trot2dtau(2,2)
      double precision lTS112(2,2),lBS112(2,2) 
      double precision lTauS112(2,2)
      double precision a0mb,a0mtau
      
      double precision a0MZ,a0MW 

      double precision ON(4,4),OCL(2,2),OCR(2,2),ONdag(4,4),OCLdag(2,2)
      double precision OCRdag(2,2)
      double precision ht,htsq,hb,htau,hbsq,htausq
      double precision modmu,mh0,mHpm,mhu0,mA0

      double precision b0mbmb,b0mtaumtau,a0mstop1,a0mstop2
      double precision a0msbot1,a0msbot2,a0mstau1,a0mstau2,a0msnu(3)

      double precision costhtsq,sinthtsq,costhbsq,sinthbsq,costhtausq
      double precision sinthtausq,lSnuS1LL,cossqbeta,sinsqbeta,sinbeta

      double precision a0muL(2),a0muR(2),a0mdL(2),a0mdR(2)
      double precision a0meL(2),a0meR(2),b0msnu(3)
      double precision ls1tt(2,2),ls1bb(2,2)
      double precision ls1tautau(2,2)

      double precision b0mt(2,2),b0mb(2,2),b0mtau(2,2),b0muL(2),b0muR(2)
      double precision b0mdL(2),b0mdR(2),b0meL(2),b0meR(2)

      double precision fmHpmMW,fmAMZ,fMWMW,fMZMZ
      double precision b0MWMW,b0MZMZ !,a0MW,a0MZ 

      double precision ls1eeLL,ls1eeRR,ls1ddLL,ls1ddRR,ls1uuLL,ls1uuRR

      double precision hhs1(4,4),hhs1s1(4),hphps1(2,2),hphps1s1(2)
      data hhs1/ 16 * 0.d0/, hhs1s1/ 4 * 0.d0/,hphps1/ 4 * 0.d0/
      data hphps1s1/ 2 * 0.d0/

      double precision b0HH(4,4),a0H(4),mG0,mGp,b0HpHp(2,2),a0Hp(2)

      double precision fChiChis1s1(4, 4), gChiChis1s1(4, 4),gmneut(4,4)
      double precision b0mneut(4,4)
      DOUBLE PRECISION aPsi(4,4),aChi(4,4),bChi(4,4)
      data aPsi/ 16 * 0.d0/, aChi/ 16 * 0.d0/, bChi/ 16 * 0.d0/

      double precision aPsic(2, 2),aPsicdag(2,2),gmch(2,2),b0mch(2,2)
      data aPsic/ 4 * 0.d0/,aPsicdag/ 4 * 0.d0/,gmch/ 4 * 0.d0/

      DOUBLE PRECISION aChic(2, 2),bChic(2, 2)
      data aChic/ 4 * 0.d0/,bChic/ 4 * 0.d0/

      double precision pis1s1p,AURG(3,3),ADRG(3,3)
      DOUBLE PRECISION AERG(3,3),SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq

      DOUBLE PRECISION muL,muR,mcL,mcR,mt1,mt2
      DOUBLE PRECISION mdL,mdR,msL,msR,mb1,mb2
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtau1,mtau2      
      DOUBLE PRECISION mneut(4),msnu(3),mchargino(2)

      double precision cossqthw,sinsqthw,cosbeta,costhw,sinthw
      double precision alpha,thw,tan2beta
      double precision thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae,mtL,mtR
      double precision yuRg(3,3),ydRG(3,3),yeRg(3,3)

      double precision beta
      double precision vev1,vev2,cos2beta,sin2beta

      DOUBLE PRECISION sinsqthw_susy,MZrun, MWrun

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

!-----------------------------------------------------------------------

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,
     $     thetamu,thetau,thetad,thetae

      common/gbrunning/ MZrun, MWrun
!-----------------------------------------------------------------------

      external dag2d,dag4d,mat3prod2d,a0,theta,b0,f,funcg

      include 'stdinputs.h'

!-------------------------------------------------
C     Nomenclature
!-------------------------------------------------


      pi = 4.d0 * datan(1.d0)

!      print("(G,ES15.4,ES15.4)"), "MZrun, MWrun = ", MZrun, MWrun

      MZ = MZrun                !MZpole
      MW = MWrun                !MWpole

      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mt1 = dsqrt(SUegg(2))
      mt2 = dsqrt(SUegg(1))

      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mb1 = dsqrt(SDegg(2))
      mb2 = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtau1 = dsqrt(SLegg(2))
      mtau2 = dsqrt(SLegg(1))

       
      msnu(1) = dsqrt(SNegg(1))
      msnu(2) = dsqrt(SNegg(2))
      msnu(3) = dsqrt(SNegg(3))

      mneut(1) = Neg(1) 
      mneut(2) = Neg(2)
      mneut(3) = Neg(3)
      mneut(4) = Neg(4)

      mchargino(1) = Ceg(1)
      mchargino(2) = Ceg(2)

      mh0  = dsqrt((mh0sq))
      mHu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      mG0 = MZ
      mGp = MW

!---------------------------------------------------------
C     Definitions.
!--------------------------------------------------------

      beta     = datan(tanbeta)
      cosbeta  = dcos(beta)
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0*beta)
      cos2beta = dcos(2.d0*beta)
      tan2beta = dtan(2.d0*beta)

      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/
     $     (mA0sq - MZ*MZ))*tan2beta)
      
      sinsqthw = sinsqthw_susy !(1.d0 - (MW/MZ)**2.d0)
      sinthw = dsqrt(sinsqthw)
      thw = dasin(dsqrt(sinsqthw))
      costhw = dcos(thw)
      cossqthw = dcos(thw)*dcos(thw)

      gnuL =   0.5d0
      guL  =   0.5d0 - 2.d0*sinsqthw/ 3.d0 
      gdL  =  -0.5d0 + sinsqthw/3.d0
      geL  =  -0.5d0 + sinsqthw
      guR  =   2.d0*sinsqthw/3.d0 
      gdR  =  -sinsqthw/3.d0
      geR  =  -sinsqthw 
      yuL  =   1.d0/3.d0 
      yuR  =  -4.d0/3.d0 
      ydL  =   1.d0/3.d0 
      ydR  =   2.d0/3.d0 
      yeL  =  -1.d0 
      yeR  =   2.d0 
      ynuL =  -1.d0 

      ht   =  yuRG(3,3) 
      hb   =  ydRg(3,3) 
      htau =  yeRG(3,3) 

      htsq =   ht*ht
      hbsq =   hb*hb
      htausq = htau*htau

!-----------------------------------
     
      rot2dt(1,1) =   dcos(thetat)
      rot2dt(1,2) =   dsin(thetat)
      rot2dt(2,1) = - rot2dt(1,2)
      rot2dt(2,2) =   rot2dt(1,1)

      rot2db(1,1) =   dcos(thetab)
      rot2db(1,2) =   dsin(thetab)
      rot2db(2,1) = - rot2db(1,2)
      rot2db(2,2) =   rot2db(1,1)

      rot2dtau(1,1) =   dcos(thetatau)
      rot2dtau(1,2) =   dsin(thetatau)
      rot2dtau(2,1) = - rot2dtau(1,2)
      rot2dtau(2,2) =   rot2dtau(1,1)

!-------------------------------------
      
      call dag2d(rot2dt,trot2dt)
      call dag2d(rot2db,trot2db)
      call dag2d(rot2dtau,trot2dtau)

!------------------------------------------------------------------------------
C     3rd family fermions 
!------------------------------------------------------------------------------

      costhtsq = (dcos(thetat))**2.d0
      sinthtsq = (dsin(thetat))**2.d0
      costhbsq = (dcos(thetab))**2.d0
      sinthbsq = (dsin(thetab))**2.d0
      costhtausq = (dcos(thetatau))**2.d0
      sinthtausq = (dsin(thetatau))**2.d0


!---------------------------------------------------  
    
      call a0(mb,q,a0mb)
      call a0(mtau,q,a0mtau)        
      call a0(mt1,q,a0mstop1)
      call a0(mt2,q,a0mstop2)
      call a0(mb1,q,a0msbot1)
      call a0(mb2,q,a0msbot2)
      call a0(mtau1,q,a0mstau1)
      call a0(mtau2,q,a0mstau2)
      
      call b0(p,mb,mb,q,b0mbmb)

!      print*,"b0mbmb in pis1s1 = ", b0mbmb
      call b0(p,mtau,mtau,q,b0mtaumtau)      
!      print*,"b0mtaumtau in pis1s1 = ", b0mtaumtau

!---------------------------------------------------

      fermions = 3.d0*hbsq*((p*p - 4.d0 * mb*mb) * b0mbmb - 
     $     2.d0 * a0mb)
      
      fermions = fermions + htausq*((p*p - 4.d0 * mtau*mtau) *
     $     b0mtaumtau - 2.d0 * a0mtau)
 
!      print*,"hb, htau = ", hb, htau
!      print*,"mb, mtau = ", mb, mtau

      sbots = 3.d0 * hbsq * (a0msbot1 + a0msbot2)

      staus =  htausq * (a0mstau1 + a0mstau2)

c$$$      print*,"stau1 = ", staus

      stops = (3.d0 * g*g / (2.d0 * cossqthW)) *
     $     (guR *(costhtsq * a0mstop2 + sinthtsq * a0mstop1) +
     $     guL * (sinthtsq * a0mstop2 + costhtsq * a0mstop1))

      sbots = sbots + (3.d0 * g*g / (2.d0 * cossqthW)) *
     $     (gdR *(costhbsq * a0msbot2 + sinthbsq * a0msbot1) +
     $     gdL * (sinthbsq * a0msbot2 + costhbsq * a0msbot1))

      staus = staus + (1.d0 * g*g / (2.d0 * cossqthW)) *
     $     (geR *(costhtausq * a0mstau2 + sinthtausq * a0mstau1) +
     $     geL * (sinthtausq * a0mstau2 + costhtausq * a0mstau1))

c$$$      print*,"a0mstau1 = ", a0mstau1, "a0mstau2 = ", a0mstau2
c$$$      print*,"ctau = ", dsqrt(costhtausq), ", stau = ", 
c$$$     $     dsqrt(sinthtausq),dsin(thetatau),dcos(thetatau)
c$$$
c$$$      print*,"stau2 = ", staus

      sups = 0.d0
      sdowns = 0.d0
      sneutrinos = 0.d0
      sleps = 0.d0

!------------------------------------
      call a0(muL,q,a0muL(1))
      call a0(muR,q,a0muR(1))
      call a0(mcL,q,a0muL(2))
      call a0(mcR,q,a0muR(2))
      call a0(mdL,q,a0mdL(1))
      call a0(mdR,q,a0mdR(1))
      call a0(msL,q,a0mdL(2))
      call a0(msR,q,a0mdR(2))
      call a0(meL,q,a0meL(1))
      call a0(meR,q,a0meR(1))
      call a0(mmuL,q,a0meL(2))
      call a0(mmuR,q,a0meR(2))
      call a0(msnu(1),q,a0msnu(1))
      call a0(msnu(2),q,a0msnu(2))
      call a0(msnu(3),q,a0msnu(3))
!-----------------------------------------

      loop1: do i = 1, 2

      sups = sups + 3.d0 * (g*g / (2.d0 * cossqthW)) *
     $     (guL * a0muL(i) + guR * a0muR(i))
      sdowns = sdowns + 3.d0 * (g*g / (2.d0 * cossqthW)) *
     $     (gdL * a0mdL(i) + gdR * a0mdR(i))
      sleps = sleps + (g*g / (2.d0 * cossqthW)) * 
     $     (geL * a0meL(i) + geR * a0meR(i))

      enddo loop1

!-------------------------------------------------------------------------------
C     stop couplings
          
      ls1tt(1, 1) = g * MZ * guL * cosbeta / costhw
      ls1tt(1, 2) = -1.d0 * ht * sgnmu * modmu / dsqrt(2.d0)
      ls1tt(2, 1) = ls1tt(1, 2)
      ls1tt(2, 2) = g * MZ * guR * cosbeta / costhW

!-----------------------------------------------------------------------
C     sbottom couplings

      ls1bb(1, 1) = (g * MZ * gdL * cosbeta / costhW) + hbsq*vev1
      ls1bb(1, 2) = hb*ADRG(3,3)/dsqrt(2.d0) 
      ls1bb(2, 1) = ls1bb(1, 2)
      ls1bb(2, 2) = (g * MZ * gdR * cosbeta / costhW) + hbsq*vev1

!-------------------------------------------------------------------------
!     stau couplings

      ls1tautau(1, 1) = (g * MZ * geL * cosbeta / costhW) + htausq*vev1
      ls1tautau(1, 2) = htau*AERG(3,3) / dsqrt(2.d0) 
      ls1tautau(2, 1) = ls1tautau(1, 2)
      ls1tautau(2, 2) = (g * MZ * geR * cosbeta / costhW) + htausq*vev1

c$$$      print*,"htau*AERG(3,3) = ", htau*AERG(3,3)

!----------------------------------------------------------------------      
C     3rd family fermions mixings

      call mat3prod2d(rot2dt,ls1tt,trot2dt,lTS112)
      call mat3prod2d(rot2db,ls1bb,trot2db,lBS112)      
      call mat3prod2d(rot2dtau,ls1tautau,trot2dtau,lTauS112)      

!---------------------------------------------------------------------
      
      call b0(p,mt1,mt1,q,b0mt(1,1))
      call b0(p,mt1,mt2,q,b0mt(1,2))
      call b0(p,mt2,mt1,q,b0mt(2,1))
      call b0(p,mt2,mt2,q,b0mt(2,2))

      call b0(p,mb1,mb1,q,b0mb(1,1))
      call b0(p,mb1,mb2,q,b0mb(1,2))
      call b0(p,mb2,mb1,q,b0mb(2,1))
      call b0(p,mb2,mb2,q,b0mb(2,2))

      call b0(p,mtau1,mtau1,q,b0mtau(1,1))
      call b0(p,mtau1,mtau2,q,b0mtau(1,2))
      call b0(p,mtau2,mtau1,q,b0mtau(2,1))
      call b0(p,mtau2,mtau2,q,b0mtau(2,2))
      
      call b0(p,muL,muL,q,b0muL(1))      
      call b0(p,mcL,mcL,q,b0muL(2))      
      call b0(p,muR,muR,q,b0muR(1))      
      call b0(p,mcR,mcR,q,b0muR(2))      
      call b0(p,mdL,mdL,q,b0mdL(1))      
      call b0(p,msL,msL,q,b0mdL(2))      
      call b0(p,mdR,mdR,q,b0mdR(1))      
      call b0(p,msR,msR,q,b0mdR(2))      
      call b0(p,meL,meL,q,b0meL(1))      
      call b0(p,mmuL,mmuL,q,b0meL(2))      
      call b0(p,meR,meR,q,b0meR(1))      
      call b0(p,mmuR,mmuR,q,b0meR(2))      

      call b0(p,msnu(1),msnu(1),q,b0msnu(1))      
      call b0(p,msnu(2),msnu(2),q,b0msnu(2))      
      call b0(p,msnu(3),msnu(3),q,b0msnu(3))      
!----------------------------------------------

      loop2: do i = 1, 2
      loop3: do j = 1, 2

      stops = stops + 3.d0 * lTS112(i, j)* lTS112(i, j) * b0mt(i,j)
      sbots = sbots + 3.d0 * lBS112(i, j)* lBS112(i, j) * b0mb(i,j)
      staus = staus + lTauS112(i, j) * lTauS112(i, j) * b0mtau(i,j)

      enddo loop3
      enddo loop2

c$$$      print*,"stau3 = ", staus

!-------------------------------------------------------------------------
C     1st and 2nd generation couplings
!-------------------------------------------------------------------------

      ls1eeLL = (g * MZ * geL * cosbeta) / costhW
      ls1eeRR = (g * MZ * geR * cosbeta) / costhW
      ls1ddLL = (g * MZ * gdL * cosbeta) / costhW
      ls1ddRR = (g * MZ * gdR * cosbeta) / costhW
      ls1uuLL = (g * MZ * guL * cosbeta) / costhW
      ls1uuRR = (g * MZ * guR * cosbeta) / costhW
      
      lSnuS1LL = (g*MZ*gnuL*cosbeta/costhW)**2.d0


      loop4: do i = 1, 2

      sups = sups + 3.d0 * ls1uuLL * ls1uuLL * b0muL(i) + 
     $     3.d0 * ls1uuRR * ls1uuRR * b0muR(i)

      sdowns = sdowns + 3.d0 * ls1ddLL * ls1ddLL * b0mdL(i) 
     $     + 3.d0 * ls1ddRR * ls1ddRR * b0mdR(i)

      sleps = sleps + ls1eeLL * ls1eeLL * b0meL(i) + 
     $     ls1eeRR * ls1eeRR * b0meR(i)
      
      enddo loop4

      sneutrinos = sneutrinos + (g*g / (2.d0 * cossqthW)) * gnuL * 
     $     (a0msnu(1) + a0msnu(2) + a0msnu(3)) +
     $     lSnuS1LL * (b0msnu(1) + b0msnu(2) + b0msnu(3))

!------------------------------------     
      call f(p,mHpm,MW,q,fmHpmMW)
      call f(p,mA0,MZ,q,fmAMZ)
      call f(p,MW,MW,q,fMWMW)
      call f(p,MZ,MZ,q,fMZMZ)

      call b0(p,MW,MW,q,b0MWMW)      
      call b0(p,MZ,MZ,q,b0MZMZ)

      call a0(MW,q,a0MW)
      call a0(MZ,q,a0MZ)
!----------------------------------

      
      cossqbeta = cosbeta*cosbeta
      sinsqbeta = sinbeta*sinbeta

      higgs = 0.d0
      
      higgs = g *g * 0.25d0 * (sinsqbeta * (2.d0 * fmHpmMW + 
     $     (fmAMZ/ cossqthW)) + cossqbeta * (2.d0 * fMwMW  + 
     $     (fMZMZ / cossqthW))) + 1.75d0 * g * g * cossqbeta * 
     $     (2.d0 * MW*MW * b0MWMW + (MZ*MZ*b0MZMZ / cossqthW)) +
     $     g*g * (2.d0 * a0MW + (a0MZ / cossqthW))

!------------------------------------------------------------------------------
      
      do i = 1,4
         do j = 1,4
            hhs1(i,j) = 0.d0
         enddo
         hhs1s1(i) = 0.d0
      enddo
      
      hhs1(1, 1) =  cosbeta * ((3.d0 * cos(alpha)**2.d0) - 
     $      (sin(alpha)**2.d0)) - sinbeta * sin(2.d0*alpha)

      hhs1(2, 2) =  cosbeta * ((3.d0 * sin(alpha)**2.d0) -
     $     (cos(alpha)**2.d0)) + sinbeta * sin(2.d0*alpha)
      hhs1(1, 2) = - 2.d0 * cosbeta * sin(2.d0*alpha) - 
     $     (sinbeta * cos(2.d0*alpha))
      hhs1(2, 1) =  hhs1(1, 2)

      hhs1(3, 3) =  cos2beta * cosbeta
      hhs1(3, 4) = -sin2beta * cosbeta
      hhs1(4, 3) =  hhs1(3, 4)
      hhs1(4, 4) = -cos2beta * cosbeta

      loop5: do i = 1, 4
      loop6: do j = 1, 4

      hhs1(i,j) = hhs1(i,j) * (g * MZ / (2.d0 * costhW))
      
      enddo loop6
      enddo loop5

!--------------------------------------------------------------------

      hhs1s1(1) = (3.d0 * cos(alpha)**2.d0) - (sin(alpha)**2.d0)
      hhs1s1(2) = (3.d0 * sin(alpha)**2.d0) - (cos(alpha)**2.d0)
      hhs1s1(3) =   cos2beta
      hhs1s1(4) = - cos2beta
      
      loop7: do i = 1, 4
      
      hhs1s1(i) = hhs1s1(i) * ((g*g) * 0.25d0 / cossqthW)

      enddo loop7

!---------------------------------------      
      call b0(p,mHu0,mHu0,q,b0HH(1,1))   
      call b0(p,mHu0,mh0,q,b0HH(1,2))      
      call b0(p,mHu0,mG0,q,b0HH(1,3))      
      call b0(p,mHu0,mA0,q,b0HH(1,4))      

      call b0(p,mh0,mHu0,q,b0HH(2,1))      
      call b0(p,mh0,mh0,q,b0HH(2,2))      
      call b0(p,mh0,mG0,q,b0HH(2,3))      
      call b0(p,mh0,mA0,q,b0HH(2,4))      

      call b0(p,mG0,mHu0,q,b0HH(3,1))      
      call b0(p,mG0,mh0,q,b0HH(3,2))      
      call b0(p,mG0,mG0,q,b0HH(3,3))      
      call b0(p,mG0,mA0,q,b0HH(3,4))      

      call b0(p,mA0,mHu0,q,b0HH(4,1))      
      call b0(p,mA0,mh0,q,b0HH(4,2))      
      call b0(p,mA0,mG0,q,b0HH(4,3))      
      call b0(p,mA0,mA0,q,b0HH(4,4)) 
      
      call a0(mHu0,q,a0H(1))
      call a0(mh0,q,a0H(2))
      call a0(mG0,q,a0H(3))
      call a0(mA0,q,a0H(4))
!------------------------------------
      
      loop8: do i = 1, 4
      loop9: do j = 1, 4
      
      higgs = higgs + 0.5d0 * hhs1(i, j)* hhs1(i, j) * b0HH(i,j)

      enddo loop9

      higgs = higgs + 0.5d0 * hhs1s1(i) * a0H(i)

      enddo loop8

!-----------------------------------------------------------

      loophphps1i: do i = 1, 2
      loophphps1j: do j = 1, 2
      
      hphps1(i,j) = 0.d0

      enddo loophphps1j

      hphps1s1(i) = 0.d0

      enddo loophphps1i

 
      hphps1(1, 1) =  cos2beta * cosbeta
      hphps1(2, 2) = -cos2beta * cosbeta + 2.d0 * cossqthW * cosbeta
      hphps1(1, 2) = -sin2beta * cosbeta + cossqthW * sinbeta 
      hphps1(2, 1) =  hphps1(1, 2)
      
      loop10: do i = 1, 2
      loop11: do j = 1, 2
      
      hphps1(i,j) = hphps1(i,j) * (g * MZ * 0.5d0 / costhW)
      
      enddo loop11
      enddo loop10

!-----------------------------------------------      
 
      hphps1s1(1) = cossqthW + sinsqthW * cos2beta
      hphps1s1(2) = cossqthW - sinsqthW * cos2beta
      
      loop12: do i = 1, 2
      hphps1s1(i) = hphps1s1(i) * (g*g* 0.25d0 / cossqthW)
      enddo loop12

!--------------------------------------
      call b0(p,mGp,mGp,q,b0HpHp(1,1))      
      call b0(p,mGp,mHpm,q,b0HpHp(1,2))      
      call b0(p,mHpm,mGp,q,b0HpHp(2,1))      
      call b0(p,mHpm,mHpm,q,b0HpHp(2,2))      
      
      call a0(mGp,q,a0Hp(1))
      call a0(mHpm,q,a0Hp(2))

!--------------------------------------
            
      loop13: do i = 1, 2
      loop14: do j = 1, 2
      
      higgs = higgs + hphps1(i, j) * hphps1(i, j)* b0HpHp(i,j)

      enddo loop14

      higgs = higgs + hphps1s1(i)  * a0Hp(i)

      enddo loop13

!----------------------------------------------------------------------------
      
      aPsi(1, 3) = - gp * 0.5d0
      aPsi(3, 1) =   aPsi(1, 3) 
      aPsi(2, 3) =   g  * 0.5d0 
      aPsi(3, 2) =   aPsi(2, 3)

      call dag4d(ON,ONdag)
      call mat3prod4d(ON,aPsi,ONdag,aChi)
      call mat3prod4d(ON,aPsi,ONdag,bChi)  

      loopni: do i = 1, 4     
      loopnj: do j = 1, 4

      call b0(p,dabs(mneut(i)),dabs(mneut(j)),q,b0mneut(i,j))
      call funcg(p,dabs(mneut(i)),dabs(mneut(j)),q,gmneut(i,j))
            
      enddo loopnj
      enddo loopni

      neutralinos = 0.d0

      loop15: do i = 1, 4
      loop16: do j = 1, 4
      
      fChiChis1s1(i, j) = real(aChi(i, j) * aChi(i,j)) + 
     $     real(bChi(i, j) * bChi(i,j))

      gChiChis1s1(i, j) = real(bChi(i, j) * aChi(i, j)) + 
     $     real(aChi(i, j) * bChi(i, j))

      neutralinos = neutralinos + 0.5d0 * 
     $     ((fChiChis1s1(i, j) * gmneut(i,j)) - 2.d0 *
     $     (gChiChis1s1(i, j) * mneut(i) * mneut(j) * b0mneut(i,j)))

      enddo loop16
      enddo loop15
      
!-------------------------------------------------------------------------

      aPsic(1, 2) = g / dsqrt(2.d0)
      
      call dag2d(aPsic,aPsicdag)
      call dag2d(OCL,OCLdag)
      call dag2d(OCR,OCRdag)
      call mat3prod2d(OCR,aPsic,OCLdag,aChic)            !<--------------------CHECKED!!
      call mat3prod2d(OCL,aPsicdag,OCRdag,bChic)


      loop19: do i = 1, 2
      loop18: do j = 1, 2
      
      call funcg(p,mchargino(i),mchargino(j),q,gmch(i,j))
      call b0(p,mchargino(i),mchargino(j),q,b0mch(i,j))
        
      enddo loop18
      enddo loop19

      charginos = 0.d0   
      
      loop23: do i = 1, 2
      loop24: do j = 1, 2

      fChiChis1s1(i, j) = real(aChic(i, j) * (aChic(i,j)) +
     $     bChic(i, j) * (bChic(i,j)))

      gChiChis1s1(i, j) = real(bChic(i, j) * aChic(i, j) + 
     $     (aChic(i, j)) * bChic(i, j))

      charginos = charginos + 
     $     (fChiChis1s1(i, j) * gmch(i,j) - 2.d0 *
     $     gChiChis1s1(i, j) * mchargino(i) * mchargino(j) * b0mch(i,j))

      enddo loop24
      enddo loop23

!----------------------------------------------------------------------
      
c$$$      print*,"in pis1s1 p = ", p, " q = ",q
c$$$
c$$$      print*,"sups = ", sups, " sdowns = ", sdowns, " sleps = ", sleps,
c$$$     $     " stops = ", stops, " sbots = ", sbots, " staus = ", staus, 
c$$$     $     " sneutrinos = ", sneutrinos, " fermions = ", fermions, 
c$$$     $     " higgs = ", higgs, " neutralinos = ", neutralinos, 
c$$$     $     " charginos = ", charginos

      pis1s1p = (sups  + sdowns + sleps  + stops  + sbots  + staus  + 
     $     sneutrinos + fermions + higgs + neutralinos + charginos)
      
      
      
      pis1s1ans =  pis1s1p/ (16.d0*pi*pi)

!------------------------------------------------------------------------------

      RETURN      
      end subroutine pis1s1

!==============================================================================
 
!------------------------------------------------------------------------------

      subroutine pis2s2(p,q,g,gp,mt,mb,mtau,tanbeta,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,ON,OCL,
     $     OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,mA0sq,vev1,vev2,
     $     pis2s2ans)

      implicit none
      integer i,j
      double precision p,q,g,gp,pis2s2ans,stops,sbots,staus
      double precision sneutrinos
      double precision sleps,sups,sdowns,fermions,higgs,charginos
      double precision neutralinos

      double precision mT, mB, mTau, tanbeta, sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      double precision rot2dt(2,2),trot2dt(2,2),rot2db(2,2),trot2db(2,2)
      double precision rot2dtau(2,2),trot2dtau(2,2)
      double precision lTS212(2,2),lBS212(2,2) 
      double precision lTauS212(2,2)
      double precision a0mtau,a0mt



      double precision a0MZ,a0MW 

      double precision ON(4,4),OCL(2,2),OCR(2,2),ONdag(4,4),OCLdag(2,2)
      double precision OCRdag(2,2)
      double precision ht,htsq,hb,htau,hbsq,htausq
      double precision modmu,mh0,mHpm,mhu0,mA0

      double precision b0mtaumtau,a0mstop1,a0mstop2,b0mtmt
      double precision a0msbot1,a0msbot2,a0mstau1,a0mstau2,a0msnu(3)

      double precision costhtsq,sinthtsq,costhbsq,sinthbsq,costhtausq
      double precision sinthtausq,lSnuS2LL,cossqbeta,sinsqbeta,sinbeta

      double precision a0muL(2),a0muR(2),a0mdL(2),a0mdR(2)
      double precision a0meL(2),a0meR(2),b0msnu(3)
      double precision ls2tt(2,2),ls2bb(2,2)
      double precision ls2tautau(2,2)

      double precision b0mt(2,2),b0mb(2,2),b0mtau(2,2),b0muL(2),b0muR(2)
      double precision b0mdL(2),b0mdR(2),b0meL(2),b0meR(2)

      double precision fmHpmMW,fmAMZ,fMWMW,fMZMZ
      double precision b0MWMW,b0MZMZ !,a0MW,a0MZ 

      double precision ls2eeLL,ls2eeRR,ls2ddLL,ls2ddRR,ls2uuLL,ls2uuRR

      double precision hhs2(4,4),hhs2s2(4),hphps2(2,2),hphps2s2(2)
      data hhs2/ 16 * 0.d0/, hhs2s2/ 4 * 0.d0/, hphps2/ 4 * 0.d0/
      data hphps2s2/ 2 * 0.d0/

      double precision b0HH(4,4),a0H(4),mG0,mGp,b0HpHp(2,2),a0Hp(2)

      double precision fChiChis2s2(4, 4),gChiChis2s2(4, 4),gmneut(4,4)
      double precision b0mneut(4,4)

      DOUBLE PRECISION aPsi(4,4),aChi(4,4),bChi(4,4)
      data aPsi/ 16 * 0.d0/, aChi/ 16 * 0.d0/, bChi/ 16 * 0.d0/

      double precision aPsic(2, 2),aPsicdag(2,2),gmch(2,2),b0mch(2,2)
      data aPsic/ 4 * 0.d0/,aPsicdag/ 4 * 0.d0/,gmch/ 4 * 0.d0/

      DOUBLE PRECISION aChic(2, 2),bChic(2, 2)
      data aChic/ 4 * 0.d0/,bChic/ 4 * 0.d0/

      DOUBLE PRECISION AERG(3,3),AURG(3,3),ADRG(3,3)
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq

      DOUBLE PRECISION muL,muR,mcL,mcR,mt1,mt2
      DOUBLE PRECISION mdL,mdR,msL,msR,mb1,mb2
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtau1,mtau2      
      DOUBLE PRECISION mneut(4),msnu(3),mchargino(2)

      double precision cossqthw,sinsqthw,cosbeta,costhw,sinthw
      double precision alpha,thw,tan2beta
      double precision thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae,mtL, mtR
      double precision yuRg(3,3),ydRG(3,3),yeRg(3,3)

      double precision beta, checkhiggs1, checkhiggs2
      double precision checkhiggs3, checkhiggs4
      double precision vev1,vev2,cos2beta,sin2beta
      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

!------------------------------------------------------------------------

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

!-------------------------------------------------------------------------

      external dag2d,dag4d,mat3prod2d,a0,theta,b0,f,funcg

      include 'stdinputs.h'

!----------------------------------------
C     Nomenclature
!----------------------------------------

      pi = 4.d0 * datan(1.d0)

      MW = MWpole
      MZ = MZpole

      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mt1 = dsqrt(SUegg(2))
      mt2 = dsqrt(SUegg(1))

      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mb1 = dsqrt(SDegg(2))
      mb2 = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtau1 = dsqrt(SLegg(2))
      mtau2 = dsqrt(SLegg(1))

      msnu(1) = dsqrt(SNegg(1))
      msnu(2) = dsqrt(SNegg(2))
      msnu(3) = dsqrt(SNegg(3))

      mneut(1) = Neg(1) 
      mneut(2) = Neg(2)
      mneut(3) = Neg(3)
      mneut(4) = Neg(4)

      mchargino(1) = Ceg(1)
      mchargino(2) = Ceg(2)

      mh0  = dsqrt((mh0sq))
      mHu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = sqrt((mA0sq))

      mG0 = MZ
      mGp = MW

!---------------------------------------------
C     Definitions.
!---------------------------------------------      

      beta = datan(tanbeta)
      cosbeta = dcos(beta)
      sinbeta = dsin(beta)
      sin2beta = dsin(2.d0*beta)
      cos2beta = dcos(2.d0*beta)
      tan2beta = dtan(2.d0*beta)

      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/
     $     (mA0sq - MZ*MZ))*tan2beta)
      
      sinsqthw = sinsqthw_susy !(1.d0 - (MW/MZ)**2.d0)
      sinthw = dsqrt(sinsqthw)
      thw = dasin(dsqrt(sinsqthw))
      costhw = dcos(thw)
      cossqthw = dcos(thw)*dcos(thw)

      gnuL =   0.5d0
      guL  =   0.5d0 - 2.d0*sinsqthw/ 3.d0 
      gdL  =  -0.5d0 + sinsqthw/3.d0
      geL  =  -0.5d0 + sinsqthw
      guR  =   2.d0*sinsqthw/3.d0 
      gdR  =  -sinsqthw/3.d0
      geR  =  -sinsqthw 
      yuL  =   1.d0/3.d0 
      yuR  =  -4.d0/3.d0 
      ydL  =   1.d0/3.d0 
      ydR  =   2.d0/3.d0 
      yeL  =  -1.d0 
      yeR  =   2.d0 
      ynuL =  -1.d0 

      ht   =  yuRG(3,3)
      hb   =  ydRg(3,3) 
      htau =  yeRG(3,3) 
      
!      print*,"ht, hb, htau = ", ht, hb, htau

      htsq =  ht*ht
      hbsq =  hb*hb
      htausq = htau*htau


!----------------------------------------------------------------
     
      rot2dt(1,1) =   dcos(thetat)
      rot2dt(1,2) =   dsin(thetat)
      rot2dt(2,1) = - rot2dt(1,2)
      rot2dt(2,2) =   rot2dt(1,1)

      rot2db(1,1) =   dcos(thetab)
      rot2db(1,2) =   dsin(thetab)
      rot2db(2,1) = - rot2db(1,2)
      rot2db(2,2) =   rot2db(1,1)

      rot2dtau(1,1) =   dcos(thetatau)
      rot2dtau(1,2) =   dsin(thetatau)
      rot2dtau(2,1) = - rot2dtau(1,2)
      rot2dtau(2,2) =   rot2dtau(1,1)

!----------------------------------------------

      call dag2d(rot2dt,trot2dt)
      call dag2d(rot2db,trot2db)
      call dag2d(rot2dtau,trot2dtau)

!---------------------------------------------
C     3rd family fermions
!---------------------------------------------

      costhtsq = (dcos(thetat))**2.d0
      sinthtsq = (dsin(thetat))**2.d0
      costhbsq = (dcos(thetab))**2.d0
      sinthbsq = (dsin(thetab))**2.d0
      costhtausq = (dcos(thetatau))**2.d0
      sinthtausq = (dsin(thetatau))**2.d0
      
      call a0(mt,q,a0mt)
      call a0(mtau,q,a0mtau)     
      call a0(mt1,q,a0mstop1)
      call a0(mt2,q,a0mstop2)
      call a0(mb1,q,a0msbot1)
      call a0(mb2,q,a0msbot2)
      call a0(mtau1,q,a0mstau1)
      call a0(mtau2,q,a0mstau2)
      
      call b0(p,mt,mt,q,b0mtmt)
      call b0(p,mtau,mtau,q,b0mtaumtau)      

!---------------------------------------------------------------------------
C     Corrections begin
!---------------------------------------------------------------------------

      sbots = 0.d0
      staus = 0.d0
      stops = 0.d0
      fermions = 0.d0

c$$$      print*,"MT", mt
c$$$      print*,"in pis2s2, ht = ", ht
c$$$      print*,"p = ", p, "q = ", q
c$$$      print*,"cossqthW = ", cossqthW

c$$$      print*,"in pis2s2, htsq = ", htsq
c$$$      print*,"B0mtmt, a0mt  = ", b0mtmt, a0mt  

      fermions = 3.d0 * htsq *((p*p - 4.d0 *mt*mt) * b0mtmt - 
     $     2.d0 * a0mt)
       
      stops = 3.d0 * htsq * (a0mstop1 + a0mstop2)
      
c$$$      print*,"a0mstop1", a0mstop1
c$$$      print*,"a0mstop2", a0mstop2
c$$$
c$$$      print*,"stops - step 1", stops/(16.d0*pi*pi)


      stops = stops - (3.d0 * g*g / (2.d0 * cossqthW)) *
     $     (guR *(costhtsq * a0mstop2 + sinthtsq * a0mstop1) +
     $     guL * (sinthtsq * a0mstop2 + costhtsq * a0mstop1))

c$$$      print*,"stops - step 2", stops/(16.d0*pi*pi)

      sbots = sbots - (3.d0 * g*g / (2.d0 * cossqthW)) *
     $     (gdR *(costhbsq * a0msbot2 + sinthbsq * a0msbot1) +
     $     gdL * (sinthbsq * a0msbot2 + costhbsq * a0msbot1))

      staus = staus - (1.d0 * g*g / (2.d0 * cossqthW)) *
     $     (geR *(costhtausq * a0mstau2 + sinthtausq * a0mstau1) +
     $     geL * (sinthtausq * a0mstau2 + costhtausq * a0mstau1))

      sups = 0.d0
      sdowns = 0.d0
      sneutrinos = 0.d0
      sleps = 0.d0

!-------------------------------------------
      call a0(muL,q,a0muL(1))
      call a0(muR,q,a0muR(1))
      call a0(mcL,q,a0muL(2))
      call a0(mcR,q,a0muR(2))
      call a0(mdL,q,a0mdL(1))
      call a0(mdR,q,a0mdR(1))
      call a0(msL,q,a0mdL(2))
      call a0(msR,q,a0mdR(2))
      call a0(meL,q,a0meL(1))
      call a0(meR,q,a0meR(1))
      call a0(mmuL,q,a0meL(2))
      call a0(mmuR,q,a0meR(2))
      call a0(msnu(1),q,a0msnu(1))
      call a0(msnu(2),q,a0msnu(2))
      call a0(msnu(3),q,a0msnu(3))
!-------------------------------------------

      loop1: do i = 1, 2

      sups = sups - 3.d0 * (g*g / (2.d0 * cossqthW)) *
     $     (guL * a0muL(i) + guR * a0muR(i))
      sdowns = sdowns - 3.d0 * (g*g / (2.d0 * cossqthW)) *
     $     (gdL * a0mdL(i) + gdR * a0mdR(i))
      sleps = sleps - (g*g / (2.d0 * cossqthW)) * 
     $     (geL * a0meL(i) + geR * a0meR(i))

      enddo loop1

!-------------------------------------------

      ls2tt(1, 1) = - (g * MZ / costhW) * guL * sinbeta + htsq * vev2 
      ls2tt(1, 2) =   (ht * AURG(3,3)) / dsqrt(2.d0) 
      ls2tt(2, 1) =   ls2tt(1, 2)
      ls2tt(2, 2) = - (g * MZ / costhW )* guR * sinbeta + htsq * vev2 

!---------------------------------------------

      ls2bb(1, 1) = - (g * MZ / costhW) * gdL * sinbeta 
      ls2bb(1, 2) = - 1.d0 * hb * sgnmu * modmu / dsqrt(2.d0)  
      ls2bb(2, 1) =   ls2bb(1, 2)
      ls2bb(2, 2) = - (g * MZ / costhW ) * gdR * sinbeta 

!---------------------------------------------

      ls2tautau(1, 1) = - (g * MZ / costhW) * geL * sinbeta 
      ls2tautau(1, 2) = - 1.d0 *  htau * sgnmu * modmu / dsqrt(2.d0)
      ls2tautau(2, 1) =   ls2tautau(1, 2)
      ls2tautau(2, 2) = - (g * MZ / costhW) * geR * sinbeta

!----------------------------------------------

      call mat3prod2d(rot2dt,ls2tt,trot2dt,lTS212)
      call mat3prod2d(rot2db,ls2bb,trot2db,lBS212)      
      call mat3prod2d(rot2dtau,ls2tautau,trot2dtau,lTauS212)      

!-----------------------------------------------      
      call b0(p,mt1,mt1,q,b0mt(1,1))
      call b0(p,mt1,mt2,q,b0mt(1,2))
      call b0(p,mt2,mt1,q,b0mt(2,1))
      call b0(p,mt2,mt2,q,b0mt(2,2))

      call b0(p,mb1,mb1,q,b0mb(1,1))
      call b0(p,mb1,mb2,q,b0mb(1,2))
      call b0(p,mb2,mb1,q,b0mb(2,1))
      call b0(p,mb2,mb2,q,b0mb(2,2))

      call b0(p,mtau1,mtau1,q,b0mtau(1,1))
      call b0(p,mtau1,mtau2,q,b0mtau(1,2))
      call b0(p,mtau2,mtau1,q,b0mtau(2,1))
      call b0(p,mtau2,mtau2,q,b0mtau(2,2))
      
      call b0(p,muL,muL,q,b0muL(1))      
      call b0(p,mcL,mcL,q,b0muL(2))      
      call b0(p,muR,muR,q,b0muR(1))      
      call b0(p,mcR,mcR,q,b0muR(2))      
      call b0(p,mdL,mdL,q,b0mdL(1))      
      call b0(p,msL,msL,q,b0mdL(2))      
      call b0(p,mdR,mdR,q,b0mdR(1))      
      call b0(p,msR,msR,q,b0mdR(2))      
      call b0(p,meL,meL,q,b0meL(1))      
      call b0(p,mmuL,mmuL,q,b0meL(2))      
      call b0(p,meR,meR,q,b0meR(1))      
      call b0(p,mmuR,mmuR,q,b0meR(2))      

      call b0(p,msnu(1),msnu(1),q,b0msnu(1))      
      call b0(p,msnu(2),msnu(2),q,b0msnu(2))      
      call b0(p,msnu(3),msnu(3),q,b0msnu(3))
      
!--------------------------------------------------

      loop2: do i = 1, 2
      loop3: do j = 1, 2

      stops = stops + 3.d0 * lTS212(i, j) * lTS212(i, j) * b0mt(i,j)
      sbots = sbots + 3.d0 * lBS212(i, j) * lBS212(i, j) * b0mb(i,j)
      staus = staus + lTauS212(i, j) * lTauS212(i, j) * b0mtau(i,j)

      enddo loop3
      enddo loop2

c$$$      print*,"stops - step 3", stops/(16.d0*pi*pi)

!---------------------------------------------------------------------------

      ls2eeLL = -(g * MZ * geL * sinbeta) / costhW
      ls2eeRR = -(g * MZ * geR * sinbeta) / costhW
      ls2ddLL = -(g * MZ * gdL * sinbeta) / costhW
      ls2ddRR = -(g * MZ * gdR * sinbeta) / costhW
      ls2uuLL = -(g * MZ * guL * sinbeta) / costhW
      ls2uuRR = -(g * MZ * guR * sinbeta) / costhW
      
      lSnuS2LL =  (g*MZ*gnuL*sinbeta/costhW)**2.d0


      loop4: do i = 1, 2

      sups = sups + 3.d0 * ls2uuLL * ls2uuLL * b0muL(i) + 
     $     3.d0 * ls2uuRR * ls2uuRR * b0muR(i)

      sdowns = sdowns + 3.d0 * ls2ddLL * ls2ddLL * b0mdL(i) 
     $     + 3.d0 * ls2ddRR * ls2ddRR * b0mdR(i)

      sleps = sleps + ls2eeLL * ls2eeLL * b0meL(i) + 
     $     ls2eeRR * ls2eeRR * b0meR(i)
      
      enddo loop4

      sneutrinos = sneutrinos - (g*g / (2.d0 * cossqthW)) * gnuL * 
     $     (a0msnu(1) + a0msnu(2) + a0msnu(3)) +
     $     lSnuS2LL * (b0msnu(1) + b0msnu(2) + b0msnu(3))

!------------------------------------------------------------------------------
       
      call f(p,mHpm,MW,q,fmHpmMW)
      call f(p,mA0,MZ,q,fmAMZ)
      call f(p,MW,MW,q,fMWMW)
      call f(p,MZ,MZ,q,fMZMZ)
      call b0(p,MW,MW,q,b0MWMW)      
      call b0(p,MZ,MZ,q,b0MZMZ)      
      call a0(MW,q,a0MW)
      call a0(MZ,q,a0MZ)

!-------------------------------------------

c$$$      print*,"fmhphmMW", fmhpmMW,p,mHpm,q
c$$$      print*,"fmA0MZ", fmAMZ,mA0
c$$$      print*,"fMWMW", fMWMW,MW
c$$$      print*,"fMZMZ", fMZMZ,MZ

!     print*,"fmhphmMW", fmhpmMW

      
      cossqbeta = cosbeta*cosbeta
      sinsqbeta = sinbeta*sinbeta
      
      higgs = 0.d0
      
      higgs = g *g * 0.25d0 * (cossqbeta * (2.d0 * fmHpmMW + 
     $     (fmAMZ/ cossqthW)) + sinsqbeta * (2.d0 * fMwMW  + 
     $     (fMZMZ / cossqthW))) + 1.75d0 * g * g * sinsqbeta * 
     $     (2.d0 * MW * MW * b0MWMW + 
     $     (MZ*MZ*b0MZMZ / cossqthW)) +
     $     g*g * (2.d0 * a0MW + (a0MZ / cossqthW))

c$$$      print*,"b0MWMW = ", b0MWMW
c$$$      print*,"b0MZMZ = ", b0MZMZ
c$$$      print*,"a0MZ = ", a0MZ
c$$$      print*,"a0MW = ", a0MW

c$$$
c$$$      print*,"g, cosb2,sinb2 = ", g, cossqbeta, sinsqbeta
c$$$      print*,"costhwsq = ", cossqthW
c$$$
c$$$      print*,"higgs1 in pis2s2", higgs !/ (16.d0 * pi *pi)

!------------------------------------------------------------------------------
      
      loophhs2i: do i = 1, 4
      loophhs2j: do j = 1, 4
      
      hhs2(i,j)   = 0.d0

      enddo loophhs2j

      hhs2s2(i) = 0.d0

      enddo loophhs2i

      hhs2(1, 1) =  sinbeta * ((3.d0 * dsin(alpha)**2.d0) - 
     $      (dcos(alpha)**2.d0)) - cosbeta * dsin(2.d0*alpha)

      hhs2(2, 2) =  sinbeta * ((3.d0 * dcos(alpha)**2.d0) -
     $     (dsin(alpha)**2.d0)) + cosbeta * dsin(2.d0*alpha)

      hhs2(1, 2) =  2.d0 * sinbeta * dsin(2.d0*alpha) - 
     $     (cosbeta * dcos(2.d0*alpha))

      hhs2(2, 1) =    hhs2(1, 2)

      hhs2(3, 3) =  - cos2beta * sinbeta
      hhs2(3, 4) =    sin2beta * sinbeta
      hhs2(4, 3) =    hhs2(3, 4)
      hhs2(4, 4) =    cos2beta * sinbeta

      loop5: do i = 1, 4
      loop6: do j = 1, 4

      hhs2(i,j) = hhs2(i,j) * (g * MZ / (2.d0 * costhW))
      
      enddo loop6
      enddo loop5

!--------------------------------------------------

      hhs2s2(1) = (3.d0 * dsin(alpha)**2.d0) - (dcos(alpha)**2.d0)
      hhs2s2(2) = (3.d0 * dcos(alpha)**2.d0) - (dsin(alpha)**2.d0)
      hhs2s2(3) = - cos2beta
      hhs2s2(4) =   cos2beta
      
      loop7: do i = 1, 4
      
      hhs2s2(i) = hhs2s2(i) * ((g*g) * 0.25d0 / cossqthW)

c$$$      print*,"hhs2s2(",i,")", hhs2s2(i)

      enddo loop7

!-------------------------------------------------- 
c$$$      print*,"alpha = ", alpha
c$$$      print*,"mHu0 = ", mHu0
c$$$      print*,"mh0 = ", mh0
c$$$      print*,"mG0 = ", mG0
c$$$      print*,"mA0 = ", mA0

      call b0(p,mHu0,mHu0,q,b0HH(1,1))      
      call b0(p,mHu0,mh0,q,b0HH(1,2))      
      call b0(p,mHu0,mG0,q,b0HH(1,3))      
      call b0(p,mHu0,mA0,q,b0HH(1,4))      

      call b0(p,mh0,mHu0,q,b0HH(2,1))      
      call b0(p,mh0,mh0,q,b0HH(2,2))      
      call b0(p,mh0,mG0,q,b0HH(2,3))      
      call b0(p,mh0,mA0,q,b0HH(2,4))      

      call b0(p,mG0,mHu0,q,b0HH(3,1))      
      call b0(p,mG0,mh0,q,b0HH(3,2))      
      call b0(p,mG0,mG0,q,b0HH(3,3))      
      call b0(p,mG0,mA0,q,b0HH(3,4))      

      call b0(p,mA0,mHu0,q,b0HH(4,1))      
      call b0(p,mA0,mh0,q,b0HH(4,2))      
      call b0(p,mA0,mG0,q,b0HH(4,3))      
      call b0(p,mA0,mA0,q,b0HH(4,4)) 
      
      call a0(mHu0,q,a0H(1))
      call a0(mh0,q,a0H(2))
      call a0(mG0,q,a0H(3))
      call a0(mA0,q,a0H(4))

c$$$      print*,"mh0", mh0
c$$$      print*,"mHu0", mHu0
c$$$      print*,"mG0", mG0
c$$$      print*,"mA0", mA0
c$$$
c$$$
c$$$      print*,"a0H", a0H(1)
c$$$      print*,"aamh", a0H(2)
c$$$      print*,"a0mG0", a0H(3)
c$$$      print*,"a0mA0", a0H(4)

!---------------------------------------------------      
      
      checkhiggs1 = 0.d0
      checkhiggs2 = 0.d0

      loop8: do i = 1, 4
      loop9: do j = 1, 4
      
      higgs = higgs + 0.5d0 * hhs2(i, j)* hhs2(i, j) * b0HH(i,j)
      checkhiggs1 = checkhiggs1 + 0.5d0 * hhs2(i, j)* hhs2(i, j) * 
     $     b0HH(i,j)

      enddo loop9

      higgs = higgs + 0.5d0 * hhs2s2(i) * a0H(i)
      checkhiggs2 = checkhiggs2 + 0.5d0 * hhs2s2(i) * a0H(i)

      enddo loop8


c$$$      print*,"higgs 1b = ", checkhiggs1  / (16.d0 * pi *pi)
c$$$      print*,"higgs 1c = ", checkhiggs2 / (16.d0 * pi *pi)
c$$$!      print*,"higgs 2 = ", higgs !/ (16.d0 * pi *pi)

!------------------------------------------------------------------------------

      hphps2(1, 1) =  - cos2beta * sinbeta
      hphps2(2, 2) =    cos2beta * sinbeta + 2.d0 * cossqthW * sinbeta
      hphps2(1, 2) =    sin2beta * sinbeta - cossqthW * cosbeta 
      hphps2(2, 1) =    hphps2(1, 2)
      
      loop10: do i = 1, 2
      loop11: do j = 1, 2
      
      hphps2(i,j) = hphps2(i,j) * (g * MZ * 0.5d0 / costhW)
c$$$      print*,"hphps2(",i,",",j,")", hphps2(i,j)
      enddo loop11
      enddo loop10

!------------------------------------------------      
 
      hphps2s2(1) = cossqthW - sinsqthW * cos2beta
      hphps2s2(2) = cossqthW + sinsqthW * cos2beta
      
      loop12: do i = 1, 2
      hphps2s2(i) = hphps2s2(i) * (g*g* 0.25d0 / cossqthW)
c$$$      print*,"hphps2s2(",i,")", hphps2s2(i)

      enddo loop12

!-------------------------------------------
      call b0(p,mGp,mGp,q,b0HpHp(1,1))      
      call b0(p,mGp,mHpm,q,b0HpHp(1,2))      
      call b0(p,mHpm,mGp,q,b0HpHp(2,1))      
      call b0(p,mHpm,mHpm,q,b0HpHp(2,2))      
      
      call a0(mGp,q,a0Hp(1))
      call a0(mHpm,q,a0Hp(2))

c$$$      print*,"mGp = ", mGp
c$$$      print*,"mHpm = ", mHpm

c$$$      print*,"a0HP(1)", a0Hp(1)
c$$$      print*,"a0HP(2)", a0Hp(2)

!----------------------------------------------
            
      checkhiggs3 = 0.d0
      checkhiggs4 = 0.d0

      loop13: do i = 1, 2
      loop14: do j = 1, 2
      
      higgs = higgs + hphps2(i, j) * hphps2(i, j)* b0HpHp(i,j)
      checkhiggs3 = checkhiggs3 + hphps2(i, j) * hphps2(i, j)* 
     $     b0HpHp(i,j)

      enddo loop14

      higgs = higgs + hphps2s2(i) * a0Hp(i)
      checkhiggs4 = checkhiggs4 + hphps2s2(i) * a0Hp(i)

      enddo loop13

c$$$      print*,"check higgs 3", checkhiggs3 / (16.d0 * pi *pi)
c$$$      print*,"check higgs 4", checkhiggs4 / (16.d0 * pi *pi)
c$$$!      print*,"check higgs 3 ", higgs !/ (16.d0 * pi *pi)
c$$$      print*,"quartic +  trilinear", (checkhiggs1 + checkhiggs2 + 
c$$$     $     checkhiggs3 + checkhiggs4) 
c$$$     $     / (16.d0 * pi *pi)

!----------------------------------------------------------------------------
      
      aPsi(1, 4) =   gp * 0.5d0
      aPsi(4, 1) =   aPsi(1, 4)
      aPsi(2, 4) = - g  * 0.5d0
      aPsi(4, 2) =   aPsi(2, 4)
    
      call dag4d(ON,ONdag)
      call mat3prod4d(ON,aPsi,ONdag,aChi)
      call mat3prod4d(ON,aPsi,ONdag,bChi)  

      loop17: do i = 1, 4
      loop181: do j = 1, 4
      
      call funcg(p,dabs(mneut(i)),dabs(mneut(j)),q,gmneut(i,j))
      call b0(p,dabs(mneut(i)),dabs(mneut(j)),q,b0mneut(i,j))
      
      enddo loop181
      enddo loop17


      neutralinos = 0.d0

      loop15: do i = 1, 4
      loop16: do j = 1, 4
      
      fChiChis2s2(i, j) = real(aChi(i, j)*(aChi(i,j)) + 
     $     bChi(i, j)*(bChi(i,j)))

      gChiChis2s2(i, j) = real(bChi(i, j) * aChi(i, j) + 
     $     (aChi(i, j)) * bChi(i, j))

      neutralinos = neutralinos + 0.5d0 * 
     $     ((fChiChis2s2(i, j) * gmneut(i,j)) - 2.d0 *
     $     (gChiChis2s2(i, j) * mneut(i) * mneut(j) * b0mneut(i,j)))

      enddo loop16
      enddo loop15

!------------------------------------------------------------------------------

      aPsic(2, 1) = g / dsqrt(2.d0)
      
      call dag2d(aPsic,aPsicdag)
      call dag2d(OCL,OCLdag)
      call dag2d(OCR,OCRdag)
      call mat3prod2d(OCR,aPsic,OCLdag,aChic) !<--------------------CHECKED!!
      call mat3prod2d(OCL,aPsicdag,OCRdag,bChic)

      loop19: do i = 1, 2
      loop18: do j = 1, 2
      
      call funcg(p,mchargino(i),mchargino(j),q,gmch(i,j))
      call b0(p,mchargino(i),mchargino(j),q,b0mch(i,j))
      
      enddo loop18
      enddo loop19

      charginos = 0.d0
      
      
      loop23: do i = 1, 2
      loop24: do j = 1, 2

      fChiChis2s2(i, j) = real(aChic(i, j) * (aChic(i,j)) +
     $     bChic(i, j) * (bChic(i,j)))

      gChiChis2s2(i, j) = real(bChic(i, j) * aChic(i, j) + 
     $     (aChic(i, j)) * bChic(i, j))

      charginos = charginos + 
     $     (fChiChis2s2(i, j) * gmch(i,j) - 2.d0 *
     $     gChiChis2s2(i, j) * mchargino(i) * mchargino(j) * b0mch(i,j))

      enddo loop24
      enddo loop23

!-------------------------------------------------------------------------------

      pis2s2ans = (sups  + sdowns  + sleps  + sbots  + staus +  stops +
     $     sneutrinos + fermions + higgs + neutralinos + charginos) 
     $     / (16.d0*pi*pi)
      
c$$$      print*,"in pis2s2"
c$$$      print*,"sups = ", sups    ! /(16.d0*pi*pi) 
c$$$      print*,"sdowns = ", sdowns !/(16.d0*pi*pi)
c$$$      print*,"sleps = ", sleps  !/(16.d0*pi*pi)
c$$$      print*,"stops = ", stops  !/(16.d0*pi*pi)
c$$$      print*,"sbots = ", sbots  !/(16.d0*pi*pi)
c$$$      print*,"staus = ", staus  !/(16.d0*pi*pi)
c$$$      print*,"sneutrinos = ", sneutrinos !/(16.d0*pi*pi)
c$$$      print*,"fermions = ", fermions !/(16.d0*pi*pi)
c$$$      print*,"higgs = ", higgs  !/(16.d0*pi*pi)
c$$$      print*,"neutralinos = ", neutralinos !/(16.d0*pi*pi)
c$$$      print*,"charginos = ", charginos !/(16.d0*pi*pi)

!-----------------------------------------------------------------------------
      RETURN      
      end subroutine pis2s2

!-----------------------------------------------------------------------------
!=============================================================================

!=============================================================================
!-----------------------------------------------------------------------------

      subroutine pis1s2(p,q,g,gp,mt,mb,mtau,tanbeta,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,ON,OCL,
     $     OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,mA0sq,vev1,vev2,
     $     pis1s2ans)

      implicit none
      integer i,j
      double precision p,q,g,gp,pis1s2ans,stops,sbots,staus
      double precision sneutrinos
      double precision sleps,sups,sdowns,fermions,higgs,charginos
      double precision neutralinos

      double precision mT, mB, mTau, tanbeta, sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      double precision rot2dt(2,2),trot2dt(2,2),rot2db(2,2),trot2db(2,2)
      double precision rot2dtau(2,2),trot2dtau(2,2)

      double precision ls2tt(2,2),ls2bb(2,2)
      double precision ls2tautau(2,2)
      double precision ls1tt(2,2),ls1bb(2,2)
      double precision ls1tautau(2,2)
      double precision lTS212(2,2),lBS212(2,2),lTauS212(2,2)
      double precision lTS112(2,2),lBS112(2,2),lTauS112(2,2)
      double precision lSnuS1LL,lSnuS2LL
      
      double precision a0mtau,a0mt

      double precision ON(4,4),OCL(2,2),OCR(2,2),ONdag(4,4),OCLdag(2,2)
      double precision OCRdag(2,2)
      double precision ht,htsq,hb,htau,hbsq,htausq
      double precision modmu,mh0,mHpm,mhu0,mA0


      double precision b0mtaumtau,a0mstop1,a0mstop2,b0mtmt
      double precision a0msbot1,a0msbot2,a0mstau1,a0mstau2,a0msnu(3)

      double precision costhtsq,sinthtsq,costhbsq,sinthbsq,costhtausq
      double precision sinthtausq,cossqbeta,sinsqbeta,sinbeta

      double precision a0muL(2),a0muR(2),a0mdL(2),a0mdR(2)
      double precision a0meL(2),a0meR(2),b0msnu(3)


      double precision b0mt(2,2),b0mb(2,2),b0mtau(2,2),b0muL(2),b0muR(2)
      double precision b0mdL(2),b0mdR(2),b0meL(2),b0meR(2)

      double precision fmHpmMW,fmAMZ,fMWMW,fMZMZ
      double precision b0MWMW,b0MZMZ !,a0MW,a0MZ 


      double precision ls2eeLL,ls2eeRR,ls2ddLL,ls2ddRR,ls2uuLL,ls2uuRR
      double precision ls1eeLL,ls1eeRR,ls1ddLL,ls1ddRR,ls1uuLL,ls1uuRR

      double precision hhs2(4,4),hphps2(2,2)
      data hhs2/ 16 * 0.d0/,hphps2/ 4 * 0.d0/

      double precision hhs1(4,4),hhs1s2(4),hphps1(2,2),hphps1s2(2)
      data hhs1/ 16 * 0.d0/,hhs1s2/ 4 * 0.d0/,hphps1/ 4 * 0.d0/
      data hphps1s2/ 2 * 0.d0/

      double precision b0HH(4,4),a0H(4),mG0,mGp,b0HpHp(2,2),a0Hp(2)

      double precision fChiChis1s2(4, 4),gChiChis1s2(4, 4),gmneut(4,4)
      double precision b0mneut(4,4)

      real*8 aPsi1(4,4),aChi1(4,4),bChi1(4,4)
      data aPsi1/ 16 * 0.d0/, aChi1/ 16 * 0.d0/, bChi1/ 16 * 0.d0/ 

      real*8 aPsi2(4,4),aChi2(4,4),bChi2(4,4)
      data aPsi2/ 16 * 0.d0/, aChi2/ 16 * 0.d0/, bChi2/ 16 * 0.d0/ 

      double precision aPsic1(2, 2),aPsic1dag(2,2),gmch(2,2),b0mch(2,2)
      data aPsic1/ 4 * 0.d0/,aPsic1dag/ 4 * 0.d0/, gmch/ 4 * 0.d0/
  
      double precision aPsic2(2, 2),aPsic2dag(2,2)
      data aPsic2/ 4 * 0.d0/,aPsic2dag/ 4 * 0.d0/

      DOUBLE PRECISION aChic1(2, 2),bChic1(2, 2),aChic2(2, 2)
      data aChic1/ 4 * 0.d0/,bChic1/ 4 * 0.d0/, aChic2/ 4 * 0.d0/

      DOUBLE PRECISION bChic2(2, 2)
      data bChic2/ 4 * 0.d0/

      DOUBLE PRECISION AERG(3,3),AURG(3,3),ADRG(3,3)
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq

      DOUBLE PRECISION muL,muR,mcL,mcR,mt1,mt2
      DOUBLE PRECISION mdL,mdR,msL,msR,mb1,mb2
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtau1,mtau2      
      DOUBLE PRECISION mneut(4),msnu(3),mchargino(2)

      double precision cossqthw,sinsqthw,cosbeta,costhw,sinthw
      double precision alpha,thw,tan2beta
      double precision thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae,mtL,mtR
      double precision yuRg(3,3),ydRG(3,3),yeRg(3,3)

      double precision beta
      double precision vev1,vev2,cos2beta,sin2beta

      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole


!-------------------------------------------------------------------------

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

!-------------------------------------------------------------------------
      external dag2d,dag4d,mat3prod2d,a0,theta,b0,f,funcg

      include 'stdinputs.h'

!-------------------------------------------------------------------------

      pi = 4.d0 * datan(1.d0)

      MW = MWpole
      MZ = MZpole

      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mt1 = dsqrt(SUegg(2))
      mt2 = dsqrt(SUegg(1))

      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mb1 = dsqrt(SDegg(2))
      mb2 = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtau1 = dsqrt(SLegg(2))
      mtau2 = dsqrt(SLegg(1))
       
      msnu(1) = dsqrt(SNegg(1))
      msnu(2) = dsqrt(SNegg(2))
      msnu(3) = dsqrt(SNegg(3))

      mneut(1) = Neg(1) 
      mneut(2) = Neg(2)
      mneut(3) = Neg(3)
      mneut(4) = Neg(4)

      mchargino(1) = Ceg(1)
      mchargino(2) = Ceg(2)

      mh0  = dsqrt((mh0sq))
      mHu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      mG0 = MZ
      mGp = MW

!----------------------------------------------
C     Definitions.
!----------------------------------------------
      
      beta     = datan(tanbeta)
      cosbeta  = dcos(beta)
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0*beta)
      cos2beta = dcos(2.d0*beta)
      tan2beta = dtan(2.d0*beta)

      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/
     $     (mA0sq - MZ*MZ))*tan2beta)
      
      sinsqthw = sinsqthw_susy !(1.d0 - (MW/MZ)**2.d0)
      sinthw = dsqrt(sinsqthw)
      thw = dasin(dsqrt(sinsqthw))
      costhw = dcos(thw)
      cossqthw = dcos(thw)*dcos(thw)

      gnuL =   0.5d0
      guL  =   0.5d0 - 2.d0*sinsqthw/ 3.d0 
      gdL  =  -0.5d0 + sinsqthw/3.d0
      geL  =  -0.5d0 + sinsqthw
      guR  =   2.d0*sinsqthw/3.d0 
      gdR  =  -sinsqthw/3.d0
      geR  =  -sinsqthw 
      yuL  =   1.d0/3.d0 
      yuR  =  -4.d0/3.d0 
      ydL  =   1.d0/3.d0 
      ydR  =   2.d0/3.d0 
      yeL  =  -1.d0 
      yeR  =   2.d0 
      ynuL =  -1.d0 

      ht   =  yuRG(3,3) 
      hb   =  ydRg(3,3) 
      htau =  yeRG(3,3) 

      htsq =   ht*ht
      hbsq =   hb*hb
      htausq = htau*htau

!----------------------------------------------------------------
      
      rot2dt(1,1) =   dcos(thetat)
      rot2dt(1,2) =   dsin(thetat)
      rot2dt(2,1) = - rot2dt(1,2)
      rot2dt(2,2) =   rot2dt(1,1)

      rot2db(1,1) =   dcos(thetab)
      rot2db(1,2) =   dsin(thetab)
      rot2db(2,1) = - rot2db(1,2)
      rot2db(2,2) =   rot2db(1,1)

      rot2dtau(1,1) =   dcos(thetatau)
      rot2dtau(1,2) =   dsin(thetatau)
      rot2dtau(2,1) = - rot2dtau(1,2)
      rot2dtau(2,2) =   rot2dtau(1,1)

!-------------------------------------------------
      
      call dag2d(rot2dt,trot2dt)
      call dag2d(rot2db,trot2db)
      call dag2d(rot2dtau,trot2dtau)

!-------------------------------------------------
C     3rd family fermions
!-------------------------------------------------

      costhtsq = (dcos(thetat))**2.d0
      sinthtsq = (dsin(thetat))**2.d0
      costhbsq = (dcos(thetab))**2.d0
      sinthbsq = (dsin(thetab))**2.d0
      costhtausq = (dcos(thetatau))**2.d0
      sinthtausq = (dsin(thetatau))**2.d0

     
      call a0(mt,q,a0mt)
      call b0(p,mt,mt,q,b0mtmt)
      call a0(mtau,q,a0mtau)          
      call a0(mt1,q,a0mstop1)
      call a0(mt2,q,a0mstop2)
      call a0(mb1,q,a0msbot1)
      call a0(mb2,q,a0msbot2)
      call a0(mtau1,q,a0mstau1)
      call a0(mtau2,q,a0mstau2)
      call b0(p,mtau,mtau,q,b0mtaumtau)      

!----------------------------------------------------------------
      
      fermions = 0.d0
      sbots = 0.d0
      staus = 0.d0
      sups = 0.d0
      sdowns = 0.d0
      sneutrinos = 0.d0
      sleps = 0.d0

      call a0(muL,q,a0muL(1))
      call a0(muR,q,a0muR(1))
      call a0(mcL,q,a0muL(2))
      call a0(mcR,q,a0muR(2))
      call a0(mdL,q,a0mdL(1))
      call a0(mdR,q,a0mdR(1))
      call a0(msL,q,a0mdL(2))
      call a0(msR,q,a0mdR(2))
      call a0(meL,q,a0meL(1))
      call a0(meR,q,a0meR(1))
      call a0(mmuL,q,a0meL(2))
      call a0(mmuR,q,a0meR(2))
      call a0(msnu(1),q,a0msnu(1))
      call a0(msnu(2),q,a0msnu(2))
      call a0(msnu(3),q,a0msnu(3))

!----------------------------------------------------------------------------

      ls1tt(1, 1) = g * MZ * guL * cosbeta / costhw
      ls1tt(1, 2) = -1.d0 * ht * sgnmu * modmu / dsqrt(2.d0)
      ls1tt(2, 1) = ls1tt(1, 2)
      ls1tt(2, 2) = g * MZ * guR * cosbeta / costhW

!-----------------------------------------------

      ls1bb(1, 1) = (g * MZ * gdL * cosbeta / costhW) + hbsq*vev1
      ls1bb(1, 2) = hb*ADRG(3,3)/dsqrt(2.d0) 
      ls1bb(2, 1) = ls1bb(1, 2)
      ls1bb(2, 2) = (g * MZ * gdR * cosbeta / costhW) + hbsq*vev1

!---------------------------------------------

      ls1tautau(1, 1) = (g * MZ * geL * cosbeta / costhW) + htausq*vev1
      ls1tautau(1, 2) = htau*AERG(3,3) / dsqrt(2.d0) 
      ls1tautau(2, 1) = ls1tautau(1, 2)
      ls1tautau(2, 2) = (g * MZ * geR * cosbeta / costhW) + htausq*vev1

!----------------------------------------------------------------------------

      call mat3prod2d(rot2dt,ls1tt,trot2dt,lTS112)
      call mat3prod2d(rot2db,ls1bb,trot2db,lBS112)      
      call mat3prod2d(rot2dtau,ls1tautau,trot2dtau,lTauS112)      

!----------------------------------------------------------------------------

      ls2tt(1, 1) = - (g * MZ / costhW) * guL * sinbeta + htsq * vev2 
      ls2tt(1, 2) =   (ht * AURG(3,3)) / dsqrt(2.d0) 
      ls2tt(2, 1) =    ls2tt(1, 2)
      ls2tt(2, 2) = - (g * MZ / costhW )* guR * sinbeta + htsq * vev2 

!------------------------------------------------------------------------

      ls2bb(1, 1) = - (g * MZ / costhW) * gdL * sinbeta 
      ls2bb(1, 2) =  -1.d0 * hb * sgnmu * modmu / dsqrt(2.d0) 
      ls2bb(2, 1) =   ls2bb(1, 2)
      ls2bb(2, 2) = - (g * MZ / costhW ) * gdR * sinbeta 

!-------------------------------------------------------------------------

      ls2tautau(1, 1) = - (g * MZ / costhW) * geL * sinbeta 
      ls2tautau(1, 2) = - 1.d0 * htau * sgnmu * modmu / dsqrt(2.d0) 
      ls2tautau(2, 1) =   ls2tautau(1, 2)
      ls2tautau(2, 2) = - (g * MZ / costhW) * geR * sinbeta 

!------------------------------------------------------------------------      

      call mat3prod2d(rot2dt,ls2tt,trot2dt,lTS212)
      call mat3prod2d(rot2db,ls2bb,trot2db,lBS212)      
      call mat3prod2d(rot2dtau,ls2tautau,trot2dtau,lTauS212)      

!-----------------------------------------
C     PV functions 
!-----------------------------------------

      call b0(p,mt1,mt1,q,b0mt(1,1))
      call b0(p,mt1,mt2,q,b0mt(1,2))
      call b0(p,mt2,mt1,q,b0mt(2,1))
      call b0(p,mt2,mt2,q,b0mt(2,2))

      call b0(p,mb1,mb1,q,b0mb(1,1))
      call b0(p,mb1,mb2,q,b0mb(1,2))
      call b0(p,mb2,mb1,q,b0mb(2,1))
      call b0(p,mb2,mb2,q,b0mb(2,2))

      call b0(p,mtau1,mtau1,q,b0mtau(1,1))
      call b0(p,mtau1,mtau2,q,b0mtau(1,2))
      call b0(p,mtau2,mtau1,q,b0mtau(2,1))
      call b0(p,mtau2,mtau2,q,b0mtau(2,2))
      
      call b0(p,muL,muL,q,b0muL(1))      
      call b0(p,mcL,mcL,q,b0muL(2))      
      call b0(p,muR,muR,q,b0muR(1))      
      call b0(p,mcR,mcR,q,b0muR(2))      
      call b0(p,mdL,mdL,q,b0mdL(1))      
      call b0(p,msL,msL,q,b0mdL(2))      
      call b0(p,mdR,mdR,q,b0mdR(1))      
      call b0(p,msR,msR,q,b0mdR(2))      
      call b0(p,meL,meL,q,b0meL(1))      
      call b0(p,mmuL,mmuL,q,b0meL(2))      
      call b0(p,meR,meR,q,b0meR(1))      
      call b0(p,mmuR,mmuR,q,b0meR(2))      

      call b0(p,msnu(1),msnu(1),q,b0msnu(1))      
      call b0(p,msnu(2),msnu(2),q,b0msnu(2))      
      call b0(p,msnu(3),msnu(3),q,b0msnu(3)) 

      stops = 0.d0
      sbots = 0.d0
      staus = 0.d0
     
!--------------------------------------------------------
C     Corrections begin
!--------------------------------------------------------

      loop2: do i = 1, 2
      loop3: do j = 1, 2

      stops = stops + 3.d0 * lTS112(i, j) * lTS212(i, j) * b0mt(i,j)
      sbots = sbots + 3.d0 * lBS112(i, j) * lBS212(i, j) * b0mb(i,j)
      staus = staus + lTauS112(i, j) * lTauS212(i, j) * b0mtau(i,j)

      enddo loop3
      enddo loop2

!----------------------------------------------------------------------------

      ls1eeLL = (g * MZ * geL * cosbeta) / costhW
      ls1eeRR = (g * MZ * geR * cosbeta) / costhW
      ls1ddLL = (g * MZ * gdL * cosbeta) / costhW
      ls1ddRR = (g * MZ * gdR * cosbeta) / costhW
      ls1uuLL = (g * MZ * guL * cosbeta) / costhW
      ls1uuRR = (g * MZ * guR * cosbeta) / costhW


      ls2eeLL = -(g * MZ * geL * sinbeta) / costhW
      ls2eeRR = -(g * MZ * geR * sinbeta) / costhW
      ls2ddLL = -(g * MZ * gdL * sinbeta) / costhW
      ls2ddRR = -(g * MZ * gdR * sinbeta) / costhW
      ls2uuLL = -(g * MZ * guL * sinbeta) / costhW
      ls2uuRR = -(g * MZ * guR * sinbeta) / costhW
   
      lSnuS1LL =   (g*MZ*gnuL*cosbeta/costhW)
      lSnuS2LL = - (g*MZ*gnuL*sinbeta/costhW)

!-----------------------------------------------
      loop4: do i = 1, 2

      sups = sups + 3.d0 * ls1uuLL * ls2uuLL * b0muL(i) + 
     $     3.d0 * ls1uuRR * ls2uuRR * b0muR(i)

      sdowns = sdowns + 3.d0 * ls1ddLL * ls2ddLL * b0mdL(i) 
     $     + 3.d0 * ls1ddRR * ls2ddRR * b0mdR(i)

      sleps = sleps + ls1eeLL * ls2eeLL * b0meL(i) + 
     $     ls1eeRR * ls2eeRR * b0meR(i)
      
      enddo loop4

      sneutrinos = sneutrinos +
     $     lSnuS1LL * lSnuS2LL * (b0msnu(1) + b0msnu(2) + b0msnu(3))

!-----------------------------------------------------------------------

      call f(p,mHpm,MW,q,fmHpmMW)
      call f(p,mA0,MZ,q,fmAMZ)
      call f(p,MW,MW,q,fMWMW)
      call f(p,MZ,MZ,q,fMZMZ)
      call b0(p,MW,MW,q,b0MWMW)      
      call b0(p,MZ,MZ,q,b0MZMZ)      
      
      cossqbeta = cosbeta*cosbeta
      sinsqbeta = sinbeta*sinbeta
      
      higgs = 0.d0
      
      higgs = g *g * 0.25d0 * sinbeta*cosbeta *((2.d0 * fMwMW -
     $     2.d0 * fmHpmMW) + ((fMZMZ  - fmAMZ)/ cossqthW) + 
     $     7.d0 * (2.d0 * MW*MW * b0MWMW + (MZ*MZ*b0MZMZ / cossqthW)))
      

!-----------------------------------------------------------------------------

      hhs1(1, 1) =  cosbeta * ((3.d0 * cos(alpha)**2.d0) - 
     $      (sin(alpha)**2.d0)) - sinbeta * sin(2.d0*alpha)

      hhs1(2, 2) =  cosbeta * ((3.d0 * sin(alpha)**2.d0) -
     $     (cos(alpha)**2.d0)) + sinbeta * sin(2.d0*alpha)

      hhs1(1, 2) = - 2.d0 * cosbeta * sin(2.d0*alpha) - 
     $     (sinbeta * cos(2.d0*alpha))

      hhs1(2, 1) =  hhs1(1, 2)

      hhs1(3, 3) =  cos2beta * cosbeta
      hhs1(3, 4) = -sin2beta * cosbeta
      hhs1(4, 3) =  hhs1(3, 4)
      hhs1(4, 4) = -cos2beta * cosbeta

      hhs2(1, 1) =  sinbeta * ((3.d0 * sin(alpha)**2.d0) - 
     $      (cos(alpha)**2.d0)) - cosbeta * sin(2.d0*alpha)

      hhs2(2, 2) =  sinbeta * ((3.d0 * cos(alpha)**2.d0) -
     $     (sin(alpha)**2.d0)) + cosbeta * sin(2.d0*alpha)

      hhs2(1, 2) =  2.d0 * sinbeta * sin(2.d0*alpha) - 
     $     (cosbeta * cos(2.d0*alpha))

      hhs2(2, 1) =  hhs2(1, 2)

      hhs2(3, 3) = - cos2beta * sinbeta
      hhs2(3, 4) =  sin2beta * sinbeta
      hhs2(4, 3) =  hhs2(3, 4)
      hhs2(4, 4) =  cos2beta * sinbeta

      loop5: do i = 1, 4
      loop6: do j = 1, 4

      hhs2(i,j) = hhs2(i,j) * (g * MZ / (2.d0 * costhW))
      hhs1(i,j) = hhs1(i,j) * (g * MZ / (2.d0 * costhW))
      
      enddo loop6
      enddo loop5

!-----------------------------------------

      hhs1s2(1) = - sin(2.d0*alpha)
      hhs1s2(2) =   sin(2.d0*alpha)
      hhs1s2(3) = 0.d0
      hhs1s2(4) = 0.d0
      
      loop7: do i = 1, 4
      
      hhs1s2(i) = hhs1s2(i) * ((g*g) * 0.25d0 / cossqthW)

      enddo loop7

!-----------------------------------------
      
      call b0(p,mHu0,mHu0,q,b0HH(1,1))      
      call b0(p,mHu0,mh0,q,b0HH(1,2))      
      call b0(p,mHu0,mG0,q,b0HH(1,3))      
      call b0(p,mHu0,mA0,q,b0HH(1,4))      

      call b0(p,mh0,mHu0,q,b0HH(2,1))      
      call b0(p,mh0,mh0,q,b0HH(2,2))      
      call b0(p,mh0,mG0,q,b0HH(2,3))      
      call b0(p,mh0,mA0,q,b0HH(2,4))      

      call b0(p,mG0,mHu0,q,b0HH(3,1))      
      call b0(p,mG0,mh0,q,b0HH(3,2))      
      call b0(p,mG0,mG0,q,b0HH(3,3))      
      call b0(p,mG0,mA0,q,b0HH(3,4))      

      call b0(p,mA0,mHu0,q,b0HH(4,1))      
      call b0(p,mA0,mh0,q,b0HH(4,2))      
      call b0(p,mA0,mG0,q,b0HH(4,3))      
      call b0(p,mA0,mA0,q,b0HH(4,4)) 
      
      call a0(mHu0,q,a0H(1))
      call a0(mh0,q,a0H(2))
      call a0(mG0,q,a0H(3))
      call a0(mA0,q,a0H(4))

!------------------------------------------      
      
      loop8: do i = 1, 4
      loop9: do j = 1, 4
      
      higgs = higgs + 0.5d0 * hhs1(i, j)* hhs2(i, j) * b0HH(i,j)

      enddo loop9

      higgs = higgs + 0.5d0 * hhs1s2(i) * a0H(i)

      enddo loop8


!-------------------------------------------------------------

      hphps1(1, 1) =  cos2beta * cosbeta
      hphps1(2, 2) = -cos2beta * cosbeta + 2.d0 * cossqthW * cosbeta
      hphps1(1, 2) = -sin2beta * cosbeta + cossqthW * sinbeta 
      hphps1(2, 1) =  hphps1(1, 2)
      
      hphps2(1, 1) =  - cos2beta * sinbeta
      hphps2(2, 2) =    cos2beta * sinbeta + 2.d0 * cossqthW * sinbeta
      hphps2(1, 2) =    sin2beta * sinbeta - cossqthW * cosbeta 
      hphps2(2, 1) =    hphps2(1, 2)
      
      loop10: do i = 1, 2
      loop11: do j = 1, 2

      hphps1(i,j) = hphps1(i,j) * (g * MZ * 0.5d0 / costhW)
      hphps2(i,j) = hphps2(i,j) * (g * MZ * 0.5d0 / costhW)
      
      enddo loop11
      enddo loop10

!------------------------------      
 
      hphps1s2(1) = - cossqthW * sin2beta
      hphps1s2(2) =   cossqthW * sin2beta 
      
      loop12: do i = 1, 2
      hphps1s2(i) = hphps1s2(i) * (g*g* 0.25d0 / cossqthW)
      enddo loop12

      call b0(p,mGp,mGp,q,b0HpHp(1,1))      
      call b0(p,mGp,mHpm,q,b0HpHp(1,2))      
      call b0(p,mHpm,mGp,q,b0HpHp(2,1))      
      call b0(p,mHpm,mHpm,q,b0HpHp(2,2))      
      
      call a0(mGp,q,a0Hp(1))
      call a0(mHpm,q,a0Hp(2))
            
      loop13: do i = 1, 2
      loop14: do j = 1, 2
      
      higgs = higgs + hphps1(i, j) * hphps2(i, j)* b0HpHp(i,j)

      enddo loop14

      higgs = higgs + hphps1s2(i) * a0Hp(i)

      enddo loop13

!----------------------------------------------------------------------------
      
      aPsi1(1, 3) = - gp * 0.5d0
      aPsi1(3, 1) =   aPsi1(1, 3) 
      aPsi1(2, 3) =   g  * 0.5d0
      aPsi1(3, 2) =   aPsi1(2, 3) 

      aPsi2(1, 4) =   gp * 0.5d0 
      aPsi2(4, 1) =   aPsi2(1, 4)
      aPsi2(2, 4) = - g  * 0.5d0 
      aPsi2(4, 2) =   aPsi2(2, 4) 
    
      call dag4d(ON,ONdag)
      call mat3prod4d(ON,aPsi1,ONdag,aChi1)
      call mat3prod4d(ON,aPsi1,ONdag,bChi1)  
      call mat3prod4d(ON,aPsi2,ONdag,aChi2)
      call mat3prod4d(ON,aPsi2,ONdag,bChi2)  

      loop17: do i = 1, 4
      loop181: do j = 1, 4
      
      call funcg(p,dabs(mneut(i)),dabs(mneut(j)),q,gmneut(i,j))
      call b0(p,dabs(mneut(i)),dabs(mneut(j)),q,b0mneut(i,j))
      
      enddo loop181
      enddo loop17


      neutralinos = 0.d0

      loop15: do i = 1, 4
      loop16: do j = 1, 4
      
      fChiChis1s2(i, j) = real(aChi1(i, j)*(aChi2(i,j)) + 
     $     bChi1(i, j)*(bChi2(i,j)))

      gChiChis1s2(i, j) = real(bChi1(i, j)  * aChi2(i, j) + 
     $     (aChi1(i, j)) * bChi2(i, j))

      neutralinos = neutralinos + 0.5d0 * 
     $     ((fChiChis1s2(i, j) * gmneut(i,j)) - 2.d0 *
     $     (gChiChis1s2(i, j) * mneut(i) * mneut(j) * b0mneut(i,j)))

      enddo loop16
      enddo loop15

!-------------------------------------------------------------------------------

      aPsic1(1, 2) = g / dsqrt(2.d0)
      aPsic2(2, 1) = g / dsqrt(2.d0)
      
      call dag2d(aPsic1,aPsic1dag)
      call dag2d(aPsic2,aPsic2dag)
      call dag2d(OCL,OCLdag)
      call dag2d(OCR,OCRdag)
      call mat3prod2d(OCR,aPsic1,OCLdag,aChic1)    !<-------------CHECKED!!
      call mat3prod2d(OCL,aPsic1dag,OCRdag,bChic1)
      call mat3prod2d(OCR,aPsic2,OCLdag,aChic2)    !<-------------CHECKED!!
      call mat3prod2d(OCL,aPsic2dag,OCRdag,bChic2)

      loop19: do i = 1, 2
      loop18: do j = 1, 2
      
      call funcg(p,mchargino(i),mchargino(j),q,gmch(i,j))
      call b0(p,mchargino(i),mchargino(j),q,b0mch(i,j))
      
      enddo loop18
      enddo loop19

      charginos = 0.d0
      
      
      loop23: do i = 1, 2
      loop24: do j = 1, 2

      fChiChis1s2(i, j) = real(aChic1(i, j) * (aChic2(i,j)) +
     $     bChic1(i, j) * (bChic2(i,j)))

      gChiChis1s2(i, j) = real(bChic1(i, j) * aChic2(i, j) + 
     $     (aChic1(i, j)) * bChic2(i, j))

      charginos = charginos + 
     $     (fChiChis1s2(i, j) * gmch(i,j) - 2.d0 *
     $     gChiChis1s2(i, j) * mchargino(i) * mchargino(j) * b0mch(i,j))

      enddo loop24
      enddo loop23

!-------------------------------------------------------------------------------

      pis1s2ans = (sups  + sdowns  + sleps  + stops  + sbots  + staus  + 
     $     sneutrinos + fermions + higgs + neutralinos + charginos)
     $     / (16.d0*pi*pi)

!------------------------------------------------------------------------------
      return 
      end subroutine pis1s2

!------------------------------------------------------------------------------
!==============================================================================
