****f* SuSeFLAV/onelooptadpole.f 
*  NAME
*    subroutine tadpole1,tadpole2
*  SYNOPSIS
*    Computes One loop tadpoles. 
*
*  FUNCTION
*     The routine calcultes one loop tadpoles
*     From hep-ph/9606211's appendix. It should be done at MSusy to minimize the
*     1-loop contributions. 
*
*  INPUTS
*     p                              - External momentum  
*     q                              - Energy scale
*     g,gp,g3                        - Gauge couplings(g = em, gp = weak, 
*                                      g3 = strong)                                     
*     mt,mb,mtau                     - pole masses of top, bottom and tau
*     mSQRG,mSDRG,mSURG,mSLRG,mSERG  - (3 X 3) mass matrix definition
*     yuRG,ydRG,yeRG                 - (3 X 3) Yukawas
*     AURG,ADRG,AERG                 - (3 X 3) Trilinear couplings
*     pizzT,piwwT                    - self energy of W and Z bosons at M_z
*     mh0sq,mhu0sq,mhpmsq,mA0sq      - physical higgs mass squared 
*     modmu                          - modulus of the \mu paramter 
*     vev1,vev2                      - vacuum expectation values of the two 
*                                      higgs doublet fields
*     M3t                            - Gaugino mass at msusy
*     tanbeta                        - the ratio of the vevs of the 
*                                      two Higgs doublet fields.
*     SUegg   =  6 eigenvalues of UP-Squark mass matrix.
*     SDegg   =  6 eigenvalues of Down-Squark mass matrix.
*     Slegg   =  6 eigenvalues of slepton mass matrix.
*     SNegg   =  3 eigenvalues of sneutrino mass matrix.
*     ON      =  (4 X 4) orthogonal matrix such that 
*                ON.MNeut.Transpose[ON] = Diag[MNeut] 
*     Neg     =  4 singular values (descending order) of the 
*                Neutralino mass matrix. 
*
*     OCR, OCL =  (2 X 2) orthogonal matrices such that 
*                 MChar = Transpose[OCR].Diag[MChar].OCL
*     Ceg      =   2 singular values of the chargino Mass Matrix
*
*  RESULT
*     delta1,delta2   =  tadpoles at a given scale(generally msusy).
*  
*  EXAMPLE
*
*      subroutine tadpole1(q,g,mt,mb,mtau,tanbeta,mSQRG,mSDRG,mSURG,
*     $     mSLRG,mSERG,AURG,ADRG,AERG,yuRG,ydRG,yeRG,SUegg,SDegg,SLegg,
*     $     SNegg,Neg,Ceg,ON,OCL,OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,
*     $     mA0sq,vev1,vev2,pizzT,piwwT,delta1)
*
*      subroutine tadpole2(q,g,mt,mb,mtau,tanbeta,mSQRG,mSDRG,mSURG,
*     $     mSLRG,mSERG,AURG,ADRG,AERG,yuRG,ydRG,yeRG,SUegg,SDegg,SLegg,
*     $     SNegg,Neg,Ceg,ON,OCL,OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,
*     $     mA0sq,vev1,vev2,pizzT,piwwT,delta2)
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
!================================================================================================

!                                    TADPOLE 1 STARTS HERE

!=================================================================================================
C//////////////////////////////////////////////////////////////////////////////////
C  1. checked: 22nd may, 2010 
C
C//////////////////////////////////////////////////////////////////////////////////      
      
      subroutine tadpole1(q,g,mb,mtau,tanbeta,
     $     ADRG,AERG,yuRG,ydRG,yeRG,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,ON,OCL,OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,
     $     mA0sq,vev1,delta1)

      implicit none
      integer i,rhn,lopt
      double precision q,g,delta1,sfermions,stops,sbots,staus,sneuts
      double precision sleps,sups,sdowns,fermions,higgs,charginos
      double precision neutralinos,gaugeBosons,gb2mwcosb,gsqmzb2mwcw

      double precision mB, mTau
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      double precision rot2dt(2,2),trot2dt(2,2),rot2db(2,2),trot2db(2,2)
      double precision rot2dtau(2,2),trot2dtau(2,2)
      double precision ltops1lr(2,2),ltops112(2,2),lbots1lr(2,2) 
      double precision lbots112(2,2),ltaus1lr(2,2),ltaus112(2,2),lsnul
      double precision a0mb,a0mtau,a0mt1,a0mt2,a0mb1,a0mb2,a0mtau1
      double precision a0mtau2,a0msnu1,a0msnu2,a0msnu3,a0muR,a0muL,a0mcR
      double precision a0mcL,a0mdR,a0mdL,a0msR,a0msL,a0meR,a0meL,a0mmuR
      double precision a0mmuL,a0mA,a0mHpm,a0mh,a0mHu
      double precision a0mneut(4),a0mch(2),a0MZ,a0MW 

      double precision ON(4,4),OCL(2,2),OCR(2,2)
      double precision ht,hb,htau,hbsq,htausq,modmu,mh0,mHpm,mhu0,mA0
      
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION ADRG(3,3)
      DOUBLE PRECISION AERG(3,3),SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq

      DOUBLE PRECISION muL,muR,mcL,mcR,mt1,mt2
      DOUBLE PRECISION mdL,mdR,msL,msR,mb1,mb2
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtau1,mtau2      
      DOUBLE PRECISION mneut(4),msnu(3),mchargino(2)

      double precision cossqthw,sinsqthw,cosbeta,costhw,sinthw,tanthW
      double precision sin2a,sinasq,cosasq,alpha,thw,tan2beta
      double precision thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae
      double precision yuRG(3,3),ydRG(3,3),yeRG(3,3)

      double precision beta,sbeta,cbeta
      double precision vev1,cos2beta

      DOUBLE PRECISION thetatz,thetabz,thetatauz
      DOUBLE PRECISION thetacz,thetasz,thetamuz
      DOUBLE PRECISION thetauz,thetadz,thetaez
      DOUBLE PRECISION thetat_DUM,thetab_DUM,thetatau_DUM
      DOUBLE PRECISION sinsqthw_susy,sgnmu,MZrun, MWrun

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole


      common/sfmixing_mz/thetatz,thetabz,thetatauz,thetacz,thetasz,
     $     thetamuz,thetauz,thetadz,thetaez

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae
      common/loops/ lopt,rhn
      common/gbrunning/ MZrun, MWrun

!-----------------------------------------------------------------
      external dag2d,mat3prod2d,a0,theta

      include 'stdinputs.h'
!--------------------------------------------------------------------------


c$$$      muR = dsqrt(SUegg(6)) 
c$$$      muL = dsqrt(SUegg(5))
c$$$      mcR = dsqrt(SUegg(4))
c$$$      mcL = dsqrt(SUegg(3))
c$$$      mt2 = dsqrt(SUegg(2))
c$$$      mt1 = dsqrt(SUegg(1))
c$$$
c$$$
c$$$      
c$$$      mdR = dsqrt(SDegg(6)) 
c$$$      mdL = dsqrt(SDegg(5))
c$$$      msR = dsqrt(SDegg(4))
c$$$      msL = dsqrt(SDegg(3))
c$$$      mb2 = dsqrt(SDegg(2))
c$$$      mb1 = dsqrt(SDegg(1))
c$$$
c$$$      meR   = dsqrt(SLegg(6))
c$$$      meL   = dsqrt(SLegg(5))
c$$$      mmuR  = dsqrt(SLegg(4))
c$$$      mmuL  = dsqrt(SLegg(3))
c$$$      mtau2 = dsqrt(SLegg(2))
c$$$      mtau1 = dsqrt(SLegg(1))



      muR = dsqrt(SUegg(5)) 
      muL = dsqrt(SUegg(6))
      mcR = dsqrt(SUegg(3))
      mcL = dsqrt(SUegg(4))
c$$$      mt2 = dsqrt(SUegg(2))
c$$$      mt1 = dsqrt(SUegg(1))

      mt1 = dsqrt(SUegg(2))
      mt2 = dsqrt(SUegg(1))


      mdR = dsqrt(SDegg(5)) 
      mdL = dsqrt(SDegg(6))
      msR = dsqrt(SDegg(3))
      msL = dsqrt(SDegg(4))
c$$$      mb2 = dsqrt(SDegg(2))
c$$$      mb1 = dsqrt(SDegg(1))

      mb1 = dsqrt(SDegg(2))
      mb2 = dsqrt(SDegg(1))

      meR   = dsqrt(SLegg(5))
      meL   = dsqrt(SLegg(6))
      mmuR  = dsqrt(SLegg(3))
      mmuL  = dsqrt(SLegg(4))
c$$$      mtau2 = dsqrt(SLegg(2))
c$$$      mtau1 = dsqrt(SLegg(1))

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

!----------------

      if(q.lt.100.d0)then

         thetat_dum = thetatz      
         thetab_dum = thetabz                
         thetatau_dum = thetatauz

      else

         thetat_dum = thetat      
         thetab_dum = thetab                
         thetatau_dum = thetatau

      endif

!--------------------------------------------------- 

!      MZ = 90.2120d0
!      MW = 79.5527d0

      MZ = MZrun
      MW = MWrun

!      print*,"MZrun, MWrun = ", MZrun, MWrun

      pi = 4.d0 * datan(1.d0)

      cosbeta = dcos(datan(tanbeta))
      beta = datan(tanbeta)
      sbeta = dsin(beta)
      cbeta = dcos(beta)
      cos2beta = cbeta*cbeta - sbeta*sbeta
      tan2beta = dtan(2.d0*datan(tanbeta))

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

!----------------------------------------------------------------
      
      rot2dt(1,1) =   dcos(thetat_dum)
      rot2dt(1,2) =   dsin(thetat_dum)
      rot2dt(2,1) = - rot2dt(1,2)
      rot2dt(2,2) =   rot2dt(1,1)

      rot2db(1,1) =   dcos(thetab_dum)
      rot2db(1,2) =   dsin(thetab_dum)
      rot2db(2,1) = - rot2db(1,2)
      rot2db(2,2) =   rot2db(1,1)

      rot2dtau(1,1) =   dcos(thetatau_dum)
      rot2dtau(1,2) =   dsin(thetatau_dum)
      rot2dtau(2,1) = - rot2dtau(1,2)
      rot2dtau(2,2) =   rot2dtau(1,1)


c$$$      print*,"thetat, thetab, thetatau = " , thetat, thetab, thetatau
c$$$      print*,"q, costhW in tadpole 1 = ", q, costhW
!----------------------------------------------------------------
      
      call dag2d(rot2dt,trot2dt)
      call dag2d(rot2db,trot2db)
      call dag2d(rot2dtau,trot2dtau)

!---------------------------------------------------------------      
C     sneutrino coupling

      lsnul = g*MZ*gnuL*cosbeta/costhW 

      ltops1lr(1, 1) = (g*MZ/costhW) * guL * cosbeta
      ltops1lr(1, 2) = - (ht * sgnmu*modmu) / dsqrt(2.d0)
      ltops1lr(2, 1) = ltops1lr(1, 2)
      ltops1lr(2, 2) = (g*MZ/costhW) * guR * cosbeta

      call mat3prod2d(rot2dt,ltops1lr,trot2dt,ltops112)

!------------------------------------------------------------------  
    
      hbsq = hb*hb

      lbots1lr(1, 1) = (g * MZ / costhW ) * gdL * cosbeta + hbsq * vev1
      lbots1lr(1, 2) =  hb * ADRG(3,3) / dsqrt(2.d0)     
      lbots1lr(2, 1) = lbots1lr(1, 2)
      lbots1lr(2, 2) = (g * MZ / costhW ) * gdR * cosbeta + hbsq * vev1
      
c$$$      print*,"hb * ADRG(3,3) = ", hb * ADRG(3,3)
c$$$      print*,"ht = ", ht, "hb = ", hb
c$$$      print*,"MW, MZ = ", MW,MZ


      call mat3prod2d(rot2db,lbots1lr,trot2db,lbots112)

!--------------------------------------------------------------------
!     stau couplings

      htausq = htau*htau
      
      ltaus1lr(1, 1) = (g * MZ / costhW) * geL * cosbeta + htausq * vev1
      ltaus1lr(1, 2) =  htau * AERG(3,3) / dsqrt(2.d0)       
      ltaus1lr(2, 1) =  ltaus1lr(1, 2)
      ltaus1lr(2, 2) = (g * MZ / costhW) * geR * cosbeta + htausq * vev1
      
      call mat3prod2d(rot2dtau,ltaus1lr,trot2dtau,ltaus112)

!----------------------------------------------------------------------------------- 
C     bottom quark and tau, ignoring others
!-----------------------------------------------------------------------------------  

      call a0(mb,q,a0mb)
      call a0(mtau,q,a0mtau)

      fermions = - (6.d0 * hbsq * a0mb) - 
     $     (2.d0 * htausq * a0mtau)

c$$$      print*,"hb = ", hb, "mb = ", mb
c$$$      print*,"htau = ", htau, "mtau = ", mtau
c$$$      print*,"vev1 = ", vev1



!---------------------------------------------------------------------------
      
      gb2mwcosb = g / (2.d0*MW*cosbeta)

      stops = 0.d0
      sbots = 0.d0
      staus = 0.d0

!---------------------------------------------------------------------------

!      print*, "mt1, mt2 = ", mt1, mt2

c$$$      print*,"ltops112(1,1), ltops112(2,2) = ", ltops112(1,1), 
c$$$     $     ltops112(2,2)

      call a0(mt1,q,a0mt1)
      call a0(mt2,q,a0mt2)
      call a0(mb1,q,a0mb1)
      call a0(mb2,q,a0mb2)

!     third generation squarks

      stops = stops + 3.d0 * gb2mwcosb * ltops112(1, 1) * a0mt1
      stops = stops + 3.d0 * gb2mwcosb * ltops112(2, 2) * a0mt2

      sbots = sbots + 3.d0 * gb2mwcosb * lbots112(1, 1) * a0mb1
      sbots = sbots + 3.d0 * gb2mwcosb * lbots112(2, 2) * a0mb2
!---------------------------------------------------------------------------
      call a0(mtau1,q,a0mtau1)
      call a0(mtau2,q,a0mtau2)

!     third generation sleptons

      staus = staus + gb2mwcosb * ltaus112(1, 1) * a0mtau1
      staus = staus + gb2mwcosb * ltaus112(2, 2) * a0mtau2
!---------------------------------------------------------------------------  
      call a0(msnu(1),q,a0msnu1)
      call a0(msnu(2),q,a0msnu2)
      call a0(msnu(3),q,a0msnu3)
    
      sneuts = 0.d0
      sneuts = sneuts + gb2mwcosb*lsnul*(a0msnu1+a0msnu2+a0msnu3)
!---------------------------------------------------------------------------

      gsqmzb2mwcw = g*g* MZ * 0.5d0 / (MW * costhW)
      sups = 0.d0
      sdowns = 0.d0

      call a0(muR,q,a0muR)
      call a0(muL,q,a0muL)
      call a0(mcR,q,a0mcR)
      call a0(mcL,q,a0mcL)
      call a0(mdR,q,a0mdR)
      call a0(mdL,q,a0mdL)
      call a0(msR,q,a0msR)
      call a0(msL,q,a0msL)

      
!     first two families of squarks
!-----------------------------------
      
      sups = sups + 3.0 * gsqmzb2mwcw * ( guL * a0muL + guR * a0muR ) 
      sups = sups + 3.0 * gsqmzb2mwcw * ( guL * a0mcL + guR * a0mcR ) 
      sdowns = sdowns + 3.0 * gsqmzb2mwcw * ( gdL * a0mdL + gdR * a0mdR)
      sdowns = sdowns + 3.0 * gsqmzb2mwcw * ( gdL * a0msL + gdR * a0msR)
      
      
      sleps = 0.d0

      call a0(meR,q,a0meR)
      call a0(meL,q,a0meL)
      call a0(mmuR,q,a0mmuR)
      call a0(mmuL,q,a0mmuL)

!     sleptons

      sleps = sleps + gsqmzb2mwcw * (geL * a0meL + geR * a0meR)
      sleps = sleps + gsqmzb2mwcw * (geL * a0mmuL + geR * a0mmuR)
 
!---------------------
C      HIGGS
!---------------------

      sinasq = dsin(alpha)*dsin(alpha)
      cosasq = dcos(alpha)*dcos(alpha)
      sin2a  = dsin(2.d0*alpha)


      higgs = 0.d0
      call a0(mA0,q,a0mA)
      call a0(mHpm,q,a0mHpm)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)

c$$$      print*, "mA0 = ", mA0, "mHpm = ", mHpm
c$$$      print*,"mh0 = ", mh0, "mHu0 = ", mHu0
c$$$      print*,"q = ", q

      higgs = higgs - ((g*g*cos2beta/(8.d0*cossqthW)) *
     $     (a0mA + 2.d0*a0mHpm)) +
     $     g*g*a0mHpm*0.5d0 + (g*g/(8.d0 * cossqthW))*a0mh*
     $     (3.d0*sinasq - cosasq + sin2a * tanbeta) +
     $     (g*g/(8.d0 * cossqthW))*a0mHu*
     $     (3.d0 * cosasq - sinasq - sin2a * tanbeta)

!-------------------------      
!     Neutralinos
!-------------------------
 
      call a0(dabs(mneut(1)),q,a0mneut(1))
      call a0(dabs(mneut(2)),q,a0mneut(2))
      call a0(dabs(mneut(3)),q,a0mneut(3))
      call a0(dabs(mneut(4)),q,a0mneut(4))

      neutralinos = 0.d0
      tanthW = dtan(thW)

      loop1: do i = 1, 4

      neutralinos = neutralinos - ((g*g*mneut(i) / (MW * cosbeta)) *
     $     (ON(i, 3) * (ON(i, 2) - ON(i, 1) * tanthW)) * a0mneut(i))
      
      enddo loop1


!---------------------
!     Charginos
!---------------------
      charginos = 0.d0 

      call a0(dabs(mchargino(1)),q,a0mch(1))
      call a0(dabs(mchargino(2)),q,a0mch(2))

      loop2: do i = 1, 2

      charginos = charginos - dsqrt(2.d0)*(g*g / (MW * cosbeta))
     $     * mchargino(i) * (OCR(i, 1) * OCL(i, 2)) * a0mch(i)

      enddo loop2


!---------------------      
C      gauge bosons
!--------------------
      
      gaugeBosons = 0.d0

      call a0(MZ,q,a0MZ)
      call a0(MW,q,a0MW)
      
      gaugeBosons = gaugeBosons + ((3.d0*g*g/4.d0)*(2.d0*a0MW + 
     $     (a0MZ/cossqthW))) + ((g*g*cos2beta/(8.d0*cossqthW))*
     $     (2.d0*a0MW + a0MZ))

!---------------------------------------------------------------------------

c$$$      if(rhn.eq.1)then
c$$$
c$$$         staus = 0.d0
c$$$         sleps = 0.d0
c$$$         sneuts = 0.d0
c$$$
c$$$      endif

      sfermions = stops + sbots + staus + sneuts + sleps + sups + sdowns

      delta1 = 0.d0
      
      delta1 = (fermions + sfermions + higgs + charginos + neutralinos + 
     $     gaugeBosons)/(16.d0*pi*pi)

c$$$      print*,'In tadpole1 '
c$$$      print*,'fermions ', fermions
c$$$      print*,'sfermions ', sfermions
c$$$      print*,'higgs ', higgs
c$$$      print*,'charginos ', charginos
c$$$      print*,'neutralinos ', neutralinos
c$$$      print*,'gaugebosons ', gaugebosons
c$$$
c$$$      print*,"stops = ", stops
c$$$      print*,"sbots = ", sbots
c$$$      print*,"staus = ", staus
c$$$      print*,"sneuts = ", sneuts
c$$$      print*,"sleps = ", sleps
c$$$      print*,"sups = ", sups
c$$$      print*,"sdowns = ", sdowns

!----------------------------------------------------------------------------

      return 
      end subroutine tadpole1

!===============================================================================

!                                  TADPOLE 2 STARTS HERE

!=================================================================================================

      subroutine tadpole2(q,g,mt,tanbeta,
     $     AURG,yuRG,ydRG,yeRG,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,ON,OCL,OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,
     $     mA0sq,vev2,delta2)

      
      implicit none
      integer i,lopt,rhn
      double precision q,g,delta2,sfermions,stops,sbots,staus,sneuts
      double precision sleps,sups,sdowns,fermions,higgs,charginos
      double precision neutralinos,gaugeBosons,gb2mwsinb,gsqmzb2mwcw


      double precision mT
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      double precision rot2dt(2,2),trot2dt(2,2),rot2db(2,2),trot2db(2,2)
      double precision rot2dtau(2,2),trot2dtau(2,2)
      double precision ltops1lr(2,2),ltops212(2,2),lbots1lr(2,2)
      double precision lbots212(2,2),ltaus1lr(2,2),ltaus212(2,2),lsnul
      double precision a0mtop,a0mt1,a0mt2,a0mb1,a0mb2,a0mtau1
      double precision a0mtau2,a0msnu1,a0msnu2,a0msnu3,a0muR,a0muL,a0mcR
      double precision a0mcL,a0mdR,a0mdL,a0msR,a0msL,a0meR,a0meL,a0mmuR
      double precision a0mmuL,a0mA,a0mHpm,a0mh,a0mHu
      double precision a0mneut(4),a0mch(2),a0MZ,a0MW 

      double precision ON(4,4),OCL(2,2),OCR(2,2)
      double precision ht,hb,htau,htsq,modmu,mh0,mHpm,mhu0,mA0 !<------
      
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION AURG(3,3)
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq

      DOUBLE PRECISION muL,muR,mcL,mcR,mt1,mt2
      DOUBLE PRECISION mdL,mdR,msL,msR,mb1,mb2
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtau1,mtau2      
      DOUBLE PRECISION mneut(4),msnu(3),mchargino(2)

      double precision cossqthw,sinsqthw,cosbeta,costhw,sinthw,tanthw
      double precision sin2a,sinasq,cosasq,alpha,thw,tan2beta,sinbeta
      double precision thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae
      double precision yuRG(3,3),ydRG(3,3),yeRG(3,3)

      double precision beta,sbeta,cbeta
      double precision vev2,cos2beta

      DOUBLE PRECISION thetatz,thetabz,thetatauz
      DOUBLE PRECISION thetacz,thetasz,thetamuz
      DOUBLE PRECISION thetauz,thetadz,thetaez
      DOUBLE PRECISION thetat_DUM,thetab_DUM,thetatau_DUM

      DOUBLE PRECISION sinsqthw_susy,sgnmu,MZrun, MWrun

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole


      common/sfmixing_mz/thetatz,thetabz,thetatauz,thetacz,thetasz,
     $     thetamuz,thetauz,thetadz,thetaez

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

      common/loops/ lopt,rhn
      common/gbrunning/ MZrun, MWrun

!------------------------------------------------------------
      external dag2d,mat3prod2d,a0,theta

      include 'stdinputs.h'

!-----------------------------------      

c$$$      muR = dsqrt(SUegg(6)) 
c$$$      muL = dsqrt(SUegg(5))
c$$$      mcR = dsqrt(SUegg(4))
c$$$      mcL = dsqrt(SUegg(3))
c$$$      mt2 = dsqrt(SUegg(2))
c$$$      mt1 = dsqrt(SUegg(1))
c$$$
c$$$
c$$$      
c$$$      mdR = dsqrt(SDegg(6)) 
c$$$      mdL = dsqrt(SDegg(5))
c$$$      msR = dsqrt(SDegg(4))
c$$$      msL = dsqrt(SDegg(3))
c$$$      mb2 = dsqrt(SDegg(2))
c$$$      mb1 = dsqrt(SDegg(1))
c$$$
c$$$      meR   = dsqrt(SLegg(6))
c$$$      meL   = dsqrt(SLegg(5))
c$$$      mmuR  = dsqrt(SLegg(4))
c$$$      mmuL  = dsqrt(SLegg(3))
c$$$      mtau2 = dsqrt(SLegg(2))
c$$$      mtau1 = dsqrt(SLegg(1))


      muR = dsqrt(SUegg(5)) 
      muL = dsqrt(SUegg(6))
      mcR = dsqrt(SUegg(3))
      mcL = dsqrt(SUegg(4))
c$$$      mt2 = dsqrt(SUegg(2))
c$$$      mt1 = dsqrt(SUegg(1))

      mt1 = dsqrt(SUegg(2))
      mt2 = dsqrt(SUegg(1))

      mdR = dsqrt(SDegg(5)) 
      mdL = dsqrt(SDegg(6))
      msR = dsqrt(SDegg(3))
      msL = dsqrt(SDegg(4))
c$$$      mb2 = dsqrt(SDegg(2))
c$$$      mb1 = dsqrt(SDegg(1))

      mb1 = dsqrt(SDegg(2))
      mb2 = dsqrt(SDegg(1))

      meR   = dsqrt(SLegg(5))
      meL   = dsqrt(SLegg(6))
      mmuR  = dsqrt(SLegg(3))
      mmuL  = dsqrt(SLegg(4))
c$$$      mtau2 = dsqrt(SLegg(2))
c$$$      mtau1 = dsqrt(SLegg(1))

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

!----------------


      if(q.lt.100.d0)then
         thetat_dum = thetatz      
         thetab_dum = thetabz                
         thetatau_dum = thetatauz

      else

         thetat_dum = thetat      
         thetab_dum = thetab                
         thetatau_dum = thetatau

      endif
      
!-----------------------------------------
!     definitions.
!-----------------------------------------

      pi = 4.d0 * datan(1.d0)
      
      MZ = MZrun                !90.074d0             !MZpole
      MW = MWrun                !78.512d0             !MWpole

c$$$      print*,"MZrun, MWrun = ", MZrun, MWrun

      cosbeta = dcos(datan(tanbeta))
      sinbeta = dsin(datan(tanbeta))
      beta = datan(tanbeta)
      sbeta = dsin(beta)
      cbeta = dcos(beta)
      cos2beta = cbeta*cbeta - sbeta*sbeta
      tan2beta = dtan(2.d0*datan(tanbeta))

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

!      print*,"ht in tagpole 2 = ", ht

!------------------------------------------------
      rot2dt(1,1) =   dcos(thetat_dum)
      rot2dt(1,2) =   dsin(thetat_dum)
      rot2dt(2,1) = - rot2dt(1,2)
      rot2dt(2,2) =   rot2dt(1,1)

      rot2db(1,1) =   dcos(thetab_dum)
      rot2db(1,2) =   dsin(thetab_dum)
      rot2db(2,1) = - rot2db(1,2)
      rot2db(2,2) =   rot2db(1,1)

      rot2dtau(1,1) =   dcos(thetatau_dum)
      rot2dtau(1,2) =   dsin(thetatau_dum)
      rot2dtau(2,1) = - rot2dtau(1,2)
      rot2dtau(2,2) =   rot2dtau(1,1)

c$$$      print*,"thetat, thetab, thetatau = ", thetat, thetab, thetatau
c$$$      print*,"costhW = ", costhW
!----------------------------------------------------------------
      
      call dag2d(rot2dt,trot2dt)
      call dag2d(rot2db,trot2db)
      call dag2d(rot2dtau,trot2dtau)

!---------------------------------------------------------------      

      lsnul = - (g * MZ / costhW) * gnuL * sinbeta 
!---------------------------------------------------------------

      htsq = ht*ht

      ltops1lr(1, 1) = - (g * MZ / costhW) * guL * sinbeta + htsq * vev2 
      ltops1lr(1, 2) = (ht * AURG(3,3)) / dsqrt(2.d0)
      ltops1lr(2, 1) = ltops1lr(1, 2)
      ltops1lr(2, 2) = - (g * MZ / costhW )* guR * sinbeta + htsq * vev2 

      call mat3prod2d(rot2dt,ltops1lr,trot2dt,ltops212)


!      print*,"ht * AURG(3,3) = ", ht * AURG(3,3)

!------------------------------------------------------------------      

      lbots1lr(1, 1) = - (g * MZ / costhW) * gdL * sinbeta     
      lbots1lr(1, 2) =  -  hb * sgnmu*modmu / dsqrt(2.d0)            
      lbots1lr(2, 1) =    lbots1lr(1, 2)
      lbots1lr(2, 2) = - (g * MZ / costhW ) * gdR * sinbeta    
      
      call mat3prod2d(rot2db,lbots1lr,trot2db,lbots212)
!--------------------------------------------------------------------
      
      ltaus1lr(1, 1) = - (g * MZ / costhW) * geL * sinbeta      
      ltaus1lr(1, 2) = - htau * sgnmu*modmu / dsqrt(2.d0)             
      ltaus1lr(2, 1) =  ltaus1lr(1, 2)
      ltaus1lr(2, 2) = - (g * MZ / costhW) * geR * sinbeta      
      
      call mat3prod2d(rot2dtau,ltaus1lr,trot2dtau,ltaus212)

!--------------------------------------------------------------------------
C     top quark, ignoring others 
!--------------------------------------------------------------------------

      call a0(mT,q,a0mtop)
      
      fermions = 0.d0
     

      fermions =  - 6.0 * htsq * a0mtop       

c$$$      print*,"ht = ", ht, "mt = ", mt, "q = ", q 
c$$$      print*,"vev2 = ", vev2
!--------------------------------------------------------------------------
      

      gb2mwsinb = g / (2.d0*MW*sinbeta)
      stops = 0.d0
      sbots = 0.d0
      staus = 0.d0
!---------------------------------------------------------------------------
      call a0(mt1,q,a0mt1)
      call a0(mt2,q,a0mt2)
      call a0(mb1,q,a0mb1)
      call a0(mb2,q,a0mb2)

!---------------------------------------------------------------------------
      stops = stops + 3.d0 * gb2mwsinb * ltops212(1, 1) * a0mt1
      stops = stops + 3.d0 * gb2mwsinb * ltops212(2, 2) * a0mt2
 
      sbots = sbots + 3.d0 * gb2mwsinb * lbots212(1, 1) * a0mb1
      sbots = sbots + 3.d0 * gb2mwsinb * lbots212(2, 2) * a0mb2
!---------------------------------------------------------------------------
      call a0(mtau1,q,a0mtau1)
      call a0(mtau2,q,a0mtau2)

!---------------------------------------------------------------------------
      staus = staus + gb2mwsinb * ltaus212(1, 1) * a0mtau1
      staus = staus + gb2mwsinb * ltaus212(2, 2) * a0mtau2
!---------------------------------------------------------------------------  
      call a0(msnu(1),q,a0msnu1)
      call a0(msnu(2),q,a0msnu2)
      call a0(msnu(3),q,a0msnu3)
    
      sneuts = 0.d0
      sneuts = sneuts + gb2mwsinb*lsnul*(a0msnu1+a0msnu2+a0msnu3)
!---------------------------------------------------------------------------

      gsqmzb2mwcw = - g*g* MZ * 0.5d0 / (MW * costhW)
      sups = 0.d0
      sdowns = 0.d0

      call a0(muR,q,a0muR)
      call a0(muL,q,a0muL)
      call a0(mcR,q,a0mcR)
      call a0(mcL,q,a0mcL)
      call a0(mdR,q,a0mdR)
      call a0(mdL,q,a0mdL)
      call a0(msR,q,a0msR)
      call a0(msL,q,a0msL)

      
!---------------------------------------------------------------------------
      
      sups = sups + 3.0 * gsqmzb2mwcw * ( guL * a0muL + guR * a0muR ) 
      sups = sups + 3.0 * gsqmzb2mwcw * ( guL * a0mcL + guR * a0mcR ) 

      sdowns = sdowns + 3.0 * gsqmzb2mwcw * ( gdL * a0mdL + gdR * a0mdR)
      sdowns = sdowns + 3.0 * gsqmzb2mwcw * ( gdL * a0msL + gdR * a0msR)
      
     
!------------------------------------

      sleps = 0.d0

      call a0(meR,q,a0meR)
      call a0(meL,q,a0meL)
      call a0(mmuR,q,a0mmuR)
      call a0(mmuL,q,a0mmuL)

      sleps = sleps + gsqmzb2mwcw * (geL * a0meL + geR * a0meR)
      sleps = sleps + gsqmzb2mwcw * (geL * a0mmuL + geR * a0mmuR)

 

!----------------
C     HIGGS
!----------------

      sinasq = dsin(alpha)*dsin(alpha)
      cosasq = dcos(alpha)*dcos(alpha)
      sin2a  = dsin(2.d0*alpha)


      higgs = 0.d0
      call a0(mA0,q,a0mA)
      call a0(mHpm,q,a0mHpm)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)

c$$$      print*, "mA0 = ", mA0, "mHpm = ", mHpm
c$$$      print*,"mh0 = ", mh0, "mHu0 = ", mHu0
c$$$      print*,"q = ", q, "g = ", g
c$$$      print*," cos2beta = ", cos2beta
c$$$      print*,"cossqthW = ", cossqthW
c$$$      print*,"cosasq = ", cosasq

      higgs = higgs + ((g*g*cos2beta/(8.d0*cossqthW)) *
     $     (a0mA + 2.d0*a0mHpm)) +
     $     g*g*a0mHpm*0.5d0 + ((g*g/(8.d0 * cossqthW))*a0mh*
     $     (3.d0*cosasq - sinasq + (sin2a/tanbeta))) +
     $     ((g*g/(8.d0 * cossqthW))*a0mHu *
     $     (3.d0 * sinasq - cosasq - (sin2a/tanbeta)))
      
!-------------------------------------------      

      call a0(dabs(mneut(1)),q,a0mneut(1))
      call a0(dabs(mneut(2)),q,a0mneut(2))
      call a0(dabs(mneut(3)),q,a0mneut(3))
      call a0(dabs(mneut(4)),q,a0mneut(4))

      neutralinos = 0.d0
      
      tanthW = dtan(thW)

      loop1: do i = 1, 4

      neutralinos = neutralinos + (g*g*mneut(i) / (MW * sinbeta)) *
     $     (ON(i, 4) * (ON(i, 2) - ON(i, 1) * tanthW)) * a0mneut(i)
      
      enddo loop1

!--------------------

      charginos = 0.d0 

      call a0(dabs(mchargino(1)),q,a0mch(1))
      call a0(dabs(mchargino(2)),q,a0mch(2))

      loop2: do i = 1, 2

      charginos = charginos - dsqrt(2.d0)*(g*g / (MW * sinbeta)) * 
     $     mchargino(i) * (OCR(i, 2) * OCL(i, 1)) * a0mch(i)

      enddo loop2

!---------------------      
C      Weak bosons
!---------------------      
      gaugeBosons = 0.d0

      call a0(MZ,q,a0MZ)
      call a0(MW,q,a0MW)
      
      gaugeBosons = gaugeBosons + ((3.d0*g*g/4.d0)*(2.d0*a0MW + 
     $     (a0MZ/cossqthW))) - ((g*g*cos2beta/(8.d0*cossqthW)) *
     $     (2.d0*a0MW + a0MZ))

!------------------------------------------------------------------------------------
c$$$      if(rhn.eq.1)then
c$$$         staus = 0.d0
c$$$         sleps = 0.d0
c$$$         sneuts = 0.d0
c$$$         endif

      sfermions = stops + sbots + staus + sneuts + sleps + sups + sdowns
      
      delta2 = 0.d0

      delta2 = (fermions + sfermions + higgs + charginos + neutralinos + 
     $     gaugeBosons)/(16.d0*pi*pi)


c$$$      print*,'In tadpole2 '
c$$$      print*,'fermions ', fermions
c$$$      print*,'sfermions ', sfermions
c$$$      print*,'higgs ', higgs
c$$$      print*,' charginos ', charginos
c$$$      print*,'neutralinos ', neutralinos
c$$$      print*,'gaugebosons ', gaugebosons

!----------------------------------------------------------------------------

      return 
      end subroutine tadpole2
C============================================================================
