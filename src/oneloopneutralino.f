****f* SuSeFLAV/oneloopneutralino.f 
*  NAME
*    subroutine neutralino
*  SYNOPSIS
*    One loop correction to neutralino. 
*  FUNCTION
*     Computes 1-loop correction to neutralinos at a given energy scale and 
*     external momenta
*
*  INPUTS
*     p                              - External momentum  
*     q                              - Energy scale
*     g,gp,g3                        - Gauge couplings(g = em, gp = weak, 
*                                      g3 = strong)                                     
*     mt,mb,mtau                     - pole masses of top, bottom and tau
*     yuRG,ydRG,yeRG                 - (3 X 3) Yukawas
*     mh0sq,mhu0sq,mhpmsq,mA0sq      - physical higgs mass squared 
*     tanbeta                        - the ratio of the vevs of the 
*                                      two Higgs doublet fields.
*     SUegg   =  6 eigenvalues of UP-Squark mass matrix.
*     SDegg   =  6 eigenvalues of Down-Squark mass matrix.
*     SLegg   =  6 eigenvalues of slepton mass matrix.
*     SNegg   =  3 eigenvalues of sneutrino mass matrix.
*     ON      =  (4 X 4) orthogonal matrix such that 
*                ON.MNeut.Transpose[ON] = Diag[MNeut] 
*     Neg     =  4 singular values (descending order) of the 
*                Neutralino mass matrix. 
*
*     OCR, OCL =  (2 X 2) orthogonal matrices such that 
*                 MChar = Transpose[OCR].Diag[MChar].OCL
*     Ceg     =   2 singular values of the chargino Mass Matrix
*
*  RESULT
*     Negm            =   4 eigenvalues of the 1-loop corrected neutralino
*                             mass matrix.
*     ONew            =   (4 X 4) corrected orthogonal matrix such that 
*                       MChar = Transpose[Onew].Diag[MNeut].ONew
*     neutmassmasstot = (4x4) corrected neutralino mass matrix
*
*  EXAMPLE
*
*      subroutine neutralino(p,q,g,gp,mt,mb,mtau,tanbeta,yuRG,
*     $     ydRG,yeRG,SUegg,SDegg,SLegg,SNegg,MNeut1,Neg,Ceg,
*     $     ON,OCL,OCR,mh0sq,mhu0sq,mHpmsq,mA0sq,ONew,Negm,neutmasstot)
*
*
*  NOTES
*    1. q, the energy scale at which the corrections are added = msusy.
*    2. Conventions and notations followed are that of BPMZ.
*    3. p = msusy
*
*  BUGS
*    ---
*  SEE ALSO
*
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C  1. couplings checked.
C  3. all expressions checked - 27th May, 2010 
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


      subroutine neutralino(p,q,g,gp,mt,mb,mtau,tanbeta,yuRG,
     $     ydRG,yeRG,SUegg,SDegg,SLegg,SNegg,MNeut1,Neg,Ceg,
     $     ON,OCL,OCR,mh0sq,mhu0sq,mHpmsq,mA0sq,ONew,Negm,
     $     neutmasstot)

      
      implicit none

      integer i,j,k,a,m
      integer AOK,INFO,LWORK

      parameter(lwork = 35)
      
      double precision work(lwork)
      double precision p,q,g,gp
      double precision mT, mB, mTau
      double precision yuRG(3,3),ydRG(3,3),yeRG(3,3)
      
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR

      double precision ON(4,4),OCL(2,2),OCR(2,2),ONdag(4,4),OCLdag(2,2)
      double precision OCRdag(2,2)
      double precision ht,htsq,hb,htau,hbsq,htausq
      double precision mh0,mHpm,mhu0,mA0
      double precision MNeut1(4,4)
      double precision sinbeta

      double precision mG0,mGp,ONew(4,4)
   
      double precision rthetat(2,2),rthetab(2,2),rthetatau(2,2)
      data rthetat/ 4 * 0.d0/, rthetab/ 4 * 0.d0/, rthetatau/ 4 * 0.d0/

      
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq

      DOUBLE PRECISION msup(2),mscharm(2),mstop(2)
      DOUBLE PRECISION msdown(2),msstrange(2),msbot(2)
      DOUBLE PRECISION msel(2),msmu(2),mstau(2),msnu(3)
      DOUBLE PRECISION mneut(4),mchargino(2),mneutm(4)


      double precision cossqthw,sinsqthw,cosbeta,costhw,sinthw
      double precision alpha,thw,tan2beta
      double precision thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae

      double precision beta
      double precision cos2beta,sin2beta


      double precision Negm(4),Neuevim(4)
      double precision ONLm(4,4),ONRm(4,4),neutmass1(4,4)
      data neutmass1/ 16 * 0.d0/
!----------------------------------------------------------------------------
      double precision apsiuq(4, 2), apsidq(4, 2), apsiel(4, 2)
      double precision apsinul(4, 2)
      data apsiuq/ 8 * 0.d0/, apsidq/ 8 * 0.d0/
      data apsiel/ 8 * 0.d0/, apsinul/ 8 * 0.d0/
      double precision bpsiuq(4, 2), bpsidq(4, 2), bpsiel(4, 2)
      double precision bpsinul(4, 2)
      data bpsiuq/ 8 * 0.d0/, bpsidq/ 8 * 0.d0/
      data bpsiel/ 8 * 0.d0/, bpsinul/ 8 * 0.d0/
      double precision apsitq(4, 2), apsibq(4, 2), apsitaul(4, 2)
      double precision apsinutl(4, 2)
      data apsitq/ 8 * 0.d0/, apsibq/ 8 * 0.d0/, apsitaul/ 8 * 0.d0/
      data apsinutl/ 8 * 0.d0/
      double precision bpsitq(4, 2), bpsibq(4, 2), bpsitaul(4, 2)
      double precision bPsiNutnut(4, 2)
      data bpsitq/ 8 * 0.d0/, bpsibq/ 8 * 0.d0/, bpsitaul/ 8 * 0.d0/
      data bPsiNutnut/ 8 * 0.d0/
      double precision aPsi0PsicW(4,2),bPsi0PsicW(4,2),aPsi0ChicW(4,2)
      double precision bPsi0ChicW(4, 2)
      data aPsi0PsicW/ 8 * 0.d0/,bPsi0PsicW/ 8 * 0.d0/
      data aPsi0ChicW/ 8 * 0.d0/, bPsi0ChicW/ 8 * 0.d0/
      double precision aPsiPsiZ(4, 4), bPsiPsiZ(4, 4), aPsiChiZ(4, 4)
      double precision bPsiChiZ(4, 4)
      data aPsiPsiZ/ 16 * 0.d0/, bPsiPsiZ/ 16 * 0.d0/
      data aPsiChiZ/ 16 * 0.d0/, bPsiChiZ/ 16 * 0.d0/
      double precision b1mneut(4),b0mneut(4)

      double precision t1(2),t2(2),sigmaL(4,4),sigmaR(4,4),sigmaS(4,4)
      
      double precision b0msup(2),sinsqthw_susy
      double precision b1msup(2),b1msdown(2),b1msel(2),b1mscharm(2)
      double precision b1msstrange(2),b1msmu(2),b1mstop(2),b1msbot(2)
      double precision b1mstau(2),b0mstop(2),b0msbot(2),b0mstau(2)
      double precision b1msnue,b1msnumu,b1msnutau,b1mchar(2),b0mchar(2)
      double precision deltaM(4,4),neutmass(4,4),delta1(4,4),delta2(4,4)
      double precision neutmasstot(4,4),deltatot(4,4),tdeltatot(4,4)
!      data neutmasstot/ 16 * 0.d0/
      data deltaM/ 16 * 0.d0/, neutmass/ 16 * 0.d0/, delta1/ 16 * 0.d0/
      data delta2/ 16 * 0.d0/,deltatot/ 16 * 0.d0/,tdeltatot/ 16 * 0.d0/
      double precision aPsiPsiHc1(4, 2), bPsiPsiHc1(4, 2)
      data aPsiPsiHc1/ 8 * 0.d0/, bPsiPsiHc1/ 8 * 0.d0/
      double precision aPsiPsiHc2(4, 2), bPsiPsiHc2(4, 2)
      data aPsiPsiHc2/ 8 * 0.d0/, bPsiPsiHc2/ 8 * 0.d0/
      double precision aPsiChiHc1(4, 2), bPsiChiHc1(4, 2)
      data aPsiChiHc1/ 8 * 0.d0/, bPsiChiHc1/ 8 * 0.d0/
      double precision aPsiChiHc2(4, 2), bPsiChiHc2(4, 2)
      data aPsiChiHc2/ 8 * 0.d0/, bPsiChiHc2/ 8 * 0.d0/
      double precision aPsiChiHHp(4, 2), bPsiChiHHp(4, 2)
      data aPsiChiHHp/ 8 * 0.d0/, bPsiChiHHp/ 8 * 0.d0/
      double precision aPsiChiHGp(4, 2), bPsiChiHGp(4, 2)
      data aPsiChiHGp/ 8 * 0.d0/, bPsiChiHGp/ 8 * 0.d0/
      double precision b1mGpmchar(2), b0mGpmchar(2)
      double precision b1mHpmchar(2), b0mHpmchar(2)
      double precision aPsiPsis1(4, 4), aPsiPsis2(4, 4)
      data aPsiPsis1/ 16 * 0.d0/, aPsiPsis2/ 16 * 0.d0/ 
      double precision aPsiPsip1(4, 4), aPsiPsip2(4, 4)
      data aPsiPsip1/ 16 * 0.d0/, aPsiPsip2/ 16 * 0.d0/
      double precision bPsiPsis1(4, 4), bPsiPsis2(4, 4) 
      data bPsiPsis1/ 16 * 0.d0/, bPsiPsis2/ 16 * 0.d0/    
      double precision bPsiPsip1(4, 4), bPsiPsip2(4, 4)
      data bPsiPsip1/ 16 * 0.d0/, bPsiPsip2/ 16 * 0.d0/
      double precision aPsiChis1(4, 4), aPsiChis2(4, 4)
      data aPsiChis1/ 16 * 0.d0/, aPsiChis2/ 16 * 0.d0/
      double precision aPsiChip1(4, 4), aPsiChip2(4, 4)
      data aPsiChip1/ 16 * 0.d0/, aPsiChip2/ 16 * 0.d0/
      double precision bPsiChis1(4, 4), bPsiChis2(4, 4)
      data bPsiChis1/ 16 * 0.d0/, bPsiChis2/ 16 * 0.d0/
      double precision bPsiChip1(4, 4), bPsiChip2(4, 4)
      data bPsiChip1/ 16 * 0.d0/, bPsiChip2/ 16 * 0.d0/
      double precision aPsiChiHu(4, 4), aPsiChih(4, 4), aPsiChiG(4, 4)
      data aPsiChiHu/ 16 * 0.d0/, aPsiChih/ 16 * 0.d0/
      data aPsiChiG/ 16 * 0.d0/
      double precision aPsiChiA(4, 4), bPsiChiHu(4, 4), bPsiChih(4, 4)
      data aPsiChiA/ 16 * 0.d0/, bPsiChiHu/ 16 * 0.d0/
      data bPsiChih/ 16 * 0.d0/
      double precision bPsiChiG(4, 4), bPsiChiA(4, 4)
      data bPsiChiG/ 16 * 0.d0/, bPsiChiA/ 16 * 0.d0/
      double precision b1mHu0mneut(4),b0mHu0mneut(4)
      double precision b1mh0mneut(4),b0mh0mneut(4)
      double precision b1mG0mneut(4),b0mG0mneut(4)
      double precision b1mA0mneut(4),b0mA0mneut(4)

      data b0msup/ 2 * 0.d0/,b1msup/ 2 * 0.d0/,b1msdown/ 2 * 0.d0/
      data b1msstrange/ 2 * 0.d0/,b1msel/ 2 * 0.d0/,b1mscharm/ 2 * 0.d0/
      data b1msmu/ 2 * 0.d0/,b1mstop/ 2 * 0.d0/,b1msbot/ 2 * 0.d0/
      data b1mstau/ 2 * 0.d0/,b0mstop/ 2 * 0.d0/,b0msbot/ 2 * 0.d0/
      data b0mstau/ 2 * 0.d0/,b1mchar/ 2 * 0.d0/,b0mchar/ 2 * 0.d0/
      data b1mGpmchar/ 2 * 0.d0/,b0mGpmchar/ 2 * 0.d0/
      data b1mHpmchar/ 2 * 0.d0/,b0mHpmchar/ 2 * 0.d0/

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole
      
      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

!-------------------------------------------------------------------------------------
      external dag2d,dag4d,mat3prod2d,a0,theta,b0,f,funcg,rmat2d,b1
      external dsyev

      include 'stdinputs.h' 

!--------------------------------------------------------------------------------
!    Nomenclature
!-------------------------------------------

      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      MZ = MZpole 

      msup(1)    = dsqrt(SUegg(6)) 
      msup(2)    = dsqrt(SUegg(5))
      mscharm(1) = dsqrt(SUegg(4))
      mscharm(2) = dsqrt(SUegg(3))
      mstop(1)   = dsqrt(SUegg(2))
      mstop(2)   = dsqrt(SUegg(1))
      
      msdown(1)    = dsqrt(SDegg(6)) 
      msdown(2)    = dsqrt(SDegg(5))
      msstrange(1) = dsqrt(SDegg(4))
      msstrange(2) = dsqrt(SDegg(3))
      msbot(1)     = dsqrt(SDegg(2))
      msbot(2)     = dsqrt(SDegg(1))

      msel(1)   = dsqrt(SLegg(6))
      msel(2)   = dsqrt(SLegg(5))
      msmu(1)   = dsqrt(SLegg(4))
      msmu(2)   = dsqrt(SLegg(3))
      mstau(1)  = dsqrt(SLegg(2))
      mstau(2)  = dsqrt(SLegg(1))

c$$$      msup(1)    = dsqrt(SUegg(5)) 
c$$$      msup(2)    = dsqrt(SUegg(6))
c$$$      mscharm(1) = dsqrt(SUegg(3))
c$$$      mscharm(2) = dsqrt(SUegg(4))
c$$$      mstop(1)   = dsqrt(SUegg(1))
c$$$      mstop(2)   = dsqrt(SUegg(2))
c$$$      
c$$$      msdown(1)    = dsqrt(SDegg(5)) 
c$$$      msdown(2)    = dsqrt(SDegg(6))
c$$$      msstrange(1) = dsqrt(SDegg(3))
c$$$      msstrange(2) = dsqrt(SDegg(4))
c$$$      msbot(1)     = dsqrt(SDegg(1))
c$$$      msbot(2)     = dsqrt(SDegg(2))
c$$$
c$$$      msel(1)   = dsqrt(SLegg(5))
c$$$      msel(2)   = dsqrt(SLegg(6))
c$$$      msmu(1)   = dsqrt(SLegg(3))
c$$$      msmu(2)   = dsqrt(SLegg(4))
c$$$      mstau(1)  = dsqrt(SLegg(1))
c$$$      mstau(2)  = dsqrt(SLegg(2))

       
      msnu(1)  = dsqrt(SNegg(1))
      msnu(2)  = dsqrt(SNegg(2))
      msnu(3)  = dsqrt(SNegg(3))

      mneut(1) = Neg(1) 
      mneut(2) = Neg(2)
      mneut(3) = Neg(3)
      mneut(4) = Neg(4)

      mneutm(1) = dabs(Neg(1)) 
      mneutm(2) = dabs(Neg(2))
      mneutm(3) = dabs(Neg(3))
      mneutm(4) = dabs(Neg(4))

      mchargino(1) = Ceg(1)
      mchargino(2) = Ceg(2)

      mh0  = dsqrt((mh0sq))
      mHu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      mG0 = MZ
      mGp = MW
!----------------------------------------------
!    General Definitions.
!----------------------------------------------
      
      beta     = datan(tanbeta)
      cosbeta  = dcos(beta)
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0*datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      tan2beta = dtan(2.d0*datan(tanbeta))

      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      sinsqthw = sinsqthw_susy
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

      ht   = yuRG(3,3)
      hb   = ydRG(3,3)
      htau = yeRG(3,3)

      htsq =   ht*ht
      hbsq =   hb*hb
      htausq = htau*htau
!-------------------------------------

      call dag2d(OCL,OCLdag)
      call dag2d(OCR,OCRdag)
!----------------------------------------------------------------
C     fermion-sfermion couplings 
!----------------------------------------------------------------



      apsiuq(1, 2) = (gp / dsqrt(2.d0)) * yuR
      bpsiuq(1, 1) = (gp / dsqrt(2.d0)) * yuL

      bpsiuq(2, 1) = g * dsqrt(2.d0) * 0.5d0

      apsitq(1, 2) = (gp / dsqrt(2.d0)) * yuR
      bpsitq(1, 1) = (gp / dsqrt(2.d0)) * yuL
      bpsitq(2, 1) = g * dsqrt(2.d0) * 0.5d0
      apsitq(4, 1) = ht
      bpsitq(4, 2) = ht
      
      apsidq(1, 2) = (gp / dsqrt(2.d0)) * ydR
      bpsidq(1, 1) = (gp / dsqrt(2.d0)) * ydL

      bpsidq(2, 1) = -g * dsqrt(2.d0) * 0.5d0

      apsibq(1, 2) = (gp / dsqrt(2.d0)) * ydR
      bpsibq(1, 1) = (gp / dsqrt(2.d0)) * ydL
      bpsibq(2, 1) = -g * dsqrt(2.d0) * 0.5d0
      apsibq(3, 1) = hb
      bpsibq(3, 2) = hb

      apsiel(1, 2) = (gp / dsqrt(2.d0)) * yeR
      bpsiel(1, 1) = (gp / dsqrt(2.d0)) * yeL

      bpsiel(2, 1) = -g * dsqrt(2.d0) * 0.5d0

      apsitaul(1, 2) = (gp / dsqrt(2.d0)) * yeR
      bpsitaul(1, 1) = (gp / dsqrt(2.d0)) * yeL
      bpsitaul(2, 1) = -g * dsqrt(2.d0) * 0.5d0
      apsitaul(3, 1) = htau
      bpsitaul(3, 2) = htau


      bpsinul(1, 1) = (gp / dsqrt(2.d0)) * ynuL
      bpsinul(2, 1) = g * dsqrt(2.d0) * 0.5d0

      bPsiNutnut(1, 1) = (gp / dsqrt(2.d0)) * ynuL
      bPsiNutnut(2, 1) = g * dsqrt(2.d0) * 0.5d0

!--------------------------------------------------------------------
C     /// mixing-  third family sfermions
!--------------------------------------------------------------------

      call rmat2d(thetat,rthetat)
      call rmat2d(thetab,rthetab)
      call rmat2d(thetatau,rthetatau)


      looptht: do i = 1, 4

      t1(1) = rthetat(1,1)*apsitq(i,1) + rthetat(1,2)*apsitq(i,2)   
      t1(2) = rthetat(2,1)*apsitq(i,1) + rthetat(2,2)*apsitq(i,2)   
      t2(1) = rthetat(1,1)*bpsitq(i,1) + rthetat(1,2)*bpsitq(i,2)   
      t2(2) = rthetat(2,1)*bpsitq(i,1) + rthetat(2,2)*bpsitq(i,2)   

      apsitq(i, 1) = t1(1)
      apsitq(i, 2) = t1(2) 
      bpsitq(i, 1) = t2(1)
      bpsitq(i, 2) = t2(2) 
      
      enddo looptht
    
      loopthb: do i = 1, 4

      t1(1) = rthetab(1,1)*apsibq(i,1) + rthetab(1,2)*apsibq(i,2)   
      t1(2) = rthetab(2,1)*apsibq(i,1) + rthetab(2,2)*apsibq(i,2)   
      t2(1) = rthetab(1,1)*bpsibq(i,1) + rthetab(1,2)*bpsibq(i,2)   
      t2(2) = rthetab(2,1)*bpsibq(i,1) + rthetab(2,2)*bpsibq(i,2)   

      apsibq(i, 1) = t1(1)
      apsibq(i, 2) = t1(2) 
      bpsibq(i, 1) = t2(1)
      bpsibq(i, 2) = t2(2) 
      
      enddo loopthb


      loopthtau: do i = 1, 4

      t1(1) = rthetatau(1,1)*apsitaul(i,1) + 
     $     rthetab(1,2)*apsitaul(i,2)   

      t1(2) = rthetatau(2,1)*apsitaul(i,1) + 
     $     rthetatau(2,2)*apsitaul(i,2)   

      t2(1) = rthetatau(1,1)*bpsitaul(i,1) + 
     $     rthetatau(1,2)*bpsitaul(i,2)   

      t2(2) = rthetatau(2,1)*bpsitaul(i,1) + 
     $     rthetatau(2,2)*bpsitaul(i,2)   

      apsitaul(i, 1) = t1(1)
      apsitaul(i, 2) = t1(2) 
      bpsitaul(i, 1) = t2(1)
      bpsitaul(i, 2) = t2(2) 
      
      enddo loopthtau

!----------------------------------------------------------
C    Neutralino-Higgs couplings
!----------------------------------------------------------    
      aPsi0PsicW(2, 1) = - g
      bPsi0PsicW(2, 1) = - g
      aPsi0PsicW(4, 2) =   g / dsqrt(2.d0)		     
      bPsi0PsicW(3, 2) = - g / dsqrt(2.d0)	     

      loop4by2i: DO i = 1, 4
      loop4by2j: DO j = 1, 2
      aPsi0ChicW(i,j) = 0.d0
      bPsi0ChicW(i,j) = 0.d0
      loop4by2k: DO k = 1, 2
      
      aPsi0ChicW(i,j) = aPsi0ChicW(i,j) + aPsi0PsicW(i,k)*OCRdag(k,j)
      bPsi0ChicW(i,j) = bPsi0ChicW(i,j) + bPsi0PsicW(i,k)*OCLdag(k,j)
      
      ENDDO loop4by2k
      ENDDO loop4by2j

      ENDDO loop4by2i
      

      loop1: do i = 1, 4
      loop2: do j = 1, 4
      
      aPsiPsiZ(i, j) = 0.d0
      bPsiPsiZ(i, j) = 0.d0
      
      enddo loop2
      enddo loop1


      aPsiPsiZ(3, 3) =   g / (2.d0 * costhW)
      aPsiPsiZ(4, 4) = - g / (2.d0 * costhW)
      
      bPsiPsiZ(3, 3) = -1.d0 * aPsiPsiZ(3, 3)
      bPsiPsiZ(4, 4) = -1.d0 * aPsiPsiZ(4, 4)

      call dag4d(ON,ONdag)
      call matmult4d(aPsiPsiZ,ONdag,aPsiChiZ)
      call matmult4d(bPsiPsiZ,ONdag,bPsiChiZ)
      
      aPsiPsiHc1(1, 2) =  gp / dsqrt(2.d0)
      bPsiPsiHc2(1, 2) =  aPsiPsiHc1(1, 2)
      aPsiPsiHc1(2, 2) =  g / dsqrt(2.d0)
      bPsiPsiHc2(2, 2) =  g / dsqrt(2.d0)
      aPsiPsiHc1(3, 1) = -g
      bPsiPsiHc2(4, 1) =  g
      
      loop5by2i: DO i = 1, 4
      loop5by2j: DO j = 1, 2
      loop5by2k: DO k = 1, 2
      
      aPsiChiHc1(i,j) = aPsiChiHc1(i,j) + aPsiPsiHc1(i,k)*OCLdag(k,j)
      bPsiChiHc1(i,j) = bPsiChiHc1(i,j) + bPsiPsiHc1(i,k)*OCRdag(k,j)
      aPsiChiHc2(i,j) = aPsiChiHc2(i,j) + aPsiPsiHc2(i,k)*OCLdag(k,j)
      bPsiChiHc2(i,j) = bPsiChiHc2(i,j) + bPsiPsiHc2(i,k)*OCRdag(k,j)
      
      ENDDO loop5by2k
      ENDDO loop5by2j
      ENDDO loop5by2i


      loop6: do i = 1, 4
      loop7: do j = 1, 2

      aPsiChiHGp(i, j) =   cosbeta * aPsiChiHc1(i, j) + 
     $     sinbeta * aPsiChiHc2(i, j)

      aPsiChiHHp(i, j) = - sinbeta * aPsiChiHc1(i, j) + 
     $     cosbeta * aPsiChiHc2(i, j)
      
      bPsiChiHGp(i, j) =   cosbeta * bPsiChiHc1(i, j) + 
     $     sinbeta * bPsiChiHc2(i, j)
      
      bPsiChiHHp(i, j) = - sinbeta * bPsiChiHc1(i, j) + 
     $     cosbeta * bPsiChiHc2(i, j)
      
      enddo loop7
      enddo loop6
      

      aPsiPsis1(1, 3) =  - gp * 0.5d0
      aPsiPsis1(3, 1) =    aPsiPsis1(1, 3)
      aPsiPsis1(2, 3) =    g * 0.5d0
      aPsiPsis1(3, 2) =    aPsiPsis1(2, 3) 
      aPsiPsis2(2, 4) =  - g * 0.5d0
      aPsiPsis2(4, 2) =    aPsiPsis2(2, 4)
      aPsiPsis2(1, 4) =    gp * 0.5d0
      aPsiPsis2(4, 1) =    aPsiPsis2(1, 4) 
      aPsiPsip1(1, 3) =  - gp * 0.5d0
      aPsiPsip1(3, 1) =    aPsiPsip1(1, 3) 
      aPsiPsip1(2, 3) =    g * 0.5d0
      aPsiPsip1(3, 2) =    aPsiPsip1(2, 3) 
      aPsiPsip2(2, 4) =    g * 0.5d0
      aPsiPsip2(4, 2) =    aPsiPsip2(2, 4) 
      aPsiPsip2(1, 4) =  - gp * 0.5d0
      aPsiPsip2(4, 1) =    aPsiPsip2(1, 4) 

      loop9:  do i = 1, 4
      loop10: do j = 1, 4
      
      bPsiPsis1(i, j) =   aPsiPsis1(i, j)
      bPsiPsis2(i, j) =   aPsiPsis2(i, j)
      bPsiPsip1(i, j) = - aPsiPsip1(i, j)
      bPsiPsip2(i, j) = - aPsiPsip2(i, j)

      enddo loop10
      enddo loop9

      call matmult4d(aPsiPsis1,ONdag,aPsiChis1)
      call matmult4d(aPsiPsis2,ONdag,aPsiChis2)
      call matmult4d(aPsiPsip1,ONdag,aPsiChip1)
      call matmult4d(aPsiPsip2,ONdag,aPsiChip2)
      call matmult4d(bPsiPsis1,ONdag,bPsiChis1)
      call matmult4d(bPsiPsis2,ONdag,bPsiChis2)
      call matmult4d(bPsiPsip1,ONdag,bPsiChip1)
      call matmult4d(bPsiPsip2,ONdag,bPsiChip2)

      
      loop11: do i = 1, 4
      loop12: do j = 1, 4

      aPsiChiHu(i,j) =   cos(alpha) * aPsiChis1(i,j) + 
     $     sin(alpha) * aPsiChis2(i,j)

      aPsiChih(i,j) = - sin(alpha) * aPsiChis1(i,j) + 
     $     cos(alpha) * aPsiChis2(i,j)

      aPsiChiG(i,j) =   cosbeta * aPsiChip1(i,j) + 
     $     sinbeta * aPsiChip2(i,j)
      
      aPsiChiA(i,j) = - sinbeta * aPsiChip1(i,j) + 
     $     cosbeta * aPsiChip2(i,j)


      bPsiChiHu(i,j) =   cos(alpha) * bPsiChis1(i,j) + 
     $     sin(alpha) * bPsiChis2(i,j)

      bPsiChih(i,j) = - sin(alpha) * bPsiChis1(i,j) + 
     $     cos(alpha) * bPsiChis2(i,j)
      
      bPsiChiG(i,j) =   cosbeta * bPsiChip1(i,j) + 
     $     sinbeta * bPsiChip2(i,j)
      
      bPsiChiA(i,j) = - sinbeta * bPsiChip1(i,j) + 
     $     cosbeta * bPsiChip2(i,j)

      enddo loop12
      enddo loop11

!-------------------------------------------------------------------
C     Computing PV functions for one loop correction to neutralino
!-------------------------------------------------------------------


      loopb1: do i = 1, 2

      call b0(p,muq,msup(i),q,b0msup(i))
      call b1(p,muq,msup(i),q,b1msup(i))

      call b1(p,md,msdown(i),q,b1msdown(i))
      call b1(p,me,msel(i),q,b1msel(i))

      call b1(p,mc,mscharm(i),q,b1mscharm(i))
      call b1(p,ms,msstrange(i),q,b1msstrange(i))
      call b1(p,mmu,msmu(i),q,b1msmu(i))
      
      call b1(p,mt,mstop(i),q,b1mstop(i))
      call b1(p,mb,msbot(i),q,b1msbot(i))
      call b1(p,mtau,mstau(i),q,b1mstau(i))

      call b0(p,mt,mstop(i),q,b0mstop(i))
      call b0(p,mb,msbot(i),q,b0msbot(i))
      call b0(p,mtau,mstau(i),q,b0mstau(i))

      call b1(p, mchargino(i), MW, q, b1mchar(i))
      call b0(p, mchargino(i), MW, q, b0mchar(i))

      call b1(p, mchargino(i), mGp, q, b1mGpmchar(i))
      call b0(p, mchargino(i), mGp, q, b0mGpmchar(i))

      call b1(p, mchargino(i), mHpm, q, b1mHpmchar(i))
      call b0(p, mchargino(i), mHpm, q, b0mHpmchar(i))

      enddo loopb1

      loopb14: do i = 1, 4

      call b1(p, mneutm(i), MZ, q, b1mneut(i))
      call b0(p, mneutm(i), MZ, q, b0mneut(i))
    
      call b1(p, mneutm(i), mHu0, q, b1mHu0mneut(i))
      call b0(p, mneutm(i), mHu0, q, b0mHu0mneut(i))

      call b1(p, mneutm(i), mh0, q, b1mh0mneut(i))
      call b0(p, mneutm(i), mh0, q, b0mh0mneut(i))

      call b1(p, mneutm(i), mG0, q, b1mG0mneut(i))
      call b0(p, mneutm(i), mG0, q, b0mG0mneut(i))

      call b1(p, mneutm(i), mA0, q, b1mA0mneut(i))
      call b0(p, mneutm(i), mA0, q, b0mA0mneut(i))

      enddo loopb14

      
      call b1(p,2.d-5,msnu(3),q,b1msnue)
      call b1(p,2.d-5,msnu(2),q,b1msnumu)
      call b1(p,2.d-5,msnu(1),q,b1msnutau)

!-------------------------------------------------------------------
C    Corrections begin
!-------------------------------------------------------------------
C    fermion-sfermion contributions
!------------------------------------------------------------------
      loopsig1i: do i = 1, 4
      loopsig1j: do j = 1, 4
            sigmaL(i,j) = 0.d0
            sigmaR(i,j) = 0.d0
            sigmaS(i,j) = 0.d0
      loopsig1k: do k = 1, 2

      sigmaL(i, j) = sigmaL(i, j) + 
     $     (3.d0 * apsiuq(i, k) * apsiuq(j, k) * real(b1msup(k)) +
     $     3.d0 * apsidq(i, k) * apsidq(j, k) * real(b1msdown(k)) +
     $     apsiel(i, k) * apsiel(j, k) * real(b1msel(k)) +
     $     apsinul(i, k) * apsinul(j, k) * real(b1msnue))

      sigmaL(i, j) = sigmaL(i, j) + 
     $     (3.d0 * apsiuq(i, k) * apsiuq(j, k) * real(b1mscharm(k))  +
     $     3.d0 * apsidq(i, k) * apsidq(j, k) * real(b1msstrange(k)) +
     $     apsiel(i, k) * apsiel(j, k) * real(b1msmu(k)) +
     $     apsinul(i, k) * apsinul(j, k) * real(b1msnumu))

      sigmaL(i, j) = sigmaL(i, j) + 
     $     (3.d0 * apsitq(i, k) * apsitq(j, k) * real(b1mstop(k)) +
     $     3.d0 * apsibq(i, k) * apsibq(j, k) * real(b1msbot(k)) +
     $     apsitaul(i, k) * apsitaul(j, k) * real(b1mstau(k)) +
     $     apsinutl(i, k) * apsinutl(j, k) * real(b1msnutau))

      
      sigmaR(i, j) = sigmaR(i, j) + 
     $     (3.d0 * bpsiuq(i, k) * bpsiuq(j, k) * real(b1msup(k))  +
     $     3.d0 * bpsidq(i, k) * bpsidq(j, k) * real(b1msdown(k)) +
     $     bpsiel(i, k) * bpsiel(j, k) * real(b1msel(k)) +
     $     bpsinul(i, k) * bpsinul(j, k) * real(b1msnue)) 

      sigmaR(i, j) = sigmaR(i, j) + 
     $     (3.d0 * bpsiuq(i, k) * bpsiuq(j, k) * real(b1mscharm(k)) +
     $     3.d0 * bpsidq(i, k) * bpsidq(j, k) * real(b1msstrange(k)) +
     $     bpsiel(i, k) * bpsiel(j, k) * real(b1msmu(k)) +
     $     bpsinul(i, k) * bpsinul(j, k) * real(b1msnumu))

      sigmaR(i, j) = sigmaR(i, j) + 
     $     (3.d0 * bpsitq(i, k) * bpsitq(j, k) * real(b1mstop(k)) +
     $     3.d0 * bpsibq(i, k) * bpsibq(j, k) * real(b1msbot(k)) +
     $     bpsitaul(i, k) * bpsitaul(j, k) * real(b1mstau(k)) +
     $     bPsiNutnut(i, k) * bPsiNutnut(j, k) * real(b1msnutau))

      sigmaS(i, j) = sigmaS(i, j) + 2.d0 * 
     $     (3.d0 * bpsitq(i, k) * apsitq(j, k) * mt * real(b0mstop(k)) +
     $     3.d0 * bpsibq(i, k) * apsibq(j, k) * mb * real(b0msbot(k))  +
     $    bpsitaul(i, k) * apsitaul(j, k) * mtau * real(b0mstau(k)))


      enddo loopsig1k
      enddo loopsig1j
      enddo loopsig1i
      
!-----------------------------------------------------------------
!     Higgs contribution
!------------------------------------------------------------------
      loopsig2i: do i = 1, 4
      loopsig2j: do j = 1, 4
      loopsig2k: do k = 1, 2

      sigmaL(i, j) = sigmaL(i, j) + 2.d0 * 
     $     (aPsi0ChicW(i, k) * aPsi0ChicW(j, k) * b1mchar(k))

      sigmaR(i, j) = sigmaR(i, j) + 2.d0 * 
     $     real(bPsi0ChicW(i, k) * bPsi0ChicW(j, k) * b1mchar(k))

      sigmaS(i, j) = sigmaS(i, j) - 8.d0 * mchargino(k) * 
     $     real(bPsi0ChicW(i, k) * aPsi0ChicW(j, k) * b0mchar(k))

      enddo loopsig2k
      enddo loopsig2j
      enddo loopsig2i 


      loopsig3i: do i = 1, 4
      loopsig3j: do j = 1, 4
      loopsig3k: do k = 1, 4
      
      sigmaL(i, j) = sigmaL(i, j) +
     $     real(aPsiChiZ(i, k) * aPsiChiZ(j, k) * b1mneut(k))

      sigmaR(i, j) = sigmaR(i, j) +
     $     real(bPsiChiZ(i, k) * bPsiChiZ(j, k) * b1mneut(k))

      sigmaS(i, j) = sigmaS(i, j) - 4.d0 * mneut(k) * 
     $     real(bPsiChiZ(i, k) * aPsiChiZ(j, k) * b0mneut(k))
      
      enddo loopsig3k
      enddo loopsig3j
      enddo loopsig3i 

!---------------------------------------------------------
      loopsig4i: do i = 1, 4
      loopsig4j: do j = 1, 4
      loopsig4k: do k = 1, 2
      
!     G+

      sigmaL(i, j) = sigmaL(i, j) + 
     $     real(aPsiChiHGp(i, k) * aPsiChiHGp(j, k) * b1mGpmchar(k))

      sigmaR(i, j) = sigmaR(i, j) + 
     $     real(bPsiChiHGp(i, k) * bPsiChiHGp(j, k) * b1mGpmchar(k))

      sigmaS(i, j) = sigmaS(i, j) + 2.d0 * mchargino(k) * 
     $     real(bPsiChiHGp(i, k) * aPsiChiHGp(j, k) * b0mGpmchar(k))
      
!------------------------------------------------------------------
!     H+

      sigmaL(i, j) = sigmaL(i, j) + 
     $     real(aPsiChiHHp(i, k) * aPsiChiHHp(j, k) * b1mHpmchar(k))

      sigmaR(i, j) = sigmaR(i, j) + 
     $     real(bPsiChiHHp(i, k) * bPsiChiHHp(j, k) * b1mHpmchar(k))

      sigmaS(i, j) = sigmaS(i, j) + 2.d0 * mchargino(k) * 
     $     real(bPsiChiHHp(i, k) * aPsiChiHHp(j, k) * b0mHpmchar(k))
      
      enddo loopsig4k
      enddo loopsig4j
      enddo loopsig4i

!--------------------------------------------------------------------  
      loopsig5i: do i = 1, 4
      loopsig5j: do j = 1, 4
      loopsig5k: do k = 1, 4

C     H


      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 *
     $     real(aPsiChiHu(i, k) * aPsiChiHu(j, k) * b1mHu0mneut(k))
      
      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 * 
     $     real(bPsiChiHu(i, k) * bPsiChiHu(j, k) * b1mHu0mneut(k))
      
      sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
     $     real(bPsiChiHu(i, k) * aPsiChiHu(j, k) * b0mHu0mneut(k))
 
!--------------------------------------------------------------------
C     h
      
      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 *
     $     real(aPsiChih(i, k) * aPsiChih(j, k) * b1mh0mneut(k))
      
      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 * 
     $     real(bPsiChih(i, k) * bPsiChih(j, k) * b1mh0mneut(k))
      
      sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
     $     real(bPsiChih(i, k) * aPsiChih(j, k) * b0mh0mneut(k))

!--------------------------------------------------------------------
C     G0
      
      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 *
     $     real(aPsiChiG(i, k) * aPsiChiG(j, k) * b1mG0mneut(k))
      
      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 * 
     $     real(bPsiChiG(i, k) * bPsiChiG(j, k) * b1mG0mneut(k))
      
      sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
     $     real(bPsiChiG(i, k) * aPsiChiG(j, k) * b0mG0mneut(k))
      
!--------------------------------------------------------------------
C     A0
      
      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 *
     $     real(aPsiChiA(i, k) * aPsiChiA(j, k) * b1mA0mneut(k))
      
      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 * 
     $     real(bPsiChiA(i, k) * bPsiChiA(j, k) * b1mA0mneut(k))
      
      sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
     $     real(bPsiChiA(i, k) * aPsiChiA(j, k) * b0mA0mneut(k))
      
      
      enddo loopsig5k
      enddo loopsig5j
      enddo loopsig5i
!-------------------------------------------------------------------      
      
      loopneui: do i = 1, 4
      loopneuj: do j = 1, 4
      
      neutmass(i,j) = MNeut1(i,j)
      
      enddo loopneuj
      enddo loopneui
      
      call matmult4d(sigmaR,neutmass,delta1)
      call matmult4d(neutmass,sigmaL,delta2)      
      call add4dmat(delta1,delta2,sigmaS,deltatot)
      call dag4d(deltatot,tdeltatot)
      call add4dmat2(deltatot,tdeltatot,deltaM)



      loopaddi: do i = 1, 4
      loopaddj: do j = 1, 4

      deltaM(i,j) = ( - deltaM(i,j)/(32.d0 * pi * pi))
      
      enddo loopaddj
      enddo loopaddi
      
      call add4dmat2(neutmass,deltaM,neutmasstot)


c$$$      print*,"one loop correction neutralino---------------"
c$$$
c$$$ 991  format( '(', f10.5, f10.5, f10.5, f10.5  ')')
c$$$c$$$  write(*,*) '********* s-electrons 6x6 matrixin test ***********'
c$$$      verbsoft: do a=1,4
c$$$      write(*,991) neutmasstot(a,1),neutmasstot(a,2),
c$$$     $     neutmasstot(a,3),neutmasstot(a,4)
c$$$      
c$$$      enddo verbsoft

      
      
      do a = 1, 4
         do m = 1, 4 

            neutmass1(a,m) = 0.d0          
            neutmass1(a,m) = neutmasstot(a,m)
            
         enddo
      enddo
 
!----------------------------------------------------------------------
C Computes eigenvalues for the given 4x4 nonsymmetric matrix neutmass1
!-----------------------------------------------------------------------      


!       call dgesvd('A','A',4,4,neutmass1,4,Negm,ONLm,4,ONRm,4,
!     $     work,lwork,info)

      call SVD(4, 4, neutmass1,4, Negm, ONRm,4, ONew,4, 1)

!      Call CEigensystem(4,neutmass1,4,Negm,ONew,4,1)

!      call DGEEV('V', 'V', 4, neutmass1, 4, Negm, Neuevim, ONLm, 4, 
!     $     ONRm, 4, WORK, LWORK, INFO)

c$$$      do a = 1, 4
c$$$         do m = 1, 4 
c$$$
c$$$            ONew(a,m) = 0.d0          
c$$$            ONew(a,m) = ONRm(a,m)
c$$$            
c$$$         enddo
c$$$      enddo
      

c$$$      if(info.eq.0) then
c$$$         AOK = AOK + 1
c$$$      endif
!----------------------------------------------------------------------
          
      return
      end subroutine neutralino
      
!=======================================================================================

