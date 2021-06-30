****f* SuSeFLAV/oneloopchargino.f 
*  NAME
*    subroutine chargino
*  SYNOPSIS
*    One loop correction to chargino. 
*  FUNCTION
*     Computes 1loop correction to charginos at a given energy scale and 
*     external momenta
*
*  INPUTS
*     p                              - External momentum  
*     q                              - Energy scale
*     g,gp,g3                        - Gauge couplings(g = em, gp = weak, 
*                                      g3 = strong)                                     
*     mt,mb,mtau                     - pole masses of top, bottom and tau
*     yuRG,ydRG,yeRG                 - (3 X 3)Yukawas
*     mh0sq,mhu0sq,mhpmsq,mA0sq      - physicsal higgs mass squared 
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
*     Ceg     =   2 singular values of the chargino Mass Matrix
*
*  RESULT
*     Cegm        =  2 eigenvalues of the 1-loop corrected chargino mass matrix.
*     OCRm, OCLTm =  (2 X 2) orthogonal matrices such that 
*                 MChar = Transpose[OCRm].Diag[MChar].OCLTm
*     charmasstot = (2x2) corrected chargino mass matrix
*
*  EXAMPLE
*
*      subroutine chargino(p,q,g,gp,mt,mb,mtau,tanbeta,yuRG,ydRG,yeRG,
*     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,MChar,ON,OCL,OCR,mh0sq,
*     $     mhu0sq,mHpmsq,mA0sq,charmasstot,OCRm,OCLTm,Cegm)
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
!---------------------------------------------------------------------------------    
      subroutine chargino(p,q,g,gp,mt,mb,mtau,tanbeta,yuRG,ydRG,yeRG,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,MChar,ON,OCL,OCR,mh0sq,
     $     mhu0sq,mHpmsq,mA0sq,charmasstot,OCRm,OCLm,Cegm)
      
      implicit none
      
      integer lwork1,info
      parameter(lwork1=15)
      double precision work1(lwork1)
      
      integer i,j,k,a,m
      double precision p,q,g,gp !,staus
      double precision mT, mB, mTau     
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR   
      double precision ON(4,4),OCL(2,2),OCR(2,2),OCLdag(2,2)
      double precision OCRdag(2,2),MChar(2,2)
      double precision ht,htsq,hb,htau,hbsq,htausq
      double precision mh0,mHpm,mhu0,mA0
      double precision sinbeta
      double precision mG0,mGp     
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
    
      double precision yuRG(3,3),ydRG(3,3),yeRG(3,3)      
      double precision beta,cos2beta,sin2beta
      double precision aPsi0PsicW(4,2),bPsi0PsicW(4,2),aChi0PsicW(4,2)
      double precision bChi0PsicW(4, 2)
      data aPsi0PsicW/ 8 * 0.d0/,bPsi0PsicW/ 8 * 0.d0/
      data aChi0PsicW/ 8 * 0.d0/, bChi0PsicW/ 8 * 0.d0/
      double precision aPsiPsiZ(2, 2), bPsiPsiZ(2, 2), aPsiChiZ(2, 2)
      double precision bPsiChiZ(2, 2)
      data aPsiPsiZ/ 4 * 0.d0/, bPsiPsiZ/ 4 * 0.d0/
      data aPsiChiZ/ 4 * 0.d0/, bPsiChiZ/ 4 * 0.d0/
      double precision aPsiChiGam(2, 2), bPsiChiGam(2, 2) 
      data aPsiChiGam/ 4 * 0.d0/, bPsiChiGam/ 4 * 0.d0/

!-------------------------------------------------------------------------

      double precision b1mneut(4),b0mneut(4),e

      double precision t1(2),t2(2),sigmaL(2,2),sigmaR(2,2),sigmaS(2,2)
      
      double precision b0msup(2)
      
      double precision b1msup(2),b1msdown(2),b1msel(2),b1mscharm(2)
      double precision b1msstrange(2),b1msmu(2),b1mstop(2),b1msbot(2)
      double precision b1mstau(2),b0mstop(2),b0msbot(2),b0mstau(2)
      double precision b1msnue,b1msnumu,b1msnutau,b1mchar(2),b0mchar(2)
      double precision b1mchargm(2),b0mchargm(2)
      double precision deltaM(2,2),charmass(2,2),delta1(2,2),delta2(2,2)
      double precision charmasstot(2,2),deltatot(2,2),tdeltatot(2,2)

!      data deltatot/ 4 * 0.d0/, tdeltatot/ 4 * 0.d0/

      data deltaM/ 4 * 0.d0/, charmass/ 4 * 0.d0/, delta1/ 4 * 0.d0/
      data delta2/ 4 * 0.d0/,deltatot/ 4 * 0.d0/,tdeltatot/ 4 * 0.d0/
      double precision aPsiPsiHc1(4, 2), bPsiPsiHc1(4, 2)
      data aPsiPsiHc1/ 8 * 0.d0/, bPsiPsiHc1/ 8 * 0.d0/
      double precision aPsiPsiHc2(4, 2), bPsiPsiHc2(4, 2)
      data aPsiPsiHc2/ 8 * 0.d0/, bPsiPsiHc2/ 8 * 0.d0/
      double precision aChiPsiHc1(4, 2), bChiPsiHc1(4, 2)
      data aChiPsiHc1/ 8 * 0.d0/, bChiPsiHc1/ 8 * 0.d0/
      double precision aChiPsiHc2(4, 2), bChiPsiHc2(4, 2)
      data aChiPsiHc2/ 8 * 0.d0/, bChiPsiHc2/ 8 * 0.d0/
      double precision aChiPsiHHp(4, 2), bChiPsiHHp(4, 2)
      data aChiPsiHHp/ 8 * 0.d0/, bChiPsiHHp/ 8 * 0.d0/
      double precision aChiPsiHGp(4, 2), bChiPsiHGp(4, 2)
      data aChiPsiHGp/ 8 * 0.d0/, bChiPsiHGp/ 8 * 0.d0/

      double precision b1mGpmchar(2), b0mGpmchar(2)
      double precision b1mHpmchar(2), b0mHpmchar(2)

      double precision b1mGpmneut(4), b0mGpmneut(4)
      double precision b1mHpmneut(4), b0mHpmneut(4)

      double precision aPsiPsis1(2, 2), aPsiPsis2(2, 2)
      data aPsiPsis1/ 4 * 0.d0/, aPsiPsis2/ 4 * 0.d0/ 
      double precision aPsiPsip1(2, 2), aPsiPsip2(2, 2)
      data aPsiPsip1/ 4 * 0.d0/, aPsiPsip2/ 4 * 0.d0/
      double precision bPsiPsis1(2, 2), bPsiPsis2(2, 2) 
      data bPsiPsis1/ 4 * 0.d0/, bPsiPsis2/ 4 * 0.d0/    
      double precision bPsiPsip1(2, 2), bPsiPsip2(2, 2)
      data bPsiPsip1/ 4 * 0.d0/, bPsiPsip2/ 4 * 0.d0/
      double precision aPsiChis1(2, 2), aPsiChis2(2, 2)
      data aPsiChis1/ 4 * 0.d0/, aPsiChis2/ 4 * 0.d0/
      double precision aPsiChip1(2, 2), aPsiChip2(2, 2)
      data aPsiChip1/ 4 * 0.d0/, aPsiChip2/ 4 * 0.d0/
      double precision bPsiChis1(2, 2), bPsiChis2(2, 2)
      data bPsiChis1/ 4 * 0.d0/, bPsiChis2/ 4 * 0.d0/
      double precision bPsiChip1(2, 2), bPsiChip2(2, 2)
      data bPsiChip1/ 4 * 0.d0/, bPsiChip2/ 4 * 0.d0/
      double precision aPsiChiHu(2, 2), aPsiChih(2, 2), aPsiChiG(2, 2)
      data aPsiChiHu/ 4 * 0.d0/, aPsiChih/ 4 * 0.d0/
      data aPsiChiG/ 4 * 0.d0/
      double precision aPsiChiA(2, 2), bPsiChiHu(2, 2), bPsiChih(2, 2)
      data aPsiChiA/ 4 * 0.d0/, bPsiChiHu/ 4 * 0.d0/
      data bPsiChih/ 4 * 0.d0/
      double precision bPsiChiG(2, 2), bPsiChiA(2, 2)
      data bPsiChiG/ 4 * 0.d0/, bPsiChiA/ 4 * 0.d0/

      double precision b1mHu0mneut(4),b0mHu0mneut(4)
      double precision b1mh0mneut(4),b0mh0mneut(4)
      double precision b1mG0mneut(4),b0mG0mneut(4)
      double precision b1mA0mneut(4),b0mA0mneut(4)

      double precision b1mHu0mchar(2),b0mHu0mchar(2)
      double precision b1mh0mchar(2),b0mh0mchar(2)
      double precision b1mG0mchar(2),b0mG0mchar(2)
      double precision b1mA0mchar(2),b0mA0mchar(2)

!---------------------------------------------------------------------------

      double precision apsicdq(2, 2), apsicuq(2, 2), apsicenul(2, 2)
      double precision apsicnuel(2, 2)
      double precision apsicbtq(2, 2),apsictbq(2, 2),apsictaunul(2, 2)
      double precision apsicnutaul(2, 2)
      double precision apsicbtqr(2, 2), apsictbqr(2, 2), 
     $     apsictaunulr(2, 2)
      double precision apsicnutaulr(2, 2)
      
      data apsicdq/ 4 * 0.d0/, apsicuq/ 4 * 0.d0/, apsicenul/ 4 * 0.d0/
     $     apsicnuel/ 4 * 0.d0/, apsicbtq/ 4 * 0.d0/, apsictbq/4*0.d0/, 
     $     apsictaunul/ 4 * 0.d0/, apsicnutaul/ 4 * 0.d0/,
     $     apsicbtqr/ 4 * 0.d0/, apsictbqr/ 4 * 0.d0/, 
     $     apsictaunulr/ 4 * 0.d0/, apsicnutaulr/ 4 * 0.d0/
      
      double precision bpsicdq(2, 2), bPsicuq(2, 2), bpsicenul(2, 2)
      double precision bpsicnuel(2, 2)
      double precision bpsicbtq(2, 2),bpsictbq(2, 2),bpsictaunul(2, 2) 
      double precision bpsicnutaul(2, 2)
      double precision bpsicbtqr(2, 2), bpsictbqr(2, 2),
     $     bpsictaunulr(2, 2)
      double precision bpsicnutaulr(2, 2)

      data bpsicdq/ 4 * 0.d0/, bPsicuq/ 4 * 0.d0/, bpsicenul/ 4 * 0.d0/
     $     bpsicnuel/ 4 * 0.d0/, bpsicbtq/ 4 * 0.d0/, 
     $     bpsictbq/ 4 * 0.d0/, 
     $     bpsictaunul/ 4 * 0.d0/, bpsicnutaul/ 4 * 0.d0/,
     $     bpsicbtqr/ 4 * 0.d0/, bpsictbqr/ 4 * 0.d0/, 
     $     bpsictaunulr/ 4 * 0.d0/, bpsicnutaulr/ 4 * 0.d0/


      double precision OCRm(2,2),OCLTm(2,2),OCRTm(2,2),OCLm(2,2)
      
      DOUBLE PRECISION charmasstotT(2,2), thetaocr,thetaocl,
     $     sinsqthw_susy
      DOUBLE PRECISION mmdag(2,2),mdagm(2,2)
      double precision Cegm(2),trchar,detchar,Mchard(2,2)
      data trchar/ 0.d0/, detchar/ 0.d0/

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole
      
      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae
            
!-------------------------------------------------------------------------------------
      external dag2d,dag4d,mat3prod2d,a0,theta,b0,f,funcg,rmat2d,b1

      include 'stdinputs.h'
       
!--------------------------------------------------------------------------------
C     Nomenclature
C------------------------------------------

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

!--------------------------------------      
!   general definitions.
!--------------------------------------
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

C---------------------------------------------------------
C    fermion-sfermion couplings
C---------------------------------------------------------

      apsicdq(1, 1) = g       
      bPsicuq(1, 1) = g       

      apsicbtq(1, 1) = g       
      apsicbtq(2, 2) = - ht	     

      bpsicbtq(2, 1) = - hb    

      bpsictbq(1, 1) = g       
      bpsictbq(2, 2) = - hb    

      apsictbq(2, 1) = - ht

      apsicenul(1, 1)   = g
      bpsicnuel(1, 1)   = g

      apsictaunul(1, 1) = g

      bpsicnutaul(1, 1) = g
      bpsicnutaul(2, 2) = - htau

      bpsictaunul(2, 1) = - htau
      

!--------------------------------------------------------------------
C     mixing- third family sfermions
!--------------------------------------------------------------------

      call rmat2d(thetat,rthetat)
      call rmat2d(thetab,rthetab)
      call rmat2d(thetatau,rthetatau)

C     stops

      looptht: do i = 1, 2

      t1(1) = rthetat(1,1)*apsicbtq(i,1) + rthetat(1,2)*apsicbtq(i,2)   
      t1(2) = rthetat(2,1)*apsicbtq(i,1) + rthetat(2,2)*apsicbtq(i,2)   
      t2(1) = rthetat(1,1)*bpsicbtq(i,1) + rthetat(1,2)*bpsicbtq(i,2)   
      t2(2) = rthetat(2,1)*bpsicbtq(i,1) + rthetat(2,2)*bpsicbtq(i,2)   

      apsicbtqr(i, 1) = t1(1)
      apsicbtqr(i, 2) = t1(2) 
      bpsicbtqr(i, 1) = t2(1)
      bpsicbtqr(i, 2) = t2(2) 
      
      enddo looptht
   
C     sbots

      loopthb: do i = 1, 2

      t1(1) = rthetab(1,1)*apsictbq(i,1) + rthetab(1,2)*apsictbq(i,2)   
      t1(2) = rthetab(2,1)*apsictbq(i,1) + rthetab(2,2)*apsictbq(i,2)   
      t2(1) = rthetab(1,1)*bpsictbq(i,1) + rthetab(1,2)*bpsictbq(i,2)   
      t2(2) = rthetab(2,1)*bpsictbq(i,1) + rthetab(2,2)*bpsictbq(i,2)   

      apsictbqr(i, 1) = t1(1)
      apsictbqr(i, 2) = t1(2) 
      bpsictbqr(i, 1) = t2(1)
      bpsictbqr(i, 2) = t2(2) 
      
      enddo loopthb

C     staus

      loopthtau: do i = 1, 2

      t1(1) = rthetatau(1,1)*apsicnutaul(i,1) + 
     $     rthetab(1,2)*apsicnutaul(i,2)   

      t1(2) = rthetatau(2,1)*apsicnutaul(i,1) + 
     $     rthetatau(2,2)*apsicnutaul(i,2)   

      t2(1) = rthetatau(1,1)*bpsicnutaul(i,1) + 
     $     rthetatau(1,2)*bpsicnutaul(i,2)   

      t2(2) = rthetatau(2,1)*bpsicnutaul(i,1) + 
     $     rthetatau(2,2)*bpsicnutaul(i,2)   

      apsicnutaulr(i, 1) = t1(1)
      apsicnutaulr(i, 2) = t1(2) 
      bpsicnutaulr(i, 1) = t2(1)
      bpsicnutaulr(i, 2) = t2(2) 
      
      enddo loopthtau
!-----------------------------------------------------
      
      aPsi0PsicW(2, 1) = - g
      bPsi0PsicW(2, 1) = - g
      aPsi0PsicW(4, 2) =   g / dsqrt(2.d0)		     
      bPsi0PsicW(3, 2) = - g / dsqrt(2.d0)	     

      loop4by2i: DO i = 1, 4
      loop4by2j: DO j = 1, 2
      aChi0PsicW(i,j) = 0.d0
      bChi0PsicW(i,j) = 0.d0
      loop4by2k: DO k = 1, 4  
            
      aChi0PsicW(i,j) = aChi0PsicW(i,j) + ON(i,k)*aPsi0PsicW(k,j)
      bChi0PsicW(i,j) = bChi0PsicW(i,j) + ON(i,k)*bPsi0PsicW(k,j)
      
      ENDDO loop4by2k
      ENDDO loop4by2j
      ENDDO loop4by2i
      

      loop1: do i = 1, 2
      loop2: do j = 1, 2
      
      aPsiPsiZ(i, j) = 0.d0
      bPsiPsiZ(i, j) = 0.d0
      
      enddo loop2
      enddo loop1


      aPsiPsiZ(1, 1) =   g * costhW
      aPsiPsiZ(2, 2) =   g * cos(2.d0 * thw) / (2.d0 * costhW)
      
      bPsiPsiZ(1, 1) =  aPsiPsiZ(1, 1)
      bPsiPsiZ(2, 2) =  aPsiPsiZ(2, 2)

      call dag2d(OCR,OCRdag)
      call dag2d(OCL,OCLdag)

      call matmult2d(aPsiPsiZ,OCRdag,aPsiChiZ)
      call matmult2d(bPsiPsiZ,OCLdag,bPsiChiZ)


      e = g * sinthw      
      
      call scmul2d(OCR,e,aPsiChiGam)
      call scmul2d(OCLdag,e,bPsiChiGam)     


      aPsiPsiHc1(1, 2) =  gp / dsqrt(2.d0)
      bPsiPsiHc2(1, 2) =  aPsiPsiHc1(1, 2)
      aPsiPsiHc1(2, 2) =  g / dsqrt(2.d0)
      bPsiPsiHc2(2, 2) =  g / dsqrt(2.d0)
      aPsiPsiHc1(3, 1) = -g
      bPsiPsiHc2(4, 1) =  g
      
      loop5by2i: DO i = 1, 4
      loop5by2j: DO j = 1, 2
      loop5by2k: DO k = 1, 4 
      
      aChiPsiHc1(i,j) = aChiPsiHc1(i,j) + ON(i,k) * bPsiPsiHc1(k,j)
      bChiPsiHc1(i,j) = bChiPsiHc1(i,j) + ON(i,k) * aPsiPsiHc1(k,j)
      aChiPsiHc2(i,j) = aChiPsiHc2(i,j) + ON(i,k) * bPsiPsiHc2(k,j)
      bChiPsiHc2(i,j) = bChiPsiHc2(i,j) + ON(i,k) * aPsiPsiHc2(k,j)
      
      ENDDO loop5by2k
      ENDDO loop5by2j
      ENDDO loop5by2i


      loop6: do i = 1, 4
      loop7: do j = 1, 2

      aChiPsiHGp(i, j) =   cosbeta * aChiPsiHc1(i, j) + 
     $     sinbeta * aChiPsiHc2(i, j)

      aChiPsiHHp(i, j) = - sinbeta * aChiPsiHc1(i, j) + 
     $     cosbeta * aChiPsiHc2(i, j)
      
      bChiPsiHGp(i, j) =   cosbeta * bChiPsiHc1(i, j) + 
     $     sinbeta * bChiPsiHc2(i, j)
      
      bChiPsiHHp(i, j) = - sinbeta * bChiPsiHc1(i, j) + 
     $     cosbeta * bChiPsiHc2(i, j)
      
      enddo loop7
      enddo loop6


      aPsiPsis1(1, 2) =   g / dsqrt(2.d0)
      aPsiPsis2(2, 1) =   g / dsqrt(2.d0)
      aPsiPsip1(1, 2) =   g / dsqrt(2.d0)
      aPsiPsip2(2, 1) = - g / dsqrt(2.d0)

!---------------------------------------------------------
      call dag2d(aPsiPsis1,bPsiPsis1)
      call dag2d(aPsiPsis2,bPsiPsis2)

      call dag2d(aPsiPsip1,bPsiPsip1)
      call dag2d(aPsiPsip2,bPsiPsip2)

      call scmul2d(bPsiPsip1, - 1.d0,bPsiPsip1)
      call scmul2d(bPsiPsip2, - 1.d0,bPsiPsip2)
      

      call matmult2d(aPsiPsis1,OCLdag,aPsiChis1)
      call matmult2d(aPsiPsis2,OCLdag,aPsiChis2)
      call matmult2d(aPsiPsip1,OCLdag,aPsiChip1)
      call matmult2d(aPsiPsip2,OCLdag,aPsiChip2)

      call matmult2d(bPsiPsis1,OCRdag,bPsiChis1)
      call matmult2d(bPsiPsis2,OCRdag,bPsiChis2)
      call matmult2d(bPsiPsip1,OCRdag,bPsiChip1)
      call matmult2d(bPsiPsip2,OCRdag,bPsiChip2)
!-----------------------------------------------------
      
      loop11: do i = 1, 2
      loop12: do j = 1, 2

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

!------------------------------------------------------------------------
C   PV functions required for one loop corrections to chargino
C------------------------------------------------------------------------
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

      call b1(p, mchargino(i), MZ, q, b1mchar(i))
      call b0(p, mchargino(i), MZ, q, b0mchar(i))

      call b1(p, mchargino(i), 0.d0, q, b1mchargm(i))
      call b0(p, mchargino(i), 0.d0, q, b0mchargm(i))

      call b1(p, mchargino(i), mGp, q, b1mGpmchar(i))
      call b0(p, mchargino(i), mGp, q, b0mGpmchar(i))

      call b1(p, mchargino(i), mHpm, q, b1mHpmchar(i))
      call b0(p, mchargino(i), mHpm, q, b0mHpmchar(i))

      call b1(p, mchargino(i), mHu0, q, b1mHu0mchar(i))
      call b0(p, mchargino(i), mHu0, q, b0mHu0mchar(i))

      call b1(p, mchargino(i), mh0, q, b1mh0mchar(i))
      call b0(p, mchargino(i), mh0, q, b0mh0mchar(i))

      call b1(p, mchargino(i), mG0, q, b1mG0mchar(i))
      call b0(p, mchargino(i), mG0, q, b0mG0mchar(i))

      call b1(p, mchargino(i), mA0, q, b1mA0mchar(i))
      call b0(p, mchargino(i), mA0, q, b0mA0mchar(i))

      enddo loopb1

!----------------------------------
      loopb14: do i = 1, 4

      call b1(p, mneutm(i), MW, q, b1mneut(i))
      call b0(p, mneutm(i), MW, q, b0mneut(i))

      call b1(p, mneutm(i), mGp, q, b1mGpmneut(i))
      call b0(p, mneutm(i), mGp, q, b0mGpmneut(i))

      call b1(p, mneutm(i), mHpm, q, b1mHpmneut(i))
      call b0(p, mneutm(i), mHpm, q, b0mHpmneut(i))

      
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

C---------------------------------------------------------------
C  corrections  begin
C--------------------------------
C     sfermion contribution
C---------------------------------------------------------------

      loopsig1i: do i = 1, 2
      loopsig1j: do j = 1, 2
            sigmaL(i,j) = 0.d0
            sigmaR(i,j) = 0.d0
            sigmaS(i,j) = 0.d0
      loopsig1k: do k = 1, 2
      

      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 * 
     $     (3.d0 * apsicdq(i, k) * apsicdq(j, k) * real(b1msup(k)) +
     $     3.d0 * apsicuq(i, k) * apsicuq(j, k) * real(b1msdown(k)) +
     $     apsicnuel(i, k) * apsicnuel(j, k) * real(b1msel(k)) +
     $     apsicenul(i, k) * apsicenul(j, k) * real(b1msnue))

      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 * 
     $     (3.d0 * apsicdq(i, k) * apsicdq(j, k) * real(b1mscharm(k)) +
     $     3.d0 * apsicuq(i, k) * apsicuq(j, k) * real(b1msstrange(k)) +
     $     apsicnuel(i, k) * apsicnuel(j, k) * real(b1msmu(k)) +
     $     apsicenul(i, k) * apsicenul(j, k) * real(b1msnumu))
      
      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 * 
     $     (3.d0 * apsicbtqr(i, k) * apsicbtqr(j, k) * real(b1mstop(k))+
     $     3.d0 * apsictbqr(i, k) * apsictbqr(j, k) * real(b1msbot(k)) +
     $     apsicnutaulr(i, k) * apsicnutaulr(j, k) * real(b1mstau(k)) +
     $     apsictaunul(i, k) * apsictaunul(j, k) * real(b1msnutau))

      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 * 
     $     (3.d0 * bpsicdq(i, k) * bpsicdq(j, k) * real(b1msup(k)) +
     $     3.d0 * bPsicuq(i, k) * bPsicuq(j, k) * real(b1msdown(k)) +
     $     bpsicnuel(i, k) * bpsicnuel(j, k) * real(b1msel(k)) +
     $     bpsicenul(i, k) * bpsicenul(j, k) * real(b1msnue))

      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 * 
     $     (3.d0 * bpsicdq(i, k) * bpsicdq(j, k) * real(b1mscharm(k)) +
     $     3.0 * bPsicuq(i, k) * bPsicuq(j, k) * real(b1msstrange(k)) +
     $     bpsicnuel(i, k) * bpsicnuel(j, k) * real(b1msmu(k)) +
     $     bpsicenul(i, k) * bpsicenul(j, k) * real(b1msnumu))
      
      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 * 
     $     (3.d0 * bpsicbtqr(i, k) * bpsicbtqr(j, k) * real(b1mstop(k))+
     $     3.d0 * bpsictbqr(i, k) * bpsictbqr(j, k) * real(b1msbot(k)) +
     $     bpsicnutaulr(i, k) * bpsicnutaulr(j, k) * real(b1mstau(k)) +
     $     bpsictaunul(i, k) * bpsictaunul(j, k) * real(b1msnutau))
      
      sigmaS(i, j) = sigmaS(i, j) + 
     $     (3.d0 * bpsicbtqr(i, k) * apsicbtqr(j, k) * mt * 
     $     real(b0mstop(k)) +
     $     3.d0 * bpsictbqr(i, k) * apsictbqr(j, k) * mb * 
     $     real(b0msbot(k)) +
     $     bpsicnutaulr(i, k) * apsicnutaulr(j, k) * mtau *
     $     real(b0mstau(k)))


      enddo loopsig1k
      enddo loopsig1j
      enddo loopsig1i

C-----------------------------------------------------------------------
C     Higgs-boson Contribution
C------------------------------------------------------------------------

      loopsig2i: do i = 1, 2
      loopsig2j: do j = 1, 2
      loopsig2k: do k = 1, 4
      
      sigmaL(i, j) = sigmaL(i, j) +
     $     (aChi0PsicW(k, i) * aChi0PsicW(k, j) * real(b1mneut(k)))
      
      sigmaR(i, j) = sigmaR(i, j) +
     $     (bChi0PsicW(k, i) * bChi0PsicW(k, j) * real(b1mneut(k)))
      
      sigmaS(i, j) = sigmaS(i, j) - 4.d0 * mneut(k) * 
     $     (bChi0PsicW(k, i) * aChi0PsicW(k, j) * real(b0mneut(k)))
      
      
      enddo loopsig2k
      enddo loopsig2j
      enddo loopsig2i 

C-----------------------------------------------------------------------
C     Z
      
      loopsig3i: do i = 1, 2
      loopsig3j: do j = 1, 2
      loopsig3k: do k = 1, 2
      
      
      sigmaL(i, j) = sigmaL(i, j) +
     $     (aPsiChiZ(i, k) * aPsiChiZ(j, k) * real(b1mchar(k)))
      
      sigmaR(i, j) = sigmaR(i, j) +
     $     (bPsiChiZ(i, k) * bPsiChiZ(j, k) * real(b1mchar(k)))
      
      sigmaS(i, j) = sigmaS(i, j) - 4.d0 * mchargino(k) * 
     $     (bPsiChiZ(i, k) * aPsiChiZ(j, k) * real(b0mchar(k)))
      
      
      sigmaL(i, j) = sigmaL(i, j) +
     $     (aPsiChiGam(i, k) * aPsiChiGam(j, k) * real(b1mchargm(k)))
      
      sigmaR(i, j) = sigmaR(i, j) +
     $     (bPsiChiGam(i, k) * bPsiChiGam(j, k) * real(b1mchargm(k)))
      
      sigmaS(i, j) = sigmaS(i, j) - 4.d0 * mchargino(k) *
     $     (bPsiChiGam(i, k) * aPsiChiGam(j, k) * real(b0mchargm(k)))
      
      
      enddo loopsig3k
      enddo loopsig3j
      enddo loopsig3i 
      
C-------------------------------------------------------------------------      
    
      loopsig4i: do i = 1, 2
      loopsig4j: do j = 1, 2
      loopsig4k: do k = 1, 4
!     G+ 
           
      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 *
     $     (aChiPsiHGp(k, i) * aChiPsiHGp(k, j) * real(b1mGpmneut(k)))
      
      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 *
     $     (bChiPsiHGp(k, i) * bChiPsiHGp(k, j) * real(b1mGpmneut(k)))
      
      sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
     $     (bChiPsiHGp(k, i) * aChiPsiHGp(k, j) * real(b0mGpmneut(k)))

C--------------------------------------------------------------------------      
!     H+	
      
      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 *
     $     (aChiPsiHHp(k, i) * aChiPsiHHp(k, j) * real(b1mHpmneut(k)))
      
      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 *
     $     (bChiPsiHHp(k, i) * bChiPsiHHp(k, j) * real(b1mHpmneut(k)))

      sigmaS(i, j) = sigmaS(i, j) + mneut(k) * 
     $     (bChiPsiHHp(k, i) * aChiPsiHHp(k, j) * real(b0mHpmneut(k)))


      enddo loopsig4k
      enddo loopsig4j
      enddo loopsig4i

C---------------------------------------------------------------------------

      loopsig5i: do i = 1, 2
      loopsig5j: do j = 1, 2
      sigmal(i,j) = 0.d0
      sigmar(i,j) = 0.d0
      sigmas(i,j) = 0.d0

      loopsig5k: do k = 1, 2
      
!     H

      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 *
     $     (aPsiChiHu(i, k) * aPsiChiHu(j, k) * real(b1mHu0mchar(k)))
      
      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 *
     $     (bPsiChiH(i, k) * bPsiChiH(j, k) * real(b1mHu0mchar(k)))
      
      sigmaS(i, j) = sigmaS(i, j) + mchargino(k) * 
     $     (bPsiChiH(i, k) * aPsiChiH(j, k) * real(b0mHu0mchar(k)))

C------------------------------------------------            
!     h
      
      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 *
     $     (aPsiChih(i, k) * aPsiChih(j, k) * real(b1mh0mchar(k)))

      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 *
     $     (bPsiChih(i, k) * bPsiChih(j, k) * real(b1mh0mchar(k)))
      
      sigmaS(i, j) = sigmaS(i, j) + mchargino(k) * 
     $     (bPsiChih(i, k) * aPsiChih(j, k) * real(b0mh0mchar(k)))

C------------------------------------------------      
!     G0	  

      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 *
     $     (aPsiChiG(i, k) * aPsiChiG(j, k) * real(b1mG0mchar(k)))

      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 *
     $     (bPsiChiG(i, k) * bPsiChiG(j, k) * real(b1mG0mchar(k)))

      sigmaS(i, j) = sigmaS(i, j) + mchargino(k) * 
     $     (bPsiChiG(i, k) * aPsiChiG(j, k) * real(b0mG0mchar(k)))

C-------------------------------------------------	
!     A0


      sigmaL(i, j) = sigmaL(i, j) + 0.5d0 *
     $     (aPsiChiA(i, k) * aPsiChiA(j, k) * real(b1mA0mchar(k)))

      sigmaR(i, j) = sigmaR(i, j) + 0.5d0 *
     $     (bPsiChiA(i, k) * bPsiChiA(j, k) * real(b1mA0mchar(k)))

      sigmaS(i, j) = sigmaS(i, j) + mchargino(k) * 
     $     (bPsiChiA(i, k) * aPsiChiA(j, k) * real(b0mA0mchar(k)))

      
      enddo loopsig5k
      enddo loopsig5j
      enddo loopsig5i
      
C------------------------------------------------------------------------

      loopcharmi: do i = 1, 2
      loopcharmj: do j = 1, 2
      charmass(i,j) = 0.d0
      charmass(i,j) = MChar(i,j)
      
      enddo loopcharmj
      enddo loopcharmi

      call matmult2d(sigmaR,charmass,delta1)
      call matmult2d(charmass,sigmaL,delta2)
     
!--------------------------------------------           
      
      loopaddi: do i = 1, 2
      loopaddj: do j = 1, 2
      
      deltatot(i,j) = delta1(i,j) + delta2(i,j) + sigmaS(i,j)

      deltaM(i,j) = ( - deltatot(i,j)/(16.d0 * pi * pi))
      
      charmasstot(i,j) = charmass(i,j) + deltaM(i,j)

c      print*," charmasstot(i,j) = ", charmasstot(i,j) 
      
      charmasstotT(j,i) =  charmasstot(i,j)

      enddo loopaddj
      enddo loopaddi

      
      
C----------------------------------------------------------
C     Diagonalizing subroutine for the corrected mass matrix
C----------------------------------------------------------

 
!      call dgesvd('A','A',2,2,charmasstot,2,Cegm,OCLTm,2,OCRm,2,
!     $     work1,lwork1,info)

      
      call SVD(2, 2, charmasstot,2, Cegm, OCLm,2, OCRm,2, 1)
      
c$$$  do a = 1, 2
c$$$  do m = 1,2 
c$$$  OCRm(a,m) = OCRTm(a,m)
c$$$  OCLTm(m,a) = OCLTm(a,m)
c$$$  enddo
c$$$  enddo 


c$$$      call matmult2d(charmasstot, charmasstotT, mmdag)
c$$$      call matmult2d(charmasstotT, charmasstot, mdagm)
c$$$
c$$$      thetaocr = datan((2.d0 * mdagm(1,2)/
c$$$     $     (mdagm(1,1) - mdagm(2,2))))/2.d0
c$$$      thetaocl = datan((2.d0 * mmdag(1,2)/
c$$$     $     (mmdag(1,1) - mmdag(2,2))))/2.d0
c$$$
c$$$
c$$$      OCLm(1,1) = dcos(thetaocL)
c$$$      OCLm(1,2) = dsin(thetaocL)
c$$$      OCLm(2,1) = -dsin(thetaocL)
c$$$      OCLm(2,2) = dcos(thetaocL)
c$$$
c$$$      OCRm(1,1) = dcos(thetaocR)
c$$$      OCRm(1,2) = dsin(thetaocR)
c$$$      OCRm(2,1) = -dsin(thetaocR)
c$$$      OCRm(2,2) = dcos(thetaocR)
c$$$
c$$$      do a = 1, 2
c$$$         do m = 1, 2 
c$$$            
c$$$            OCRTm(m,a) = OCRm(a,m)
c$$$            OCLTm(m,a) = OCLm(a,m)
c$$$
c$$$         enddo
c$$$      enddo 
c$$$
c$$$
c$$$      call mat3prod2d(OCLm, charmasstot, OCRTm, Mchard)
c$$$
c$$$
c$$$      if(dabs(Mchard(1,1)).lt.1.d-6) then
c$$$         
c$$$         Cegm(1) = Mchard(1,2)
c$$$         Cegm(2) = Mchard(2,1)
c$$$         
c$$$      else
c$$$         
c$$$         Cegm(1) = Mchard(1,1)
c$$$         Cegm(2) = Mchard(2,2)
c$$$         
c$$$      endif
c$$$
c$$$      



      return

      end subroutine chargino

C==========================================================================================
