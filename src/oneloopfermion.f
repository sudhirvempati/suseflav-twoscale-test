****f* SuSeFLAV/oneloopfermion.f 
*  NAME
*    Subroutine topcor, bottomcor and taucor
*  SYNOPSIS
*    One loop correction to fermions. 
*  FUNCTION
*     Computes self energy for third generation fermions at a given energy scale and external momenta
*  INPUTS
*     p                              - External momentum  
*     q                              - Energy scale
*     g,gp,g3                        - Gauge couplings(g = em, gp = weak, g3 = strong)                                     
*     mt,mb,mtau                     - pole masses of top, bottom and tau
*     yuRG,ydRG,yeRG                 - Yukawas
*     mh0sq,mhu0sq,mhpmsq,mA0sq      - physicsal higgs mass squared 
*     M3t                            - Gluino mass at msusy
*     tanbeta                        - the ratio of the vevs of the two–Higgs doublet fields.
*     SUegg   =  6 eigenvalues (ascending order) of UP-Squark mass matrix.
*     SDegg   =  6 eigenvalues (ascending order) of Down-Squark mass matrix.
*     Slegg   =  6 eigenvalues (ascending order) of slepton mass matrix.
*     SNegg   =  3 eigenvalues (ascending order) of sneutrino mass matrix.
*     ON      =  (4 X 4) orthogonal matrix such that 
*                ON.MNeut.Transpose[ON] = Diag[MNeut] 
*     Neg     =  4 singular values (descending order) of the Neutralino mass matrix. 
*
*     OCR, OCL =  (2 X 2) orthogonal matrices such that 
*                 MChar = Transpose[OCR].Diag[MChar].OCL
*     Ceg     =   2 singular values of the Neutralino Mass Matrix
*
*  RESULT
*   correctiom,bcor     =  correction to top, bottom and tau
*   mtaucor
*
*  EXAMPLE
*
*      SUBROUTINE topcor(p,q,g,gp,g3,M3t,mt,mb,tanbeta,
*     $     yuRG,ydRG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
*     $     mh0sq,mhu0sq,mhpmsq,mA0sq,ON,OCL,
*     $     OCR,sinsqtheff,correction)
*
*      SUBROUTINE bottomcor(p,q,g,gp,g3,M3t,mt,mb,tanbeta,
*     $     yuRG,ydRG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
*     $     mh0sq,mhu0sq,mhpmsq,mA0sq,ON,OCL,
*     $     OCR,M2tz,sinsqtheff,mbcor,mbdrbar)
*
*      SUBROUTINE taucor(p,q,g,gp,M3t,mtau,tanbeta,yeRG,
*     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,
*     $     mhu0sq,mhpmsq,mA0sq,ON,OCL,
*     $     OCR,sinsqtheff,mtaucor)
*
*  NOTES
*    1. q, the energy scale at which the corrections are added = MZ
*    2. Conventions and notations followed are that of BPMZ.
*  BUGS
*    ---
*  SEE ALSO
*   ---
*
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!================================================================================================
C
C      Corrections to fermions
C     
C---------------------------------------------------------------------------      
      SUBROUTINE topcor(p,q,g,gp,g3,M3t,mt,mb,tanbeta,
     $     yuRG,ydRG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,ON,OCL,
     $     OCR,sinsqtheff,correction)

      IMPLICIT NONE 
      
      INTEGER i,j
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),correction
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,g3,M3t
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      DOUBLE PRECISION yuRG(3,3),ydRG(3,3)
      
      double precision mT, mB, sinsqtheff
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2

      DOUBLE PRECISION ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta,alpha

      DOUBLE PRECISION b1mneut1mt1,b1mneut1mt2,b1mneut2mt1,b1mneut2mt2
      DOUBLE PRECISION b1mneut3mt1,b1mneut3mt2,b1mneut4mt1,b1mneut4mt2
      DOUBLE PRECISION b0mneut1mt1,b0mneut1mt2,b0mneut2mt1,b0mneut2mt2
      DOUBLE PRECISION b0mneut3mt1,b0mneut3mt2,b0mneut4mt1,b0mneut4mt2

      DOUBLE PRECISION sinsqthw,g
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION beta
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetatz,thetabz,thetatauz
      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetacz,thetasz,thetamuz
      DOUBLE PRECISION thetauz,thetadz,thetaez
    
      DOUBLE PRECISION mu,rthetat(2,2)
      DOUBLE PRECISION stopterm, mchargino(2)

      DOUBLE PRECISION b1M3mt1,b1M3mt2,b0M3mt1,b0M3mt2
      DOUBLE PRECISION toptermt,btmtermt,tneutterm,charginoterm
      
      DOUBLE PRECISION b1mbMw,b1mbmHpm,b0mbMw,b0mbmHpm
      DOUBLE PRECISION b1mtmHu,b1mtmh,b1mtMz,b1mtmA
      DOUBLE PRECISION b0mtmHu,b0mtmh,b0mtMz,b0mtmA
      DOUBLE PRECISION b1mch1mb1,b1mch1mb2,b1mch2mb1,b1mch2mb2
      DOUBLE PRECISION b0mch1mb1,b0mch1mb2,b0mch2mb1,b0mch2mb2

      DOUBLE PRECISION fnotf(4,2),gnotf(4,2)
      DOUBLE PRECISION nmneut(4)

      DOUBLE PRECISION hxnot(4,2),bxnot(4,2),hxpos(2,2),bxpos(2,2)
      DOUBLE PRECISION tt(2), t1(2), t2(2), rthetab(2,2)
      DATA tt/ 2 * 0.d0/, t1/ 2 * 0.d0/, t2/ 2 *0.d0/,
     $     rthetab/ 4 *0.d0/
!-----------------------------------------------------------------
      DOUBLE PRECISION  bPsicsbl(2), bPsicsbr(2), aPsicsbl(2)
      data bPsicsbl/ 2 * 0.d0/, bPsicsbr/ 2 * 0.d0/,
     $     aPsicsbl/ 2 * 0.d0/

      DOUBLE PRECISION  aPsicsbr(2) !, aPsicCSbotl(2)
      DATA aPsicsbr/ 2 * 0.d0/  !, aPsicCSbotl/ 2 * 0.d0/
      
      DOUBLE PRECISION aChicsbr(2), aChicsbl(2), bChicsbl(2) !<---------SHOULD BE COMPLEX
      DATA aChicsbr/ 2 * 0.d0/,aChicsbl/ 2 * 0.d0/, 
     $     bChicsbl/ 2 * 0.d0/
      
      DOUBLE PRECISION bChicsbr(2), aChTsb(2, 2), bChTsb(2, 2)
      DATA bChicsbr/ 2 * 0.d0/,aChTsb/ 4 * 0.d0/, 
     $     bChTsb/ 4 * 0.d0/
      
      DOUBLE PRECISION  fChTsb(2, 2), gChTsb(2, 2)
      DATA fChTsb/ 4 * 0.d0/, gChTsb/ 4 * 0.d0/
      
      DOUBLE PRECISION aPsi0str(4), aPsi0stl(4), bPsi0str(4),
     $     bPsi0stl(4), achif(4,2), bchif(4,2)
      
      DOUBLE PRECISION aChi0str(4), aChi0stl(4), bChi0str(4),
     $     bChi0stl(4)
      
      DATA aChi0str/ 4 * 0.d0/, aChi0stl/ 4 * 0.d0/,
     $     bChi0str/ 4 * 0.d0/, bChi0stl/ 4 * 0.d0/
!-----------------------------------------------------
      DOUBLE PRECISION mbpole, mtaupole, Mtpole, MWpole, MZpole
      double precision MW,MZ,pi
!----------------

      DOUBLE PRECISION  sinsqthw_mz
      common/sinsq_mz/sinsqthw_mz
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole
      common/sfmixing_mz/thetatz,thetabz,thetatauz,thetacz,thetasz,
     $      thetamuz,thetauz,thetadz,thetaez

      
      EXTERNAL b1,b0,mat3prod4d,dag4d,dag2d,mat3prod2d,theta,rmat2d
      

      include 'stdinputs.h'
!---------------------------------------------------------------------------

      pi = 4.d0*datan(1.d0)
      MW = MWpole
      MZ = MZpole

      mu = mUQ

      sinsqthw = sinsqtheff
      sinthw = dsqrt(sinsqtheff)
      thw = dasin(sinthw)
      costhw = dcos(thw)
      cos2thw = dcos(2.d0*thw)
      
      gnuL = 0.5d0
      guL = 0.5d0 - 2.d0*sinsqthw/ 3.d0 
      gdL = -0.5d0 + sinsqthw/3.d0
      geL = -0.5d0 + sinsqthw
      guR =  2.d0*sinsqthw/3.d0 
      gdR = -sinsqthw/3.d0
      geR = -sinsqthw 
      yuL = 1.d0/3.d0 
      yuR = -4.d0/3.d0 
      ydL = 1.d0/3.d0 
      ydR = 2.d0/3.d0 
      yeL = -1.d0 
      yeR = 2.d0 
      ynuL = -1.d0       

C------------------------------------------------------------

      muR = dsqrt(SUegg(5)) 
      muL = dsqrt(SUegg(6))
      mcR = dsqrt(SUegg(3))
      mcL = dsqrt(SUegg(4))
      mtR = dsqrt(SUegg(1))
      mtL = dsqrt(SUegg(2))
      
      mdR = dsqrt(SDegg(5)) 
      mdL = dsqrt(SDegg(6))
      msR = dsqrt(SDegg(3))
      msL = dsqrt(SDegg(4))
      mbR = dsqrt(SDegg(1))
      mbL = dsqrt(SDegg(2))

      meR   = dsqrt(SLegg(5))
      meL   = dsqrt(SLegg(6))
      mmuR  = dsqrt(SLegg(3))
      mmuL  = dsqrt(SLegg(4))
      mtauR = dsqrt(SLegg(1))
      mtauL = dsqrt(SLegg(2))

      
      snu(1) = dsqrt(SNegg(1))
      snu(2) = dsqrt(SNegg(2))
      snu(3) = dsqrt(SNegg(3))

      mneut(1) = Neg(1) 
      mneut(2) = Neg(2)
      mneut(3) = Neg(3)
      mneut(4) = Neg(4)


      nmneut(1) = dabs(Neg(1)) 
      nmneut(2) = dabs(Neg(2))
      nmneut(3) = dabs(Neg(3))
      nmneut(4) = dabs(Neg(4))

      mchargino(1) = Ceg(1)
      mchargino(2) = Ceg(2)
 
!-----------------------------------------------------------------------


      costhetat = dcos(thetatz)
      sinthetat = dsin(thetatz)
      costhetab = dcos(thetabz)
      sinthetab = dsin(thetabz)
      costhetatau = dcos(thetatauz)
      sinthetatau = dsin(thetatauz)


      mt1 = mtL
      mt2 = mtR

      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR

c$$$      mt1 = mtR
c$$$      mt2 = mtL
c$$$
c$$$      mb1 = mbR
c$$$      mb2 = mbL
c$$$
c$$$      mtau1 = mtauR
c$$$      mtau2 = mtauL

C----------------------------------------------------------------

      mh0 = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0= dsqrt((mA0sq))


      beta = datan(tanbeta)
      tan2beta = dtan(2.d0 * beta)
      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta = dcos(beta)
      cos2beta = dcos(2.d0 * beta)
      sinbeta = dsin(beta)
      sin2beta = dsin(2.d0 * beta)
      correction = 0.d0
C--------------------------------------------------------------------


      call b1(p,M3t,mt1,q,b1M3mt1)
      call b1(p,M3t,mt2,q,b1M3mt2)

      call b0(p,M3t,mt1,q,b0M3mt1)
      call b0(p,M3t,mt2,q,b0M3mt2)


      call b1(p,mt,mHu0,q,b1mtmHu)
      call b1(p,mt,mh0,q,b1mtmh)

      call b1(p,mt,mA0,q,b1mtmA)
      call b1(p,mt,MZ,q,b1mtMz)

      call b1(p,mb,mHpm,q,b1mbmHpm)
      call b1(p,mb,MW,q,b1mbMw)


      call b0(p,mt,mHu0,q,b0mtmHu)
      call b0(p,mt,mh0,q,b0mtmh)

      call b0(p,mt,mA0,q,b0mtmA)
      call b0(p,mt,MZ,q,b0mtMz)

      call b0(p,mb,mHpm,q,b0mbmHpm)
      call b0(p,mb,MW,q,b0mbMw)

      
C---------------------------------------------------------------------
      stopterm = 0.d0
      toptermt = 0.d0
      btmtermt = 0.d0

      stopterm = g3**2.d0*(4.d0/3.d0)*((b1M3mt1 + b1M3mt2) - 
     $     (dsin(2.d0*thetatz)*M3t/mtpole*(b0M3mt1 - b0M3mt2))) 


      toptermt = 0.5d0*yuRG(3,3)**2.d0*(((salpha**2.d0)*(b1mtmHu + 
     $     b0mtmHu)) +    
     $     ((calpha**2.d0)*(b1mtmh + b0mtmh)) + 
     $     ((cosbeta**2.d0)*(b1mtmA - b0mtmA)) +
     $     ((sinbeta**2.d0)*(b1mtMz - b0mtMz))) +
     $     ((g/costhw)**2.d0*((guL*guL + guR*guR)*b1mtMz + !<---------check, gtl,gtr: def!!---------------
     $     4.d0*guL*guR*b0mtMz))

      btmtermt = 0.5d0*((ydRG(3,3)*sinbeta)**2.d0 + 
     $     (yuRG(3,3)*cosbeta)**2.d0)*b1mbmHpm +
     $     0.5d0*(g*g + (ydRG(3,3)*cosbeta)**2.d0 + 
     $     (yuRG(3,3)*sinbeta)**2.d0)*b1mbMw +
     $     (ydRG(3,3)*cosbeta)**2.d0*(b0mbmHpm - b0mbMw) -
     $     (g*sinthw*2.d0/3.d0)**2.d0*(5.d0 + 6.d0*dlog(q/mt))

C----------------------------------------------------------------------------
      
C     Neutralinos

      loopsii: DO i = 1,4
      loopsij: DO j = 1,2

      fnotf(i,j) = 0.d0
      gnotf(i,j) = 0.d0
      
      achif(i,j) = 0.d0
      bchif(i,j) = 0.d0

      ENDDO loopsij
      ENDDO loopsii


      aPsi0str(1)  = -4.d0*gp/(3.d0*dsqrt(2.d0))
      aPsi0str(2)  = 0.d0
      aPsi0str(3)  = 0.d0
      aPsi0str(4)  = 0.d0

      bPsi0stl(1)  =  gp/(3.d0*dsqrt(2.d0))
      bPsi0stl(2)  =  g/dsqrt(2.d0)
      bPsi0stl(3)  = 0.d0
      bPsi0stl(4)  = 0.d0

      aPsi0stl(1) = 0.d0
      aPsi0stl(2) = 0.d0
      aPsi0stl(3) = 0.d0
      aPsi0stl(4) =  yuRG(3,3)

      bPsi0str(1) = 0.d0
      bPsi0str(2) = 0.d0
      bPsi0str(3) = 0.d0     
      bPsi0str(4) = yuRG(3,3)


C-----------------------------------------------------------------

      loop4by2i: DO i =1,4

      aChi0stl(i) = 0.d0
      bChi0stl(i) = 0.d0
      aChi0str(i) = 0.d0
      bChi0str(i) = 0.d0

      loop4by2j: DO j = 1,4
      


      aChi0stl(i) = aChi0stl(i) + ON(i,j) * aPsi0stl(j)
      bChi0stl(i) = bChi0stl(i) + ON(i,j) * bPsi0stl(j)
      aChi0str(i) = aChi0str(i) + ON(i,j) * aPsi0str(j)
      bChi0str(i) = bChi0str(i) + ON(i,j) * bPsi0str(j)


      ENDDO loop4by2j
      ENDDO loop4by2i
      

      tneutterm = 0.d0



      call b1(p,nmneut(1),mt1,q,b1mneut1mt1)
      call b1(p,nmneut(1),mt2,q,b1mneut1mt2)
      call b1(p,nmneut(2),mt1,q,b1mneut2mt1)
      call b1(p,nmneut(2),mt2,q,b1mneut2mt2)
      call b1(p,nmneut(3),mt1,q,b1mneut3mt1)
      call b1(p,nmneut(3),mt2,q,b1mneut3mt2)
      call b1(p,nmneut(4),mt1,q,b1mneut4mt1)
      call b1(p,nmneut(4),mt2,q,b1mneut4mt2)

      call b0(p,nmneut(1),mt1,q,b0mneut1mt1)
      call b0(p,nmneut(1),mt2,q,b0mneut1mt2)
      call b0(p,nmneut(2),mt1,q,b0mneut2mt1)
      call b0(p,nmneut(2),mt2,q,b0mneut2mt2)
      call b0(p,nmneut(3),mt1,q,b0mneut3mt1)
      call b0(p,nmneut(3),mt2,q,b0mneut3mt2)
      call b0(p,nmneut(4),mt1,q,b0mneut4mt1)
      call b0(p,nmneut(4),mt2,q,b0mneut4mt2)


      hxnot(1,1) = b1mneut1mt1
      hxnot(1,2) = b1mneut1mt2
      hxnot(2,1) = b1mneut2mt1
      hxnot(2,2) = b1mneut2mt2
      hxnot(3,1) = b1mneut3mt1
      hxnot(3,2) = b1mneut3mt2
      hxnot(4,1) = b1mneut4mt1
      hxnot(4,2) = b1mneut4mt2
      
      bxnot(1,1) = mneut(1)*b0mneut1mt1/mtpole
      bxnot(1,2) = mneut(1)*b0mneut1mt2/mtpole
      bxnot(2,1) = mneut(2)*b0mneut2mt1/mtpole
      bxnot(2,2) = mneut(2)*b0mneut2mt2/mtpole
      bxnot(3,1) = mneut(3)*b0mneut3mt1/mtpole
      bxnot(3,2) = mneut(3)*b0mneut3mt2/mtpole
      bxnot(4,1) = mneut(4)*b0mneut4mt1/mtpole
      bxnot(4,2) = mneut(4)*b0mneut4mt2/mtpole

      call rmat2d(thetatz,rthetat)


      loopnlinoi: DO i = 1,4

      t1(1) = rthetat(1,1) * aChi0stl(i) + 
     $     rthetat(1,2) * aChi0str(i)  

      t1(2) = rthetat(2,1) * aChi0stl(i) + 
     $     rthetat(2,2) * aChi0str(i)  


      t2(1) = rthetat(1,1) * bChi0stl(i) +
     $     rthetat(1,2) * bChi0str(i)

      t2(2) = rthetat(2,1) * bChi0stl(i) +
     $     rthetat(2,2) * bChi0str(i)

      loopnlinoj: DO j = 1,2

      achif(i,j) = t1(j)
      bchif(i,j) = t2(j)

      fnotf(i,j)  = achif(i,j)**2.d0 + bchif(i,j)**2.d0
      gnotf(i,j)  = 2.d0 * (achif(i,j) * bchif(i,j))

      tneutterm = tneutterm + 0.5d0*(fnotf(i,j)*
     $     hxnot(i,j) + gnotf(i,j)*bxnot(i,j))
      
      ENDDO loopnlinoj

      ENDDO loopnlinoi 

C--------------------------------
C     Chargino Terms

      charginoterm = 0.d0
      
      call b1(p,dabs(mchargino(1)),mb1,q,b1mch1mb1)
      call b1(p,dabs(mchargino(1)),mb2,q,b1mch1mb2)
      call b1(p,dabs(mchargino(2)),mb1,q,b1mch2mb1)
      call b1(p,dabs(mchargino(2)),mb2,q,b1mch2mb2)

      call b0(p,dabs(mchargino(1)),mb1,q,b0mch1mb1)
      call b0(p,dabs(mchargino(1)),mb2,q,b0mch1mb2)
      call b0(p,dabs(mchargino(2)),mb1,q,b0mch2mb1)
      call b0(p,dabs(mchargino(2)),mb2,q,b0mch2mb2)



      hxpos(1,1) = b1mch1mb1
      hxpos(1,2) = b1mch1mb2
      hxpos(2,1) = b1mch2mb1
      hxpos(2,2) = b1mch2mb2


      bxpos(1,1) = mchargino(1)*b0mch1mb1/mtpole
      bxpos(1,2) = mchargino(1)*b0mch1mb2/mtpole
      bxpos(2,1) = mchargino(2)*b0mch2mb1/mtpole
      bxpos(2,2) = mchargino(2)*b0mch2mb2/mtpole
!----------------------------------------------
C     couplings

       bPsicsbl(1) = g
       bPsicsbr(2) = -ydRG(3,3)
       aPsicsbl(2) = -yuRG(3,3)

!-----------------------------------

       aChicsbl(1) = OCR(1,1) * aPsicsbl(1) + 
     $      OCR(1,2) * aPsicsbl(2)

       aChicsbl(2) = OCR(2,1) * aPsicsbl(1) + 
     $      OCR(2,2) * aPsicsbl(2)

       bChicsbl(1) = OCL(1,1) * bPsicsbl(1) +
     $      OCL(1,2) * bPsicsbl(2)

       bChicsbl(2) = OCL(2,1) * bPsicsbl(1) +
     $      OCL(2,2) * bPsicsbl(2)

       aChicsbr(1) = OCR(1,1) * aPsicsbr(1) +
     $      OCR(1,2) * aPsicsbr(2) 

       aChicsbr(2) = OCR(2,1) * aPsicsbr(1) +
     $      OCR(2,2) * aPsicsbr(2) 


       bChicsbr(1) = OCL(1,1) * bPsicsbr(1) +
     $      OCL(1,2) * bPsicsbr(2)

       bChicsbr(2) = OCL(2,1) * bPsicsbr(1) +
     $      OCL(2,2) * bPsicsbr(2)

!---------------------------------------------------       
       call rmat2d(thetabz,rthetab)


       loopchmti: DO i = 1, 2
       
       tt(1) = aChicsbl(i)
       tt(2) = aChicsbr(i)
       
       t1(1) = rthetab(1,1) *  tt(1) + rthetab(1,2) *  tt(2) 
       t1(2) = rthetab(2,1) *  tt(1) + rthetab(2,2) *  tt(2) 

       tt(1) = bChicsbl(i) 
       tt(2) = bChicsbr(i)
       
       t2(1) = rthetab(1,1) *  tt(1) + rthetab(1,2) *  tt(2) 
       t2(2) = rthetab(2,1) *  tt(1) + rthetab(2,2) *  tt(2) 

       
       loopchmtj: DO j = 1, 2
       
       aChTsb(i, j) = t1(j)
       bChTsb(i, j) = t2(j)


       fChTsb(i, j) = (aChTsb(i, j))**2.d0 + 
     $      (bChTsb(i, j))**2.d0

       gChTsb(i, j) = 2.d0 * (aChTsb(i, j) * 
     $      bChTsb(i, j)) 
       
       charginoterm = charginoterm +
     $      ((fChTsb(i, j) * hxpos(i,j)) +
     $      (gChTsb(i, j) * bxpos(i,j))) * 0.5d0
       
      ENDDO loopchmtj
      ENDDO loopchmti



!----------------------------------------------------------------

      correction = (stopterm + btmtermt + toptermt + tneutterm + 
     $     charginoterm )*
     $     mtpole/(16.d0*pi*pi)

c$$$      print*,"stopterm, btmtermt+toptermt, tneutterm, charginoterm, 
c$$$     $mt, mtpole, alphas = ", stopterm, btmtermt+toptermt, tneutterm, 
c$$$     $     charginoterm, mt, mtpole, g3*g3/(4.d0*pi)

      RETURN

      END SUBROUTINE topcor
C===========================================================================


      SUBROUTINE bottomcor(p,q,g,gp,g3,M3t,mt,mb,tanbeta,
     $     yuRG,ydRG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,ON,OCL,
     $     OCR,sinsqtheff,mbcor,mbdrbar)


      IMPLICIT NONE 
      
      INTEGER i,j
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),mbcor,correction
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,g3,M3t
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      DOUBLE PRECISION yuRG(3,3),ydRG(3,3)
      
      double precision mT, mB, mbdrbar
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2

      DOUBLE PRECISION sinsqtheff
      DOUBLE PRECISION ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta,alpha

      DOUBLE PRECISION sinsqthw,g
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION beta
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetatz,thetabz,thetatauz

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetacz,thetasz,thetamuz
      DOUBLE PRECISION thetauz,thetadz,thetaez
      DOUBLE PRECISION mu,rthetat(2,2),t1(2),t2(2)
      DOUBLE PRECISION stopterm
      DOUBLE PRECISION mchargino(2),tt(2),rthetab(2,2)

      DOUBLE PRECISION b1M3mb1,b1M3mb2,b0M3mb1,b0M3mb2
      DOUBLE PRECISION toptermt,btmtermt,tneutterm,charginoterm
      
      DOUBLE PRECISION b1mtMw,b1mtmHpm,b0mtMw,b0mtmHpm
      DOUBLE PRECISION b1mbmHu,b1mbmh,b1mbMz,b1mbmA,b0mchmt(2,2)
      DOUBLE PRECISION b0mbmHu,b0mbmh,b0mbMz,b0mbmA,b1mchmt(2,2)
      DOUBLE PRECISION b0mneutmb(4,2),b1mneutmb(4,2)
!----------------------------------------------------------------------
      DOUBLE PRECISION  bPsicstl(2), bPsicstr(2)
      data bPsicstl/ 2 * 0.d0/,bPsicstr/ 2 * 0.d0/

      DOUBLE PRECISION aPsicstl(2), aPsicstr(2) 
      data aPsicstl/ 2 * 0.d0/,aPsicstr/ 2 * 0.d0/

      DOUBLE PRECISION aChicstr(2) !,aPsicCStopl(2)
      data aChicstr/ 2 * 0.d0/ !,aPsicCStopl/ 2 * 0.d0/

      DOUBLE PRECISION aChicstl(2), bChicstl(2),bChicstr(2)
      data aChicstl/ 2 * 0.d0/,bChicstl/ 2 * 0.d0/,
     $     bChicstr/ 2 * 0.d0/

      DOUBLE PRECISION aChBstop(2, 2), bChBstop(2, 2)
      data aChBstop/ 4 * 0.d0/,bChBstop/ 4 * 0.d0/

      DOUBLE PRECISION fChBstop(2, 2), gChBstop(2, 2)
      data fChBstop/ 4 * 0.d0/,gChBstop/ 4 * 0.d0/

      DOUBLE PRECISION aPsi0sbr(4), bPsi0sbr(4)
      data aPsi0sbr/ 4 * 0.d0/,bPsi0sbr/ 4 * 0.d0/

      DOUBLE PRECISION aPsi0sbl(4),bPsi0sbl(4) 
      data aPsi0sbl/ 4 * 0.d0/,bPsi0sbl/ 4 * 0.d0/

      DOUBLE PRECISION aChi0sbl(4), bChi0sbl(4)
      data aChi0sbl/ 4 * 0.d0/, bChi0sbl/ 4 * 0.d0/

      DOUBLE PRECISION aChi0sbr(4), bChi0sbr(4)
      data aChi0sbr/ 4 * 0.d0/, bChi0sbr/ 4 * 0.d0/

      DOUBLE PRECISION aNeutsb(4, 2), bNeutsb(4, 2)
      data aNeutsb/ 8 * 0.d0/,bNeutsb/ 8 * 0.d0/
   
      DOUBLE PRECISION fNeutsb(4, 2), gNeutsb(4, 2)
      data fNeutsb/ 8 * 0.d0/,gNeutsb/ 8 * 0.d0/

      DOUBLE PRECISION neutm(4, 2),nmneut(4)

      DOUBLE PRECISION mbpole, mtaupole, Mtpole,mbmzdrbar,mbMZmsbar
      double precision thetabz1, thetatz1
      double precision MWpole, MZpole,MW,MZ,pi

!----------------

      common/sfmixing_mz/thetatz,thetabz,thetatauz,thetacz,thetasz,
     $     thetamuz,thetauz,thetadz,thetaez
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole
      common/qcd_cor/mbmzdrbar,mbMZmsbar
      

!-----
      EXTERNAL b1,b0,mat3prod4d,dag4d,dag2d,mat3prod2d,theta,rmat2d


      include 'stdinputs.h'
!---------------------------------------------------------------------------


      MW = Mwpole
      MZ = MZpole
      mu = mUQ

      sinsqthw = sinsqtheff 
      sinthw = dsqrt(sinsqtheff)
      thw = dasin(sinthw)
      costhw = dcos(thw)
      cos2thw = dcos(2.d0*thw)
      pi = 4.d0*datan(1.d0)   

      gnuL = 0.5d0
      guL = 0.5d0 - 2.d0*sinsqthw/ 3.d0 
      gdL = -0.5d0 + sinsqthw/3.d0
      geL = -0.5d0 + sinsqthw
      guR =  2.d0*sinsqthw/3.d0 
      gdR = -sinsqthw/3.d0
      geR = -sinsqthw 
      yuL = 1.d0/3.d0 
      yuR = -4.d0/3.d0 
      ydL = 1.d0/3.d0 
      ydR = 2.d0/3.d0 
      yeL = -1.d0 
      yeR = 2.d0 
      ynuL = -1.d0 

!----------------------------------

      thetabz1 = thetabz
      thetatz1 = thetatz

C------------------------------------------------------------

c$$$      muR = dsqrt(SUegg(6)) 
c$$$      muL = dsqrt(SUegg(5))
c$$$      mcR = dsqrt(SUegg(4))
c$$$      mcL = dsqrt(SUegg(3))
c$$$      mtR = dsqrt(SUegg(2))
c$$$      mtL = dsqrt(SUegg(1))
c$$$      
c$$$      mdR = dsqrt(SDegg(6)) 
c$$$      mdL = dsqrt(SDegg(5))
c$$$      msR = dsqrt(SDegg(4))
c$$$      msL = dsqrt(SDegg(3))
c$$$      mbR = dsqrt(SDegg(2))
c$$$      mbL = dsqrt(SDegg(1))
c$$$
c$$$      meR   = dsqrt(SLegg(6))
c$$$      meL   = dsqrt(SLegg(5))
c$$$      mmuR  = dsqrt(SLegg(4))
c$$$      mmuL  = dsqrt(SLegg(3))
c$$$      mtauR = dsqrt(SLegg(2))
c$$$      mtauL = dsqrt(SLegg(1))


      muR = dsqrt(SUegg(5)) 
      muL = dsqrt(SUegg(6))
      mcR = dsqrt(SUegg(3))
      mcL = dsqrt(SUegg(4))
      mtR = dsqrt(SUegg(1))
      mtL = dsqrt(SUegg(2))
      
      mdR = dsqrt(SDegg(5)) 
      mdL = dsqrt(SDegg(6))
      msR = dsqrt(SDegg(3))
      msL = dsqrt(SDegg(4))
      mbR = dsqrt(SDegg(1))
      mbL = dsqrt(SDegg(2))

      meR   = dsqrt(SLegg(5))
      meL   = dsqrt(SLegg(6))
      mmuR  = dsqrt(SLegg(3))
      mmuL  = dsqrt(SLegg(4))
      mtauR = dsqrt(SLegg(1))
      mtauL = dsqrt(SLegg(2))

      
      snu(1) = dsqrt(SNegg(1))
      snu(2) = dsqrt(SNegg(2))
      snu(3) = dsqrt(SNegg(3))

      mneut(1) = Neg(1) 
      mneut(2) = Neg(2)
      mneut(3) = Neg(3)
      mneut(4) = Neg(4)

      nmneut(1) = dabs(Neg(1)) 
      nmneut(2) = dabs(Neg(2))
      nmneut(3) = dabs(Neg(3))
      nmneut(4) = dabs(Neg(4))

      mchargino(1) = Ceg(1)
      mchargino(2) = Ceg(2)

      
!-----------------------------------------------------------------------

      costhetat = dcos(thetatz1)
      sinthetat = dsin(thetatz1)
      costhetab = dcos(thetabz1)
      sinthetab = dsin(thetabz1)
      costhetatau = dcos(thetatauz)
      sinthetatau = dsin(thetatauz)

      mt1 = mtL
      mt2 = mtR

      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR

c$$$      mt1 = mtR
c$$$      mt2 = mtL
c$$$
c$$$      mb1 = mbR
c$$$      mb2 = mbL
c$$$
c$$$      mtau1 = mtauR
c$$$      mtau2 = mtauL

C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))


      beta = datan(tanbeta)
      tan2beta = dtan(2.d0 * beta)
      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/
     $     (mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta = dcos(beta)
      cos2beta = dcos(2.d0 * beta)
      sinbeta = dsin(beta)
      sin2beta = dsin(2.d0 * beta)
      correction = 0.d0

C--------------------------------------------------------------------
      

c$$$      call b1(2.d-5,M3t,mb1,q,b1M3mb1)
c$$$      call b1(2.d-5,M3t,mb2,q,b1M3mb2)
c$$$      call b0(2.d-5,M3t,mb1,q,b0M3mb1)
c$$$      call b0(2.d-5,M3t,mb2,q,b0M3mb2)

      call b1(p,M3t,mb1,q,b1M3mb1)
      call b1(p,M3t,mb2,q,b1M3mb2)
      call b0(p,M3t,mb1,q,b0M3mb1)
      call b0(p,M3t,mb2,q,b0M3mb2)

      call b1(p,mb,mHu0,q,b1mbmHu)
      call b1(p,mb,mh0,q,b1mbmh)

      call b1(p,mb,mA0,q,b1mbmA)
      call b1(p,mb,MZ,q,b1mbMz)

      call b1(p,mt,mHpm,q,b1mtmHpm)
      call b1(p,mt,MW,q,b1mtMw)
      call b0(p,mb,mHu0,q,b0mbmHu)
      call b0(p,mb,mh0,q,b0mbmh)

      call b0(p,mb,mA0,q,b0mbmA)
      call b0(p,mb,MZ,q,b0mbMz)

      call b0(p,mt,mHpm,q,b0mtmHpm)
      call b0(p,mt,MW,q,b0mtMw)

c$$$      call b1(2.d-5,nmneut(1),mb1,q,b1mneutmb(1,1))
c$$$      call b1(2.d-5,nmneut(1),mb2,q,b1mneutmb(1,2))
c$$$      call b1(2.d-5,nmneut(2),mb1,q,b1mneutmb(2,1))
c$$$      call b1(2.d-5,nmneut(2),mb2,q,b1mneutmb(2,2))
c$$$      call b1(2.d-5,nmneut(3),mb1,q,b1mneutmb(3,1))
c$$$      call b1(2.d-5,nmneut(3),mb2,q,b1mneutmb(3,2))
c$$$      call b1(2.d-5,nmneut(4),mb1,q,b1mneutmb(4,1))
c$$$      call b1(2.d-5,nmneut(4),mb2,q,b1mneutmb(4,2))
c$$$
c$$$      call b0(2.d-5,nmneut(1),mb1,q,b0mneutmb(1,1))
c$$$      call b0(2.d-5,nmneut(1),mb2,q,b0mneutmb(1,2))
c$$$      call b0(2.d-5,nmneut(2),mb1,q,b0mneutmb(2,1))
c$$$      call b0(2.d-5,nmneut(2),mb2,q,b0mneutmb(2,2))
c$$$      call b0(2.d-5,nmneut(3),mb1,q,b0mneutmb(3,1))
c$$$      call b0(2.d-5,nmneut(3),mb2,q,b0mneutmb(3,2))
c$$$      call b0(2.d-5,nmneut(4),mb1,q,b0mneutmb(4,1))
c$$$      call b0(2.d-5,nmneut(4),mb2,q,b0mneutmb(4,2))
c$$$      
c$$$
c$$$      call b1(2.d-5,dabs(mchargino(1)),dabs(mt1),q,b1mchmt(1,1))
c$$$      call b1(2.d-5,dabs(mchargino(1)),dabs(mt2),q,b1mchmt(1,2))
c$$$      call b1(2.d-5,dabs(mchargino(2)),dabs(mt1),q,b1mchmt(2,1))
c$$$      call b1(2.d-5,dabs(mchargino(2)),dabs(mt2),q,b1mchmt(2,2))
c$$$
c$$$      call b0(2.d-5,dabs(mchargino(1)),mt1,q,b0mchmt(1,1))
c$$$      call b0(2.d-5,dabs(mchargino(1)),mt2,q,b0mchmt(1,2))
c$$$      call b0(2.d-5,dabs(mchargino(2)),mt1,q,b0mchmt(2,1))
c$$$      call b0(2.d-5,dabs(mchargino(2)),mt2,q,b0mchmt(2,2))

      call b1(p,nmneut(1),mb1,q,b1mneutmb(1,1))
      call b1(p,nmneut(1),mb2,q,b1mneutmb(1,2))
      call b1(p,nmneut(2),mb1,q,b1mneutmb(2,1))
      call b1(p,nmneut(2),mb2,q,b1mneutmb(2,2))
      call b1(p,nmneut(3),mb1,q,b1mneutmb(3,1))
      call b1(p,nmneut(3),mb2,q,b1mneutmb(3,2))
      call b1(p,nmneut(4),mb1,q,b1mneutmb(4,1))
      call b1(p,nmneut(4),mb2,q,b1mneutmb(4,2))

      call b0(p,nmneut(1),mb1,q,b0mneutmb(1,1))
      call b0(p,nmneut(1),mb2,q,b0mneutmb(1,2))
      call b0(p,nmneut(2),mb1,q,b0mneutmb(2,1))
      call b0(p,nmneut(2),mb2,q,b0mneutmb(2,2))
      call b0(p,nmneut(3),mb1,q,b0mneutmb(3,1))
      call b0(p,nmneut(3),mb2,q,b0mneutmb(3,2))
      call b0(p,nmneut(4),mb1,q,b0mneutmb(4,1))
      call b0(p,nmneut(4),mb2,q,b0mneutmb(4,2))
      

      call b1(p,dabs(mchargino(1)),dabs(mt1),q,b1mchmt(1,1))
      call b1(p,dabs(mchargino(1)),dabs(mt2),q,b1mchmt(1,2))
      call b1(p,dabs(mchargino(2)),dabs(mt1),q,b1mchmt(2,1))
      call b1(p,dabs(mchargino(2)),dabs(mt2),q,b1mchmt(2,2))

      call b0(p,dabs(mchargino(1)),mt1,q,b0mchmt(1,1))
      call b0(p,dabs(mchargino(1)),mt2,q,b0mchmt(1,2))
      call b0(p,dabs(mchargino(2)),mt1,q,b0mchmt(2,1))
      call b0(p,dabs(mchargino(2)),mt2,q,b0mchmt(2,2))

      
C---------------------------------------------------------------------

      stopterm = 0.d0
      btmtermt = 0.d0
      toptermt = 0.d0
      charginoterm = 0.d0

      stopterm =  ((g3**2.d0)/(12.d0 * pi*pi)) * (b1M3mb1 + b1M3mb2 - 
     $     (dsin(2.d0*thetabz1) * (M3t/mb) * (b0M3mb1 - b0M3mb2))) 

c$$$      print*,"b1m3b1 = ", b1M3mb1
c$$$      print*,"b1m3b2 = ", b1M3mb2
c$$$
c$$$      print*,"b0m3b1 = ", b0M3mb1
c$$$      print*,"b0m3b2 = ", b0M3mb2
c$$$
c$$$      print*,"p = ", p, "q = ", q
c$$$      print*,"g3 = ",g3,"M3t = ", M3t
c$$$      print*,"stopterm, mb = ", stopterm, mb, thetabz1

      btmtermt = 0.5d0 * (ydRG(3,3))**2.d0 *
     $     (((calpha**2.d0)*(b1mbmHu + b0mbmHu)) +
     $     ((salpha**2.d0)*(b1mbmh + b0mbmh)) + 
     $     ((sinbeta**2.d0)*(b1mbmA - b0mbmA)) +
     $     ((cosbeta**2.d0)*(b1mbMz - b0mbMz))) +
     $     (g/costhw)**2.d0 * ((gdL*gdL + gdR*gdR)*b1mbMz + 
     $     4.d0*gdL*gdR*b0mbMz) + 
     $     ((yuRG(3,3)*sinbeta)**2.d0 * (b0mtmhpm - b0mtmw)) 


      toptermt = 0.5d0*((yuRG(3,3)*cosbeta)**2.d0 + 
     $     (ydRG(3,3)*sinbeta)**2.d0)*b1mtmHpm +
     $     0.5d0*(g*g + (yuRG(3,3)*sinbeta)**2.d0 + 
     $     (ydRG(3,3)*cosbeta)**2.d0)*b1mtMw 


C----------------------------------------------------------------------------

      aPsicstl(1) = g
      aPsicstl(2) = 0.d0
      
      aPsicstr(1) = 0.d0
      aPsicstr(2) = -yuRG(3,3)

      bPsicstr(1) = 0.d0
      bPsicstr(2) = 0.d0

      bPsicstl(1) = 0.d0
      bPsicstl(2) = -ydRG(3,3)
!---
  
      aChicstl(1) = OCR(1,1)*aPsicstl(1) + OCR(1,2)*aPsicstl(2)
      aChicstl(2) = OCR(2,1)*aPsicstl(1) + OCR(2,2)*aPsicstl(2)

      bChicstl(1) = OCL(1,1)*bPsicstl(1) + OCL(1,2)*bPsicstl(2)
      bChicstl(2) = OCL(2,1)*bPsicstl(1) + OCL(2,2)*bPsicstl(2)
 
      aChicstr(1) = OCR(1,1)*aPsicstr(1) + OCR(1,2)*aPsicstr(2)
      aChicstr(2) = OCR(2,1)*aPsicstr(1) + OCR(2,2)*aPsicstr(2)

      bChicstr(1) = OCL(1,1)*bPsicstr(1) + OCL(1,2)*bPsicstr(2)
      bChicstr(2) = OCL(2,1)*bPsicstr(1) + OCL(2,2)*bPsicstr(2)

      
        
      call rmat2d(thetatz1,rthetat)

      

      loopchi: DO i = 1, 2
      


      tt(1) = aChicstl(i)
      tt(2) = aChicstr(i)

      t1(1) = rthetat(1,1)*tt(1) + rthetat(1,2)*tt(2)
      t1(2) = rthetat(2,1)*tt(1) + rthetat(2,2)*tt(2)

      tt(1) = bChicstl(i)
      tt(2) = bChicstr(i)      
      
      t2(1) = rthetat(1,1)*tt(1) + rthetat(1,2)*tt(2)
      t2(2) = rthetat(2,1)*tt(1) + rthetat(2,2)*tt(2)


!----------------------------------------------------------------------------

      loopchj: DO j =1, 2 

      fChBstop(i, j) = 0.d0
      gChBstop(i, j) = 0.d0

      aChBstop(i, j) = t1(j)
      bChBstop(i, j) = t2(j)

      fChBstop(i, j) = (aChBstop(i, j))**2.d0 + (bChBstop(i, j))**2.d0
      gChBstop(i, j) = 2.d0*(aChBstop(i, j)*bChBstop(i, j)) 
      

      charginoterm = charginoterm + (fChBstop(i, j)*b1mchmt(i, j) +
     $     gChBstop(i, j) * (mchargino(i)/mb) * b0mchmt(i, j))*0.5d0


      

      enddo loopchj
      enddo loopchi




!-------------------------------------------------------------------
C          Neutralinos



      aPsi0sbr(1) = gp/(dsqrt(2.d0)*3.d0) * 2.d0
      aPsi0sbr(2) = 0.d0
      aPsi0sbr(3) = 0.d0
      aPsi0sbr(4) = 0.d0

      bPsi0sbl(1) = gp/ (dsqrt(2.d0)* 3.d0)
      bPsi0sbl(2) = -dsqrt(2.d0)*g*0.5d0
      bPsi0sbl(3) = 0.d0
      bPsi0sbl(4) = 0.d0

      aPsi0sbl(1) = 0.d0
      aPsi0sbl(2) = 0.d0
      aPsi0sbl(3) = ydRG(3,3)
      aPsi0sbl(4) = 0.d0

      bPsi0sbr(1) = 0.d0
      bPsi0sbr(2) = 0.d0
      bPsi0sbr(3) = ydRG(3,3)
      bPsi0sbr(4) = 0.d0


!---------------------------------------------------------------
      looponi: DO i = 1, 4

      aChi0sbl(i) = 0.d0
      bChi0sbl(i) = 0.d0
      aChi0sbr(i) = 0.d0
      bChi0sbr(i) = 0.d0

      looponj: DO j = 1, 4


      aChi0sbl(i) = aChi0sbl(i) + ON(i,j)*aPsi0sbl(j)
      bChi0sbl(i) = bChi0sbl(i) + ON(i,j)*bPsi0sbl(j)
      aChi0sbr(i) = aChi0sbr(i) + ON(i,j)*aPsi0sbr(j)
      bChi0sbr(i) = bChi0sbr(i) + ON(i,j)*bPsi0sbr(j)

      ENDDO looponj
      ENDDO looponi
!--------------------------------------------------------------

      call rmat2d(thetabz1,rthetab)
      tneutterm = 0.d0
      loopneui: DO i = 1, 4

      tt(1) = aChi0sbl(i)
      tt(2) = aChi0sbr(i)      
      
      t1(1) = 0.d0
      t1(2) = 0.d0

      t1(1) = rthetab(1,1)*tt(1) + rthetab(1,2)*tt(2)
      t1(2) = rthetab(2,1)*tt(1) + rthetab(2,2)*tt(2)

      tt(1) = bChi0sbl(i)
      tt(2) = bChi0sbr(i)

      t2(1) = 0.d0
      t2(2) = 0.d0

      t2(1) = rthetab(1,1)*tt(1) + rthetab(1,2)*tt(2)
      t2(2) = rthetab(2,1)*tt(1) + rthetab(2,2)*tt(2)
     
      loopneuj: DO j = 1, 2

      neutm(i,j) = 0.d0
      aNeutsb(i, j) = t1(j)
      bNeutsb(i, j) = t2(j)

      fNeutsb(i, j) = 0.d0
      gNeutsb(i, j) = 0.d0

      fNeutsb(i, j) = (aNeutsb(i, j))**2.d0 + 
     $     (bNeutsb(i, j))**2.d0
      
      gNeutsb(i, j) = 2.d0*(aNeutsb(i, j)*bNeutsb(i, j)) 
      
      neutm(i, j) = (fNeutsb(i, j) * b1mneutmb(i, j) + 
     $     gNeutsb(i, j) * (mneut(i)/mb) * b0mneutmb(i, j))*0.5d0

      tneutterm = neutm(i,j) + tneutterm

    
      ENDDO loopneuj
      ENDDO loopneui


!--------------MSbar---> DRbar

        mbdrbar = mbMZmsbar*(
     $        (1.d0 - (g3**2.d0/(12.d0 * pi * pi)) - 
     $     (23.d0 * (g3**4.d0)/(72.d0 *16.d0 *pi*pi* pi * pi)) + 
     $     (3.d0 * g**2.d0/(128.d0 * pi * pi)) + 
     $     (13.d0 * gp**2.d0/(1152.d0 * pi *pi) * (5.d0/3.d0))))

!------------------------------------------------------------------------------

      correction = (- stopterm - (toptermt + btmtermt +
     $     tneutterm + charginoterm) /(16.d0*pi*pi))

c$$$      print*,"in btmcr stopterm, topterm, btmterm, neutterm, 
c$$$     $     charterm, mb = ", stopterm, toptermt, btmtermt, tneutterm, 
c$$$     $     charginoterm, mb
c$$$      
c$$$      print*,"mbmzmsbar = ", mbmzmsbar,mbdrbar,mbmzdrbar
c$$$      print*,"mb cor = ", mbdrbar/(1 + correction), 
c$$$     $     mbmzdrbar/(1 + correction)

      mbcor = correction


      RETURN

      END SUBROUTINE bottomcor

C====================================================================================

      SUBROUTINE taucor(p,q,g,gp,M3t,mtau,tanbeta,yeRG,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,
     $     mhu0sq,mhpmsq,mA0sq,ON,OCL,
     $     OCR,sinsqtheff,mtaucor)

      IMPLICIT NONE 
      
      INTEGER i,j
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),correction,mtaucor
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,M3t
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      DOUBLE PRECISION yeRG(3,3)
      
      double precision mTau
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2

      DOUBLE PRECISION sinsqtheff
      DOUBLE PRECISION ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta,alpha

      DOUBLE PRECISION b1mneutmtau(4,2),b0mneutmtau(4,2)

      DOUBLE PRECISION sinsqthw,g
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION beta
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetatz,thetabz,thetatauz

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetacz,thetasz,thetamuz
      DOUBLE PRECISION thetauz,thetadz,thetaez

      DOUBLE PRECISION mu,rthetatau(2,2)
      DOUBLE PRECISION stopterm
      DOUBLE PRECISION mchargino(2),tt(2),t1(2),t2(2),nmneut(4)

      DOUBLE PRECISION b1M3mtau1,b1M3mtau2,b0M3mtau1,b0M3mtau2
      DOUBLE PRECISION toptermt,btmtermt,tneutterm,charginoterm
      
      DOUBLE PRECISION b10Mw,b10mHpm
      DOUBLE PRECISION b1mtaumHu,b1mtaumh,b1mtauMz,b1mtaumA
      DOUBLE PRECISION b0mtaumHu,b0mtaumh,b0mtauMz,b0mtaumA
      DOUBLE PRECISION b1mch10,b1mch20
      DOUBLE PRECISION b0mch10,b0mch20

      DOUBLE PRECISION aPsicTauSnul(2), bPsicTauSnul(2)
      data aPsicTauSnul/ 2 * 0.d0/, bPsicTauSnul/ 2 * 0.d0/

      DOUBLE PRECISION aChicTauSnul(2), bChicTauSnul(2)
      data aChicTauSnul/ 2 * 0.d0/, bChicTauSnul/ 2 * 0.d0/

      DOUBLE PRECISION fChiTauSnu(2), gChiTauSnu(2) 
      data fChiTauSnu/ 2 * 0.d0/, gChiTauSnu/ 2 * 0.d0/

      DOUBLE PRECISION aPsi0TauStaur(4), bPsi0TauStaur(4)
      data aPsi0TauStaur/ 4 *0.d0/,bPsi0TauStaur/ 4 * 0.d0/

      DOUBLE PRECISION aPsi0TauStaul(4), bPsi0TauStaul(4)
      data aPsi0TauStaul/ 4 *0.d0/,bPsi0TauStaul/ 4 * 0.d0/

      DOUBLE PRECISION aChi0TauStaul(4), bChi0TauStaul(4)
      data aChi0TauStaul/ 4 *0.d0/,bChi0TauStaul/ 4 * 0.d0/

      DOUBLE PRECISION aChi0TauStaur(4), bChi0TauStaur(4)
      data aChi0TauStaur/ 4 *0.d0/,bChi0TauStaur/ 4 * 0.d0/

      DOUBLE PRECISION aNeutTauStau(4,2), bNeutTauStau(4,2)
      data aNeutTauStau/ 8 * 0.d0/, bNeutTauStau/ 8 * 0.d0/

      DOUBLE PRECISION fNeutTauStau(4, 2), gNeutTauStau(4, 2)
      data fNeutTauStau/ 8 * 0.d0/, gNeutTauStau/ 8 * 0.d0/

      double precision mtaupoledrbar,mTauMZmsbar,mTauMZdrbar,mtaup

      DOUBLE PRECISION mbpole, mtaupole, Mtpole, MWpole, MZpole
      double precision MW,MZ,pi
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole
      common/sfmixing_mz/thetatz,thetabz,thetatauz,thetacz,thetasz,
     $     thetamuz,thetauz,thetadz,thetaez
      common/mtaurmass/ mtaupoledrbar,mTauMZmsbar,mTauMZdrbar
      
      EXTERNAL b1,b0,mat3prod4d,dag4d,dag2d,mat3prod2d,theta,rmat2d,
     $     b0lt

      include 'stdinputs.h'
!---------------------------------------------------------------------------
      
      pi = 4.d0*datan(1.d0)
      MW = MWpole
      MZ = MZpole
      mu = mUQ

      sinsqthw = sinsqtheff 
      sinthw = dsqrt(sinsqtheff)
      thw = dasin(sinthw)
      costhw = dcos(thw)
      cos2thw = dcos(2.d0*thw)
      

      gnuL = 0.5d0
      guL = 0.5d0 - 2.d0*sinsqthw/ 3.d0 
      gdL = -0.5d0 + sinsqthw/3.d0
      geL = -0.5d0 + sinsqthw
      guR =  2.d0*sinsqthw/3.d0 
      gdR = -sinsqthw/3.d0
      geR = -sinsqthw 
      yuL = 1.d0/3.d0 
      yuR = -4.d0/3.d0 
      ydL = 1.d0/3.d0 
      ydR = 2.d0/3.d0 
      yeL = -1.d0 
      yeR = 2.d0 
      ynuL = -1.d0 

C------------------------------------------------------------

      mTaup = mTauMZmsbar * (1.d0 - (3.d0/1.d0) * 
     $     (gp*gp - g*g)/128.d0)

!      print*,"mtaumzdrbar = ", mtaup, mtaumzdrbar

C------------------------------------------------------------
c$$$      muR = dsqrt(SUegg(6)) 
c$$$      muL = dsqrt(SUegg(5))
c$$$      mcR = dsqrt(SUegg(4))
c$$$      mcL = dsqrt(SUegg(3))
c$$$      mtR = dsqrt(SUegg(2))
c$$$      mtL = dsqrt(SUegg(1))
c$$$      
c$$$      mdR = dsqrt(SDegg(6)) 
c$$$      mdL = dsqrt(SDegg(5))
c$$$      msR = dsqrt(SDegg(4))
c$$$      msL = dsqrt(SDegg(3))
c$$$      mbR = dsqrt(SDegg(2))
c$$$      mbL = dsqrt(SDegg(1))
c$$$
c$$$      meR   = dsqrt(SLegg(6))
c$$$      meL   = dsqrt(SLegg(5))
c$$$      mmuR  = dsqrt(SLegg(4))
c$$$      mmuL  = dsqrt(SLegg(3))
c$$$      mtauR = dsqrt(SLegg(2))
c$$$      mtauL = dsqrt(SLegg(1))

      muR = dsqrt(SUegg(5)) 
      muL = dsqrt(SUegg(6))
      mcR = dsqrt(SUegg(3))
      mcL = dsqrt(SUegg(4))
      mtR = dsqrt(SUegg(1))
      mtL = dsqrt(SUegg(2))
      
      mdR = dsqrt(SDegg(5)) 
      mdL = dsqrt(SDegg(6))
      msR = dsqrt(SDegg(3))
      msL = dsqrt(SDegg(4))
      mbR = dsqrt(SDegg(1))
      mbL = dsqrt(SDegg(2))

      meR   = dsqrt(SLegg(5))
      meL   = dsqrt(SLegg(6))
      mmuR  = dsqrt(SLegg(3))
      mmuL  = dsqrt(SLegg(4))
      mtauR = dsqrt(SLegg(1))
      mtauL = dsqrt(SLegg(2))
      
      snu(1) = dsqrt(SNegg(1))
      snu(2) = dsqrt(SNegg(2))
      snu(3) = dsqrt(SNegg(3))

      nmneut(1) = dabs(Neg(1)) 
      nmneut(2) = dabs(Neg(2))
      nmneut(3) = dabs(Neg(3))
      nmneut(4) = dabs(Neg(4))

      mneut(1) = Neg(1) 
      mneut(2) = Neg(2)
      mneut(3) = Neg(3)
      mneut(4) = Neg(4)

      mchargino(1) = Ceg(1)
      mchargino(2) = Ceg(2)
 

!-----------------------------------------------------------------------

      costhetat = dcos(thetatz)
      sinthetat = dsin(thetatz)
      costhetab = dcos(thetabz)
      sinthetab = dsin(thetabz)
      costhetatau = dcos(thetatauz)
      sinthetatau = dsin(thetatauz)

      
      mt1 = mtL
      mt2 = mtR

      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR

c$$$      mt1 = mtR
c$$$      mt2 = mtL
c$$$
c$$$      mb1 = mbR
c$$$      mb2 = mbL
c$$$
c$$$      mtau1 = mtauR
c$$$      mtau2 = mtauL

C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))


      beta = datan(tanbeta)
      tan2beta = dtan(2.d0 * beta)
      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/
     $     (mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta = dcos(beta)
      cos2beta = dcos(2.d0 * beta)
      sinbeta = dsin(beta)
      sin2beta = dsin(2.d0 * beta)
      correction = 0.d0

C--------------------------------------------------------------------


      call b1(p,M3t,mtau1,q,b1M3mtau1)
      call b1(p,M3t,mtau2,q,b1M3mtau2)

      call b1(p,mtau,mHu0,q,b1mtaumHu)
      call b1(p,mtau,mh0,q,b1mtaumh)

      call b1(p,mtau,mA0,q,b1mtaumA)
      call b1(p,mtau,MZ,q,b1mtauMz)

!      call b1(p,mtau,mHpm,q,b10mHpm)
!      call b1(p,mtau,MW,q,b10Mw)

      call b1(p,2.d-5,mHpm,q,b10mHpm)
      call b1(p,2.d-5,MW,q,b10Mw)

      call b0(p,M3t,mtau1,q,b0M3mtau1)
      call b0(p,M3t,mtau2,q,b0M3mtau2)

      call b0(p,mtau,mHu0,q,b0mtaumHu)
      call b0(p,mtau,mh0,q,b0mtaumh)

      call b0(p,mtau,mA0,q,b0mtaumA)
      call b0(p,mtau,MZ,q,b0mtauMz)


C---------------------------------------------------------------------

c$$$      print*,"p = ", p, "q = ", q 
c$$$      print*,"mtau = ", mtau, "mtaupole = ", mtaupole
c$$$      print*,"thetatau = ", thetatauz

      stopterm = 0.d0
      btmtermt = 0.d0
      toptermt = 0.d0


      btmtermt =  0.5d0*yeRG(3,3)**2.d0*(((calpha**2.d0)*
     $     (b1mtaumHu + b0mtaumHu)) +    
     $     ((salpha**2.d0)*(b1mtaumh + b0mtaumh)) + 
     $     ((sinbeta**2.d0)*(b1mtaumA - b0mtaumA)) +
     $     ((cosbeta**2.d0)*(b1mtauMz - b0mtauMz))) +
     $     ((g/costhw)**2.d0*((geL*geL + geR*geR)*b1mtauMz + !<---------check, gtl,gtr: def!!---------------
     $     4.d0*geL*geR*b0mtauMz))

      

      toptermt =  0.5d0*(((yeRG(3,3)*sinbeta)**2.d0)*b10mHpm +
     $     (g*g + (yeRG(3,3)*cosbeta)**2.d0)*b10Mw)


C----------------------------------------------------------------------------
      
C     Neutralinos


      aPsi0TauStaur(1) = gp*dsqrt(2.d0)
      aPsi0TauStaur(2) = 0.d0
      aPsi0TauStaur(3) = 0.d0
      aPsi0TauStaur(4) = 0.d0

      aPsi0TauStaul(1) = 0.d0
      aPsi0TauStaul(2) = 0.d0
      aPsi0TauStaul(3) = yeRG(3,3)
      aPsi0TauStaul(4) = 0.d0

      bPsi0TauStaul(1) = -gp/dsqrt(2.d0)
      bPsi0TauStaul(2) = -g/dsqrt(2.d0)
      bPsi0TauStaul(3) = 0.d0
      bPsi0TauStaul(4) = 0.d0

      bPsi0TauStaur(1) = 0.d0
      bPsi0TauStaur(2) = 0.d0
      bPsi0TauStaur(3) = yeRG(3,3)
      bPsi0TauStaur(4) = 0.d0

      tneutterm = 0.d0
!--------------------------------------------------------------

      call b1(p,nmneut(1),mtau1,q,b1mneutmtau(1,1))
      call b1(p,nmneut(1),mtau2,q,b1mneutmtau(1,2))
      call b1(p,nmneut(2),mtau1,q,b1mneutmtau(2,1))
      call b1(p,nmneut(2),mtau2,q,b1mneutmtau(2,2))
      call b1(p,nmneut(3),mtau1,q,b1mneutmtau(3,1))
      call b1(p,nmneut(3),mtau2,q,b1mneutmtau(3,2))
      call b1(p,nmneut(4),mtau1,q,b1mneutmtau(4,1))
      call b1(p,nmneut(4),mtau2,q,b1mneutmtau(4,2))

      call b0(p,nmneut(1),mtau1,q,b0mneutmtau(1,1))
      call b0(p,nmneut(1),mtau2,q,b0mneutmtau(1,2))
      call b0(p,nmneut(2),mtau1,q,b0mneutmtau(2,1))
      call b0(p,nmneut(2),mtau2,q,b0mneutmtau(2,2))
      call b0(p,nmneut(3),mtau1,q,b0mneutmtau(3,1))
      call b0(p,nmneut(3),mtau2,q,b0mneutmtau(3,2))
      call b0(p,nmneut(4),mtau1,q,b0mneutmtau(4,1))
      call b0(p,nmneut(4),mtau2,q,b0mneutmtau(4,2))

!---------------------------------------------------------------

      looponi: DO i = 1, 4

      aChi0TauStaul(i) = 0.d0
      bChi0TauStaul(i) = 0.d0
      aChi0TauStaur(i) = 0.d0
      bChi0TauStaur(i) = 0.d0

      looponj: DO j = 1, 4

      aChi0TauStaul(i) = aChi0TauStaul(i) + ON(i,j)*aPsi0TauStaul(j)
      bChi0TauStaul(i) = bChi0TauStaul(i) + ON(i,j)*bPsi0TauStaul(j)
      aChi0TauStaur(i) = aChi0TauStaur(i) + ON(i,j)*aPsi0TauStaur(j)
      bChi0TauStaur(i) = bChi0TauStaur(i) + ON(i,j)*bPsi0TauStaur(j)


      ENDDO looponj
      ENDDO looponi
!--------------------------------------------------------------

      call rmat2d(thetatauz,rthetatau)

      loopneui: DO i = 1, 4

      tt(1) = aChi0TauStaul(i)
      tt(2) = aChi0TauStaur(i)      

      t1(1) = rthetatau(1,1)*tt(1) + rthetatau(1,2)*tt(2)
      t1(2) = rthetatau(2,1)*tt(1) + rthetatau(2,2)*tt(2)

      tt(1) = bChi0TauStaul(i)
      tt(2) = bChi0TauStaur(i)

      t2(1) = rthetatau(1,1)*tt(1) + rthetatau(1,2)*tt(2)
      t2(2) = rthetatau(2,1)*tt(1) + rthetatau(2,2)*tt(2)
     
      loopneuj: DO j = 1, 2

      aNeutTauStau(i,j) = t1(j)
      bNeutTauStau(i,j) = t2(j)
 
      fNeutTauStau(i,j) = (aNeutTauStau(i,j))**2.d0 + 
     $                   (bNeutTauStau(i,j))**2.d0

      gNeutTauStau(i,j) = 2.d0*(aNeutTauStau(i,j)*bNeutTauStau(i,j))  

     
      tneutterm = tneutterm + (fNeutTauStau(i,j)*b1mneutmtau(i,j)+ 
     $        gNeutTauStau(i,j)*(mneut(i)/mtaup)*b0mneutmtau(i,j))*
     $        0.5d0

      ENDDO loopneuj
      ENDDO loopneui


C------------------------------------------------------------------------------
      
C     Chargino Terms

      aPsicTauSnul(1) = g
      aPsicTauSnul(2) = 0.d0

      bPsicTauSnul(1) = 0.d0
      bPsicTauSnul(2) = -yeRG(3,3)

C        Mass eignebasis of charginos

      aChicTauSnul(1) = OCR(1,1)*aPsicTauSnul(1) + 
     $                  OCR(1,2)*aPsicTauSnul(2)

      aChicTauSnul(2) = OCR(2,1)*aPsicTauSnul(1) + 
     $                  OCR(2,2)*aPsicTauSnul(2)

      bChicTauSnul(1) = OCL(1,1)*bPsicTauSnul(1) + 
     $                  OCL(1,2)*bPsicTauSnul(2)

      bChicTauSnul(2) = OCL(2,1)*bPsicTauSnul(1) +
     $                   OCL(2,2)*bPsicTauSnul(2)

!---------------------------------------------------------------------

      call b1(p,dabs(mchargino(1)),snu(1),q,b1mch10)
      call b1(p,dabs(mchargino(2)),snu(1),q,b1mch20)

      call b0(p,dabs(mchargino(1)),snu(1),q,b0mch10)
      call b0(p,dabs(mchargino(2)),snu(1),q,b0mch20)
          
!----------------------------------------------------------------------

      charginoterm = 0.d0
    
      loopchginoi: DO i = 1,2
  
      fChiTauSnu(i) = 0.d0
      gChiTauSnu(i) = 0.d0      
      fChiTauSnu(i) = (aChicTauSnul(i))**2.d0 + (bChicTauSnul(i))**2.d0
      gChiTauSnu(i) = 2.d0*(aChicTauSnul(i)*bChicTauSnul(i))
      
      ENDDO loopchginoi
      
      charginoterm = (fChiTauSnu(1)*b1mch10 + 
     $     gChiTauSnu(1)*(mchargino(1)/mtaup)*b0mch10)* 0.5d0


      charginoterm  = charginoterm + (fChiTauSnu(2)*b1mch20 + 
     $     gChiTauSnu(2)*(mchargino(2)/mtaup)*b0mch20) * 0.5d0

!----------------------------------------------------------

      correction = (charginoterm + btmtermt + toptermt + tneutterm)
     $     /((16.d0*pi*pi))  

c$$$      print*,"neut, char, higgs = ", tneutterm ,
c$$$     $     charginoterm ,(btmtermt+toptermt)
c$$$
c$$$      print*,"mtaupole, mtau_mz = ", mtaupole,mtau,mtaup
      
      mtaucor =  correction    
      
      RETURN

      END SUBROUTINE taucor

!============================================================================================================
