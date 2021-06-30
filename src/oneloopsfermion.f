****f* SuSeFLAV/oneloopsfermion.f 
*  NAME
*    oneloopsfermion
*  SYNOPSIS
*    One loop correction to sfermions. 
*  FUNCTION
*     Computes self energy for all sfermions at a given energy scale and external momenta
*  INPUTS
*     p                              - External momentum  
*     q                              - Energy scale
*     g,gp,g3                        - Gauge couplings(g = em, gp = weak, 
*                                      g3 = strong)                                     
*     mt,mb,mtau                     - pole masses of top, bottom and tau
*     mSQRG,mSDRG,mSURG,mSLRG,mSERG  - (3 X 3) mass matrix definition
*     yuRG,ydRG,yeRG                 - (3 X 3)Yukawas
*     AURG,ADRG,AERG                 - (3 X 3)Trilinear couplings
*     pizzT,piwwT                    - self energy of W and Z bosons at M_z
*     mh0sq,mhu0sq,mhpmsq,mA0sq      - physicsal higgs mass squared 
*     modmu                          - modulus of the \mu paramter 
*     vev1,vev2                      - vacuum expectation values of the two higgs 
*                                      doublet fields
*     M3t                            - Gaugino mass at msusy
*     tanbeta                        - the ratio of the vevs of the two–Higgs 
*                                      doublet fields.
*     SUegg   =  6 eigenvalues of UP-Squark mass matrix.
*     SDegg   =  6 eigenvalues of Down-Squark mass matrix.
*     Slegg   =  6 eigenvalues of slepton mass matrix.
*     SNegg   =  3 eigenvalues of sneutrino mass matrix.
*     ON      =  (4 X 4) orthogonal matrix such that 
*                ON.MNeut.Transpose[ON] = Diag[MNeut] 
*     Neg     =  4 singular values (descending order) of the 
*                 Neutralino mass matrix. 
*
*     OCR, OCL =  (2 X 2) orthogonal matrices such that 
*                 MChar = OCL.Diag[MChar].Transpose[OCR]
*     Ceg     =   2 singular values of the chargino Mass Matrix
*
*  RESULT
*   S_eg     =  2 eigenvalues of the corrected mass matrix( L and R components).
*
*  EXAMPLE
*
*      SUBROUTINE pistop(p,q,g,gp,g3,mt,mb,tanbeta,mSQRG,
*     $     mSURG,yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
*     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
*     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3t,delpist,STeg)
*    
*      SUBROUTINE pisbottom(p,q,g,gp,g3,mt,mb,tanbeta,mSQRG,mSDRG,
*     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
*     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
*     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3t,delpisbtm,SBeg)
*      
*      SUBROUTINE pisdown(p,q,g,gp,g3,tanbeta,mSQRG,mSDRG,
*     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
*     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
*     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3t,SDeg)
*
*      SUBROUTINE pisstrange(p,q,g,gp,g3,tanbeta,mSQRG,mSDRG,
*     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
*     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
*     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3t,delpisst,SSTeg)
* 
*       
*  NOTES
*    1. q, the energy scale at which the corrections are added = msusy.
*    2. Running values of gauge couplings( Rge output) are used. 
*    3. Pole masses: DRbar scheme is followed.
*    4. Conventions and notations followed are that of BPMZ.
*  BUGS
*    ---
*  SEE ALSO
*    DSYEV - Diagonalizing Routine.
*    LAPACK
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C

C=====================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     1. Removed dependency on stddef.dat from pistop,pisbottom,pisdown,pisstrange,pistaul,
C     pismul,pisel,pischarm,pisupq,pisnutaul,pisnumul,pisnuel.
C     2. Removed unused variables from pistop,pisbottom,pisdown,pisstrange,pistaul,
C     pismul,pisel,pischarm,pisupq,pisnutaul,pisnumul,pisnuel.
C     3. Corrected some terms in pistop,pisbottom,pisdown,pisstrange,pistaul,
C     pismul,pisel,pischarm,pisupq,pisnutaul,pisnumul,pisnuel.
C     4. The couplings nomenclature changing is not yet done.
C     5. All the variables in all subroutines has not yet initialized to zero.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

C======================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C     pistop is checked on 20/05/2010 @ 14:30.
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      
      SUBROUTINE pistop(p,q,g,gp,g3,mt,mb,tanbeta,mSQRG,
     $     mSURG,yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3t,delpist,STeg)
      
      
      IMPLICIT NONE 

      INTEGER i,j,info,AOK,lwork
      parameter(lwork = 35)
      DOUBLE PRECISION work(lwork)
      DOUBLE PRECISION tanbeta,mSQRG(3,3),mSURG(3,3)
      DOUBLE PRECISION AURG(3,3),ADRG(3,3)
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yuRG(3,3),ydRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      
      double precision mT, mB
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR,MSQU3(2,2)      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g,g3,M3t
      
      double precision beta,sgnmu,vev1,vev2      
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION pist(2,2)
      
      DOUBLE PRECISION higgsterm(2,2),higgsc(2)
      data higgsterm/ 4 * 0.d0/,higgsc/ 2 * 0.d0/
      DOUBLE PRECISION dnu(4),dnd(4),cn(4)


      DOUBLE PRECISION lstlstlR(4,2),lstlst12(4,2)
      data lstlstlR/ 8 * 0.d0/,lstlst12/ 8 * 0.d0/
      
      DOUBLE PRECISION lstrstlR(4,2),lstRst12(4,2)
      data lstrstlR/ 8 * 0.d0/,lstRst12/ 8 * 0.d0/

      
      DOUBLE PRECISION aPsi0str(4), bPsi0str(4), aPsi0stl(4)
      data aPsi0str/ 4 * 0.d0/, bPsi0str/ 4 * 0.d0/,
     $     aPsi0stl/ 4 * 0.d0/


      DOUBLE PRECISION bPsi0stl(4)
      data bPsi0stl/ 4 * 0.d0/
      
      DOUBLE PRECISION aChi0stl(4), bChi0stl(4), aChi0str(4)
      data aChi0stl/ 4 * 0.d0/, bChi0stl/ 4 * 0.d0/, 
     $     aChi0str/ 4 * 0.d0/
      
      DOUBLE PRECISION bChi0str(4)
      data bChi0str/ 4 * 0.d0/

      
      DOUBLE PRECISION gChi0tstll(4), fChi0tstll(4)
      data gChi0tstll/ 4 * 0.d0/, fChi0tstll/ 4 * 0.d0/

      DOUBLE PRECISION gChi0tstlr(4), fChi0tstlr(4)
      data gChi0tstlr/ 4 * 0.d0/, fChi0tstlr/ 4 * 0.d0/


      DOUBLE PRECISION gChi0tstrr(4), fChi0tstrr(4)
      data gChi0tstrr/ 4 * 0.d0/, fChi0tstrr/ 4 * 0.d0/

      DOUBLE PRECISION bPsicbstl(2), bPsicbstr(2), aPsicbstl(2)
      data bPsicbstl/ 2 * 0.d0/,bPsicbstr/ 2 * 0.d0/,
     $     aPsicbstl/ 2 * 0.d0/
      
      DOUBLE PRECISION aPsicbstr(2)
      data aPsicbstr/ 2 * 0.d0/


      DOUBLE PRECISION aChicbstr(2), aChicbstl(2)
      data aChicbstr/ 2 * 0.d0/,
     $     aChicbstl/ 2 * 0.d0/

      DOUBLE PRECISION bChicbstl(2),bChicbstr(2)
      data bChicbstl/ 2 * 0.d0/,bChicbstr/ 2 * 0.d0/

      DOUBLE PRECISION aChBStop(2,2),bChBStop(2,2)
      data aChBStop/ 4 * 0.d0/,bChBStop/ 4 * 0.d0/

      DOUBLE PRECISION gtterm(2,2), cstop(2,2), csbtm(2,2) 
      data gtterm/ 4 * 0.d0/,cstop/ 4 * 0.d0/,csbtm/ 4 * 0.d0/

      DOUBLE PRECISION chargino(2,2), neutralino(2,2)
      data chargino/ 4 * 0.d0/, neutralino/ 4 * 0.d0/

      DOUBLE PRECISION fChbstll(2), gChbstll(2) 
      data fChbstll/ 2 * 0.d0/,gChbstll/ 2 * 0.d0/

      DOUBLE PRECISION fChbstlr(2),gChbstlr(2) 
      data fChbstlr/ 2 * 0.d0/,gChbstlr/ 2 * 0.d0/

      DOUBLE PRECISION fChbstrr(2), gChbstrr(2)
      data fChbstrr/ 2 * 0.d0/,gChbstrr/ 2 * 0.d0/

      DOUBLE PRECISION intm1(2)  , intm2(2)
      data intm1/ 2 * 0.d0/, intm2/ 2 * 0.d0/

      DOUBLE PRECISION lHcstrsblr(2,2),lHcstrsb12(2,2)
      data lHcstrsblr/ 4 * 0.d0/, lHcstrsb12/ 4 * 0.d0/

      DOUBLE PRECISION lHcstlsblR(2,2),lHcstlsb12(2,2)
      data lHcstlsblR/ 4 * 0.d0/, lHcstlsb12/ 4 * 0.d0/

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION lHstlst12(4,2),lHstrst12(4,2)
      data lHstlst12/ 8 * 0.d0/, lHstrst12/ 8 * 0.d0/

      
      DOUBLE PRECISION rthetat(2,2),rthetab(2,2),ralpha(2,2)

      DOUBLE PRECISION ggt,ft10,ft20,a0t1,a0t2,b0M3mt,a0mb1,a0mb2
      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz

      DOUBLE PRECISION b0h0mstop(4,2),b0hcmbtm(2,2)
      data b0h0mstop/ 8 * 0.d0/,b0hcmbtm/ 4 * 0.d0/

      DOUBLE PRECISION a0md1,a0md2,a0ms1,a0ms2,a0mtau1,a0mtau2,a0snu2
      DOUBLE PRECISION a0me1,a0me2,a0mmu1,a0mmu2,a0snu1,a0snu3
      DOUBLE PRECISION a0mc1,a0mc2,a0mu1,a0mu2
      
      DOUBLE PRECISION fmb1MW,fmb2MW,fmt1MZ,fmt2MZ
      DOUBLE PRECISION gmchargino1mb,gmchargino2mb
      DOUBLE PRECISION b0mchargino1mb,b0mchargino2mb
      DOUBLE PRECISION b0mneut1mt,b0mneut2mt,b0mneut3mt,b0mneut4mt
      DOUBLE PRECISION gmneut1mt,gmneut2mt,gmneut3mt,gmneut4mt
      
      DOUBLE PRECISION STeg(2),stmix(2,2)

      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2
      DOUBLE PRECISION nmneut(4),sinsqthw_susy, delpist(2,2)

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW

      double precision :: mwrun,mzrun
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole
      
!------------
      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

      common/sinsq_susy/sinsqthw_susy

!------------
      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d
      

      include 'stdinputs.h'
!      include 'stddef.dat'

C------------------------------------------------------------
      
      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      MZ = MZpole

      MZrun = dsqrt(((g*g) + 
     $     (gp*gp) * (3.d0/5.d0)) *
     $     (vev1**2.d0 + vev2**2.d0))/2.d0

      MWrun = dsqrt((g*g)*(vev1**2.d0 + vev2**2.d0))/2.d0

c$$$      MW = MWrun
c$$$      MZ = MZrun

c$$$      print*,"mwrun,mzrun = ", mwrun,mzrun
c$$$      print*,"mw,mz = ", mw,mz

      STeg(1) = 0.d0
      STeg(2) = 0.d0

      
      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw = dsqrt(sinsqthw)
      thw = dasin(sinthw)
      costhw = dcos(thw)
      cos2thw = dcos(2.d0*thw)

      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))
      
      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))

       
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

      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)

      mt1 = mtL 
      mt2 = mtR

      mb1 = mbL 
      mb2 = mbR

      mtau1 = mtauL 
      mtau2 = mtauR

C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta     = datan(tanbeta)
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/
     $     (mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta  = dcos(datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)

!-----------------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0

!----------------------------------------------------------------
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

!---------------------------------------------------------------    

      pist(1,1) = 0.d0
      pist(2,1) = 0.d0
      pist(1,2) = 0.d0
      pist(2,2) = 0.d0
      delpist(1,1) = 0.d0
      delpist(2,1) = 0.d0
      delpist(1,2) = 0.d0
      delpist(2,2) = 0.d0

C---------------------------------------------------------------    
      
      lstlstlR(1,1) = g*MZ*guL*cosbeta/costhw
      lstlstlR(1,2) = - yuRG(3,3)*sgnmu*modmu/dsqrt(2.d0)
      
      lstlstlR(2,1) = (- g*MZ*guL*sinbeta/costhw) + 
     $     (yuRG(3,3)**2.d0) * vev2
      lstlstlR(2,2) = (yuRG(3,3) * AURG(3,3))/dsqrt(2.d0)
      
      lstlstlR(3,1) = 0.d0
      lstlstlR(3,2) = ((-sgnmu*modmu*cosbeta*yuRG(3,3) + 
     $     yuRG(3,3) * AURG(3,3) * sinbeta))/dsqrt(2.d0)
      
      lstlstlR(4,1) = 0.d0
      lstlstlR(4,2) = -((-sgnmu*modmu*sinbeta*yuRG(3,3) - 
     $     yuRG(3,3) * AURG(3,3) * cosbeta))/dsqrt(2.d0)
      

      lstrstlR(1,1) = lstlstlR(1,2)
      lstrstlR(1,2) = g*MZ*guR*cosbeta/costhw

      lstrstlR(2,1) = lstlstlR(2,2)
      lstrstlR(2,2) = -g*MZ*guR*sinbeta/costhw + 
     $     (yuRG(3,3)**2.d0) * vev2

      lstrstlR(3,1) = -lstlstlR(3,2)
      lstrstlR(3,2) = 0.d0
      
      lstrstlR(4,1) = -lstlstlR(4,2)
      lstrstlR(4,2) = 0.d0

!--------------------------------------------------------------------------------
      
      loopmixst: do i = 1, 4

      intm1(1) = lstlstlR(i,1)
      intm1(2) = lstlstlR(i,2)

      call rmat2d(thetat,rthetat)
      
      intm2(1) = rthetat(1,1)*intm1(1) + rthetat(1,2)*intm1(2)
      intm2(2) = rthetat(2,1)*intm1(1) + rthetat(2,2)*intm1(2)

      lstlst12(i,1) = intm2(1)
      lstlst12(i,2) = intm2(2)

      intm1(1) = lstrstlR(i,1)
      intm1(2) = lstrstlR(i,2)

      intm2(1) = rthetat(1,1)*intm1(1) + rthetat(1,2)*intm1(2)
      intm2(2) = rthetat(2,1)*intm1(1) + rthetat(2,2)*intm1(2)

      lstRst12(i,1) = intm2(1)
      lstRst12(i,2) = intm2(2)

      intm1(1)  = 0.d0
      intm1(2)  = 0.d0 
      intm2(1) = 0.d0
      intm2(2) = 0.d0 
      
      enddo loopmixst


!--------------------------------------------------------------------------
C     Mixing CP-even Higgs 

      loopmixcpeven: do i =1,2
      
      intm1(1) = lstlst12(1,i)
      intm1(2) = lstlst12(2,i)

      call rmat2d(alpha,ralpha)
      
      intm2(1) = ralpha(1,1)*intm1(1) + ralpha(1,2)*intm1(2)
      intm2(2) = ralpha(2,1)*intm1(1) + ralpha(2,2)*intm1(2)

      lHstlst12(1,i) = intm2(1)
      lHstlst12(2,i) = intm2(2)

      intm1(1) = lstRst12(1,i)
      intm1(2) = lstRst12(2,i)

      intm2(1) = ralpha(1,1)*intm1(1) + ralpha(1,2)*intm1(2)
      intm2(2) = ralpha(2,1)*intm1(1) + ralpha(2,2)*intm1(2)

      lHstrst12(1, i) = intm2(1)
      lHstrst12(2, i) = intm2(2)
      
      ENDDO loopmixcpeven

!----------------------------------------------------------------------
C     Feynman rules for Charged Higgs
!----------------------------------------------------------------------
!     (Hpm Gpm, L R) basis

      lHcstlsblR(1,1) = (g*MW*sin2beta - (yuRG(3,3)**2.d0) * vev2 *
     $     cosbeta - (ydRG(3,3)**2.d0) * vev1 * sinbeta)/dsqrt(2.d0)
      
      lHcstlsblR(1,2) = (-sgnmu*modmu * ydRG(3,3) * cosbeta -
     $     ADRG(3,3) * ydRG(3,3) * sinbeta)

      lHcstlsblR(2,1) = (-g*MW*cos2beta - (yuRG(3,3)**2.d0) * vev2*
     $     sinbeta + (ydRG(3,3)**2.d0) * vev1 * cosbeta)/dsqrt(2.d0)

      lHcstlsblR(2,2) = ydRG(3,3) * (-sgnmu*modmu*sinbeta + 
     $     ADRG(3,3)*cosbeta)

!-------------------------------------------------------------------------

      intm1(1) = lHcstlsblR(1,1)
      intm1(2) = lHcstlsblR(1,2)

      call rmat2d(thetab,rthetab)
      
      intm2(1) = rthetab(1,1)*intm1(1) + rthetab(1,2)*intm1(2)
      intm2(2) = rthetab(2,1)*intm1(1) + rthetab(2,2)*intm1(2)

      
      lHcstlsb12(1, 1) = intm2(1)
      lHcstlsb12(1, 2) = intm2(2)

      intm1(1) = lHcstlsblR(2, 1)
      intm1(2) = lHcstlsblR(2, 2)


      call rmat2d(thetab,rthetab)
      
      intm2(1) = rthetab(1,1)*intm1(1) + rthetab(1,2)*intm1(2)
      intm2(2) = rthetab(2,1)*intm1(1) + rthetab(2,2)*intm1(2)

      lHcstlsb12(2, 1) = intm2(1)
      lHcstlsb12(2, 2) = intm2(2)

!---------------------------------------------------------------------------
!     (Hpm Gpm, L R) basis

      lHcstrsblr(1,1) = yuRG(3,3)*(-sgnmu*modmu*sinbeta - 
     $     AURG(3,3)*cosbeta)

      lHcstrsblr(1,2) = -1.d0 *(yuRG(3,3)*ydRG(3,3)*(vev1*cosbeta + 
     $     vev2*sinbeta)/dsqrt(2.d0))

      lHcstrsblr(2,1) = -yuRG(3,3)*(-sgnmu*modmu*cosbeta + 
     $     AURG(3,3)*sinbeta)  

      lHcstrsblr(2,2) = 0.d0

!-------------------------------------------------------------------------

      intm1(1) = lHcstrsblr(1,1)
      intm1(2) = lHcstrsblr(1,2)

      call rmat2d(thetab,rthetab)
      
      intm2(1) = rthetab(1,1)*intm1(1) + rthetab(1,2)*intm1(2)
      intm2(2) = rthetab(2,1)*intm1(1) + rthetab(2,2)*intm1(2)

      lHcstrsb12(1,1) = intm2(1)
      lHcstrsb12(1,2) = intm2(2)

      intm1(1) = lHcstrsblr(2,1)
      intm1(2) = lHcstrsblr(2,2)

      call rmat2d(thetab,rthetab)
      
      intm2(1) = rthetab(1,1)*intm1(1) + rthetab(1,2)*intm1(2)
      intm2(2) = rthetab(2,1)*intm1(1) + rthetab(2,2)*intm1(2)
      
      lHcstrsb12(2, 1) = intm2(1)
      lHcstrsb12(2, 2) = intm2(2)

!------------------------------------------------------------------------
C     Feynman rules - neutralino
!------------------------------------------------------------------------

      aPsi0str(1) = -4.d0*gp/(3.d0*dsqrt(2.d0))
      aPsi0str(2) = 0.d0
      aPsi0str(3) = 0.d0
      aPsi0str(4) = 0.d0

      bPsi0stl(1) = gp/(3.d0*dsqrt(2.d0))
      bPsi0stl(2) = g/dsqrt(2.d0)
      bPsi0stl(3) = 0.d0
      bPsi0stl(4) = 0.d0

      aPsi0stl(1) = 0.d0
      aPsi0stl(2) = 0.d0
      aPsi0stl(3) = 0.d0 
      aPsi0stl(4) = yuRG(3,3)

      bPsi0str(1) = 0.d0
      bPsi0str(2) = 0.d0
      bPsi0str(3) = 0.d0
      bPsi0str(4) = yuRG(3,3)


      aChi0stl(1) = ON(1,1)*aPsi0stl(1) + ON(1,2)*aPsi0stl(2) +
     $     ON(1,3)*aPsi0stl(1) + ON(1,4)*aPsi0stl(4) 
      aChi0stl(2) = ON(2,1)*aPsi0stl(1) + ON(2,2)*aPsi0stl(2)+
     $     ON(2,3)*aPsi0stl(1) + ON(2,4)*aPsi0stl(4) 
      aChi0stl(3) = ON(3,1)*aPsi0stl(1) + ON(3,2)*aPsi0stl(2)+
     $     ON(3,3)*aPsi0stl(1) + ON(3,4)*aPsi0stl(4) 
      aChi0stl(4) = ON(4,1)*aPsi0stl(1) + ON(4,2)*aPsi0stl(2)+
     $     ON(4,3)*aPsi0stl(1) + ON(4,4)*aPsi0stl(4) 



      bChi0stl(1) = ON(1,1)*bPsi0stl(1) + ON(1,2)*bPsi0stl(2)+
     $     ON(1,3)*bPsi0stl(1) + ON(1,4)*bPsi0stl(4) 
      bChi0stl(2) = ON(2,1)*bPsi0stl(1) + ON(2,2)*bPsi0stl(2)+
     $     ON(2,3)*bPsi0stl(1) + ON(2,4)*bPsi0stl(4) 
      bChi0stl(3) = ON(3,1)*bPsi0stl(1) + ON(3,2)*bPsi0stl(2)+
     $     ON(3,3)*bPsi0stl(1) + ON(3,4)*bPsi0stl(4) 
      bChi0stl(4) = ON(4,1)*bPsi0stl(1) + ON(4,2)*bPsi0stl(2)+
     $     ON(4,3)*bPsi0stl(1) + ON(4,4)*bPsi0stl(4) 



      aChi0str(1) = ON(1,1)*aPsi0str(1) + ON(1,2)*aPsi0str(2) +
     $     ON(1,3)*aPsi0str(3) + ON(1,4)*aPsi0str(4)
      aChi0str(2) = ON(2,1)*aPsi0str(1) + ON(2,2)*aPsi0str(2) +
     $     ON(2,3)*aPsi0str(3) + ON(2,4)*aPsi0str(4)
      aChi0str(3) = ON(3,1)*aPsi0str(1) + ON(3,2)*aPsi0str(2) +
     $     ON(3,3)*aPsi0str(3) + ON(3,4)*aPsi0str(4)
      aChi0str(4) = ON(4,1)*aPsi0str(1) + ON(4,2)*aPsi0str(2) +
     $     ON(4,3)*aPsi0str(3) + ON(4,4)*aPsi0str(4)



      bChi0str(1) = ON(1,1)*bPsi0str(1) + ON(1,2)*bPsi0str(2) +
     $     ON(1,3)*bPsi0str(3) + ON(1,4)*bPsi0str(4) 
      bChi0str(2) = ON(2,1)*bPsi0str(1) + ON(2,2)*bPsi0str(2) +
     $     ON(2,3)*bPsi0str(3) + ON(2,4)*bPsi0str(4) 
      bChi0str(3) = ON(3,1)*bPsi0str(1) + ON(3,2)*bPsi0str(2) +
     $     ON(3,3)*bPsi0str(3) + ON(3,4)*bPsi0str(4) 
      bChi0str(4) = ON(4,1)*bPsi0str(1) + ON(4,2)*bPsi0str(2) +
     $     ON(4,3)*bPsi0str(3) + ON(4,4)*bPsi0str(4) 


      loopchtllrr: DO i = 1, 4

      fChi0tstll(i) = (aChi0stl(i)*aChi0stl(i) + 
     $     bChi0stl(i)*bChi0stl(i))
      gChi0tstll(i) = (bChi0stl(i)*aChi0stl(i) + 
     $     bChi0stl(i)*aChi0stl(i))
      fChi0tstrr(i) = (aChi0str(i)*aChi0str(i) + 
     $     bChi0str(i)*bChi0str(i))
      gChi0tstrr(i) = (bChi0str(i)*aChi0str(i) + 
     $     bChi0str(i)*aChi0str(i))
      fChi0tstlr(i) = (aChi0str(i)*aChi0stl(i) + 
     $     bChi0str(i)*bChi0stl(i))
      gChi0tstlr(i) = (bChi0stl(i)*aChi0str(i) + 
     $     bChi0str(i)*aChi0stl(i))
      
      ENDDO loopchtllrr

!---------------------------------------------------------------------------
C     Feynman Rules for chargino
!---------------------------------------------------------------------------
      
      aPsicbstl(1) = g
      aPsicbstl(2) = 0.d0
      aPsicbstr(1) = 0.d0
      aPsicbstr(2) = -yuRG(3,3)
      bPsicbstl(1) = 0.d0
      bPsicbstl(2) = -ydRG(3,3)
      bPsicbstr(1) = 0.d0
      bPsicbstr(2) = 0.d0
      

c$$$      aChicbstl(1) = OCL(1,1)*aPsicbstl(1) + OCL(1,2)*aPsicbstl(2)
c$$$      aChicbstl(2) = OCL(2,1)*aPsicbstl(1) + OCL(2,2)*aPsicbstl(2)
c$$$
c$$$      bChicbstl(1) = OCR(1,1)*bPsicbstl(1) + OCR(1,2)*bPsicbstl(2)
c$$$      bChicbstl(2) = OCR(2,1)*bPsicbstl(1) + OCR(2,2)*bPsicbstl(2)
c$$$
c$$$      aChicbstr(1) = OCL(1,1)*aPsicbstr(1) + OCL(1,2)*aPsicbstr(2)
c$$$      aChicbstr(2) = OCL(2,1)*aPsicbstr(1) + OCL(2,2)*aPsicbstr(2)
c$$$
c$$$      bChicbstr(1) = OCR(1,1)*bPsicbstr(1) + OCR(1,2)*bPsicbstr(2)
c$$$      bChicbstr(2) = OCR(2,1)*bPsicbstr(1) + OCR(2,2)*bPsicbstr(2)

      aChicbstl(1) = OCR(1,1)*aPsicbstl(1) + OCR(1,2)*aPsicbstl(2)
      aChicbstl(2) = OCR(2,1)*aPsicbstl(1) + OCR(2,2)*aPsicbstl(2)

      bChicbstl(1) = OCL(1,1)*bPsicbstl(1) + OCL(1,2)*bPsicbstl(2)
      bChicbstl(2) = OCL(2,1)*bPsicbstl(1) + OCL(2,2)*bPsicbstl(2)

      aChicbstr(1) = OCR(1,1)*aPsicbstr(1) + OCR(1,2)*aPsicbstr(2)
      aChicbstr(2) = OCR(2,1)*aPsicbstr(1) + OCR(2,2)*aPsicbstr(2)

      bChicbstr(1) = OCL(1,1)*bPsicbstr(1) + OCL(1,2)*bPsicbstr(2)
      bChicbstr(2) = OCL(2,1)*bPsicbstr(1) + OCL(2,2)*bPsicbstr(2)
      
      
      loopchst: do i = 1, 2
      
      fChbstll(i) = (aChicbstl(i)*aChicbstl(i) +
     $     bChicbstl(i)*bChicbstl(i))
      gChbstll(i) = (bChicbstl(i)*aChicbstl(i) +
     $     aChicbstl(i)*bChicbstl(i))
      fChbstlr(i) = (aChicbstl(i)*aChicbstr(i) +
     $     bChicbstl(i)*bChicbstr(i))
      gChbstlr(i) = (bChicbstl(i)*aChicbstr(i) +
     $     aChicbstl(i)*bChicbstr(i))
      fChbstrr(i) = (aChicbstr(i)*aChicbstr(i) +
     $     bChicbstr(i)*bChicbstr(i))
      gChbstrr(i) = (bChicbstr(i)*aChicbstr(i) +
     $     aChicbstr(i)*bChicbstr(i))
      
      ENDDO loopchst

!-------------------------------------------
C     Corrections Begins From Here
!-------------------------------------------
      
      call funcg(p,M3t,mt,q,ggt)
      
!     print*,"ggt = ",ggt
      
      call f(p,mt1,0.d0,q,ft10)
      call f(p,mt2,0.d0,q,ft20)
      
c$$$  print*,"ft10 = ",ft10
c$$$  print* ,"q inside pistop = ", q
      

      call a0(muR,q,a0mu2)
      call a0(muL,q,a0mu1)
      call a0(mcR,q,a0mc2)
      call a0(mcL,q,a0mc1)
      call a0(mt2,q,a0t2)
      call a0(mt1,q,a0t1)

      call a0(msR,q,a0ms2)
      call a0(msL,q,a0ms1)
      call a0(mdR,q,a0md2)
      call a0(mdL,q,a0md1)
      call a0(mb2,q,a0mb2)
      call a0(mb1,q,a0mb1)

      call a0(meR,q,a0me2)
      call a0(meL,q,a0me1)
      call a0(mmuR,q,a0mmu2)
      call a0(mmuL,q,a0mmu1)
      call a0(mtau2,q,a0mtau2)
      call a0(mtau1,q,a0mtau1)

      call a0(snu(1),q,a0snu1)
      call a0(snu(2),q,a0snu2)
      call a0(snu(3),q,a0snu3)


      call a0(mA0,q,a0mA)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)
      call a0(mHpm,q,a0mHpm)
      call a0(MW,q,a0Mw)
      call a0(MZ,q,a0Mz)
      
      call f(p,mb1,MW,q,fmb1Mw)
      call f(p,mb2,MW,q,fmb2Mw)
      call f(p,mt1,MZ,q,fmt1Mz)
      call f(p,mt2,MZ,q,fmt2Mz)
      

      call b0(p,M3t,mt,q,b0M3mt)
      
!----------------------------------------------------------------
      
      gtterm(1, 1) = 4.d0 * ((g3**2.d0)/3.d0) *
     $     (((costhetat)**2.d0) * (ft10 + a0t1) + 
     $     ((sinthetat)**2.d0) * (ft20 + a0t2) + 2.d0 * ggt)
      
c$$$  print*,"gtterm(1,1) in pistop", gtterm(1,1)
c$$$  print*,"a0t2", a0t2
c$$$  print*,"a0t1", a0t1

      gtterm(2, 2) = 4.d0 * ((g3**2.d0)/3.d0) *
     $     (2.d0*ggt + ((sinthetat)**2.d0) *(ft10 + a0t1) + 
     $     ((costhetat)**2.d0) * (ft20 + a0t2))

      gtterm(1, 2) = 4.d0 * ((g3**2.d0)/3.d0) *
     $     (4.d0 * M3t*mt * b0M3mt + 
     $     (sinthetat*costhetat) * (ft10 - a0t1 - ft20 + a0t2))

!---------------------------------------------------------------------

      cstop(1, 1) = (yuRG(3,3)**2.d0) * ((sinthetat)**2.d0 * a0t1 + 
     $   (costhetat)**2.d0 * a0t2)

      cstop(2, 2) = (yuRG(3,3)**2.d0) * ((costhetat)**2.d0 * a0t1 + 
     $     (sinthetat)**2.d0 * a0t2)

      cstop(1, 2) = (yuRG(3,3)**2.d0) * costhetat*sinthetat * 3.d0 *
     $     (a0t1 - a0t2)
!---------------------------------------------------------------------

      csbtm(1, 1) = (ydRG(3,3)**2.d0) * ((sinthetab)**2.d0 * a0mb1 + 
     $     (costhetab)**2.d0 * a0mb2)

      csbtm(2, 2) = (yuRG(3,3)**2.d0) * ((costhetab)**2.d0 * a0mb1 + 
     $     (sinthetab)**2.d0 * a0mb2)
      
!---------------------------------------------------------------------
      
      higgsterm(1,1) = 0.5d0 * ((yuRG(3,3)**2.d0) * dnu(1) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((yuRG(3,3)**2.d0) * dnu(2) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((yuRG(3,3)**2.d0) * dnu(3) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((yuRG(3,3)**2.d0) * dnu(4) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(4)) * a0mA
      

      higgsterm(2,2) = 0.5d0 * ((yuRG(3,3)**2.d0) * dnu(1) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((yuRG(3,3)**2.d0) * dnu(2) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((yuRG(3,3)**2.d0) * dnu(3) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((yuRG(3,3)**2.d0) * dnu(4) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(4)) * a0mA


      higgsterm(1,1) = higgsterm(1,1) + 
     $     ((ydRG(3,3)**2.d0) * dnu(3) + ((g*g) * 
     $     ((guL*0.5d0/(costhw**2.d0)) - 0.5d0) * cn(3))) * a0mhpm +
     $     ((ydRG(3,3)**2.d0) * dnu(4) + ((g*g) * 
     $     ((guL*0.5d0/(costhw**2.d0)) - 0.5d0) * cn(4))) * a0MW


      higgsterm(2,2) = higgsterm(2,2) + 
     $     ((yuRG(3,3)**2.d0) * dnd(3) + 
     $     ((g*g) * (guR * 0.5d0/(costhw**2.d0)) * cn(3))) * a0mhpm +
     $     ((yuRG(3,3)**2.d0) * dnd(4) + 
     $     ((g*g) * (guR * 0.5d0/(costhw**2.d0)) * cn(4))) * a0MW

!----------------------------------------------------------------------------

      call b0(p,mHu0,mt1,q,b0h0mstop(1,1))
      call b0(p,mHu0,mt2,q,b0h0mstop(1,2))
      call b0(p,mh0,mt1,q,b0h0mstop(2,1))
      call b0(p,mh0,mt2,q,b0h0mstop(2,2))
      call b0(p,MZ,mt1,q,b0h0mstop(3,1))
      call b0(p,MZ,mt2,q,b0h0mstop(3,2))
      call b0(p,mA0,mt1,q,b0h0mstop(4,1))
      call b0(p,mA0,mt2,q,b0h0mstop(4,2))
      

      higgsterm(1,2) = 0.d0

      loophi: DO i = 1, 4        !<-------------neutral higgs terms
      loophj: DO j = 1, 2

      higgsterm(1,1) = higgsterm(1,1) + (lHstlst12(i,j)**2.d0) *
     $     b0h0mstop(i,j)

      higgsterm(1,2) = higgsterm(1,2) + (lHstlst12(i,j) *
     $     lHstrst12(i,j)) * b0h0mstop(i,j)

      higgsterm(2,2) = higgsterm(2,2) + (lHstrst12(i,j)**2.d0) *
     $     b0h0mstop(i,j)

      ENDDO loophj
      ENDDO loophi
!-----------------------------------------------------------------------

c$$$      call b0(p,mb1,MW,q,b0hcmbtm(2,1))
c$$$      call b0(p,mb2,MW,q,b0hcmbtm(2,2))
c$$$      call b0(p,mb1,mHpm,q,b0hcmbtm(1,1))
c$$$      call b0(p,mb2,mHpm,q,b0hcmbtm(1,2))
 
      call b0(p,MW,mb1,q,b0hcmbtm(2,1))
      call b0(p,MW,mb2,q,b0hcmbtm(2,2))
      call b0(p,mHpm,mb1,q,b0hcmbtm(1,1))
      call b0(p,mHpm,mb2,q,b0hcmbtm(1,2))
      
      
      loophci: DO i = 1, 2      !<-------------charged higgs terms
      loophcj: DO j = 1, 2

      higgsterm(1, 1) = higgsterm(1,1) + (lHcstlsb12(i,j)**2.d0) *
     $     b0hcmbtm(i,j)

      higgsterm(1, 2) = higgsterm(1,2) + 
     $     (lHcstlsb12(i, j)*lHcstrsb12(i, j)) * b0hcmbtm(i,j)

      higgsterm(2, 2) = higgsterm(2, 2) + (lHcstrsb12(i,j)**2.d0) *
     $     b0hcmbtm(i,j)
      
      ENDDO loophcj
      ENDDO loophci
      
!------------------------------------------------------------------------------
      
      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (4.d0*(g*g)/costhw**2.d0) * (guL*guL) * a0Mz  + 
     $     2.d0*(g*g) * a0Mw + 
     $     ((4.d0/9.d0)*g*g*sinsqthw)*((costhetat**2.d0) * ft10 + 
     $     (sinthetat**2.d0) * ft20) +
     $     (g*guL/costhw)**2.d0 * ((costhetat**2.d0) * fmt1MZ  +
     $     (sinthetat**2.d0) * fmt2MZ) +
     $     (g*g)*0.5d0*((costhetab**2.d0)*fmb1Mw + (sinthetab**2.d0)*
     $     fmb2Mw) + 
     $     (g*g)*0.25d0*((costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*
     $     a0t2 + 2.d0 * ((costhetab**2.d0)*a0mb1 + (sinthetab**2.d0)*
     $     a0mb2)) +
     $     (g*g)*0.5d0*(1.5d0*a0mu1 + !1.5d0*a0mu2 +
     $     1.5d0*a0mc1 +
     $     1.5d0*((costhetat**2.d0)*a0t1 + 
     $     (sinthetat**2.d0)*a0t2) -
     $     1.5d0*a0md1 - 1.5d0*a0ms1 - !1.5d0*a0md2 -
     $     1.5d0*((costhetab**2.d0)*a0mb1 +
     $     (sinthetab**2.d0)*a0mb2) +
     $     0.5d0*(a0snu1 + a0snu2 + a0snu3) -
     $     0.5d0*(a0me1 + a0mmu1 + !a0me2 +
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2)) +
     $     (gp**2.d0)*0.25d0*(yuL**2.d0)*((costhetat**2.d0)*a0t1 +
     $     (sinthetat**2.d0)*a0t2) +
     $     (gp**2.d0)*0.25d0*yuL*(3.d0*yuL*(a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 +(sinthetat**2.d0)*a0t2) +
     $     3.d0*yuR*(a0mu2 + a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 +(costhetat**2.d0)*a0t2) +
     $     3.d0*ydL*(a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0mb1 + (sinthetab**2.d0)*a0mb1) +
     $     3.d0*ydR*(a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0mb1 + (costhetab**2.d0)*a0mb2) +
     $     yeL*(a0me1 + a0mmu1 + 
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2) +
     $     yeR*(a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 +(costhetatau**2.d0)*a0mtau2) +
     $     ynuL*(a0snu1 + a0snu2 + a0snu3))
      
!-----------------------------------------------------------------------------

      higgsterm(2, 2) = higgsterm(2,2) +
     $     (4.d0*(g*g)/costhw**2.d0) * (guR**2.d0) * a0MZ + 
     $     ((4.d0/9.d0)*g*g*sinsqthw) * ((sinthetat)*ft10 + 
     $     (costhetat**2.d0)*ft20) +
     $     (g*guR/costhw)**2.d0 * ((sinthetat**2.d0)*fmt1MZ + 
     $     (costhetat**2.d0)*fmt2MZ) +
     $     (gp**2.d0) * 0.25d0 * (yuR**2.d0) *
     $     ((sinthetat**2.d0)*a0t1 + (costhetat**2.d0)*a0t2) +
     $     (gp**2.d0) * 0.25d0 * yuR * 
     $     (3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t2) +
     $     3.d0 * yuR * (a0mu2 +  a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 + (costhetat**2.d0)*a0t2) +
     $     3.d0 * ydL * (a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0mb1 + (sinthetab**2.d0)*a0mb2) +
     $     3.d0 * ydR * (a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0mb1 + (costhetab**2.d0)*a0mb2) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2) +
     $     yeR*(a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL*(a0snu1 + a0snu2 + a0snu3))

!----------------------------------------------------------------------------

      higgsterm(1,2) = higgsterm(1,2) + 
     $     (gp**2.d0) * 0.25d0 * 4.d0 * yuL*yuR * sinthetat*costhetat *
     $     (a0t1 - a0t2) +
     $     ((4.d0/9.d0)*g*g*sinsqthw) * sinthetat*costhetat * 
     $     (ft10 - ft20) -
     $     ((g*g)/costhw**2.d0) * guL*guR * sinthetat*costhetat *
     $     (fmt1MZ - fmt2MZ)

!----------------------------------------------------------------------
C     Chargino term
!----------------------------------------------------------------------

      call funcg(p,mchargino(1),mb,q,gmchargino1mb)
      call funcg(p,mchargino(2),mb,q,gmchargino2mb)
      call b0(p,mchargino(1),mb,q,b0mchargino1mb)
      call b0(p,mchargino(2),mb,q,b0mchargino2mb)
      
      chargino(1, 1) =  fChbstll(1)*gmchargino1mb -
     $     gChbstll(1)*mchargino(1)*mb*b0mchargino1mb*2.d0 +
     $     fChbstll(2)*gmchargino2mb -
     $     gChbstll(2)*mchargino(2)*mb*b0mchargino2mb*2.d0

      chargino(1, 2) =  fChbstlr(1)*gmchargino1mb -
     $    gChbstlr(1)*mchargino(1)*mb*b0mchargino1mb*2.d0 +
     $     fChbstlr(2)*gmchargino2mb -
     $    gChbstlr(2)*mchargino(2)*mb*b0mchargino2mb*2.d0

      chargino(2, 2) =  fChbstrr(1)*gmchargino1mb - 
     $     gChbstrr(1)*mchargino(1)*mb*b0mchargino1mb*2.d0 +
     $     fChbstrr(2)*gmchargino2mb - 
     $     gChbstrr(2)*mchargino(2)*mb*b0mchargino2mb*2.d0

!---------------------------------------------------------------------
C     neutralino terms
!---------------------------------------------------------------------

      call b0(p,nmneut(1),mt,q,b0mneut1mt)
      call b0(p,nmneut(2),mt,q,b0mneut2mt)
      call b0(p,nmneut(3),mt,q,b0mneut3mt)
      call b0(p,nmneut(4),mt,q,b0mneut4mt)
      
      call funcg(p,nmneut(1),mt,q,gmneut1mt)
      call funcg(p,nmneut(2),mt,q,gmneut2mt)
      call funcg(p,nmneut(3),mt,q,gmneut3mt)
      call funcg(p,nmneut(4),mt,q,gmneut4mt)

      neutralino(1, 1) = 
     $     fChi0tstll(1)*gmneut1mt - gChi0tstll(1)*2.d0*
     $     mneut(1)*mt*b0mneut1mt +
     $     fChi0tstll(2)*gmneut2mt - gChi0tstll(2)*2.d0*
     $     mneut(2)*mt*b0mneut2mt +
     $     fChi0tstll(3)*gmneut3mt - gChi0tstll(3)*2.d0*
     $     mneut(3)*mt*b0mneut3mt +
     $     fChi0tstll(4)*gmneut4mt - gChi0tstll(4)*2.d0*
     $     mneut(4)*mt*b0mneut4mt


      neutralino(2, 2) =
     $     fChi0tstrr(i)*gmneut1mt - gChi0tstrr(1)*2.d0*
     $     mneut(1)*mt*b0mneut1mt +
     $     fChi0tstrr(2)*gmneut2mt - gChi0tstrr(2)*2.d0*
     $     mneut(2)*mt*b0mneut2mt +
     $     fChi0tstrr(3)*gmneut3mt - gChi0tstrr(3)*2.d0*
     $     mneut(3)*mt*b0mneut3mt +
     $     fChi0tstrr(4)*gmneut4mt - gChi0tstrr(4)*2.d0*
     $     mneut(4)*mt*b0mneut4mt

      neutralino(1, 2) = 
     $     fChi0tstlr(1)*gmneut1mt - gChi0tstlr(1)*2.d0*
     $     mneut(1)*mt*b0mneut1mt +
     $     fChi0tstlr(2)*gmneut2mt - gChi0tstlr(2)*2.d0*
     $     mneut(2)*mt*b0mneut2mt +
     $     fChi0tstlr(3)*gmneut3mt - gChi0tstlr(3)*2.d0*
     $     mneut(3)*mt*b0mneut3mt +
     $     fChi0tstlr(4)*gmneut4mt - gChi0tstlr(4)*2.d0*
     $     mneut(4)*mt*b0mneut4mt 
      
!-------------------------------------------------------------------

      MSQU3(1,1) = mSQRG(3,3) + MT**2.d0 + guL*MZ*MZ*dcos(2.d0*beta)
      MSQU3(1,2) = MT * (AURG(3,3) - (sgnmu*modmu/dtan(beta)))
      MSQU3(2,1) = MSQU3(1,2)
      MSQU3(2,2) = mSURG(3,3) + MT**2.d0 + guR*MZ*MZ*dcos(2.d0*beta)

c$$$      print*,"mSQRG(3,3) = ", mSQRG(3,3)
c$$$      print*,"mSURG(3,3) = ", mSURG(3,3)
c$$$      print*,"AURG(3,3) = ", AURG(3,3)
c$$$      print*,"mt = ", mt
c$$$      print*,"mz = ", mz
c$$$      print*,"sinsqthw = ", sinsqthw
c$$$
c$$$      print*,"MSQU3(1,1) = ", MSQU3(1,1)
c$$$      print*,"MSQU3(1,2) = ", MSQU3(1,2)
c$$$      print*,"MSQU3(2,1) = ", MSQU3(2,1)
c$$$      print*,"MSQU3(2,2) = ", MSQU3(2,2)

!----------------------------------------------------------------------

      pist(1,1) = (cstop(1,1) + csbtm(1,1) + gtterm(1,1) +
     $     chargino(1,1) + neutralino(1,1) + higgsterm(1,1))/
     $     (16.d0*pi*pi)

      pist(1,2) = (cstop(1,2) + csbtm(1,2) + gtterm(1,2)  +
     $     chargino(1,2)+ neutralino(1,2) + higgsterm(1,2))/
     $     (16.d0*pi*pi)

      pist(2,2) = (cstop(2,2) + csbtm(2,2) + gtterm(2,2) +
     $     chargino(2,2)+ neutralino(2,2) +  higgsterm(2,2))/
     $     (16.d0*pi*pi)
      
      pist(2,1) = pist(1,2)

      delpist(1,1) = pist(1,1)
      delpist(1,2) = pist(1,2)
      delpist(2,1) = pist(2,1)
      delpist(2,2) = pist(2,2)


c$$$      print*,"cstop(1,1) = ", cstop(1,1)
c$$$      print*,"csbtm(1,1) = ", csbtm(1,1)
c$$$      print*,"gtterm(1,1) = ", gtterm(1,1)
c$$$      print*,"chargino(1,1) = ", chargino(1,1)
c$$$      print*,"neutralino(1,1) = ", neutralino(1,1)
c$$$      print*,"higgsterm(1,1) = ", higgsterm(1,1)
c$$$
c$$$      print*,"pisbtm(1,1) = ", pisbtm(1,1)
c$$$      print*,"pisbtm(1,2) = ", pisbtm(1,2)
c$$$      print*,"pisbtm(2,1) = ", pisbtm(2,1)
c$$$      print*,"pisbtm(2,2) = ", pisbtm(2,2)

      pist(1,1) = MSQU3(1,1) - pist(1,1)
      pist(1,2) = MSQU3(1,2) - pist(1,2)
      pist(2,2) = MSQU3(2,2) - pist(2,2)
      pist(2,1) = pist(1,2) 

!----------------------------------------------------------------------
C     find the singular values and the diagonalising matrices
!----------------------------------------------------------------------

      info  = 10
      AOK = 0
      
!     call dsyev('V','U',2,pist,2,STeg,work,lwork,info)

      Call CEigensystem(2,pist,2,STeg,stmix,2,0)
      
      if(info.eq.0) then
         AOK = AOK + 1
      endif
      
      
      RETURN

      END SUBROUTINE pistop

C=======================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C     pisbottom is checked on 22/05/2010 @ 14:30.
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      
      SUBROUTINE pisbottom(p,q,g,gp,g3,mt,mb,tanbeta,mSQRG,mSDRG,
     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3t,delpisbtm,SBeg)
      
      IMPLICIT NONE 

      INTEGER i,j,info,AOK,lwork
      parameter(lwork = 35)
      DOUBLE PRECISION work(lwork)
      DOUBLE PRECISION tanbeta,mSQRG(3,3),mSDRG(3,3)
      DOUBLE PRECISION AURG(3,3),ADRG(3,3)
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yuRG(3,3),ydRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      
      double precision mT, mB
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR,MSQD3(2,2)      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g,g3,M3t

      double precision beta
      double precision vev1,vev2
      
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae
      DOUBLE PRECISION pisbtm(2,2)
      
      DOUBLE PRECISION higgsterm(2,2)
      DOUBLE PRECISION dnu(4),dnd(4),cn(4)

      DOUBLE PRECISION lsblsblr(4,2),lsblsb12(4,2)
      DOUBLE PRECISION lsbrcblr(4,2),lsbrsb12(4,2)
      
      DOUBLE PRECISION aPsi0sbr(4), bPsi0sbr(4), aPsi0sbl(4)
      DOUBLE PRECISION bPsi0sbl(4) 
      DOUBLE PRECISION aChi0sbl(4),bChi0sbl(4), aChi0sbr(4)
      DOUBLE PRECISION bChi0sbr(4)

      DOUBLE PRECISION gChi0bsbll(4), fChi0bsbll(4)
      DOUBLE PRECISION gChi0bsblr(4), fChi0bsblr(4)
      DOUBLE PRECISION gChi0bsbrr(4), fChi0bsbrr(4)

      DOUBLE PRECISION bPsictsbl(2), bPsictsbr(2), aPsictsbl(2)
      DOUBLE PRECISION aPsictsbr(2)
      DOUBLE PRECISION aChictsbr(2), aChictsbl(2)
      DOUBLE PRECISION bChictsbl(2),bChictsbr(2)

      DOUBLE PRECISION gtterm(2,2), cstop(2,2), csbtm(2,2) 
      data gtterm/ 4 * 0.d0/, cstop/ 4 * 0.d0/, csbtm/ 4 * 0.d0/

      DOUBLE PRECISION chargino(2,2), neutralino(2,2)
      data chargino/ 4 * 0.d0/,neutralino/ 4 * 0.d0/

      DOUBLE PRECISION fChtsbLL(2), gChtsbLL(2) 
      DOUBLE PRECISION fChtsbLR(2), gChtsbLR(2) 
      DOUBLE PRECISION fChtsbRR(2), gChtsbRR(2)
      
      DOUBLE PRECISION intm1(2), intm2(2)

      DOUBLE PRECISION lHcsbrstlr(2,2),lHcsbrst12(2,2)
      DOUBLE PRECISION lHcsblst12(2,2)
      DOUBLE PRECISION lHcsblstlr(2,2) 

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION lHsblsb12(4,2),lHsbrsb12(4,2)
      data lHsblsb12/ 8 * 0.d0/, lHsbrsb12/ 8 * 0.d0/

      
      DOUBLE PRECISION rthetat(2,2),rthetab(2,2),ralpha(2,2)

      DOUBLE PRECISION ggb,fb10,fb20,a0b1,a0b2,b0M3mb,a0t1,a0t2
      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz
      DOUBLE PRECISION b0h0msbtm(4,2),b0hcmstop(2,2)

      DOUBLE PRECISION a0md1,a0md2,a0ms1,a0ms2,a0mtau1,a0mtau2,a0snu2
      DOUBLE PRECISION a0me1,a0me2,a0mmu1,a0mmu2,a0snu1,a0snu3
      DOUBLE PRECISION a0mc1,a0mc2,a0mu1,a0mu2
      
      DOUBLE PRECISION fmb1MZ,fmb2MZ,fmt1MW,fmt2MW
      DOUBLE PRECISION gmchargino1mt,gmchargino2mt
      DOUBLE PRECISION b0mchargino1mt,b0mchargino2mt
      DOUBLE PRECISION b0mneut1mb,b0mneut2mb,b0mneut3mb,b0mneut4mb
      DOUBLE PRECISION gmneut1mb,gmneut2mb,gmneut3mb,gmneut4mb
      
      DOUBLE PRECISION SBeg(2),nmneut(4),sbmix(2,2)

      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2
      DOUBLE PRECISION sinsqthw_susy,delpisbtm(2,2)
      
      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae
      common/sinsq_susy/sinsqthw_susy


      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d

      include 'stdinputs.h'

!     include 'stddef.dat'

C------------------------------------------------------------

      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      MZ = MZpole
      
      SBeg(1) = 0.d0
      SBeg(2) = 0.d0
      
      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw   = dsqrt(sinsqthw)
      thw = dasin(sinthw)
      costhw  = dcos(thw)
      cos2thw = dcos(2.d0*thw)


      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))
      
      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))
       
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
      
      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)

      mt1 = mtL 
      mt2 = mtR

      mb1 = mbL 
      mb2 = mbR

      mtau1 = mtauL 
      mtau2 = mtauR

!-----------------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta     = datan(tanbeta)
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha = 0.5d0 * datan(((mA0sq + MZ*MZ)/
     $     (mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta  = dcos(datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)

!-----------------------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0

!------------------------------------------------------------------------

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

!-------------------------------------------------------------------------
    
C-------------------------------------------------------------------------
!----change matrix name from stop-->bottom here.........................


      lsblsblr(1,1) = (g*MZ*gdL*cosbeta/costhw) + 
     $     (ydRG(3,3)**2.d0) * vev1
      lsblsblr(1,2) = (ydRG(3,3)*ADRG(3,3))/dsqrt(2.d0)

      lsblsblr(2,1) = ((- g*MZ*gdL*sinbeta)/costhw)
      lsblsblr(2,2) = (- ydRG(3,3)*sgnmu*modmu)/dsqrt(2.d0)

      lsblsblr(3,1) = 0.d0
      lsblsblr(3,2) = (-1.d0/dsqrt(2.d0)) *
     $     (-sgnmu*modmu*sinbeta*ydRG(3,3) + 
     $     ydRG(3,3)*ADRG(3,3)*cosbeta)

      lsblsblr(4,1) = 0.d0
      lsblsblr(4,2) = (-1.d0/dsqrt(2.d0)) * 
     $     (-sgnmu*modmu*cosbeta*ydRG(3,3) - 
     $     ydRG(3,3)*ADRG(3,3)*sinbeta)



      lsbrcblr(1,1) = lsblsblr(1,2)
      lsbrcblr(1,2) = (g*MZ*gdR*cosbeta/costhw) + 
     $     ydRG(3,3)**2.d0 * vev1

      lsbrcblr(2,1) = lsblsblr(2,2)
      lsbrcblr(2,2) = -g*MZ*gdR*sinbeta/costhw

      lsbrcblr(3,1) = -lsblsblr(3,2)
      lsbrcblr(3,2) = 0.d0

      lsbrcblr(4,1) = -lsblsblr(4,2)
      lsbrcblr(4,2) = 0.d0

!--------------------------------------------------------------------------------  
  
      loopmixst: do i = 1, 4

      intm1(1) = lsblsblr(i,1)
      intm1(2) = lsblsblr(i,2)

      call rmat2d(thetab,rthetab)
            
      intm2(1) = rthetab(1,1)*intm1(1) + rthetab(1,2)*intm1(2)
      intm2(2) = rthetab(2,1)*intm1(1) + rthetab(2,2)*intm1(2)

      lsblsb12(i,1) = intm2(1)
      lsblsb12(i,2) = intm2(2)

      intm1(1) = lsbrcblr(i,1)
      intm1(2) = lsbrcblr(i,2)

      intm2(1) = rthetab(1,1)*intm1(1) + rthetab(1,2)*intm1(2)
      intm2(2) = rthetab(2,1)*intm1(1) + rthetab(2,2)*intm1(2)

      lsbrsb12(i,1) = intm2(1)
      lsbrsb12(i,2) = intm2(2)

      enddo loopmixst

!--------------------------------------------------------------------------
C     Mixing CP-even Higgs 
!--------------------------------------------------------------------------

      loopmixcpeven: do i = 1, 2
      
      intm1(1) = lsblsb12(1,i)
      intm1(2) = lsblsb12(2,i)

      call rmat2d(alpha,ralpha)
            
      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lHsblsb12(1,i) = intm2(1)
      lHsblsb12(2,i) = intm2(2)

      intm1(1) = lsbrsb12(1,i)
      intm1(2) = lsbrsb12(2,i)

      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lHsbrsb12(1,i) = intm2(1)
      lHsbrsb12(2,i) = intm2(2)
     
      enddo loopmixcpeven

!----------------------------------------------------------------------
C    Feynman rules - charged higgs
!<--------- (H+ G+, L R) basis

      lHcsblstlr(1, 1) = ( g*MW*sin2beta - 
     $     (yuRG(3,3)**2.d0) * vev2 * cosbeta - 
     $     (ydRG(3,3)**2.d0) * vev1 * sinbeta)/dsqrt(2.d0)

      lHcsblstlr(1,2) = (-sgnmu*modmu*yuRG(3,3)*sinbeta -
     $     AURG(3,3)*yuRG(3,3)*cosbeta)


      lHcsblstlr(2,1) = (-g*MW*cos2beta - 
     $     yuRG(3,3)**2.d0 * vev2 * sinbeta + 
     $     ydRG(3,3)**2.d0 * vev1 * cosbeta)/dsqrt(2.d0)

      lHcsblstlr(2,2) = - yuRG(3,3) * (-sgnmu*modmu*cosbeta + 
     $     AURG(3,3)*sinbeta)

!-------------------------------------------------------------------------

      intm1(1) = lHcsblstlr(1,1)
      intm1(2) = lHcsblstlr(1,2)

      call rmat2d(thetat,rthetat)
            
      intm2(1) = rthetat(1,1)*intm1(1) + rthetat(1,2)*intm1(2)
      intm2(2) = rthetat(2,1)*intm1(1) + rthetat(2,2)*intm1(2)

 
      lHcsblst12(1, 1) = intm2(1)
      lHcsblst12(1, 2) = intm2(2)

      intm1(1) = lHcsblstlr(2, 1)
      intm1(2) = lHcsblstlr(2, 2)


      call rmat2d(thetat,rthetat)
            
      intm2(1) = rthetat(1,1)*intm1(1) + rthetat(1,2)*intm1(2)
      intm2(2) = rthetat(2,1)*intm1(1) + rthetat(2,2)*intm1(2)

      lHcsblst12(2, 1) = intm2(1)
      lHcsblst12(2, 2) = intm2(2)

!---------------------------------------------------------------------------
!<------------ (H+ G+, L R) basis

      lHcsbrstlr(1, 1) = ydRG(3,3)*(-sgnmu*modmu*cosbeta - 
     $     ADRG(3,3)*sinbeta)

      lHcsbrstlr(1, 2) = - ydRG(3,3)*yuRG(3,3) * (vev2*sinbeta + 
     $     vev1*cosbeta)/dsqrt(2.d0)


      lHcsbrstlr(2, 1) = ydRG(3,3)*(-sgnmu*modmu*sinbeta + 
     $     ADRG(3,3)*cosbeta)  

      lHcsbrstlr(2, 2) = 0.d0

!-------------------------------------------------------------------------

      intm1(1) = lHcsbrstlr(1,1)
      intm1(2) = lHcsbrstlr(1,2)

      call rmat2d(thetat,rthetat)
            
      intm2(1) = rthetat(1,1)*intm1(1) + rthetat(1,2)*intm1(2)
      intm2(2) = rthetat(2,1)*intm1(1) + rthetat(2,2)*intm1(2)

      lHcsbrst12(1,1) = intm2(1)
      lHcsbrst12(1,2) = intm2(2)

      intm1(1) = lHcsbrstlr(2,1)
      intm1(2) = lHcsbrstlr(2,2)

      call rmat2d(thetat,rthetat)
            
      intm2(1) = rthetat(1,1)*intm1(1) + rthetat(1,2)*intm1(2)
      intm2(2) = rthetat(2,1)*intm1(1) + rthetat(2,2)*intm1(2)
  
      lHcsbrst12(2, 1) = intm2(1)
      lHcsbrst12(2, 2) = intm2(2)

!------------------------------------------------------------------------
C     Feynman rules - neutralino
!------------------------------------------------------------------------

      aPsi0sbr(1) = (2.d0*gp/(3.d0*dsqrt(2.d0)))
      aPsi0sbr(2) = 0.d0
      aPsi0sbr(3) = 0.d0
      aPsi0sbr(4) = 0.d0

      bPsi0sbl(1) = gp/(3.d0*dsqrt(2.d0))
      bPsi0sbl(2) = - g/dsqrt(2.d0)
      bPsi0sbl(3) = 0.d0
      bPsi0sbl(4) = 0.d0

      aPsi0sbl(1) = 0.d0
      aPsi0sbl(2) = 0.d0
      aPsi0sbl(4) = 0.d0 
      aPsi0sbl(3) = ydRG(3,3)

      bPsi0sbr(1) = 0.d0
      bPsi0sbr(2) = 0.d0
      bPsi0sbr(4) = 0.d0
      bPsi0sbr(3) = ydRG(3,3)


      aChi0sbl(1) = ON(1,1)*aPsi0sbl(1) + ON(1,2)*aPsi0sbl(2) +
     $                 ON(1,3)*aPsi0sbl(1) + ON(1,4)*aPsi0sbl(4) 
      aChi0sbl(2) = ON(2,1)*aPsi0sbl(1) + ON(2,2)*aPsi0sbl(2) +
     $                 ON(2,3)*aPsi0sbl(1) + ON(2,4)*aPsi0sbl(4) 
      aChi0sbl(3) = ON(3,1)*aPsi0sbl(1) + ON(3,2)*aPsi0sbl(2) +
     $                 ON(3,3)*aPsi0sbl(1) + ON(3,4)*aPsi0sbl(4) 
      aChi0sbl(4) = ON(4,1)*aPsi0sbl(1) + ON(4,2)*aPsi0sbl(2) +
     $                 ON(4,3)*aPsi0sbl(1) + ON(4,4)*aPsi0sbl(4) 


      bChi0sbl(1) = ON(1,1)*bPsi0sbl(1) + ON(1,2)*bPsi0sbl(2) +
     $                 ON(1,3)*bPsi0sbl(1) + ON(1,4)*bPsi0sbl(4) 
      bChi0sbl(2) = ON(2,1)*bPsi0sbl(1) + ON(2,2)*bPsi0sbl(2) +
     $                 ON(2,3)*bPsi0sbl(1) + ON(2,4)*bPsi0sbl(4) 
      bChi0sbl(3) = ON(3,1)*bPsi0sbl(1) + ON(3,2)*bPsi0sbl(2) +
     $                 ON(3,3)*bPsi0sbl(1) + ON(3,4)*bPsi0sbl(4) 
      bChi0sbl(4) = ON(4,1)*bPsi0sbl(1) + ON(4,2)*bPsi0sbl(2) +
     $                 ON(4,3)*bPsi0sbl(1) + ON(4,4)*bPsi0sbl(4) 


      aChi0sbr(1) = ON(1,1)*aPsi0sbr(1) + ON(1,2)*aPsi0sbr(2) +
     $                 ON(1,3)*aPsi0sbr(3) + ON(1,4)*aPsi0sbr(4)
      aChi0sbr(2) = ON(2,1)*aPsi0sbr(1) + ON(2,2)*aPsi0sbr(2) +
     $                 ON(2,3)*aPsi0sbr(3) + ON(2,4)*aPsi0sbr(4)
      aChi0sbr(3) = ON(3,1)*aPsi0sbr(1) + ON(3,2)*aPsi0sbr(2) +
     $                 ON(3,3)*aPsi0sbr(3) + ON(3,4)*aPsi0sbr(4)
      aChi0sbr(4) = ON(4,1)*aPsi0sbr(1) + ON(4,2)*aPsi0sbr(2) +
     $                 ON(4,3)*aPsi0sbr(3) + ON(4,4)*aPsi0sbr(4)


      bChi0sbr(1) = ON(1,1)*bPsi0sbr(1) + ON(1,2)*bPsi0sbr(2) +
     $                 ON(1,3)*bPsi0sbr(3) + ON(1,4)*bPsi0sbr(4) 
      bChi0sbr(2) = ON(2,1)*bPsi0sbr(1) + ON(2,2)*bPsi0sbr(2) +
     $                 ON(2,3)*bPsi0sbr(3) + ON(2,4)*bPsi0sbr(4) 
      bChi0sbr(3) = ON(3,1)*bPsi0sbr(1) + ON(3,2)*bPsi0sbr(2) +
     $                 ON(3,3)*bPsi0sbr(3) + ON(3,4)*bPsi0sbr(4) 
      bChi0sbr(4) = ON(4,1)*bPsi0sbr(1) + ON(4,2)*bPsi0sbr(2) +
     $                 ON(4,3)*bPsi0sbr(3) + ON(4,4)*bPsi0sbr(4) 


      loopchtllrr: do i = 1, 4

      fChi0bsbll(i) = (aChi0sbl(i)*aChi0sbl(i) + 
     $    bChi0sbl(i)*bChi0sbl(i))
      gChi0bsbll(i) = (bChi0sbl(i)*aChi0sbl(i) + 
     $    bChi0sbl(i)*aChi0sbl(i))
      fChi0bsbrr(i) = (aChi0sbr(i)*aChi0sbr(i) + 
     $     bChi0sbr(i)*bChi0sbr(i))
      gChi0bsbrr(i) = (bChi0sbr(i)*aChi0sbr(i) + 
     $     bChi0sbr(i)*aChi0sbr(i))
      fChi0bsblr(i) = (aChi0sbr(i)*aChi0sbl(i) + 
     $     bChi0sbr(i)*bChi0sbl(i))
      gChi0bsblr(i) = (bChi0sbl(i)*aChi0sbr(i) + 
     $     bChi0sbr(i)*aChi0sbl(i))
      
      ENDDO loopchtllrr

!---------------------------------------------------------------------------
C     Feynman Rules for chargino
!---------------------------------------------------------------------------
      
      aPsictsbl(1) = 0.d0
      aPsictsbl(2) = - yuRG(3,3)

      aPsictsbr(1) = 0.d0
      aPsictsbr(2) = 0.d0

      bPsictsbl(1) = g
      bPsictsbl(2) = 0.d0

      bPsictsbr(1) = 0.d0
      bPsictsbr(2) = - ydRG(3,3)
      

      aChictsbl(1) = OCR(1,1)*aPsictsbl(1) + OCR(1,2)*aPsictsbl(2)
      aChictsbl(2) = OCR(2,1)*aPsictsbl(1) + OCR(2,2)*aPsictsbl(2)

      bChictsbl(1) = OCL(1,1)*bPsictsbl(1) + OCL(1,2)*bPsictsbl(2)
      bChictsbl(2) = OCL(2,1)*bPsictsbl(1) + OCL(2,2)*bPsictsbl(2)

      aChictsbr(1) = OCR(1,1)*aPsictsbr(1) + OCR(1,2)*aPsictsbr(2)
      aChictsbr(2) = OCR(2,1)*aPsictsbr(1) + OCR(2,2)*aPsictsbr(2)

      bChictsbr(1) = OCL(1,1)*bPsictsbr(1) + OCL(1,2)*bPsictsbr(2)
      bChictsbr(2) = OCL(2,1)*bPsictsbr(1) + OCL(2,2)*bPsictsbr(2)

      
      loopchst: do i = 1, 2

      fChtsbLL(i) = (aChictsbl(i)*aChictsbl(i) +
     $     bChictsbl(i)*bChictsbl(i))
      gChtsbLL(i) = (bChictsbl(i)*aChictsbl(i) +
     $     aChictsbl(i)*bChictsbl(i))
      fChtsbLR(i) = (aChictsbl(i)*aChictsbr(i) +
     $     bChictsbl(i)*bChictsbr(i))
      gChtsbLR(i) = (bChictsbl(i)*aChictsbr(i) +
     $     aChictsbl(i)*bChictsbr(i))
      fChtsbRR(i) = (aChictsbr(i)*aChictsbr(i) +
     $     bChictsbr(i)*bChictsbr(i))
      gChtsbRR(i) = (bChictsbr(i)*aChictsbr(i) +
     $     aChictsbr(i)*bChictsbr(i))
      
      enddo loopchst

!--------------------------------------------------------------------
C     Corrections Begin
!--------------------------------------------------------------------

      call funcg(p,M3t,mb,q,ggb)

!     print*,"ggt"
      
      call f(p,mb1,0.d0,q,fb10)
      call f(p,mb2,0.d0,q,fb20)

      call a0(muR,q,a0mu2)
      call a0(muL,q,a0mu1)
      call a0(mcR,q,a0mc2)
      call a0(mcL,q,a0mc1)
      call a0(mt2,q,a0t2)
      call a0(mt1,q,a0t1)

      call a0(msR,q,a0ms2)
      call a0(msL,q,a0ms1)
      call a0(mdR,q,a0md2)
      call a0(mdL,q,a0md1)
      call a0(mb2,q,a0b2)
      call a0(mb1,q,a0b1)

      call a0(meR,q,a0me2)
      call a0(meL,q,a0me1)
      call a0(mmuR,q,a0mmu2)
      call a0(mmuL,q,a0mmu1)
      call a0(mtau2,q,a0mtau2)
      call a0(mtau1,q,a0mtau1)

      call a0(snu(1),q,a0snu1)
      call a0(snu(2),q,a0snu2)
      call a0(snu(3),q,a0snu3)


      call a0(mA0,q,a0mA)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)
      call a0(mHpm,q,a0mHpm)
      call a0(MW,q,a0Mw)
      call a0(MZ,q,a0Mz)

      call f(p,mt1,MW,q,fmt1Mw)
      call f(p,mt2,MW,q,fmt2Mw)
      call f(p,mb1,MZ,q,fmb1Mz)
      call f(p,mb2,MZ,q,fmb2Mz)


      call b0(p,M3t,mb,q,b0M3mb)
      
!--------------------------------------------------------------------
      
      gtterm(1,1) = (4.d0 * (g3**2.d0)/3.d0) *
     $     (2.d0*ggb + (costhetab)**2.d0 * (fb10 + a0b1) + 
     $     (sinthetab)**2.d0 * (fb20 + a0b2))

      gtterm(2,2) = (4.d0 * (g3**2.d0)/3.d0) *
     $     (2.d0*ggb + (sinthetab)**2.d0 * (fb10 + a0b1) + 
     $     (costhetab)**2.d0 * (fb20 + a0b2))

      gtterm(1,2) = (4.d0*(g3**2.d0)/3.d0) *
     $     (4.d0* M3t*mb * b0M3mb + 
     $     sinthetab*costhetab * (fb10 - a0b1 - fb20 + a0b2))

!--------------------------------------------------------------------

      csbtm(1, 1) = (ydRG(3,3)**2.d0) * ((sinthetab)**2.d0*a0b1 + 
     $     (costhetab)**2.d0*a0b2)

      csbtm(2, 2) = (ydRG(3,3)**2.d0) * ((costhetab)**2.d0*a0b1 + 
     $     (sinthetab)**2.d0*a0b2)

      csbtm(1, 2) = (ydRG(3,3)**2.d0) * costhetab*sinthetab * 3.d0 *
     $     (a0b1 - a0b2)

!--------------------------------------------------------------------

      cstop(1, 1) = (yuRG(3,3)**2.d0) * ((sinthetat)**2.d0 * a0t1 + 
     $     (costhetat)**2.d0 * a0t2)

      cstop(2, 2) = (ydRG(3,3)**2.d0) * ((costhetat)**2.d0 * a0t1 + 
     $     (sinthetat)**2.d0 * a0t2)

!--------------------------------------------------------------------

      higgsterm(1,1) = 0.5d0 * ((ydRG(3,3)**2.d0) * dnd(1) -
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((ydRG(3,3)**2.d0) * dnd(2) - 
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((ydRG(3,3)**2.d0) * dnd(3) - 
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((ydRG(3,3)**2.d0) * dnd(4) - 
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA

      higgsterm(1,2) = 0.d0

      higgsterm(2,2) = 0.5d0 * ((ydRG(3,3)**2.d0) * dnd(1) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((ydRG(3,3)**2.d0) * dnd(2) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((ydRG(3,3)**2.d0) * dnd(3) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((ydRG(3,3)**2.d0) * dnd(4) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA


      higgsterm(1,1) = higgsterm(1,1) + 
     $     ((yuRG(3,3)**2.d0) * dnd(3) + ((g*g) * 
     $     ((gdL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(3))) * a0Mhpm +
     $     ((yuRG(3,3)**2.d0) * dnd(4) + ((g*g) * 
     $     ((gdL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(4))) * a0MW

      higgsterm(2,2) = higgsterm(2,2) + 
     $     ((ydRG(3,3)**2.d0) * dnu(3) +  
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(3)) * a0Mhpm +
     $     ((ydRG(3,3)**2.d0) * dnu(4) + 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(4)) * a0MW  !<----------------chkd

!---------------------------------------------------------------------------

      call b0(p,mHu0,mb1,q,b0h0msbtm(1,1))
      call b0(p,mHu0,mb2,q,b0h0msbtm(1,2))
      call b0(p,mh0,mb1,q,b0h0msbtm(2,1))
      call b0(p,mh0,mb2,q,b0h0msbtm(2,2))
      call b0(p,MZ,mb1,q,b0h0msbtm(3,1))
      call b0(p,MZ,mb2,q,b0h0msbtm(3,2))
      call b0(p,mA0,mb1,q,b0h0msbtm(4,1))
      call b0(p,mA0,mb2,q,b0h0msbtm(4,2))
       

      loophi: DO i = 1, 4       !<-------------neutral higgs terms
      loophj: DO j = 1, 2
      
      higgsterm(1,1) = higgsterm(1,1) + 
     $     (lHsblsb12(i,j)**2.d0) * b0h0msbtm(i,j)
      
      higgsterm(1,2) = higgsterm(1,2) + 
     $     lHsblsb12(i,j) * lHsbrsb12(i,j) * b0h0msbtm(i,j)
      
      higgsterm(2,2) = higgsterm(2,2) + 
     $     (lHsbrsb12(i,j)**2.d0) * b0h0msbtm(i,j)
      
      ENDDO loophj
      ENDDO loophi

!----------------------------------------------------------------------------

      call b0(p,mt1,mHpm,q,b0hcmstop(1,1))
      call b0(p,mt2,mHpm,q,b0hcmstop(1,2))
      call b0(p,mt1,MW,q,b0hcmstop(2,1))
      call b0(p,mt2,MW,q,b0hcmstop(2,2))


      loophci: do i = 1, 2      !<-------------charged higgs terms
      loophcj: do j = 1, 2

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (lHcsblst12(i,j)**2.d0) * b0hcmstop(i,j)

      higgsterm(1, 2) = higgsterm(1,2) + 
     $     lHcsblst12(i, j) *lHcsbrst12(i, j) *b0hcmstop(i,j)

      higgsterm(2, 2) = higgsterm(2, 2) + 
     $     (lHcsbrst12(i,j)**2.d0) * b0hcmstop(i,j)
      
      enddo loophcj
      enddo loophci

!------------------------------------------------------------------------------

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (4.d0*(g*g)/costhw**2.d0)*(gdL*gdL)*a0MZ  + 
     $     2.d0*(g*g)*a0MW + 
     $     ((1.d0/9.d0)*g*g*sinsqthw) * ((costhetab**2.d0) * 
     $     fb10 + (sinthetab**2.d0) * fb20) +
     $     (g*gdL/costhw)**2.d0 * ((costhetab**2.d0) * fmb1Mz +
     $     (sinthetab**2.d0) * fmb2mz) +
     $     (g*g)*0.5d0*((costhetat)**2.d0 * fmt1Mw + (sinthetat**2.d0) *
     $     fmt2Mw) +
     $     (g*g) * 0.25d0 * (1.d0 * ((costhetab**2.d0) * a0b1 + 
     $     (sinthetab**2.d0) * a0b2) + 2.d0* ((costhetat**2.d0) * a0t1 + 
     $     (sinthetat**2.d0) * a0t2)) -
     $     (g*g) * 0.5d0 * (- 1.5d0 * a0md1 - 1.5d0*a0ms1 +
     $     1.5d0 * a0mu1 + 1.5d0*a0mc1 +
     $     0.5d0 * (a0snu1 + a0snu2 + a0snu3) -
     $     0.5d0 * (a0me1 + a0mmu1 + ((costhetatau**2.d0) * a0mtau1 +
     $     (sinthetatau**2.d0) * a0mtau2))) +
     $     (gp**2.d0) * 0.25d0 *(ydL**2.d0) * ((costhetab**2.d0) * a0b1+
     $     (sinthetab**2.d0)*a0b2) +
     $     (gp**2.d0) * 0.25d0 * ydL * ((3.d0*ydL*(a0md1 + a0ms1 + 
     $     ((costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2))) +
     $     (3.d0 * ydR *(a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2)) +
     $     (3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t1)) +
     $     (3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (costhetat**2.d0)*a0t2 + (sinthetat**2.d0)*a0t1)) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))
     
!----------------------------------------------------------------------------------

      higgsterm(2, 2) = higgsterm(2,2) +
     $     (4.d0*(g*g)/costhw**2.d0) * (gdR**2.d0) * a0MZ + 
     $     ((1.d0/9.d0)*g*g*sinsqthw) * ((sinthetab)*fb10 + 
     $     (costhetab**2.d0)*fb20) +
     $     (g*gdR/costhw)**2.d0 * ((sinthetab**2.d0) * fmb1MZ + 
     $     (costhetab**2.d0) * fmb2MZ) +
     $     (gp**2.d0) * 0.25d0 * (ydR**2.d0) * 
     $     ((sinthetab**2.d0) * a0b1 + (costhetab**2.d0) * a0b2) +
     $     (gp**2.d0) * 0.25d0 * ydR *  
     $     (3.d0 * ydL * (a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2) +
     $     3.d0 * ydR * (a0md2 +  a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2) +
     $     3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t2) +
     $     3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 + (costhetat**2.d0)*a0t2) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))

!----------------------------------------------------------------------

      higgsterm(1,2) = higgsterm(1,2) + 
     $     (gp**2.d0) * ydL*ydR * sinthetab*costhetab *
     $     (a0b1 - a0b2) +
     $     ((1.d0/9.d0) * g*g * sinsqthw) * sinthetab*costhetab * 
     $     (fb10 - fb20) -
     $     ((g*g)/costhw**2.d0) * gdL*gdR * sinthetab*costhetab *
     $     (fmb1Mz - fmb2Mz)    !<-----------------------chked

!----------------------------------------------------------------------
C       Chargino term
!----------------------------------------------------------------------

      call funcg(p,mchargino(1),mt,q,gmchargino1mt)
      call funcg(p,mchargino(2),mt,q,gmchargino2mt)
      call b0(p,mchargino(1),mt,q,b0mchargino1mt)
      call b0(p,mchargino(2),mt,q,b0mchargino2mt)
      
      chargino(1, 1) =  fChtsbLL(1)*gmchargino1mt -
     $     gChtsbLL(1)*mchargino(1)*mt*b0mchargino1mt*2.d0 +
     $     fChtsbLL(2)*gmchargino2mt -
     $     gChtsbLL(2)*mchargino(2)*mt*b0mchargino2mt*2.d0
      
      chargino(1, 2) =  fChtsbLR(1)*gmchargino1mt -
     $     gChtsbLR(1)*mchargino(1)*mt*b0mchargino1mt*2.d0 +
     $     fChtsbLR(2)*gmchargino2mt -
     $     gChtsbLR(2)*mchargino(2)*mt*b0mchargino2mt*2.d0

      chargino(2, 2) = fChtsbRR(1)*gmchargino1mt - 
     $     gChtsbRR(1)*mchargino(1)*mt*b0mchargino1mt*2.d0 +
     $     fChtsbRR(2)*gmchargino2mt - 
     $     gChtsbRR(2)*mchargino(2)*mt*b0mchargino2mt*2.d0
      
!---------------------------------------------------------------------
C     neutralino terms
!---------------------------------------------------------------------
      
      call b0(p,nmneut(1),mb,q,b0mneut1mb)
      call b0(p,nmneut(2),mb,q,b0mneut2mb)
      call b0(p,nmneut(3),mb,q,b0mneut3mb)
      call b0(p,nmneut(4),mb,q,b0mneut4mb)
      
      call funcg(p,nmneut(1),mb,q,gmneut1mb)
      call funcg(p,nmneut(2),mb,q,gmneut2mb)
      call funcg(p,nmneut(3),mb,q,gmneut3mb)
      call funcg(p,nmneut(4),mb,q,gmneut4mb)
      
      neutralino(1, 1) = 
     $     fChi0bsbll(1)*gmneut1mb - gChi0bsbll(1)*2.d0*
     $     mneut(1)*mb*b0mneut1mb +
     $     fChi0bsbll(2)*gmneut2mb - gChi0bsbll(2)*2.d0*
     $     mneut(2)*mb*b0mneut2mb +
     $     fChi0bsbll(3)*gmneut3mb - gChi0bsbll(3)*2.d0*
     $     mneut(3)*mb*b0mneut3mb +
     $     fChi0bsbll(4)*gmneut4mb - gChi0bsbll(4)*2.d0*
     $     mneut(4)*mb*b0mneut4mb
      
      
      neutralino(2, 2) = 
     $     fChi0bsbrr(i)*gmneut1mb - gChi0bsbrr(1)*2.d0*
     $     mneut(1)*mb*b0mneut1mb +
     $     fChi0bsbrr(2)*gmneut2mb - gChi0bsbrr(2)*2.d0*
     $     mneut(2)*mb*b0mneut2mb +
     $     fChi0bsbrr(3)*gmneut3mb - gChi0bsbrr(3)*2.d0*
     $     mneut(3)*mb*b0mneut3mb +
     $     fChi0bsbrr(4)*gmneut4mb - gChi0bsbrr(4)*2.d0*
     $     mneut(4)*mb*b0mneut4mb
      
      neutralino(1, 2) = 
     $     fChi0bsblr(1)*gmneut1mb - gChi0bsblr(1)*2.d0*
     $     mneut(1)*mb*b0mneut1mb +
     $     fChi0bsblr(2)*gmneut2mb - gChi0bsblr(2)*2.d0*
     $     mneut(2)*mb*b0mneut2mb +
     $     fChi0bsblr(3)*gmneut3mb - gChi0bsblr(3)*2.d0*
     $     mneut(3)*mb*b0mneut3mb +
     $     fChi0bsblr(4)*gmneut4mb - gChi0bsblr(4)*2.d0*
     $     mneut(4)*mb*b0mneut4mb 
      
!-------------------------------------------------------------------
      
      MSQD3(1,1) = mSQRG(3,3) + MB**2.d0 + gdL*MZ*MZ*dcos(2.d0*beta)
      MSQD3(1,2) = MB*((ADRG(3,3) - sgnmu*modmu*dtan(beta)))
      MSQD3(2,1) = MSQD3(1,2)
      MSQD3(2,2) = mSDRG(3,3) + MB**2.d0 + gdR*MZ*MZ*dcos(2.d0*beta)

c$$$      print*,"MSQD3(1,1) = ", MSQD3(1,1)
c$$$      print*,"MSQD3(1,2) = ", MSQD3(1,2)
c$$$      print*,"MSQD3(2,1) = ", MSQD3(2,1)
c$$$      print*,"MSQD3(2,2) = ", MSQD3(2,2)
      
!----------------------------------------------------------------------
      
      pisbtm(1,1) = (cstop(1,1) + csbtm(1,1) + gtterm(1,1) +
     $     chargino(1,1) + neutralino(1,1) + higgsterm(1,1))/
     $     (16.d0*pi*pi)
      
      pisbtm(1,2) = (cstop(1,2) + csbtm(1,2) + gtterm(1,2)  +
     $     chargino(1,2)+ neutralino(1,2) + higgsterm(1,2))/
     $     (16.d0*pi*pi)
      
      pisbtm(2,2) = (cstop(2,2) + csbtm(2,2) + gtterm(2,2) +
     $     chargino(2,2)+ neutralino(2,2) +  higgsterm(2,2))/
     $     (16.d0*pi*pi)
      
      pisbtm(2,1) = pisbtm(1,2)
      
      delpisbtm(1,1) = pisbtm(1,1)
      delpisbtm(1,2) = pisbtm(1,2)
      delpisbtm(2,1) = pisbtm(2,1)
      delpisbtm(2,2) = pisbtm(2,2)

c$$$      print*,"MSQD3(1,1) = ", MSQD3(1,1)
c$$$      print*,"MSQD3(1,2) = ", MSQD3(1,2)
c$$$      print*,"MSQD3(2,1) = ", MSQD3(2,1)
c$$$      print*,"MSQD3(2,2) = ", MSQD3(2,2)

c$$$      print*,"cstop(1,1) = ", cstop(1,1)
c$$$      print*,"csbtm(1,1) = ", csbtm(1,1)
c$$$      print*,"gtterm(1,1) = ", gtterm(1,1)
c$$$      print*,"chargino(1,1) = ", chargino(1,1)
c$$$      print*,"neutralino(1,1) = ", neutralino(1,1)
c$$$      print*,"higgsterm(1,1) = ", higgsterm(1,1)
c$$$
c$$$      print*,"pisbtm(1,1) = ", pisbtm(1,1)
c$$$      print*,"pisbtm(1,2) = ", pisbtm(1,2)
c$$$      print*,"pisbtm(2,1) = ", pisbtm(2,1)
c$$$      print*,"pisbtm(2,2) = ", pisbtm(2,2)


      pisbtm(1,1) = MSQD3(1,1) - pisbtm(1,1)
      pisbtm(1,2) = MSQD3(1,2) - pisbtm(1,2)
      pisbtm(2,2) = MSQD3(2,2) - pisbtm(2,2)
      pisbtm(2,1) = pisbtm(1,2) 
      
!----------------------------------------------------------------------
C     find the singular values and the diagonalising matrices
!----------------------------------------------------------------------
      
      info  = 10
      AOK = 0
      
C     call dsyev('V','U',2,pisbtm,2,SBeg,work,lwork,info)
      
      Call CEigensystem(2,pisbtm,2,SBeg,sbmix,2,0)
      if(info.eq.0) then
         AOK = AOK + 1
      endif
      
      
      RETURN
      
      END SUBROUTINE pisbottom
      
C====================================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C     pisdown is checked on 24/05/2010 @ 18:30.
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 
      SUBROUTINE pisdown(p,q,g,gp,g3,tanbeta,mSQRG,mSDRG,
     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3t,SDeg)


      IMPLICIT NONE 

      INTEGER i,j,info,AOK,lwork
      parameter(lwork = 35)
      DOUBLE PRECISION work(lwork)
      DOUBLE PRECISION tanbeta,mSQRG(3,3),mSDRG(3,3)
      DOUBLE PRECISION AURG(3,3),ADRG(3,3)
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yuRG(3,3),ydRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR,MSQD1(2,2)      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g,g3,M3t
   
      double precision vev1,vev2,beta
 
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetac,thetas,thetamu

      DOUBLE PRECISION thetau,thetad,thetae,costhetau
      DOUBLE PRECISION costhetad,sinthetad,sinthetau

      DOUBLE PRECISION pisdn(2,2)
 
      DOUBLE PRECISION higgsterm(2,2)
      DOUBLE PRECISION dnu(4),dnd(4),cn(4)

      DOUBLE PRECISION lsdnLsdnLR(4,2),lsdnLsdn12(4,2)
      DOUBLE PRECISION lsdnRsdnLR(4,2),lsdnRsdn12(4,2)
      
      DOUBLE PRECISION aPsi0sdnr(4), bPsi0sdnr(4), aPsi0sdnl(4)
      DOUBLE PRECISION bPsi0sdnl(4) 
      DOUBLE PRECISION aChi0sdnl(4), bChi0sdnl(4), aChi0sdnr(4)
      DOUBLE PRECISION bChi0sdnr(4)

      DOUBLE PRECISION gChi0dsdnLL(4), fChi0dsdnLL(4)
      DOUBLE PRECISION gChi0dsdnLR(4), fChi0dsdnLR(4)
      DOUBLE PRECISION gChi0dsdnRR(4), fChi0dsdnRR(4)

      DOUBLE PRECISION bPsicsdnl(2), bPsicsdnr(2), aPsicsdnl(2)
      DOUBLE PRECISION aPsicsdnr(2)
      DOUBLE PRECISION aChicsdnr(2), aChicsdnl(2)
      DOUBLE PRECISION bChicsdnl(2),bChicsdnr(2)


      DOUBLE PRECISION gtterm(2,2), cstop(2,2), csbtm(2,2) 
      data gtterm/ 4 * 0.d0/, cstop/ 4 * 0.d0/, csbtm/ 4 * 0.d0/

      DOUBLE PRECISION  chargino(2,2), neutralino(2,2)
      data chargino/ 4 * 0.d0/,neutralino/ 4 * 0.d0/

      DOUBLE PRECISION fChsdnLL(2), gChsdnLL(2) 
      DOUBLE PRECISION fChsdnLR(2), gChsdnLR(2) 
      DOUBLE PRECISION fChsdnRR(2), gChsdnRR(2)
 
      DOUBLE PRECISION intm1(2), intm2(2)

      DOUBLE PRECISION lHcsuprsdnlr(2,2),lHcsuprsdn12(2,2)
      DOUBLE PRECISION lHcsuplsdn12(2,2)
      DOUBLE PRECISION lHcsuplsdnlr(2,2) 

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION lHsdnLsdn12(4,2),lHsdnRsdn12(4,2)
      data lHsdnLsdn12/ 8 * 0.d0/, lHsdnRsdn12/ 8 * 0.d0/
      
      DOUBLE PRECISION rthetau(2,2),rthetad(2,2),ralpha(2,2)

      DOUBLE PRECISION ggd,fd10,fd20,a0d1,a0d2,b0M3md
      DOUBLE PRECISION a0b1,a0b2,a0t1,a0t2
      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz
      DOUBLE PRECISION b0h0msdn(4,2),b0hcmsu(2,2)

      DOUBLE PRECISION a0md1,a0md2,a0ms1,a0ms2,a0mtau1,a0mtau2,a0snu2
      DOUBLE PRECISION a0me1,a0me2,a0mmu1,a0mmu2,a0snu1,a0snu3
      DOUBLE PRECISION a0mc1,a0mc2,a0mu1,a0mu2,nmneut(4)
      
      DOUBLE PRECISION fmd1MZ,fmd2MZ,fmu1MW,fmu2MW
      DOUBLE PRECISION gmchargino1mu,gmchargino2mu
      DOUBLE PRECISION b0mchargino1mu,b0mchargino2mu
      DOUBLE PRECISION b0mneut1md,b0mneut2md,b0mneut3md,b0mneut4md
      DOUBLE PRECISION gmneut1md,gmneut2md,gmneut3md,gmneut4md
      
      DOUBLE PRECISION SDeg(2),sdmix(2,2)
      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2
      DOUBLE PRECISION mu1,mu2,md1,md2,mu,sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae
  
      common/sinsq_susy/sinsqthw_susy

      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d
  

      include 'stdinputs.h'

!-------------------------------------------

      pi = 4.d0 * datan(1.d0)
      MZ = MZpole
      MW = MWpole

      mu = MUq

      SDeg(1) = 0.d0
      SDeg(2) = 0.d0


      sinsqthw = sinsqthw_susy  !1.d0 - (MW/MZ)**2.d0
      sinthw = dsqrt(sinsqthw)
      thw = dasin(sinthw)
      costhw = dcos(thw)
      cos2thw = dcos(2.d0*thw)

      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))
      
      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))


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

      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)

      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)

      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)

      costhetau = dcos(thetau)
      sinthetau = dsin(thetau)

      costhetad = dcos(thetad)
      sinthetad = dsin(thetad)


      mt1 = mtL 
      mt2 = mtR

      mb1 = mbL 
      mb2 = mbR

      mtau1 = mtauL 
      mtau2 = mtauR

      mu1 = muL 
      mu2 = muR

      md1 = mdL 
      md2 = mdR

C-------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta = datan(tanbeta)
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta = dcos(datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)

!--------------------------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0

!-----------------------------------------------------------

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

!---------------------------------------------------------------    
!---------------------------------------------------------------

      lsdnLsdnLR(1,1) = (g*MZ*gdL*cosbeta/costhw) + 
     $     (ydRG(1,1)**2.d0) * vev1
      lsdnLsdnLR(1,2) = (ydRG(1,1)*ADRG(1,1))/dsqrt(2.d0)

      lsdnLsdnLR(2,1) = ((- g*MZ*gdL*sinbeta)/costhw)
      lsdnLsdnLR(2,2) = (-ydRG(1,1)*sgnmu*modmu)/dsqrt(2.d0)

      lsdnLsdnLR(3,1) = 0.d0
      lsdnLsdnLR(3,2) = (-1.d0/dsqrt(2.d0)) *
     $     (-sgnmu*modmu*sinbeta*ydRG(1,1) + 
     $     ydRG(1,1)*ADRG(1,1)*cosbeta)

      lsdnLsdnLR(4,1) = 0.d0
      lsdnLsdnLR(4,2) = (-1.d0/dsqrt(2.d0)) * 
     $     (-sgnmu*modmu*cosbeta*ydRG(1,1) - 
     $     ydRG(1,1)*ADRG(1,1)*sinbeta)


      lsdnRsdnLR(1,1) = lsdnLsdnLR(1,2)
      lsdnRsdnLR(1,2) = (g*MZ*gdR*cosbeta/costhw) + 
     $     ydRG(1,1)**2.d0 * vev1

      lsdnRsdnLR(2,1) = lsdnLsdnLR(2,2)
      lsdnRsdnLR(2,2) = -g*MZ*gdR*sinbeta/costhw

      lsdnRsdnLR(3,1) = -lsdnLsdnLR(3,2)
      lsdnRsdnLR(3,2) = 0.d0

      lsdnRsdnLR(4,1) = -lsdnLsdnLR(4,2)
      lsdnRsdnLR(4,2) = 0.d0

!--------------------------------------------------------------------------------  
  
      loopmixst: DO i = 1, 4

      intm1(1) = lsdnLsdnLR(i,1)
      intm1(2) = lsdnLsdnLR(i,2)

      call rmat2d(thetad,rthetad)
            
      intm2(1) = rthetad(1,1)*intm1(1) + rthetad(1,2)*intm1(2)
      intm2(2) = rthetad(2,1)*intm1(1) + rthetad(2,2)*intm1(2)

      lsdnLsdn12(i,1) = intm2(1)
      lsdnLsdn12(i,2) = intm2(2)

      intm1(1) = lsdnRsdnLR(i,1)
      intm1(2) = lsdnRsdnLR(i,2)

      intm2(1) = rthetad(1,1)*intm1(1) + rthetad(1,2)*intm1(2)
      intm2(2) = rthetad(2,1)*intm1(1) + rthetad(2,2)*intm1(2)

      lsdnRsdn12(i,1) = intm2(1)
      lsdnRsdn12(i,2) = intm2(2)

      ENDDO loopmixst
!--------------------------------------------------------------------------
C       Mixing CP-even Higgs 
!----------------------------------------------------------------------

      loopmixcpeven: DO i = 1, 2
      
      intm1(1) = lsdnLsdn12(1,i)
      intm1(2) = lsdnLsdn12(2,i)

      call rmat2d(alpha,ralpha)
            
      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lHsdnLsdn12(1,i) = intm2(1)
      lHsdnLsdn12(2,i) = intm2(2)

      intm1(1) = lsdnRsdn12(1,i)
      intm1(2) = lsdnRsdn12(2,i)

      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lHsdnRsdn12(1,i) = intm2(1)
      lHsdnRsdn12(2,i) = intm2(2)
     
      ENDDO loopmixcpeven

!----------------------------------------------------------------------
C    Feynman rules - charged higgs
!<--------- (H+ G+, L R) basis


      lHcsuplsdnlr(1, 1) = (g*MW*sin2beta - 
     $     (yuRG(1,1)**2.d0) * vev2 * cosbeta - 
     $     (ydRG(1,1)**2.d0) * vev1 * sinbeta)/dsqrt(2.d0)

      lHcsuplsdnlr(1,2) = (-sgnmu*modmu*yuRG(1,1)*sinbeta -
     $     AURG(1,1)*yuRG(1,1)*cosbeta)


      lHcsuplsdnlr(2,1) = (-g*MW*cos2beta - 
     $     yuRG(1,1)**2.d0 * vev2 * sinbeta + 
     $     ydRG(1,1)**2.d0 * vev1 * cosbeta)/dsqrt(2.d0)

      lHcsuplsdnlr(2,2) = - yuRG(1,1) * (-sgnmu*modmu*cosbeta + 
     $     AURG(1,1)*sinbeta)

!-------------------------------------------------------------------------

      intm1(1) = lHcsuplsdnlr(1,1)
      intm1(2) = lHcsuplsdnlr(1,2)

      call rmat2d(thetau,rthetau)
            
      intm2(1) = rthetau(1,1)*intm1(1) + rthetau(1,2)*intm1(2)
      intm2(2) = rthetau(2,1)*intm1(1) + rthetau(2,2)*intm1(2)

 
      lHcsuplsdn12(1, 1) = intm2(1)
      lHcsuplsdn12(1, 2) = intm2(2)

      intm1(1) = lHcsuplsdnlr(2, 1)
      intm1(2) = lHcsuplsdnlr(2, 2)

           
      intm2(1) = rthetau(1,1)*intm1(1) + rthetau(1,2)*intm1(2)
      intm2(2) = rthetau(2,1)*intm1(1) + rthetau(2,2)*intm1(2)

      lHcsuplsdn12(2, 1) = intm2(1)
      lHcsuplsdn12(2, 2) = intm2(2)

!---------------------------------------------------------------------------
!<------------ (H+ G+, L R) basis

      lHcsuprsdnlr(1, 1) = ydRG(1,1)*(-sgnmu*modmu*cosbeta - 
     $     ADRG(1,1)*sinbeta)

      lHcsuprsdnlr(1, 2) = ydRG(1,1)*yuRG(1,1)*(- vev2*sinbeta - 
     $     vev1*cosbeta)/dsqrt(2.d0)

      lHcsuprsdnlr(2, 1) = ydRG(1,1)*(-sgnmu*modmu*sinbeta +
     $     ADRG(1,1)*cosbeta)  

      lHcsuprsdnlr(2, 2) = 0.d0

!-------------------------------------------------------------------------

      intm1(1) = lHcsuprsdnlr(1,1)
      intm1(2) = lHcsuprsdnlr(1,2)

      call rmat2d(thetau,rthetau)
            
      intm2(1) = rthetau(1,1)*intm1(1) + rthetau(1,2)*intm1(2)
      intm2(2) = rthetau(2,1)*intm1(1) + rthetau(2,2)*intm1(2)

      lHcsuprsdn12(1,1) = intm2(1)
      lHcsuprsdn12(1,2) = intm2(2)

      intm1(1) = lHcsuprsdnlr(2,1)
      intm1(2) = lHcsuprsdnlr(2,2)

            
      intm2(1) = rthetau(1,1)*intm1(1) + rthetau(1,2)*intm1(2)
      intm2(2) = rthetau(2,1)*intm1(1) + rthetau(2,2)*intm1(2)
  
      lHcsuprsdn12(2, 1) = intm2(1)
      lHcsuprsdn12(2, 2) = intm2(2)

!------------------------------------------------------------------------
C     Feynman rules - neutralino
!------------------------------------------------------------------------

      aPsi0sdnr(1) = (2.d0*gp/(3.d0*dsqrt(2.d0)))
      aPsi0sdnr(2) = 0.d0
      aPsi0sdnr(3) = 0.d0
      aPsi0sdnr(4) = 0.d0

      bPsi0sdnl(1) = gp/(3.d0*dsqrt(2.d0))
      bPsi0sdnl(2) = - g/dsqrt(2.d0)
      bPsi0sdnl(3) = 0.d0
      bPsi0sdnl(4) = 0.d0

      aPsi0sdnl(1) = 0.d0
      aPsi0sdnl(2) = 0.d0
      aPsi0sdnl(3) = 0.d0 
      aPsi0sdnl(4) = ydRG(1,1)

      bPsi0sdnr(1) = 0.d0
      bPsi0sdnr(2) = 0.d0
      bPsi0sdnr(3) = 0.d0
      bPsi0sdnr(4) = ydRG(1,1)


      aChi0sdnl(1) = ON(1,1)*aPsi0sdnl(1) + ON(1,2)*aPsi0sdnl(2) +
     $                 ON(1,3)*aPsi0sdnl(1) + ON(1,4)*aPsi0sdnl(4) 
      aChi0sdnl(2) = ON(2,1)*aPsi0sdnl(1) + ON(2,2)*aPsi0sdnl(2) +
     $                 ON(2,3)*aPsi0sdnl(1) + ON(2,4)*aPsi0sdnl(4) 
      aChi0sdnl(3) = ON(3,1)*aPsi0sdnl(1) + ON(3,2)*aPsi0sdnl(2) +
     $                 ON(3,3)*aPsi0sdnl(1) + ON(3,4)*aPsi0sdnl(4) 
      aChi0sdnl(4) = ON(4,1)*aPsi0sdnl(1) + ON(4,2)*aPsi0sdnl(2) +
     $                 ON(4,3)*aPsi0sdnl(1) + ON(4,4)*aPsi0sdnl(4) 


      bChi0sdnl(1) = ON(1,1)*bPsi0sdnl(1) + ON(1,2)*bPsi0sdnl(2) +
     $                 ON(1,3)*bPsi0sdnl(1) + ON(1,4)*bPsi0sdnl(4) 
      bChi0sdnl(2) = ON(2,1)*bPsi0sdnl(1) + ON(2,2)*bPsi0sdnl(2) +
     $                 ON(2,3)*bPsi0sdnl(1) + ON(2,4)*bPsi0sdnl(4) 
      bChi0sdnl(3) = ON(3,1)*bPsi0sdnl(1) + ON(3,2)*bPsi0sdnl(2) +
     $                 ON(3,3)*bPsi0sdnl(1) + ON(3,4)*bPsi0sdnl(4)  
      bChi0sdnl(4) = ON(4,1)*bPsi0sdnl(1) + ON(4,2)*bPsi0sdnl(2) +
     $                 ON(4,3)*bPsi0sdnl(1) + ON(4,4)*bPsi0sdnl(4) 




      aChi0sdnr(1) = ON(1,1)*aPsi0sdnr(1) + ON(1,2)*aPsi0sdnr(2) +
     $                 ON(1,3)*aPsi0sdnr(3) + ON(1,4)*aPsi0sdnr(4)
      aChi0sdnr(2) = ON(2,1)*aPsi0sdnr(1) + ON(2,2)*aPsi0sdnr(2) +
     $                 ON(2,3)*aPsi0sdnr(3) + ON(2,4)*aPsi0sdnr(4)
      aChi0sdnr(3) = ON(3,1)*aPsi0sdnr(1) + ON(3,2)*aPsi0sdnr(2) +
     $                 ON(3,3)*aPsi0sdnr(3) + ON(3,4)*aPsi0sdnr(4)
      aChi0sdnr(4) = ON(4,1)*aPsi0sdnr(1) + ON(4,2)*aPsi0sdnr(2) +
     $                 ON(4,3)*aPsi0sdnr(3) + ON(4,4)*aPsi0sdnr(4)


      bChi0sdnr(1) = ON(1,1)*bPsi0sdnr(1) + ON(1,2)*bPsi0sdnr(2) +
     $                 ON(1,3)*bPsi0sdnr(3) + ON(1,4)*bPsi0sdnr(4) 
      bChi0sdnr(2) = ON(2,1)*bPsi0sdnr(1) + ON(2,2)*bPsi0sdnr(2) +
     $                 ON(2,3)*bPsi0sdnr(3) + ON(2,4)*bPsi0sdnr(4) 
      bChi0sdnr(3) = ON(3,1)*bPsi0sdnr(1) + ON(3,2)*bPsi0sdnr(2) +
     $                 ON(3,3)*bPsi0sdnr(3) + ON(3,4)*bPsi0sdnr(4) 
      bChi0sdnr(4) = ON(4,1)*bPsi0sdnr(1) + ON(4,2)*bPsi0sdnr(2) +
     $                 ON(4,3)*bPsi0sdnr(3) + ON(4,4)*bPsi0sdnr(4) 


      loopchtllrr: DO i = 1, 4

      fChi0dsdnLL(i) = (aChi0sdnl(i)*aChi0sdnl(i) + 
     $     bChi0sdnl(i)*bChi0sdnl(i))
      gChi0dsdnLL(i) = (bChi0sdnl(i)*aChi0sdnl(i) + 
     $     bChi0sdnl(i)*aChi0sdnl(i))
      fChi0dsdnRR(i) = (aChi0sdnr(i)*aChi0sdnr(i) + 
     $     bChi0sdnr(i)*bChi0sdnr(i))
      gChi0dsdnRR(i) = (bChi0sdnr(i)*aChi0sdnr(i) + 
     $     bChi0sdnr(i)*aChi0sdnr(i))
      fChi0dsdnLR(i) = (aChi0sdnr(i)*aChi0sdnl(i) + 
     $     bChi0sdnr(i)*bChi0sdnl(i))
      gChi0dsdnLR(i) = (bChi0sdnl(i)*aChi0sdnr(i) + 
     $     bChi0sdnr(i)*aChi0sdnl(i))
      
      ENDDO loopchtllrr

!---------------------------------------------------------------------------
c     Feynman Rules for chargino
!---------------------------------------------------------------------------

      aPsicsdnl(1) = 0.d0
      aPsicsdnl(2) = - yuRG(1,1)

      aPsicsdnr(1) = 0.d0
      aPsicsdnr(2) = 0.d0

      bPsicsdnl(1) = g
      bPsicsdnl(2) = 0.d0

      bPsicsdnr(1) = 0.d0
      bPsicsdnr(2) = - ydRG(1,1)


      aChicsdnl(1) = OCR(1,1)*aPsicsdnl(1) + OCR(1,2)*aPsicsdnl(2)
      aChicsdnl(2) = OCR(2,1)*aPsicsdnl(1) + OCR(2,2)*aPsicsdnl(2)

      bChicsdnl(1) = OCL(1,1)*bPsicsdnl(1) + OCL(1,2)*bPsicsdnl(2)
      bChicsdnl(2) = OCL(2,1)*bPsicsdnl(1) + OCL(2,2)*bPsicsdnl(2)

      aChicsdnr(1) = OCR(1,1)*aPsicsdnr(1) + OCR(1,2)*aPsicsdnr(2)
      aChicsdnr(2) = OCR(2,1)*aPsicsdnr(1) + OCR(2,2)*aPsicsdnr(2)

      bChicsdnr(1) = OCL(1,1)*bPsicsdnr(1) + OCL(1,2)*bPsicsdnr(2)
      bChicsdnr(2) = OCL(2,1)*bPsicsdnr(1) + OCL(2,2)*bPsicsdnr(2)

      
      loopchst: DO i = 1, 2

      fChsdnLL(i) = (aChicsdnl(i)*aChicsdnl(i) +
     $     bChicsdnl(i)*bChicsdnl(i))
      gChsdnLL(i) = (bChicsdnl(i)*aChicsdnl(i) +
     $     aChicsdnl(i)*bChicsdnl(i))
      fChsdnLR(i) = (aChicsdnl(i)*aChicsdnr(i) +
     $     bChicsdnl(i)*bChicsdnr(i))
      gChsdnLR(i) = (bChicsdnl(i)*aChicsdnr(i) +
     $     aChicsdnl(i)*bChicsdnr(i))
      fChsdnRR(i) = (aChicsdnr(i)*aChicsdnr(i) +
     $     bChicsdnr(i)*bChicsdnr(i))
      gChsdnRR(i) = (bChicsdnr(i)*aChicsdnr(i) +
     $     aChicsdnr(i)*bChicsdnr(i))

      ENDDO loopchst

!--------------------------------------------------------------------  
C      Corrections Begin
!--------------------------------------------------------------------

       call funcg(p,M3t,md,q,ggd)

  !     print*,"ggt"
       
       call f(p,md1,0.d0,q,fd10)
       call f(p,md2,0.d0,q,fd20)

       call a0(mu2,q,a0mu2)
       call a0(mu1,q,a0mu1)
       call a0(mcR,q,a0mc2)
       call a0(mcL,q,a0mc1)
       call a0(mt2,q,a0t2)
       call a0(mt1,q,a0t1)

       call a0(msR,q,a0ms2)
       call a0(msL,q,a0ms1)
       call a0(mdR,q,a0md2)
       call a0(mdL,q,a0md1)
       call a0(mb2,q,a0b2)
       call a0(mb1,q,a0b1)

       call a0(meR,q,a0me2)
       call a0(meL,q,a0me1)
       call a0(mmuR,q,a0mmu2)
       call a0(mmuL,q,a0mmu1)
       call a0(mtau2,q,a0mtau2)
       call a0(mtau1,q,a0mtau1)

       call a0(snu(1),q,a0snu1)
       call a0(snu(2),q,a0snu2)
       call a0(snu(3),q,a0snu3)


       call a0(mA0,q,a0mA)
       call a0(mh0,q,a0mh)
       call a0(mHu0,q,a0mHu)
       call a0(mHpm,q,a0mHpm)
       call a0(MW,q,a0Mw)
       call a0(MZ,q,a0Mz)

       call f(p,mu1,MW,q,fmu1Mw)
       call f(p,mu2,MW,q,fmu2Mw)
       call f(p,md1,MZ,q,fmd1Mz)
       call f(p,md2,MZ,q,fmd2Mz)

       call a0(md1,q,a0d1)
       call a0(md2,q,a0d2)

       call b0(p,M3t,md,q,b0M3md)
    
!---initialize to zero
       gtterm(1,1) = 0.d0
       gtterm(1,2) = 0.d0
       gtterm(2,1) = 0.d0
       gtterm(2,2) = 0.d0

       csbtm(1,1) = 0.d0
       csbtm(1,2) = 0.d0
       csbtm(2,1) = 0.d0
       csbtm(2,2) = 0.d0
       
       cstop(1,1) = 0.d0
       cstop(1,2) = 0.d0
       cstop(2,1) = 0.d0
       cstop(2,2) = 0.d0

       higgsterm(1,1) = 0.d0
       higgsterm(1,2) = 0.d0
       higgsterm(2,1) = 0.d0
       higgsterm(2,2) = 0.d0

!--------------------------------------------------------------
      gtterm(1, 1) = (4.d0 * (g3**2.d0)/3.d0) *
     $      (2.d0 * ggd + (costhetad)**2.d0 * (fd10 + a0md1) + 
     $      (sinthetad)**2.d0 * (fd20 + a0md2))

      gtterm(2, 2) = (4.d0 * (g3**2.d0)/3.d0) *
     $     (2.d0 * ggd + (sinthetad)**2.d0 * (fd10 + a0d1) + 
     $     (costhetad)**2.d0 * (fd20 + a0md2))

      gtterm(1, 2) = (4.d0 * (g3**2.d0)/3.d0) *
     $     (4.d0 * M3t*md * b0M3md + 
     $     sinthetad*costhetad * (fd10 - a0md1 - fd20 + a0md2))

!----------------------------------------------------------------

      csbtm(1, 1) = (ydRG(1,1)**2.d0) * ((sinthetad)**2.d0 * a0d1 + 
     $   (costhetad)**2.d0 * a0d2)

      csbtm(2, 2) = (ydRG(1,1)**2.d0) * ((costhetad)**2.d0 * a0d1 + 
     $     (sinthetad)**2.d0 * a0d2)

      csbtm(1, 2) = (ydRG(1,1)**2.d0) * costhetad*sinthetad * 3.d0 *
     $     (a0d1 - a0d2)
!------------------------------------------------------------------

      cstop(1, 1) = (yuRG(1,1)**2.d0) * ((sinthetau)**2.d0 * a0mu1 + 
     $     (costhetau)**2.d0 * a0mu2)

      cstop(2, 2) = (yuRG(1,1)**2.d0) * ((costhetau)**2.d0 * a0mu1 + 
     $     (sinthetau)**2.d0 * a0mu2)
!--------------------------------------------------------------------

      higgsterm(1,1) = 0.5d0 * ((ydRG(1,1)**2.d0) * dnd(1) - 
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((ydRG(1,1)**2.d0) * dnd(2) - 
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((ydRG(1,1)**2.d0) * dnd(3) - 
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((ydRG(1,1)**2.d0) * dnd(4) -
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA

      higgsterm(2,2) = 0.5d0 * ((ydRG(1,1)**2.d0) * dnd(1) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((ydRG(1,1)**2.d0) * dnd(2) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((ydRG(1,1)**2.d0) * dnd(3) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((ydRG(1,1)**2.d0) * dnd(4) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA

      higgsterm(1,2) = 0.d0

      higgsterm(1,1) = higgsterm(1,1) + 
     $     ((yuRG(1,1)**2.d0) * dnd(3) + ((g*g) * 
     $     ((gdL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(3))) * a0mHpm +
     $     ((yuRG(1,1)**2.d0) * dnd(4) + ((g*g) * 
     $     ((gdL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(4))) * a0MW 

      higgsterm(2,2) = higgsterm(2,2) + 
     $     ((ydRG(1,1)**2.d0) * dnu(3) + 
     $     ((g*g) * gdR*0.5d0/(costhw**2.d0)) * cn(3)) * a0mHpm +
     $     ((ydRG(1,1)**2.d0) * dnu(4) +  
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(4)) * a0MW !<---- checked 

!------------------------------------------------------------------------------

      call b0(p,mHu0,md1,q,b0h0msdn(1,1))
      call b0(p,mHu0,md2,q,b0h0msdn(1,2))
      call b0(p,mh0,md1,q,b0h0msdn(2,1))
      call b0(p,mh0,md2,q,b0h0msdn(2,2))
      call b0(p,MZ,md1,q,b0h0msdn(3,1))
      call b0(p,MZ,md2,q,b0h0msdn(3,2))
      call b0(p,mA0,md1,q,b0h0msdn(4,1))
      call b0(p,mA0,md2,q,b0h0msdn(4,2))
      

      loophi: DO i = 1, 4       !<-------------neutral higgs terms
      loophj: DO j = 1, 2

      higgsterm(1,1) = higgsterm(1,1) + 
     $     (lHsdnLsdn12(i,j)**2.d0) * b0h0msdn(i,j)

      higgsterm(1,2) = higgsterm(1,2) + 
     $     lHsdnLsdn12(i,j) * lHsdnRsdn12(i,j) * b0h0msdn(i,j)

      higgsterm(2,2) = higgsterm(2,2) + 
     $     (lHsdnRsdn12(i,j)**2.d0) * b0h0msdn(i,j)
      

      ENDDO loophj
      ENDDO loophi
!-------------------------------------------------------------------------------

      call b0(p,mu1,mHpm,q,b0hcmsu(1,1))
      call b0(p,mu2,mHpm,q,b0hcmsu(1,2))
      call b0(p,mu1,MW,q,b0hcmsu(2,1))
      call b0(p,mu2,MW,q,b0hcmsu(2,2))
      

      loophci: DO i = 1, 2      !<-------------charged higgs terms
      loophcj: DO j = 1, 2

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (lHcsuplsdn12(i,j)**2.d0) * b0hcmsu(i,j)        

      higgsterm(1, 2) = higgsterm(1,2) + 
     $     lHcsuplsdn12(i, j) * lHcsuprsdn12(i, j) * b0hcmsu(i,j)

      higgsterm(2, 2) = higgsterm(2, 2) +
     $     (lHcsuprsdn12(i,j)**2.d0) * b0hcmsu(i,j)
      

      ENDDO loophcj
      ENDDO loophci


!---------------------------------------------------

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (4.d0*(g*g)/costhw**2.d0)*(gdL*gdL)*a0MZ  + 
     $     2.d0*(g*g)*a0MW + 
     $     ((1.d0/9.d0)*g*g*sinsqthw) * ((costhetad**2.d0) * 
     $     fd10 + (sinthetad**2.d0) * fd20) +
     $     (g*gdL/costhw)**2.d0 * ((costhetad**2.d0) * fmd1Mz +
     $     (sinthetad**2.d0) * fmd2mz) +
     $     (g*g)*0.5d0*((costhetau)**2.d0 * fmu1Mw + (sinthetau**2.d0) *
     $     fmu2Mw) +
     $     (g*g) * 0.25d0 * (1.d0 * ((costhetad**2.d0) * a0md1 + 
     $     (sinthetad**2.d0) * a0md2) + 2.d0* ((costhetau**2.d0)*a0mu1 + 
     $     (sinthetau**2.d0) * a0mu2)) -
     $     (g*g) * 0.5d0 * (- 1.5d0 * a0md1 - 1.5d0*a0ms1 +
     $     1.5d0 * a0mu1 + 1.5d0*a0mc1 +
     $     0.5d0 * (a0snu1 + a0snu2 + a0snu3) -
     $     0.5d0 * (a0me1 + a0mmu1 + ((costhetatau**2.d0) * a0mtau1 +
     $     (sinthetatau**2.d0) * a0mtau2))) +
     $     (gp**2.d0) * 0.25d0 *(ydL**2.d0) * ((costhetad**2.d0) *a0md1+
     $     (sinthetad**2.d0)*a0md2) +
     $     (gp**2.d0) * 0.25d0 * ydL * ((3.d0*ydL*(a0md1 + a0ms1 + 
     $     ((costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2))) +
     $     (3.d0 * ydR *(a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2)) +
     $     (3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t1)) +
     $     (3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (costhetat**2.d0)*a0t2 + (sinthetat**2.d0)*a0t1)) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))



      higgsterm(2, 2) = higgsterm(2,2) +
     $     (4.d0*(g*g)/costhw**2.d0) * (gdR**2.d0) * a0MZ + 
     $     ((1.d0/9.d0)*g*g*sinsqthw) * ((sinthetad)*fd10 + 
     $     (costhetad**2.d0)*fd20) +
     $     (g*gdR/costhw)**2.d0 * ((sinthetad**2.d0) * fmd1MZ + 
     $     (costhetad**2.d0) * fmd2MZ) +
     $     (gp**2.d0) * 0.25d0 * (ydR**2.d0) * 
     $     ((sinthetad**2.d0) * a0d1 + (costhetad**2.d0) * a0d2) +
     $     (gp**2.d0) * 0.25d0 * ydR *  
     $     (3.d0 * ydL * (a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2) +
     $     3.d0 * ydR * (a0md2 +  a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2) +
     $     3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t2) +
     $     3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 + (costhetat**2.d0)*a0t2) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))



      higgsterm(1,2) = higgsterm(1,2) + 
     $     (gp**2.d0) * ydL*ydR * sinthetad*costhetad *
     $     (a0d1 - a0d2) +
     $     ((1.d0/9.d0) * g*g * sinsqthw) * sinthetad*costhetad * 
     $     (fd10 - fd20) -
     $     ((g*g)/costhw**2.d0) * gdL*gdR * sinthetad*costhetad *
     $     (fmd1Mz - fmd2Mz)    !<-----------------------chked

!----------------------------------------------------------------------
C     Chargino term
!----------------------------------------------------------------------

      call funcg(p,mchargino(1),mu,q,gmchargino1mu)
      call funcg(p,mchargino(2),mu,q,gmchargino2mu)
      call b0(p,mchargino(1),mu,q,b0mchargino1mu)
      call b0(p,mchargino(2),mu,q,b0mchargino2mu)
      
      chargino(1, 1) =  fChsdnLL(1)*gmchargino1mu -
     $     gChsdnLL(1)*mchargino(1)*mu*b0mchargino1mu*2.d0 +
     $     fChsdnLL(2)*gmchargino2mu -
     $     gChsdnLL(2)*mchargino(2)*mu*b0mchargino2mu*2.d0

      chargino(1, 2) =  fChsdnLR(1)*gmchargino1mu -
     $     gChsdnLR(1)*mchargino(1)*mu*b0mchargino1mu*2.d0 +
     $     fChsdnLR(2)*gmchargino2mu -
     $     gChsdnLR(2)*mchargino(2)*mu*b0mchargino2mu*2.d0

      chargino(2, 2) =  fChsdnRR(1)*gmchargino1mu - 
     $     gChsdnRR(1)*mchargino(1)*mu*b0mchargino1mu*2.d0 +
     $     fChsdnRR(2)*gmchargino2mu - 
     $     gChsdnRR(2)*mchargino(2)*mu*b0mchargino2mu*2.d0

!---------------------------------------------------------------------
C     neutralino terms
!---------------------------------------------------------------------

      call b0(p,nmneut(1),md,q,b0mneut1md)
      call b0(p,nmneut(2),md,q,b0mneut2md)
      call b0(p,nmneut(3),md,q,b0mneut3md)
      call b0(p,nmneut(4),md,q,b0mneut4md)
      
      call funcg(p,nmneut(1),md,q,gmneut1md)
      call funcg(p,nmneut(2),md,q,gmneut2md)
      call funcg(p,nmneut(3),md,q,gmneut3md)
      call funcg(p,nmneut(4),md,q,gmneut4md)

      neutralino(1, 1) = 
     $     fChi0dsdnLL(1)*gmneut1md - gChi0dsdnLL(1)*2.d0*
     $     mneut(1)*md*b0mneut1md +
     $     fChi0dsdnLL(2)*gmneut2md - gChi0dsdnLL(2)*2.d0*
     $     mneut(2)*md*b0mneut2md +
     $     fChi0dsdnLL(3)*gmneut3md - gChi0dsdnLL(3)*2.d0*
     $     mneut(3)*md*b0mneut3md +
     $     fChi0dsdnLL(4)*gmneut4md - gChi0dsdnLL(4)*2.d0*
     $     mneut(4)*md*b0mneut4md


      neutralino(2, 2) = 
     $     fChi0dsdnRR(i)*gmneut1md - gChi0dsdnRR(1)*2.d0*
     $     mneut(1)*md*b0mneut1md +
     $     fChi0dsdnRR(2)*gmneut2md - gChi0dsdnRR(2)*2.d0*
     $     mneut(2)*md*b0mneut2md +
     $     fChi0dsdnRR(3)*gmneut3md - gChi0dsdnRR(3)*2.d0*
     $     mneut(3)*md*b0mneut3md +
     $     fChi0dsdnRR(4)*gmneut4md - gChi0dsdnRR(4)*2.d0*
     $     mneut(4)*md*b0mneut4md

      neutralino(1, 2) = 
     $     fChi0dsdnLR(1)*gmneut1md - gChi0dsdnLR(1)*2.d0*
     $     mneut(1)*md*b0mneut1md +
     $     fChi0dsdnLR(2)*gmneut2md - gChi0dsdnLR(2)*2.d0*
     $     mneut(2)*md*b0mneut2md +
     $     fChi0dsdnLR(3)*gmneut3md - gChi0dsdnLR(3)*2.d0*
     $     mneut(3)*md*b0mneut3md +
     $     fChi0dsdnLR(4)*gmneut4md - gChi0dsdnLR(4)*2.d0*
     $     mneut(4)*md*b0mneut4md 
      
!-------------------------------------------------------------------

      MSQD1(1,1) = mSQRG(1,1) + MD**2.d0 + gdL*MZ*MZ*dcos(2.d0*beta)
      MSQD1(1,2) = MD*((ADRG(1,1) - sgnmu*modmu*dtan(beta)))
      MSQD1(2,1) = MSQD1(1,2)
      MSQD1(2,2) = mSDRG(1,1) + MD**2.d0 + gdR*MZ*MZ*dcos(2.d0*beta)

!----------------------------------------------------------------------

      pisdn(1,1) = (cstop(1,1) + csbtm(1,1) + gtterm(1,1) +
     $     chargino(1,1) + neutralino(1,1) + higgsterm(1,1))/
     $     (16.d0*pi*pi)

      pisdn(1,2) = (cstop(1,2) + csbtm(1,2) + gtterm(1,2)  +
     $     chargino(1,2)+ neutralino(1,2) + higgsterm(1,2))/
     $     (16.d0*pi*pi)

      pisdn(2,2) = (cstop(2,2) + csbtm(2,2) + gtterm(2,2) +
     $     chargino(2,2)+ neutralino(2,2) +  higgsterm(2,2))/
     $     (16.d0*pi*pi)
      
      pisdn(2,1) = pisdn(1,2)

!----------------------------------------------------------------------

      pisdn(1,1) = MSQD1(1,1) - pisdn(1,1)
      pisdn(1,2) = MSQD1(1,2) - pisdn(1,2)
      pisdn(2,2) = MSQD1(2,2) - pisdn(2,2)
      pisdn(2,1) = pisdn(1,2) 

!----------------------------------------------------------------------

C     find the singular values and the diagonalising matrices

      info = 10
      AOK = 0
      
    !  call dsyev('V','U',2,pisdn,2,SDeg,work,lwork,info)
      
      Call CEigensystem(2,pisdn,2,SDeg,sdmix,2,1)
     
      if(info.eq.0) then
         AOK = AOK + 1
      endif

      
      RETURN

      END SUBROUTINE pisdown

C==================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C     pisstrange is checked on 24/05/2010 @ 18:30.
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      SUBROUTINE pisstrange(p,q,g,gp,g3,tanbeta,mSQRG,mSDRG,
     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3t,delpisst,SSTeg)

 
      IMPLICIT NONE 

      INTEGER i,j,info,AOK,lwork
      parameter(lwork = 35)
      DOUBLE PRECISION work(lwork)
      DOUBLE PRECISION tanbeta,mSQRG(3,3),mSDRG(3,3)
      DOUBLE PRECISION AURG(3,3),ADRG(3,3)
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yuRG(3,3),ydRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 

      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR,MSQD3(2,2)      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g,g3,M3t

      double precision beta
      double precision vev1,vev2
 
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION costhetac,sinthetac,costhetas,sinthetas

      DOUBLE PRECISION thetau,thetad,thetae

      DOUBLE PRECISION pisstr(2,2)
 
      DOUBLE PRECISION higgsterm(2,2)
      DOUBLE PRECISION dnu(4),dnd(4),cn(4)

      DOUBLE PRECISION lsstrLstrLR(4,2),lsstrLstr12(4,2)
      DOUBLE PRECISION lsstrRsstrLR(4,2),lsstrRsstr12(4,2)
      
      DOUBLE PRECISION aPsi0sstrr(4), bPsi0sstrr(4), aPsi0sstrl(4)
      DOUBLE PRECISION bPsi0sstrl(4) 
      DOUBLE PRECISION aChi0sstrl(4), bChi0sstrl(4), aChi0sstrr(4)
      DOUBLE PRECISION bChi0sstrr(4)

      DOUBLE PRECISION gChi0strsstrLL(4), fChi0strsstrLL(4)
      DOUBLE PRECISION gChi0strsstrLR(4), fChi0strsstrLR(4)
      DOUBLE PRECISION gChi0strsstrRR(4), fChi0strsstrRR(4)

      DOUBLE PRECISION bPsicstrl(2), bPsicstrr(2), aPsicstrl(2)
      DOUBLE PRECISION aPsicstrr(2)
      DOUBLE PRECISION aChicstrr(2), aChicstrl(2)
      DOUBLE PRECISION bChicstrl(2),bChicstrr(2)

      DOUBLE PRECISION gtterm(2,2), cstop(2,2), csbtm(2,2) 
      data gtterm/ 4 * 0.d0/, cstop/ 4 * 0.d0/, csbtm/ 4 * 0.d0/

      DOUBLE PRECISION  chargino(2,2), neutralino(2,2)
      data chargino/ 4 * 0.d0/,neutralino/ 4 * 0.d0/

      DOUBLE PRECISION fChstrLL(2), gChstrLL(2) 
      DOUBLE PRECISION fChstrLR(2), gChstrLR(2) 
      DOUBLE PRECISION fChstrRR(2), gChstrRR(2)
 
      DOUBLE PRECISION intm1(2), intm2(2)

      DOUBLE PRECISION lHcstrrsclr(2,2),lHcstrrsc12(2,2)
      DOUBLE PRECISION lHcstrlsc12(2,2)
      DOUBLE PRECISION lHcstrlsclr(2,2) 

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION lHsstrLstr12(4,2),lHsstrRsstr12(4,2)
      data lHsstrLstr12/ 8 * 0.d0/, lHsstrRsstr12/ 8 * 0.d0/
      
      DOUBLE PRECISION rthetas(2,2),rthetac(2,2),ralpha(2,2)

      DOUBLE PRECISION ggs,fms10,fms20,a0b1,a0b2,b0M3ms,a0t1,a0t2
      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz
      DOUBLE PRECISION b0h0mstr(4,2),b0hcmscharm(2,2)

      DOUBLE PRECISION a0md1,a0md2,a0ms1,a0ms2,a0mtau1,a0mtau2,a0snu2
      DOUBLE PRECISION a0me1,a0me2,a0mmu1,a0mmu2,a0snu1,a0snu3
      DOUBLE PRECISION a0mc1,a0mc2,a0mu1,a0mu2,nmneut(4)
      
      DOUBLE PRECISION fms1MZ,fms2MZ,fmc1MW,fmc2MW
      DOUBLE PRECISION gmchargino1mc,gmchargino2mc
      DOUBLE PRECISION b0mchargino1mc,b0mchargino2mc
      DOUBLE PRECISION b0mneut1ms,b0mneut2ms,b0mneut3ms,b0mneut4ms
      DOUBLE PRECISION gmneut1ms,gmneut2ms,gmneut3ms,gmneut4ms
      
      DOUBLE PRECISION SSTeg(2),sstmix(2,2)

      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2,mc1,mc2,
     $     ms1,ms2
      DOUBLE PRECISION delpisst(2,2)

      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae


      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d
  

      include 'stdinputs.h'
C------------------------------------------------------------

      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      MZ = MZpole

      SSTeg(1) = 0.d0
      SSTeg(2) = 0.d0

      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw = dsqrt(sinsqthw)
      thw = dasin(sinthw)
      costhw = dcos(thw)
      cos2thw = dcos(2.d0*thw)
      

      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))

      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))     

       
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
      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)

      costhetas = dcos(thetas)
      sinthetas = dsin(thetas)

      costhetac = dcos(thetac)
      sinthetac = dsin(thetac)

      mt1 = mtL
      mt2 = mtR
      
      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR

      mc1 = mcL
      mc2 = mcR
      
      ms1 = msL
      ms2 = msR

C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta     = datan(tanbeta)
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha  = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta  = dcos(datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)

!-----------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0

!----------------------------------------------------------

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

!---------------------------------------------------------------    

      lsstrLstrLR(1,1) = (g*MZ*gdL*cosbeta/costhw) + 
     $     (ydRG(2,2)**2.d0) * vev1
      lsstrLstrLR(1,2) = (ydRG(2,2)*ADRG(2,2))/dsqrt(2.d0)

      lsstrLstrLR(2,1) = ((- g*MZ*gdL*sinbeta)/costhw)
      lsstrLstrLR(2,2) = (-ydRG(2,2)*sgnmu*modmu)/dsqrt(2.d0)

      lsstrLstrLR(3,1) = 0.d0
      lsstrLstrLR(3,2) = (-1.d0/dsqrt(2.d0)) *
     $     (-sgnmu*modmu*sinbeta*ydRG(2,2) + 
     $     ydRG(2,2)*ADRG(2,2)*cosbeta)

      lsstrLstrLR(4,1) = 0.d0
      lsstrLstrLR(4,2) = (-1.d0/dsqrt(2.d0)) * 
     $     (-sgnmu*modmu*cosbeta*ydRG(2,2) - 
     $     ydRG(2,2)*ADRG(2,2)*sinbeta)



      lsstrRsstrLR(1,1) = lsstrLstrLR(1,2)
      lsstrRsstrLR(1,2) = (g*MZ*gdR*cosbeta/costhw) + 
     $     ydRG(2,2)**2.d0 * vev1

      lsstrRsstrLR(2,1) = lsstrLstrLR(2,2)
      lsstrRsstrLR(2,2) = -g*MZ*gdR*sinbeta/costhw

      lsstrRsstrLR(3,1) = -lsstrLstrLR(3,2)
      lsstrRsstrLR(3,2) = 0.d0

      lsstrRsstrLR(4,1) = -lsstrLstrLR(4,2)
      lsstrRsstrLR(4,2) = 0.d0

!--------------------------------------------------------  
  
      loopmixst: DO i = 1, 4

      intm1(1) = lsstrLstrLR(i,1)
      intm1(2) = lsstrLstrLR(i,2)

      call rmat2d(thetas,rthetas)
            
      intm2(1) = rthetas(1,1)*intm1(1) + rthetas(1,2)*intm1(2)
      intm2(2) = rthetas(2,1)*intm1(1) + rthetas(2,2)*intm1(2)

      lsstrLstr12(i,1) = intm2(1)
      lsstrLstr12(i,2) = intm2(2)

      intm1(1) = lsstrRsstrLR(i,1)
      intm1(2) = lsstrRsstrLR(i,2)

      intm2(1) = rthetas(1,1)*intm1(1) + rthetas(1,2)*intm1(2)
      intm2(2) = rthetas(2,1)*intm1(1) + rthetas(2,2)*intm1(2)

      lsstrRsstr12(i,1) = intm2(1)
      lsstrRsstr12(i,2) = intm2(2)

      ENDDO loopmixst
!--------------------------------------------------------------------------
C       Mixing CP-even Higgs 
!--------------------------------------------------------------------------

      loopmixcpeven: DO i = 1, 2
      
      intm1(1) = lsstrLstr12(1,i)
      intm1(2) = lsstrLstr12(2,i)

      call rmat2d(alpha,ralpha)
            
      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lHsstrLstr12(1,i) = intm2(1)
      lHsstrLstr12(2,i) = intm2(2)

      intm1(1) = lsstrRsstr12(1,i)
      intm1(2) = lsstrRsstr12(2,i)

      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lHsstrRsstr12(1,i) = intm2(1)
      lHsstrRsstr12(2,i) = intm2(2)
     
      ENDDO loopmixcpeven
!----------------------------------------------------------------------
C    Feynman rules - charged higgs
!<--------- (H+ G+, L R) basis

      lHcstrlsclr(1, 1) = (g*MW*sin2beta - 
     $     (yuRG(2,2)**2.d0) * vev2 * cosbeta - 
     $     (ydRG(2,2)**2.d0) * vev1 * sinbeta)/dsqrt(2.d0)

      lHcstrlsclr(1,2) = (-sgnmu*modmu*yuRG(2,2)*sinbeta -
     $     AURG(2,2)*yuRG(2,2)*cosbeta)


      lHcstrlsclr(2,1) = (-g*MW*cos2beta - 
     $     yuRG(2,2)**2.d0 * vev2 * sinbeta + 
     $     ydRG(2,2)**2.d0 * vev1 * cosbeta)/dsqrt(2.d0)

      lHcstrlsclr(2,2) = - yuRG(2,2) * (-sgnmu*modmu*cosbeta + 
     $     AURG(2,2)*sinbeta)

!-----------------------------------------------------------------

      intm1(1) = lHcstrlsclr(1,1)
      intm1(2) = lHcstrlsclr(1,2)

      call rmat2d(thetac,rthetac)
            
      intm2(1) = rthetac(1,1)*intm1(1) + rthetac(1,2)*intm1(2)
      intm2(2) = rthetac(2,1)*intm1(1) + rthetac(2,2)*intm1(2)

 
      lHcstrlsc12(1, 1) = intm2(1)
      lHcstrlsc12(1, 2) = intm2(2)

      intm1(1) = lHcstrlsclr(2, 1)
      intm1(2) = lHcstrlsclr(2, 2)


            
      intm2(1) = rthetac(1,1)*intm1(1) + rthetac(1,2)*intm1(2)
      intm2(2) = rthetac(2,1)*intm1(1) + rthetac(2,2)*intm1(2)

      lHcstrlsc12(2, 1) = intm2(1)
      lHcstrlsc12(2, 2) = intm2(2)

!---------------------------------------------------------------------------
!<------------ (H+ G+, L R) basis


      lHcstrrsclr(1, 1) = ydRG(2,2)*(-sgnmu*modmu*cosbeta - 
     $     ADRG(2,2)*sinbeta)

      lHcstrrsclr(1, 2) = - ydRG(2,2)*yuRG(2,2) * (vev2*sinbeta + 
     $     vev1*cosbeta)/dsqrt(2.d0)


      lHcstrrsclr(2, 1) = ydRG(2,2)*(-sgnmu*modmu*sinbeta + 
     $     ADRG(2,2)*cosbeta)  

      lHcstrrsclr(2, 2) = 0.d0


!-------------------------------------------------------------------------

      intm1(1) = lHcstrrsclr(1,1)
      intm1(2) = lHcstrrsclr(1,2)

      call rmat2d(thetac,rthetac)
            
      intm2(1) = rthetac(1,1)*intm1(1) + rthetac(1,2)*intm1(2)
      intm2(2) = rthetac(2,1)*intm1(1) + rthetac(2,2)*intm1(2)

      lHcstrrsc12(1,1) = intm2(1)
      lHcstrrsc12(1,2) = intm2(2)

      intm1(1) = lHcstrrsclr(2,1)
      intm1(2) = lHcstrrsclr(2,2)

            
      intm2(1) = rthetac(1,1)*intm1(1) + rthetac(1,2)*intm1(2)
      intm2(2) = rthetac(2,1)*intm1(1) + rthetac(2,2)*intm1(2)
  
      lHcstrrsc12(2, 1) = intm2(1)
      lHcstrrsc12(2, 2) = intm2(2)

!------------------------------------------------------------------------
C     Feynman rules - neutralino
!------------------------------------------------------------------------

      aPsi0sstrr(1) = (2.d0*gp/(3.d0*dsqrt(2.d0)))
      aPsi0sstrr(2) = 0.d0
      aPsi0sstrr(3) = 0.d0
      aPsi0sstrr(4) = 0.d0

      bPsi0sstrl(1) = gp/(3.d0*dsqrt(2.d0))
      bPsi0sstrl(2) = - g/dsqrt(2.d0)
      bPsi0sstrl(3) = 0.d0
      bPsi0sstrl(4) = 0.d0

      aPsi0sstrl(1) = 0.d0
      aPsi0sstrl(2) = 0.d0
      aPsi0sstrl(3) = 0.d0 
      aPsi0sstrl(4) = ydRG(2,2)

      bPsi0sstrr(1) = 0.d0
      bPsi0sstrr(2) = 0.d0
      bPsi0sstrr(3) = 0.d0
      bPsi0sstrr(4) = ydRG(2,2)


      aChi0sstrl(1) = ON(1,1)*aPsi0sstrl(1) + ON(1,2)*aPsi0sstrl(2) +
     $                 ON(1,3)*aPsi0sstrl(1) + ON(1,4)*aPsi0sstrl(4) 
      aChi0sstrl(2) = ON(2,1)*aPsi0sstrl(1) + ON(2,2)*aPsi0sstrl(2) +
     $                 ON(2,3)*aPsi0sstrl(1) + ON(2,4)*aPsi0sstrl(4) 
      aChi0sstrl(3) = ON(3,1)*aPsi0sstrl(1) + ON(3,2)*aPsi0sstrl(2) +
     $                 ON(3,3)*aPsi0sstrl(1) + ON(3,4)*aPsi0sstrl(4) 
      aChi0sstrl(4) = ON(4,1)*aPsi0sstrl(1) + ON(4,2)*aPsi0sstrl(2) +
     $                 ON(4,3)*aPsi0sstrl(1) + ON(4,4)*aPsi0sstrl(4) 


      bChi0sstrl(1) = ON(1,1)*bPsi0sstrl(1) + ON(1,2)*bPsi0sstrl(2) +
     $                 ON(1,3)*bPsi0sstrl(1) + ON(1,4)*bPsi0sstrl(4) 
      bChi0sstrl(2) = ON(2,1)*bPsi0sstrl(1) + ON(2,2)*bPsi0sstrl(2) +
     $                 ON(2,3)*bPsi0sstrl(1) + ON(2,4)*bPsi0sstrl(4) 
      bChi0sstrl(3) = ON(3,1)*bPsi0sstrl(1) + ON(3,2)*bPsi0sstrl(2) +
     $                 ON(3,3)*bPsi0sstrl(1) + ON(3,4)*bPsi0sstrl(4) 
      bChi0sstrl(4) = ON(4,1)*bPsi0sstrl(1) + ON(4,2)*bPsi0sstrl(2) +
     $                 ON(4,3)*bPsi0sstrl(1) + ON(4,4)*bPsi0sstrl(4) 




      aChi0sstrr(1) = ON(1,1)*aPsi0sstrr(1) + ON(1,2)*aPsi0sstrr(2) +
     $                 ON(1,3)*aPsi0sstrr(3) + ON(1,4)*aPsi0sstrr(4)
      aChi0sstrr(2) = ON(2,1)*aPsi0sstrr(1) + ON(2,2)*aPsi0sstrr(2) +
     $                 ON(2,3)*aPsi0sstrr(3) + ON(2,4)*aPsi0sstrr(4)
      aChi0sstrr(3) = ON(3,1)*aPsi0sstrr(1) + ON(3,2)*aPsi0sstrr(2) +
     $                 ON(3,3)*aPsi0sstrr(3) + ON(3,4)*aPsi0sstrr(4)
      aChi0sstrr(4) = ON(4,1)*aPsi0sstrr(1) + ON(4,2)*aPsi0sstrr(2) +
     $                 ON(4,3)*aPsi0sstrr(3) + ON(4,4)*aPsi0sstrr(4)


      bChi0sstrr(1) = ON(1,1)*bPsi0sstrr(1) + ON(1,2)*bPsi0sstrr(2) +
     $                 ON(1,3)*bPsi0sstrr(3) + ON(1,4)*bPsi0sstrr(4) 
      bChi0sstrr(2) = ON(2,1)*bPsi0sstrr(1) + ON(2,2)*bPsi0sstrr(2) +
     $                 ON(2,3)*bPsi0sstrr(3) + ON(2,4)*bPsi0sstrr(4) 
      bChi0sstrr(3) = ON(3,1)*bPsi0sstrr(1) + ON(3,2)*bPsi0sstrr(2) +
     $                 ON(3,3)*bPsi0sstrr(3) + ON(3,4)*bPsi0sstrr(4) 
      bChi0sstrr(4) = ON(4,1)*bPsi0sstrr(1) + ON(4,2)*bPsi0sstrr(2) +
     $                 ON(4,3)*bPsi0sstrr(3) + ON(4,4)*bPsi0sstrr(4) 


      loopchtllrr: DO i = 1, 4

      fChi0strsstrLL(i) = (aChi0sstrl(i)*aChi0sstrl(i) + 
     $     bChi0sstrl(i)*bChi0sstrl(i))
      gChi0strsstrLL(i) = (bChi0sstrl(i)*aChi0sstrl(i) + 
     $     bChi0sstrl(i)*aChi0sstrl(i))
      fChi0strsstrRR(i) = (aChi0sstrr(i)*aChi0sstrr(i) + 
     $     bChi0sstrr(i)*bChi0sstrr(i))
      gChi0strsstrRR(i) = (bChi0sstrr(i)*aChi0sstrr(i) + 
     $     bChi0sstrr(i)*aChi0sstrr(i))
      fChi0strsstrLR(i) = (aChi0sstrr(i)*aChi0sstrl(i) + 
     $     bChi0sstrr(i)*bChi0sstrl(i))
      gChi0strsstrLR(i) = (bChi0sstrl(i)*aChi0sstrr(i) + 
     $     bChi0sstrr(i)*aChi0sstrl(i))
      
      ENDDO loopchtllrr

!---------------------------------------------------------------------------
c     Feynman Rules for chargino
!---------------------------------------------------------------------------

      aPsicstrl(1) = 0.d0
      aPsicstrl(2) = - yuRG(2,2)

      aPsicstrr(1) = 0.d0
      aPsicstrr(2) = 0.d0

      bPsicstrl(1) = g
      bPsicstrl(2) = 0.d0

      bPsicstrr(1) = 0.d0
      bPsicstrr(2) = - ydRG(2,2)
      

      aChicstrl(1) = OCR(1,1)*aPsicstrl(1) + OCR(1,2)*aPsicstrl(2)
      aChicstrl(2) = OCR(2,1)*aPsicstrl(1) + OCR(2,2)*aPsicstrl(2)

      bChicstrl(1) = OCL(1,1)*bPsicstrl(1) + OCL(1,2)*bPsicstrl(2)
      bChicstrl(2) = OCL(2,1)*bPsicstrl(1) + OCL(2,2)*bPsicstrl(2)

      aChicstrr(1) = OCR(1,1)*aPsicstrr(1) + OCR(1,2)*aPsicstrr(2)
      aChicstrr(2) = OCR(2,1)*aPsicstrr(1) + OCR(2,2)*aPsicstrr(2)

      bChicstrr(1) = OCL(1,1)*bPsicstrr(1) + OCL(1,2)*bPsicstrr(2)
      bChicstrr(2) = OCL(2,1)*bPsicstrr(1) + OCL(2,2)*bPsicstrr(2)

      
      loopchst: DO i =1,2

      fChstrLL(i) = (aChicstrl(i)*aChicstrl(i) +
     $     bChicstrl(i)*bChicstrl(i))
      gChstrLL(i) = (bChicstrl(i)*aChicstrl(i) +
     $     aChicstrl(i)*bChicstrl(i))
      fChstrLR(i) = (aChicstrl(i)*aChicstrr(i) +
     $     bChicstrl(i)*bChicstrr(i))
      gChstrLR(i) = (bChicstrl(i)*aChicstrr(i) +
     $     aChicstrl(i)*bChicstrr(i))
      fChstrRR(i) = (aChicstrr(i)*aChicstrr(i) +
     $     bChicstrr(i)*bChicstrr(i))
      gChstrRR(i) = (bChicstrr(i)*aChicstrr(i) +
     $     aChicstrr(i)*bChicstrr(i))

      ENDDO loopchst

!--------------------------------------------------------------------
C     Corrections Begin
!------------------------------------------------------------------

      call funcg(p,M3t,ms,q,ggs)
  
      call f(p,ms1,0.d0,q,fms10)
      call f(p,ms2,0.d0,q,fms20)

      call a0(muR,q,a0mu2)
      call a0(muL,q,a0mu1)
      call a0(mc2,q,a0mc2)
      call a0(mc1,q,a0mc1)
      call a0(mt2,q,a0t2)
      call a0(mt1,q,a0t1)

      call a0(ms2,q,a0ms2)
      call a0(ms1,q,a0ms1)
      call a0(mdR,q,a0md2)
      call a0(mdL,q,a0md1)
      call a0(mb2,q,a0b2)
      call a0(mb1,q,a0b1)

      call a0(meR,q,a0me2)
      call a0(meL,q,a0me1)
      call a0(mmuR,q,a0mmu2)
      call a0(mmuL,q,a0mmu1)
      call a0(mtau2,q,a0mtau2)
      call a0(mtau1,q,a0mtau1)

      call a0(snu(1),q,a0snu1)
      call a0(snu(2),q,a0snu2)
      call a0(snu(3),q,a0snu3)


      call a0(mA0,q,a0mA)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)
      call a0(mHpm,q,a0mHpm)
      call a0(MW,q,a0Mw)
      call a0(MZ,q,a0Mz)

      call f(p,mc1,MW,q,fmc1Mw)
      call f(p,mc2,MW,q,fmc2Mw)
      call f(p,ms1,MZ,q,fms1Mz)
      call f(p,ms2,MZ,q,fms2Mz)


      call b0(p,M3t,ms,q,b0M3ms)
      

!---initialize to zero
       gtterm(1,1) = 0.d0
       gtterm(1,2) = 0.d0
       gtterm(2,1) = 0.d0
       gtterm(2,2) = 0.d0

       csbtm(1,1) = 0.d0
       csbtm(1,2) = 0.d0
       csbtm(2,1) = 0.d0
       csbtm(2,2) = 0.d0
       
       cstop(1,1) = 0.d0
       cstop(1,2) = 0.d0
       cstop(2,1) = 0.d0
       cstop(2,2) = 0.d0

       higgsterm(1,1) = 0.d0
       higgsterm(1,2) = 0.d0
       higgsterm(2,1) = 0.d0
       higgsterm(2,2) = 0.d0


!--------------------------------------------------------------
      gtterm(1,1) = (4.d0 * (g3**2.d0)/3.d0) *
     $     (2.d0*ggs + (costhetas)**2.d0 * (fms10 + a0ms1) + 
     $     (sinthetas)**2.d0 * (fms20 + a0ms2))

      gtterm(2,2) = (4.d0 * (g3**2.d0)/3.d0) *
     $     (2.d0*ggs + (sinthetas)**2.d0 * (fms10 + a0ms1) + 
     $     (costhetas)**2.d0 * (fms20 + a0ms2))

      gtterm(1,2) = (4.d0*(g3**2.d0)/3.d0) *
     $     (4.d0* M3t*ms * b0M3ms + 
     $     sinthetas*costhetas * (fms10 - a0ms1 - fms20 + a0ms2))


!-------------------------------------------------------
      
      csbtm(1, 1) = (ydRG(2,2)**2.d0) * ((sinthetas)**2.d0 * a0ms1 + 
     $     (costhetas)**2.d0 * a0ms2)

      csbtm(2, 2) = (ydRG(2,2)**2.d0) * ((costhetas)**2.d0 * a0ms1 + 
     $     (sinthetas)**2.d0 * a0ms2)

      csbtm(1, 2) = (ydRG(2,2)**2.d0) * costhetas*sinthetas * 3.d0 *
     $     (a0ms1 - a0ms2)

!------------------------------------------------------------------

      cstop(1, 1) = (yuRG(2,2)**2.d0) * ((sinthetac)**2.d0 * a0mc1 + 
     $     (costhetac)**2.d0 * a0mc2)

      cstop(2, 2) = (yuRG(2,2)**2.d0) * ((costhetac)**2.d0 * a0mc1 + 
     $     (sinthetac)**2.d0 * a0mc2)

!--------------------------------------------------------------------

      higgsterm(1,1) = 0.5d0 * ((ydRG(2,2)**2.d0) * dnd(1) -
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((ydRG(2,2)**2.d0) * dnd(2) - 
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((ydRG(2,2)**2.d0) * dnd(3) - 
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((ydRG(2,2)**2.d0) * dnd(4) - 
     $     ((g*g)*gdL*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA



      higgsterm(1,2) = 0.d0

      higgsterm(2,2) = 0.5d0 * ((ydRG(2,2)**2.d0) * dnd(1) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((ydRG(2,2)**2.d0) * dnd(2) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((ydRG(2,2)**2.d0) * dnd(3) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((ydRG(2,2)**2.d0) * dnd(4) - 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA



      higgsterm(1,1) = higgsterm(1,1) + 
     $     ((yuRG(2,2)**2.d0) * dnd(3) + ((g*g) * 
     $     ((gdL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(3))) * a0Mhpm +
     $     ((yuRG(2,2)**2.d0) * dnd(4) + ((g*g) * 
     $     ((gdL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(4))) * a0MW

      higgsterm(2,2) = higgsterm(2,2) + 
     $     ((ydRG(2,2)**2.d0) * dnu(3) +  
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(3)) * a0Mhpm +
     $     ((ydRG(2,2)**2.d0) * dnu(4) + 
     $     ((g*g)*gdR*0.5d0/(costhw**2.d0)) * cn(4)) * a0MW  !<----------------chkd

!--------------------------

      call b0(p,mHu0,ms1,q,b0h0mstr(1,1))
      call b0(p,mHu0,ms2,q,b0h0mstr(1,2))
      call b0(p,mh0,ms1,q,b0h0mstr(2,1))
      call b0(p,mh0,ms2,q,b0h0mstr(2,2))
      call b0(p,MZ,ms1,q,b0h0mstr(3,1))
      call b0(p,MZ,ms2,q,b0h0mstr(3,2))
      call b0(p,mA0,ms1,q,b0h0mstr(4,1))
      call b0(p,mA0,ms2,q,b0h0mstr(4,2))
      

      loophi: DO i = 1, 4       !<-------------neutral higgs terms
      loophj: DO j = 1, 2

      higgsterm(1,1) = higgsterm(1,1) + 
     $     (lHsstrLstr12(i,j)**2.d0) * b0h0mstr(i,j)

      higgsterm(1,2) = higgsterm(1,2) + 
     $     lHsstrLstr12(i,j) * lHsstrRsstr12(i,j) * b0h0mstr(i,j)

      higgsterm(2,2) = higgsterm(2,2) + 
     $     (lHsstrRsstr12(i,j)**2.d0) * b0h0mstr(i,j)
      

      ENDDO loophj
      ENDDO loophi
!------------------------------------------------------------------

      call b0(p,mc1,mHpm,q,b0hcmscharm(1,1))
      call b0(p,mc2,mHpm,q,b0hcmscharm(1,2))
      call b0(p,mc1,MW,q,b0hcmscharm(2,1))
      call b0(p,mc2,MW,q,b0hcmscharm(2,2))
      

      loophci: DO i = 1, 2       !<-------------charged higgs terms
      loophcj: DO j = 1, 2

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (lHcstrlsc12(i,j)**2.d0) * b0hcmscharm(i,j)     

      higgsterm(1, 2) = higgsterm(1,2) + 
     $     lHcstrlsc12(i, j) *lHcstrrsc12(i, j)*b0hcmscharm(i,j)

      higgsterm(2, 2) = higgsterm(2, 2) + 
     $     (lHcstrrsc12(i,j)**2.d0) * b0hcmscharm(i,j)
      
   

      ENDDO loophcj
      ENDDO loophci
!------------------------------------------------------------------------------

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (4.d0*(g*g)/costhw**2.d0)*(gdL*gdL)*a0MZ  + 
     $     2.d0*(g*g)*a0MW + 
     $     ((1.d0/9.d0)*g*g*sinsqthw) * ((costhetas**2.d0) * 
     $     fms10 + (sinthetas**2.d0) * fms20) +
     $     (g*gdL/costhw)**2.d0 * ((costhetas**2.d0) * fms1MZ +
     $     (sinthetas**2.d0) * fms2MZ) +
     $     (g*g)*0.5d0*((costhetac)**2.d0 * fmc1Mw + (sinthetac**2.d0) *
     $     fmc2Mw) +
     $     (g*g) * 0.25d0 * (1.d0 * ((costhetas**2.d0) * a0ms1 + 
     $     (sinthetas**2.d0) * a0ms2) + 2.d0*((costhetac**2.d0) *a0mc1 + 
     $     (sinthetac**2.d0) * a0mc2)) -
     $     (g*g) * 0.5d0 * (- 1.5d0 * a0md1 - 1.5d0*a0ms1 +
     $     1.5d0 * a0mu1 + 1.5d0*a0mc1 +
     $     0.5d0 * (a0snu1 + a0snu2 + a0snu3) -
     $     0.5d0 * (a0me1 + a0mmu1 + ((costhetatau**2.d0) * a0mtau1 +
     $     (sinthetatau**2.d0) * a0mtau2))) +
     $     (gp**2.d0) * 0.25d0 *(ydL**2.d0) * ((costhetas**2.d0)*a0ms1 +
     $     (sinthetas**2.d0)*a0ms2) +
     $     (gp**2.d0) * 0.25d0 * ydL * ((3.d0*ydL*(a0md1 + a0ms1 + 
     $     ((costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2))) +
     $     (3.d0 * ydR *(a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2)) +
     $     (3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t1)) +
     $     (3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (costhetat**2.d0)*a0t2 + (sinthetat**2.d0)*a0t1)) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))

     
!----------------------------------------------------------------------------------

      higgsterm(2, 2) = higgsterm(2,2) +
     $     (4.d0*(g*g)/costhw**2.d0) * (gdR**2.d0) * a0MZ + 
     $     ((1.d0/9.d0)*g*g*sinsqthw) * ((sinthetas)* fms10 + 
     $     (costhetas**2.d0)*fms20) +
     $     (g*gdR/costhw)**2.d0 * ((sinthetas**2.d0) * fms1MZ + 
     $     (costhetas**2.d0) * fms2MZ) +
     $     (gp**2.d0) * 0.25d0 * (ydR**2.d0) * 
     $     ((sinthetas**2.d0) * a0ms1 + (costhetas**2.d0) * a0ms2) +
     $     (gp**2.d0) * 0.25d0 * ydR *  
     $     (3.d0 * ydL * (a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2) +
     $     3.d0 * ydR * (a0md2 +  a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2) +
     $     3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t2) +
     $     3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 + (costhetat**2.d0)*a0t2) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))


      higgsterm(1,2) = higgsterm(1,2) + 
     $     (gp**2.d0) * ydL*ydR * sinthetas*costhetas *
     $     (a0ms1 - a0ms2) +
     $     ((1.d0/9.d0) * g*g * sinsqthw) * sinthetas*costhetas * 
     $     (fms10 - fms20) -
     $     ((g*g)/costhw**2.d0) * gdL*gdR * sinthetas*costhetas *
     $     (fms1Mz - fms2Mz)    !<-----------------------chked


!----------------------------------------------------------------------
C     Chargino term
!----------------------------------------------------------------------

      call funcg(p,mchargino(1),mc,q,gmchargino1mc)
      call funcg(p,mchargino(2),mc,q,gmchargino2mc)
      call b0(p,mchargino(1),mc,q,b0mchargino1mc)
      call b0(p,mchargino(2),mc,q,b0mchargino2mc)
      
      chargino(1, 1) =  fChstrLL(1)*gmchargino1mc -
     $     gChstrLL(1)*mchargino(1)*mc*b0mchargino1mc*2.d0 +
     $     fChstrLL(2)*gmchargino2mc -
     $     gChstrLL(2)*mchargino(2)*mc*b0mchargino2mc*2.d0

      chargino(1, 2) =  fChstrLR(1)*gmchargino1mc -
     $     gChstrLR(1)*mchargino(1)*mc*b0mchargino1mc*2.d0 +
     $     fChstrLR(2)*gmchargino2mc -
     $     gChstrLR(2)*mchargino(2)*mc*b0mchargino2mc*2.d0

      chargino(2, 2) =  fChstrRR(1)*gmchargino1mc - 
     $     gChstrRR(1)*mchargino(1)*mc*b0mchargino1mc*2.d0 +
     $     fChstrRR(2)*gmchargino2mc - 
     $     gChstrRR(2)*mchargino(2)*mc*b0mchargino2mc*2.d0

!---------------------------------------------------------------------
C     neutralino terms
!---------------------------------------------------------------------

      call b0(p,nmneut(1),ms,q,b0mneut1ms)
      call b0(p,nmneut(2),ms,q,b0mneut2ms)
      call b0(p,nmneut(3),ms,q,b0mneut3ms)
      call b0(p,nmneut(4),ms,q,b0mneut4ms)
      
      call funcg(p,nmneut(1),ms,q,gmneut1ms)
      call funcg(p,nmneut(2),ms,q,gmneut2ms)
      call funcg(p,nmneut(3),ms,q,gmneut3ms)
      call funcg(p,nmneut(4),ms,q,gmneut4ms)

      neutralino(1, 1) = 
     $     fChi0strsstrLL(1)*gmneut1ms - gChi0strsstrLL(1)*2.d0*
     $     mneut(1)*ms*b0mneut1ms +
     $     fChi0strsstrLL(2)*gmneut2ms - gChi0strsstrLL(2)*2.d0*
     $     mneut(2)*ms*b0mneut2ms +
     $     fChi0strsstrLL(3)*gmneut3ms - gChi0strsstrLL(3)*2.d0*
     $     mneut(3)*ms*b0mneut3ms +
     $     fChi0strsstrLL(4)*gmneut4ms - gChi0strsstrLL(4)*2.d0*
     $     mneut(4)*ms*b0mneut4ms


      neutralino(2, 2) = 
     $     fChi0strsstrRR(i)*gmneut1ms - gChi0strsstrRR(1)*2.d0*
     $     mneut(1)*ms*b0mneut1ms +
     $     fChi0strsstrRR(2)*gmneut2ms - gChi0strsstrRR(2)*2.d0*
     $     mneut(2)*ms*b0mneut2ms +
     $     fChi0strsstrRR(3)*gmneut3ms - gChi0strsstrRR(3)*2.d0*
     $     mneut(3)*ms*b0mneut3ms +
     $     fChi0strsstrRR(4)*gmneut4ms - gChi0strsstrRR(4)*2.d0*
     $     mneut(4)*ms*b0mneut4ms

      neutralino(1, 2) = 
     $     fChi0strsstrLR(1)*gmneut1ms - gChi0strsstrLR(1)*2.d0*
     $     mneut(1)*ms*b0mneut1ms +
     $     fChi0strsstrLR(2)*gmneut2ms - gChi0strsstrLR(2)*2.d0*
     $     mneut(2)*ms*b0mneut2ms +
     $     fChi0strsstrLR(3)*gmneut3ms - gChi0strsstrLR(3)*2.d0*
     $     mneut(3)*ms*b0mneut3ms +
     $     fChi0strsstrLR(4)*gmneut4ms - gChi0strsstrLR(4)*2.d0*
     $     mneut(4)*ms*b0mneut4ms 
      
!-------------------------------------------------------------------

      MSQD3(1,1) = mSQRG(2,2) + MS**2.d0 + gdL*MZ*MZ*dcos(2.d0*beta)
      MSQD3(1,2) = MS*((ADRG(2,2) - sgnmu*modmu*dtan(beta)))
      MSQD3(2,1) = MSQD3(1,2)
      MSQD3(2,2) = mSDRG(2,2) + MS**2.d0 + gdR*MZ*MZ*dcos(2.d0*beta)

!----------------------------------------------------------------------

      pisstr(1,1) = (cstop(1,1) + csbtm(1,1) + gtterm(1,1) +
     $     chargino(1,1) + neutralino(1,1)+ higgsterm(1,1))/
     $     (16.d0*pi*pi)

      pisstr(1,2) = (cstop(1,2) + csbtm(1,2) + gtterm(1,2)  +
     $     chargino(1,2)+ neutralino(1,2) + higgsterm(1,2))/
     $     (16.d0*pi*pi)

      pisstr(2,2) = (cstop(2,2) + csbtm(2,2) + gtterm(2,2) +
     $     chargino(2,2)+ neutralino(2,2) +  higgsterm(2,2))/
     $     (16.d0*pi*pi)
      
      pisstr(2,1) = pisstr(1,2)


      delpisst(1,1) = pisstr(1,1)
      delpisst(1,2) = pisstr(1,2)
      delpisst(2,1) = pisstr(2,1)
      delpisst(2,2) = pisstr(2,2)

!--------------------------------------------------------------------


      pisstr(1,1) = MSQD3(1,1) - pisstr(1,1)
      pisstr(1,2) = MSQD3(1,2) - pisstr(1,2)
      pisstr(2,2) = MSQD3(2,2) - pisstr(2,2)
      pisstr(2,1) = pisstr(1,2) 

!----------------------------------------------------------------------

C     find the singular values and the diagonalising matrices

      info  = 10
      AOK = 0
      
C     call dsyev('V','U',2,pisstr,2,SSTeg,work,lwork,info)
      
      if(info.eq.0) then
         AOK = AOK + 1
      endif

      Call CEigensystem(2,pisstr,2,SSTeg,sstmix,2,1)
      
      RETURN

      END SUBROUTINE pisstrange

C==================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C     pischarm is checked on 30/05/2010 @ 11:30.
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      SUBROUTINE pischarm(p,q,g,gp,g3,tanbeta,mSQRG,
     $     mSURG,yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3t,delpisc,SCeg)

      IMPLICIT NONE 

      INTEGER i,j,info,AOK,lwork
      parameter(lwork=35)
      DOUBLE PRECISION work(lwork)
      DOUBLE PRECISION tanbeta,mSQRG(3,3),mSURG(3,3)
      DOUBLE PRECISION AURG(3,3),ADRG(3,3)
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yuRG(3,3),ydRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR,MSQU3(2,2)      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3),nmneut(4)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g,g3,M3t

      double precision vev1,vev2,beta
 
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION costhetac,sinthetac,costhetas,sinthetas


      DOUBLE PRECISION thetau,thetad,thetae

      DOUBLE PRECISION pisc(2,2)
 
      DOUBLE PRECISION higgsterm(2,2)
      DOUBLE PRECISION dnu(4),dnd(4),cn(4)

      DOUBLE PRECISION lscLscLR(4,2),lscLsc12(4,2)
      DOUBLE PRECISION lscRscLR(4,2),lscRsc12(4,2)
      
      DOUBLE PRECISION aPsi0cscr(4), bPsi0cscr(4), aPsi0cscl(4)
      DOUBLE PRECISION bPsi0cscl(4) 
      DOUBLE PRECISION aChi0cscl(4), bChi0cscl(4), aChi0cscr(4)
      DOUBLE PRECISION bChi0cscr(4)

      DOUBLE PRECISION gChi0cscLL(4), fChi0cscLL(4)
      DOUBLE PRECISION gChi0cscLR(4), fChi0cscLR(4)
      DOUBLE PRECISION gChi0cscRR(4), fChi0cscRR(4)

      DOUBLE PRECISION bPsicstrscl(2), bPsicstrscr(2), aPsicstrscl(2)
      DOUBLE PRECISION aPsicstrscr(2)
      DOUBLE PRECISION aChicstrscr(2), aChicstrscl(2)
      DOUBLE PRECISION bChicstrscl(2),bChicstrscr(2)
      DOUBLE PRECISION gtterm(2,2), cstop(2,2), csbtm(2,2) 
      DOUBLE PRECISION chargino(2,2), neutralino(2,2)
      DOUBLE PRECISION fChstrscLL(2), gChstrscLL(2) 
      DOUBLE PRECISION fChstrscLR(2), gChstrscLR(2) 
      DOUBLE PRECISION fChstrscRR(2), gChstrscRR(2)
 
      DOUBLE PRECISION intm1(2),intm2(2)

      DOUBLE PRECISION lHcscrsstlr(2,2),lHcscrsstr12(2,2)
      DOUBLE PRECISION lHcsclsst12(2,2)
      DOUBLE PRECISION lHcsclsstlr(2,2) 

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION lHscLsc12(4,2),lHscRsc12(4,2)
      data lHscLsc12/ 8 * 0.d0/, lHscRsc12/ 8 * 0.d0/

      
      DOUBLE PRECISION ralpha(2,2)

      DOUBLE PRECISION ggc,fc10,fc20,a0c1,a0c2,b0M3mc,a0mb1,a0mb2
      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz,mc1,mc2
      DOUBLE PRECISION b0h0msch(4,2),b0hcmstr(2,2),a0t1,a0t2

      DOUBLE PRECISION a0md1,a0md2,a0ms1,a0ms2,a0mtau1,a0mtau2,a0snu2
      DOUBLE PRECISION a0me1,a0me2,a0mmu1,a0mmu2,a0snu1,a0snu3
      DOUBLE PRECISION a0mc1,a0mc2,a0mu1,a0mu2
      
      DOUBLE PRECISION fms1MW,fms2MW,fmc1MZ,fmc2MZ
      DOUBLE PRECISION gmchargino1ms,gmchargino2ms
      DOUBLE PRECISION b0mchargino1ms,b0mchargino2ms
      DOUBLE PRECISION b0mneut1mc,b0mneut2mc,b0mneut3mc,b0mneut4mc
      DOUBLE PRECISION gmneut1mc,gmneut2mc,gmneut3mc,gmneut4mc

        
      DOUBLE PRECISION SCeg(2),scmix(2,2)

      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2
      DOUBLE PRECISION ms1,ms2, delpisc(2,2)
      
      double precision rthetac(2,2),rthetas(2,2)

      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae


      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d
  

      include 'stdinputs.h'

C---------------------------------------------------------------------

      pi = 4.d0 * datan(1.d0)
      MZ = MZpole
      MW = MWpole

      SCeg(1) = 0.d0
      SCeg(2) = 0.d0


      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw = dsqrt(sinsqthw)
      thw = dasin(sinthw)
      costhw = dcos(thw)
      cos2thw = dcos(2.d0*thw)
 
 
      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))
      
      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))

       
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

!----------------------------------     

      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)

      costhetac = dcos(thetac)
      sinthetac = dsin(thetac)

      costhetas = dcos(thetas)
      sinthetas = dsin(thetas)
     
!----------------------------------

      mt1 = mtL
      mt2 = mtR

      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR

      mc1 = mcL
      mc2 = mcR

      ms1 = msL
      ms2 = msR


C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))
      
      beta     = datan(tanbeta)
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha  = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta  = dcos(datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)

!----------------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0

!----------------------------------------------------------
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

!---------------------------------------------------------------

      pisc(1,1) = 0.d0
      pisc(2,1) = 0.d0
      pisc(1,2) = 0.d0
      pisc(2,2) = 0.d0
C---------------------------------------------------------------

      lscLscLR(1,1) = g*MZ*guL*cosbeta/costhw
      lscLscLR(1,2) = -yuRG(2,2)*sgnmu*modmu/dsqrt(2.d0)
      
      lscLscLR(2,1) = (- g*MZ*guL*sinbeta/costhw) + 
     $     (yuRG(2,2)**2.d0) * vev2
      lscLscLR(2,2) = (yuRG(2,2) * AURG(2,2))/dsqrt(2.d0)
      
      lscLscLR(3,1) = 0.d0
      lscLscLR(3,2) = ((-sgnmu*modmu*cosbeta*yuRG(2,2) + 
     $     yuRG(2,2) * AURG(2,2) * sinbeta))/dsqrt(2.d0)
      
      lscLscLR(4,1) = 0.d0
      lscLscLR(4,2) = -((-sgnmu*modmu*sinbeta*yuRG(2,2) - 
     $     yuRG(2,2) * AURG(2,2) * cosbeta))/dsqrt(2.d0)
      

      lscRscLR(1,1) = lscLscLR(1,2)
      lscRscLR(1,2) = g*MZ*guR*cosbeta/costhw

      lscRscLR(2,1) = lscLscLR(2,2)
      lscRscLR(2,2) = -g*MZ*guR*sinbeta/costhw + 
     $     (yuRG(2,2)**2.d0) * vev2

      lscRscLR(3,1) = -lscLscLR(3,2)
      lscRscLR(3,2) = 0.d0
      
      lscRscLR(4,1) = -lscLscLR(4,2)
      lscRscLR(4,2) = 0.d0

!----------------------------------------------
      
      loopmixst: DO i = 1, 4

      intm1(1) = lscLscLR(i,1)
      intm1(2) = lscLscLR(i,2)

      call rmat2d(thetac,rthetac)
      
      intm2(1) = rthetac(1,1)*intm1(1) + rthetac(1,2)*intm1(2)
      intm2(2) = rthetac(2,1)*intm1(1) + rthetac(2,2)*intm1(2)

      lscLsc12(i,1) = intm2(1)
      lscLsc12(i,2) = intm2(2)

      intm1(1) = lscRscLR(i,1)
      intm1(2) = lscRscLR(i,2)

      intm2(1) = rthetac(1,1)*intm1(1) + rthetac(1,2)*intm1(2)
      intm2(2) = rthetac(2,1)*intm1(1) + rthetac(2,2)*intm1(2)

      lscRsc12(i,1) = intm2(1)
      lscRsc12(i,2) = intm2(2)

      ENDDO loopmixst

!--------------------------------------------------------------------------
C     Mixing CP-even Higgs 
!--------------------------------------------------------------------------

      loopmixcpeven: DO i = 1, 2
      
      intm1(1) = lscLsc12(1,i)
      intm1(2) = lscLsc12(2,i)

      call rmat2d(alpha,ralpha)
      
      intm2(1) = ralpha(1,1)*intm1(1) + ralpha(1,2)*intm1(2)
      intm2(2) = ralpha(2,1)*intm1(1) + ralpha(2,2)*intm1(2)

      lHscLsc12(1,i) = intm2(1)
      lHscLsc12(2,i) = intm2(2)

      intm1(1) = lscRsc12(1,i)
      intm1(2) = lscRsc12(2,i)

      intm2(1) = ralpha(1,1)*intm1(1) + ralpha(1,2)*intm1(2)
      intm2(2) = ralpha(2,1)*intm1(1) + ralpha(2,2)*intm1(2)

      lHscRsc12(1, i) = intm2(1)
      lHscRsc12(2, i) = intm2(2)
      
      ENDDO loopmixcpeven
!----------------------------------------------------------------------
C     Feynman rules for charged higgs
!----------------------------------------------------------------------
!     <--------- (H+ G+, L R) basis

      lHcsclsstlr(1, 1) = (g*MW*sin2beta - yuRG(2,2)**2.d0 * vev2 *
     $     cosbeta - ydRG(2,2)**2.d0 * vev1 * sinbeta)/dsqrt(2.d0)

      lHcsclsstlr(1,2) = (-sgnmu*modmu*ydRG(2,2)*cosbeta -
     $     ADRG(2,2)*ydRG(2,2)*sinbeta)

      lHcsclsstlr(2,1) = (-g*MW*cos2beta - yuRG(2,2)**2.d0 * vev2 *
     $     sinbeta + ydRG(2,2)**2.d0 * vev1 * cosbeta)/dsqrt(2.d0)

      lHcsclsstlr(2,2) = ydRG(2,2) * (-sgnmu*modmu*sinbeta + 
     $     ADRG(2,2)*cosbeta)

!-------------------------------------------------------------------------
      
      intm1(1) = lHcsclsstlr(1, 1)
      intm1(2) = lHcsclsstlr(1, 2)

      call rmat2d(thetas,rthetas)
      
      intm2(1) = rthetas(1,1)*intm1(1) + rthetas(1,2)*intm1(2)
      intm2(2) = rthetas(2,1)*intm1(1) + rthetas(2,2)*intm1(2)

      
      lHcsclsst12(1, 1) = intm2(1)
      lHcsclsst12(1, 2) = intm2(2)

      intm1(1) = lHcsclsstlr(2, 1)
      intm1(2) = lHcsclsstlr(2, 2)


      call rmat2d(thetas,rthetas)
      
      intm2(1) = rthetas(1,1)*intm1(1) + rthetas(1,2)*intm1(2)
      intm2(2) = rthetas(2,1)*intm1(1) + rthetas(2,2)*intm1(2)

      lHcsclsst12(2, 1) = intm2(1)
      lHcsclsst12(2, 2) = intm2(2)

!---------------------------------------------------------------------------
!     <------------ (H+ G+, L R) basis

      lHcscrsstlr(1,1) = yuRG(2,2)*(-sgnmu*modmu*sinbeta - 
     $     AURG(2,2)*cosbeta)

      lHcscrsstlr(1,2) = yuRG(2,2)*ydRG(2,2)*(- vev1*cosbeta - 
     $     vev2*sinbeta)/dsqrt(2.d0)

      lHcscrsstlr(2,1) = -yuRG(2,2)*(-sgnmu*modmu*cosbeta + 
     $     AURG(2,2)*sinbeta)  

      lHcscrsstlr(2,2) = 0.d0

!-------------------------------------------------------------------------

      intm1(1) = lHcscrsstlr(1,1)
      intm1(2) = lHcscrsstlr(1,2)

      call rmat2d(thetas,rthetas)
      
      intm2(1) = rthetas(1,1)*intm1(1) + rthetas(1,2)*intm1(2)
      intm2(2) = rthetas(2,1)*intm1(1) + rthetas(2,2)*intm1(2)

      lHcscrsstr12(1,1) = intm2(1)
      lHcscrsstr12(1,2) = intm2(2)

      intm1(1) = lHcscrsstlr(2,1)
      intm1(2) = lHcscrsstlr(2,2)

      call rmat2d(thetas,rthetas)
      
      intm2(1) = rthetas(1,1)*intm1(1) + rthetas(1,2)*intm1(2)
      intm2(2) = rthetas(2,1)*intm1(1) + rthetas(2,2)*intm1(2)
      
      lHcscrsstr12(2, 1) = intm2(1)
      lHcscrsstr12(2, 2) = intm2(2)

!------------------------------------------------------------------------
C     Feynman rules - neutralino
!------------------------------------------------------------------------


      aPsi0cscr(1) = -4.d0*gp/(3.d0*dsqrt(2.d0))
      aPsi0cscr(2) = 0.d0
      aPsi0cscr(3) = 0.d0
      aPsi0cscr(4) = 0.d0

      bPsi0cscl(1) = gp/(3.d0*dsqrt(2.d0))
      bPsi0cscl(2) = g/dsqrt(2.d0)
      bPsi0cscl(3) = 0.d0
      bPsi0cscl(4) = 0.d0

      aPsi0cscl(1) = 0.d0
      aPsi0cscl(2) = 0.d0
      aPsi0cscl(3) = 0.d0 
      aPsi0cscl(4) = yuRG(2,2)

      bPsi0cscr(1) = 0.d0
      bPsi0cscr(2) = 0.d0
      bPsi0cscr(3) = 0.d0
      bPsi0cscr(4) = yuRG(2,2)


      aChi0cscl(1) = ON(1,1)*aPsi0cscl(1) + ON(1,2)*aPsi0cscl(2) +
     $                 ON(1,3)*aPsi0cscl(1) + ON(1,4)*aPsi0cscl(4) 
      aChi0cscl(2) = ON(2,1)*aPsi0cscl(1) + ON(2,2)*aPsi0cscl(2) +
     $                 ON(2,3)*aPsi0cscl(1) + ON(2,4)*aPsi0cscl(4) 
      aChi0cscl(3) = ON(3,1)*aPsi0cscl(1) + ON(3,2)*aPsi0cscl(2) +
     $                 ON(3,3)*aPsi0cscl(1) + ON(3,4)*aPsi0cscl(4) 
      aChi0cscl(4) = ON(4,1)*aPsi0cscl(1) + ON(4,2)*aPsi0cscl(2) +
     $                 ON(4,3)*aPsi0cscl(1) + ON(4,4)*aPsi0cscl(4) 


      bChi0cscl(1) = ON(1,1)*bPsi0cscl(1) + ON(1,2)*bPsi0cscl(2) +
     $                 ON(1,3)*bPsi0cscl(1) + ON(1,4)*bPsi0cscl(4) 
      bChi0cscl(2) = ON(2,1)*bPsi0cscl(1) + ON(2,2)*bPsi0cscl(2) +
     $                 ON(2,3)*bPsi0cscl(1) + ON(2,4)*bPsi0cscl(4) 
      bChi0cscl(3) = ON(3,1)*bPsi0cscl(1) + ON(3,2)*bPsi0cscl(2) +
     $                 ON(3,3)*bPsi0cscl(1) + ON(3,4)*bPsi0cscl(4) 
      bChi0cscl(4) = ON(4,1)*bPsi0cscl(1) + ON(4,2)*bPsi0cscl(2) +
     $                 ON(4,3)*bPsi0cscl(1) + ON(4,4)*bPsi0cscl(4) 




      aChi0cscr(1) = ON(1,1)*aPsi0cscr(1) + ON(1,2)*aPsi0cscr(2) +
     $                 ON(1,3)*aPsi0cscr(3) + ON(1,4)*aPsi0cscr(4)
      aChi0cscr(2) = ON(2,1)*aPsi0cscr(1) + ON(2,2)*aPsi0cscr(2) +
     $                 ON(2,3)*aPsi0cscr(3) + ON(2,4)*aPsi0cscr(4)
      aChi0cscr(3) = ON(3,1)*aPsi0cscr(1) + ON(3,2)*aPsi0cscr(2) +
     $                 ON(3,3)*aPsi0cscr(3) + ON(3,4)*aPsi0cscr(4)
      aChi0cscr(4) = ON(4,1)*aPsi0cscr(1) + ON(4,2)*aPsi0cscr(2) +
     $                 ON(4,3)*aPsi0cscr(3) + ON(4,4)*aPsi0cscr(4)


      bChi0cscr(1) = ON(1,1)*bPsi0cscr(1) + ON(1,2)*bPsi0cscr(2) +
     $                 ON(1,3)*bPsi0cscr(3) + ON(1,4)*bPsi0cscr(4) 
      bChi0cscr(2) = ON(2,1)*bPsi0cscr(1) + ON(2,2)*bPsi0cscr(2) +
     $                 ON(2,3)*bPsi0cscr(3) + ON(2,4)*bPsi0cscr(4) 
      bChi0cscr(3) = ON(3,1)*bPsi0cscr(1) + ON(3,2)*bPsi0cscr(2) +
     $                 ON(3,3)*bPsi0cscr(3) + ON(3,4)*bPsi0cscr(4) 
      bChi0cscr(4) = ON(4,1)*bPsi0cscr(1) + ON(4,2)*bPsi0cscr(2) +
     $                 ON(4,3)*bPsi0cscr(3) + ON(4,4)*bPsi0cscr(4) 


      loopchtllrr: DO i = 1, 4

      fChi0cscLL(i) = (aChi0cscl(i)*aChi0cscl(i) + 
     $     bChi0cscl(i)*bChi0cscl(i))
      gChi0cscLL(i) = (bChi0cscl(i)*aChi0cscl(i) + 
     $     bChi0cscl(i)*aChi0cscl(i))
      fChi0cscRR(i) = (aChi0cscr(i)*aChi0cscr(i) + 
     $     bChi0cscr(i)*bChi0cscr(i))
      gChi0cscRR(i) = (bChi0cscr(i)*aChi0cscr(i) + 
     $     bChi0cscr(i)*aChi0cscr(i))
      fChi0cscLR(i) = (aChi0cscr(i)*aChi0cscl(i) + 
     $     bChi0cscr(i)*bChi0cscl(i))
      gChi0cscLR(i) = (bChi0cscl(i)*aChi0cscr(i) + 
     $     bChi0cscr(i)*aChi0cscl(i))
      
      ENDDO loopchtllrr

!---------------------------------------------------------------------------
c     Feynman Rules for chargino
      

      aPsicstrscl(1) = g
      aPsicstrscl(2) = 0.d0
      aPsicstrscr(1) = 0.d0
      aPsicstrscr(2) = -yuRG(2,2)
      bPsicstrscl(1) = 0.d0
      bPsicstrscl(2) = -ydRG(2,2)
      bPsicstrscr(1) = 0.d0
      bPsicstrscr(2) = 0.d0
      

      aChicstrscl(1) = OCR(1,1)*aPsicstrscl(1) + OCR(1,2)*aPsicstrscl(2)
      aChicstrscl(2) = OCR(2,1)*aPsicstrscl(1) + OCR(2,2)*aPsicstrscl(2)

      bChicstrscl(1) = OCL(1,1)*bPsicstrscl(1) + OCL(1,2)*bPsicstrscl(2)
      bChicstrscl(2) = OCL(2,1)*bPsicstrscl(1) + OCL(2,2)*bPsicstrscl(2)

      aChicstrscr(1) = OCR(1,1)*aPsicstrscr(1) + OCR(1,2)*aPsicstrscr(2)
      aChicstrscr(2) = OCR(2,1)*aPsicstrscr(1) + OCR(2,2)*aPsicstrscr(2)

      bChicstrscr(1) = OCL(1,1)*bPsicstrscr(1) + OCL(1,2)*bPsicstrscr(2)
      bChicstrscr(2) = OCL(2,1)*bPsicstrscr(1) + OCL(2,2)*bPsicstrscr(2)

      
      loopchst: DO i = 1, 2

      fChstrscLL(i) = (aChicstrscl(i)*aChicstrscl(i) +
     $     bChicstrscl(i)*bChicstrscl(i))
      gChstrscLL(i) = (bChicstrscl(i)*aChicstrscl(i) +
     $     aChicstrscl(i)*bChicstrscl(i))
      fChstrscLR(i) = (aChicstrscl(i)*aChicstrscr(i) +
     $     bChicstrscl(i)*bChicstrscr(i))
      gChstrscLR(i) = (bChicstrscl(i)*aChicstrscr(i) +
     $     aChicstrscl(i)*bChicstrscr(i))
      fChstrscRR(i) = (aChicstrscr(i)*aChicstrscr(i) +
     $     bChicstrscr(i)*bChicstrscr(i))
      gChstrscRR(i) = (bChicstrscr(i)*aChicstrscr(i) +
     $     aChicstrscr(i)*bChicstrscr(i))

      ENDDO loopchst

!--------------------------------------------------------------------
C     Corrections Begin
!------------------------------------------------------------------

      call funcg(p,M3t,mc,q,ggc)
      
      call f(p,mc1,0.d0,q,fc10)
      call f(p,mc2,0.d0,q,fc20)

      call a0(muR,q,a0mu2)
      call a0(muL,q,a0mu1)
      call a0(mc2,q,a0c2)
      call a0(mc1,q,a0c1)
      call a0(mc2,q,a0mc2)
      call a0(mc1,q,a0mc1)
      call a0(mt2,q,a0t2)
      call a0(mt1,q,a0t1)

      call a0(ms2,q,a0ms2)
      call a0(ms1,q,a0ms1)
      call a0(mdR,q,a0md2)
      call a0(mdL,q,a0md1)
      call a0(mb2,q,a0mb2)
      call a0(mb1,q,a0mb1)

      call a0(meR,q,a0me2)
      call a0(meL,q,a0me1)
      call a0(mmuR,q,a0mmu2)
      call a0(mmuL,q,a0mmu1)
      call a0(mtau2,q,a0mtau2)
      call a0(mtau1,q,a0mtau1)

      call a0(snu(1),q,a0snu1)
      call a0(snu(2),q,a0snu2)
      call a0(snu(3),q,a0snu3)


      call a0(mA0,q,a0mA)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)
      call a0(mHpm,q,a0mHpm)
      call a0(MW,q,a0Mw)
      call a0(MZ,q,a0Mz)

      call f(p,ms1,MW,q,fms1Mw)
      call f(p,ms2,MW,q,fms2Mw)
      call f(p,mc1,MZ,q,fmc1Mz)
      call f(p,mc2,MZ,q,fmc2Mz)


      call b0(p,M3t,mc,q,b0M3mc)
      

      initi:do i = 1,2
      initj:do j= 1,2

      gtterm(i,j) = 0.d0
      cstop(i,j) = 0.d0
      csbtm(i,j) = 0.d0
      higgsterm(i,j) = 0.d0

      enddo initj
      enddo initi
!--------------------------------------------------------------
      gtterm(1, 1) = 4.d0 * ((g3**2.d0)/3.d0) *
     $     (2.d0*ggc + (costhetac)**2.d0 * (fc10 + a0c1) + 
     $     (sinthetac)**2.d0 * (fc20 + a0c2))

      gtterm(2, 2) = 4.d0 * ((g3**2.d0)/3.d0) * 
     $     (2.d0*ggc + (sinthetac)**2.d0*(fc10 + a0c1) + 
     $     (costhetac)**2.d0*(fc20 + a0c2))

      gtterm(1, 2) = 4.d0 * ((g3**2.d0)/3.d0) *
     $     (4.d0 * M3t*mc * b0M3mc + 
     $     sinthetac*costhetac * (fc10 - a0c1 - fc20 + a0c2))

!---------------------------------------------------------------------

      cstop(1, 1) = (yuRG(2,2)**2.d0) * ((sinthetac)**2.d0 * a0c1 + 
     $     (costhetac)**2.d0 * a0c2)

      cstop(2, 2) = (yuRG(2,2)**2.d0) * ((costhetac)**2.d0 * a0c1 + 
     $     (sinthetac)**2.d0 * a0c2)

      cstop(1, 2) = (yuRG(2,2)**2.d0) * costhetac*sinthetac * 3.d0 *
     $     (a0c1 - a0c2)

!----------------------------------------------------------------------

      csbtm(1, 1) = (ydRG(2,2)**2.d0) * ((sinthetas)**2.d0 * a0ms1 + 
     $     (costhetas)**2.d0 * a0ms2)

      csbtm(2, 2) = (yuRG(2,2)**2.d0) * ((costhetas)**2.d0 * a0ms1 + 
     $     (sinthetas)**2.d0 * a0ms2)

!--------------------------------------------------------------------

      higgsterm(1,1) = 0.5d0 * ((yuRG(2,2)**2.d0) * dnu(1) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((yuRG(2,2)**2.d0) * dnu(2) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((yuRG(2,2)**2.d0) * dnu(3) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((yuRG(2,2)**2.d0) * dnu(4) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(4)) * a0mA


      higgsterm(2,2) = 0.5d0 * ((yuRG(2,2)**2.d0) * dnu(1) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((yuRG(2,2)**2.d0) * dnu(2) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((yuRG(2,2)**2.d0) * dnu(3) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((yuRG(2,2)**2.d0) * dnu(4) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(4)) * a0mA


      higgsterm(1,1) = higgsterm(1,1) + 
     $     ((ydRG(2,2)**2.d0) * dnu(3) + ((g*g) * 
     $     ((guL*0.5d0/(costhw**2.d0)) - 0.5d0) * cn(3))) * a0mhpm +
     $     ((ydRG(2,2)**2.d0) * dnu(4) + ((g*g) * 
     $     ((guL*0.5d0/(costhw**2.d0)) - 0.5d0) * cn(4))) * a0MW


      higgsterm(2,2) = higgsterm(2,2) + 
     $     ((yuRG(2,2)**2.d0) * dnd(3) + 
     $     ((g*g) * (guR * 0.5d0/(costhw**2.d0)) * cn(3))) * a0mhpm +
     $     ((yuRG(2,2)**2.d0) * dnd(4) + 
     $     ((g*g) * (guR * 0.5d0/(costhw**2.d0)) * cn(4))) * a0MW


!---------------------------------------------------------------

      call b0(p,mHu0,mc1,q,b0h0msch(1,1))
      call b0(p,mHu0,mc2,q,b0h0msch(1,2))
      call b0(p,mh0,mc1,q,b0h0msch(2,1))
      call b0(p,mh0,mc2,q,b0h0msch(2,2))
      call b0(p,MZ,mc1,q,b0h0msch(3,1))
      call b0(p,MZ,mc2,q,b0h0msch(3,2))
      call b0(p,mA0,mc1,q,b0h0msch(4,1))
      call b0(p,mA0,mc2,q,b0h0msch(4,2))
      

      loophi: DO i = 1, 4       !<-------------neutral higgs terms
      loophj: DO j = 1, 2

      higgsterm(1,1) = higgsterm(1,1) + 
     $     (lHscLsc12(i,j)**2.d0) * b0h0msch(i,j)

      higgsterm(1,2) = higgsterm(1,2) + 
     $     lHscLsc12(i,j) * lHscRsc12(i,j) * b0h0msch(i,j)

      higgsterm(2,2) = higgsterm(2,2) + 
     $     (lHscRsc12(i,j)**2.d0) * b0h0msch(i,j)
      
      ENDDO loophj
      ENDDO loophi

!-------------------------------------

      call b0(p,ms1,mHpm,q,b0hcmstr(1,1))
      call b0(p,ms2,mHpm,q,b0hcmstr(1,2))
      call b0(p,ms1,MW,q,b0hcmstr(2,1))
      call b0(p,ms2,MW,q,b0hcmstr(2,2))
      

      loophci: DO i = 1, 2      !<-------------charged higgs terms
      loophcj: DO j = 1, 2

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (lHcsclsst12(i,j)**2.d0) * b0hcmstr(i,j)        

      higgsterm(1, 2) = higgsterm(1,2) + 
     $     lHcsclsst12(i, j) * lHcscrsstr12(i, j) * b0hcmstr(i,j)

      higgsterm(2, 2) = higgsterm(2, 2) + 
     $     (lHcscrsstr12(i,j)**2.d0) * b0hcmstr(i,j)
      
      ENDDO loophcj
      ENDDO loophci

!----------------------------------------------


      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (4.d0*(g*g)/costhw**2.d0) * (guL*guL) * a0Mz  + 
     $     2.d0*(g*g) * a0Mw + 
     $     ((4.d0/9.d0)*g*g*sinsqthw)*((costhetac**2.d0) * fc10 + 
     $     (sinthetac**2.d0) * fc20) +
     $     (g*guL/costhw)**2.d0 * ((costhetac**2.d0) * fmc1MZ  +
     $     (sinthetac**2.d0) * fmc2MZ) +
     $     (g*g)*0.5d0*((costhetas**2.d0)*fms1Mw + (sinthetas**2.d0)*
     $     fms2Mw) + 
     $     (g*g)*0.25d0*((costhetac**2.d0)*a0c1 + (sinthetac**2.d0)*
     $     a0c2 + 2.d0 * ((costhetas**2.d0)*a0ms1 + (sinthetas**2.d0)*
     $     a0ms2)) +
     $     (g*g) * 0.5d0 * (1.5d0*a0mu1 + 1.5d0*a0mc1 +
     $     1.5d0*((costhetat**2.d0)*a0t1 + 
     $     (sinthetat**2.d0)*a0t2) -
     $     1.5d0*a0md1 - 1.5d0*a0ms1 -
     $     1.5d0*((costhetab**2.d0)*a0mb1 +
     $     (sinthetab**2.d0)*a0mb2) +
     $     0.5d0*(a0snu1 + a0snu2 + a0snu3) -
     $     0.5d0*(a0me1 + a0mmu1 +
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2)) +
     $     (gp**2.d0)*0.25d0*(yuL**2.d0) * ((costhetac**2.d0) * a0c1 +
     $     (sinthetac**2.d0) * a0c2) +
     $     (gp**2.d0)*0.25d0 * yuL * (3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 +(sinthetat**2.d0)*a0t2) +
     $     3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 +(costhetat**2.d0)*a0t2) +
     $     3.d0 * ydL * (a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0mb1 + (sinthetab**2.d0)*a0mb1) +
     $     3.d0 * ydR * (a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0mb1 + (costhetab**2.d0)*a0mb2) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 +(costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))


!----------------------------------

      higgsterm(2, 2) = higgsterm(2,2) +
     $     (4.d0*(g*g)/costhw**2.d0) * (guR**2.d0) * a0MZ + 
     $     ((4.d0/9.d0)*g*g*sinsqthw) * ((sinthetac)*fc10 + 
     $     (costhetac**2.d0)*fc20) +
     $     (g*guR/costhw)**2.d0 * ((sinthetac**2.d0)*fmc1MZ + 
     $     (costhetac**2.d0)*fmc2MZ) +
     $     (gp**2.d0) * 0.25d0 * (yuR**2.d0) *
     $     ((sinthetac**2.d0)*a0c1 + (costhetac**2.d0)*a0c2) +
     $     (gp**2.d0) * 0.25d0 * yuR * 
     $     (3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t2) +
     $     3.d0 * yuR * (a0mu2 +  a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 + (costhetat**2.d0)*a0t2) +
     $     3.d0 * ydL * (a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0mb1 + (sinthetab**2.d0)*a0mb2) +
     $     3.d0 * ydR * (a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0mb1 + (costhetab**2.d0)*a0mb2) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))
!----------------------------------------------------------------------

      higgsterm(1,2) = higgsterm(1,2) + 
     $     (gp**2.d0) * yuL*yuR * sinthetac*costhetac *
     $     (a0c1 - a0c2) +
     $     ((4.d0/9.d0) * g*g * sinsqthw) * sinthetac*costhetac * 
     $     (fc10 - fc20) -
     $     ((g*g)/costhw**2.d0) * guL*guR * sinthetac*costhetac *
     $     (fmc1Mz - fmc2Mz)

!----------------------------------------------------------------------
C     Chargino term
!----------------------------------------------------------------------

      chargino(1,1) = 0.d0
      chargino(1,2) = 0.d0
      chargino(2,1) = 0.d0
      chargino(2,2) = 0.d0

      call funcg(p,mchargino(1),ms,q,gmchargino1ms)
      call funcg(p,mchargino(2),ms,q,gmchargino2ms)
      call b0(p,mchargino(1),ms,q,b0mchargino1ms)
      call b0(p,mchargino(2),ms,q,b0mchargino2ms)
      
      chargino(1, 1) = chargino(1, 1) + fChstrscLL(1)*gmchargino1ms -
     $     gChstrscLL(1)*mchargino(1)*ms*b0mchargino1ms*2.d0 +
     $     fChstrscLL(2)*gmchargino2ms -
     $     gChstrscLL(2)*mchargino(2)*ms*b0mchargino2ms*2.d0

      chargino(1, 2) = chargino(1, 2) + fChstrscLR(1)*gmchargino1ms -
     $     gChstrscLR(1)*mchargino(1)*ms*b0mchargino1ms*2.d0 +
     $     fChstrscLR(2)*gmchargino2ms -
     $     gChstrscLR(2)*mchargino(2)*ms*b0mchargino2ms*2.d0

      chargino(2, 2) = chargino(2, 2) + fChstrscRR(1)*gmchargino1ms - 
     $     gChstrscRR(1)*mchargino(1)*ms*b0mchargino1ms*2.d0 +
     $     fChstrscRR(2)*gmchargino2ms - 
     $     gChstrscRR(2)*mchargino(2)*ms*b0mchargino2ms*2.d0

!---------------------------------------------------------------------
C     neutralino terms
!---------------------------------------------------------------------

      neutralino(1,1) = 0.d0
      neutralino(1,2) = 0.d0
      neutralino(2,1) = 0.d0
      neutralino(2,2) = 0.d0

      call b0(p,nmneut(1),mc,q,b0mneut1mc)
      call b0(p,nmneut(2),mc,q,b0mneut2mc)
      call b0(p,nmneut(3),mc,q,b0mneut3mc)
      call b0(p,nmneut(4),mc,q,b0mneut4mc)
      
      call funcg(p,nmneut(1),mc,q,gmneut1mc)
      call funcg(p,nmneut(2),mc,q,gmneut2mc)
      call funcg(p,nmneut(3),mc,q,gmneut3mc)
      call funcg(p,nmneut(4),mc,q,gmneut4mc)

      neutralino(1, 1) = neutralino(1, 1) +
     $     fChi0cscLL(1)*gmneut1mc - gChi0cscLL(1)*2.d0*
     $     mneut(1)*mc*b0mneut1mc +
     $     fChi0cscLL(2)*gmneut2mc - gChi0cscLL(2)*2.d0*
     $     mneut(2)*mc*b0mneut2mc +
     $     fChi0cscLL(3)*gmneut3mc - gChi0cscLL(3)*2.d0*
     $     mneut(3)*mc*b0mneut3mc +
     $     fChi0cscLL(4)*gmneut4mc - gChi0cscLL(4)*2.d0*
     $     mneut(4)*mc*b0mneut4mc


      neutralino(2, 2) = neutralino(2, 2) +
     $     fChi0cscRR(i)*gmneut1mc - gChi0cscRR(1)*2.d0*
     $     mneut(1)*mc*b0mneut1mc +
     $     fChi0cscRR(2)*gmneut2mc - gChi0cscRR(2)*2.d0*
     $     mneut(2)*mc*b0mneut2mc +
     $     fChi0cscRR(3)*gmneut3mc - gChi0cscRR(3)*2.d0*
     $     mneut(3)*mc*b0mneut3mc +
     $     fChi0cscRR(4)*gmneut4mc - gChi0cscRR(4)*2.d0*
     $     mneut(4)*mc*b0mneut4mc

      neutralino(1, 2) = neutralino(1, 2) +
     $     fChi0cscLR(1)*gmneut1mc - gChi0cscLR(1)*2.d0*
     $     mneut(1)*mc*b0mneut1mc +
     $     fChi0cscLR(2)*gmneut2mc - gChi0cscLR(2)*2.d0*
     $     mneut(2)*mc*b0mneut2mc +
     $     fChi0cscLR(3)*gmneut3mc - gChi0cscLR(3)*2.d0*
     $     mneut(3)*mc*b0mneut3mc +
     $     fChi0cscLR(4)*gmneut4mc - gChi0cscLR(4)*2.d0*
     $     mneut(4)*mc*b0mneut4mc 
      
!-------------------------------------------------------------------

      MSQU3(1,1) = mSQRG(2,2) + mc**2.d0 + guL*MZ*MZ*dcos(2.d0*beta)
      MSQU3(1,2) = mc * ((AURG(2,2) - sgnmu*modmu/dtan(beta)))
      MSQU3(2,1) = MSQU3(1,2)
      MSQU3(2,2) = mSURG(2,2) + mc**2.d0 + guR*MZ*MZ*dcos(2.d0*beta)

!----------------------------------------------------------------------
      
      pisc(1,1) = (cstop(1,1) + csbtm(1,1) + gtterm(1,1) +
     $     chargino(1,1) + neutralino(1,1) + higgsterm(1,1))/
     $     (16.d0*pi*pi)

      pisc(1,2) = (cstop(1,2) + csbtm(1,2) + gtterm(1,2) +
     $     chargino(1,2) + neutralino(1,2) + higgsterm(1,2))/
     $     (16.d0*pi*pi)

      pisc(2,2) = (cstop(2,2) + csbtm(2,2) + gtterm(2,2) +
     $     chargino(2,2) + neutralino(2,2) +  higgsterm(2,2))/
     $     (16.d0*pi*pi)
      
      pisc(2,1) = pisc(1,2)

      delpisc(1,1) = pisc(1,1)
      delpisc(1,2) = pisc(1,2)
      delpisc(2,1) = pisc(2,1)
      delpisc(2,2) = pisc(2,2)


      pisc(1,1) = MSQU3(1,1) - pisc(1,1)
      pisc(1,2) = MSQU3(1,2) - pisc(1,2)
      pisc(2,2) = MSQU3(2,2) - pisc(2,2)
      pisc(2,1) = pisc(1,2) 

!---------------------------------------------

C     find the singular values and the diagonalising matrices

      info  = 10
      AOK = 0
      
    !  call dsyev('V','U',2,pisc,2,SCeg,work,lwork,info)
      
      if(info.eq.0) then
         AOK = AOK + 1
      endif     
      
      Call CEigensystem(2,pisc,2,SCeg,scmix,2,1)

      RETURN

      END SUBROUTINE pischarm

C==============================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     
C     pisupq  is checked on 31/05/2010 @ 21:30.
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      SUBROUTINE pisupq(p,q,g,gp,g3,tanbeta,mSQRG,
     $     mSURG,yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3t,SUeg)
 
      IMPLICIT NONE 

      INTEGER i,j,info,AOK,lwork
      parameter(lwork=35)
      DOUBLE PRECISION work(lwork)
      DOUBLE PRECISION tanbeta,mSQRG(3,3),mSURG(3,3)
      DOUBLE PRECISION AURG(3,3),ADRG(3,3)
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yuRG(3,3),ydRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR,MSQU3(2,2)      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3),nmneut(4)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g,g3,M3t

      double precision vev1,vev2,beta
 
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetac,thetas,thetamu

      DOUBLE PRECISION thetau,thetad,thetae
      DOUBLE PRECISION costhetau,sinthetau,costhetad,sinthetad
  
      DOUBLE PRECISION pisup(2,2)
 
      DOUBLE PRECISION higgsterm(2,2)
      DOUBLE PRECISION dnu(4),dnd(4),cn(4)

      DOUBLE PRECISION lsupLsupLR(4,2),lsupLsup12(4,2)
      DOUBLE PRECISION lsupRsupLR(4,2),lsupRsup12(4,2)
      
      DOUBLE PRECISION aPsi0usupr(4), bPsi0usupr(4), aPsi0usupl(4)
      DOUBLE PRECISION bPsi0usupl(4) 
      DOUBLE PRECISION aChi0usupl(4), bChi0usupl(4), aChi0usupr(4)
      DOUBLE PRECISION bChi0usupr(4)

      DOUBLE PRECISION gChi0usupLL(4), fChi0usupLL(4)
      DOUBLE PRECISION gChi0usupLR(4), fChi0usupLR(4)
      DOUBLE PRECISION gChi0usupRR(4), fChi0usupRR(4)

      DOUBLE PRECISION bPsicdsupl(2), bPsicdsupr(2), aPsicdsupl(2)
      DOUBLE PRECISION aPsicdsupr(2)
      DOUBLE PRECISION aChicdsupr(2), aChicdsupl(2)
      DOUBLE PRECISION bChicdsupl(2),bChicdsupr(2)
      DOUBLE PRECISION gtterm(2,2), cstop(2,2), csbtm(2,2) 
      DOUBLE PRECISION  chargino(2,2), neutralino(2,2)
      DOUBLE PRECISION fChdsupLL(2), gChdsupLL(2) 
      DOUBLE PRECISION fChdsupLR(2), gChdsupLR(2) 
      DOUBLE PRECISION fChdsupRR(2), gChdsupRR(2)
 
      DOUBLE PRECISION intm1(2), intm2(2)

      DOUBLE PRECISION lHcsuprsdnlr(2,2),lHcsuprsdn12(2,2)
      DOUBLE PRECISION lHcsuplsdn12(2,2)
      DOUBLE PRECISION lHcsuplsdnlr(2,2) 

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION lHsupLsup12(4,2),lHsupRsup12(4,2)
      data lHsupLsup12/ 8 * 0.d0/, lHsupRsup12/ 8 * 0.d0/
      
      DOUBLE PRECISION rthetau(2,2),rthetad(2,2),ralpha(2,2)

!     double precision rthetac(2,2),rthetas(2,2)

      DOUBLE PRECISION ggu,fu10,fu20,a0u1,a0u2,b0M3mu,a0mb1,a0mb2
      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz,a0c1,a0c2
      DOUBLE PRECISION b0h0msu(4,2),b0hcmsd(2,2),a0t1,a0t2,mc1,mc2

      DOUBLE PRECISION a0md1,a0md2,a0ms1,a0ms2,a0mtau1,a0mtau2,a0snu2
      DOUBLE PRECISION a0me1,a0me2,a0mmu1,a0mmu2,a0snu1,a0snu3
      DOUBLE PRECISION a0mc1,a0mc2,a0mu1,a0mu2
      
      DOUBLE PRECISION fmd1MW,fmd2MW,fmu1MZ,fmu2MZ
      DOUBLE PRECISION gmchargino1ms,gmchargino2ms
      DOUBLE PRECISION b0mchargino1ms,b0mchargino2ms
      DOUBLE PRECISION b0mneut1mu,b0mneut2mu,b0mneut3mu,b0mneut4mu
      DOUBLE PRECISION gmneut1mu,gmneut2mu,gmneut3mu,gmneut4mu

        
      DOUBLE PRECISION SUeg(2),sumix(2,2)

      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2
      DOUBLE PRECISION ms1,ms2,mu1,mu2,mu,md1,md2
  
      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae


      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d
  

      include 'stdinputs.h'


C-----------------------------

      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      MZ = MZpole

      mu = mUQ

      SUeg(1) = 0.d0
      SUeg(2) = 0.d0

      
      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw   = dsqrt(sinsqthw)
      thw      = dasin(sinthw)
      costhw   = dcos(thw)
      cos2thw  = dcos(2.d0*thw)


      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))
      
      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))

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

!--------------------------------------------------
      
      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)

      costhetau = dcos(thetau)
      sinthetau = dsin(thetau)

      costhetad = dcos(thetad)
      sinthetad = dsin(thetad)

     
!-----------------------------------------------

      mt1 = mtL
      mt2 = mtR

      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR

      mc1 = mcL
      mc2 = mcR

      ms1 = msL
      ms2 = msR

      mu1 = muL
      mu2 = muR

      md1 = mdL
      md2 = mdR


C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta     = datan(tanbeta)
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta  = dcos(datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)

!-----------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0

!----------------------------------------------------------
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

!---------------------------------------------------------------    
      pisup(1,1) = 0.d0
      pisup(2,1) = 0.d0
      pisup(1,2) = 0.d0
      pisup(2,2) = 0.d0
C---------------------------------------------------------------    

      lsupLsupLR(1,1) = g*MZ*guL*cosbeta/costhw
      lsupLsupLR(1,2) = -yuRG(1,1)*sgnmu*modmu/dsqrt(2.d0)
      
      lsupLsupLR(2,1) = (- g*MZ*guL*sinbeta/costhw) + 
     $     (yuRG(1,1)**2.d0) * vev2
      lsupLsupLR(2,2) = (yuRG(1,1) * AURG(1,1))/dsqrt(2.d0)
      
      lsupLsupLR(3,1) = 0.d0
      lsupLsupLR(3,2) = ((-sgnmu*modmu*cosbeta*yuRG(1,1) + 
     $     yuRG(1,1) * AURG(1,1) * sinbeta))/dsqrt(2.d0)
      
      lsupLsupLR(4,1) = 0.d0
      lsupLsupLR(4,2) = -((-sgnmu*modmu*sinbeta*yuRG(1,1) - 
     $     yuRG(1,1) * AURG(1,1) * cosbeta))/dsqrt(2.d0)
      

      lsupRsupLR(1,1) = lsupLsupLR(1,2)
      lsupRsupLR(1,2) = g*MZ*guR*cosbeta/costhw

      lsupRsupLR(2,1) = lsupLsupLR(2,2)
      lsupRsupLR(2,2) = -g*MZ*guR*sinbeta/costhw + 
     $     (yuRG(1,1)**2.d0) * vev2

      lsupRsupLR(3,1) = -lsupLsupLR(3,2)
      lsupRsupLR(3,2) = 0.d0
      
      lsupRsupLR(4,1) = -lsupLsupLR(4,2)
      lsupRsupLR(4,2) = 0.d0

!--------------------------------------------------------------------------------  
  
      loopmixst: DO i = 1, 4

      intm1(1) = lsupLsupLR(i,1)
      intm1(2) = lsupLsupLR(i,2)

      call rmat2d(thetau,rthetau)
            
      intm2(1) = rthetau(1,1)*intm1(1) + rthetau(1,2)*intm1(2)
      intm2(2) = rthetau(2,1)*intm1(1) + rthetau(2,2)*intm1(2)

      lsupLsup12(i,1) = intm2(1)
      lsupLsup12(i,2) = intm2(2)

      intm1(1) = lsupRsupLR(i,1)
      intm1(2) = lsupRsupLR(i,2)

      intm2(1) = rthetau(1,1)*intm1(1) + rthetau(1,2)*intm1(2)
      intm2(2) = rthetau(2,1)*intm1(1) + rthetau(2,2)*intm1(2)

      lsupRsup12(i,1) = intm2(1)
      lsupRsup12(i,2) = intm2(2)

      ENDDO loopmixst

!--------------------------------------------------------------------------
C       Mixing CP-even Higgs 
!--------------------------------------------------------------------------

      loopmixcpeven: DO i = 1, 2
      
      intm1(1) = lsupLsup12(1,i)
      intm1(2) = lsupLsup12(2,i)

      call rmat2d(alpha,ralpha)
            
      intm2(1) = ralpha(1,1)*intm1(1) + ralpha(1,2)*intm1(2)
      intm2(2) = ralpha(2,1)*intm1(1) + ralpha(2,2)*intm1(2)

      lHsupLsup12(1,i) = intm2(1)
      lHsupLsup12(2,i) = intm2(2)

      intm1(1) = lsupRsup12(1,i)
      intm1(2) = lsupRsup12(2,i)

      intm2(1) = ralpha(1,1)*intm1(1) + ralpha(1,2)*intm1(2)
      intm2(2) = ralpha(2,1)*intm1(1) + ralpha(2,2)*intm1(2)

      lHsupRsup12(1, i) = intm2(1)
      lHsupRsup12(2, i) = intm2(2)
     
      ENDDO loopmixcpeven

!----------------------------------------------------------------------
C    Feynman rules - charged higgs
!----------------------------------------------------------------------
!     <--------- (H+ G+, L R) basis

      lHcsuplsdnlr(1, 1) = (g*MW*sin2beta - yuRG(1,1)**2.d0 * vev2 *
     $     cosbeta - ydRG(1,1)**2.d0 * vev1 * sinbeta)/dsqrt(2.d0)

      lHcsuplsdnlr(1,2) = (-sgnmu*modmu*ydRG(1,1)*cosbeta -
     $     ADRG(1,1)*ydRG(1,1)*sinbeta)

      lHcsuplsdnlr(2,1) = (-g*MW*cos2beta - yuRG(1,1)**2.d0 * vev2 *
     $     sinbeta + ydRG(1,1)**2.d0 * vev1 * cosbeta)/dsqrt(2.d0)

      lHcsuplsdnlr(2,2) = ydRG(1,1) * (-sgnmu*modmu*sinbeta + 
     $     ADRG(1,1)*cosbeta)

!-------------------------------------------------------------------------

      intm1(1) = lHcsuplsdnlr(1,1)
      intm1(2) = lHcsuplsdnlr(1,2)

      call rmat2d(thetad,rthetad)
            
      intm2(1) = rthetad(1,1)*intm1(1) + rthetad(1,2)*intm1(2)
      intm2(2) = rthetad(2,1)*intm1(1) + rthetad(2,2)*intm1(2)

 
      lHcsuplsdn12(1, 1) = intm2(1)
      lHcsuplsdn12(1, 2) = intm2(2)

      intm1(1) = lHcsuplsdnlr(2, 1)
      intm1(2) = lHcsuplsdnlr(2, 2)

            
      intm2(1) = rthetad(1,1)*intm1(1) + rthetad(1,2)*intm1(2)
      intm2(2) = rthetad(2,1)*intm1(1) + rthetad(2,2)*intm1(2)

      lHcsuplsdn12(2, 1) = intm2(1)
      lHcsuplsdn12(2, 2) = intm2(2)

!---------------------------------------------------------------------------
!     <------------ (H+ G+, L R) basis

      lHcsuprsdnlr(1,1) = yuRG(1,1) * (-sgnmu*modmu*sinbeta - 
     $     AURG(1,1)*cosbeta)

      lHcsuprsdnlr(1,2) = yuRG(1,1)*ydRG(1,1) * (- vev1*cosbeta - 
     $     vev2*sinbeta)/dsqrt(2.d0)

      lHcsuprsdnlr(2,1) = -yuRG(1,1) * (-sgnmu*modmu*cosbeta + 
     $     AURG(1,1)*sinbeta)  

      lHcsuprsdnlr(2,2) = 0.d0

!-------------------------------------------------------------------------

      intm1(1) = lHcsuprsdnlr(1,1)
      intm1(2) = lHcsuprsdnlr(1,2)

      call rmat2d(thetad,rthetad)
             
      intm2(1) = rthetad(1,1)*intm1(1) + rthetad(1,2)*intm1(2)
      intm2(2) = rthetad(2,1)*intm1(1) + rthetad(2,2)*intm1(2)

      lHcsuprsdn12(1,1) = intm2(1)
      lHcsuprsdn12(1,2) = intm2(2)

      intm1(1) = lHcsuprsdnlr(2,1)
      intm1(2) = lHcsuprsdnlr(2,2)

              
      intm2(1) = rthetad(1,1)*intm1(1) + rthetad(1,2)*intm1(2)
      intm2(2) = rthetad(2,1)*intm1(1) + rthetad(2,2)*intm1(2)
  
      lHcsuprsdn12(2, 1) = intm2(1)
      lHcsuprsdn12(2, 2) = intm2(2)

!------------------------------------------------------------------------
C     Feynman rules - neutralino
!------------------------------------------------------------------------


      aPsi0usupr(1) = -4.d0*gp/(3.d0*dsqrt(2.d0))
      aPsi0usupr(2) = 0.d0
      aPsi0usupr(3) = 0.d0
      aPsi0usupr(4) = 0.d0

      bPsi0usupl(1) = gp/(3.d0*dsqrt(2.d0))
      bPsi0usupl(2) = g/dsqrt(2.d0)
      bPsi0usupl(3) = 0.d0
      bPsi0usupl(4) = 0.d0

      aPsi0usupl(1) = 0.d0
      aPsi0usupl(2) = 0.d0
      aPsi0usupl(3) = 0.d0 
      aPsi0usupl(4) = yuRG(1,1)

      bPsi0usupr(1) = 0.d0
      bPsi0usupr(2) = 0.d0
      bPsi0usupr(3) = 0.d0
      bPsi0usupr(4) = yuRG(1,1)


      aChi0usupl(1) = ON(1,1)*aPsi0usupl(1) + ON(1,2)*aPsi0usupl(2) +
     $                 ON(1,3)*aPsi0usupl(1) + ON(1,4)*aPsi0usupl(4) 
      aChi0usupl(2) = ON(2,1)*aPsi0usupl(1) + ON(2,2)*aPsi0usupl(2) +
     $                 ON(2,3)*aPsi0usupl(1) + ON(2,4)*aPsi0usupl(4) 
      aChi0usupl(3) = ON(3,1)*aPsi0usupl(1) + ON(3,2)*aPsi0usupl(2) +
     $                 ON(3,3)*aPsi0usupl(1) + ON(3,4)*aPsi0usupl(4) 
      aChi0usupl(4) = ON(4,1)*aPsi0usupl(1) + ON(4,2)*aPsi0usupl(2) +
     $                 ON(4,3)*aPsi0usupl(1) + ON(4,4)*aPsi0usupl(4) 


      bChi0usupl(1) = ON(1,1)*bPsi0usupl(1) + ON(1,2)*bPsi0usupl(2) +
     $                 ON(1,3)*bPsi0usupl(1) + ON(1,4)*bPsi0usupl(4) 
      bChi0usupl(2) = ON(2,1)*bPsi0usupl(1) + ON(2,2)*bPsi0usupl(2) +
     $                 ON(2,3)*bPsi0usupl(1) + ON(2,4)*bPsi0usupl(4) 
      bChi0usupl(3) = ON(3,1)*bPsi0usupl(1) + ON(3,2)*bPsi0usupl(2) +
     $                 ON(3,3)*bPsi0usupl(1) + ON(3,4)*bPsi0usupl(4) 
      bChi0usupl(4) = ON(4,1)*bPsi0usupl(1) + ON(4,2)*bPsi0usupl(2) +
     $                 ON(4,3)*bPsi0usupl(1) + ON(4,4)*bPsi0usupl(4) 




      aChi0usupr(1) = ON(1,1)*aPsi0usupr(1) + ON(1,2)*aPsi0usupr(2) +
     $                 ON(1,3)*aPsi0usupr(3) + ON(1,4)*aPsi0usupr(4)
      aChi0usupr(2) = ON(2,1)*aPsi0usupr(1) + ON(2,2)*aPsi0usupr(2) +
     $                 ON(2,3)*aPsi0usupr(3) + ON(2,4)*aPsi0usupr(4)
      aChi0usupr(3) = ON(3,1)*aPsi0usupr(1) + ON(3,2)*aPsi0usupr(2) +
     $                 ON(3,3)*aPsi0usupr(3) + ON(3,4)*aPsi0usupr(4)
      aChi0usupr(4) = ON(4,1)*aPsi0usupr(1) + ON(4,2)*aPsi0usupr(2) +
     $                 ON(4,3)*aPsi0usupr(3) + ON(4,4)*aPsi0usupr(4)


      bChi0usupr(1) = ON(1,1)*bPsi0usupr(1) + ON(1,2)*bPsi0usupr(2) +
     $                 ON(1,3)*bPsi0usupr(3) + ON(1,4)*bPsi0usupr(4) 
      bChi0usupr(2) = ON(2,1)*bPsi0usupr(1) + ON(2,2)*bPsi0usupr(2) +
     $                 ON(2,3)*bPsi0usupr(3) + ON(2,4)*bPsi0usupr(4) 
      bChi0usupr(3) = ON(3,1)*bPsi0usupr(1) + ON(3,2)*bPsi0usupr(2) +
     $                 ON(3,3)*bPsi0usupr(3) + ON(3,4)*bPsi0usupr(4) 
      bChi0usupr(4) = ON(4,1)*bPsi0usupr(1) + ON(4,2)*bPsi0usupr(2) +
     $                 ON(4,3)*bPsi0usupr(3) + ON(4,4)*bPsi0usupr(4) 


      loopchtllrr: DO i = 1, 4

      fChi0usupLL(i) = (aChi0usupl(i)*aChi0usupl(i) + 
     $     bChi0usupl(i)*bChi0usupl(i))
      gChi0usupLL(i) = (bChi0usupl(i)*aChi0usupl(i) + 
     $     bChi0usupl(i)*aChi0usupl(i))
      fChi0usupRR(i) = (aChi0usupr(i)*aChi0usupr(i) + 
     $     bChi0usupr(i)*bChi0usupr(i))
      gChi0usupRR(i) = (bChi0usupr(i)*aChi0usupr(i) + 
     $     bChi0usupr(i)*aChi0usupr(i))
      fChi0usupLR(i) = (aChi0usupr(i)*aChi0usupl(i) + 
     $     bChi0usupr(i)*bChi0usupl(i))
      gChi0usupLR(i) = (bChi0usupl(i)*aChi0usupr(i) + 
     $     bChi0usupr(i)*aChi0usupl(i))
      
      ENDDO loopchtllrr

!---------------------------------------------------------------------------
c     Feynman Rules for chargino
!---------------------------------------------------------------------------

      aPsicdsupl(1) = g
      aPsicdsupl(2) = 0.d0
      aPsicdsupr(1) = 0.d0
      aPsicdsupr(2) = -yuRG(1,1)
      bPsicdsupl(1) = 0.d0
      bPsicdsupl(2) = -ydRG(1,1)
      bPsicdsupr(1) = 0.d0
      bPsicdsupr(2) = 0.d0
     
     

      aChicdsupl(1) = OCR(1,1)*aPsicdsupl(1) + OCR(1,2)*aPsicdsupl(2)
      aChicdsupl(2) = OCR(2,1)*aPsicdsupl(1) + OCR(2,2)*aPsicdsupl(2)

      bChicdsupl(1) = OCL(1,1)*bPsicdsupl(1) + OCL(1,2)*bPsicdsupl(2)
      bChicdsupl(2) = OCL(2,1)*bPsicdsupl(1) + OCL(2,2)*bPsicdsupl(2)

      aChicdsupr(1) = OCR(1,1)*aPsicdsupr(1) + OCR(1,2)*aPsicdsupr(2)
      aChicdsupr(2) = OCR(2,1)*aPsicdsupr(1) + OCR(2,2)*aPsicdsupr(2)

      bChicdsupr(1) = OCL(1,1)*bPsicdsupr(1) + OCL(1,2)*bPsicdsupr(2)
      bChicdsupr(2) = OCL(2,1)*bPsicdsupr(1) + OCL(2,2)*bPsicdsupr(2)

      
      loopchst: DO i = 1, 2

      fChdsupLL(i) = (aChicdsupl(i)*aChicdsupl(i) +
     $     bChicdsupl(i)*bChicdsupl(i))
      gChdsupLL(i) = (bChicdsupl(i)*aChicdsupl(i) +
     $     aChicdsupl(i)*bChicdsupl(i))
      fChdsupLR(i) = (aChicdsupl(i)*aChicdsupr(i) +
     $     bChicdsupl(i)*bChicdsupr(i))
      gChdsupLR(i) = (bChicdsupl(i)*aChicdsupr(i) +
     $     aChicdsupl(i)*bChicdsupr(i))
      fChdsupRR(i) = (aChicdsupr(i)*aChicdsupr(i) +
     $     bChicdsupr(i)*bChicdsupr(i))
      gChdsupRR(i) = (bChicdsupr(i)*aChicdsupr(i) +
     $     aChicdsupr(i)*bChicdsupr(i))

      ENDDO loopchst

!--------------------------------------------------------------------
C      Corrections Begin
!------------------------------------------------------------------

       call funcg(p,M3t,mu,q,ggu)

  !     print*,"ggt"
       
       call f(p,mu1,0.d0,q,fu10)
       call f(p,mu2,0.d0,q,fu20)

       call a0(mu2,q,a0mu2)
       call a0(mu1,q,a0mu1)
       call a0(mu2,q,a0u2)
       call a0(mu1,q,a0u1)

       call a0(mc2,q,a0mc2)
       call a0(mc1,q,a0mc1)
       call a0(mc2,q,a0c2)
       call a0(mc1,q,a0c1)
       call a0(mt2,q,a0t2)
       call a0(mt1,q,a0t1)

       call a0(ms2,q,a0ms2)
       call a0(ms1,q,a0ms1)
       call a0(md2,q,a0md2)
       call a0(md1,q,a0md1)
       call a0(mb2,q,a0mb2)
       call a0(mb1,q,a0mb1)

       call a0(meR,q,a0me2)
       call a0(meL,q,a0me1)
       call a0(mmuR,q,a0mmu2)
       call a0(mmuL,q,a0mmu1)
       call a0(mtau2,q,a0mtau2)
       call a0(mtau1,q,a0mtau1)

       call a0(snu(1),q,a0snu1)
       call a0(snu(2),q,a0snu2)
       call a0(snu(3),q,a0snu3)


       call a0(mA0,q,a0mA)
       call a0(mh0,q,a0mh)
       call a0(mHu0,q,a0mHu)
       call a0(mHpm,q,a0mHpm)
       call a0(MW,q,a0Mw)
       call a0(MZ,q,a0Mz)

       call f(p,md1,MW,q,fmd1Mw)
       call f(p,md2,MW,q,fmd2Mw)
       call f(p,mu1,MZ,q,fmu1Mz)
       call f(p,mu2,MZ,q,fmu2Mz)


       call b0(p,M3t,mu,q,b0M3mu)
!-------------------------
      initi:do i = 1,2
      initj:do j= 1,2

      gtterm(i,j) = 0.d0
      cstop(i,j) = 0.d0
      csbtm(i,j) = 0.d0
      higgsterm(i,j) = 0.d0

      enddo initj
      enddo initi
    
!--------------------------------------------------------------
      gtterm(1, 1) = 4.d0 * ((g3**2.d0)/3.d0) *
     $      (2.d0*ggu + (costhetau)**2.d0 * (fu10 + a0mu1) + 
     $      (sinthetau)**2.d0*(fu20 + a0mu2))

      gtterm(2, 2) = 4.d0 * ((g3**2.d0)/3.d0) *
     $     (2.d0*ggu + (sinthetau)**2.d0*(fu10 + a0mu1) + 
     $     (costhetau)**2.d0 * (fu20 + a0mu2))

      gtterm(1, 2) = 4.d0 * ((g3**2.d0)/3.d0) *
     $     (4.d0 * M3t*mu * b0M3mu + 
     $     sinthetau*costhetau * (fu10 - a0mu1 - fu20 + a0mu2))

!----------------------------------------------------------------

      cstop(1, 1) = (yuRG(1,1)**2.d0) * ((sinthetau)**2.d0 * a0u1 + 
     $   (costhetau)**2.d0 * a0u2)

      cstop(2, 2) = (yuRG(1,1)**2.d0) * ((costhetau)**2.d0 * a0u1 + 
     $     (sinthetau)**2.d0 * a0u2)

      cstop(1, 2) = (yuRG(1,1)**2.d0) * costhetau*sinthetau * 3.d0 *
     $     (a0u1 - a0u2)
!------------------------------------------------------------------

      csbtm(1, 1) = (ydRG(1,1)**2.d0) * ((sinthetad)**2.d0 * a0md1 + 
     $     (costhetad)**2.d0 * a0md2)
      
      csbtm(2, 2) = (yuRG(1,1)**2.d0) * ((costhetad)**2.d0 * a0md1 + 
     $     (sinthetad)**2.d0 * a0md2)

!--------------------------------------------------------------------

      higgsterm(1,1) = 0.5d0 * ((yuRG(1,1)**2.d0) * dnu(1) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((yuRG(1,1)**2.d0) * dnu(2) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((yuRG(1,1)**2.d0) * dnu(3) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((yuRG(1,1)**2.d0) * dnu(4) - 
     $     ((g*g) * guL * 0.5d0/(costhw**2.d0)) * cn(4)) * a0mA

      
      higgsterm(2,2) = 0.5d0 * ((yuRG(1,1)**2.d0) * dnu(1) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((yuRG(1,1)**2.d0) * dnu(2) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((yuRG(1,1)**2.d0) * dnu(3) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((yuRG(1,1)**2.d0) * dnu(4) - 
     $     ((g*g) * guR * 0.5d0/(costhw**2.d0)) * cn(4)) * a0mA


      higgsterm(1,1) = higgsterm(1,1) + 
     $     ((ydRG(1,1)**2.d0) * dnu(3) + ((g*g) * 
     $     ((guL*0.5d0/(costhw**2.d0)) - 0.5d0) * cn(3))) * a0mhpm +
     $     ((ydRG(1,1)**2.d0) * dnu(4) + ((g*g) * 
     $     ((guL*0.5d0/(costhw**2.d0)) - 0.5d0) * cn(4))) * a0MW


      higgsterm(2,2) = higgsterm(2,2) + 
     $     ((yuRG(1,1)**2.d0) * dnd(3) + 
     $     ((g*g) * (guR * 0.5d0/(costhw**2.d0)) * cn(3))) * a0mhpm +
     $     ((yuRG(1,1)**2.d0) * dnd(4) + 
     $     ((g*g) * (guR * 0.5d0/(costhw**2.d0)) * cn(4))) * a0MW
     

!--------------------------

      call b0(p,mHu0,mu1,q,b0h0msu(1,1))
      call b0(p,mHu0,mu2,q,b0h0msu(1,2))
      call b0(p,mh0,mu1,q,b0h0msu(2,1))
      call b0(p,mh0,mu2,q,b0h0msu(2,2))
      call b0(p,MZ,mu1,q,b0h0msu(3,1))
      call b0(p,MZ,mu2,q,b0h0msu(3,2))
      call b0(p,mA0,mu1,q,b0h0msu(4,1))
      call b0(p,mA0,mu2,q,b0h0msu(4,2))
      

      loophi: DO i = 1, 4       !<-------------neutral higgs terms
      loophj: DO j = 1, 2

      higgsterm(1,1) = higgsterm(1,1) + 
     $     (lHsupLsup12(i,j)**2.d0) * b0h0msu(i,j)

      higgsterm(1,2) = higgsterm(1,2) + 
     $     lHsupLsup12(i,j) * lHsupRsup12(i,j) * b0h0msu(i,j)

      higgsterm(2,2) = higgsterm(2,2) + 
     $     (lHsupRsup12(i,j)**2.d0) * b0h0msu(i,j)
      
      ENDDO loophj
      ENDDO loophi

!-------------------------------------

      call b0(p,md1,mHpm,q,b0hcmsd(1,1))
      call b0(p,md2,mHpm,q,b0hcmsd(1,2))
      call b0(p,md1,MW,q,b0hcmsd(2,1))
      call b0(p,md2,MW,q,b0hcmsd(2,2))
      

      loophci: DO i = 1, 2      !<-------------charged higgs terms
      loophcj: DO j = 1, 2

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (lHcsuplsdn12(i,j)**2.d0) * b0hcmsd(i,j)        

      higgsterm(1, 2) = higgsterm(1,2) + 
     $     lHcsuplsdn12(i, j) * lHcsuprsdn12(i, j) * b0hcmsd(i,j)

      higgsterm(2, 2) = higgsterm(2, 2) + 
     $     (lHcsuprsdn12(i,j)**2.d0) * b0hcmsd(i,j)
      
      ENDDO loophcj
      ENDDO loophci

!-----------------------------------------------

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (4.d0*(g*g)/costhw**2.d0) * (guL*guL) * a0Mz  + 
     $     2.d0*(g*g) * a0Mw + 
     $     ((4.d0/9.d0)*g*g*sinsqthw)*((costhetau**2.d0) * fu10 + 
     $     (sinthetau**2.d0) * fu20) +
     $     (g*guL/costhw)**2.d0 * ((costhetau**2.d0) * fmu1MZ  +
     $     (sinthetau**2.d0) * fmu2MZ) +
     $     (g*g)*0.5d0*((costhetad**2.d0)*fmd1Mw + (sinthetad**2.d0) *
     $     fmd2Mw) + 
     $     (g*g)*0.25d0*((costhetau**2.d0)*a0u1 + (sinthetau**2.d0) *
     $     a0u2 + 2.d0 * ((costhetad**2.d0)*a0md1 + (sinthetad**2.d0) *
     $     a0md2)) +
     $     (g*g) * 0.5d0 * (1.5d0*a0mu1 + 1.5d0*a0mc1 +
     $     1.5d0*((costhetat**2.d0)*a0t1 + 
     $     (sinthetat**2.d0)*a0t2) -
     $     1.5d0*a0md1 - 1.5d0*a0ms1 -
     $     1.5d0*((costhetab**2.d0)*a0mb1 +
     $     (sinthetab**2.d0)*a0mb2) +
     $     0.5d0*(a0snu1 + a0snu2 + a0snu3) -
     $     0.5d0*(a0me1 + a0mmu1 +
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2)) +
     $     (gp**2.d0)*0.25d0*(yuL**2.d0) * ((costhetau**2.d0) * a0u1 +
     $     (sinthetau**2.d0) * a0u2) +
     $     (gp**2.d0)*0.25d0 * yuL * (3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 +(sinthetat**2.d0)*a0t2) +
     $     3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 +(costhetat**2.d0)*a0t2) +
     $     3.d0 * ydL * (a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0mb1 + (sinthetab**2.d0)*a0mb1) +
     $     3.d0 * ydR * (a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0mb1 + (costhetab**2.d0)*a0mb2) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 +(costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))

    
!-----------------------------------------------------

      higgsterm(2, 2) = higgsterm(2,2) +
     $     (4.d0*(g*g)/costhw**2.d0) * (guR**2.d0) * a0MZ + 
     $     ((4.d0/9.d0)*g*g*sinsqthw) * ((sinthetau)*fu10 + 
     $     (costhetau**2.d0)*fu20) +
     $     (g*guR/costhw)**2.d0 * ((sinthetau**2.d0)*fmu1MZ + 
     $     (costhetau**2.d0)*fmu2MZ) +
     $     (gp**2.d0) * 0.25d0 * (yuR**2.d0) *
     $     ((sinthetau**2.d0)*a0u1 + (costhetau**2.d0)*a0u2) +
     $     (gp**2.d0) * 0.25d0 * yuR * 
     $     (3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t2) +
     $     3.d0 * yuR * (a0mu2 +  a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 + (costhetat**2.d0)*a0t2) +
     $     3.d0 * ydL * (a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0mb1 + (sinthetab**2.d0)*a0mb2) +
     $     3.d0 * ydR * (a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0mb1 + (costhetab**2.d0)*a0mb2) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))


!----------------------------------------------

      higgsterm(1,2) = higgsterm(1,2) + 
     $     (gp**2.d0) * yuL*yuR * sinthetau*costhetau *
     $     (a0u1 - a0u2) +
     $     ((4.d0/9.d0) * g*g * sinsqthw) * sinthetau*costhetau * 
     $     (fu10 - fu20) -
     $     ((g*g)/costhw**2.d0) * guL*guR * sinthetau*costhetau *
     $     (fmu1Mz - fmu2Mz)

!----------------------------------------------------------------------
C     Chargino term
!----------------------------------------------------------------------
      
      chargino(1,1) = 0.d0
      chargino(1,2) = 0.d0
      chargino(2,1) = 0.d0
      chargino(2,2) = 0.d0

      call funcg(p,mchargino(1),ms,q,gmchargino1ms)
      call funcg(p,mchargino(2),ms,q,gmchargino2ms)
      call b0(p,mchargino(1),ms,q,b0mchargino1ms)
      call b0(p,mchargino(2),ms,q,b0mchargino2ms)
      
      chargino(1, 1) = chargino(1, 1) + fChdsupLL(1)*gmchargino1ms -
     $     gChdsupLL(1)*mchargino(1)*md*b0mchargino1ms*2.d0 +
     $     fChdsupLL(2)*gmchargino2ms -
     $     gChdsupLL(2)*mchargino(2)*md*b0mchargino2ms*2.d0

      chargino(1, 2) = chargino(1, 2) + fChdsupLR(1)*gmchargino1ms -
     $     gChdsupLR(1)*mchargino(1)*md*b0mchargino1ms*2.d0 +
     $     fChdsupLR(2)*gmchargino2ms -
     $     gChdsupLR(2)*mchargino(2)*md*b0mchargino2ms*2.d0

      chargino(2, 2) = chargino(2, 2) + fChdsupRR(1)*gmchargino1ms - 
     $     gChdsupRR(1)*mchargino(1)*md*b0mchargino1ms*2.d0 +
     $     fChdsupRR(2)*gmchargino2ms - 
     $     gChdsupRR(2)*mchargino(2)*md*b0mchargino2ms*2.d0

!---------------------------------------------------------------------
C     neutralino terms
!---------------------------------------------------------------------

      neutralino(1,1) = 0.d0
      neutralino(1,2) = 0.d0
      neutralino(2,1) = 0.d0
      neutralino(2,2) = 0.d0

      call b0(p,nmneut(1),mu,q,b0mneut1mu)
      call b0(p,nmneut(2),mu,q,b0mneut2mu)
      call b0(p,nmneut(3),mu,q,b0mneut3mu)
      call b0(p,nmneut(4),mu,q,b0mneut4mu)
      
      call funcg(p,nmneut(1),mu,q,gmneut1mu)
      call funcg(p,nmneut(2),mu,q,gmneut2mu)
      call funcg(p,nmneut(3),mu,q,gmneut3mu)
      call funcg(p,nmneut(4),mu,q,gmneut4mu)

      neutralino(1, 1) = neutralino(1, 1) +
     $     fChi0usupLL(1)*gmneut1mu - gChi0usupLL(1)*2.d0*
     $     mneut(1)*mu*b0mneut1mu +
     $     fChi0usupLL(2)*gmneut2mu - gChi0usupLL(2)*2.d0*
     $     mneut(2)*mu*b0mneut2mu +
     $     fChi0usupLL(3)*gmneut3mu - gChi0usupLL(3)*2.d0*
     $     mneut(3)*mu*b0mneut3mu +
     $     fChi0usupLL(4)*gmneut4mu - gChi0usupLL(4)*2.d0*
     $     mneut(4)*mu*b0mneut4mu


      neutralino(2, 2) = neutralino(2, 2) +
     $     fChi0usupRR(i)*gmneut1mu - gChi0usupRR(1)*2.d0*
     $     mneut(1)*mu*b0mneut1mu +
     $     fChi0usupRR(2)*gmneut2mu - gChi0usupRR(2)*2.d0*
     $     mneut(2)*mu*b0mneut2mu +
     $     fChi0usupRR(3)*gmneut3mu - gChi0usupRR(3)*2.d0*
     $     mneut(3)*mu*b0mneut3mu +
     $     fChi0usupRR(4)*gmneut4mu - gChi0usupRR(4)*2.d0*
     $     mneut(4)*mu*b0mneut4mu

      neutralino(1, 2) = neutralino(1, 2) +
     $     fChi0usupLR(1)*gmneut1mu - gChi0usupLR(1)*2.d0*
     $     mneut(1)*mu*b0mneut1mu +
     $     fChi0usupLR(2)*gmneut2mu - gChi0usupLR(2)*2.d0*
     $     mneut(2)*mu*b0mneut2mu +
     $     fChi0usupLR(3)*gmneut3mu - gChi0usupLR(3)*2.d0*
     $     mneut(3)*mu*b0mneut3mu +
     $     fChi0usupLR(4)*gmneut4mu - gChi0usupLR(4)*2.d0*
     $     mneut(4)*mu*b0mneut4mu 
      
!-------------------------------------------------------------------

      MSQU3(1,1) = mSQRG(1,1) + mu**2.d0 + guL*MZ*MZ*dcos(2.d0*beta)
      MSQU3(1,2) = mu * ((AURG(1,1) - sgnmu*modmu/dtan(beta)))
      MSQU3(2,1) = MSQU3(1,2)
      MSQU3(2,2) = mSURG(2,2) + mu**2.d0 + guR*MZ*MZ*dcos(2.d0*beta)

!----------------------------------------------------------------------

      pisup(1,1) = (cstop(1,1) + csbtm(1,1) + gtterm(1,1) +
     $     chargino(1,1) + neutralino(1,1) + higgsterm(1,1))/
     $     (16.d0*pi*pi)

      pisup(1,2) = (cstop(1,2) + csbtm(1,2) + gtterm(1,2) +
     $     chargino(1,2) + neutralino(1,2) + higgsterm(1,2))/
     $     (16.d0*pi*pi)

      pisup(2,2) = (cstop(2,2) + csbtm(2,2) + gtterm(2,2) +
     $     chargino(2,2) + neutralino(2,2) + higgsterm(2,2))/
     $     (16.d0*pi*pi)
      
      pisup(2,1) = pisup(1,2)


      pisup(1,1) = MSQU3(1,1) - pisup(1,1)
      pisup(1,2) = MSQU3(1,2) - pisup(1,2)
      pisup(2,2) = MSQU3(2,2) - pisup(2,2)
      pisup(2,1) = pisup(1,2) 

!----------------------------------------------------------------------

C     find the singular values and the diagonalising matrices

      info  = 10
      AOK = 0
      
   !   call dsyev('V','U',2,pisup,2,SUeg,work,lwork,info)
      
      if(info.eq.0) then
         AOK = AOK + 1
      endif
      
      Call CEigensystem(2,pisup,2,SUeg,sumix,2,1)
      
      RETURN

      END SUBROUTINE pisupq

C==============================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     
C     pistaul is checked on 22/05/2010 @ 17:30.
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


      SUBROUTINE pistaul(p,q,g,gp,mtau,tanbeta,
     $     mSLRG,mSERG,yeRG,AERG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,M3t,delpistau,STaueg)


 
      IMPLICIT NONE 

      INTEGER i,j,info,AOK,lwork
      parameter(lwork = 35)
      DOUBLE PRECISION work(lwork)
      DOUBLE PRECISION tanbeta,mSLRG(3,3),mSERG(3,3)
      DOUBLE PRECISION AERG(3,3),SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yeRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      double precision  mTau
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3),nmneut(4)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g,M3t  
      double precision vev1,beta
 
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetac,thetas,thetamu

      DOUBLE PRECISION thetau,thetad,thetae

      DOUBLE PRECISION pistau(2,2)
      data pistau/ 4 * 0.d0/

      DOUBLE PRECISION higgsterm(2,2)
      data higgsterm/ 4 * 0.d0/
      
      DOUBLE PRECISION dnu(4),dnd(4),cn(4)

      DOUBLE PRECISION lstaulstaulr(4,2),lstaulstau12(4,2)
      DOUBLE PRECISION lstaurstaulr(4,2),lstaurstau12(4,2)
      
      DOUBLE PRECISION aPsi0tstaur(4), bPsi0tstaur(4), aPsi0tstaul(4)
      DOUBLE PRECISION bPsi0tstaul(4) 
      DOUBLE PRECISION aChi0tstaul(4), bChi0tstaul(4), aChi0tstaur(4)
      DOUBLE PRECISION bChi0tstaur(4)

      DOUBLE PRECISION gChi0taustauLL(4), fChi0taustauLL(4)
      DOUBLE PRECISION gChi0taustauLR(4), fChi0taustauLR(4)
      DOUBLE PRECISION gChi0taustauRR(4), fChi0taustauRR(4)

      DOUBLE PRECISION bPsicnustaul(2), bPsicnustaur(2), aPsicnustaul(2)
      DOUBLE PRECISION aPsicnustaur(2)
      DOUBLE PRECISION aChicnustaur(2), aChicnustaul(2)
      DOUBLE PRECISION bChicnustaul(2),bChicnustaur(2)

      DOUBLE PRECISION gtterm(2,2), cstop(2,2), csbtm(2,2) 
      data gtterm/ 4 * 0.d0/, cstop/ 4 * 0.d0/, csbtm/ 4 * 0.d0/

      DOUBLE PRECISION  chargino(2,2), neutralino(2,2), mSLD3(2,2)
      data chargino/ 4 * 0.d0/,neutralino/ 4 * 0.d0/, mSLD3/ 4 * 0.d0/

      DOUBLE PRECISION fChnustauLL(2), gChnustauLL(2) 
      DOUBLE PRECISION fChnustauLR(2), gChnustauLR(2) 
      DOUBLE PRECISION fChnustauRR(2), gChnustauRR(2)
 
      DOUBLE PRECISION intm1(2), intm2(2)
      data intm1/ 2 * 0.d0/, intm2/ 2 * 0.d0/

      DOUBLE PRECISION lHcstaursnulr(2,2),lHcstaursnu12(2,2)
      DOUBLE PRECISION lHcstaulsnu12(2,2)
      DOUBLE PRECISION lHcstaulsnulr(2,2) 

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION lhstauslstau12(4,2),lhstaurstau12(4,2)
      data lhstauslstau12/ 8 * 0.d0/, lhstaurstau12/ 8 * 0.d0/

      
      DOUBLE PRECISION rthetatau(2,2),ralpha(2,2)

      DOUBLE PRECISION ggtau,ftau10,ftau20,a0b1,a0b2,b0M3mtau,a0t1,a0t2
      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz
      DOUBLE PRECISION b0h0mstau(4,2),b0hcmstop(2,2),fsnu1MW
      data b0h0mstau/ 8 * 0.d0/,b0hcmstop/ 4 * 0.d0/

      DOUBLE PRECISION a0md1,a0md2,a0ms1,a0ms2,a0mtau1,a0mtau2,a0snu2
      DOUBLE PRECISION a0me1,a0me2,a0mmu1,a0mmu2,a0snu1,a0snu3
      DOUBLE PRECISION a0mc1,a0mc2,a0mu1,a0mu2
      
      DOUBLE PRECISION fmtau1MZ,fmtau2MZ,f0MW
      DOUBLE PRECISION gmchargino1mtau,gmchargino2mtau
      DOUBLE PRECISION b0mchargino1mtau,b0mchargino2mtau,b0mneut4mtau
      DOUBLE PRECISION b0mneut1mtau,b0mneut2mtau,b0mneut3mtau
      DOUBLE PRECISION gmneut1mtau,gmneut2mtau,gmneut3mtau,gmneut4mtau
      
      DOUBLE PRECISION STaueg(2),staumix(2,2)
      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2
      DOUBLE PRECISION delpistau(2,2)
  
      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae


      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d
  

      include 'stdinputs.h'


C------------------------------------------------------------

      pi = 4.d0 * datan(1.d0)
      MZ = MZpole
      MW = MWpole

      STaueg(1) = 0.d0
      STaueg(2) = 0.d0
     
      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw = dsqrt(sinsqthw)
      thw = dasin(sinthw)
      costhw = dcos(thw)
      cos2thw = dcos(2.d0*thw)
 

      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))
      
      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))
      
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

      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)

      mt1 = mtL
      mt2 = mtR

      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR

C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))
      
      beta     = datan(tanbeta)
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha  = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta  = dcos(datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)

!-----------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0
!----------------------------------------------------------
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

!---------------------------------------------------------------    

      lstaulstaulr(1,1) = (g*MZ*geL*cosbeta/costhw) + 
     $     (yeRG(3,3)**2.d0) * vev1
      lstaulstaulr(1,2) = (yeRG(3,3)*AERG(3,3))/dsqrt(2.d0)

      lstaulstaulr(2,1) = ((- g*MZ*geL*sinbeta)/costhw)
      lstaulstaulr(2,2) = (-yeRG(3,3)*sgnmu*modmu)/dsqrt(2.d0)

      lstaulstaulr(3,1) = 0.d0
      lstaulstaulr(3,2) = (-1.d0/dsqrt(2.d0)) *
     $     (-sgnmu*modmu*sinbeta*yeRG(3,3) + 
     $     yeRG(3,3)*AERG(3,3)*cosbeta)

      lstaulstaulr(4,1) = 0.d0
      lstaulstaulr(4,2) = (-1.d0/dsqrt(2.d0)) * 
     $     (-sgnmu*modmu*cosbeta*yeRG(3,3) - 
     $     yeRG(3,3)*AERG(3,3)*sinbeta)



      lstaurstaulr(1,1) = lstaulstaulr(1,2)
      lstaurstaulr(1,2) = (g*MZ*geR*cosbeta/costhw) + 
     $     yeRG(3,3)**2.d0 * vev1

      lstaurstaulr(2,1) = lstaulstaulr(2,2)
      lstaurstaulr(2,2) = -g*MZ*geR*sinbeta/costhw

      lstaurstaulr(3,1) = -lstaulstaulr(3,2)
      lstaurstaulr(3,2) = 0.d0

      lstaurstaulr(4,1) = -lstaulstaulr(4,2)
      lstaurstaulr(4,2) = 0.d0


!--------------------------------------------------------------  
  
      loopmixst: DO i = 1, 4

      intm1(1) = lstaulstaulr(i,1)
      intm1(2) = lstaulstaulr(i,2)

      call rmat2d(thetatau,rthetatau)
            
      intm2(1) = rthetatau(1,1)*intm1(1) + rthetatau(1,2)*intm1(2)
      intm2(2) = rthetatau(2,1)*intm1(1) + rthetatau(2,2)*intm1(2)

      lstaulstau12(i,1) = intm2(1)
      lstaulstau12(i,2) = intm2(2)

      intm1(1) = lstaurstaulr(i,1)
      intm1(2) = lstaurstaulr(i,2)

      intm2(1) = rthetatau(1,1)*intm1(1) + rthetatau(1,2)*intm1(2)
      intm2(2) = rthetatau(2,1)*intm1(1) + rthetatau(2,2)*intm1(2)

      lstaurstau12(i,1) = intm2(1)
      lstaurstau12(i,2) = intm2(2)

      ENDDO loopmixst

!--------------------------------------------------------------------------
C       Mixing CP-even Higgs 
!--------------------------------------------------------------------------

      loopmixcpeven: DO i = 1, 2
      
      intm1(1) = lstaulstau12(1,i)
      intm1(2) = lstaulstau12(2,i)

      call rmat2d(alpha,ralpha)
            
      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lhstauslstau12(1,i) = intm2(1)
      lhstauslstau12(2,i) = intm2(2)

      intm1(1) = lstaurstau12(1,i)
      intm1(2) = lstaurstau12(2,i)

      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lhstaurstau12(1,i) = intm2(1)
      lhstaurstau12(2,i) = intm2(2)
     
      ENDDO loopmixcpeven
!----------------------------------------------------------------------
C     Feynman rules - charged higgs
!<--------- (H+ G+, L R) basis

      lHcstaulsnulr(1, 1) = ( g*MW*sin2beta - 
     $     (yeRG(3,3)**2.d0) * vev1 * sinbeta)/dsqrt(2.d0)

      lHcstaulsnulr(1,2) = 0.d0 

      lHcstaulsnulr(2,1) = (-g*MW*cos2beta - 
     $     yeRG(3,3)**2.d0 * vev1 * cosbeta)/dsqrt(2.d0)

      lHcstaulsnulr(2,2) = 0.d0


!-----------------------------------------------

      intm1(1) = lHcstaulsnulr(1,1)
      intm1(2) = lHcstaulsnulr(1,2)

      call rmat2d(thetatau,rthetatau)
            
      intm2(1) = rthetatau(1,1)*intm1(1) + rthetatau(1,2)*intm1(2)
      intm2(2) = rthetatau(2,1)*intm1(1) + rthetatau(2,2)*intm1(2)

 
      lHcstaulsnu12(1, 1) = intm2(1)
      lHcstaulsnu12(1, 2) = intm2(2)

      intm1(1) = lHcstaulsnulr(2, 1)
      intm1(2) = lHcstaulsnulr(2, 2)


                  
      intm2(1) = rthetatau(1,1)*intm1(1) + rthetatau(1,2)*intm1(2)
      intm2(2) = rthetatau(2,1)*intm1(1) + rthetatau(2,2)*intm1(2)

      lHcstaulsnu12(2, 1) = intm2(1)
      lHcstaulsnu12(2, 2) = intm2(2)

!---------------------------------------------------------------------------
!<------------ (H+ G+, L R) basis

      lHcstaursnulr(1, 1) = yeRG(3,3)*(-sgnmu*modmu*cosbeta - 
     $     AERG(3,3)*sinbeta)

      lHcstaursnulr(1, 2) = 0.d0

      lHcstaursnulr(2, 1) = yeRG(3,3) * (-sgnmu*modmu*sinbeta + 
     $     AERG(3,3)*cosbeta)  

      lHcstaursnulr(2, 2) = 0.d0

!-------------------------------------------------------------------------

      intm1(1) = lHcstaursnulr(1,1)
      intm1(2) = lHcstaursnulr(1,2)

      call rmat2d(thetatau,rthetatau)
            
      intm2(1) = rthetatau(1,1)*intm1(1) + rthetatau(1,2)*intm1(2)
      intm2(2) = rthetatau(2,1)*intm1(1) + rthetatau(2,2)*intm1(2)

      lHcstaursnu12(1,1) = intm2(1)
      lHcstaursnu12(1,2) = intm2(2)

      intm1(1) = lHcstaursnulr(2,1)
      intm1(2) = lHcstaursnulr(2,2)

                 
      intm2(1) = rthetatau(1,1)*intm1(1) + rthetatau(1,2)*intm1(2)
      intm2(2) = rthetatau(2,1)*intm1(1) + rthetatau(2,2)*intm1(2)
  
      lHcstaursnu12(2, 1) = intm2(1)
      lHcstaursnu12(2, 2) = intm2(2)

!------------------------------------------------------------------------
C     Feynman rules - neutralino
!------------------------------------------------------------------------

      aPsi0tstaur(1) = yeR*gp/(dsqrt(2.d0))
      aPsi0tstaur(2) = 0.d0
      aPsi0tstaur(3) = 0.d0
      aPsi0tstaur(4) = 0.d0

      bPsi0tstaul(1) = yeL*gp/(dsqrt(2.d0))
      bPsi0tstaul(2) = g/dsqrt(2.d0)
      bPsi0tstaul(3) = 0.d0
      bPsi0tstaul(4) = 0.d0

      aPsi0tstaul(1) = 0.d0
      aPsi0tstaul(2) = 0.d0
      aPsi0tstaul(3) = 0.d0 
      aPsi0tstaul(4) = yeRG(3,3)

      bPsi0tstaur(1) = 0.d0
      bPsi0tstaur(2) = 0.d0
      bPsi0tstaur(3) = 0.d0
      bPsi0tstaur(4) = yeRG(3,3)


      aChi0tstaul(1) = ON(1,1)*aPsi0tstaul(1) + ON(1,2)*aPsi0tstaul(2) +
     $                 ON(1,3)*aPsi0tstaul(1) + ON(1,4)*aPsi0tstaul(4) 
      aChi0tstaul(2) = ON(2,1)*aPsi0tstaul(1) + ON(2,2)*aPsi0tstaul(2) +
     $                 ON(2,3)*aPsi0tstaul(1) + ON(2,4)*aPsi0tstaul(4) 
      aChi0tstaul(3) = ON(3,1)*aPsi0tstaul(1) + ON(3,2)*aPsi0tstaul(2) +
     $                 ON(3,3)*aPsi0tstaul(1) + ON(3,4)*aPsi0tstaul(4) 
      aChi0tstaul(4) = ON(4,1)*aPsi0tstaul(1) + ON(4,2)*aPsi0tstaul(2) +
     $                 ON(4,3)*aPsi0tstaul(1) + ON(4,4)*aPsi0tstaul(4) 


      bChi0tstaul(1) = ON(1,1)*bPsi0tstaul(1) + ON(1,2)*bPsi0tstaul(2) +
     $                 ON(1,3)*bPsi0tstaul(1) + ON(1,4)*bPsi0tstaul(4) 
      bChi0tstaul(2) = ON(2,1)*bPsi0tstaul(1) + ON(2,2)*bPsi0tstaul(2) +
     $                 ON(2,3)*bPsi0tstaul(1) + ON(2,4)*bPsi0tstaul(4) 
      bChi0tstaul(3) = ON(3,1)*bPsi0tstaul(1) + ON(3,2)*bPsi0tstaul(2) +
     $                 ON(3,3)*bPsi0tstaul(1) + ON(3,4)*bPsi0tstaul(4) 
      bChi0tstaul(4) = ON(4,1)*bPsi0tstaul(1) + ON(4,2)*bPsi0tstaul(2) +
     $                 ON(4,3)*bPsi0tstaul(1) + ON(4,4)*bPsi0tstaul(4) 




      aChi0tstaur(1) = ON(1,1)*aPsi0tstaur(1) + ON(1,2)*aPsi0tstaur(2) +
     $                 ON(1,3)*aPsi0tstaur(3) + ON(1,4)*aPsi0tstaur(4)
      aChi0tstaur(2) = ON(2,1)*aPsi0tstaur(1) + ON(2,2)*aPsi0tstaur(2) +
     $                 ON(2,3)*aPsi0tstaur(3) + ON(2,4)*aPsi0tstaur(4)
      aChi0tstaur(3) = ON(3,1)*aPsi0tstaur(1) + ON(3,2)*aPsi0tstaur(2) +
     $                 ON(3,3)*aPsi0tstaur(3) + ON(3,4)*aPsi0tstaur(4)
      aChi0tstaur(4) = ON(4,1)*aPsi0tstaur(1) + ON(4,2)*aPsi0tstaur(2) +
     $                 ON(4,3)*aPsi0tstaur(3) + ON(4,4)*aPsi0tstaur(4)


      bChi0tstaur(1) = ON(1,1)*bPsi0tstaur(1) + ON(1,2)*bPsi0tstaur(2) +
     $                 ON(1,3)*bPsi0tstaur(3) + ON(1,4)*bPsi0tstaur(4) 
      bChi0tstaur(2) = ON(2,1)*bPsi0tstaur(1) + ON(2,2)*bPsi0tstaur(2) +
     $                 ON(2,3)*bPsi0tstaur(3) + ON(2,4)*bPsi0tstaur(4) 
      bChi0tstaur(3) = ON(3,1)*bPsi0tstaur(1) + ON(3,2)*bPsi0tstaur(2) +
     $                 ON(3,3)*bPsi0tstaur(3) + ON(3,4)*bPsi0tstaur(4) 
      bChi0tstaur(4) = ON(4,1)*bPsi0tstaur(1) + ON(4,2)*bPsi0tstaur(2) +
     $                 ON(4,3)*bPsi0tstaur(3) + ON(4,4)*bPsi0tstaur(4) 


      loopchtllrr: DO i = 1, 4

      fChi0taustauLL(i) = (aChi0tstaul(i)*aChi0tstaul(i) + 
     $     bChi0tstaul(i)*bChi0tstaul(i))
      gChi0taustauLL(i) = (bChi0tstaul(i)*aChi0tstaul(i) + 
     $     bChi0tstaul(i)*aChi0tstaul(i))
      fChi0taustauRR(i) = (aChi0tstaur(i)*aChi0tstaur(i) + 
     $     bChi0tstaur(i)*bChi0tstaur(i))
      gChi0taustauRR(i) = (bChi0tstaur(i)*aChi0tstaur(i) + 
     $     bChi0tstaur(i)*aChi0tstaur(i))
      fChi0taustauLR(i) = (aChi0tstaur(i)*aChi0tstaul(i) + 
     $     bChi0tstaur(i)*bChi0tstaul(i))
      gChi0taustauLR(i) = (bChi0tstaul(i)*aChi0tstaur(i) + 
     $     bChi0tstaur(i)*aChi0tstaul(i))
      
      ENDDO loopchtllrr

!---------------------------------------------------------------------------
c     Feynman Rules for chargino

      aPsicnustaul(1) = 0.d0
      aPsicnustaul(2) = 0.d0

      aPsicnustaur(1) = 0.d0
      aPsicnustaur(2) = 0.d0

      bPsicnustaul(1) = g
      bPsicnustaul(2) = 0.d0

      bPsicnustaur(1) = 0.d0
      bPsicnustaur(2) = - yeRG(3,3)
      

      aChicnustaul(1) = OCR(1,1)*aPsicnustaul(1) + 
     $     OCR(1,2)*aPsicnustaul(2)
      aChicnustaul(2) = OCR(2,1)*aPsicnustaul(1) + 
     $     OCR(2,2)*aPsicnustaul(2)

      bChicnustaul(1) = OCL(1,1)*bPsicnustaul(1) + 
     $     OCL(1,2)*bPsicnustaul(2)
      bChicnustaul(2) = OCL(2,1)*bPsicnustaul(1) + 
     $     OCL(2,2)*bPsicnustaul(2)

      aChicnustaur(1) = OCR(1,1)*aPsicnustaur(1) + 
     $     OCR(1,2)*aPsicnustaur(2)
      aChicnustaur(2) = OCR(2,1)*aPsicnustaur(1) + 
     $     OCR(2,2)*aPsicnustaur(2)

      bChicnustaur(1) = OCL(1,1)*bPsicnustaur(1) + 
     $     OCL(1,2)*bPsicnustaur(2)
      bChicnustaur(2) = OCL(2,1)*bPsicnustaur(1) + 
     $     OCL(2,2)*bPsicnustaur(2)

      
      loopchst: DO i = 1, 2

      fChnustauLL(i) = (aChicnustaul(i)*aChicnustaul(i) +
     $     bChicnustaul(i)*bChicnustaul(i))
      gChnustauLL(i) = (bChicnustaul(i)*aChicnustaul(i) +
     $     aChicnustaul(i)*bChicnustaul(i))
      fChnustauLR(i) = (aChicnustaul(i)*aChicnustaur(i) +
     $     bChicnustaul(i)*bChicnustaur(i))
      gChnustauLR(i) = (bChicnustaul(i)*aChicnustaur(i) +
     $     aChicnustaul(i)*bChicnustaur(i))
      fChnustauRR(i) = (aChicnustaur(i)*aChicnustaur(i) +
     $     bChicnustaur(i)*bChicnustaur(i))
      gChnustauRR(i) = (bChicnustaur(i)*aChicnustaur(i) +
     $     aChicnustaur(i)*bChicnustaur(i))

      ENDDO loopchst

!--------------------------------------------------------------------
C     Corrections Begin
!------------------------------------------------------------------

      call funcg(p,M3t,mtau,q,ggtau)

!     print*,"ggt"
      
      call f(p,mtau1,0.d0,q,ftau10)
      call f(p,mtau2,0.d0,q,ftau20)

      call a0(muR,q,a0mu2)
      call a0(muL,q,a0mu1)
      call a0(mcR,q,a0mc2)
      call a0(mcL,q,a0mc1)
      call a0(mt2,q,a0t2)
      call a0(mt1,q,a0t1)

      call a0(msR,q,a0ms2)
      call a0(msL,q,a0ms1)
      call a0(mdR,q,a0md2)
      call a0(mdL,q,a0md1)
      call a0(mb2,q,a0b2)
      call a0(mb1,q,a0b1)

      call a0(meR,q,a0me2)
      call a0(meL,q,a0me1)
      call a0(mmuR,q,a0mmu2)
      call a0(mmuL,q,a0mmu1)
      call a0(mtau2,q,a0mtau2)
      call a0(mtau1,q,a0mtau1)

      call a0(snu(1),q,a0snu1)
      call a0(snu(2),q,a0snu2)
      call a0(snu(3),q,a0snu3)


      call a0(mA0,q,a0mA)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)
      call a0(mHpm,q,a0mHpm)
      call a0(MW,q,a0Mw)
      call a0(MZ,q,a0Mz)

      call f(p,0.d0,MW,q,f0Mw)
      call f(p,0.d0,MW,q,f0Mw)
      call f(p,mtau1,MZ,q,fmtau1Mz)
      call f(p,mtau2,MZ,q,fmtau2Mz)

      call f(p,snu(1),MW,q,fsnu1MW)
      call b0(p,M3t,mtau,q,b0M3mtau)
      
!--------------------------------------------------------------
      gtterm(1, 1) = 0.d0

      gtterm(2, 2) = 0.d0

      gtterm(1, 2) = 0.d0
!----------------------------------------------------------------

      csbtm(1, 1) = yeRG(3,3)**2.d0*((sinthetatau)**2.d0*a0mtau1 + 
     $     (costhetatau)**2.d0*a0mtau2)

      csbtm(2, 2) = (yeRG(3,3)**2.d0) * ((costhetatau)**2.d0*a0mtau1 + 
     $     (sinthetatau)**2.d0*a0mtau2)

      csbtm(1, 2) = (yeRG(3,3)**2.d0) * costhetatau*sinthetatau* 1.d0 *
     $     (a0mtau1 - a0mtau2)
!------------------------------------------------------------------

      cstop(1, 1) = 0.d0

      cstop(2, 2) = yeRG(3,3)**2.d0 * a0snu1
!--------------------------------------------------------------------

      higgsterm(1,1) = 0.5d0 * ((yeRG(3,3)**2.d0) * dnd(1) -
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((yeRG(3,3)**2.d0) * dnd(2) - 
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((yeRG(3,3)**2.d0) * dnd(3) - 
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((yeRG(3,3)**2.d0) * dnd(4) - 
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA


      higgsterm(1,2) = 0.d0

      higgsterm(2,2) = 0.5d0*(yeRG(3,3)**2.d0 * dnd(1) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0*(yeRG(3,3)**2.d0 * dnd(2) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0*(yeRG(3,3)**2.d0 * dnd(3) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(3)) * a0Mz +
     $     0.5d0*(yeRG(3,3)**2.d0 * dnd(4) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA
      

      higgsterm(1,1) = higgsterm(1,1) + 
     $     (0.d0 + ((g*g) * 
     $     ((geL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(3))) * a0Mhpm +
     $     (0.d0 + ((g*g) * 
     $     ((geL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(4))) * a0MW

      higgsterm(2,2) = higgsterm(2,2) + 
     $     ((yeRG(3,3)**2.d0) * dnu(3) +  
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(3)) * a0Mhpm +
     $     ((yeRG(3,3)**2.d0) * dnu(4) + 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(4)) * a0MW !<----------------chkd
      
!--------------------------
      
      call b0(p,mHu0,mtau1,q,b0h0mstau(1,1))
      call b0(p,mHu0,mtau2,q,b0h0mstau(1,2))
      call b0(p,mh0,mtau1,q,b0h0mstau(2,1))
      call b0(p,mh0,mtau2,q,b0h0mstau(2,2))
      call b0(p,MZ,mtau1,q,b0h0mstau(3,1))
      call b0(p,MZ,mtau2,q,b0h0mstau(3,2))
      call b0(p,mA0,mtau1,q,b0h0mstau(4,1))
      call b0(p,mA0,mtau2,q,b0h0mstau(4,2))
      
      
      loophi: DO i = 1, 4       !<-------------neutral higgs terms
      loophj: DO j = 1, 2
      
      higgsterm(1,1) = higgsterm(1,1) + 
     $     (lhstauslstau12(i,j)**2.d0) * b0h0mstau(i,j)
      
      higgsterm(1,2) = higgsterm(1,2) + 
     $     lhstauslstau12(i,j) * lhstaurstau12(i,j) * b0h0mstau(i,j)
      
      higgsterm(2,2) = higgsterm(2,2) + 
     $     (lhstaurstau12(i,j)**2.d0) * b0h0mstau(i,j)
      
      ENDDO loophj
      ENDDO loophi
!-------------------------------------
      
      call b0(p,0.d0,mHpm,q,b0hcmstop(1,1))
      call b0(p,0.d0,mHpm,q,b0hcmstop(1,2))
      call b0(p,0.d0,MW,q,b0hcmstop(2,1))
      call b0(p,0.d0,MW,q,b0hcmstop(2,2))
      
      
      loophci: DO i = 1, 2      !<-------------charged higgs terms
      loophcj: DO j = 1, 2
      
      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (lHcstaulsnu12(i,j)**2.d0) * b0hcmstop(i,j)       
      
      higgsterm(1, 2) = higgsterm(1,2) + 
     $     lHcstaulsnu12(i, j) * lHcstaursnu12(i, j) *b0hcmstop(i,j)
      
      higgsterm(2, 2) = higgsterm(2, 2) + 
     $     (lHcstaursnu12(i,j)**2.d0) * b0hcmstop(i,j)
      

      ENDDO loophcj
      ENDDO loophci

!------------------------------------------------------------------------------

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (4.d0*(g*g)/costhw**2.d0)*(geL*geL)*a0MZ  + 
     $     2.d0*(g*g)*a0MW + 
     $     ((1.d0 *g*g*sinsqthw) * ((costhetatau**2.d0) * 
     $     ftau10 + (sinthetatau**2.d0) * ftau20)) +
     $     (g*geL/costhw)**2.d0 * ((costhetatau**2.d0) * fmtau1Mz +
     $     (sinthetab**2.d0) * fmtau2mz) +
     $     (g*g)*0.5d0 * fsnu1MW  +
     $     (g*g) * 0.25d0 * (1.d0 * ((costhetatau**2.d0) * a0mtau1 + 
     $     (sinthetatau**2.d0) * a0mtau2) + 2.d0 * a0snu1 ) -
     $     (g*g) * 0.5d0 * (- 1.5d0 * a0md1 - 1.5d0*a0ms1 +
     $     1.5d0 * a0mu1 + 1.5d0*a0mc1 +
     $     0.5d0 * (a0snu1 + a0snu2 + a0snu3) -
     $     0.5d0 * (a0me1 + a0mmu1 + ((costhetatau**2.d0) * a0mtau1 +
     $     (sinthetatau**2.d0) * a0mtau2))) +
     $     (gp**2.d0) * 0.25d0 *(yeL**2.d0) * 
     $     ((costhetatau**2.d0) * a0mtau1 + 
     $     (sinthetatau**2.d0) * a0mtau2) +
     $     (gp**2.d0) * 0.25d0 * yeL * ((3.d0*ydL*(a0md1 + a0ms1 + 
     $     ((costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2))) +
     $     (3.d0 * ydR *(a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2)) +
     $     (3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t1)) +
     $     (3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (costhetat**2.d0)*a0t2 + (sinthetat**2.d0)*a0t1)) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))

      
!----------------------------------------------------------------------------------

      higgsterm(2, 2) = higgsterm(2,2) +
     $     (4.d0*(g*g)/costhw**2.d0) * (geR**2.d0) * a0MZ + 
     $     ((1.d0)*g*g*sinsqthw) * ((sinthetatau)*ftau10 + 
     $     (costhetatau**2.d0)*ftau20) +
     $     (g*geR/costhw)**2.d0 * ((sinthetatau**2.d0) * fmtau1MZ + 
     $     (costhetatau**2.d0) * fmtau2MZ) +
     $     (gp**2.d0) * 0.25d0 * (yeR**2.d0) * 
     $     ((sinthetatau**2.d0) * a0mtau1 + 
     $     (costhetatau**2.d0) * a0mtau2) +
     $     (gp**2.d0) * 0.25d0 * yeR *  
     $     (3.d0 * ydL * (a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2) +
     $     3.d0 * ydR * (a0md2 +  a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2) +
     $     3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t2) +
     $     3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 + (costhetat**2.d0)*a0t2) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))

      
!----------------------------------------------------------------------
      
      higgsterm(1,2) = higgsterm(1,2) + 
     $     (gp**2.d0)* yeL*yeR * sinthetatau*costhetatau *
     $     (a0mtau1 - a0mtau2) +
     $     (1.d0 * g*g* sinsqthw) * sinthetatau*costhetatau * 
     $     (ftau10 - ftau20) -
     $     ((g*g)/costhw**2.d0) * geL*geR * sinthetatau*costhetatau *
     $     (fmtau1Mz - fmtau2Mz)

!----------------------------------------------------------------------
!     C       Chargino term
!----------------------------------------------------------------------
      
      call funcg(p,mchargino(1),0.d0,q,gmchargino1mtau)
      call funcg(p,mchargino(2),0.d0,q,gmchargino2mtau)
      call b0(p,mchargino(1),0.d0,q,b0mchargino1mtau)
      call b0(p,mchargino(2),0.d0,q,b0mchargino2mtau)
      
      chargino(1, 1) =  fChnustauLL(1)*gmchargino1mtau -
     $     gChnustauLL(1)*mchargino(1)*0.d0*b0mchargino1mtau*2.d0 +
     $     fChnustauLL(2)*gmchargino2mtau -
     $     gChnustauLL(2)*mchargino(2)*0.d0*b0mchargino2mtau*2.d0
      
      chargino(1, 2) =  fChnustauLR(1)*gmchargino1mtau -
     $     gChnustauLR(1)*mchargino(1)*0.d0*b0mchargino1mtau*2.d0 +
     $     fChnustauLR(2)*gmchargino2mtau -
     $     gChnustauLR(2)*mchargino(2)*0.d0*b0mchargino2mtau*2.d0
      
      chargino(2, 2) =  fChnustauRR(1)*gmchargino1mtau - 
     $     gChnustauRR(1)*mchargino(1)*0.d0*b0mchargino1mtau*2.d0 +
     $     fChnustauRR(2)*gmchargino2mtau - 
     $     gChnustauRR(2)*mchargino(2)*0.d0*b0mchargino2mtau*2.d0
      
!---------------------------------------------------------------------
!     C     neutralino terms
      
      call b0(p,nmneut(1),mtau,q,b0mneut1mtau)
      call b0(p,nmneut(2),mtau,q,b0mneut2mtau)
      call b0(p,nmneut(3),mtau,q,b0mneut3mtau)
      call b0(p,nmneut(4),mtau,q,b0mneut4mtau)
      
      call funcg(p,nmneut(1),mtau,q,gmneut1mtau)
      call funcg(p,nmneut(2),mtau,q,gmneut2mtau)
      call funcg(p,nmneut(3),mtau,q,gmneut3mtau)
      call funcg(p,nmneut(4),mtau,q,gmneut4mtau)
      
      neutralino(1, 1) =  
     $     fChi0taustauLL(1)*gmneut1mtau - gChi0taustauLL(1)*2.d0*
     $     mneut(1)*mtau*b0mneut1mtau +
     $     fChi0taustauLL(2)*gmneut2mtau - gChi0taustauLL(2)*2.d0*
     $     mneut(2)*mtau*b0mneut2mtau +
     $     fChi0taustauLL(3)*gmneut3mtau - gChi0taustauLL(3)*2.d0*
     $     mneut(3)*mtau*b0mneut3mtau +
     $     fChi0taustauLL(4)*gmneut4mtau - gChi0taustauLL(4)*2.d0*
     $     mneut(4)*mtau*b0mneut4mtau
      
      
      neutralino(2, 2) =  
     $     fChi0taustauRR(i)*gmneut1mtau - gChi0taustauRR(1)*2.d0*
     $     mneut(1)*mtau*b0mneut1mtau +
     $     fChi0taustauRR(2)*gmneut2mtau - gChi0taustauRR(2)*2.d0*
     $     mneut(2)*mtau*b0mneut2mtau +
     $     fChi0taustauRR(3)*gmneut3mtau - gChi0taustauRR(3)*2.d0*
     $     mneut(3)*mtau*b0mneut3mtau +
     $     fChi0taustauRR(4)*gmneut4mtau - gChi0taustauRR(4)*2.d0*
     $     mneut(4)*mtau*b0mneut4mtau
      
      neutralino(1, 2) = 
     $     fChi0taustauLR(1)*gmneut1mtau - gChi0taustauLR(1)*2.d0*
     $     mneut(1)*mtau*b0mneut1mtau +
     $     fChi0taustauLR(2)*gmneut2mtau - gChi0taustauLR(2)*2.d0*
     $     mneut(2)*mtau*b0mneut2mtau +
     $     fChi0taustauLR(3)*gmneut3mtau - gChi0taustauLR(3)*2.d0*
     $     mneut(3)*mtau*b0mneut3mtau +
     $     fChi0taustauLR(4)*gmneut4mtau - gChi0taustauLR(4)*2.d0*
     $     mneut(4)*mtau*b0mneut4mtau 
      
!-------------------------------------------------------------------

      MSLD3(1,1) = mSLRG(3,3) + mtau**2.d0 + geL*MZ*MZ*dcos(2.d0*beta)
      MSLD3(1,2) = mtau*((AERG(3,3) - sgnmu*modmu*dtan(beta)))
      MSLD3(2,1) = MSLD3(1,2)
      MSLD3(2,2) = mSERG(3,3) + mtau**2.d0 + geR*MZ*MZ*dcos(2.d0*beta)

!----------------------------------------------------------------------

      pistau(1,1) = (cstop(1,1) + csbtm(1,1) + gtterm(1,1) +
     $     chargino(1,1) + neutralino(1,1) + higgsterm(1,1))/
     $     (16.d0*pi*pi)

      pistau(1,2) = (cstop(1,2) + csbtm(1,2) + gtterm(1,2)  +
     $     chargino(1,2) + neutralino(1,2) + higgsterm(1,2))/
     $     (16.d0*pi*pi)

      pistau(2,2) = (cstop(2,2) + csbtm(2,2) + gtterm(2,2) +
     $     chargino(2,2) + neutralino(2,2) + higgsterm(2,2))/
     $     (16.d0*pi*pi)
      
      pistau(2,1) = pistau(1,2)



      delpistau(1,1) = pistau(1,1)
      delpistau(1,2) = pistau(1,2)
      delpistau(2,1) = pistau(2,1)
      delpistau(2,2) = pistau(2,2)

      pistau(1,1) = MSLD3(1,1) - pistau(1,1)
      pistau(1,2) = MSLD3(1,2) - pistau(1,2)
      pistau(2,2) = MSLD3(2,2) - pistau(2,2)
      pistau(2,1) = pistau(1,2) 


!----------------------------------------------------------------------
C     find the singular values and the diagonalising matrices
!----------------------------------------------------------------------

      info  = 10
      AOK = 0
      
C     call dsyev('V','U',2,pistau,2,STaueg,work,lwork,info)
      
      if(info.eq.0) then
         AOK = AOK + 1
      endif
  
      Call CEigensystem(2,pistau,2,STaueg,staumix,2,1)
      
      RETURN

      END SUBROUTINE pistaul

C==================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     
C     pismul is checked on 28/05/2010 @ 12:30.
C     
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      SUBROUTINE pisMUl(p,q,g,gp,tanbeta,mSLRG,mSERG,
     $     yeRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,M3t,SMueg)
 
      IMPLICIT NONE 

      INTEGER i,j,info,AOK,lwork
      parameter(lwork=35)
      DOUBLE PRECISION work(lwork)
      DOUBLE PRECISION tanbeta,mSLRG(3,3),mSERG(3,3)
      DOUBLE PRECISION AERG(3,3),SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yeRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3),nmneut(4)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g,M3t
   
      double precision vev1,beta
 
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION costhetamu,sinthetamu

      DOUBLE PRECISION thetau,thetad,thetae

      DOUBLE PRECISION pisMu(2,2)
      data pisMu/ 4 * 0.d0/
 
      DOUBLE PRECISION higgsterm(2,2)
      data higgsterm/ 4 * 0.d0/

      DOUBLE PRECISION dnu(4),dnd(4),cn(4)

      DOUBLE PRECISION lsmuLsmuLR(4,2),lsmuLsmu12(4,2)
      data lsmuLsmuLR/ 8 * 0.d0/,lsmuLsmu12/ 8 * 0.d0/
      
      DOUBLE PRECISION lsmuRsmuLR(4,2),lsmuRsmu12(4,2)
      data lsmuRsmuLR/ 8 * 0.d0/,lsmuRsmu12/ 8 * 0.d0/

      
      DOUBLE PRECISION aPsi0smur(4), bPsi0smur(4), aPsi0smul(4)
      data aPsi0smur/ 4 * 0.d0/, bPsi0smur/ 4 * 0.d0/,
     $     aPsi0smul/ 4 * 0.d0/


      DOUBLE PRECISION bPsi0smul(4)
      data bPsi0smul/ 4 * 0.d0/
      
      DOUBLE PRECISION aChi0smul(4), bChi0smul(4), aChi0smur(4)
      data aChi0smul/ 4 * 0.d0/, bChi0smul/ 4 * 0.d0/, 
     $     aChi0smur/ 4 * 0.d0/
      
      DOUBLE PRECISION bChi0smur(4)
      data bChi0smur/ 4 * 0.d0/

      
      DOUBLE PRECISION gChi0musmuLL(4), fChi0musmuLL(4)
      data gChi0musmuLL/ 4 * 0.d0/, fChi0musmuLL/ 4 * 0.d0/

      DOUBLE PRECISION gChi0musmuLR(4), fChi0musmuLR(4)
      data gChi0musmuLR/ 4 * 0.d0/, fChi0musmuLR/ 4 * 0.d0/


      DOUBLE PRECISION gChi0musmuRR(4), fChi0musmuRR(4)
      data gChi0musmuRR/ 4 * 0.d0/, fChi0musmuRR/ 4 * 0.d0/

      DOUBLE PRECISION bPsicsmul(2), bPsicsmur(2), aPsicsmul(2)
      data bPsicsmul/ 2 * 0.d0/,bPsicsmur/ 2 * 0.d0/,
     $     aPsicsmul/ 2 * 0.d0/
      
      DOUBLE PRECISION aPsicsmur(2)
      data aPsicsmur/ 2 * 0.d0/

      DOUBLE PRECISION aChicsmur(2), aChicsmul(2)
      data aChicsmur/ 2 * 0.d0/,
     $     aChicsmul/ 2 * 0.d0/

      DOUBLE PRECISION bChicsmul(2),bChicsmur(2)
      data bChicsmul/ 2 * 0.d0/,bChicsmur/ 2 * 0.d0/

      DOUBLE PRECISION aChsmu(2,2),bChsmu(2,2)
      data aChsmu/ 4 * 0.d0/,bChsmu/ 4 * 0.d0/

      DOUBLE PRECISION gtterm(2,2), cstop(2,2), csbtm(2,2) 
      data gtterm/ 4 * 0.d0/,cstop/ 4 * 0.d0/,csbtm/ 4 * 0.d0/

      DOUBLE PRECISION  chargino(2,2), neutralino(2,2), mSLD3(2,2)
      data chargino/ 4 * 0.d0/,neutralino/ 4 * 0.d0/, mSLD3/ 4 * 0.d0/


      DOUBLE PRECISION fChsmuLL(2), gChsmuLL(2) 
      data fChsmuLL/ 2 * 0.d0/,gChsmuLL/ 2 * 0.d0/

      DOUBLE PRECISION fChsmuLR(2), gChsmuLR(2) 
      data fChsmuLR/ 2 * 0.d0/,gChsmuLR/ 2 * 0.d0/

      DOUBLE PRECISION fChsmuRR(2), gChsmuRR(2)
      data fChsmuRR/ 2 * 0.d0/,gChsmuRR/ 2 * 0.d0/

      DOUBLE PRECISION intm1(2),intm2(2)
      data intm1/ 2 * 0.d0/, intm2/ 2 * 0.d0/

      DOUBLE PRECISION lHcsmursnulr(2,2),lHcsmursnu12(2,2)
      data lHcsmursnulr/ 4 * 0.d0/, lHcsmursnu12/ 4 * 0.d0/

      DOUBLE PRECISION lHcsmulsnulr(2,2),lHcsmulsnu12(2,2)
      data lHcsmulsnulr/ 4 * 0.d0/, lHcsmulsnu12/ 4 * 0.d0/

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION lHsmuLsmu12(4,2),lHsmuRsmu12(4,2)
      data lHsmuLsmu12/ 8 * 0.d0/, lHsmuRsmu12/ 8 * 0.d0/

!-----------------------------------------------------------------------------

      DOUBLE PRECISION rthetamu(2,2),ralpha(2,2)

      DOUBLE PRECISION ggmu,fmu10,fmu20,a0b1,a0b2,b0M3mmu,a0t1,a0t2
      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz
      DOUBLE PRECISION b0h0msmu(4,2),b0hcmstop(2,2)

      DOUBLE PRECISION a0md1,a0md2,a0ms1,a0ms2,a0mtau1,a0mtau2,a0snu2
      DOUBLE PRECISION a0me1,a0me2,a0mmu1,a0mmu2,a0snu1,a0snu3
      DOUBLE PRECISION a0mc1,a0mc2,a0mu1,a0mu2
      
      DOUBLE PRECISION fmmu1MZ,fmmu2MZ,f0MW,fsnu2MW
      DOUBLE PRECISION gmchargino1mmu,gmchargino2mmu
      DOUBLE PRECISION b0mchargino1mmu,b0mchargino2mmu,b0mneut4mmu
      DOUBLE PRECISION b0mneut1mmu,b0mneut2mmu,b0mneut3mmu
      DOUBLE PRECISION gmneut1mmu,gmneut2mmu,gmneut3mmu,gmneut4mmu
      
      DOUBLE PRECISION SMueg(2),smumix(2,2)

      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2
      DOUBLE PRECISION mmu1,mmu2
  
      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae


      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d
  

      include 'stdinputs.h'

C--------------------------

      pi = 4.d0 * datan(1.d0)

      MW = MWpole
      MZ = MZpole

      SMueg(1) = 0.d0
      SMueg(2) = 0.d0

      
      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw = dsqrt(sinsqthw)
      thw = dasin(sinthw)
      costhw = dcos(thw)
      cos2thw = dcos(2.d0*thw)


      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))
      
      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))
       
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
!------------------------

      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)
      costhetamu = dcos(thetamu)
      sinthetamu = dsin(thetamu)

!----------------------------------------

      mt1 = mtL
      mt2 = mtR

      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR

      mmu1 = mmuL
      mmu2 = mmuR

C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta     = datan(tanbeta) 
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha  = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta  = dcos(datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)

!-----------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0

!----------------------------------------------------------
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

!-----------------------------------------------------------    

      lsmuLsmuLR(1,1) = (g*MZ*geL*cosbeta/costhw) + 
     $     (yeRG(2,2)**2.d0) * vev1
      lsmuLsmuLR(1,2) = (yeRG(2,2)*AERG(2,2))/dsqrt(2.d0)

      lsmuLsmuLR(2,1) = ((- g*MZ*geL*sinbeta)/costhw)
      lsmuLsmuLR(2,2) = -1.d0 * (yeRG(2,2)*sgnmu*modmu)/dsqrt(2.d0)

      lsmuLsmuLR(3,1) = 0.d0
      lsmuLsmuLR(3,2) = (-1.d0/dsqrt(2.d0)) *
     $     (-1.d0 * sgnmu*modmu*sinbeta*yeRG(2,2) + 
     $     yeRG(2,2)*AERG(2,2)*cosbeta)

      lsmuLsmuLR(4,1) = 0.d0
      lsmuLsmuLR(4,2) = (-1.d0/dsqrt(2.d0)) * 
     $     (-1.d0 * sgnmu*modmu*cosbeta*yeRG(2,2) - 
     $     yeRG(2,2)*AERG(2,2)*sinbeta)



      lsmuRsmuLR(1,1) = lsmuLsmuLR(1,2)
      lsmuRsmuLR(1,2) = (g*MZ*geR*cosbeta/costhw) + 
     $     yeRG(2,2)**2.d0 * vev1

      lsmuRsmuLR(2,1) = lsmuLsmuLR(2,2)
      lsmuRsmuLR(2,2) = -g*MZ*geR*sinbeta/costhw

      lsmuRsmuLR(3,1) = -lsmuLsmuLR(3,2)
      lsmuRsmuLR(3,2) = 0.d0

      lsmuRsmuLR(4,1) = -lsmuLsmuLR(4,2)
      lsmuRsmuLR(4,2) = 0.d0

!--------------------------------------------------------------------------------  
  
      loopmixst: DO i = 1, 4

      intm1(1) = lsmuLsmuLR(i,1)
      intm1(2) = lsmuLsmuLR(i,2)

      call rmat2d(thetamu,rthetamu)
            
      intm2(1) = rthetamu(1,1)*intm1(1) + rthetamu(1,2)*intm1(2)
      intm2(2) = rthetamu(2,1)*intm1(1) + rthetamu(2,2)*intm1(2)

      lsmuLsmu12(i,1) = intm2(1)
      lsmuLsmu12(i,2) = intm2(2)

      intm1(1) = lsmuRsmuLR(i,1)
      intm1(2) = lsmuRsmuLR(i,2)

      intm2(1) = rthetamu(1,1)*intm1(1) + rthetamu(1,2)*intm1(2)
      intm2(2) = rthetamu(2,1)*intm1(1) + rthetamu(2,2)*intm1(2)

      lsmuRsmu12(i,1) = intm2(1)
      lsmuRsmu12(i,2) = intm2(2)

      ENDDO loopmixst
!--------------------------------------------------------------------------
C       Mixing CP-even Higgs 

      loopmixcpeven: DO i = 1, 2
      
      intm1(1) = lsmuLsmu12(1,i)
      intm1(2) = lsmuLsmu12(2,i)

      call rmat2d(alpha,ralpha)
            
      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lHsmuLsmu12(1,i) = intm2(1)
      lHsmuLsmu12(2,i) = intm2(2)

      intm1(1) = lsmuRsmu12(1,i)
      intm1(2) = lsmuRsmu12(2,i)

      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lHsmuRsmu12(1,i) = intm2(1)
      lHsmuRsmu12(2,i) = intm2(2)
     
      ENDDO loopmixcpeven
!----------------------------------------------------------------------
C    Feynman rules - charged higgs
!<--------- (H+ G+, L R) basis
      
      lHcsmulsnulr(1, 1) = ( g*MW*sin2beta - 
     $     (yeRG(2,2)**2.d0) * vev1 * sinbeta)/dsqrt(2.d0)

      lHcsmulsnulr(1,2) = 0.d0 

      lHcsmulsnulr(2,1) = (-g*MW*cos2beta - 
     $     yeRG(2,2)**2.d0 * vev1 * cosbeta)/dsqrt(2.d0)

      lHcsmulsnulr(2,2) = 0.d0


!----------------------------------------

      intm1(1) = lHcsmulsnulr(1,1)
      intm1(2) = lHcsmulsnulr(1,2)

      call rmat2d(thetamu,rthetamu)
            
      intm2(1) = rthetamu(1,1)*intm1(1) + rthetamu(1,2)*intm1(2)
      intm2(2) = rthetamu(2,1)*intm1(1) + rthetamu(2,2)*intm1(2)

 
      lHcsmulsnu12(1, 1) = intm2(1)
      lHcsmulsnu12(1, 2) = intm2(2)

      intm1(1) = lHcsmulsnulr(2, 1)
      intm1(2) = lHcsmulsnulr(2, 2)


                  
      intm2(1) = rthetamu(1,1)*intm1(1) + rthetamu(1,2)*intm1(2)
      intm2(2) = rthetamu(2,1)*intm1(1) + rthetamu(2,2)*intm1(2)

      lHcsmulsnu12(2, 1) = intm2(1)
      lHcsmulsnu12(2, 2) = intm2(2)

!---------------------------------------------------------------------------
!<------------ (H+ G+, L R) basis

      lHcsmursnulr(1, 1) = yeRG(2,2)*(-1.d0 * sgnmu*modmu*cosbeta - 
     $     AERG(2,2)*sinbeta)

      lHcsmursnulr(1, 2) = 0.d0

      lHcsmursnulr(2, 2) = yeRG(2,2)*(-1.d0 * sgnmu*modmu*sinbeta + 
     $     AERG(2,2)*cosbeta)  

      lHcsmursnulr(2, 1) = 0.d0

!-------------------------------------------------------------------------

      intm1(1) = lHcsmursnulr(1,1)
      intm1(2) = lHcsmursnulr(1,2)

      call rmat2d(thetamu,rthetamu)
            
      intm2(1) = rthetamu(1,1)*intm1(1) + rthetamu(1,2)*intm1(2)
      intm2(2) = rthetamu(2,1)*intm1(1) + rthetamu(2,2)*intm1(2)

      lHcsmursnu12(1,1) = intm2(1)
      lHcsmursnu12(1,2) = intm2(2)

      intm1(1) = lHcsmursnulr(2,1)
      intm1(2) = lHcsmursnulr(2,2)

                 
      intm2(1) = rthetamu(1,1)*intm1(1) + rthetamu(1,2)*intm1(2)
      intm2(2) = rthetamu(2,1)*intm1(1) + rthetamu(2,2)*intm1(2)
  
      lHcsmursnu12(2, 1) = intm2(1)
      lHcsmursnu12(2, 2) = intm2(2)

!------------------------------------------------------------------------
C     Feynman rules - neutralino
!------------------------------------------------------------------------

      aPsi0smur(1) = yeR*gp/(dsqrt(2.d0))
      aPsi0smur(2) = 0.d0
      aPsi0smur(3) = 0.d0
      aPsi0smur(4) = 0.d0

      bPsi0smul(1) = yeL*gp/(dsqrt(2.d0))
      bPsi0smul(2) = g/dsqrt(2.d0)
      bPsi0smul(3) = 0.d0
      bPsi0smul(4) = 0.d0

      aPsi0smul(1) = 0.d0
      aPsi0smul(2) = 0.d0
      aPsi0smul(3) = 0.d0 
      aPsi0smul(4) = yeRG(2,2)

      bPsi0smur(1) = 0.d0
      bPsi0smur(2) = 0.d0
      bPsi0smur(3) = 0.d0
      bPsi0smur(4) = yeRG(2,2)


      aChi0smul(1) = ON(1,1)*aPsi0smul(1) + ON(1,2)*aPsi0smul(2) +
     $                 ON(1,3)*aPsi0smul(1) + ON(1,4)*aPsi0smul(4) 
      aChi0smul(2) = ON(2,1)*aPsi0smul(1) + ON(2,2)*aPsi0smul(2) +
     $                 ON(2,3)*aPsi0smul(1) + ON(2,4)*aPsi0smul(4) 
      aChi0smul(3) = ON(3,1)*aPsi0smul(1) + ON(3,2)*aPsi0smul(2) +
     $                 ON(3,3)*aPsi0smul(1) + ON(3,4)*aPsi0smul(4) 
      aChi0smul(4) = ON(4,1)*aPsi0smul(1) + ON(4,2)*aPsi0smul(2) +
     $                 ON(4,3)*aPsi0smul(1) + ON(4,4)*aPsi0smul(4) 


      bChi0smul(1) = ON(1,1)*bPsi0smul(1) + ON(1,2)*bPsi0smul(2) +
     $                 ON(1,3)*bPsi0smul(1) + ON(1,4)*bPsi0smul(4) 
      bChi0smul(2) = ON(2,1)*bPsi0smul(1) + ON(2,2)*bPsi0smul(2) +
     $                 ON(2,3)*bPsi0smul(1) + ON(2,4)*bPsi0smul(4) 
      bChi0smul(3) = ON(3,1)*bPsi0smul(1) + ON(3,2)*bPsi0smul(2) +
     $                 ON(3,3)*bPsi0smul(1) + ON(3,4)*bPsi0smul(4) 
      bChi0smul(4) = ON(4,1)*bPsi0smul(1) + ON(4,2)*bPsi0smul(2) +
     $                 ON(4,3)*bPsi0smul(1) + ON(4,4)*bPsi0smul(4) 




      aChi0smur(1) = ON(1,1)*aPsi0smur(1) + ON(1,2)*aPsi0smur(2) +
     $                 ON(1,3)*aPsi0smur(3) + ON(1,4)*aPsi0smur(4)
      aChi0smur(2) = ON(2,1)*aPsi0smur(1) + ON(2,2)*aPsi0smur(2) +
     $                 ON(2,3)*aPsi0smur(3) + ON(2,4)*aPsi0smur(4)
      aChi0smur(3) = ON(3,1)*aPsi0smur(1) + ON(3,2)*aPsi0smur(2) +
     $                 ON(3,3)*aPsi0smur(3) + ON(3,4)*aPsi0smur(4)
      aChi0smur(4) = ON(4,1)*aPsi0smur(1) + ON(4,2)*aPsi0smur(2) +
     $                 ON(4,3)*aPsi0smur(3) + ON(4,4)*aPsi0smur(4)


      bChi0smur(1) = ON(1,1)*bPsi0smur(1) + ON(1,2)*bPsi0smur(2) +
     $                 ON(1,3)*bPsi0smur(3) + ON(1,4)*bPsi0smur(4) 
      bChi0smur(2) = ON(2,1)*bPsi0smur(1) + ON(2,2)*bPsi0smur(2) +
     $                 ON(2,3)*bPsi0smur(3) + ON(2,4)*bPsi0smur(4) 
      bChi0smur(3) = ON(3,1)*bPsi0smur(1) + ON(3,2)*bPsi0smur(2) +
     $                 ON(3,3)*bPsi0smur(3) + ON(3,4)*bPsi0smur(4) 
      bChi0smur(4) = ON(4,1)*bPsi0smur(1) + ON(4,2)*bPsi0smur(2) +
     $                 ON(4,3)*bPsi0smur(3) + ON(4,4)*bPsi0smur(4) 


      loopchtllrr: DO i = 1, 4

      fChi0musmuLL(i) = (aChi0smul(i)*aChi0smul(i) + 
     $     bChi0smul(i)*bChi0smul(i))
      gChi0musmuLL(i) = (bChi0smul(i)*aChi0smul(i) + 
     $     bChi0smul(i)*aChi0smul(i))
      fChi0musmuRR(i) = (aChi0smur(i)*aChi0smur(i) + 
     $     bChi0smur(i)*bChi0smur(i))
      gChi0musmuRR(i) = (bChi0smur(i)*aChi0smur(i) + 
     $     bChi0smur(i)*aChi0smur(i))
      fChi0musmuLR(i) = (aChi0smur(i)*aChi0smul(i) + 
     $     bChi0smur(i)*bChi0smul(i))
      gChi0musmuLR(i) = (bChi0smul(i)*aChi0smur(i) + 
     $     bChi0smur(i)*aChi0smul(i))
      
      ENDDO loopchtllrr

!---------------------------------------------------------------------------
c     Feynman Rules for chargino

      aPsicsmul(1) = 0.d0
      aPsicsmul(2) = 0.d0

      aPsicsmur(1) = 0.d0
      aPsicsmur(2) = 0.d0

      bPsicsmul(1) = g
      bPsicsmul(2) = 0.d0

      bPsicsmur(1) = 0.d0
      bPsicsmur(2) = - yeRG(2,2)

      aChicsmul(1) = OCR(1,1)*aPsicsmul(1) + OCR(1,2)*aPsicsmul(2)
      aChicsmul(2) = OCR(2,1)*aPsicsmul(1) + OCR(2,2)*aPsicsmul(2)

      bChicsmul(1) = OCL(1,1)*bPsicsmul(1) + OCL(1,2)*bPsicsmul(2)
      bChicsmul(2) = OCL(2,1)*bPsicsmul(1) + OCL(2,2)*bPsicsmul(2)

      aChicsmur(1) = OCR(1,1)*aPsicsmur(1) + OCR(1,2)*aPsicsmur(2)
      aChicsmur(2) = OCR(2,1)*aPsicsmur(1) + OCR(2,2)*aPsicsmur(2)

      bChicsmur(1) = OCL(1,1)*bPsicsmur(1) + OCL(1,2)*bPsicsmur(2)
      bChicsmur(2) = OCL(2,1)*bPsicsmur(1) + OCL(2,2)*bPsicsmur(2)

      
      loopchst: DO i = 1, 2

      fChsmuLL(i) = (aChicsmul(i)*aChicsmul(i) +
     $     bChicsmul(i)*bChicsmul(i))
      gChsmuLL(i) = (bChicsmul(i)*aChicsmul(i) +
     $     aChicsmul(i)*bChicsmul(i))
      fChsmuLR(i) = (aChicsmul(i)*aChicsmur(i) +
     $     bChicsmul(i)*bChicsmur(i))
      gChsmuLR(i) = (bChicsmul(i)*aChicsmur(i) +
     $     aChicsmul(i)*bChicsmur(i))
      fChsmuRR(i) = (aChicsmur(i)*aChicsmur(i) +
     $     bChicsmur(i)*bChicsmur(i))
      gChsmuRR(i) = (bChicsmur(i)*aChicsmur(i) +
     $     aChicsmur(i)*bChicsmur(i))

      ENDDO loopchst

!--------------------------------------------------------------------
      
C     Corrections Begin
!------------------------------------------------------------------

      call funcg(p,M3t,mmu,q,ggmu)

!     print*,"ggt"
      
      call f(p,mmu1,0.d0,q,fmu10)
      call f(p,mmu2,0.d0,q,fmu20)

      call a0(muR,q,a0mu2)
      call a0(muL,q,a0mu1)
      call a0(mcR,q,a0mc2)
      call a0(mcL,q,a0mc1)
      call a0(mt2,q,a0t2)
      call a0(mt1,q,a0t1)

      call a0(msR,q,a0ms2)
      call a0(msL,q,a0ms1)
      call a0(mdR,q,a0md2)
      call a0(mdL,q,a0md1)
      call a0(mb2,q,a0b2)
      call a0(mb1,q,a0b1)

      call a0(meR,q,a0me2)
      call a0(meL,q,a0me1)
      call a0(mmu2,q,a0mmu2)
      call a0(mmu1,q,a0mmu1)
      call a0(mtau2,q,a0mtau2)
      call a0(mtau1,q,a0mtau1)

      call a0(snu(1),q,a0snu1)
      call a0(snu(2),q,a0snu2)
      call a0(snu(3),q,a0snu3)


      call a0(mA0,q,a0mA)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)
      call a0(mHpm,q,a0mHpm)
      call a0(MW,q,a0Mw)
      call a0(MZ,q,a0Mz)

      call f(p,0.d0,MW,q,f0Mw)
      call f(p,0.d0,MW,q,f0Mw)
      call f(p,mmu1,MZ,q,fmmu1Mz)
      call f(p,mmu2,MZ,q,fmmu2Mz)
      
      call f(p,snu(2),MW,q,fsnu2MW)

      call b0(p,M3t,mmu,q,b0M3mmu)
    
!--------------------------------------------------------------
      gtterm(1, 1) = 0.d0

      gtterm(2, 2) = 0.d0

      gtterm(1, 2) = 0.d0
!----------------------------------------------------------------

      csbtm(1, 1) = yeRG(2,2)**2.d0 * ((sinthetamu)**2.d0*a0mmu1 + 
     $     (costhetamu)**2.d0*a0mmu2)

      csbtm(2, 2) = yeRG(2,2)**2.d0 * ((costhetamu)**2.d0*a0mmu1 + 
     $     (sinthetamu)**2.d0*a0mmu2)

      csbtm(1, 2) = yeRG(2,2)**2.d0 * costhetamu*sinthetamu * 1.d0 *
     $     (a0mmu1 - a0mmu2)

!------------------------------------------------------------------
      cstop(1, 1) = 0.d0

      cstop(2, 2) = yeRG(2,2)**2.d0 * a0snu2
!-------------------------------------


      higgsterm(1,1) = 0.5d0 * ((yeRG(2,2)**2.d0) * dnd(1) -
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((yeRG(2,2)**2.d0) * dnd(2) - 
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((yeRG(2,2)**2.d0) * dnd(3) - 
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((yeRG(2,2)**2.d0) * dnd(4) - 
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA
      
      
      higgsterm(1,2) = 0.d0
      

      higgsterm(2,2) = 0.5d0*(yeRG(2,2)**2.d0 * dnd(1) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0*(yeRG(2,2)**2.d0 * dnd(2) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0*(yeRG(2,2)**2.d0 * dnd(3) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(3)) * a0Mz +
     $     0.5d0*(yeRG(2,2)**2.d0 * dnd(4) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA


      higgsterm(1,1) = higgsterm(1,1) + 
     $     (0.d0 + ((g*g) * 
     $     ((geL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(3))) * a0Mhpm +
     $     (0.d0 + ((g*g) * 
     $     ((geL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(4))) * a0MW

      higgsterm(2,2) = higgsterm(2,2) + 
     $     ((yeRG(2,2)**2.d0) * dnu(3) +  
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(3)) * a0Mhpm +
     $     ((yeRG(2,2)**2.d0) * dnu(4) + 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(4)) * a0MW !<----------------chkd

     
!--------------------------
      
      call b0(p,mHu0,mmu1,q,b0h0msmu(1,1))
      call b0(p,mHu0,mmu2,q,b0h0msmu(1,2))
      call b0(p,mh0,mmu1,q,b0h0msmu(2,1))
      call b0(p,mh0,mmu2,q,b0h0msmu(2,2))
      call b0(p,MZ,mmu1,q,b0h0msmu(3,1))
      call b0(p,MZ,mmu2,q,b0h0msmu(3,2))
      call b0(p,mA0,mmu1,q,b0h0msmu(4,1))
      call b0(p,mA0,mmu2,q,b0h0msmu(4,2))
      
      
      loophi: DO i = 1, 4        !<-------------neutral higgs terms 
      loophj: DO j = 1, 2
      
      higgsterm(1,1) = higgsterm(1,1) + 
     $     (lHsmuLsmu12(i,j)**2.d0) * b0h0msmu(i,j)
      
      higgsterm(1,2) = higgsterm(1,2) + 
     $     lHsmuLsmu12(i,j) * lHsmuRsmu12(i,j) * b0h0msmu(i,j)
      
      higgsterm(2,2) = higgsterm(2,2) + 
     $     (lHsmuRsmu12(i,j)**2.d0) * b0h0msmu(i,j)
      
      ENDDO loophj
      ENDDO loophi
!-------------------------------------
      
      call b0(p,0.d0,mHpm,q,b0hcmstop(1,1))
      call b0(p,0.d0,mHpm,q,b0hcmstop(1,2))
      call b0(p,0.d0,MW,q,b0hcmstop(2,1))
      call b0(p,0.d0,MW,q,b0hcmstop(2,2))
      
      
      loophci: DO i = 1, 2      !<-------------charged higgs terms
      loophcj: DO j = 1, 2
      
      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (lHcsmulsnu12(i,j)**2.d0) * b0hcmstop(i,j)       
      
      higgsterm(1, 2) = higgsterm(1,2) + 
     $     lHcsmulsnu12(i, j)*lHcsmursnu12(i, j)*b0hcmstop(i,j)
      
      higgsterm(2, 2) = higgsterm(2, 2) + 
     $     (lHcsmursnu12(i,j)**2.d0) * b0hcmstop(i,j)
      
      ENDDO loophcj
      ENDDO loophci

!---------------------

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (4.d0*(g*g)/costhw**2.d0)*(geL*geL)*a0MZ  + 
     $     2.d0*(g*g)*a0MW + 
     $     ((1.d0 *g*g*sinsqthw) * ((costhetamu**2.d0) * 
     $     fmu10 + (sinthetamu**2.d0) * fmu20)) +
     $     (g*geL/costhw)**2.d0 * ((costhetamu**2.d0) * fmmu1Mz +
     $     (sinthetamu**2.d0) * fmmu2mz) +
     $     (g*g)*0.5d0 * fsnu2MW  +
     $     (g*g) * 0.25d0 * (1.d0 * ((costhetamu**2.d0) * a0mmu1 + 
     $     (sinthetamu**2.d0) * a0mmu2) + 2.d0 * a0snu2 ) -
     $     (g*g) * 0.5d0 * (- 1.5d0 * a0md1 - 1.5d0*a0ms1 +
     $     1.5d0 * a0mu1 + 1.5d0*a0mc1 +
     $     0.5d0 * (a0snu1 + a0snu2 + a0snu3) -
     $     0.5d0 * (a0me1 + a0mmu1 + ((costhetatau**2.d0) * a0mtau1 +
     $     (sinthetatau**2.d0) * a0mtau2))) +
     $     (gp**2.d0) * 0.25d0 *(yeL**2.d0) * 
     $     ((costhetamu**2.d0) * a0mmu1 + 
     $     (sinthetamu**2.d0) * a0mmu2) +
     $     (gp**2.d0) * 0.25d0 * yeL * ((3.d0*ydL*(a0md1 + a0ms1 + 
     $     ((costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2))) +
     $     (3.d0 * ydR *(a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2)) +
     $     (3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t1)) +
     $     (3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (costhetat**2.d0)*a0t2 + (sinthetat**2.d0)*a0t1)) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))


      higgsterm(2, 2) = higgsterm(2,2) +
     $     (4.d0*(g*g)/costhw**2.d0) * (geR**2.d0) * a0MZ + 
     $     ((1.d0)*g*g*sinsqthw) * ((sinthetamu)*fmu10 + 
     $     (costhetamu**2.d0)*fmu20) +
     $     (g*geR/costhw)**2.d0 * ((sinthetamu**2.d0) * fmmu1MZ + 
     $     (costhetamu**2.d0) * fmmu2MZ) +
     $     (gp**2.d0) * 0.25d0 * (yeR**2.d0) * 
     $     ((sinthetamu**2.d0) * a0mmu1 + 
     $     (costhetamu**2.d0) * a0mmu2) +
     $     (gp**2.d0) * 0.25d0 * yeR *  
     $     (3.d0 * ydL * (a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2) +
     $     3.d0 * ydR * (a0md2 +  a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2) +
     $     3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t2) +
     $     3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 + (costhetat**2.d0)*a0t2) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))

      
!----------------------------------------------
      
      higgsterm(1,2) = higgsterm(1,2) + 
     $     (gp**2.d0) * yeL*yeR * sinthetamu*costhetamu *
     $     (a0mmu1 - a0mmu2) +
     $     (1.d0 * g*g* sinsqthw) * sinthetamu*costhetamu * 
     $     (fmu10 - fmu20) -
     $     ((g*g)/costhw**2.d0) * geL*geR * sinthetamu*costhetamu *
     $     (fmmu1Mz - fmmu2Mz)    

!----------------------------------------------------------------------
!     Chargino term
!----------------------------------------------------------------------
      
      call funcg(p,mchargino(1),0.d0,q,gmchargino1mmu)
      call funcg(p,mchargino(2),0.d0,q,gmchargino2mmu)
      call b0(p,mchargino(1),0.d0,q,b0mchargino1mmu)
      call b0(p,mchargino(2),0.d0,q,b0mchargino2mmu)
      
      chargino(1, 1) =  fChsmuLL(1)*gmchargino1mmu -
     $     gChsmuLL(1)*mchargino(1)*0.d0*b0mchargino1mmu*2.d0 +
     $     fChsmuLL(2)*gmchargino2mmu -
     $     gChsmuLL(2)*mchargino(2)*0.d0*b0mchargino2mmu*2.d0
      
      chargino(1, 2) =  fChsmuLR(1)*gmchargino1mmu -
     $     gChsmuLR(1)*mchargino(1)*0.d0*b0mchargino1mmu*2.d0 +
     $     fChsmuLR(2)*gmchargino2mmu -
     $     gChsmuLR(2)*mchargino(2)*0.d0*b0mchargino2mmu*2.d0
      
      chargino(2, 2) =  fChsmuRR(1)*gmchargino1mmu - 
     $     gChsmuRR(1)*mchargino(1)*0.d0*b0mchargino1mmu*2.d0 +
     $     fChsmuRR(2)*gmchargino2mmu - 
     $     gChsmuRR(2)*mchargino(2)*0.d0*b0mchargino2mmu*2.d0
      
!---------------------------------------------------------------------
!     neutralino terms
!---------------------------------------------------------------------
      
      call b0(p,nmneut(1),mmu,q,b0mneut1mmu)
      call b0(p,nmneut(2),mmu,q,b0mneut2mmu)
      call b0(p,nmneut(3),mmu,q,b0mneut3mmu)
      call b0(p,nmneut(4),mmu,q,b0mneut4mmu)
      
      call funcg(p,nmneut(1),mmu,q,gmneut1mmu)
      call funcg(p,nmneut(2),mmu,q,gmneut2mmu)
      call funcg(p,nmneut(3),mmu,q,gmneut3mmu)
      call funcg(p,nmneut(4),mmu,q,gmneut4mmu)
      
      neutralino(1, 1) = 
     $     fChi0musmuLL(1)*gmneut1mmu - gChi0musmuLL(1)*2.d0*
     $     mneut(1)*mmu*b0mneut1mmu +
     $     fChi0musmuLL(2)*gmneut2mmu - gChi0musmuLL(2)*2.d0*
     $     mneut(2)*mmu*b0mneut2mmu +
     $     fChi0musmuLL(3)*gmneut3mmu - gChi0musmuLL(3)*2.d0*
     $     mneut(3)*mmu*b0mneut3mmu +
     $     fChi0musmuLL(4)*gmneut4mmu - gChi0musmuLL(4)*2.d0*
     $     mneut(4)*mmu*b0mneut4mmu
      
      
      neutralino(2, 2) = 
     $     fChi0musmuRR(i)*gmneut1mmu - gChi0musmuRR(1)*2.d0*
     $     mneut(1)*mmu*b0mneut1mmu +
     $     fChi0musmuRR(2)*gmneut2mmu - gChi0musmuRR(2)*2.d0*
     $     mneut(2)*mmu*b0mneut2mmu +
     $     fChi0musmuRR(3)*gmneut3mmu - gChi0musmuRR(3)*2.d0*
     $     mneut(3)*mmu*b0mneut3mmu +
     $     fChi0musmuRR(4)*gmneut4mmu - gChi0musmuRR(4)*2.d0*
     $     mneut(4)*mmu*b0mneut4mmu
      
      neutralino(1, 2) = 
     $     fChi0musmuLR(1)*gmneut1mmu - gChi0musmuLR(1)*2.d0*
     $     mneut(1)*mmu*b0mneut1mmu +
     $     fChi0musmuLR(2)*gmneut2mmu - gChi0musmuLR(2)*2.d0*
     $     mneut(2)*mmu*b0mneut2mmu +
     $     fChi0musmuLR(3)*gmneut3mmu - gChi0musmuLR(3)*2.d0*
     $     mneut(3)*mmu*b0mneut3mmu +
     $     fChi0musmuLR(4)*gmneut4mmu - gChi0musmuLR(4)*2.d0*
     $     mneut(4)*mmu*b0mneut4mmu 
      
!-------------------------------------------------------------------
      
      MSLD3(1,1) = mSLRG(2,2) + mmu**2.d0 + geL*MZ*MZ*dcos(2.d0*beta)
      MSLD3(1,2) = mmu*((AERG(2,2) - (sgnmu*modmu)*dtan(beta)))
      MSLD3(2,1) = MSLD3(1,2)
      MSLD3(2,2) = mSERG(2,2) + mmu**2.d0 + geR*MZ*MZ*dcos(2.d0*beta)

!----------------------------------------------------------------------

      pisMu(1,1) = (cstop(1,1) + csbtm(1,1) + gtterm(1,1) +
     $     chargino(1,1) + neutralino(1,1) + higgsterm(1,1))/
     $     (16.d0*pi*pi)

      pisMu(1,2) = (cstop(1,2) + csbtm(1,2) + gtterm(1,2)  +
     $     chargino(1,2) + neutralino(1,2) + higgsterm(1,2))/
     $     (16.d0*pi*pi)

      pisMu(2,2) = (cstop(2,2) + csbtm(2,2) + gtterm(2,2) +
     $     chargino(2,2) + neutralino(2,2) +  higgsterm(2,2))/
     $     (16.d0*pi*pi)
      
      pisMu(2,1) = pisMu(1,2)


      pisMu(1,1) = MSLD3(1,1) - pisMu(1,1)
      pisMu(1,2) = MSLD3(1,2) - pisMu(1,2)
      pisMu(2,2) = MSLD3(2,2) - pisMu(2,2)
      pisMu(2,1) = pisMu(1,2) 

!----------------------------------------------------------------------
C     find the singular values and the diagonalising matrices
!----------------------------------------------------------------------

      info = 10
      AOK = 0
      
    !  call dsyev('V','U',2,pisMu,2,SMueg,work,lwork,info)
      
      Call CEigensystem(2,pisMu,2,SMueg,smumix,2,1)
      if(info.eq.0) then
         AOK = AOK + 1
      endif

      RETURN

      END SUBROUTINE pisMul

C==================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     
C     pisel is checked on 30/05/2010 @ 16:30.
C     
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      
      SUBROUTINE pisEl(p,q,g,gp,tanbeta,mSLRG,mSERG,yeRG,
     $     AERG,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu,
     $     modmu,ON,OCL,OCR,vev1,M3t,SEeg)

 
      IMPLICIT NONE 

      INTEGER i,j,info,AOK,lwork
      parameter(lwork=35)
      DOUBLE PRECISION work(lwork)
      DOUBLE PRECISION tanbeta, mSLRG(3,3),mSERG(3,3)
      DOUBLE PRECISION AERG(3,3),SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yeRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g,M3t

      double precision vev1,beta
 
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetac,thetas,thetamu

      DOUBLE PRECISION thetau,thetad,thetae
      DOUBLE PRECISION costhetae,sinthetae

      DOUBLE PRECISION pisE(2,2)
 
      DOUBLE PRECISION higgsterm(2,2)
      DOUBLE PRECISION dnu(4),dnd(4),cn(4)

      DOUBLE PRECISION lselLselLR(4,2),lselLsel12(4,2)
      DOUBLE PRECISION lselRselLR(4,2),lselRsel12(4,2)
      
      DOUBLE PRECISION aPsi0selr(4), bPsi0selr(4), aPsi0sell(4)
      DOUBLE PRECISION bPsi0sell(4) 
      DOUBLE PRECISION aChi0sell(4), bChi0sell(4), aChi0selr(4)
      DOUBLE PRECISION bChi0selr(4)

      DOUBLE PRECISION gChi0elselLL(4), fChi0elselLL(4)
      DOUBLE PRECISION gChi0elselLR(4), fChi0elselLR(4)
      DOUBLE PRECISION gChi0elselRR(4), fChi0elselRR(4)

      DOUBLE PRECISION bPsicsell(2), bPsicselr(2), aPsicsell(2)
      DOUBLE PRECISION aPsicselr(2)
      DOUBLE PRECISION aChicselr(2), aChicsell(2)
      DOUBLE PRECISION bChicsell(2),bChicselr(2)

      DOUBLE PRECISION gtterm(2,2), cstop(2,2), csbtm(2,2) 
      data gtterm/ 4 * 0.d0/, cstop/ 4 * 0.d0/, csbtm/ 4 * 0.d0/

      DOUBLE PRECISION  chargino(2,2), neutralino(2,2), mSLD3(2,2)
      data chargino/ 4 * 0.d0/,neutralino/ 4 * 0.d0/, mSLD3/ 4 * 0.d0/

      DOUBLE PRECISION fChselLL(2), gChselLL(2) 
      DOUBLE PRECISION fChselLR(2), gChselLR(2) 
      DOUBLE PRECISION fChselRR(2), gChselRR(2)
 
      DOUBLE PRECISION intm1(2), intm2(2)

      DOUBLE PRECISION lHcselrsnulr(2,2),lHcselrsnu12(2,2)
      DOUBLE PRECISION lHcsellsnu12(2,2)
      DOUBLE PRECISION lHcsellsnulr(2,2) 

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION lHselLsel12(4,2),lHselRsel12(4,2)
      data lHselLsel12/ 8 * 0.d0/, lHselRsel12/ 8 * 0.d0/

      
      DOUBLE PRECISION rthetae(2,2),ralpha(2,2)

      DOUBLE PRECISION ggme,fme10,fme20,a0b1,a0b2,b0M3me,a0t1,a0t2
      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz
      DOUBLE PRECISION b0h0mse(4,2),b0hcmstop(2,2)

      DOUBLE PRECISION a0md1,a0md2,a0ms1,a0ms2,a0mtau1,a0mtau2,a0snu2
      DOUBLE PRECISION a0me1,a0me2,a0mmu1,a0mmu2,a0snu1,a0snu3
      DOUBLE PRECISION a0mc1,a0mc2,a0mu1,a0mu2
      
      DOUBLE PRECISION fme1MZ,fme2MZ,f0MW,fsnu1MW
      DOUBLE PRECISION gmchargino1me,gmchargino2me
      DOUBLE PRECISION b0mchargino1me,b0mchargino2me,b0mneut4me
      DOUBLE PRECISION b0mneut1me,b0mneut2me,b0mneut3me
      DOUBLE PRECISION gmneut1me,gmneut2me,gmneut3me,gmneut4me
      
      DOUBLE PRECISION SEeg(2),semix(2,2)

      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2
      DOUBLE PRECISION me1,me2,nmneut(4)

      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae
  

      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d
  

      include 'stdinputs.h'

C------------------------------------
      
      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      MZ = MZpole

      SEeg(1) = 0.d0
      SEeg(2) = 0.d0

      
      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw = dsqrt(sinsqthw)
      thw = dasin(sinthw)
      costhw = dcos(thw)
      cos2thw = dcos(2.d0*thw)

 
      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))
      
      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))

       
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

!-------------------------------------------------------

      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)
      costhetae   = dcos(thetae)
      sinthetae   = dsin(thetae)


      mt1 = mtL
      mt2 = mtR

      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR

      me1 = meL
      me2 = meR

C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta     = datan(tanbeta) 
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha  = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      salpha = dsin(alpha)
      calpha = dcos(alpha)
      cosbeta  = dcos(datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)
!-----------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0
!----------------------------------------------------------
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

!---------------------------------------------------------------    

      lselLselLR(1,1) = (g*MZ*geL*cosbeta/costhw) + 
     $     (yeRG(1,1)**2.d0) * vev1
      lselLselLR(1,2) = (yeRG(1,1)*AERG(1,1))/dsqrt(2.d0)

      lselLselLR(2,1) = ((- g*MZ*geL*sinbeta)/costhw)
      lselLselLR(2,2) = (-yeRG(1,1)*sgnmu*modmu)/dsqrt(2.d0)

      lselLselLR(3,1) = 0.d0
      lselLselLR(3,2) = (-1.d0/dsqrt(2.d0)) *
     $     (-sgnmu*modmu*sinbeta*yeRG(1,1) + 
     $     yeRG(1,1)*AERG(1,1)*cosbeta)

      lselLselLR(4,1) = 0.d0
      lselLselLR(4,2) = (-1.d0/dsqrt(2.d0)) * 
     $     (-sgnmu*modmu*cosbeta*yeRG(1,1) - 
     $     yeRG(1,1)*AERG(1,1)*sinbeta)



      lselRselLR(1,1) = lselLselLR(1,2)
      lselRselLR(1,2) = (g*MZ*geR*cosbeta/costhw) + 
     $     yeRG(1,1)**2.d0 * vev1

      lselRselLR(2,1) = lselLselLR(2,2)
      lselRselLR(2,2) = -g*MZ*geR*sinbeta/costhw

      lselRselLR(3,1) = -lselLselLR(3,2)
      lselRselLR(3,2) = 0.d0

      lselRselLR(4,1) = -lselLselLR(4,2)
      lselRselLR(4,2) = 0.d0

!------------------------------------  
  
      loopmixst: DO i = 1, 4

      intm1(1) = lselLselLR(i,1)
      intm1(2) = lselLselLR(i,2)

      call rmat2d(thetae,rthetae)
            
      intm2(1) = rthetae(1,1)*intm1(1) + rthetae(1,2)*intm1(2)
      intm2(2) = rthetae(2,1)*intm1(1) + rthetae(2,2)*intm1(2)

      lselLsel12(i,1) = intm2(1)
      lselLsel12(i,2) = intm2(2)

      intm1(1) = lselRselLR(i,1)
      intm1(2) = lselRselLR(i,2)

      intm2(1) = rthetae(1,1)*intm1(1) + rthetae(1,2)*intm1(2)
      intm2(2) = rthetae(2,1)*intm1(1) + rthetae(2,2)*intm1(2)

      lselRsel12(i,1) = intm2(1)
      lselRsel12(i,2) = intm2(2)

      ENDDO loopmixst

!--------------------------------------------------------------------------
C       Mixing CP-even Higgs 
!--------------------------------------------------------------------------

      loopmixcpeven: DO i = 1, 2
      
      intm1(1) = lselLsel12(1,i)
      intm1(2) = lselLsel12(2,i)

      call rmat2d(alpha,ralpha)
            
      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lHselLsel12(1,i) = intm2(1)
      lHselLsel12(2,i) = intm2(2)

      intm1(1) = lselRsel12(1,i)
      intm1(2) = lselRsel12(2,i)

      intm2(1) = ralpha(2,2)*intm1(1) + ralpha(2,1)*intm1(2)
      intm2(2) = ralpha(1,2)*intm1(1) + ralpha(1,1)*intm1(2)

      lHselRsel12(1,i) = intm2(1)
      lHselRsel12(2,i) = intm2(2)
     
      ENDDO loopmixcpeven
!----------------------------------------------------------------------
C     Feynman rules - charged higgs
!<--------- (H+ G+, L R) basis


      lHcsellsnulr(1, 1) = (g*MW*sin2beta - 
     $     (yeRG(1,1)**2.d0) * vev1 * sinbeta)/dsqrt(2.d0)

      lHcsellsnulr(1,2) = 0.d0 

      lHcsellsnulr(2,1) = (-g*MW*cos2beta - 
     $     yeRG(1,1)**2.d0 * vev1 * cosbeta)/dsqrt(2.d0)

      lHcsellsnulr(2,2) = 0.d0

!-------------------------------------

      intm1(1) = lHcsellsnulr(1,1)
      intm1(2) = lHcsellsnulr(1,2)

      call rmat2d(thetae,rthetae)
            
      intm2(1) = rthetae(1,1)*intm1(1) + rthetae(1,2)*intm1(2)
      intm2(2) = rthetae(2,1)*intm1(1) + rthetae(2,2)*intm1(2)

 
      lHcsellsnu12(1, 1) = intm2(1)
      lHcsellsnu12(1, 2) = intm2(2)

      intm1(1) = lHcsellsnulr(2, 1)
      intm1(2) = lHcsellsnulr(2, 2)


                  
      intm2(1) = rthetae(1,1)*intm1(1) + rthetae(1,2)*intm1(2)
      intm2(2) = rthetae(2,1)*intm1(1) + rthetae(2,2)*intm1(2)

      lHcsellsnu12(2, 1) = intm2(1)
      lHcsellsnu12(2, 2) = intm2(2)

!-------------------------------------------
!<------------ (H+ G+, L R) basis

      lHcselrsnulr(1, 1) = yeRG(1,1)*(-sgnmu*modmu*cosbeta - 
     $     AERG(1,1)*sinbeta)

      lHcselrsnulr(1, 2) = 0.d0

      lHcselrsnulr(2, 2) = yeRG(1,1)*(-sgnmu*modmu*sinbeta + 
     $     AERG(1,1)*cosbeta)  

      lHcselrsnulr(2, 1) = 0.d0
!-------------------------------------------

      intm1(1) = lHcselrsnulr(1,1)
      intm1(2) = lHcselrsnulr(1,2)

      call rmat2d(thetae,rthetae)
            
      intm2(1) = rthetae(1,1)*intm1(1) + rthetae(1,2)*intm1(2)
      intm2(2) = rthetae(2,1)*intm1(1) + rthetae(2,2)*intm1(2)

      lHcselrsnu12(1,1) = intm2(1)
      lHcselrsnu12(1,2) = intm2(2)

      intm1(1) = lHcselrsnulr(2,1)
      intm1(2) = lHcselrsnulr(2,2)

                 
      intm2(1) = rthetae(1,1)*intm1(1) + rthetae(1,2)*intm1(2)
      intm2(2) = rthetae(2,1)*intm1(1) + rthetae(2,2)*intm1(2)
  
      lHcselrsnu12(2, 1) = intm2(1)
      lHcselrsnu12(2, 2) = intm2(2)

!------------------------------------------------------------------------
C     Feynman rules - neutralino
!------------------------------------------------------------------------


      aPsi0selr(1) = yeR*gp/(dsqrt(2.d0))
      aPsi0selr(2) = 0.d0
      aPsi0selr(3) = 0.d0
      aPsi0selr(4) = 0.d0

      bPsi0sell(1) = yeL*gp/(dsqrt(2.d0))
      bPsi0sell(2) = g/dsqrt(2.d0)
      bPsi0sell(3) = 0.d0
      bPsi0sell(4) = 0.d0

      aPsi0sell(1) = 0.d0
      aPsi0sell(2) = 0.d0
      aPsi0sell(3) = 0.d0 
      aPsi0sell(4) = yeRG(1,1)

      bPsi0selr(1) = 0.d0
      bPsi0selr(2) = 0.d0
      bPsi0selr(3) = 0.d0
      bPsi0selr(4) = yeRG(1,1)


      aChi0sell(1) = ON(1,1)*aPsi0sell(1) + ON(1,2)*aPsi0sell(2) +
     $     ON(1,3)*aPsi0sell(1) + ON(1,4)*aPsi0sell(4) 
      aChi0sell(2) = ON(2,1)*aPsi0sell(1) + ON(2,2)*aPsi0sell(2) +
     $     ON(2,3)*aPsi0sell(1) + ON(2,4)*aPsi0sell(4) 
      aChi0sell(3) = ON(3,1)*aPsi0sell(1) + ON(3,2)*aPsi0sell(2) +
     $     ON(3,3)*aPsi0sell(1) + ON(3,4)*aPsi0sell(4) 
      aChi0sell(4) = ON(4,1)*aPsi0sell(1) + ON(4,2)*aPsi0sell(2) +
     $     ON(4,3)*aPsi0sell(1) + ON(4,4)*aPsi0sell(4) 


      bChi0sell(1) = ON(1,1)*bPsi0sell(1) + ON(1,2)*bPsi0sell(2) +
     $     ON(1,3)*bPsi0sell(1) + ON(1,4)*bPsi0sell(4) 
      bChi0sell(2) = ON(2,1)*bPsi0sell(1) + ON(2,2)*bPsi0sell(2) +
     $     ON(2,3)*bPsi0sell(1) + ON(2,4)*bPsi0sell(4) 
      bChi0sell(3) = ON(3,1)*bPsi0sell(1) + ON(3,2)*bPsi0sell(2) +
     $     ON(3,3)*bPsi0sell(1) + ON(3,4)*bPsi0sell(4) 
      bChi0sell(4) = ON(4,1)*bPsi0sell(1) + ON(4,2)*bPsi0sell(2) +
     $     ON(4,3)*bPsi0sell(1) + ON(4,4)*bPsi0sell(4) 




      aChi0selr(1) = ON(1,1)*aPsi0selr(1) + ON(1,2)*aPsi0selr(2) +
     $     ON(1,3)*aPsi0selr(3) + ON(1,4)*aPsi0selr(4)
      aChi0selr(2) = ON(2,1)*aPsi0selr(1) + ON(2,2)*aPsi0selr(2) +
     $     ON(2,3)*aPsi0selr(3) + ON(2,4)*aPsi0selr(4)
      aChi0selr(3) = ON(3,1)*aPsi0selr(1) + ON(3,2)*aPsi0selr(2) +
     $     ON(3,3)*aPsi0selr(3) + ON(3,4)*aPsi0selr(4)
      aChi0selr(4) = ON(4,1)*aPsi0selr(1) + ON(4,2)*aPsi0selr(2) +
     $     ON(4,3)*aPsi0selr(3) + ON(4,4)*aPsi0selr(4)


      bChi0selr(1) = ON(1,1)*bPsi0selr(1) + ON(1,2)*bPsi0selr(2) +
     $     ON(1,3)*bPsi0selr(3) + ON(1,4)*bPsi0selr(4) 
      bChi0selr(2) = ON(2,1)*bPsi0selr(1) + ON(2,2)*bPsi0selr(2) +
     $     ON(2,3)*bPsi0selr(3) + ON(2,4)*bPsi0selr(4) 
      bChi0selr(3) = ON(3,1)*bPsi0selr(1) + ON(3,2)*bPsi0selr(2) +
     $     ON(3,3)*bPsi0selr(3) + ON(3,4)*bPsi0selr(4) 
      bChi0selr(4) = ON(4,1)*bPsi0selr(1) + ON(4,2)*bPsi0selr(2) +
     $     ON(4,3)*bPsi0selr(3) + ON(4,4)*bPsi0selr(4) 


      loopchtllrr: DO i = 1, 4

      fChi0elselLL(i) = (aChi0sell(i)*aChi0sell(i) + 
     $     bChi0sell(i)*bChi0sell(i))
      gChi0elselLL(i) = (bChi0sell(i)*aChi0sell(i) + 
     $     bChi0sell(i)*aChi0sell(i))
      fChi0elselRR(i) = (aChi0selr(i)*aChi0selr(i) + 
     $     bChi0selr(i)*bChi0selr(i))
      gChi0elselRR(i) = (bChi0selr(i)*aChi0selr(i) + 
     $     bChi0selr(i)*aChi0selr(i))
      fChi0elselLR(i) = (aChi0selr(i)*aChi0sell(i) + 
     $     bChi0selr(i)*bChi0sell(i))
      gChi0elselLR(i) = (bChi0sell(i)*aChi0selr(i) + 
     $     bChi0selr(i)*aChi0sell(i))
      
      ENDDO loopchtllrr

!---------------------------------------------------------------------------
c     Feynman Rules for chargino
!---------------------------------------------------------------------------
      
      aPsicsell(1) = 0.d0
      aPsicsell(2) = 0.d0

      aPsicselr(1) = 0.d0
      aPsicselr(2) = 0.d0

      bPsicsell(1) = g
      bPsicsell(2) = 0.d0

      bPsicselr(1) = 0.d0
      bPsicselr(2) = - yeRG(1,1)
      

      aChicsell(1) = OCR(1,1)*aPsicsell(1) + OCR(1,2)*aPsicsell(2)
      aChicsell(2) = OCR(2,1)*aPsicsell(1) + OCR(2,2)*aPsicsell(2)

      bChicsell(1) = OCL(1,1)*bPsicsell(1) + OCL(1,2)*bPsicsell(2)
      bChicsell(2) = OCL(2,1)*bPsicsell(1) + OCL(2,2)*bPsicsell(2)

      aChicselr(1) = OCR(1,1)*aPsicselr(1) + OCR(1,2)*aPsicselr(2)
      aChicselr(2) = OCR(2,1)*aPsicselr(1) + OCR(2,2)*aPsicselr(2)

      bChicselr(1) = OCL(1,1)*bPsicselr(1) + OCL(1,2)*bPsicselr(2)
      bChicselr(2) = OCL(2,1)*bPsicselr(1) + OCL(2,2)*bPsicselr(2)

      
      loopchst: DO i = 1, 2

      fChselLL(i) = (aChicsell(i)*aChicsell(i) +
     $     bChicsell(i)*bChicsell(i))
      gChselLL(i) = (bChicsell(i)*aChicsell(i) +
     $     aChicsell(i)*bChicsell(i))
      fChselLR(i) = (aChicsell(i)*aChicselr(i) +
     $     bChicsell(i)*bChicselr(i))
      gChselLR(i) = (bChicsell(i)*aChicselr(i) +
     $     aChicsell(i)*bChicselr(i))
      fChselRR(i) = (aChicselr(i)*aChicselr(i) +
     $     bChicselr(i)*bChicselr(i))
      gChselRR(i) = (bChicselr(i)*aChicselr(i) +
     $     aChicselr(i)*bChicselr(i))

      ENDDO loopchst

!--------------------------------------------------------------------
C     Corrections Begin
!------------------------------------------------------------------

      call funcg(p,M3t,me,q,ggme)
     
      call f(p,me1,0.d0,q,fme10)
      call f(p,me2,0.d0,q,fme20)

      call a0(muR,q,a0mu2)
      call a0(muL,q,a0mu1)
      call a0(mcR,q,a0mc2)
      call a0(mcL,q,a0mc1)
      call a0(mt2,q,a0t2)
      call a0(mt1,q,a0t1)

      call a0(msR,q,a0ms2)
      call a0(msL,q,a0ms1)
      call a0(mdR,q,a0md2)
      call a0(mdL,q,a0md1)
      call a0(mb2,q,a0b2)
      call a0(mb1,q,a0b1)

      call a0(me2,q,a0me2)
      call a0(me1,q,a0me1)
      call a0(mmuR,q,a0mmu2)
      call a0(mmuL,q,a0mmu1)
      call a0(mtau2,q,a0mtau2)
      call a0(mtau1,q,a0mtau1)

      call a0(snu(1),q,a0snu1)
      call a0(snu(2),q,a0snu2)
      call a0(snu(3),q,a0snu3)


      call a0(mA0,q,a0mA)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)
      call a0(mHpm,q,a0mHpm)
      call a0(MW,q,a0Mw)
      call a0(MZ,q,a0Mz)

      call f(p,0.d0,MW,q,f0Mw)
      call f(p,0.d0,MW,q,f0Mw)
      call f(p,me1,MZ,q,fme1Mz)
      call f(p,me2,MZ,q,fme2Mz)

      call f(p,snu(1),MW,q,fsnu1MW)

      call b0(p,M3t,me,q,b0M3me)
      
!----------------------------------------------------------------
      gtterm(1, 1) = 0.d0

      gtterm(2, 2) = 0.d0

      gtterm(1, 2) = 0.d0
!----------------------------------------------------------------

      csbtm(1, 1) = yeRG(1,1)**2.d0 * ((sinthetae)**2.d0*a0me1 + 
     $     (costhetae)**2.d0*a0me2)

      csbtm(2, 2) = yeRG(1,1)**2.d0 * ((costhetae)**2.d0*a0me1 + 
     $     (sinthetae)**2.d0*a0me2)

      csbtm(1, 2) = yeRG(1,1)**2.d0 * costhetae*sinthetae * 1.d0 *
     $     (a0me1 - a0me2)
!------------------------------------------------------------------

      cstop(1, 1) = 0.d0

      cstop(2, 2) =  yeRG(1,1)**2.d0 * a0snu1
!--------------------------------------------------------------------

      higgsterm(1,1) = 0.5d0 * ((yeRG(1,1)**2.d0) * dnd(1) -
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0 * ((yeRG(1,1)**2.d0) * dnd(2) - 
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0 * ((yeRG(1,1)**2.d0) * dnd(3) - 
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(3)) * a0MZ +
     $     0.5d0 * ((yeRG(1,1)**2.d0) * dnd(4) - 
     $     ((g*g)*geL*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA
      
      
      higgsterm(1,2) = 0.d0

      higgsterm(2,2) = 0.5d0*(yeRG(1,1)**2.d0 * dnd(1) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(1)) * a0mHu +
     $     0.5d0*(yeRG(1,1)**2.d0 * dnd(2) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(2)) * a0mh +
     $     0.5d0*(yeRG(1,1)**2.d0 * dnd(3) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(3)) * a0Mz +
     $     0.5d0*(yeRG(1,1)**2.d0 * dnd(4) - 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(4)) * a0mA


      higgsterm(1,1) = higgsterm(1,1) + 
     $     (0.d0 + ((g*g) * 
     $     ((geL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(3))) * a0Mhpm +
     $     (0.d0 + ((g*g) * 
     $     ((geL*0.5d0/(costhw**2.d0)) + 0.5d0) * cn(4))) * a0MW

      higgsterm(2,2) = higgsterm(2,2) + 
     $     ((yeRG(1,1)**2.d0) * dnu(3) +  
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(3)) * a0Mhpm +
     $     ((yeRG(1,1)**2.d0) * dnu(4) + 
     $     ((g*g)*geR*0.5d0/(costhw**2.d0)) * cn(4)) * a0MW !<----------------chkd

!--------------------------
      
      call b0(p,mHu0,me1,q,b0h0mse(1,1))
      call b0(p,mHu0,me2,q,b0h0mse(1,2))
      call b0(p,mh0,me1,q,b0h0mse(2,1))
      call b0(p,mh0,me2,q,b0h0mse(2,2))
      call b0(p,MZ,me1,q,b0h0mse(3,1))
      call b0(p,MZ,me2,q,b0h0mse(3,2))
      call b0(p,mA0,me1,q,b0h0mse(4,1))
      call b0(p,mA0,me2,q,b0h0mse(4,2))
      
      
      loophi: DO i = 1, 4       !<-------------neutral higgs terms
      loophj: DO j = 1, 2
      
      higgsterm(1,1) = higgsterm(1,1) + 
     $     (lHselLsel12(i,j)**2.d0) * b0h0mse(i,j)
      
      higgsterm(1,2) = higgsterm(1,2) + 
     $     lHselLsel12(i,j) * lHselRsel12(i,j) * b0h0mse(i,j)
      
      higgsterm(2,2) = higgsterm(2,2) + 
     $     (lHselRsel12(i,j)**2.d0) * b0h0mse(i,j)
      
      ENDDO loophj
      ENDDO loophi
!-------------------------------------
      
      call b0(p,0.d0,mHpm,q,b0hcmstop(1,1))
      call b0(p,0.d0,mHpm,q,b0hcmstop(1,2))
      call b0(p,0.d0,MW,q,b0hcmstop(2,1))
      call b0(p,0.d0,MW,q,b0hcmstop(2,2))
      
      
      loophci: DO i = 1, 2      !<-------------charged higgs terms
      loophcj: DO j = 1, 2
      
      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (lHcsellsnu12(i,j)**2.d0) * b0hcmstop(i,j)       
      
      higgsterm(1, 2) = higgsterm(1,2) + 
     $     lHcsellsnu12(i, j)*lHcselrsnu12(i, j)*b0hcmstop(i,j)
      
      higgsterm(2, 2) = higgsterm(2, 2) + 
     $     (lHcselrsnu12(i,j)**2.d0) * b0hcmstop(i,j)
      
      ENDDO loophcj
      ENDDO loophci

!--------------------------------------

      higgsterm(1, 1) = higgsterm(1,1) + 
     $     (4.d0*(g*g)/costhw**2.d0)*(geL*geL)*a0MZ  + 
     $     2.d0*(g*g)*a0MW + 
     $     ((1.d0 *g*g*sinsqthw) * ((costhetae**2.d0) * 
     $     fme10 + (sinthetae**2.d0) * fme20)) +
     $     (g*geL/costhw)**2.d0 * ((costhetae**2.d0) * fme1Mz +
     $     (sinthetae**2.d0) * fme2mz) +
     $     (g*g)*0.5d0 * fsnu1MW  +
     $     (g*g) * 0.25d0 * (1.d0 * ((costhetae**2.d0) * a0me1 + 
     $     (sinthetae**2.d0) * a0me2) + 2.d0 * a0snu1 ) -
     $     (g*g) * 0.5d0 * (- 1.5d0 * a0md1 - 1.5d0*a0ms1 +
     $     1.5d0 * a0mu1 + 1.5d0*a0mc1 +
     $     0.5d0 * (a0snu1 + a0snu2 + a0snu3) -
     $     0.5d0 * (a0me1 + a0mmu1 + ((costhetatau**2.d0) * a0mtau1 +
     $     (sinthetatau**2.d0) * a0mtau2))) +
     $     (gp**2.d0) * 0.25d0 *(yeL**2.d0) * 
     $     ((costhetae**2.d0) * a0me1 + 
     $     (sinthetae**2.d0) * a0me2) +
     $     (gp**2.d0) * 0.25d0 * yeL * ((3.d0*ydL*(a0md1 + a0ms1 + 
     $     ((costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2))) +
     $     (3.d0 * ydR *(a0md2 + a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2)) +
     $     (3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t1)) +
     $     (3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (costhetat**2.d0)*a0t2 + (sinthetat**2.d0)*a0t1)) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))
     
!--------------------------------------------
 
      higgsterm(2, 2) = higgsterm(2,2) +
     $     (4.d0*(g*g)/costhw**2.d0) * (geR**2.d0) * a0MZ + 
     $     ((1.d0)*g*g*sinsqthw) * ((sinthetae)*fme10 + 
     $     (costhetae**2.d0)*fme20) +
     $     (g*geR/costhw)**2.d0 * ((sinthetae**2.d0) * fme1MZ + 
     $     (costhetae**2.d0) * fme2MZ) +
     $     (gp**2.d0) * 0.25d0 * (yeR**2.d0) * 
     $     ((sinthetae**2.d0) * a0me1 + 
     $     (costhetae**2.d0) * a0me2) +
     $     (gp**2.d0) * 0.25d0 * yeR *  
     $     (3.d0 * ydL * (a0md1 + a0ms1 + 
     $     (costhetab**2.d0)*a0b1 + (sinthetab**2.d0)*a0b2) +
     $     3.d0 * ydR * (a0md2 +  a0ms2 + 
     $     (sinthetab**2.d0)*a0b1 + (costhetab**2.d0)*a0b2) +
     $     3.d0 * yuL * (a0mu1 + a0mc1 + 
     $     (costhetat**2.d0)*a0t1 + (sinthetat**2.d0)*a0t2) +
     $     3.d0 * yuR * (a0mu2 + a0mc2 + 
     $     (sinthetat**2.d0)*a0t1 + (costhetat**2.d0)*a0t2) +
     $     yeL * (a0me1 + a0mmu1 + 
     $     (sinthetatau**2.d0)*a0mtau2 + (costhetatau**2.d0)*a0mtau1) +
     $     yeR * (a0me2 + a0mmu2 + 
     $     (sinthetatau**2.d0)*a0mtau1 + (costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0snu1 + a0snu2 + a0snu3))
      
!--------------------
      
      higgsterm(1,2) = higgsterm(1,2) + 
     $     (gp**2.d0) * yeL*yeR * sinthetae*costhetae *
     $     (a0me1 - a0me2) +
     $     (1.d0 * g*g * sinsqthw) * sinthetae*costhetae * 
     $     (fme10 - fme20) -
     $     ((g*g)/costhw**2.d0) * geL*geR * sinthetae*costhetae *
     $     (fme1Mz - fme2Mz)

!----------------------------------------------------------------------
!     Chargino term
!----------------------------------------------------------------------
      
      call funcg(p,mchargino(1),0.d0,q,gmchargino1me)
      call funcg(p,mchargino(2),0.d0,q,gmchargino2me)
      call b0(p,mchargino(1),0.d0,q,b0mchargino1me)
      call b0(p,mchargino(2),0.d0,q,b0mchargino2me)
      
      chargino(1, 1) =  fChselLL(1)*gmchargino1me -
     $     gChselLL(1)*mchargino(1)*0.d0*b0mchargino1me*2.d0 +
     $     fChselLL(2)*gmchargino2me -
     $     gChselLL(2)*mchargino(2)*0.d0*b0mchargino2me*2.d0
      
      chargino(1, 2) =  fChselLR(1)*gmchargino1me -
     $     gChselLR(1)*mchargino(1)*0.d0*b0mchargino1me*2.d0 +
     $     fChselLR(2)*gmchargino2me -
     $     gChselLR(2)*mchargino(2)*0.d0*b0mchargino2me*2.d0
      
      chargino(2, 2) =  fChselRR(1)*gmchargino1me - 
     $     gChselRR(1)*mchargino(1)*0.d0*b0mchargino1me*2.d0 +
     $     fChselRR(2)*gmchargino2me - 
     $     gChselRR(2)*mchargino(2)*0.d0*b0mchargino2me*2.d0
      
!---------------------------------------------------------------------
!     neutralino terms
!---------------------------------------------------------------------
      
      call b0(p,nmneut(1),me,q,b0mneut1me)
      call b0(p,nmneut(2),me,q,b0mneut2me)
      call b0(p,nmneut(3),me,q,b0mneut3me)
      call b0(p,nmneut(4),me,q,b0mneut4me)
      
      call funcg(p,nmneut(1),me,q,gmneut1me)
      call funcg(p,nmneut(2),me,q,gmneut2me)
      call funcg(p,nmneut(3),me,q,gmneut3me)
      call funcg(p,nmneut(4),me,q,gmneut4me)
      
      neutralino(1, 1) = 
     $     fChi0elselLL(1)*gmneut1me - gChi0elselLL(1)*2.d0*
     $     mneut(1)*me*b0mneut1me +
     $     fChi0elselLL(2)*gmneut2me - gChi0elselLL(2)*2.d0*
     $     mneut(2)*me*b0mneut2me +
     $     fChi0elselLL(3)*gmneut3me - gChi0elselLL(3)*2.d0*
     $     mneut(3)*me*b0mneut3me +
     $     fChi0elselLL(4)*gmneut4me - gChi0elselLL(4)*2.d0*
     $     mneut(4)*me*b0mneut4me
      
      
      neutralino(2, 2) = 
     $     fChi0elselRR(i)*gmneut1me - gChi0elselRR(1)*2.d0*
     $     mneut(1)*me*b0mneut1me +
     $     fChi0elselRR(2)*gmneut2me - gChi0elselRR(2)*2.d0*
     $     mneut(2)*me*b0mneut2me +
     $     fChi0elselRR(3)*gmneut3me - gChi0elselRR(3)*2.d0*
     $     mneut(3)*me*b0mneut3me +
     $     fChi0elselRR(4)*gmneut4me - gChi0elselRR(4)*2.d0*
     $     mneut(4)*me*b0mneut4me
      
      neutralino(1, 2) = 
     $     fChi0elselLR(1)*gmneut1me - gChi0elselLR(1)*2.d0*
     $     mneut(1)*me*b0mneut1me +
     $     fChi0elselLR(2)*gmneut2me - gChi0elselLR(2)*2.d0*
     $     mneut(2)*me*b0mneut2me +
     $     fChi0elselLR(3)*gmneut3me - gChi0elselLR(3)*2.d0*
     $     mneut(3)*me*b0mneut3me +
     $     fChi0elselLR(4)*gmneut4me - gChi0elselLR(4)*2.d0*
     $     mneut(4)*me*b0mneut4me 
      
!-------------------------------------------------------------------

      MSLD3(1,1) = mSLRG(1,1) + me**2.d0 + geL*MZ*MZ*dcos(2.d0*beta)
      MSLD3(1,2) = me * ((AERG(1,1) - sgnmu*modmu*dtan(beta)))
      MSLD3(2,1) = MSLD3(1,2)
      MSLD3(2,2) = mSERG(1,1) + me**2.d0 + geR*MZ*MZ*dcos(2.d0*beta)

!----------------------------------------------------------------------

      pisE(1,1) = (cstop(1,1) + csbtm(1,1) + gtterm(1,1) +
     $     chargino(1,1) + neutralino(1,1) + higgsterm(1,1))/
     $     (16.d0*pi*pi)

      pisE(1,2) = (cstop(1,2) + csbtm(1,2) + gtterm(1,2)  +
     $     chargino(1,2) + neutralino(1,2) + higgsterm(1,2))/
     $     (16.d0*pi*pi)

      pisE(2,2) = (cstop(2,2) + csbtm(2,2) + gtterm(2,2) +
     $     chargino(2,2) + neutralino(2,2) +  higgsterm(2,2))/
     $     (16.d0*pi*pi)
      
      pisE(2,1) = pisE(1,2)


      pisE(1,1) = MSLD3(1,1) - pisE(1,1)
      pisE(1,2) = MSLD3(1,2) - pisE(1,2)
      pisE(2,2) = MSLD3(2,2) - pisE(2,2)
      pisE(2,1) = pisE(1,2) 

!----------------------------------------------------------------------
C     find the singular values and the diagonalising matrices
!----------------------------------------------------------------------
      
      info  = 10
      AOK = 0
      
    !  call dsyev('V','U',2,pisE,2,SEeg,work,lwork,info)
      
      if(info.eq.0) then
         AOK = AOK + 1
      endif
            
      Call CEigensystem(2,pisE,2,SEeg,semix,2,1)
      RETURN

      END SUBROUTINE pisEl

C==================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     
C     pisnutau is checked on 01/06/2010 @ 13:30.
C     
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      
      SUBROUTINE pitausnu(p,q,g,gp,mtau,tanbeta,yeRG,AERG,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,
     $     mhpmsq,mA0sq,sgnmu,modmu,ON,OCL,OCR,vev1,tsnu)

 
      IMPLICIT NONE 

      INTEGER i,j
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION AERG(3,3),SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yeRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      double precision mTau
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g

      double precision vev1,beta
 
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION pitsnu,tsnu

      DOUBLE PRECISION dnu(4),dnd(4),cn(4)

      
      DOUBLE PRECISION  charginoterm, neutralinoterm
      DOUBLE PRECISION stops,sbotterm,higgsterm
      data stops/ 1 * 0.d0/,sbotterm/ 1 * 0.d0/,higgsterm/ 1 * 0.d0/
      data charginoterm/ 1 * 0.d0/, neutralinoterm/1 * 0.d0/
      
      DOUBLE PRECISION intm1(2), intm2(2)

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      
      DOUBLE PRECISION rthetatau(2,2),ralpha(2,2)

      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz

      DOUBLE PRECISION a0mdL,a0mdR,a0msL,a0msR,a0msnu2
      DOUBLE PRECISION a0meL,a0meR,a0mmuL,a0mmuR,a0msnu1,a0msnu3
      DOUBLE PRECISION a0mcR,a0mcL,a0muL,a0muR,a0mbL,a0mbR,a0mtL,a0mtR
      DOUBLE PRECISION a0mtau1,a0mtau2

      DOUBLE PRECISION fmtau1MW,fmtau2MW,fsnu1MZ
      DOUBLE PRECISION gmchargino10,gmchargino20
      DOUBLE PRECISION gmneut10,gmneut20,gmneut30,gmneut40

      DOUBLE PRECISION aPsi0snul(4), bPsi0snul(4) 
      data aPsi0snul/ 4 * 0.d0/,bPsi0snul/ 4 * 0.d0/

      DOUBLE PRECISION aChi0snul(4), bChi0snul(4)
      data aChi0snul/ 4 * 0.d0/,bChi0snul/ 4 * 0.d0/
      
      DOUBLE PRECISION gChi0nusnuLL(4), fChi0nusnuLL(4)
      data gChi0nusnuLL/ 4 * 0.d0/,fChi0nusnuLL/ 4 * 0.d0/


      DOUBLE PRECISION bPsicsnul(2), aPsicsnul(2)
      data bPsicsnul/ 2 * 0.d0/,aPsicsnul/ 2 * 0.d0/

      DOUBLE PRECISION aPsicCStaul(2),aChicsnur(2), aChicsnul(2)
      data aPsicCStaul/ 2 * 0.d0/,aChicsnur/ 2 * 0.d0/,
     $     aChicsnul/ 2 * 0.d0/

      DOUBLE PRECISION bChicsnul(2),bChicsnur(2)
      data bChicsnul/ 2 * 0.d0/,bChicsnur/ 2 * 0.d0/

      DOUBLE PRECISION aChsnu(2, 2), bChsnu(2, 2)
      data aChsnu/ 4 * 0.d0/,bChsnu/ 4 * 0.d0/

      DOUBLE PRECISION fChsnuLL(2), gChsnuLL(2)
      data fChsnuLL/ 2 * 0.d0/,gChsnuLL/ 2 * 0.d0/
     
      DOUBLE PRECISION  lsnulsnulr(4),lhsnulsnu12(4)
      DOUBLE PRECISION lHcsnulstaulr(2, 2),lHcsnulstau12(2, 2)
      DOUBLE PRECISION b0mHusnu1,b0mhsnu1,b0mzsnu1,b0mAsnu1
      DOUBLE PRECISION b0mHpmtau1,b0mHpmtau2,b0mwmtau1,b0mwmtau2
      DOUBLE PRECISION b0mch1mtau,b0mch2mtau,nmneut(4)
      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2

      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

      common/sinsq_susy/sinsqthw_susy


      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d
  

      include 'stdinputs.h'

C----------------------------------------------------

      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      MZ = MZpole

      pitsnu = 0.d0
      tsnu = 0.d0     
      
      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw   = dsqrt(sinsqthw)
      thw      = dasin(sinthw)
      costhw   = dcos(thw)
      cos2thw  = dcos(2.d0*thw)
     

      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))
      
      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))
      
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

!-------------------------------------

      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)

      mt1 = mtL
      mt2 = mtR

      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR


C---------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta     = datan(tanbeta)
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      salpha   = dsin(alpha)
      calpha   = dcos(alpha)
      cosbeta  = dcos(datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)
!-----------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0
!-------------------------------------------------------------
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

!-----------------------------------------------------------

      lsnulsnulr(1) =  g * mz * gnuL * cosbeta/costhw
      lsnulsnulr(2) = -g * mz * gnuL * sinbeta/costhw
      lsnulsnulr(3) = 0.d0
      lsnulsnulr(4) = 0.d0

C     Mix CP-even Higgses up

      intm1(1) = lsnulsnulr(1)
      intm1(2) = lsnulsnulr(2)
      
      call rmat2d(alpha,ralpha)

      intm2(1) = ralpha(1,1)*intm1(1) + ralpha(1,2)*intm1(2)
      intm2(2) = ralpha(2,1)*intm1(1) + ralpha(2,2)*intm1(2)
      
      lhsnulsnu12(1) = intm2(1)
      lhsnulsnu12(2) = intm2(2)
      lhsnulsnu12(3) = 0.d0
      lhsnulsnu12(4) = 0.d0

!-------------------------------------------------------------------
C     Charged Higgs Feynman rules
!     (H+ G+, L R) basis

      lHcsnulstaulr(1, 1) = ((g*MW*dsin(2.d0*beta)
     $     - yeRG(3,3)**2.d0 * vev1 *sinbeta)/dsqrt(2.d0))

      lHcsnulstaulr(1, 2) = (- sgnmu*modmu*yeRG(3,3)*cosbeta - 
     $     yeRG(3,3) * AERG(3,3) * sinbeta)

      lHcsnulstaulr(2, 1) = (-g*MW*dcos(2.d0*beta) 
     $     + yeRG(3,3)**2.d0*vev1*cosbeta)/dsqrt(2.d0)

      lHcsnulstaulr(2, 2) = (-yeRG(3,3)*sgnmu*modmu*sinbeta + 
     $     yeRG(3,3) * AERG(3,3) * cosbeta) 


      intm1(1) = lHcsnulstaulr(1, 1)
      intm1(2) = lHcsnulstaulr(1, 2)

      call rmat2d(thetatau,rthetatau)

      intm2(1) = rthetatau(1,1)*intm1(1) + rthetatau(1,2)*intm1(2)
      intm2(2) = rthetatau(2,1)*intm1(1) + rthetatau(2,2)*intm1(2)
      
      lHcsnulstau12(1, 1) = intm2(1)
      lHcsnulstau12(1, 2) = intm2(2)

      intm1(1) = lHcsnulstaulr(2, 1)
      intm1(2) = lHcsnulstaulr(2, 2)

      intm2(1) = rthetatau(1,1)*intm1(1) + rthetatau(1,2)*intm1(2)
      intm2(2) = rthetatau(2,1)*intm1(1) + rthetatau(2,2)*intm1(2)

      lHcsnulstau12(2, 1) = intm2(1)
      lHcsnulstau12(2, 2) = intm2(2)

!-------------------------------------------------------------------------------
C     Neutralino Feynman rule
!-------------------------------------------------------------------------------

      bPsi0snul(1) = gp*ynuL/dsqrt(2.d0)
      bPsi0snul(2) = g / dsqrt(2.d0)
      bPsi0snul(3) = 0.d0
      bPsi0snul(4) = 0.d0

      aPsi0snul(1) = 0.d0
      aPsi0snul(2) = 0.d0
      aPsi0snul(3) = 0.d0
      aPsi0snul(4) = yeRG(3,3) 

!------------

      aChi0snul(1) = 0.d0
      aChi0snul(2) = 0.d0
      aChi0snul(3) = 0.d0
      aChi0snul(4) = 0.d0

      bChi0snul(1) = 0.d0
      bChi0snul(2) = 0.d0
      bChi0snul(3) = 0.d0
      bChi0snul(4) = 0.d0


      loopi: DO i = 1, 4
      loopj: DO j = 1, 4

      aChi0snul(i) = aChi0snul(i) + ON(i,j) * aPsi0snul(j)
      bChi0snul(i) = bChi0snul(i) + ON(i,j) * bPsi0snul(j)

      ENDDO loopj
      ENDDO loopi
!-------------


      lcosnui: DO i = 1, 4

      fChi0nusnuLL(i) = (aChi0snul(i) * aChi0snul(i) + 
     $     bChi0snul(i) * bChi0snul(i))
      gChi0nusnuLL(i) = (bChi0snul(i) * aChi0snul(i) + 
     $     bChi0snul(i) * aChi0snul(i))


      ENDDO lcosnui
      
!----------------------------------------------------------------------------
C     Chargino Feynman Rules


      aPsicsnul(1) = g
      aPsicsnul(2) = 0.d0

      bPsicsnul(1) =  0.d0
      bPsicsnul(2) = - yeRG(3,3)
      
!--------------
      aChicsnul(1) = 0.d0
      aChicsnul(2) = 0.d0

      bChicsnul(1) = 0.d0
      bChicsnul(2) = 0.d0

      loopchi: DO i = 1, 2
      loopchj: DO j = 1, 2

      aChicsnul(i) = aChicsnul(i) + OCR(i,j)*aPsicsnul(j)
      bChicsnul(i) = bChicsnul(i) + OCL(i,j)*bPsicsnul(j)

      ENDDO loopchj
      ENDDO loopchi
!---------------------

      lcosnuchi: DO i = 1, 2
      
      fChsnuLL(i) = (aChicsnul(i) * aChicsnul(i) +
     $     bChicsnul(i) * bChicsnul(i))
      gChsnuLL(i) = (bChicsnul(i)* aChicsnul(i) +
     $     aChicsnul(i)* bChicsnul(i))
      
      ENDDO lcosnuchi
!-----------------------------------------------------------

      call a0(mA0,q,a0mA)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)
      call a0(mHpm,q,a0mHpm)
      call a0(MW,q,a0Mw)
      call a0(MZ,q,a0Mz)

      call a0(muR,q,a0muR)
      call a0(muL,q,a0muL)
      call a0(mcR,q,a0mcR)
      call a0(mcL,q,a0mcL)
      call a0(mtL,q,a0mtL)
      call a0(mtR,q,a0mtR)

      call a0(msR,q,a0msR)
      call a0(msL,q,a0msL)
      call a0(mdR,q,a0mdR)
      call a0(mdL,q,a0mdL)
      call a0(mbL,q,a0mbL)
      call a0(mbR,q,a0mbR)

      call a0(meR,q,a0meR)
      call a0(meL,q,a0meL)
      call a0(mmuR,q,a0mmuR)
      call a0(mmuL,q,a0mmuL)
      call a0(mtau1,q,a0mtau1)
      call a0(mtau2,q,a0mtau2)

      call a0(snu(1),q,a0msnu1)
      call a0(snu(2),q,a0msnu2)
      call a0(snu(3),q,a0msnu3)

      call f(p,mtau1,MW,q,fmtau1Mw)
      call f(p,mtau2,MW,q,fmtau2Mw)
      call f(p,snu(1),Mz,q,fsnu1mz)


      call b0(p,mHu0,snu(1),q,b0mHusnu1)
      call b0(p,mh0,snu(1),q,b0mhsnu1)
      call b0(p,MZ,snu(1),q,b0mzsnu1)
      call b0(p,mA0,snu(1),q,b0mAsnu1)
      call b0(p,mtau1,mHpm,q,b0mHpmtau1)
      call b0(p,mtau2,mHpm,q,b0mHpmtau2)
      call b0(p,mtau1,MW,q,b0mwmtau1)
      call b0(p,mtau2,MW,q,b0mwmtau2)

!-----------------------------------------------------------------------------------
! Corrections begin

      sbotterm =  yeRG(3,3)**2.d0 * ((sinthetatau)**2.d0*a0mtau1 +
     $     (costhetatau)**2.d0*a0mtau2)


      higgsterm = 0.d0 +  
     $   0.5d0*((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0) * cn(1)) * a0mHu+
     $     0.5*((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0) * cn(2)) * a0mh +
     $     0.5*((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0) * cn(3)) * a0mz +
     $     0.5*((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0) * cn(4)) * a0mA 

 
      higgsterm = higgsterm + ((yeRG(3,3)**2.d0) * dnu(3) + (g)**2.d0 *
     $     ((gnuL*0.5d0/costhw**2.d0) - 0.5d0) * cn(3)) * a0mHpm +
     $     ((yeRG(3,3)**2.d0) * dnu(4) + (g)**2.d0 *
     $     ((gnuL*0.5d0/costhw**2.d0 - 0.5d0)) * cn(4)) * a0mw


      higgsterm = higgsterm + (lhsnulsnu12(1))**2.d0 * b0mHusnu1 +
     $     (lhsnulsnu12(2))**2.d0*b0mhsnu1 +
     $     (lhsnulsnu12(3))**2.d0*b0mzsnu1 +
     $     (lhsnulsnu12(4))**2.d0*b0mAsnu1
      



      higgsterm = higgsterm + (lHcsnulstau12(1,1))**2.d0*b0mHpmtau1 +
     $     (lHcsnulstau12(1,2))**2.d0*b0mHpmtau2 +
     $     (lHcsnulstau12(2,1))**2.d0*b0mwmtau1 +
     $     (lHcsnulstau12(2,2))**2.d0*b0mwmtau2
      

      higgsterm = higgsterm + 4.d0 * ((g*gnuL)**2.d0/(costhw)**2.d0) * 
     $     a0mz + 
     $     2.d0*(g)**2.d0*a0mw + 
     $     (g*gnuL/costhw)**2.d0 * fsnu1mz +
     $     (g)**2.d0*0.5d0*((costhetatau)**2.d0*fmtau1MW + 
     $     (sinthetatau)**2.d0*fmtau2MW) +
     $     (g)**2.d0*0.25*(a0msnu1 + 2.d0 *
     $     ((costhetatau)**2.d0 * a0mtau1 + 
     $     (sinthetatau)**2.d0 * a0mtau2)) +
     $     (g*g) * 0.5d0 * (1.5d0*a0muL + 1.5d0*a0mcL +
     $     1.5d0*((costhetat**2.d0)*a0mtL + 
     $     (sinthetat**2.d0)*a0mtR) -
     $     1.5d0*a0mdL - 1.5d0*a0msL -
     $     1.5d0*((costhetab**2.d0)*a0mbL +
     $     (sinthetab**2.d0)*a0mbR) +
     $     0.5d0*(a0msnu1 + a0msnu2 + a0msnu3) -
     $     0.5d0*(a0meL + a0mmuR +
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2)) +
     $     (gp**2.d0)*0.25d0*(ynuL**2.d0) * (a0msnu1) +
     $     (gp**2.d0)*0.25d0 * ynuL * (3.d0 * yuL * (a0muL + a0mcL + 
     $     (costhetat**2.d0)*a0mtL +(sinthetat**2.d0)*a0mtR) +
     $     3.d0 * yuR * (a0muR + a0mcR + 
     $     (sinthetat**2.d0)*a0mtL +(costhetat**2.d0)*a0mtR) +
     $     3.d0 * ydL * (a0mdL + a0msL + 
     $     (costhetab**2.d0)*a0mbL + (sinthetab**2.d0)*a0mbR) +
     $     3.d0 * ydR * (a0mdR + a0msR + 
     $     (sinthetab**2.d0)*a0mbL + (costhetab**2.d0)*a0mbR) +
     $     yeL * (a0meL + a0mmuL + 
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2) +
     $     yeR * (a0meR + a0mmuR + 
     $     (sinthetatau**2.d0)*a0mtau1 +(costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0msnu1 + a0msnu2 + a0msnu3))


      
!----------------------------------------------------------

      call funcg(p,mchargino(1),mtau,q,gmchargino10)
      call funcg(p,mchargino(2),mtau,q,gmchargino20)
      call funcg(p,nmneut(1),0.d0,q,gmneut10)
      call funcg(p,nmneut(2),0.d0,q,gmneut20)
      call funcg(p,nmneut(3),0.d0,q,gmneut30)
      call funcg(p,nmneut(4),0.d0,q,gmneut40)

      call b0(p,mchargino(1),mtau,q,b0mch1mtau)
      call b0(p,mchargino(2),mtau,q,b0mch2mtau)

      charginoterm =  fChsnuLL(1)*gmchargino10
     $     -  2.d0*mchargino(1)*mtau*gChsnuLL(1)*b0mch1mtau
     $     + fChsnuLL(2)*gmchargino20
     $     -  2.d0*mchargino(2)*mtau*gChsnuLL(2)*b0mch2mtau


      neutralinoterm = fChi0nusnuLL(1)*gmneut10
     $     + fChi0nusnuLL(2)*gmneut20   
     $     + fChi0nusnuLL(3)*gmneut30
     $     + fChi0nusnuLL(4)*gmneut40

!-------------------------------------------------------------

      pitsnu = 0.d0
      pitsnu = (1.d0/(16.d0*(pi)**2.d0)) *  
     $     (sbotterm + charginoterm + neutralinoterm + higgsterm)



      tsnu = dsqrt((snu(1)*snu(1)) - pitsnu)	  

      RETURN

      END SUBROUTINE pitausnu

C===================================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     
C     pisnumu is checked on 01/06/2010 @ 15:30.
C     
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      SUBROUTINE pimulsnu(p,q,g,gp,tanbeta,yeRG,AERG,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,
     $     mhpmsq,mA0sq,sgnmu,modmu,ON,OCL,OCR,vev1,musnu)
      
 
      IMPLICIT NONE 

      INTEGER i,j
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION AERG(3,3),SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yeRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3),nmneut(4)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g

      double precision vev1,beta
 
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau
      DOUBLE PRECISION pimusnu,musnu
 
      DOUBLE PRECISION dnu(4),dnd(4),cn(4)

      
      DOUBLE PRECISION chargino,neutralino
      DOUBLE PRECISION stops,sbottom,higgs
      data stops/ 1 * 0.d0/,sbottom/ 1 * 0.d0/,higgs/ 1 * 0.d0/
      data chargino/ 1 * 0.d0/, neutralino/1 * 0.d0/

      
      DOUBLE PRECISION intm1(2), intm2(2)

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      
      DOUBLE PRECISION rthetamu(2,2),ralpha(2,2)

      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz

      DOUBLE PRECISION a0mdL,a0mdR,a0msL,a0msR,a0msnu2
      DOUBLE PRECISION a0meL,a0meR,a0mmuL,a0mmuR,a0msnu1,a0msnu3
      DOUBLE PRECISION a0mcR,a0mcL,a0muL,a0muR,a0mbL,a0mbR,a0mtL,a0mtR
      DOUBLE PRECISION a0mtau1,a0mtau2

      DOUBLE PRECISION fmmu1MW,fmmu2MW,fsnu2MZ
      DOUBLE PRECISION gmchargino10,gmchargino20
      DOUBLE PRECISION gmneut10,gmneut20,gmneut30,gmneut40

      DOUBLE PRECISION aPsi0snul(4), bPsi0snul(4) 
      data aPsi0snul/ 4 * 0.d0/,bPsi0snul/ 4 * 0.d0/

      DOUBLE PRECISION aChi0snul(4), bChi0snul(4)
      data aChi0snul/ 4 * 0.d0/,bChi0snul/ 4 * 0.d0/
      
      DOUBLE PRECISION gChi0nusnuLL(4), fChi0nusnuLL(4)
      data gChi0nusnuLL/ 4 * 0.d0/,fChi0nusnuLL/ 4 * 0.d0/


      DOUBLE PRECISION bPsicsnul(2), aPsicsnul(2)
      data bPsicsnul/ 2 * 0.d0/,aPsicsnul/ 2 * 0.d0/

      DOUBLE PRECISION aPsicCStaul(2),aChicsnur(2), aChicsnul(2)
      data aPsicCStaul/ 2 * 0.d0/,aChicsnur/ 2 * 0.d0/,
     $     aChicsnul/ 2 * 0.d0/

      DOUBLE PRECISION bChicsnul(2),bChicsnur(2)
      data bChicsnul/ 2 * 0.d0/,bChicsnur/ 2 * 0.d0/

      DOUBLE PRECISION aChsnu(2, 2), bChsnu(2, 2)
      data aChsnu/ 4 * 0.d0/,bChsnu/ 4 * 0.d0/

      DOUBLE PRECISION fChsnuLL(2), gChsnuLL(2)
      data fChsnuLL/ 2 * 0.d0/,gChsnuLL/ 2 * 0.d0/
     
      DOUBLE PRECISION lsnulsnulr(4),lHsnulsnu12(4)
      DOUBLE PRECISION lHcsnulsmulr(2, 2),lHcsnulsmu12(2, 2)
      DOUBLE PRECISION b0mHusnu2,b0mhsnu2,b0mzsnu2,b0mAsnu2
      DOUBLE PRECISION b0mHpmmuL,b0mHpmmuR,b0mwmmu1,b0mwmmu2
      DOUBLE PRECISION b0mch1mmu,b0mch2mmu
      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2

      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae
 

      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d
  

      include 'stdinputs.h'
C------------------------------------------------------------

      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      MZ = MZpole

      
      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw   = dsqrt(sinsqthw)
      thw      = dasin(sinthw)
      costhw   = dcos(thw)
      cos2thw  = dcos(2.d0*thw)
      
     
      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))
      
      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))
      
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

!--------------------------------------------

      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)

      mt1 = mtL
      mt2 = mtR

      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR

C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta     = datan(tanbeta)
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      salpha   = dsin(alpha)
      calpha   = dcos(alpha)
      cosbeta  = dcos(datan(tanbeta))
      cos2beta = dcos(2.d0*datan(tanbeta))
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)
!-----------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0

!-------------------------------------------------------------

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

!-----------------------------------------------------------

      lsnulsnulr(1) =  g * mz * gnuL * cosbeta/costhw
      lsnulsnulr(2) = -g * mz * gnuL * sinbeta/costhw
      lsnulsnulr(3) = 0.d0
      lsnulsnulr(4) = 0.d0

C     Mix CP-even Higgses up

      intm1(1) = lsnulsnulr(1)
      intm1(2) = lsnulsnulr(2)
      
      call rmat2d(alpha,ralpha)

      intm2(1) = ralpha(1,1)*intm1(1) + ralpha(1,2)*intm1(2)
      intm2(2) = ralpha(2,1)*intm1(1) + ralpha(2,2)*intm1(2)
      
      lHsnulsnu12(1) = intm2(1)
      lHsnulsnu12(2) = intm2(2)
      lHsnulsnu12(3) = 0.d0
      lHsnulsnu12(4) = 0.d0

!-------------------------------------------------------------------
C     Charged Higgs Feynman rules
!     (H+ G+, L R) basis

      lHcsnulsmulr(1, 1) = ((g*MW*dsin(2.d0*beta) -
     $     yeRG(2,2)**2.d0 * vev1 * sinbeta)/dsqrt(2.d0))

      lHcsnulsmulr(1, 2) = (- sgnmu*modmu*yeRG(2,2)*cosbeta - 
     $     yeRG(2,2) * AERG(2,2) * sinbeta)

      lHcsnulsmulr(2, 1) = (-g*MW*dcos(2.d0*beta) +
     $     yeRG(2,2)**2.d0 * vev1 * cosbeta)/dsqrt(2.d0)

      lHcsnulsmulr(2, 2) = (-yeRG(2,2)*sgnmu*modmu*sinbeta +
     $     yeRG(2,2) * AERG(2,2) * cosbeta)

      intm1(1) = lHcsnulsmulr(1, 1)
      intm1(2) = lHcsnulsmulr(1, 2)

      call rmat2d(thetamu,rthetamu)

      intm2(1) = rthetamu(1,1)*intm1(1) + rthetamu(1,2)*intm1(2)
      intm2(2) = rthetamu(2,1)*intm1(1) + rthetamu(2,2)*intm1(2)
      
      lHcsnulsmu12(1, 1) = intm2(1)
      lHcsnulsmu12(1, 2) = intm2(2)

      intm1(1) = lHcsnulsmulr(2, 1)
      intm1(2) = lHcsnulsmulr(2, 2)

      intm2(1) = rthetamu(1,1)*intm1(1) + rthetamu(1,2)*intm1(2)
      intm2(2) = rthetamu(2,1)*intm1(1) + rthetamu(2,2)*intm1(2)

      lHcsnulsmu12(2, 1) = intm2(1)
      lHcsnulsmu12(2, 2) = intm2(2)

!-------------------------------------------------------------------------------
C     Neutralino Feynman rules
!-------------------------------------------------------------------------------
      
      bPsi0snul(1) = gp*ynuL/dsqrt(2.d0)
      bPsi0snul(2) = g / dsqrt(2.d0)
      bPsi0snul(3) = 0.d0
      bPsi0snul(4) = 0.d0

      aPsi0snul(1) = 0.d0
      aPsi0snul(2) = 0.d0
      aPsi0snul(3) = 0.d0
      aPsi0snul(4) = yeRG(2,2) 

!------------

      aChi0snul(1) = 0.d0
      aChi0snul(2) = 0.d0
      aChi0snul(3) = 0.d0
      aChi0snul(4) = 0.d0

      bChi0snul(1) = 0.d0
      bChi0snul(2) = 0.d0
      bChi0snul(3) = 0.d0
      bChi0snul(4) = 0.d0

!------------

      loopi: DO i = 1, 4
      loopj: DO j = 1, 4

      aChi0snul(i) = aChi0snul(i) + ON(i,j)*aPsi0snul(j)
      bChi0snul(i) = bChi0snul(i) + ON(i,j)*bPsi0snul(j)

      ENDDO loopj
      ENDDO loopi
!-------------


      lcosnui: DO i = 1, 4

      fChi0nusnuLL(i) = (aChi0snul(i) * aChi0snul(i) + 
     $     bChi0snul(i) * bChi0snul(i))
      gChi0nusnuLL(i) = (bChi0snul(i) * aChi0snul(i) + 
     $     bChi0snul(i) * aChi0snul(i))

      ENDDO lcosnui
      
!----------------------------------------------------------------------------
C     Chargino Feynman Rules
!----------------------------------------------------------------------------


      aPsicsnul(1) = g
      aPsicsnul(2) = 0.d0

      bPsicsnul(1) =  0.d0
      bPsicsnul(2) = - yeRG(2,2)
      
!--------------
      aChicsnul(1) = 0.d0
      aChicsnul(2) = 0.d0

      bChicsnul(1) = 0.d0
      bChicsnul(2) = 0.d0


      
!--------------

      loopchi: DO i = 1, 2
      loopchj: DO j = 1, 2

      aChicsnul(i) = aChicsnul(i) + OCR(i,j)*aPsicsnul(j)
      bChicsnul(i) = bChicsnul(i) + OCL(i,j)*bPsicsnul(j)

      ENDDO loopchj
      ENDDO loopchi

!---------------------

      lcosnuchi: DO i = 1, 2
      
      fChsnuLL(i) = (aChicsnul(i) * aChicsnul(i) +
     $     bChicsnul(i) * bChicsnul(i))
      gChsnuLL(i) = (bChicsnul(i)* aChicsnul(i) +
     $     aChicsnul(i)* bChicsnul(i))
      
      ENDDO lcosnuchi
!-----------------------------------------------------------
      
      call a0(mA0,q,a0mA)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)
      call a0(mHpm,q,a0mHpm)
      call a0(MW,q,a0Mw)
      call a0(MZ,q,a0Mz)

      call a0(muR,q,a0muR)
      call a0(muL,q,a0muL)
      call a0(mcR,q,a0mcR)
      call a0(mcL,q,a0mcL)
      call a0(mtL,q,a0mtL)
      call a0(mtR,q,a0mtR)

      call a0(msR,q,a0msR)
      call a0(msL,q,a0msL)
      call a0(mdR,q,a0mdR)
      call a0(mdL,q,a0mdL)
      call a0(mbL,q,a0mbL)
      call a0(mbR,q,a0mbR)

      call a0(meR,q,a0meR)
      call a0(meL,q,a0meL)
      call a0(mmuR,q,a0mmuR)
      call a0(mmuL,q,a0mmuL)
      call a0(mtau1,q,a0mtau1)
      call a0(mtau2,q,a0mtau2)

      call a0(snu(1),q,a0msnu1)
      call a0(snu(2),q,a0msnu2)
      call a0(snu(3),q,a0msnu3)

      call f(p,mmuL,MW,q,fmmu1Mw)
      call f(p,mmuR,MW,q,fmmu2Mw)
      call f(p,snu(2),MZ,q,fsnu2mz)

      call b0(p,mHu0,snu(2),q,b0mHusnu2)
      call b0(p,mh0,snu(2),q,b0mhsnu2)
      call b0(p,MZ,snu(2),q,b0mzsnu2)
      call b0(p,mA0,snu(2),q,b0mAsnu2)

      call b0(p,mmuL,mHpm,q,b0mHpmmuL)
      call b0(p,mmuR,mHpm,q,b0mHpmmuR)
      call b0(p,mmuL,MW,q,b0mwmmu1)
      call b0(p,mmuR,MW,q,b0mwmmu2)

      call b0(p,mchargino(1),mmu,q,b0mch1mmu)
      call b0(p,mchargino(2),mmu,q,b0mch2mmu)

!-----------------------------------------------------------------------------------
! Corrections begin
!-----------------------------------------------------------------------------------

      sbottom = 0.d0
      sbottom = (yeRG(2,2)**2.d0 * a0mmuR)

      higgs = 0.d0
      higgs = 0.d0 + 0.5d0 * 
     $     (((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0) * cn(1)) * a0mHu +
     $     ((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0) * cn(2)) * a0mh +
     $     ((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0) * cn(3)) * a0mz +
     $     ((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0) * cn(4)) * a0mA) 


      higgs = higgs + ((yeRG(2,2)**2.d0) * dnu(3) + (g)**2.d0 * 
     $     ((gnuL*0.5d0/costhw**2.d0) - 0.5d0) * cn(3)) * a0mHpm +
     $     ((yeRG(2,2)**2.d0) * dnu(4) + (g)**2.d0 *  
     $     ((gnuL*0.5d0/costhw**2.d0) - 0.5d0) * cn(4)) * a0MW


      higgs = higgs + (lHsnulsnu12(1))**2.d0 * b0mHusnu2 +
     $     (lHsnulsnu12(2))**2.d0 * b0mhsnu2 +
     $     (lHsnulsnu12(3))**2.d0 * b0mzsnu2 +
     $     (lHsnulsnu12(4))**2.d0 * b0mAsnu2
      

      higgs = higgs + (lHcsnulsmu12(1,1))**2.d0 * b0mHpmmuL +
     $     (lHcsnulsmu12(1,2))**2.d0 * b0mHpmmuR +
     $     (lHcsnulsmu12(2,1))**2.d0 * b0mwmmu1 +
     $     (lHcsnulsmu12(2,2))**2.d0 * b0mwmmu2

!-------------------------------------------------


      higgs = higgs + 4.d0 * ((g*gnuL)**2.d0/(costhw)**2.d0) * a0mz + 
     $     2.d0*(g)**2.d0*a0mw + 
     $     (g*gnuL/costhw)**2.d0 * fsnu2mz +
     $     (g)**2.d0 * 0.5d0 * fmmu1MW + 
     $     (g)**2.d0 * 0.25 * (a0msnu2 + 2.d0 * a0mmuL) + 
     $     (g*g) * 0.5d0 * (1.5d0*a0muL + 1.5d0*a0mcL +
     $     1.5d0*((costhetat**2.d0)*a0mtL + 
     $     (sinthetat**2.d0)*a0mtR) -
     $     1.5d0*a0mdL - 1.5d0*a0msL -
     $     1.5d0*((costhetab**2.d0)*a0mbL +
     $     (sinthetab**2.d0)*a0mbR) +
     $     0.5d0*(a0msnu1 + a0msnu2 + a0msnu3) -
     $     0.5d0*(a0meL + a0mmuR +
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2)) +
     $     (gp**2.d0) * 0.25d0 * (ynuL**2.d0) * (a0msnu2) +
     $     (gp**2.d0) * 0.25d0 * ynuL * (3.d0 * yuL * (a0muL + a0mcL + 
     $     (costhetat**2.d0)*a0mtL +(sinthetat**2.d0)*a0mtR) +
     $     3.d0 * yuR * (a0muR + a0mcR + 
     $     (sinthetat**2.d0)*a0mtL +(costhetat**2.d0)*a0mtR) +
     $     3.d0 * ydL * (a0mdL + a0msL + 
     $     (costhetab**2.d0)*a0mbL + (sinthetab**2.d0)*a0mbR) +
     $     3.d0 * ydR * (a0mdR + a0msR + 
     $     (sinthetab**2.d0)*a0mbL + (costhetab**2.d0)*a0mbR) +
     $     yeL * (a0meL + a0mmuL + 
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2) +
     $     yeR * (a0meR + a0mmuR + 
     $     (sinthetatau**2.d0)*a0mtau1 +(costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0msnu1 + a0msnu2 + a0msnu3))

     
!----------------------------------------------------------

      call funcg(p,mchargino(1),mmu,q,gmchargino10)
      call funcg(p,mchargino(2),mmu,q,gmchargino20)
      call funcg(p,mneut(1),0.d0,q,gmneut10)
      call funcg(p,mneut(2),0.d0,q,gmneut20)
      call funcg(p,mneut(3),0.d0,q,gmneut30)
      call funcg(p,mneut(4),0.d0,q,gmneut40)


      chargino = 0.d0
      chargino = fChsnuLL(1)*gmchargino10 -
     $     2.d0*mchargino(1)*mmu*gChsnuLL(1)*b0mch1mmu +
     $     fChsnuLL(2)*gmchargino20 -
     $     2.d0*mchargino(2)*mmu*gChsnuLL(2)*b0mch2mmu

      neutralino = 0.d0
      neutralino = fChi0nusnuLL(1)*gmneut10 +
     $     fChi0nusnuLL(2)*gmneut20 +   
     $     fChi0nusnuLL(3)*gmneut30 + 
     $     fChi0nusnuLL(4)*gmneut40  

!-------------------------------------------------------------

      pimusnu = (1.d0/(16.d0*(pi)**2.d0)) * 
     $     (sbottom + higgs + chargino + neutralino)


      musnu = dsqrt((snu(2)*snu(2)) - pimusnu)	  

      RETURN

      END SUBROUTINE pimulsnu
      
C===================================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     
C     pisnuel is checked on 01/06/2010 @ 15:30.
C     
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
 
       SUBROUTINE pielsnu(p,q,g,gp,tanbeta,yeRG,AERG,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,
     $     mhpmsq,mA0sq,sgnmu,modmu,ON,OCL,OCR,vev1,elsnu)
      

      IMPLICIT NONE 

      INTEGER i,j
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION AERG(3,3),SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,yeRG(3,3)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3),nmneut(4)

      DOUBLE PRECISION gp,calpha,cos2beta,cosbeta,salpha,sin2beta
      DOUBLE PRECISION sinbeta,tan2beta
      DOUBLE PRECISION sinsqthw,g

      double precision vev1,beta
 
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION piesnu,elsnu
 
      DOUBLE PRECISION dnu(4),dnd(4),cn(4)

      
      DOUBLE PRECISION  chargino, neutralino
      DOUBLE PRECISION stops,sbottom,higgs
      data stops/ 1 * 0.d0/,sbottom/ 1 * 0.d0/,higgs/ 1 * 0.d0/
      data chargino/ 1 * 0.d0/, neutralino/1 * 0.d0/
      
      DOUBLE PRECISION  intm1(2), intm2(2)

      DOUBLE PRECISION mchargino(2),ON(4,4),OCL(2,2),OCR(2,2)
      
      DOUBLE PRECISION rthetae(2,2),ralpha(2,2)


      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz

      DOUBLE PRECISION a0mdL,a0mdR,a0msL,a0msR,a0msnu2
      DOUBLE PRECISION a0meL,a0meR,a0mmuL,a0mmuR,a0msnu1,a0msnu3
      DOUBLE PRECISION a0mcR,a0mcL,a0muL,a0muR,a0mbL,a0mbR,a0mtL,a0mtR
      DOUBLE PRECISION a0mtau1,a0mtau2

      DOUBLE PRECISION fme1MW,fme2MW,fsnu3MZ
      DOUBLE PRECISION gmchargino10,gmchargino20
      DOUBLE PRECISION gmneut10,gmneut20,gmneut30,gmneut40

      DOUBLE PRECISION aPsi0snul(4), bPsi0snul(4) 
      data aPsi0snul/ 4 * 0.d0/,bPsi0snul/ 4 * 0.d0/

      DOUBLE PRECISION aChi0snul(4), bChi0snul(4)
      data aChi0snul/ 4 * 0.d0/,bChi0snul/ 4 * 0.d0/
      
      DOUBLE PRECISION gChi0nusnuLL(4), fChi0nusnuLL(4)
      data gChi0nusnuLL/ 4 * 0.d0/,fChi0nusnuLL/ 4 * 0.d0/


      DOUBLE PRECISION bPsicsnul(2), aPsicsnul(2)
      data bPsicsnul/ 2 * 0.d0/,aPsicsnul/ 2 * 0.d0/

      DOUBLE PRECISION aPsicCStaul(2),aChicsnur(2), aChicsnul(2)
      data aPsicCStaul/ 2 * 0.d0/,aChicsnur/ 2 * 0.d0/,
     $     aChicsnul/ 2 * 0.d0/

      DOUBLE PRECISION bChicsnul(2),bChicsnur(2)
      data bChicsnul/ 2 * 0.d0/,bChicsnur/ 2 * 0.d0/

      DOUBLE PRECISION aChsnu(2, 2), bChsnu(2, 2)
      data aChsnu/ 4 * 0.d0/,bChsnu/ 4 * 0.d0/

      DOUBLE PRECISION fChsnuLL(2), gChsnuLL(2)
      data fChsnuLL/ 2 * 0.d0/,gChsnuLL/ 2 * 0.d0/
     
      DOUBLE PRECISION  lsnulsnulr(4),lHsnulsnu12(4)
      DOUBLE PRECISION lHcsnulseller(2, 2),lHcsnulsel12(2, 2)
      data lsnulsnulr/ 4 * 0.d0/, lHsnulsnu12/ 4 * 0.d0/
      data lHcsnulseller/ 4 * 0.d0/, lHcsnulsel12/ 4 * 0.d0/

      DOUBLE PRECISION b0mHusnu3,b0mhsnu3,b0mzsnu3,b0mAsnu3
      DOUBLE PRECISION b0mHpmeL,b0mHpmeR,b0mwme1,b0mwme2

      double precision b0mch1me,b0mch2me
 
      DOUBLE PRECISION mt1,mt2,mb1,mb2,mtau1,mtau2

      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

 
      EXTERNAL funcg,f,b0,a0,rmat2d,mat3prod2d
  

      include 'stdinputs.h'

C------------------------------------------------------------

      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      MZ = MZpole

      elsnu = 0.d0
      
      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw   = dsqrt(sinsqthw)
      thw      = dasin(sinthw)
      costhw   = dcos(thw)
      cos2thw  = dcos(2.d0*thw)
      

      muL = dsqrt(SUegg(6)) 
      muR = dsqrt(SUegg(5))
      mcL = dsqrt(SUegg(4))
      mcR = dsqrt(SUegg(3))
      mtL = dsqrt(SUegg(2))
      mtR = dsqrt(SUegg(1))
      
      mdL = dsqrt(SDegg(6)) 
      mdR = dsqrt(SDegg(5))
      msL = dsqrt(SDegg(4))
      msR = dsqrt(SDegg(3))
      mbL = dsqrt(SDegg(2))
      mbR = dsqrt(SDegg(1))

      meL   = dsqrt(SLegg(6))
      meR   = dsqrt(SLegg(5))
      mmuL  = dsqrt(SLegg(4))
      mmuR  = dsqrt(SLegg(3))
      mtauL = dsqrt(SLegg(2))
      mtauR = dsqrt(SLegg(1))

      
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

!----------------------------------------

      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)


      mt1 = mtL
      mt2 = mtR

      mb1 = mbL
      mb2 = mbR

      mtau1 = mtauL
      mtau2 = mtauR

C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta     = datan(tanbeta)
      tan2beta = dtan(2.d0*datan(tanbeta))
      alpha    = 0.5d0*datan(((mA0sq + MZ*MZ)/(mA0sq - MZ*MZ))*tan2beta)
      salpha   = dsin(alpha)
      calpha   = dcos(alpha)
      cosbeta  = dcos(beta)
      cos2beta = dcos(2.d0 * beta)
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)

!-----------------------------------------------------------

      cn(1) = -dcos(2.d0*alpha)
      cn(2) =  dcos(2.d0*alpha)
      cn(3) = -dcos(2.d0*beta)
      cn(4) =  dcos(2.d0*beta)
      
      dnu(1) = dsin(alpha)**2.d0
      dnu(2) = dcos(alpha)**2.d0
      dnu(3) = dsin(beta)**2.d0
      dnu(4) = dcos(beta)**2.d0

      dnd(1) = dcos(alpha)**2.d0
      dnd(2) = dsin(alpha)**2.d0
      dnd(3) = dcos(beta)**2.d0
      dnd(4) = dsin(beta)**2.d0

!-------------------------------------------------------------

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

!-----------------------------------------------------------

      lsnulsnulr(1) =  g * mz * gnuL * cosbeta/costhw
      lsnulsnulr(2) = -g * mz * gnuL * sinbeta/costhw
      lsnulsnulr(3) = 0.d0
      lsnulsnulr(4) = 0.d0

C     Mix CP-even Higgses up

      intm1(1) = lsnulsnulr(1)
      intm1(2) = lsnulsnulr(2)
      
      call rmat2d(alpha,ralpha)

      intm2(1) = ralpha(1,1)*intm1(1) + ralpha(1,2)*intm1(2)
      intm2(2) = ralpha(2,1)*intm1(1) + ralpha(2,2)*intm1(2)
      
      lHsnulsnu12(1) = intm2(1)
      lHsnulsnu12(2) = intm2(2)
      lHsnulsnu12(3) = 0.d0
      lHsnulsnu12(4) = 0.d0

!-------------------------------------------------------------------
C     Charged Higgs Feynman rules
!     (H+ G+, L R) basis
      
      lHcsnulseller(1, 1) = (g*MW*dsin(2.d0*beta) -
     $     yeRG(1,1)**2.d0 * vev1 * sinbeta)/dsqrt(2.d0)

      lHcsnulseller(1, 2) = (-sgnmu*modmu*yeRG(1,1)*cosbeta - 
     $     yeRG(1,1)* AERG(1,1) * sinbeta)

      lHcsnulseller(2, 1) = (-g*MW*dcos(2.d0*beta) +
     $     yeRG(1,1)**2.d0 * vev1 * cosbeta)/dsqrt(2.d0)

      lHcsnulseller(2, 2) = (-yeRG(1,1)*sgnmu*modmu*sinbeta +
     $     yeRG(1,1) * AERG(1,1) * cosbeta) 

      intm1(1) = lHcsnulseller(1, 1)
      intm1(2) = lHcsnulseller(1, 2)

      call rmat2d(thetae,rthetae)

      intm2(1) = rthetae(1,1)*intm1(1) + rthetae(1,2)*intm1(2)
      intm2(2) = rthetae(2,1)*intm1(1) + rthetae(2,2)*intm1(2)
      
      lHcsnulsel12(1, 1) = intm2(1)
      lHcsnulsel12(1, 2) = intm2(2)

      intm1(1) = lHcsnulseller(2, 1)
      intm1(2) = lHcsnulseller(2, 2)

      intm2(1) = rthetae(1,1)*intm1(1) + rthetae(1,2)*intm1(2)
      intm2(2) = rthetae(2,1)*intm1(1) + rthetae(2,2)*intm1(2)

      lHcsnulsel12(2, 1) = intm2(1)
      lHcsnulsel12(2, 2) = intm2(2)

!-------------------------------------------------------------------------------
C     Neutralino Feynman rules
!-------------------------------------------------------------------------------
      
      bPsi0snul(1) = gp*ynuL/dsqrt(2.d0)
      bPsi0snul(2) = g / dsqrt(2.d0)
      bPsi0snul(3) = 0.d0
      bPsi0snul(4) = 0.d0

      aPsi0snul(1) = 0.d0
      aPsi0snul(2) = 0.d0
      aPsi0snul(3) = 0.d0
      aPsi0snul(4) = yeRG(1,1) 

!------------

      aChi0snul(1) = 0.d0
      aChi0snul(2) = 0.d0
      aChi0snul(3) = 0.d0
      aChi0snul(4) = 0.d0

      bChi0snul(1) = 0.d0
      bChi0snul(2) = 0.d0
      bChi0snul(3) = 0.d0
      bChi0snul(4) = 0.d0



!------------

      loopi: DO i = 1, 4
      loopj: DO j = 1, 4

      aChi0snul(i) = aChi0snul(i) + ON(i,j)*aPsi0snul(j)
      bChi0snul(i) = bChi0snul(i) + ON(i,j)*bPsi0snul(j)

      ENDDO loopj
      ENDDO loopi

!-------------


      lcosnui: DO i = 1, 4

      fChi0nusnuLL(i) = (aChi0snul(i) * aChi0snul(i) + 
     $     bChi0snul(i) * bChi0snul(i))
      gChi0nusnuLL(i) = (bChi0snul(i) * aChi0snul(i) + 
     $     bChi0snul(i) * aChi0snul(i))

      ENDDO lcosnui
      
!----------------------------------------------------------------------------
C     Chargino Feynman Rules
!----------------------------------------------------------------------------

      aPsicsnul(1) = g
      aPsicsnul(2) = 0.d0

      bPsicsnul(1) =  0.d0
      bPsicsnul(2) = - yeRG(1,1)
      
!--------------
      aChicsnul(1) = 0.d0
      aChicsnul(2) = 0.d0

      bChicsnul(1) = 0.d0
      bChicsnul(2) = 0.d0


      
!--------------
      loopchi: DO i = 1, 2
      loopchj: DO j = 1, 2

      aChicsnul(i) = aChicsnul(i) + OCR(i,j)*aPsicsnul(j)
      bChicsnul(i) = bChicsnul(i) + OCL(i,j)*bPsicsnul(j)

      ENDDO loopchj
      ENDDO loopchi
!---------------------

      lcosnuchi: DO i = 1, 2
      
      fChsnuLL(i) = (aChicsnul(i) * aChicsnul(i) +
     $     bChicsnul(i) * bChicsnul(i))
      gChsnuLL(i) = (bChicsnul(i)* aChicsnul(i) +
     $     aChicsnul(i)* bChicsnul(i))
      
      ENDDO lcosnuchi
!-----------------------------------------------------------

       call a0(mA0,q,a0mA)
       call a0(mh0,q,a0mh)
       call a0(mHu0,q,a0mHu)
       call a0(mHpm,q,a0mHpm)
       call a0(MW,q,a0Mw)
       call a0(MZ,q,a0Mz)

       call a0(muR,q,a0muR)
       call a0(muL,q,a0muL)
       call a0(mcR,q,a0mcR)
       call a0(mcL,q,a0mcL)
       call a0(mtL,q,a0mtL)
       call a0(mtR,q,a0mtR)

       call a0(msR,q,a0msR)
       call a0(msL,q,a0msL)
       call a0(mdR,q,a0mdR)
       call a0(mdL,q,a0mdL)
       call a0(mbL,q,a0mbL)
       call a0(mbR,q,a0mbR)

       call a0(meR,q,a0meR)
       call a0(meL,q,a0meL)
       call a0(mmuR,q,a0mmuR)
       call a0(mmuL,q,a0mmuL)
       call a0(mtau1,q,a0mtau1)
       call a0(mtau2,q,a0mtau2)

       call a0(snu(1),q,a0msnu1)
       call a0(snu(2),q,a0msnu2)
       call a0(snu(3),q,a0msnu3)

       call f(p,meL,MW,q,fme1Mw)
       call f(p,meR,MW,q,fme2Mw)
       call f(p,snu(3),MZ,q,fsnu3Mz)

      call b0(p,mHu0,snu(3),q,b0mHusnu3)
      call b0(p,mh0,snu(3),q,b0mhsnu3)
      call b0(p,MZ,snu(3),q,b0mzsnu3)
      call b0(p,mA0,snu(3),q,b0mAsnu3)

      call b0(p,meL,mHpm,q,b0mHpmeL)
      call b0(p,meR,mHpm,q,b0mHpmeR)
      call b0(p,meL,MW,q,b0mwme1)
      call b0(p,meR,MW,q,b0mwme2)

      call b0(p,mchargino(1),me,q,b0mch1me)
      call b0(p,mchargino(2),me,q,b0mch2me)


!-----------------------------------------------------------------------------------
!     Corrections Begin
!-----------------------------------------------------------------------------------
      sbottom = 0.d0
      sbottom = yeRG(1,1)**2.d0 * (a0meR) !<==  changed from nue_L --> nue_R

      higgs = 0.d0
      higgs = 0.d0 + 0.5d0 *  
     $     (((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0) * cn(1)) * a0mHu +
     $     ((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0)  * cn(2)) * a0mh +
     $     ((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0)  * cn(3)) * a0mz +
     $     ((-(g)**2.d0*gnuL*0.5d0/(costhw)**2.d0)  * cn(4)) * a0mA) 



      higgs = higgs + ((yeRG(1,1)**2.d0) * dnu(3) + (g)**2.d0 * 
     $     ((gnuL*0.5d0/costhw**2.d0) - 0.5d0) * cn(3)) * a0mHpm +
     $     ((yeRG(1,1)**2.d0) * dnu(4) + (g)**2.d0 * 
     $     ((gnuL*0.5d0/costhw**2.d0) - 0.5d0) * cn(4)) * a0mw


      higgs = higgs + (lHsnulsnu12(1))**2.d0 * b0mHusnu3 +
     $     (lHsnulsnu12(2))**2.d0 * b0mhsnu3 +
     $     (lHsnulsnu12(3))**2.d0 * b0mzsnu3 +
     $     (lHsnulsnu12(4))**2.d0 * b0mAsnu3
      

      higgs = higgs + (lHcsnulsel12(1,1))**2.d0 * b0mHpmeL +
     $     (lHcsnulsel12(1,2))**2.d0 *b0mHpmeR +
     $     (lHcsnulsel12(2,1))**2.d0 * b0mwme1 +
     $     (lHcsnulsel12(2,2))**2.d0 * b0mwme2
      

!-----------------------------------------------------------------------------------

      higgs = higgs + 4.d0 * ((g*gnuL)**2.d0/(costhw)**2.d0) * a0mz + 
     $     2.d0*(g)**2.d0*a0mw + 
     $     (g*gnuL/costhw)**2.d0 * fsnu3mz +
     $     (g)**2.d0 * 0.5d0 * fme1MW + 
     $     (g)**2.d0 * 0.25 * (a0msnu3 + 2.d0 * a0meL) + 
     $     (g*g) * 0.5d0 * (1.5d0*a0muL + 1.5d0*a0mcL +
     $     1.5d0*((costhetat**2.d0)*a0mtL + 
     $     (sinthetat**2.d0)*a0mtR) -
     $     1.5d0*a0mdL - 1.5d0*a0msL -
     $     1.5d0*((costhetab**2.d0)*a0mbL +
     $     (sinthetab**2.d0)*a0mbR) +
     $     0.5d0 * (a0msnu1 + a0msnu2 + a0msnu3) -
     $     0.5d0 * (a0meL + a0mmuR +
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2)) +
     $     (gp**2.d0) * 0.25d0 * (ynuL**2.d0) * (a0msnu3) +
     $     (gp**2.d0) * 0.25d0 * ynuL * (3.d0 * yuL * (a0muL + a0mcL + 
     $     (costhetat**2.d0)*a0mtL +(sinthetat**2.d0)*a0mtR) +
     $     3.d0 * yuR * (a0muR + a0mcR + 
     $     (sinthetat**2.d0)*a0mtL +(costhetat**2.d0)*a0mtR) +
     $     3.d0 * ydL * (a0mdL + a0msL + 
     $     (costhetab**2.d0)*a0mbL + (sinthetab**2.d0)*a0mbR) +
     $     3.d0 * ydR * (a0mdR + a0msR + 
     $     (sinthetab**2.d0)*a0mbL + (costhetab**2.d0)*a0mbR) +
     $     yeL * (a0meL + a0mmuL + 
     $     (costhetatau**2.d0)*a0mtau1 + (sinthetatau**2.d0)*a0mtau2) +
     $     yeR * (a0meR + a0mmuR + 
     $     (sinthetatau**2.d0)*a0mtau1 +(costhetatau**2.d0)*a0mtau2) +
     $     ynuL * (a0msnu1 + a0msnu2 + a0msnu3))

     
!---------------------------------------------------

      call funcg(p,mchargino(1),me,q,gmchargino10)
      call funcg(p,mchargino(2),me,q,gmchargino20)
      call funcg(p,mneut(1),me,q,gmneut10)
      call funcg(p,mneut(2),me,q,gmneut20)
      call funcg(p,mneut(3),me,q,gmneut30)
      call funcg(p,mneut(4),me,q,gmneut40)

      chargino = 0.d0

      chargino = fChsnuLL(1)*gmchargino10 -
     $     2.d0 * mchargino(1) * me * gChsnuLL(1) * b0mch1me +
     $     fChsnuLL(2) * gmchargino20 -
     $     2.d0 * mchargino(2) * me * gChsnuLL(2) * b0mch2me

      neutralino = 0.d0
      neutralino =  fChi0nusnuLL(1)*gmneut10 +
     $     fChi0nusnuLL(2)*gmneut20 +    
     $     fChi0nusnuLL(3)*gmneut30 + 
     $     fChi0nusnuLL(4)*gmneut40

!-------------------------------------------------------------

      piesnu = 0.d0
      piesnu = (1.d0/(16.d0*(pi)**2.d0)) * 
     $     (sbottom + higgs + chargino + neutralino)


      elsnu = dsqrt((snu(3)*snu(3)) - piesnu)	  

      RETURN

      END SUBROUTINE pielsnu
C=========================================================================================
