****f* SuSeFLAV/oneloopselfenergy.f 
*  NAME
*    oneloopselfenergy
*  SYNOPSIS
*    One loop correction to Z,W,hA,Hpm. 
*  FUNCTION
*     Computes self energy for W,Z,hA and Hpm, given, energy scale and external momenta
*  INPUTS
*     p                              - External momentum  
*     q                              - Energy scale
*     g,gp,g3                        - Gauge couplings(g = em, gp = weak, g3 = strong)                                     
*     mt,mb,mtau                     - pole masses of top, bottom and tau
*     mSQRG,mSDRG,mSURG,mSLRG,mSERG  - (2 X 2) mass matrix definition
*     AURG,ADRG,AERG                 - Trilinear couplings
*     mh0sq,mhu0sq,mhpmsq,mA0sq      - physical higgs mass squared 
*     modmu                          - modulus of the \mu paramter 
*     tanbeta                        - the ratio of the vevs of the two Higgs doublet fields.
*     SUegg    =  6 eigenvalues (ascending order) of UP-Squark mass matrix.
*     SDegg    =  6 eigenvalues (ascending order) of Down-Squark mass matrix.
*     Slegg    =  6 eigenvalues (ascending order) of slepton mass matrix.
*     SNegg    =  3 eigenvalues (ascending order) of sneutrino mass matrix.
*     ON       =  (4 X 4) orthogonal matrix such that 
*                ON.MNeut.Transpose[ON] = Diag[MNeut] 
*     Neg      =  4 singular values (descending order) of the Neutralino mass matrix. 
*
*     OCR, OCL =  (2 X 2) orthogonal matrices such that 
*                 MChar = Transpose[OCR].Diag[MChar].OCL
*     Ceg      =   2 singular values of the Neutralino Mass Matrix
*
*  RESULT
*
*   pizzT      - Z boson self energy 
*   piwwT      - W boson self energy
*   piaaT      - Higgs boson A self energy  
*   piHpHm     - Charged higgs H+ self energy  
*
*  EXAMPLE
*
*          SUBROUTINE pizz(p,q,g,mt,mb,mtau,tanbeta,mSQRG,mSDRG,mSURG,mSLRG,
*     $     mSERG,AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,
*     $     mhu0sq,mhpmsq,mA0sq,modmu,ON,OCL,OCR,pizzT)
*       
*           SUBROUTINE piww(p,q,g,mt,mb,mtau,tanbeta,mSQRG,mSDRG,mSURG,mSLRG,
*     $     mSERG,AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,
*     $     mhu0sq,mhpmsq,mA0sq,modmu,ON,OCL,OCR,piwwT)
*
*           SUBROUTINE piaa(p,q,g,gp,mt,mb,mtau,tanbeta,mSQRG,mSDRG,mSURG,
*     $     mSLRG,mSERG,yuRG,ydRG,yeRG,AURG,ADRG,AERG,SUegg,SDegg,SLegg,
*     $     SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,modmu,ON,OCL,OCR,
*     $     piaaT)
*
*          SUBROUTINE pihphm(p,q,g,gp,mt,mb,mtau,tanbeta,mSQRG,mSDRG,mSURG,
*     $     mSLRG,mSERG,AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
*     $     yuRG,yeRG,ydRG,mh0sq,mhu0sq,mhpmsq,mA0sq,modmu,ON,OCL,OCR,
*     $     pihphmT)
*
*  NOTES
*    1. q, the energy scale at which the corrections are added = MZ for W and Z bosons.
*    2. Running values of gauge couplings( Rge output) are used. 
*    3. Pole masses: DRbar scheme is followed.
*    4. Conventions followed are that of BPMZ.
*  BUGS
*    ---
*  SEE ALSO
*    ----
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     1. Inside every subroutine running DRbar mass (MW,MZ etc.) should be used, but not used.
C     2. pizz, piww, piaa, pihphm is checked.
C     3. Some B0, G etc. function's argument leads to diverging value (like B0(mu,mu) in piaa 
C     etc.). This functions needs to be checked.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


C=========================================================================================
C     TRANSVERSE PART OF Z-boson self energy
C
C----------------------------------------------------------------------------------------
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C     pizz is checked @ 07:30 on 02/06/2010.
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


      SUBROUTINE pizz(p,q,g,mt,mb,mtau,tanbeta,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,
     $     mhu0sq,mhpmsq,mA0sq,ON,OCL,OCR,sinsqthw,pizzT)

      IMPLICIT NONE 

      INTEGER i,j
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      double precision mT, mB, mTau
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

      DOUBLE PRECISION up(2), dn(2), el(2)
      DOUBLE PRECISION top(2, 2), btm(2, 2), tau(2, 2)

      DOUBLE PRECISION pizzT,ans,sinsqthw,g
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION b0MwMw,b0mzmh0,b0mzmHu0
      DOUBLE PRECISION alpha,beta
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION thetatz,thetabz,thetatauz

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetacz,thetasz,thetamuz

      DOUBLE PRECISION thetau,thetad,thetae
      DOUBLE PRECISION thetauz,thetadz,thetaez

      DOUBLE PRECISION smHz,susyp,smp,sumHz,sleps

      DOUBLE PRECISION hmumu,hmcmc,hmtopmtop
      DOUBLE PRECISION hmdmd,hmsms,hmbmb,h00
      DOUBLE PRECISION hmeme,hmmummu,hmtaumtau

      DOUBLE PRECISION hmneut1mneut1,hmneut1mneut2
      DOUBLE PRECISION hmneut1mneut3,hmneut1mneut4
      DOUBLE PRECISION hmneut2mneut1,hmneut2mneut2
      DOUBLE PRECISION hmneut2mneut3,hmneut2mneut4
      DOUBLE PRECISION hmneut3mneut1,hmneut3mneut2
      DOUBLE PRECISION hmneut3mneut3,hmneut3mneut4
      DOUBLE PRECISION hmneut4mneut1,hmneut4mneut2
      DOUBLE PRECISION hmneut4mneut3,hmneut4mneut4
      DOUBLE PRECISION hmch1mch1,hmch1mch2,hmch2mch1,hmch2mch2

      DOUBLE PRECISION b0mtopmtop,b0mbmb,b0mtaumtau
      
      DOUBLE PRECISION b0mneut1mneut1,b0mneut1mneut2
      DOUBLE PRECISION b0mneut1mneut3,b0mneut1mneut4      
      DOUBLE PRECISION b0mneut2mneut1,b0mneut2mneut2
      DOUBLE PRECISION b0mneut2mneut3,b0mneut2mneut4      
      DOUBLE PRECISION b0mneut3mneut1,b0mneut3mneut2
      DOUBLE PRECISION b0mneut3mneut3,b0mneut3mneut4
      DOUBLE PRECISION b0mneut4mneut1,b0mneut4mneut2
      DOUBLE PRECISION b0mneut4mneut3,b0mneut4mneut4
      DOUBLE PRECISION b0mch1mch1,b0mch1mch2,b0mch2mch1,b0mch2mch2

      DOUBLE PRECISION b22tmzmh0,b22tMwMw,b22tmA0mHu0,b22tmzmHu0
      DOUBLE PRECISION b22tmA0mh0, b22tmHpmmHpm

      DOUBLE PRECISION b22tmuLmuL,b22tmuRmuR,b22tmcLmcL,b22tmcRmcR
      DOUBLE PRECISION b22tmtLmtL,b22tmtRmtR,b22tmtLmtR,b22tmtRmtL

      DOUBLE PRECISION b22tmdLmdL,b22tmdRmdR,b22tmsLmsL,b22tmsRmsR
      DOUBLE PRECISION b22tmbLmbL,b22tmbRmbR,b22tmbLmbR,b22tmbRmbL

      DOUBLE PRECISION b22tmeLmeL,b22tmeRmeR,b22tmmuLmmuL,b22tmmuRmmuR
      DOUBLE PRECISION b22tmtauLmtauL,b22tmtauRmtauR
      DOUBLE PRECISION b22tmtauLmtauR,b22tmtauRmtauL

      DOUBLE PRECISION b22tsnu1snu1, b22tsnu2snu2,b22tsnu3snu3

      DOUBLE PRECISION fnotz(4,4),gnotz(4,4),fposz(2,2),gposz(2,2)
      DOUBLE PRECISION mchargino(2),charginoterm

      DOUBLE PRECISION apsinotz(4,4), bpsinotz(4,4)
      DOUBLE PRECISION achinotz(4,4), bchinotz(4,4)

      DOUBLE PRECISION apsiposz(2,2), bpsiposz(2,2)
      DOUBLE PRECISION achiposz(2,2), bchiposz(2,2)
      
      DOUBLE PRECISION neutterm,hxnot(4,4),bxnot(4,4)
      DOUBLE PRECISION hxpos(2,2),bxpos(2,2),OCRdag(2,2),OCLdag(2,2)
      DOUBLE PRECISION ON(4,4),ONdag(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION mu,mtop,tan2beta,nmneut(4),sinsqthw_susy
      Integer lopt,rhn
      double precision MZrun,MWrun
      
      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole

!      DOUBLE PRECISION  sinsqthw_mz
c$$$       common/sinsq_mz/sinsqthw_mz
      
           
      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

      common/sfmixing_mz/thetatz,thetabz,thetatauz,thetacz,thetasz,
     $     thetamuz,thetauz,thetadz,thetaez

      common/sinsq_susy/sinsqthw_susy
      common/loops/ lopt,rhn

      common/gbrunning/ MZrun,MWrun


      EXTERNAL b0,b22t,h,mat3prod4d,dag4d,dag2d,mat3prod2d,theta


      include 'stdinputs.h'

!---------------------------------------------------------------------------
     
      pi = 4.d0 * datan(1.d0)

      MW = MWpole
      MZ = MZpole
      
      mu   = mUQ
      mtop = mt
      
      
      sinthw   = dsqrt(sinsqthw)
      thw      = dasin(sinthw)
      costhw   = dcos(thw)
      cos2thw  = dcos(2.d0*thw)
      


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

     
      pizzT = 0.d0

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

      if(q.lt.92.d0)then
         costhetat = dcos(thetatz)
         sinthetat = dsin(thetatz)
         costhetab = dcos(thetabz)
         sinthetab = dsin(thetabz)
         costhetatau = dcos(thetatauz)
         sinthetatau = dsin(thetatauz)

      else

         costhetat = dcos(thetat)
         sinthetat = dsin(thetat)
         costhetab = dcos(thetab)
         sinthetab = dsin(thetab)
         costhetatau = dcos(thetatau)
         sinthetatau = dsin(thetatau)

      endif



      beta = datan(tanbeta)

!---------------------------------------------------------------------------

      up(1) = guL
      up(2) = guR

      top(1, 1) = guL * costhetat**2.d0 - guR * sinthetat**2.d0
      top(1, 2) = (guL + guR) * costhetat * sinthetat
      top(2, 2) = guR * costhetat**2.d0 - guL * sinthetat**2.d0
      top(2, 1) = top(1, 2)
      
      dn(1) = gdL
      dn(2) = gdR
      
      btm(1, 1) = gdL * costhetab**2.d0 - gdR * sinthetab**2.d0
      btm(1, 2) = (gdL + gdR) * costhetab * sinthetab
      btm(2, 2) = gdR * costhetab**2.d0 - gdL * sinthetab**2.d0
      btm(2, 1) = btm(1, 2)
      
      el(1) = geL
      el(2) = geR
      
      tau(1, 1) = geL * costhetatau**2.d0 - geR * sinthetatau**2.d0
      tau(1, 2) = (geL + geR) * costhetatau * sinthetatau
      tau(2, 2) = geR * costhetatau**2.d0 - geL * sinthetatau**2.d0
      tau(2, 1) = tau(1, 2)
    
C----------------------------------------------------------------
     

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      tan2beta = dtan(2.d0*beta)
      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/
     $     (mA0sq - MZ*MZ))*tan2beta)



C-----------------------------------------------------------------------    

      call b0(p,MW,MW,q,b0MwMw)
      call b0(p,mz,mh0,q,b0mzmh0)     
      
      call b0(p,mz,mHu0,q,b0mzmHu0)
      
      
      call h(p,mu,mu,q,hmumu)
      call h(p,mc,mc,q,hmcmc)
      call h(p,mtop,mtop,q,hmtopmtop)
      call b0(p,mtop,mtop,q,b0mtopmtop)
      
      call h(p,md,md,q,hmdmd)
      call h(p,ms,ms,q,hmsms)
      call h(p,mb,mb,q,hmbmb)
      call b0(p,mb,mb,q,b0mbmb)

      call h(p,me,me,q,hmeme)
      call h(p,mmu,mmu,q,hmmummu)
      call h(p,mtau,mtau,q,hmtaumtau)
      call b0(p,mtau,mtau,q,b0mtaumtau)
!-------------------------------------------------------------------

      call b0(p,nmneut(1),nmneut(1),q,b0mneut1mneut1)
      call b0(p,nmneut(1),nmneut(2),q,b0mneut1mneut2)
      call b0(p,nmneut(1),nmneut(3),q,b0mneut1mneut3)
      call b0(p,nmneut(1),nmneut(4),q,b0mneut1mneut4)

      call b0(p,nmneut(2),nmneut(1),q,b0mneut2mneut1)
      call b0(p,nmneut(2),nmneut(2),q,b0mneut2mneut2)
      call b0(p,nmneut(2),nmneut(3),q,b0mneut2mneut3)
      call b0(p,nmneut(2),nmneut(4),q,b0mneut2mneut4)

      call b0(p,nmneut(3),nmneut(1),q,b0mneut3mneut1)
      call b0(p,nmneut(3),nmneut(2),q,b0mneut3mneut2)
      call b0(p,nmneut(3),nmneut(3),q,b0mneut3mneut3)
      call b0(p,nmneut(3),nmneut(4),q,b0mneut3mneut4)

      call b0(p,nmneut(4),nmneut(1),q,b0mneut4mneut1)
      call b0(p,nmneut(4),nmneut(2),q,b0mneut4mneut2)
      call b0(p,nmneut(4),nmneut(3),q,b0mneut4mneut3)
      call b0(p,nmneut(4),nmneut(4),q,b0mneut4mneut4)
                 
      call b0(p,dabs(mchargino(1)),dabs(mchargino(1)),q,b0mch1mch1)
      call b0(p,dabs(mchargino(1)),dabs(mchargino(2)),q,b0mch1mch2)
      call b0(p,dabs(mchargino(2)),dabs(mchargino(1)),q,b0mch2mch1)
      call b0(p,dabs(mchargino(2)),dabs(mchargino(2)),q,b0mch2mch2)
     
      call h(p,nmneut(1),nmneut(1),q,hmneut1mneut1)
      call h(p,nmneut(1),nmneut(2),q,hmneut1mneut2)
      call h(p,nmneut(1),nmneut(3),q,hmneut1mneut3)
      call h(p,nmneut(1),nmneut(4),q,hmneut1mneut4)
 
      call h(p,nmneut(2),nmneut(1),q,hmneut2mneut1)
      call h(p,nmneut(2),nmneut(2),q,hmneut2mneut2)
      call h(p,nmneut(2),nmneut(3),q,hmneut2mneut3)
      call h(p,nmneut(2),nmneut(4),q,hmneut2mneut4)
 
      call h(p,nmneut(3),nmneut(1),q,hmneut3mneut1)
      call h(p,nmneut(3),nmneut(2),q,hmneut3mneut2)
      call h(p,nmneut(3),nmneut(3),q,hmneut3mneut3)
      call h(p,nmneut(3),nmneut(4),q,hmneut3mneut4)
 
      call h(p,nmneut(4),nmneut(1),q,hmneut4mneut1)
      call h(p,nmneut(4),nmneut(2),q,hmneut4mneut2)
      call h(p,nmneut(4),nmneut(3),q,hmneut4mneut3)
      call h(p,nmneut(4),nmneut(4),q,hmneut4mneut4)
      
      call h(p,dabs(mchargino(1)),dabs(mchargino(1)),q,hmch1mch1)
      call h(p,dabs(mchargino(1)),dabs(mchargino(2)),q,hmch1mch2)
      call h(p,dabs(mchargino(2)),dabs(mchargino(1)),q,hmch2mch1)
      call h(p,dabs(mchargino(2)),dabs(mchargino(2)),q,hmch2mch2)
      
      call h(p,2.d-5,2.d-5,q,h00)

     
      call b22t(p,mz,mh0,q,b22tmzmh0)
      call b22t(p,Mw,Mw,q,b22tMwMw)
      call b22t(p,mA0,mHu0,q,b22tmA0mHu0)
      call b22t(p,mz,mHu0,q,b22tmzmHu0)
      call b22t(p,mA0,mh0,q,b22tmA0mh0)
      call b22t(p,mHpm,mHpm,q,b22tmHpmmHpm)

      call b22t(p,muL,muL,q,b22tmuLmuL)
      call b22t(p,muR,muR,q,b22tmuRmuR)
      call b22t(p,mcL,mcL,q,b22tmcLmcL)
      call b22t(p,mcR,mcR,q,b22tmcRmcR)
      call b22t(p,mtL,mtL,q,b22tmtLmtL)
      call b22t(p,mtR,mtR,q,b22tmtRmtR)
      call b22t(p,mtL,mtR,q,b22tmtLmtR)
      call b22t(p,mtR,mtL,q,b22tmtRmtL)

      call b22t(p,mdL,mdL,q,b22tmdLmdL)
      call b22t(p,mdR,mdR,q,b22tmdRmdR)
      call b22t(p,msL,msL,q,b22tmsLmsL)
      call b22t(p,msR,msR,q,b22tmsRmsR)
      call b22t(p,mbL,mbL,q,b22tmbLmbL)
      call b22t(p,mbR,mbR,q,b22tmbRmbR)
      call b22t(p,mbL,mbR,q,b22tmbLmbR)
      call b22t(p,mbR,mbL,q,b22tmbRmbL)
      
      call b22t(p,meL,meL,q,b22tmeLmeL)
      call b22t(p,meR,meR,q,b22tmeRmeR)
      call b22t(p,mmuL,mmuL,q,b22tmmuLmmuL)
      call b22t(p,mmuR,mmuR,q,b22tmmuRmmuR)
      call b22t(p,mtauL,mtauL,q,b22tmtauLmtauL)
      call b22t(p,mtauR,mtauR,q,b22tmtauRmtauR)
      call b22t(p,mtauL,mtauR,q,b22tmtauLmtauR)
      call b22t(p,mtauR,mtauL,q,b22tmtauRmtauL)
      
      

      call b22t(p,snu(1),snu(1),q,b22tsnu1snu1)
      call b22t(p,snu(2),snu(2),q,b22tsnu2snu2)
      call b22t(p,snu(3),snu(3),q,b22tsnu3snu3)

   
C----------------------------------------------------------------

      smHz = - ((dsin(alpha - beta))**2.d0) * ((b22tmzmh0) - !<-------MW, MZ running should be used....
     $     (mz*mz) * b0mzmh0 + b22tmA0mHu0) - 
     $     2.d0 * (costhw**4.d0) * (2.d0 * (p*p) + MW*MW -
     $     (((mz*mz)*(sinthw)**4.d0)/costhw**2.d0)) * b0MwMw -
     $     (8.d0 * (costhw**4.d0) + (cos2thw)**2.d0) * b22tMwMw


      sumHz = - ((dcos(alpha - beta))**2.d0) * ((b22tmzmHu0) +
     $     (b22tmA0mh0) - (mz*mz) * b0mzmHu0) -
     $     ((cos2thw)**2.d0) * b22tmHpmmHpm

C-------------------------------------------------------------------------------
      
      susyp = - 12.d0 * ((up(1)**2.d0) * b22tmuLmuL + !<------- up squarks, first generation
     $     (up(2)**2.d0) * b22tmuRmuR) -
     $     12.d0 * ((up(1)**2.d0) * b22tmcLmcL + !<------- second generation        
     $     (up(2)**2.d0) * b22tmcRmcR) -   
     $     12.d0 * ((top(1,1)**2.d0) * b22tmtLmtL + !<------- third generation 
     $     (top(1,2)**2.d0) * b22tmtLmtR +
     $     (top(2,1)**2.d0) * b22tmtRmtL +
     $     (top(2,2)**2.d0) * b22tmtRmtR) -
     $     12.d0 * ((dn(1)**2.d0) * b22tmdLmdL + !<------- down squarks, first generation
     $     (dn(2)**2.d0) * b22tmdRmdR) -
     $     12.d0 * ((dn(1)**2.d0) * b22tmsLmsL + !<------- second generation
     $     (dn(2)**2.d0) * b22tmsRmsR) -                       
     $     12.d0 * ((btm(1,1)**2.d0) * b22tmbLmbL + !<------- third generation 
     $     (btm(1,2)**2.d0) * b22tmbLmbR +
     $     (btm(2,1)**2.d0) * b22tmbRmbL +
     $     (btm(2,2)**2.d0) * b22tmbRmbR) 

      sleps = - 4.d0 * ((el(1)**2.d0) * b22tmeLmeL + !<------- sleptons, first generation
     $     (el(2)**2.d0) * b22tmeRmeR) -
     $     4.d0 * ((el(1)**2.d0) * b22tmmuLmmuL + !<------- second generation
     $     (el(2)**2.d0) * b22tmmuRmmuR) -                       
     $     4.d0 * ((tau(1,1)**2.d0) * b22tmtauLmtauL + !<------- third generation 
     $     (tau(1,2)**2.d0) * b22tmtauLmtauR  +
     $     (tau(2,1)**2.d0) * b22tmtauRmtauL +
     $     (tau(2,2)**2.d0) * b22tmtauRmtauR) -
     $     (b22tsnu1snu1 + b22tsnu2snu2 + b22tsnu3snu3) !<-------- sneutrinos
      
C--------------------------------------------------------------------------

      smp = 3.d0 * (hmumu * ((guL)**2.d0 + (guR)**2.d0) +      !<--------- up sector
     $     hmcmc * ((guL)**2.d0 + (guR)**2.d0) +
     $     hmtopmtop * ((guL)**2.d0 + (guR)**2.d0) - 
     $     4.d0 * guL*guR * (mtop**2.d0) * b0mtopmtop) +
     $     3.d0 * (hmdmd * ((gdL)**2.d0 + (gdR)**2.d0) +        !<--------- down sector
     $     hmsms * ((gdL)**2.d0 + (gdR)**2.d0) +
     $     hmbmb * ((gdL)**2.d0 + (gdR)**2.d0) - 
     $     4.d0 * gdL*gdR * (mb**2.d0) * b0mbmb) +
     $     1.d0 * (hmeme * ((geL)**2.d0 + (geR)**2.d0) +        !<--------- electron sector
     $     hmmummu * ((geL)**2.d0 + (geR)**2.d0) + 
     $     hmtaumtau * ((geL)**2.d0 + (geR)**2.d0) - 
     $     4.d0 * geL*geR * (mtau**2.d0) * b0mtaumtau) +
     $     3.d0 * h00 * 0.25d0                                  !<-------neutrino sector

      
C---------------------------------------------------------------------------
C     Neutralinos
C---------------------------------------------------------------------------

      loopsii: DO i = 1, 4
      loopsij: DO j = 1, 4

      apsinotz(i,j) = 0.d0
      bpsinotz(i,j) = 0.d0
      
      achinotz(i,j) = 0.d0
      bchinotz(i,j) = 0.d0

      ENDDO loopsij
      ENDDO loopsii


      apsinotz(3,3) =  g/(2.d0*costhw)
      apsinotz(4,4) = -1.d0*apsinotz(3,3)
      bpsinotz(3,3) = -1.d0*apsinotz(3,3)
      bpsinotz(4,4) = -1.d0*apsinotz(4,4)

      

      call dag4d(ON,ONdag)

      call mat3prod4d(ON,apsinotz,ONdag,achinotz)
      call mat3prod4d(ON,bpsinotz,ONdag,bchinotz)  

      loopnoti: DO i = 1, 4 
      loopnotj: DO j = 1, 4

      fnotz(i,j) = achinotz(i,j)**2.d0 + bchinotz(i,j)**2.d0
      gnotz(i,j) = 2.d0 * (bchinotz(i,j)*achinotz(i,j))
      
      
      ENDDO loopnotj
      ENDDO loopnoti

      neutterm = 0.d0


      hxnot(1,1) = hmneut1mneut1
      hxnot(1,2) = hmneut1mneut2
      hxnot(1,3) = hmneut1mneut3
      hxnot(1,4) = hmneut1mneut4
      hxnot(2,1) = hmneut2mneut1
      hxnot(2,2) = hmneut2mneut2
      hxnot(2,3) = hmneut2mneut3
      hxnot(2,4) = hmneut2mneut4
      hxnot(3,1) = hmneut3mneut1
      hxnot(3,2) = hmneut3mneut2
      hxnot(3,3) = hmneut3mneut3
      hxnot(3,4) = hmneut3mneut4
      hxnot(4,1) = hmneut4mneut1
      hxnot(4,2) = hmneut4mneut2
      hxnot(4,3) = hmneut4mneut3
      hxnot(4,4) = hmneut4mneut4
      
      
      bxnot(1,1) = mneut(1)*mneut(1)*b0mneut1mneut1
      bxnot(1,2) = mneut(2)*mneut(1)*b0mneut1mneut2
      bxnot(1,3) = mneut(3)*mneut(1)*b0mneut1mneut3
      bxnot(1,4) = mneut(4)*mneut(1)*b0mneut1mneut4
      bxnot(2,1) = mneut(2)*mneut(1)*b0mneut2mneut1
      bxnot(2,2) = mneut(2)*mneut(2)*b0mneut2mneut2
      bxnot(2,3) = mneut(2)*mneut(3)*b0mneut2mneut3
      bxnot(2,4) = mneut(2)*mneut(4)*b0mneut2mneut4
      bxnot(3,1) = mneut(3)*mneut(1)*b0mneut3mneut1
      bxnot(3,2) = mneut(3)*mneut(2)*b0mneut3mneut2
      bxnot(3,3) = mneut(3)*mneut(3)*b0mneut3mneut3
      bxnot(3,4) = mneut(3)*mneut(4)*b0mneut3mneut4
      bxnot(4,1) = mneut(4)*mneut(1)*b0mneut4mneut1
      bxnot(4,2) = mneut(4)*mneut(2)*b0mneut4mneut2
      bxnot(4,3) = mneut(4)*mneut(3)*b0mneut4mneut3
      bxnot(4,4) = mneut(4)*mneut(4)*b0mneut4mneut4


      loopnlinoi: DO i = 1, 4
      loopnlinoj: DO j = 1, 4


      neutterm = neutterm + 0.5d0 * ((costhw/g)**2.d0) * (fnotz(i,j) *
     $     hxnot(i,j) + 2.d0 * gnotz(i,j) * bxnot(i,j))
        
      ENDDO loopnlinoj
      ENDDO loopnlinoi 


C------------------------------------------------------------------------------
C     Chargino Terms
C------------------------------------------------------------------------------

      apsiposz(1,1) = g*costhw
      apsiposz(1,2) = 0.d0
      apsiposz(2,1) = 0.d0
      apsiposz(2,2) = (g*cos2thw)/(2.d0 * costhw)
      bpsiposz(1,1) = apsiposz(1,1)
      bpsiposz(1,2) = 0.d0
      bpsiposz(2,1) = 0.d0
      bpsiposz(2,2) = apsiposz(2,2)


      call dag2d(OCL,OCLdag)
      call dag2d(OCR,OCRdag)
      call mat3prod2d(OCR,apsiposz,OCRdag,achiposz)    !<--------------------CHECK!!
      call mat3prod2d(OCL,bpsiposz,OCLdag,bchiposz)
      
      loopposi: DO i = 1, 2 
      loopposj: DO j = 1, 2

      fposz(i,j) = achiposz(i,j)**2.d0 + bchiposz(i,j)**2.d0
      gposz(i,j) = 2.d0 * (bchiposz(i,j)*achiposz(i,j))
      
      ENDDO loopposj
      ENDDO loopposi

      charginoterm = 0.d0

      hxpos(1,1) = hmch1mch1
      hxpos(1,2) = hmch1mch2
      hxpos(2,1) = hmch2mch1
      hxpos(2,2) = hmch2mch2

      bxpos(1,1) = mchargino(1)*mchargino(1)*b0mch1mch1
      bxpos(1,2) = mchargino(1)*mchargino(2)*b0mch1mch2
      bxpos(2,1) = mchargino(2)*mchargino(1)*b0mch2mch1
      bxpos(2,2) = mchargino(2)*mchargino(2)*b0mch2mch2

      
      loopchginoi: DO i = 1, 2
      loopchginoj: DO j = 1, 2


      charginoterm = charginoterm + (costhw/g)**2.d0 * (fposz(i,j) *
     $     hxpos(i,j) + 2.d0 * gposz(i,j) * bxpos(i,j))

        
      ENDDO loopchginoj
      ENDDO loopchginoi

      
C--------------------------------------------------------------------------------
c$$$      if(rhn.eq.1)then
c$$$         sleps = 0.d0
c$$$      endif

      ans = smp + charginoterm + neutterm + sumHz + smHz + susyp 
     $     + sleps

      pizzT = ans *
     $     (g**2.d0)/((costhw**2.d0) * 16.d0 *(pi**2.d0))


      RETURN

      END SUBROUTINE pizz

C=========================================================================================
C     TRANSVERSE PART OF W-boson self energy
C
C----------------------------------------------------------------------------------------
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C     piww is checked @ 09:00 on 02/06/2010.
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      SUBROUTINE piww(p,q,g,mt,mb,mtau,tanbeta,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,mh0sq,
     $     mhu0sq,mhpmsq,mA0sq,ON,OCL,OCR,sinsqthw,piwwT)

      IMPLICIT NONE 

      INTEGER i,j,k,lopt,rhn
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2)
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      double precision mT, mB, mTau
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

c$$$      DOUBLE PRECISION up(2), dn(2), el(2)
c$$$      DOUBLE PRECISION top(2, 2), btm(2, 2), tau(2, 2)

      DOUBLE PRECISION piwwT,ans,sinsqthw,g
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION b0mhMw,b0MzMw,b0Mw0,b0mHuMw,b0mAmHpm
      DOUBLE PRECISION alpha,beta
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION thetatz,thetabz,thetatauz

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetacz,thetasz,thetamuz

c$$$      DOUBLE PRECISION costhetac,sinthetac,costhetas,sinthetas
c$$$      DOUBLE PRECISION costhetamu,sinthetamu

      DOUBLE PRECISION thetau,thetad,thetae
      DOUBLE PRECISION thetauz,thetadz,thetaez

c$$$      DOUBLE PRECISION costhetau,sinthetau,costhetad,sinthetad
c$$$      DOUBLE PRECISION sinthetae

      DOUBLE PRECISION smhz,susyp,smp,sumhz,sleps

      DOUBLE PRECISION hmumd,hmcms,hmtopmb,h0me,h0mmu,h0mtau
      DOUBLE PRECISION hmneut1ch1,hmneut1ch2,hmneut2ch1,hmneut2ch2
      DOUBLE PRECISION hmneut3ch1,hmneut3ch2,hmneut4ch1,hmneut4ch2

      DOUBLE PRECISION b0mneut1ch1,b0mneut1ch2,b0mneut2ch1,b0mneut2ch2
      DOUBLE PRECISION b0mneut3ch1,b0mneut3ch2,b0mneut4ch1,b0mneut4ch2
     
      DOUBLE PRECISION b22tmhMw,b22tmHumHpm,b22tMzMw,b22tMw0
      DOUBLE PRECISION b22tmhmHpm,b22tmHuMw,b22tmAmHpm     

      DOUBLE PRECISION b22tmuLmdL,b22tmcLmsL,b22tmtLmbL,b22tmtLmbR
      DOUBLE PRECISION b22tmtRmbL,b22tmtRmbR,b22tsnu1meL,b22tsnu2mmuL
      DOUBLE PRECISION b22tsnu3mtauL,b22tsnu3mtauR

      DOUBLE PRECISION fijw(4,2),gijw(4,2)
      DOUBLE PRECISION mchargino(2),achi(4,2),bchi(4,2)

      DOUBLE PRECISION apsinotpw(4,2), bpsinotpw(4,2)
      DOUBLE PRECISION achinotpw(4,2), bchinotpw(4,2)
    
      DOUBLE PRECISION npchino,hxnotp(4,2),bxnotp(4,2)
      DOUBLE PRECISION OCRdag(2,2),OCLdag(2,2)
      DOUBLE PRECISION ON(4,4),ONdag(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION mu,mtop,tan2beta,w(2,2),nmneut(4)

      double precision testchk

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole


      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

      common/sfmixing_mz/thetatz,thetabz,thetatauz,thetacz,thetasz,
     $     thetamuz,thetauz,thetadz,thetaez
            
      common/loops/ lopt,rhn


      EXTERNAL b0,b22t,h,dag2d


      include 'stdinputs.h'

!---------------------------------------------------------------------------

      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      MZ = MZpole

      mu = mUQ
      mtop = mt

      
      sinthw   = dsqrt((sinsqthw))
      thw      = dasin(sinthw)
      costhw   = dcos(thw)
      cos2thw  = dcos(2.d0*thw)

      beta     = datan(tanbeta)
      

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

      piwwT = 0.d0

!----------------------------------------------------------------------------

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


!----------------------------------------------------------------------------

      if(q.lt.93.d0)then

         costhetat = dcos(thetatz)
         sinthetat = dsin(thetatz)
         costhetab = dcos(thetabz)
         sinthetab = dsin(thetabz)
         costhetatau = dcos(thetatauz)
         sinthetatau = dsin(thetatauz)


      else

         costhetat = dcos(thetat)
         sinthetat = dsin(thetat)
         costhetab = dcos(thetab)
         sinthetab = dsin(thetab)
         costhetatau = dcos(thetatau)
         sinthetatau = dsin(thetatau)

      endif

C--------------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))


      tan2beta = dtan(2.d0 * beta)
      alpha    = 0.5d0*datan(((mA0sq + MZ*MZ)/
     $     (mA0sq - MZ*MZ))*tan2beta)


C---------------------------------------------------------------------    
      call b0(p,mHu0,MW,q,b0mhuMw)
      call b0(p,mh0,MW,q,b0mhMw)
      call b0(p,MZ,MW,q,b0MzMw)
      call b0(p,MW,2.d-5,q,b0Mw0)
      call b0(p,mA0,mHpm,q,b0mAmHpm)
      
      call h(p,mu,md,q,hmumd)
      call h(p,mc,ms,q,hmcms)
      call h(p,mt,mb,q,hmtopmb)
      call h(p,2.d-5,me,q,h0me)
      call h(p,2.d-5,mmu,q,h0mmu)
      call h(p,2.d-5,mtau,q,h0mtau)


      call h(p,nmneut(1),dabs(mchargino(1)),q,hmneut1ch1)
      call h(p,nmneut(1),dabs(mchargino(2)),q,hmneut1ch2)
      call h(p,nmneut(2),dabs(mchargino(1)),q,hmneut2ch1)
      call h(p,nmneut(2),dabs(mchargino(2)),q,hmneut2ch2)
      call h(p,nmneut(3),dabs(mchargino(1)),q,hmneut3ch1)
      call h(p,nmneut(3),dabs(mchargino(2)),q,hmneut3ch2)
      call h(p,nmneut(4),dabs(mchargino(1)),q,hmneut4ch1)
      call h(p,nmneut(4),dabs(mchargino(2)),q,hmneut4ch2)

      call b0(p,nmneut(1),dabs(mchargino(1)),q,b0mneut1ch1)
      call b0(p,nmneut(1),dabs(mchargino(2)),q,b0mneut1ch2)
      call b0(p,nmneut(2),dabs(mchargino(1)),q,b0mneut2ch1)
      call b0(p,nmneut(2),dabs(mchargino(2)),q,b0mneut2ch2)
      call b0(p,nmneut(3),dabs(mchargino(1)),q,b0mneut3ch1)
      call b0(p,nmneut(3),dabs(mchargino(2)),q,b0mneut3ch2)
      call b0(p,nmneut(4),dabs(mchargino(1)),q,b0mneut4ch1)
      call b0(p,nmneut(4),dabs(mchargino(2)),q,b0mneut4ch2)

!-------------------------------------------------------------------
      
      call b22t(p,mh0,MW,q,b22tmhMw)
      call b22t(p,mHu0,mHpm,q,b22tmHumHpm)
      call b22t(p,MZ,MW,q,b22tMzMw)
      call b22t(p,MW,2.d-5,q,b22tMw0)
      call b22t(p,mh0,mHpm,q,b22tmhmHpm)
      call b22t(p,mHu0,MW,q,b22tmHuMw)
      call b22t(p,mA0,mHpm,q,b22tmAmHpm)


      call b22t(p,mtL,mbL,q,b22tmtLmbL)
      call b22t(p,mtL,mbR,q,b22tmtLmbR)
      call b22t(p,mtR,mbL,q,b22tmtRmbL)
      call b22t(p,mtR,mbR,q,b22tmtRmbR)

!      call b22t(80.386d0,5002.6d0,4992.5d0,91.188d0,testchk)
            

      call b22t(p,muL,mdL,q,b22tmuLmdL)
      call b22t(p,mcL,msL,q,b22tmcLmsL)
      call b22t(p,snu(3),meL,q,b22tsnu1meL)     
      call b22t(p,snu(2),mmuL,q,b22tsnu2mmuL)

c$$$      print*,"muL, mdL = ", muL, mdL
c$$$      print*,"mcL, msL = ", mcL, msL

c$$$      call b22t(p,muR,mdR,q,b22tmuLmdL) !<=== changed from fL --> fR, one should change the variable name as well
c$$$      call b22t(p,mcR,msR,q,b22tmcLmsL)
c$$$      call b22t(p,snu(3),meR,q,b22tsnu1meL)     
c$$$      call b22t(p,snu(2),mmuR,q,b22tsnu2mmuL)

c$$$      print*,"snu(1), snu(2), snu(3) = ", snu(1), snu(2), snu(3)
c$$$      print*,"meL, meR, mmuL, mmuR = ", meL, meR, mmuL, mmuR
c$$$      print*,"mtauL, mtauR = ", mtauL, mtauR
      
      call b22t(p,snu(1),mtauL,q,b22tsnu3mtauL)
      call b22t(p,snu(1),mtauR,q,b22tsnu3mtauR)
C-------------------------------------------------------------------------------
      
      smHz = - ((dsin(alpha - beta))**2.d0 * ((b22tmhMw) - !<-------MW, MZ running mass should be used.
     $     (MW*MW) * b0mhMw + b22tmHumHpm)) - 
     $     ((costhw**2.d0) * (4.d0 * (p*p) + MW*MW + (MZ*MZ)) - 
     $     (MZ*MZ * (sinthw**4.d0))) * b0MzMw -
     $     ((8.d0 * (costhw**2.d0) + 1.d0) * b22tMzMw) - 
     $     (sinsqthw * (8.d0 * b22tMw0 + 4.d0 * p*p * b0Mw0))
      
      sumHz = - (((dcos(alpha - beta))**2.d0) * ((b22tmhmHpm) +
     $     (b22tmHuMw) - (MW*MW) * b0mHuMw)) - b22tmAmHpm
C------------------------------------------------------------------------------- !<----------- checked
      
      w(1,1) = costhetat*costhetab
      w(1,2) = costhetat*sinthetab
      w(2,1) = sinthetat*costhetab 
      w(2,2) = sinthetat*sinthetab
     
c$$$      print*,"w(1,1) = ", w(1,1)
c$$$      print*,"w(1,2) = ", w(1,2)
c$$$      print*,"w(2,1) = ", w(2,1)
c$$$      print*,"w(2,2) = ", w(2,2)

      susyp = 0.d0

c$$$      susyp = - 6.d0 * (b22tmuLmdL + b22tmcLmsL) - !<---first 2 doublets
c$$$     $     2.d0 * (b22tsnu1meL + b22tsnu2mmuL) -
c$$$     $     6.d0 * ((w(2,2)**2.d0) * b22tmtLmbL + !<----3rd doublet
c$$$     $     (w(2,1)**2.d0) * b22tmtLmbR +
c$$$     $     (w(1,2)**2.d0) * b22tmtRmbL +
c$$$     $     (w(1,1)**2.d0) * b22tmtRmbR)

      susyp = - 6.d0 * (b22tmuLmdL + b22tmcLmsL) - !<---first 2 doublets
     $     2.d0 * (b22tsnu1meL + b22tsnu2mmuL) -
     $     6.d0 * ((w(1,1)**2.d0) * b22tmtLmbL + !<----3rd doublet
     $     (w(1,2)**2.d0) * b22tmtLmbR +
     $     (w(2,1)**2.d0) * b22tmtRmbL +
     $     (w(2,2)**2.d0) * b22tmtRmbR)

c$$$      print*, "b22tmcLmsL = ", b22tmcLmsL
c$$$      print*, "b22tmuLmdL = ", b22tmuLmdL
c$$$      print*, "b22tmsnu1LmeL = ", b22tsnu1meL
c$$$      print*, "b22tmsnu2LmmuL = ", b22tsnu2mmuL
     
c$$$      sleps = -2.d0 * (((costhetatau)**2.d0) * b22tsnu3mtauR +
c$$$     $     ((sinthetatau)**2.d0) * b22tsnu3mtauL)

      sleps = -2.d0 * (((costhetatau)**2.d0) * b22tsnu3mtauL +
     $     ((sinthetatau)**2.d0) * b22tsnu3mtauR)


C---------------------------------------------------------------------------------
    
      smp = 1.5d0 * (hmumd + hmcms + hmtopmb) + 
     $     0.5d0 * (h0me + h0mmu + h0mtau)     

C---------------------------------------------------------------------------
C     Neutralinos and Charginos
C---------------------------------------------------------------------------

      loopsii: DO i = 1, 4
      loopsij: DO j = 1, 2

      apsinotpw(i,j) = 0.d0
      bpsinotpw(i,j) = 0.d0
      
      achinotpw(i,j) = 0.d0
      bchinotpw(i,j) = 0.d0

      achi(i,j) = 0.d0
      bchi(i,j) = 0.d0

      ENDDO loopsij
      ENDDO loopsii


      apsinotpw(2,1) = - 1.d0 * g
      apsinotpw(4,2) =  g/(dsqrt(2.d0))
      bpsinotpw(2,1) = - 1.d0 * g
      bpsinotpw(3,2) = - 1.d0 * (g/(dsqrt(2.d0)))
 

      call dag2d(ON,ONdag)
      call dag2d(OCL,OCLdag)
      call dag2d(OCR,OCRdag)

!--------------------------------------------------------------

      loop4by2i: DO i = 1,4
      loop4by2j: DO j = 1,2

      achi(i,j) = 0.d0
      bchi(i,j) = 0.d0

      loop4by2k: DO k = 1,4
      
      
      achi(i,j) = achi(i,j) + ON(i,k) * apsinotpw(k,j)
      bchi(i,j) = bchi(i,j) + ON(i,k) * bpsinotpw(k,j)
      
      
      ENDDO loop4by2k
      ENDDO loop4by2j
      ENDDO loop4by2i
!---------------------------------------------------------------
      loop2i: DO i = 1, 4
      loop2j: DO j = 1, 2

      achinotpw(i,j) = 0.d0
      bchinotpw(i,j) = 0.d0

      loop2k: DO k = 1,2

c$$$      achinotpw(i,j) = achinotpw(i,j) + achi(i,k) * OCR(k,j) <--- 1.1.2 
c$$$      bchinotpw(i,j) = bchinotpw(i,j) + bchi(i,k) * OCL(k,j)

      achinotpw(i,j) = achinotpw(i,j) + achi(i,k) * OCRdag(k,j)
      bchinotpw(i,j) = bchinotpw(i,j) + bchi(i,k) * OCLdag(k,j)
      
      ENDDO loop2k
      ENDDO loop2j
      ENDDO loop2i

!----------------------------------------------------------------

      loopnoti: DO i = 1, 4 
      loopnotj: DO j = 1, 2

      fijw(i,j) = (achinotpw(i,j)**2.d0) + (bchinotpw(i,j)**2.d0)
      gijw(i,j) = 2.d0 * real(bchinotpw(i,j) * achinotpw(i,j))
      
      
      ENDDO loopnotj
      ENDDO loopnoti

!--------------------------------------------------------------------------

      hxnotp(1,1) = hmneut1ch1
      hxnotp(1,2) = hmneut1ch2
      hxnotp(2,1) = hmneut2ch1
      hxnotp(2,2) = hmneut2ch2
      hxnotp(3,1) = hmneut3ch1
      hxnotp(3,2) = hmneut3ch2
      hxnotp(4,1) = hmneut4ch1
      hxnotp(4,2) = hmneut4ch2

      

      bxnotp(1,1) = mneut(1)*mchargino(1)*b0mneut1ch1
      bxnotp(1,2) = mneut(1)*mchargino(2)*b0mneut1ch2
      bxnotp(2,1) = mneut(2)*mchargino(1)*b0mneut2ch1
      bxnotp(2,2) = mneut(2)*mchargino(2)*b0mneut2ch2
      bxnotp(3,1) = mneut(3)*mchargino(1)*b0mneut3ch1
      bxnotp(3,2) = mneut(3)*mchargino(2)*b0mneut3ch2
      bxnotp(4,1) = mneut(4)*mchargino(1)*b0mneut4ch1
      bxnotp(4,2) = mneut(4)*mchargino(2)*b0mneut4ch2

!      print*,"bxnotp(1,2) = ", mneut(1),mchargino(2),b0mneut1ch2

  
      npchino = 0.d0

      loopnlinoi: DO i = 1, 4
      loopnlinoj: DO j = 1, 2


      npchino = npchino + (1.d0/(g*g)) * (fijw(i,j) *
     $     hxnotp(i,j) + 2.d0 * gijw(i,j) * bxnotp(i,j))
      
      ENDDO loopnlinoj
      ENDDO loopnlinoi 
     
c$$$      print*,"npchino terms = ", g,mchargino(1),mchargino(2)
c$$$ 10   format(/A/)
c$$$      print 10, "fijw : "
c$$$      opu22a: do i = 1, 4
c$$$ 223  format(1x,2(2x,1pe11.4))
c$$$      write(*,223) (fijw(i,j),j = 1, 2)
c$$$      enddo opu22a  
c$$$      print 10, "hxnotp : "
c$$$      opu22b: do i = 1, 4
c$$$      write(*,223) (hxnotp(i,j),j = 1, 2)
c$$$      enddo opu22b  
c$$$      print 10, "gijw : "
c$$$      opu22c: do i = 1, 4
c$$$      write(*,223) (gijw(i,j),j = 1, 2)
c$$$      enddo opu22c  
c$$$      print 10, "bxnotp : "
c$$$      opu22d: do i = 1, 4
c$$$      write(*,223) (bxnotp(i,j),j = 1, 2)
c$$$      enddo opu22d  
      

C----------------------------------------------------------------------

!      print*,"piwwt terms = ", sumHz, smHz,npchino,susyp,smp,sleps


c$$$      if(rhn.eq.1)then
c$$$         sleps = 0.d0
c$$$      endif

      ans =  sumHz + smHz  + npchino + susyp + smp + sleps

      piwwT =  ans * ((g**2.d0)/(16.d0 * (pi**2.d0))) 



      RETURN

      END SUBROUTINE piww

C=============================================================================================

C=============================================================================================
C     Higgs boson mA self energy
C
C---------------------------------------------------------------------------------------------
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C     piaa is checked @ 19:30 on 20/06/2010.
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      SUBROUTINE piaa(p,q,g,gp,mt,mb,mtau,tanbeta,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu,modmu,ON,
     $     OCL,OCR,piaaT)

      IMPLICIT NONE 

      INTEGER i,j
      DOUBLE PRECISION tanbeta
      DOUBLE PRECISION AURG(3,3),ADRG(3,3)
      DOUBLE PRECISION AERG(3,3),SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,gp,sgnmu
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR,yeRG(3,3)
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR,ydRG(3,3)
      DOUBLE PRECISION cosbeta,cos2beta,yuRG(3,3),sinbeta


      double precision mT, mB, mTau
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3)

      DOUBLE PRECISION up(2), dn(2), el(2)
      DOUBLE PRECISION top(2, 2), btm(2, 2), tau(2, 2)

!      DOUBLE PRECISION mu,mc,mtop,md,ms,mb,me,mmu,mtau

      DOUBLE PRECISION piaaT,sinsqthw,g
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION alpha,beta
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION thetatz,thetabz,thetatauz

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau

      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetacz,thetasz,thetamuz

c$$$      DOUBLE PRECISION costhetac,sinthetac,costhetas,sinthetas
c$$$      DOUBLE PRECISION costhetamu,sinthetamu

      DOUBLE PRECISION thetau,thetad,thetae
      DOUBLE PRECISION thetauz,thetadz,thetaez

c$$$      DOUBLE PRECISION costhetau,sinthetau,costhetad,sinthetad
c$$$      DOUBLE PRECISION costhetae,sinthetae

      DOUBLE PRECISION smhz,susyp,smp,susyp1,neutterm,sumHz

      DOUBLE PRECISION a0mt,a0mc,a0mu,b0mtmt,b0mcmc,b0mumu
      DOUBLE PRECISION a0muL,a0muR,a0mcL,a0mcR,a0mtL,a0mtR

      DOUBLE PRECISION b0muLmuR,b0muRmuL,b0mcLmcR,b0mcRmcL
      DOUBLE PRECISION b0mtLmtR,b0mtRmtL

      DOUBLE PRECISION b0mdLmdR,b0mdRmdL,b0msLmsR,b0msRmsL
      DOUBLE PRECISION b0mbLmbR,b0mbRmbL

      DOUBLE PRECISION fmhpmMw,fmhuMz,fmhMz,a0Mw,a0Mz

      
      DOUBLE PRECISION gmneut1mneut1,gmneut1mneut2
      DOUBLE PRECISION gmneut1mneut3,gmneut1mneut4
      DOUBLE PRECISION gmneut2mneut1,gmneut2mneut2
      DOUBLE PRECISION gmneut2mneut3,gmneut2mneut4
      DOUBLE PRECISION gmneut3mneut1,gmneut3mneut2
      DOUBLE PRECISION gmneut3mneut3,gmneut3mneut4
      DOUBLE PRECISION gmneut4mneut1,gmneut4mneut2
      DOUBLE PRECISION gmneut4mneut3,gmneut4mneut4
      DOUBLE PRECISION gmch1mch1,gmch1mch2,gmch2mch1,gmch2mch2

      DOUBLE PRECISION a0mb,a0ms,a0md,a0me,a0mmu,a0mtau,a0snu3
      DOUBLE PRECISION b0mbmb,b0mtaumtau
      DOUBLE PRECISION b0mmummu,b0msms,b0meme,b0mdmd

      DOUBLE PRECISION b0mtauLmtauR,b0mtauRmtauL,b0meLmeR,b0meRmeL
      DOUBLE PRECISION b0mmuLmmuR,b0mmuRmmuL,sin2beta
      
      DOUBLE PRECISION a0mdL,a0mdR,a0msL,a0msR,a0mtauL,a0mtauR,a0snu2
      DOUBLE PRECISION a0meL,a0meR,a0mmuL,a0mmuR,a0mbL,a0mbR,a0snu1

      DOUBLE PRECISION b0mneut1mneut1,b0mneut1mneut2
      DOUBLE PRECISION b0mneut1mneut3,b0mneut1mneut4      
      DOUBLE PRECISION b0mneut2mneut1,b0mneut2mneut2
      DOUBLE PRECISION b0mneut2mneut3,b0mneut2mneut4      
      DOUBLE PRECISION b0mneut3mneut1,b0mneut3mneut2
      DOUBLE PRECISION b0mneut3mneut3,b0mneut3mneut4
      DOUBLE PRECISION b0mneut4mneut1,b0mneut4mneut2
      DOUBLE PRECISION b0mneut4mneut3,b0mneut4mneut4
      DOUBLE PRECISION b0mch1mch1,b0mch1mch2,b0mch2mch1,b0mch2mch2

      DOUBLE PRECISION a0mA,a0mh,a0mhpm,a0mhu
      DOUBLE PRECISION b0mAmHu,b0mAmh,b0MWmHu
      DOUBLE PRECISION b0Mzmh,b0MzmHu


      DOUBLE PRECISION fnotA(4,4),gnotA(4,4)
      DOUBLE PRECISION a2chiA(4,4),b2chiA(4,4),achinotp1(4,4)
      DOUBLE PRECISION achinotp2(4,4),bchinotp1(4,4),bchinotp2(4,4)
      
      DOUBLE PRECISION mchargino(2)

      DOUBLE PRECISION apsinotp1(4,4), bpsinotp1(4,4)
      DOUBLE PRECISION apsinotp2(4,4), bpsinotp2(4,4)
      DOUBLE PRECISION achip1(4,4), bchip1(4,4),achip2(4,4), bchip2(4,4)
      DOUBLE PRECISION apsiposp1(2,2), bpsiposp1(2,2)
      DOUBLE PRECISION apsiposp2(2,2), bpsiposp2(2,2)
      DOUBLE PRECISION achiposA(2,2),bchiposA(2,2)
      DOUBLE PRECISION fposA(2,2),gposA(2,2),hxpos(2,2),bxpos(2,2)
      DOUBLE PRECISION a2chiposA(2,2),b2chiposA(2,2),charginoterm
     
      DOUBLE PRECISION hxnot(4,4),bxnot(4,4)
      DOUBLE PRECISION OCRdag(2,2),OCLdag(2,2)
      DOUBLE PRECISION ON(4,4),ONdag(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION mu,mtop,tan2beta

      DOUBLE PRECISION  lAas(2),lAah(2)
      DOUBLE PRECISION  lAgs(2), lAgh(2)
      DOUBLE PRECISION  lAa12(2,2), lAahh(2,2)
      DOUBLE PRECISION lAaaa,lAagg,malpha,nmneut(4)

      DOUBLE PRECISION ralpha(2,2),mralpha(2,2)

      DOUBLE PRECISION sinsqthw_susy

      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole
      
      
      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

      common/sfmixing_mz/thetatz,thetabz,thetatauz,thetacz,thetasz,
     $     thetamuz,thetauz,thetadz,thetaez

            
      EXTERNAL b0,a0,h,funcg,f,dag2d,theta,mat3prod4d,rmat2d
      EXTERNAL scmul4d,add4d,mat3prod2d


      include 'stdinputs.h'

!---------------------------------------------------------------------------


      pi = 4.d0 * datan(1.d0)
      
      MW = MWpole
      MZ = MZpole

      mu = mUQ
      mtop = mt

      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw = dsqrt(dabs(sinsqthw))
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

      piaaT = 0.d0

!----------------------------------------------------------------------------

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


      if(q.lt.94.d0)then

         costhetat = dcos(thetatz)
         sinthetat = dsin(thetatz)
         costhetab = dcos(thetabz)
         sinthetab = dsin(thetabz)
         costhetatau = dcos(thetatauz)
         sinthetatau = dsin(thetatauz)

      else

         costhetat = dcos(thetat)
         sinthetat = dsin(thetat)
         costhetab = dcos(thetab)
         sinthetab = dsin(thetab)
         costhetatau = dcos(thetatau)
         sinthetatau = dsin(thetatau)

      endif

      up(1) = guL
      up(2) = guR

      top(1, 1) = guL*costhetat**2.d0 - guR*sinthetat**2.d0
      top(1, 2) = (guL + guR)*costhetat*sinthetat
      top(2, 2) = guR*costhetat**2.d0 - guL*sinthetat**2.d0
      top(2, 1) = top(1, 2)
      
      dn(1) = gdL
      dn(2) = gdR
      
      btm(1, 1) = gdL*costhetab**2.d0 - gdR*sinthetab**2.d0
      btm(1, 2) = (gdL + gdR)*costhetab*sinthetab
      btm(2, 2) = gdR*costhetab**2.d0 - gdL*sinthetab**2.d0
      btm(2, 1) = btm(1, 2)
      
      
      el(1) = geL
      el(2) = geR
      
      tau(1, 1) = geL*costhetatau**2.d0 - geR*sinthetatau**2.d0
      tau(1, 2) = (geL + geR)*costhetatau*sinthetatau
      tau(2, 2) = geR*costhetatau**2.d0 - geL*sinthetatau**2.d0
      tau(2, 1) = tau(1, 2)    

C---------------------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta = datan(tanbeta)
      tan2beta = dtan(2.d0 * beta)
      alpha = 0.5d0*datan(((mA0sq + MZ*MZ)/
     $     (mA0sq - MZ*MZ))*tan2beta)
      cosbeta  = dcos(beta)
      cos2beta = dcos(2.d0 * beta)
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)

!---------------------------------------------------------------------------

      call a0(mt,q,a0mt)
      call a0(mc,q,a0mc)
      call a0(mu,q,a0mu)

      call a0(mb,q,a0mb)
      call a0(ms,q,a0ms)
      call a0(md,q,a0md)

      call a0(me,q,a0me)
      call a0(mmu,q,a0mmu)
      call a0(mtau,q,a0mtau)


      call b0(p,mt,mt,q,b0mtmt)
      call b0(p,mc,mc,q,b0mcmc)
      call b0(p,mu,mu,q,b0mumu)

      call b0(p,mb,mb,q,b0mbmb)
      call b0(p,ms,ms,q,b0msms)
      call b0(p,md,md,q,b0mdmd)

      call b0(p,me,me,q,b0meme)
      call b0(p,mmu,mmu,q,b0mmummu)
      call b0(p,mtau,mtau,q,b0mtaumtau)

      call a0(muL,q,a0muL)
      call a0(mcL,q,a0mcL)
      call a0(muR,q,a0muR)
      call a0(mcR,q,a0mcR)
      call a0(mtL,q,a0mtL)
      call a0(mtR,q,a0mtR)

      call a0(mdL,q,a0mdL)
      call a0(msL,q,a0msL)
      call a0(mdR,q,a0mdR)
      call a0(msR,q,a0msR)
      call a0(mbL,q,a0mbL)
      call a0(mbR,q,a0mbR)

      call a0(meL,q,a0meL)
      call a0(mmuL,q,a0mmuL)
      call a0(meR,q,a0meR)
      call a0(mmuR,q,a0mmuR)
      call a0(mtauL,q,a0mtauL)
      call a0(mtauR,q,a0mtauR)

      call a0(snu(1),q,a0snu1)
      call a0(snu(2),q,a0snu2)
      call a0(snu(3),q,a0snu3)

      call a0(mA0,q,a0mA)
      call a0(mh0,q,a0mh)
      call a0(mhpm,q,a0mhpm)
      call a0(mHu0,q,a0mhu)
 

      call b0(p,muL,muR,q,b0muLmuR)
      call b0(p,muR,muL,q,b0muRmuL)
      call b0(p,mcL,mcR,q,b0mcLmcR)
      call b0(p,mcR,mcL,q,b0mcRmcL)
      call b0(p,mtL,mtR,q,b0mtLmtR)
      call b0(p,mtR,mtL,q,b0mtRmtL)

      call b0(p,mdL,mdR,q,b0mdLmdR)
      call b0(p,mdR,mdL,q,b0mdRmdL)
      call b0(p,msL,msR,q,b0msLmsR)
      call b0(p,msR,msL,q,b0msRmsL)
      call b0(p,mbL,mbR,q,b0mbLmbR)
      call b0(p,mbR,mbL,q,b0mbRmbL)
      
      call b0(p,meL,meR,q,b0meLmeR)
      call b0(p,meR,meL,q,b0meRmeL)
      call b0(p,mmuL,mmuR,q,b0mmuLmmuR)
      call b0(p,mmuR,mmuL,q,b0mmuRmmuL)
      call b0(p,mtauL,mtauR,q,b0mtauLmtauR)
      call b0(p,mtauR,mtauL,q,b0mtauRmtauL)

      
      call f(p,mHpm,MW,q,fmhpmMw)
      call f(p,mhu0,MZ,q,fmhuMz)
      call f(p,mh0,MZ,q,fmhMz)
 
      call a0(MW,q,a0Mw)
      call a0(MZ,q,a0Mz)

      call funcg(p,nmneut(1),nmneut(1),q,gmneut1mneut1)
      call funcg(p,nmneut(1),nmneut(2),q,gmneut1mneut2)
      call funcg(p,nmneut(1),nmneut(3),q,gmneut1mneut3)
      call funcg(p,nmneut(1),nmneut(4),q,gmneut1mneut4)
 
      call funcg(p,nmneut(2),nmneut(1),q,gmneut2mneut1)
      call funcg(p,nmneut(2),nmneut(2),q,gmneut2mneut2)
      call funcg(p,nmneut(2),nmneut(3),q,gmneut2mneut3)
      call funcg(p,nmneut(2),nmneut(4),q,gmneut2mneut4)
 
      call funcg(p,nmneut(3),nmneut(1),q,gmneut3mneut1)
      call funcg(p,nmneut(3),nmneut(2),q,gmneut3mneut2)
      call funcg(p,nmneut(3),nmneut(3),q,gmneut3mneut3)
      call funcg(p,nmneut(3),nmneut(4),q,gmneut3mneut4)

 
      call funcg(p,nmneut(4),nmneut(1),q,gmneut4mneut1)
      call funcg(p,nmneut(4),nmneut(2),q,gmneut4mneut2)
      call funcg(p,nmneut(4),nmneut(3),q,gmneut4mneut3)
      call funcg(p,nmneut(4),nmneut(4),q,gmneut4mneut4)

      call funcg(p,dabs(mchargino(1)),dabs(mchargino(1)),q,gmch1mch1)
      call funcg(p,dabs(mchargino(1)),dabs(mchargino(2)),q,gmch1mch2)
      call funcg(p,dabs(mchargino(2)),dabs(mchargino(1)),q,gmch2mch1)
      call funcg(p,dabs(mchargino(2)),dabs(mchargino(2)),q,gmch2mch2)
            
      call b0(p,dabs(mchargino(1)),dabs(mchargino(1)),q,b0mch1mch1)
      call b0(p,dabs(mchargino(1)),dabs(mchargino(2)),q,b0mch1mch2)
      call b0(p,dabs(mchargino(2)),dabs(mchargino(1)),q,b0mch2mch1)
      call b0(p,dabs(mchargino(2)),dabs(mchargino(2)),q,b0mch2mch2)
 
      call b0(p,nmneut(1),nmneut(1),q,b0mneut1mneut1)
      call b0(p,nmneut(1),nmneut(2),q,b0mneut1mneut2)
      call b0(p,nmneut(1),nmneut(3),q,b0mneut1mneut3)
      call b0(p,nmneut(1),nmneut(4),q,b0mneut1mneut4)

      call b0(p,nmneut(2),nmneut(1),q,b0mneut2mneut1)
      call b0(p,nmneut(2),nmneut(2),q,b0mneut2mneut2)
      call b0(p,nmneut(2),nmneut(3),q,b0mneut2mneut3)
      call b0(p,nmneut(2),nmneut(4),q,b0mneut2mneut4)

      call b0(p,nmneut(3),nmneut(1),q,b0mneut3mneut1)
      call b0(p,nmneut(3),nmneut(2),q,b0mneut3mneut2)
      call b0(p,nmneut(3),nmneut(3),q,b0mneut3mneut3)
      call b0(p,nmneut(3),nmneut(4),q,b0mneut3mneut4)

      call b0(p,nmneut(4),nmneut(1),q,b0mneut4mneut1)
      call b0(p,nmneut(4),nmneut(2),q,b0mneut4mneut2)
      call b0(p,nmneut(4),nmneut(3),q,b0mneut4mneut3)
      call b0(p,nmneut(4),nmneut(4),q,b0mneut4mneut4)

      call b0(p,mA0,mHu0,q,b0mAmHu)
      call b0(p,mA0,mh0,q,b0mAmh)
      call b0(p,MW,mHu0,q,b0MwmHu)
      call b0(p,MZ,mh0,q,b0Mzmh)
      call b0(p,MZ,mHu0,q,b0MzmHu)
      
      
C------------------------------------------------------------------------------
      
      smp = (cosbeta**2.d0) * (3.d0 * ((yuRG(3,3)**2.d0) * !<----------up sector fermions
     $     (p*p * b0mtmt - 2.d0 * a0mt) + (yuRG(2,2)**2.d0) * 
     $     (p*p * b0mcmc - 2.d0 * a0mc) + (yuRG(1,1)**2.d0) * 
     $     (p*p * b0mumu - 2.d0 * a0mu))) + 
     $     (sinbeta**2.d0) * (3.d0 * ((ydRG(3,3)**2.d0) * !<-----------down sector fermions
     $     (p*p * b0mbmb - 2.d0 * a0mb) + (ydRG(2,2)**2.d0) * 
     $     (p*p * b0msms - 2.d0 * a0ms) + (ydRG(1,1)**2.d0) * 
     $     (p*p * b0mdmd - 2.d0 * a0md))) +
     $     (sinbeta**2.d0) * (1.d0 * ((yeRG(3,3)**2.d0) * !<-----------leptons
     $     (p*p * b0mtaumtau - 2.d0 * a0mtau) + (yeRG(2,2)**2.d0) * 
     $     (p*p * b0mmummu - 2.d0 * a0mmu) + (yeRG(1,1)**2.d0) * 
     $     (p*p * b0meme - 2.d0 * a0me))) 
      
 

!---------------------------------------------------------------------------------

      susyp1 = 3.d0 * ((yuRG(1,1)*cosbeta)**2.d0 - 
     $     (((g*g)/(costhw**2.d0)) * !<------up sector sfermions terms
     $     cos2beta * 0.5d0 * guL)) * a0muL +
     $     3.d0 * ((yuRG(2,2)*cosbeta)**2.d0 - 
     $     (((g*g)/(costhw**2.d0)) * !<------up sector sfermions terms
     $     cos2beta * 0.5d0 * guL)) * a0mcL +
     $     3.d0 * ((yuRG(1,1)*cosbeta)**2.d0 - (((g*g)/(costhw**2.d0)) *
     $     cos2beta * 0.5d0 * guR)) * a0muR +
     $     3.d0 * ((yuRG(2,2)*cosbeta)**2.d0 - (((g*g)/(costhw**2.d0)) *
     $     cos2beta * 0.5d0 * guR)) * a0mcR +
     $     3.d0 * ((yuRG(3,3)*cosbeta)**2.d0 - (((g*g)/(costhw**2.d0)) * !<--------stops
     $     cos2beta * 0.5d0 * guL)) * (costhetat**2.d0 * a0mtL + 
     $     (sinthetat)**2.d0 * a0mtR) +
     $     3.d0 * ((yuRG(3,3)*cosbeta)**2.d0 - (((g*g)/(costhw**2.d0)) * 
     $     cos2beta * 0.5d0 * guR)) * ((sinthetat**2.d0) * a0mtL + 
     $     (costhetat**2.d0) * a0mtR) +
     $     3.d0 * ((ydRG(1,1)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) * !<------down sector sfermions term
     $     cos2beta * 0.5d0 * gdL)) * a0mdL +
     $     3.d0 * ((ydRG(2,2)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) * !<------down sector sfermions term
     $     cos2beta * 0.5d0 * gdL)) * a0msL +
     $     3.d0 * ((ydRG(1,1)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) *
     $     cos2beta * 0.5d0 * gdR)) * a0mdR +
     $     3.d0 * ((ydRG(2,2)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) *
     $     cos2beta * 0.5d0 * gdR)) * a0msR +
     $     3.d0 * ((ydRG(3,3)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) * !<--------sbottom
     $     cos2beta * 0.5d0 * gdL)) * (costhetab**2.d0 * a0mbL + 
     $     (sinthetab)**2.d0 * a0mbR) +
     $     3.d0 * ((ydRG(3,3)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) * 
     $     cos2beta * 0.5d0 * gdR)) * ((sinthetab**2.d0) * a0mbL + 
     $     (costhetab**2.d0) * a0mbR) +
     $     1.d0 * ((yeRG(1,1)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) * !<------sleptons
     $     cos2beta * 0.5d0 * geL)) * a0meL +
     $     1.d0 * ((yeRG(2,2)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) * !<------sleptons
     $     cos2beta * 0.5d0 * geL)) * a0mmuL +
     $     1.d0 * ((yeRG(1,1)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) *
     $     cos2beta * 0.5d0 * geR)) * a0meR  +
     $     1.d0 * ((yeRG(2,2)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) *
     $     cos2beta * 0.5d0 * geR)) * a0mmuR +
     $     1.d0 * ((yeRG(3,3)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) * !<--------stau
     $     cos2beta * 0.5d0 * geL)) * ((costhetatau**2.d0) * a0mtauL + 
     $     (sinthetatau)**2.d0 * a0mtauR) +
     $     1.d0 * ((yeRG(3,3)*sinbeta)**2.d0 - (((g*g)/(costhw**2.d0)) * 
     $     cos2beta * 0.5d0 * geR)) * ((sinthetatau**2.d0) * a0mtauL + 
     $     (costhetatau**2.d0) * a0mtauR) -
     $     ((g*g)/costhw**2.d0) * 0.5d0 * gnuL * cos2beta * 
     $     (a0snu1 + a0snu2 +  a0snu3) !<--------sneutrinos



  
!------------------------------------------------------------------------------------------------


c$$$      susyp =  3.d0 * ((yuRG(1,1) * (- sgnmu*modmu * sinbeta -                    !<------------squarks
c$$$     $     AURG(1,1)*cosbeta)))**2.d0 * 0.5d0 * (b0muLmuR + b0muRmuL) + 
c$$$     $     3.d0 * ((yuRG(2,2) * (- sgnmu*modmu*sinbeta - 
c$$$     $     AURG(2,2)*cosbeta))**2.d0) * 0.5d0 * (b0mcLmcR + b0mcRmcL) +
c$$$     $     3.d0 * ((ydRG(1,1) * (- sgnmu*modmu*cosbeta - 
c$$$     $     ADRG(1,1)*sinbeta)))**2.d0 * 0.5d0 * (b0mdLmdR + b0mdRmdL) + 
c$$$     $     3.d0 * ((ydRG(2,2) * (- sgnmu*modmu*cosbeta - 
c$$$     $     ADRG(2,2)*sinbeta))**2.d0) * 0.5d0 * (b0msLmsR + b0msRmsL) +
c$$$     $     3.d0 * (yuRG(3,3) * (sgnmu*modmu*sinbeta + 
c$$$     $     AURG(3,3)*cosbeta))**2.d0 * 0.5d0 * (b0mtLmtR + b0mtRmtL) + 	
c$$$     $     3.d0 * (ydRG(3,3) * (sgnmu*modmu*cosbeta + 
c$$$     $     ADRG(3,3)*sinbeta))**2.d0 * 0.5d0 * (b0mbLmbR + b0mbRmbL) + 	
c$$$     $     1.d0 * (yeRG(1,1) * (- sgnmu*modmu*cosbeta -                          !<------------sleptons
c$$$     $     AERG(1,1)*sinbeta))**2.d0 * 0.5d0 * (b0meLmeR + b0meRmeL) + 
c$$$     $     1.d0 * (yeRG(2,2) * (- sgnmu*modmu*cosbeta - 
c$$$     $     AERG(2,2)*sinbeta))**2.d0 * 0.5d0 * (b0mmuLmmuR + 
c$$$     $     b0mmuRmmuL) +
c$$$     $     1.d0 * (yeRG(3,3) * (- sgnmu*modmu*cosbeta - 
c$$$     $     AERG(3,3)*sinbeta))**2.d0 * 0.5d0 * (b0mtauLmtauR + 
c$$$     $     b0mtauRmtauL)

      susyp =  3.d0 * ((yuRG(1,1) * (sgnmu*modmu * sinbeta +                    !<------------squarks
     $     AURG(1,1)*cosbeta)))**2.d0 * 0.5d0 * (b0muLmuR + b0muRmuL) + 
     $     3.d0 * ((yuRG(2,2) * (sgnmu*modmu*sinbeta + 
     $     AURG(2,2)*cosbeta))**2.d0) * 0.5d0 * (b0mcLmcR + b0mcRmcL) +
     $     3.d0 * ((ydRG(1,1) * (sgnmu*modmu*cosbeta +
     $     ADRG(1,1)*sinbeta)))**2.d0 * 0.5d0 * (b0mdLmdR + b0mdRmdL) + 
     $     3.d0 * ((ydRG(2,2) * (sgnmu*modmu*cosbeta + 
     $     ADRG(2,2)*sinbeta))**2.d0) * 0.5d0 * (b0msLmsR + b0msRmsL) +
     $     3.d0 * (yuRG(3,3) * (sgnmu*modmu*sinbeta + 
     $     AURG(3,3)*cosbeta))**2.d0 * 0.5d0 * (b0mtLmtR + b0mtRmtL) + 	
     $     3.d0 * (ydRG(3,3) * (sgnmu*modmu*cosbeta + 
     $     ADRG(3,3)*sinbeta))**2.d0 * 0.5d0 * (b0mbLmbR + b0mbRmbL) + 	
     $     1.d0 * (yeRG(1,1) * (sgnmu*modmu*cosbeta +                          !<------------sleptons
     $     AERG(1,1)*sinbeta))**2.d0 * 0.5d0 * (b0meLmeR + b0meRmeL) + 
     $     1.d0 * (yeRG(2,2) * (sgnmu*modmu*cosbeta +
     $     AERG(2,2)*sinbeta))**2.d0 * 0.5d0 * (b0mmuLmmuR + 
     $     b0mmuRmmuL) +
     $     1.d0 * (yeRG(3,3) * (sgnmu*modmu*cosbeta + 
     $     AERG(3,3)*sinbeta))**2.d0 * 0.5d0 * (b0mtauLmtauR + 
     $     b0mtauRmtauL)
      

!---------------------------------------------------------------------------------

      smHz = (g*g) * 0.25d0 * (2.d0 * fmHpmMw + (((sin(alpha - beta))/     !<------------check alpha
     $     costhw)**2.d0 * fmhuMz) + 
     $     (((cos(alpha - beta))/costhw)**2.d0 * fmhMz)) + 
     $     ((g*g) * (2.d0 * a0Mw + a0Mz/(costhw**2.d0)))

!--------------------------------------------------------------------------------- !<----checked


      lAas(1) = -cos2beta*cosbeta*(g*MZ/(2.d0*costhw))
      lAas(2) =  cos2beta*sinbeta*(g*MZ/(2.d0*costhw))

      call rmat2d(alpha,ralpha) 
      
      lAah(1) = ralpha(1,1)*lAas(1) +ralpha(1,2)*lAas(2)
      lAah(2) = ralpha(2,1)*lAas(1) +ralpha(2,2)*lAas(2)

      lAgs(1) = -sin2beta*cosbeta*(g*Mz/(2.d0*costhw))                 !<--- trilinear
      lAgs(2) =  sin2beta*sinbeta*(g*Mz/(2.d0*costhw))
      
      lAgh(1) = ralpha(1,1)*lAgs(1) +ralpha(1,2)*lAgs(2)               !<---quartic coupling
      lAgh(2) = ralpha(2,1)*lAgs(1) +ralpha(2,2)*lAgs(2)


      lAa12(1,1) = -cos2beta*(g*g)/(4.d0*costhw**2.d0)
      lAa12(1,2) = 0.d0
      lAa12(2,1) = 0.d0
      lAa12(2,2) =  cos2beta*(g*g)/(4.d0*costhw**2.d0)

      malpha = -1.d0*alpha

      call rmat2d(malpha,mralpha) 

      call mat3prod2d(ralpha,lAa12,mralpha,lAahh)

      lAaaa = (g*g)/(4.d0*costhw**2.d0) * 3.d0 * (cos2beta**2.d0)

      lAagg = (3.d0 * (sin2beta**2.d0) - 1.d0)*(g*g)/(4.d0*costhw**2.d0)

!-----------------------------------------------------------------------------------------------------      
 

      sumHz =  0.5d0 * (2.d0 * ((lAah(2)**2.d0) * b0mAmh + 
     $     (lAah(1)**2.d0) * b0mAmHu + (lAgh(2)**2.d0) * b0Mzmh +
     $     (lAgh(1)**2.d0) * b0MzmHu) +
     $     lAahh(1,1) * a0mHu + lAahh(2,2) * a0mh + 
     $     lAaaa * a0mA + lAagg * a0Mz) +             
     $     (g*g) * (MW*Mw) * 0.5d0 * b0MwmHu + 
     $     ((g*g)/(4.d0*costhw*costhw)) * ((costhw**2.d0 * (1.d0 + 
     $     (sin2beta**2.d0)) - sinsqthw * (cos2beta**2.d0)) * a0Mw + 
     $     (cos2beta**2.d0) * a0mhpm)
      

!----------------------------------------------------------------------------------------
C       Neutralinos
!----------------------------------------------------------------------------------------
      
      loopsii: DO i = 1, 4
      loopsij: DO j = 1, 4
      
      apsinotp1(i,j) = 0.d0
      bpsinotp1(i,j) = 0.d0
      
      apsinotp2(i,j) = 0.d0
      bpsinotp2(i,j) = 0.d0
      
      ENDDO loopsij
      ENDDO loopsii
      
      
      apsinotp1(1,3) =  -gp*0.5d0
      apsinotp2(1,4) =  -gp*0.5d0
      apsinotp1(2,3) =   g*0.5d0
      apsinotp2(2,4) =   g*0.5d0
      
      apsinotp1(3,1) =   apsinotp1(1,3) 
      apsinotp2(4,1) =   apsinotp2(1,4) 
      apsinotp1(3,2) =   apsinotp1(2,3) 
      apsinotp2(4,2) =   apsinotp2(2,4) 
      
      
      loopbii: DO i = 1, 4
      loopbij: DO j = 1, 4
      
      bpsinotp1(i,j) = -1.d0 * apsinotp1(i,j)
      bpsinotp2(i,j) = -1.d0 * apsinotp2(i,j)
      
      ENDDO loopbij
      ENDDO loopbii
      
      call dag4d(ON,ONdag)
      
      call mat3prod4d(ON,apsinotp1,ONdag,achip1)
      call scmul4d(achip1,-1.d0 * sinbeta,achinotp1)

      call mat3prod4d(ON,apsinotp2,ONdag,achip2)
      call scmul4d(achip2,cosbeta,achinotp2)


      call mat3prod4d(ON,bpsinotp1,ONdag,bchip1)
      call scmul4d(bchip1,-1.d0 * sinbeta,bchinotp1)

      call mat3prod4d(ON,bpsinotp2,ONdag,bchip2)
      call scmul4d(bchip2,cosbeta,bchinotp2)
      
      
      call add4d(achinotp1,achinotp2,a2chiA)
      call add4d(bchinotp1,bchinotp2,b2chiA)
          
!-----------------------------------------------------------
      loopnoti: DO i = 1, 4 
      loopnotj: DO j = 1, 4

      fnotA(i,j) = a2chiA(i,j)**2.d0 + b2chiA(i,j)**2.d0
      gnotA(i,j) = 2.d0 * (b2chiA(i,j)*a2chiA(i,j))
      
      ENDDO loopnotj
      ENDDO loopnoti

      neutterm = 0.d0

      hxnot(1,1) = gmneut1mneut1
      hxnot(1,2) = gmneut1mneut2
      hxnot(1,3) = gmneut1mneut3
      hxnot(1,4) = gmneut1mneut4
      hxnot(2,1) = gmneut2mneut1
      hxnot(2,2) = gmneut2mneut2
      hxnot(2,3) = gmneut2mneut3
      hxnot(2,4) = gmneut2mneut4 
      hxnot(3,1) = gmneut3mneut1
      hxnot(3,2) = gmneut3mneut2
      hxnot(3,3) = gmneut3mneut3
      hxnot(3,4) = gmneut3mneut4
      hxnot(4,1) = gmneut4mneut1
      hxnot(4,2) = gmneut4mneut2
      hxnot(4,3) = gmneut4mneut3
      hxnot(4,4) = gmneut4mneut4

      
      bxnot(1,1) = mneut(1)*mneut(1)*b0mneut1mneut1
      bxnot(1,2) = mneut(1)*mneut(2)*b0mneut1mneut2
      bxnot(1,3) = mneut(1)*mneut(3)*b0mneut1mneut3
      bxnot(1,4) = mneut(1)*mneut(4)*b0mneut1mneut4
      bxnot(2,1) = mneut(2)*mneut(1)*b0mneut2mneut1
      bxnot(2,2) = mneut(2)*mneut(2)*b0mneut2mneut2
      bxnot(2,3) = mneut(2)*mneut(3)*b0mneut2mneut3
      bxnot(2,4) = mneut(2)*mneut(4)*b0mneut2mneut4
      bxnot(3,1) = mneut(3)*mneut(1)*b0mneut3mneut1
      bxnot(3,2) = mneut(3)*mneut(2)*b0mneut3mneut2
      bxnot(3,3) = mneut(3)*mneut(3)*b0mneut3mneut3
      bxnot(3,4) = mneut(3)*mneut(4)*b0mneut3mneut4
      bxnot(4,1) = mneut(4)*mneut(1)*b0mneut4mneut1
      bxnot(4,2) = mneut(4)*mneut(2)*b0mneut4mneut2
      bxnot(4,3) = mneut(4)*mneut(3)*b0mneut4mneut3
      bxnot(4,4) = mneut(4)*mneut(4)*b0mneut4mneut4




      loopnlinoi: DO i = 1, 4
      loopnlinoj: DO j = 1, 4


      neutterm = neutterm + 0.5d0 * (fnotA(i,j) *
     $     hxnot(i,j) - 2.d0 * gnotA(i,j) * bxnot(i,j))
        
      ENDDO loopnlinoj
      ENDDO loopnlinoi 

C------------------------------------------------------------------------------
C       Chargino Terms
C------------------------------------------------------------------------------
      
      apsiposp1(1,1) = 0.d0
      apsiposp1(1,2) = g/(dsqrt(2.d0))
      apsiposp1(2,1) = 0.d0
      apsiposp1(2,2) = 0.d0

      apsiposp2(1,1) = 0.d0
      apsiposp2(2,1) = -1.d0 * (g/(dsqrt(2.d0)))
      apsiposp2(1,2) = 0.d0
      apsiposp2(2,2) = 0.d0

      bpsiposp2(1,1) = -1.d0*apsiposp2(1,1)
      bpsiposp2(1,2) = -1.d0*apsiposp2(2,1)
      bpsiposp2(2,1) = -1.d0*apsiposp2(1,2)
      bpsiposp2(2,2) = -1.d0*apsiposp2(2,2)

      bpsiposp1(1,1) = -1.d0*apsiposp1(1,1)
      bpsiposp1(1,2) = -1.d0*apsiposp1(2,1)
      bpsiposp1(2,1) = -1.d0*apsiposp1(1,2)
      bpsiposp1(2,2) = -1.d0*apsiposp1(2,2)

      loopchii: DO i = 1, 2
      loopchij: DO j = 1, 2

      achiposA(i,j) = -sinbeta*apsiposp1(i,j) + cosbeta*apsiposp2(i,j)
      bchiposA(i,j) = -sinbeta*bpsiposp1(i,j) + cosbeta*bpsiposp2(i,j)
      
      ENDDO loopchij
      ENDDO loopchii

      call dag2d(OCL,OCLdag)
      call dag2d(OCR,OCRdag)

c$$$      call mat3prod2d(OCL,achiposA,OCRdag,a2chiposA)    !<--------------------CHECK!!
c$$$      call mat3prod2d(OCR,bchiposA,OCLdag,b2chiposA)

      call mat3prod2d(OCR,achiposA,OCLdag,a2chiposA)    !<--------------------CHECK!!
      call mat3prod2d(OCL,bchiposA,OCRdag,b2chiposA)
      
      loopposi: DO i = 1, 2 
      loopposj: DO j = 1, 2

      fposA(i,j) = a2chiposA(i,j)**2.d0 + b2chiposA(i,j)**2.d0
      gposA(i,j) = 2.d0 * (b2chiposA(i,j)*a2chiposA(i,j))
      
      ENDDO loopposj
      ENDDO loopposi

      charginoterm = 0.d0

      hxpos(1,1) = gmch1mch1
      hxpos(1,2) = gmch1mch2
      hxpos(2,1) = gmch2mch1
      hxpos(2,2) = gmch2mch2

      bxpos(1,1) = mchargino(1)*mchargino(1)*b0mch1mch1
      bxpos(1,2) = mchargino(1)*mchargino(2)*b0mch1mch2
      bxpos(2,1) = mchargino(2)*mchargino(1)*b0mch2mch1
      bxpos(2,2) = mchargino(2)*mchargino(2)*b0mch2mch2
      
      loopchginoi: DO i = 1, 2
      loopchginoj: DO j = 1, 2


      charginoterm = charginoterm + (fposA(i,j) *
     $     hxpos(i,j) - 2.d0 * gposA(i,j) * bxpos(i,j))
        
      ENDDO loopchginoj
      ENDDO loopchginoi 


      piaaT = (charginoterm + neutterm + smp + susyp + smHz + susyp1 
     $            + sumHz)/(16.d0 * pi*pi)


      RETURN

      END SUBROUTINE piaa

C=========================================================================================




C=========================================================================================
C     Charged higgs (H+-) self energy
C
C----------------------------------------------------------------------------------------
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C
C     pihphm is checked @ 22:30 on 20/06/2010.
C
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

      SUBROUTINE pihphm(p,q,g,gp,mt,mb,mtau,tanbeta,AURG,ADRG,AERG,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     yuRG,yeRG,ydRG,mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu,modmu,ON,
     $     OCL,OCR,pihphmT)

      IMPLICIT NONE 

      INTEGER i,j,k
      DOUBLE PRECISION tanbeta,gp
      DOUBLE PRECISION AURG(3,3),ADRG(3,3)
      DOUBLE PRECISION AERG(3,3),SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      DOUBLE PRECISION Neg(4),Ceg(2),modmu,sgnmu
      DOUBLE PRECISION mh0sq,mhu0sq,mhpmsq,mA0sq
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      DOUBLE PRECISION cosbeta, sinbeta, cos2beta
 
      double precision mT, mB, mTau
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),snu(3),yuRG(3,3),ydRG(3,3),nmneut(4)

      DOUBLE PRECISION pihphmT,ans,sinsqthw,g
      DOUBLE PRECISION p,q,mh0,mHu0,mHpm,mA0          
      DOUBLE PRECISION b0mhMw,b0MzMw,b0Mw0,b0mHuMw,b0mAmHpm
      DOUBLE PRECISION alpha,beta
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION thetatz,thetabz,thetatauz

      DOUBLE PRECISION costhetat,sinthetat,costhetab,sinthetab
      DOUBLE PRECISION costhetatau,sinthetatau,sin2beta
      DOUBLE PRECISION yeRG(3,3)

      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetacz,thetasz,thetamuz

c$$$      DOUBLE PRECISION costhetac,sinthetac,costhetas,sinthetas
c$$$      DOUBLE PRECISION costhetamu,sinthetamu

      DOUBLE PRECISION thetau,thetad,thetae
      DOUBLE PRECISION thetauz,thetadz,thetaez

c$$$      DOUBLE PRECISION costhetau,sinthetau,costhetad,sinthetad
c$$$      DOUBLE PRECISION costhetae,sinthetae

      DOUBLE PRECISION smhz,susyp,smp
      double precision higgs

      DOUBLE PRECISION hmumd,hmcms,hmtopmb,h0me,h0mmu,h0mtau
      DOUBLE PRECISION hmneut1ch1,hmneut1ch2,hmneut2ch1,hmneut2ch2
      DOUBLE PRECISION hmneut3ch1,hmneut3ch2,hmneut4ch1,hmneut4ch2

      DOUBLE PRECISION b0mneut1ch1,b0mneut1ch2,b0mneut2ch1,b0mneut2ch2
      DOUBLE PRECISION b0mneut3ch1,b0mneut3ch2,b0mneut4ch1,b0mneut4ch2
      DOUBLE PRECISION b0mHumHpm,b0MwmA
     
      DOUBLE PRECISION gch1mneut1,gch1mneut2,gch1mneut3,gch1mneut4
      DOUBLE PRECISION gch2mneut1,gch2mneut2,gch2mneut3,gch2mneut4

      DOUBLE PRECISION b0ch1mneut1,b0ch1mneut2,b0ch1mneut3,b0ch1mneut4
      DOUBLE PRECISION b0ch2mneut1,b0ch2mneut2,b0ch2mneut3,b0ch2mneut4

      DOUBLE PRECISION a0muL,a0muR,a0mcL,a0mcR,a0msL,a0msR,a0mtauR
      DOUBLE PRECISION a0mdL,a0mdR,a0meL,a0meR,a0mmuL,a0mmuR,a0mtauL
      DOUBLE PRECISION a0snu1,a0snu2,a0mtR,a0mtL,a0mbR,a0mbL,a0snu3 
      DOUBLE PRECISION a0mA,a0mh,a0mHu,a0mHpm,a0Mw,a0Mz
      
      DOUBLE PRECISION fmAMw,fmhMw,fmHpm0,fmHpmMz,fmHuMw,rlHpH0Hm(2,2)
      DOUBLE PRECISION ralpha(2,2),mralpha(2,2),rlHpHmH0H0(2,2)
     
      DOUBLE PRECISION b22tmhMw,b22tmHumHpm,b22tMzMw,b22tMw0
      DOUBLE PRECISION b22tmhmHpm,b22tmHuMw,b22tmAmHpm     

      DOUBLE PRECISION ldHpud(2,2),ldHpen(2,2),rldHpud(2,2),rldHpen(2,2)
      DOUBLE PRECISION lHpH0Hm(2,2),lHpHmH0H0(2,2),lHpHmG0G0
      DOUBLE PRECISION lHpHmA0A0
      DOUBLE PRECISION b0muLmdL,b0mcLmsL,b0snu1meL,b0snu2mmuL,b0mtmb
      DOUBLE PRECISION b0mtLmbL,b0mtLmbR,b0mtRmbL,b0mtRmbR,b0mhmHpm
      DOUBLE PRECISION b0snu3mtauL,b0snu3mtauR,g0mtau,gmtmb
      DOUBLE PRECISION rmthetat(2,2),rmthetab(2,2),rmthetatau(2,2)

      DOUBLE PRECISION fijh(4,2),gijh(4,2)
      DOUBLE PRECISION mchargino(2)

      DOUBLE PRECISION apsiposh1(4,2),bpsiposh1(4,2),achi(4,2),bchi(4,2)
      DOUBLE PRECISION achiposh(4,2),bchiposh(4,2)
      double precision apsiposh2(4,2),bpsiposh2(4,2)

      DOUBLE PRECISION npchino,hxnoth(2,4),bxnoth(2,4)
      DOUBLE PRECISION OCRdag(2,2),OCLdag(2,2)
      DOUBLE PRECISION ON(4,4),ONdag(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION mu,mtop,tan2beta
      DOUBLE PRECISION sinsqthw_susy


      double precision mbpole,mtaupole,Mtpole,MZpole,pi,MZ,MWpole,MW
      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole


      common/sfmixing_mz/thetatz,thetabz,thetatauz,thetacz,thetasz,
     $     thetamuz,thetauz,thetadz,thetaez

      common/sinsq_susy/sinsqthw_susy
            
      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

      EXTERNAL b0,b22t,h,dag2d,theta,rmat2d,f,funcg


      include 'stdinputs.h'
!---------------------------------------------------------------------------

      pi = 4.d0 * datan(1.d0)

      MW = MWpole
      MZ = MZpole

      mu = mUQ
      mtop = mt

      sinsqthw = sinsqthw_susy !1.d0 - (MW/MZ)**2.d0
      sinthw = dsqrt(dabs(sinsqthw))
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

      piHphmT = 0.d0

!----------------------------------------------------------------------------

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

      

      costhetat = dcos(thetat)
      sinthetat = dsin(thetat)
      costhetab = dcos(thetab)
      sinthetab = dsin(thetab)
      costhetatau = dcos(thetatau)
      sinthetatau = dsin(thetatau)

C----------------------------------------------------------------

      mh0  = dsqrt((mh0sq))
      mhu0 = dsqrt((mhu0sq))
      mHpm = dsqrt((mHpmsq))
      mA0  = dsqrt((mA0sq))

      beta = datan(tanbeta)
      tan2beta = dtan(2.d0 * beta)
      alpha  = 0.5d0 * datan(((mA0sq + MZ*MZ)/
     $     (mA0sq - MZ*MZ))*tan2beta)
      cosbeta  = dcos(beta)
      cos2beta = dcos(2.d0 * beta)
      sinbeta  = dsin(beta)
      sin2beta = dsin(2.d0 * beta)


C---------------------------------------------------------------    

      call b0(p,mh0,MW,q,b0mhMw)
      call b0(p,MZ,MW,q,b0MzMw)
      call b0(p,MW,0.d0,q,b0Mw0)
      call b0(p,mA0,mHpm,q,b0mAmHpm)
      
      call h(p,mu,md,q,hmumd)
      call h(p,mc,ms,q,hmcms)
      call h(p,mtop,mb,q,hmtopmb)
      call h(p,0.d0,me,q,h0me)
      call h(p,0.d0,mmu,q,h0mmu)
      call h(p,0.d0,mtau,q,h0mtau)


      call h(p,nmneut(1),dabs(mchargino(1)),q,hmneut1ch1)
      call h(p,nmneut(1),dabs(mchargino(2)),q,hmneut1ch2)
      call h(p,nmneut(2),dabs(mchargino(1)),q,hmneut2ch1)
      call h(p,nmneut(2),dabs(mchargino(2)),q,hmneut2ch2)
      call h(p,nmneut(3),dabs(mchargino(1)),q,hmneut3ch1)
      call h(p,nmneut(3),dabs(mchargino(2)),q,hmneut3ch2)
      call h(p,nmneut(4),dabs(mchargino(1)),q,hmneut4ch1)
      call h(p,nmneut(4),dabs(mchargino(2)),q,hmneut4ch2)

      call b0(p,nmneut(1),dabs(mchargino(1)),q,b0mneut1ch1)
      call b0(p,nmneut(1),dabs(mchargino(2)),q,b0mneut1ch2)
      call b0(p,nmneut(2),dabs(mchargino(1)),q,b0mneut2ch1)
      call b0(p,nmneut(2),dabs(mchargino(2)),q,b0mneut2ch2)
      call b0(p,nmneut(3),dabs(mchargino(1)),q,b0mneut3ch1)
      call b0(p,nmneut(3),dabs(mchargino(2)),q,b0mneut3ch2)
      call b0(p,nmneut(4),dabs(mchargino(1)),q,b0mneut4ch1)
      call b0(p,nmneut(4),dabs(mchargino(2)),q,b0mneut4ch2)

!-------------------------------------------------------------------
      
      call b22t(p,mh0,MW,q,b22tmhMw)
      call b22t(p,mHu0,mHpm,q,b22tmHumHpm)
      call b22t(p,MZ,MW,q,b22tMzMw)
      call b22t(p,MW,0.d0,q,b22tMw0)
      call b22t(p,mh0,mHpm,q,b22tmhmHpm)
      call b22t(p,mHu0,MW,q,b22tmHuMw)
      call b22t(p,mA0,mHpm,q,b22tmAmHpm)


      call b0(p,muL,mdL,q,b0muLmdL)
      call b0(p,mcL,msL,q,b0mcLmsL)
      call b0(p,mtL,mbL,q,b0mtLmbL)
      call b0(p,mtL,mbR,q,b0mtLmbR)
      call b0(p,mtR,mbL,q,b0mtRmbL)
      call b0(p,mtR,mbR,q,b0mtRmbR)
     
            
      call b0(p,snu(1),meL,q,b0snu1meL)     
      call b0(p,snu(2),mmuL,q,b0snu2mmuL)

      call b0(p,snu(3),mtauL,q,b0snu3mtauL)
      call b0(p,snu(3),mtauR,q,b0snu3mtauR)

!--------------------------------------------------------------------------

      call a0(muR,q,a0muR)
      call a0(muL,q,a0muL)
      call a0(mcR,q,a0mcR)
      call a0(mcL,q,a0mcL)
      call a0(mtR,q,a0mtR)
      call a0(mtL,q,a0mtL)

      call a0(msR,q,a0msR)
      call a0(msL,q,a0msL)
      call a0(mdR,q,a0mdR)
      call a0(mdL,q,a0mdL)
      call a0(mbR,q,a0mbR)
      call a0(mbL,q,a0mbL)

      call a0(meR,q,a0meR)
      call a0(meL,q,a0meL)
      call a0(mmuR,q,a0mmuR)
      call a0(mmuL,q,a0mmuL)
      call a0(mtauR,q,a0mtauR)
      call a0(mtauL,q,a0mtauL)

      call a0(snu(1),q,a0snu1)
      call a0(snu(2),q,a0snu2)
      call a0(snu(3),q,a0snu3)

      call b0(p,mt,mb,q,b0mtmb)
      call funcg(p,0.d0,mtau,q,g0mtau)
      call funcg(p,mt,mb,q,gmtmb)
      
      
      call a0(mA0,q,a0mA)
      call a0(mh0,q,a0mh)
      call a0(mHu0,q,a0mHu)
      call a0(mHpm,q,a0mHpm)
      call a0(MW,q,a0Mw)
      call a0(MZ,q,a0Mz)

      call b0(p,mHu0,mHpm,q,b0mHumHpm)
      call b0(p,MW,mA0,q,b0MwmA)
      call b0(p,mh0,mHpm,q,b0mhmHpm)

C---------------------------------------------------------------------------    
      
      smp = 3.d0 * ((yuRG(3,3) * yuRG(3,3) * cosbeta * cosbeta) + !<-----first two generation => neglected 
     $     (ydRG(3,3) * ydRG(3,3) * sinbeta * sinbeta)) * gmtmb + 
     $     sinbeta * sinbeta * (yeRG(3,3)*yeRG(3,3)) * g0mtau -
     $     6.d0 * ydRG(3,3) * yuRG(3,3) * mt * mb * sin2beta * b0mtmb

C---------------------------------------------------------------------------
            
      ldHpud(1, 1) = (g*MW*sin2beta/dsqrt(2.d0)) - 
     $     yuRG(3,3) * mt * cosbeta - ydRG(3,3) * mb * sinbeta

      ldHpud(1, 2) = ydRG(3,3)*(- sgnmu*modmu*cosbeta - 
     $     ADRG(3,3)*sinbeta)

      ldHpud(2, 1) = yuRG(3,3)*(- sgnmu*modmu*sinbeta - 
     $     AURG(3,3)*cosbeta)

      ldHpud(2, 2) = - 1.d0 * (yuRG(3,3) * mb * cosbeta + 
     $     ydRG(3,3) * mt * sinbeta)
      
      call rmat2d(thetat,rmthetat)
      call rmat2d(-1.d0 * thetab,rmthetab)

      call mat3prod2d(rmthetat,ldHpud,rmthetab,rldHpud)
      
!-----------------------------------------------------------------------

      ldHpen(1, 1) = (g*MW*sin2beta/dsqrt(2.d0)) - 
     $     yeRG(3,3) * mtau * sinbeta
      ldHpen(1, 2) = yeRg(3,3) * (- sgnmu*modmu*cosbeta - 
     $     AERG(3,3)*sinbeta)
      ldHpen(2, 1) = 0.d0
      ldHpen(2, 2) = 0.d0

      call rmat2d(-1.d0 * thetatau,rmthetatau)

      call matmult2d(ldHpen,rmthetatau,rldHpen)

!---------------------------------------------------------------------------------------


      susyp = 3.d0 * ((g*MW*sin2beta/dsqrt(2.d0))**2.d0) * (b0muLmdL + !<---first two families
     $     b0mcLmsL) + 
     $     ((g*MW*sin2beta/dsqrt(2.d0))**2.d0) * (b0snu1meL + 
     $     b0snu2mmuL) +
     $     3.d0 * ((rldHpud(1,1)**2.d0) * b0mtLmbL +                   !<-----third family
     $     (rldHpud(1,2)**2.d0) * b0mtLmbR + 
     $     (rldHpud(2,1)**2.d0) * b0mtRmbL +
     $     (rldHpud(2,2)**2.d0) * b0mtRmbR) +
     $     ((rldHpen(1, 1)**2.d0) * b0snu3mtauL +  
     $     (rldHpen(1, 2)**2.d0) * b0snu3mtauR) +
     $     3.d0 * (g*g) * cos2beta *                                   !<----squark 1st and 2nd family 
     $     (-0.5d0 * (guL/costhw**2.d0) + 0.5d0) * (a0muL + a0mcL) +
     $     3.d0 * (g*g) * cos2beta * (-0.5d0 * (guR/costhw**2.d0)) * 
     $     (a0muR + a0mcR) +
     $     3.d0 * (g*g) * cos2beta * (-0.5d0 * (gdL/costhw**2.d0) - 
     $     0.5d0) * (a0mdL + a0msL) +
     $     3.d0 * (g*g) * cos2beta * (-0.5d0 * (gdR/costhw**2.d0)) * 
     $     (a0mdR + a0msR) +
     $     (g*g) * cos2beta * (-0.5d0 * (gnuL/costhw**2.d0) + 0.5d0 ) * !<--------slepton sector  1st and 2nd family
     $     (a0snu1 + a0snu2) +
     $     (g*g) * cos2beta * (-0.5d0 * (geL/costhw**2.d0) - 0.5d0) * 
     $     (a0meL + a0mmuL) + 
     $     (g*g) * cos2beta * (-0.5d0 * (geR/costhw**2.d0)) * 
     $     (a0meR + a0mmuR) +
     $     3.d0 * ((ydRG(3,3)**2.d0) * sinbeta**2.d0 - !<---------third family sup
     $     (g*g) * cos2beta * 0.5d0 * (guL/costhw**2.d0) + 
     $     0.5d0 * cos2beta * (g*g)) * 
     $     ((costhetat**2.d0) * a0mtL + (sinthetat**2.d0) * a0mtR) +
     $     3.d0 * ((yuRG(3,3)**2.d0) * cosbeta*cosbeta - 
     $     (g*g) * cos2beta * 0.5d0 * (guR/costhw**2.d0)) * 
     $     ((sinthetat**2.d0) * a0mtL + (costhetat**2.d0) * a0mtR) +
     $     3.d0 * ((yuRG(3,3)**2.d0) * cosbeta**2.d0 - 
     $     (g*g) * cos2beta * 0.5d0 * (gdL/costhw**2.d0) - 
     $     0.5d0 * cos2beta * (g*g)) * 
     $     ((costhetab**2.d0) * a0mbL + (sinthetab**2.d0) * a0mdR) +
     $     3.d0 * ((ydRG(3,3)**2.d0) * sinbeta**2.d0 - 
     $     (g*g) * cos2beta * 0.5d0 * (gdR/costhw**2.d0)) *
     $     ((sinthetab**2.d0) * a0mbL + (costhetab**2.d0) * a0mbR) +
     $     ((yeRG(3,3)**2.d0) * sinbeta**2.d0 - 
     $     (g*g) * cos2beta * 0.5d0 * (gnuL/costhw**2.d0) + 
     $     0.5d0 * cos2beta * (g*g)) * a0snu3 + 
     $     (-(g*g) * cos2beta * 0.5d0 * (geL/costhw**2.d0) - 
     $     0.5d0 * cos2beta * (g*g)) * 
     $     ((costhetatau**2.d0) * a0mtauL + 
     $     (sinthetatau**2.d0) * a0mtauR) +
     $     ((yeRG(3,3)**2.d0) * sinbeta**2.d0 - 
     $     (g*g) * cos2beta * 0.5d0 * (geR/costhw**2.d0)) * 
     $     ((sinthetatau**2.d0) * a0mtauL + 
     $     (costhetatau**2.d0) *a0mtauR)

!-------------------------------------------------------------------------------- <--- checked
      call f(p,mA0,MW,q,fmAMw)
      call f(p,mh0,MW,q,fmhMw)
      call f(p,mHpm,0.d0,q,fmHpm0)
      call f(p,mHpm,MZ,q,fmHpmMz)
      call f(p,mHu0,MW,q,fmHuMw)

      call b0(p,mHu0,MW,q,b0mHuMw)


      smHz = (g*g) * 0.25d0 * ((sin(alpha - beta))**2.d0 * fmHuMw +
     $     (cos(alpha - beta))**2.d0 * fmhMw + fmAMw +
     $     (cos2thw/costhw)**2.d0 * fmHpmMz) +
     $     (g*g) * sinsqthw * fmHpm0 + 2.d0 * (g*g) * a0Mw +
     $     (g*g) * (cos2thw/costhw)**2.d0 * a0Mz +
     $     (g*g) * (MW*MW) * 0.25d0 * b0MwmA

!---------------------------------------------------------------------------------

      lHpH0Hm(1, 1) = (-sin2beta*cosbeta + costhw**2.d0 * sinbeta) *
     $     (g * MZ * 0.5d0/costhw)
      lHpH0Hm(1, 2) = (sin2beta*sinbeta - costhw**2.d0 * cosbeta) *
     $     (g * MZ * 0.5d0/costhw)
      lHpH0Hm(2, 1) = (-cos2beta*cosbeta + 
     $     2.d0 * costhw**2.d0 * cosbeta) * (g * MZ * 0.5d0/costhw)
      lHpH0Hm(2, 2) =  (cos2beta*sinbeta + 
     $     2.d0 * costhw**2.d0 * sinbeta) * (g * MZ * 0.5d0/costhw)


      call rmat2d(alpha,ralpha)
      call matmult2d(ralpha,lHpH0Hm,rlHpH0Hm)

!-------------------------------------------------------------------------------

      lHpHmH0H0(1,1) = costhw**2.d0 - sinsqthw*cos2beta
      lHpHmH0H0(1,2) = costhw**2.d0 * sin2beta
      lHpHmH0H0(2,1) = lHpHmH0H0(1,2) 
      lHpHmH0H0(2,2) = costhw**2.d0 + sinsqthw*cos2beta

      call rmat2d(-1.d0 * alpha,mralpha)
      call mat3prod2d(ralpha,lHpHmH0H0,mralpha,rlHpHmH0H0)

!------------------------------------------------------------------------------

      lHpHmG0G0 = (costhw**2.d0 * (1.d0 + (sin2beta**2.d0)) - 
     $     sinsqthw * (cos2beta**2.d0)) * 0.25d0 *((g*g)/(costhw**2.d0))

      lHpHmA0A0 = (cos2beta**2.d0) * 0.25d0 * ((g*g)/(costhw**2.d0))


      higgs = 0.d0
      higgs = (rlHpH0Hm(1,1)**2.d0) * b0mHuMw +
     $     (rlHpH0Hm(1,2)**2.d0) * b0mHumHpm +
     $     (rlHpH0Hm(2,1)**2.d0) * b0mhMw +
     $     (rlHpH0Hm(2,2)**2.d0) * b0mhmHpm +
     $     (2.d0*(sin2beta**2.d0) - 1.d0) * (g*g)*(0.25d0/costhw**2.d0)*
     $     a0MW + (2.d0*(cos2beta**2.d0) * (g*g)*(0.25d0/costhw**2.d0))*
     $     a0mHpm +
     $     0.5d0 * (rlHpHmH0H0(1, 1) * (0.25d0*(g*g)/(costhw**2.d0)) *
     $     a0mHu + rlHpHmH0H0(2, 2) * (0.25d0*(g*g)/(costhw**2.d0)) * 
     $     a0mh + lHpHmG0G0 * a0mz + lHpHmA0A0 * a0mA)


!------------------------------------------------------------------------------
C     Neutralinos and Charginos
!------------------------------------------------------------------------------

      loopsii: DO i = 1, 4
      loopsij: DO j = 1, 2

      apsiposh1(i,j) = 0.d0
      bpsiposh1(i,j) = 0.d0

      apsiposh2(i,j) = 0.d0
      bpsiposh2(i,j) = 0.d0

      achiposh(i,j) = 0.d0
      bchiposh(i,j) = 0.d0

      ENDDO loopsij
      ENDDO loopsii


      apsiposh1(1,2) = gp/dsqrt(2.d0)
      apsiposh1(2,2) = g/dsqrt(2.d0)
      apsiposh1(3,1) = -1.d0*g

      bpsiposh2(1,2) = gp /dsqrt(2.d0)
      bpsiposh2(2,2) = g/dsqrt(2.d0)
      bpsiposh2(4,1) = g



      call dag2d(ON,ONdag)
      call dag2d(OCL,OCLdag)
      call dag2d(OCR,OCRdag)

!-------------------------------------------------------------- 

      loop4by2i: DO i = 1, 4
      loop4by2j: DO j = 1, 2
      achi(i,j) = 0.d0
      bchi(i,j) = 0.d0
      loop4by2k: DO k = 1, 4


      achi(i,j) = achi(i,j) + ON(i,k)*apsiposh1(k,j)
      bchi(i,j) = bchi(i,j) + ON(i,k)*bpsiposh2(k,j)


      ENDDO loop4by2k
      ENDDO loop4by2j
      ENDDO loop4by2i

!---------------------------------------------------------------

      loop2i: DO i = 1, 4
      loop2j: DO j = 1, 2
      achiposh(i,j) = 0.d0
      bchiposh(i,j) = 0.d0
      loop2k: DO k = 1, 2

      achiposh(i,j) = achiposh(i,j) + achi(i,k)*OCLdag(k,j)
      bchiposh(i,j) = bchiposh(i,j) + bchi(i,k)*OCRdag(k,j)

      ENDDO loop2k
      ENDDO loop2j
      ENDDO loop2i

!----------------------------------------------------------------

      loopnoti: DO i = 1, 4 
      loopnotj: DO j = 1, 2

      fijh(i,j) = (achiposh(i,j)*sinbeta)**2.d0 + 
     $     (bchiposh(i,j)*cosbeta)**2.d0
      gijh(i,j) = 2.d0 * (bchiposh(i,j) * achiposh(i,j) * (-1.d0) *
     $     cosbeta * sinbeta)

      ENDDO loopnotj
      ENDDO loopnoti


!--------------------------------------------------------------------------
      call funcg(p,dabs(mchargino(1)),nmneut(1),q,gch1mneut1)
      call funcg(p,dabs(mchargino(1)),nmneut(2),q,gch1mneut2)
      call funcg(p,dabs(mchargino(1)),nmneut(3),q,gch1mneut3)
      call funcg(p,dabs(mchargino(1)),nmneut(4),q,gch1mneut4)
      call funcg(p,dabs(mchargino(2)),nmneut(1),q,gch2mneut1)
      call funcg(p,dabs(mchargino(2)),nmneut(2),q,gch2mneut2)
      call funcg(p,dabs(mchargino(2)),nmneut(3),q,gch2mneut3)
      call funcg(p,dabs(mchargino(2)),nmneut(4),q,gch2mneut4)

      call b0(p,dabs(mchargino(1)),nmneut(1),q,b0ch1mneut1)
      call b0(p,dabs(mchargino(1)),nmneut(2),q,b0ch1mneut2)
      call b0(p,dabs(mchargino(1)),nmneut(3),q,b0ch1mneut3)
      call b0(p,dabs(mchargino(1)),nmneut(4),q,b0ch1mneut4)
      call b0(p,dabs(mchargino(2)),nmneut(1),q,b0ch2mneut1)
      call b0(p,dabs(mchargino(2)),nmneut(2),q,b0ch2mneut2)
      call b0(p,dabs(mchargino(2)),nmneut(3),q,b0ch2mneut3)
      call b0(p,dabs(mchargino(2)),nmneut(4),q,b0ch2mneut4)



      hxnoth(1,1) = gch1mneut1
      hxnoth(1,2) = gch1mneut2
      hxnoth(1,3) = gch1mneut3
      hxnoth(1,4) = gch1mneut4
      hxnoth(2,1) = gch2mneut1
      hxnoth(2,2) = gch2mneut2
      hxnoth(2,3) = gch2mneut3
      hxnoth(2,4) = gch2mneut4


      bxnoth(1,1) = mneut(1)*mchargino(1)*b0ch1mneut1
      bxnoth(1,2) = mneut(2)*mchargino(1)*b0ch1mneut2
      bxnoth(1,3) = mneut(3)*mchargino(1)*b0ch1mneut3
      bxnoth(1,4) = mneut(4)*mchargino(1)*b0ch1mneut4
      bxnoth(2,1) = mneut(1)*mchargino(2)*b0ch2mneut1
      bxnoth(2,2) = mneut(2)*mchargino(2)*b0ch2mneut2
      bxnoth(2,3) = mneut(3)*mchargino(2)*b0ch2mneut3
      bxnoth(2,4) = mneut(4)*mchargino(2)*b0ch2mneut4



      npchino = 0.d0

      loopnlinoi: DO i = 1, 4
      loopnlinoj: DO j = 1, 2


      npchino = npchino + (fijh(i,j) *
     $     hxnoth(j,i) - 2.d0 * gijh(i,j) * bxnoth(j,i))

      ENDDO loopnlinoj
      ENDDO loopnlinoi 
C------------------------------------------------------------------------------

      ans = ( smHz + smp  + npchino + higgs + susyp)

      pihphmT = ans/(16.d0*(pi**2.d0))



      RETURN

      END SUBROUTINE piHphm
C============================================================================================

!===============================================================================

      subroutine pizgamma(p,q,g,mt,mb,mtau,SUegg,SDegg,
     $     SLegg,Ceg,mhpmsq,OCL,OCR,alphhat,result)

      implicit none
      
      integer i

      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6)
      DOUBLE PRECISION Ceg(2),mhpmsq
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
 
      double precision mT, mB, mTau
      DOUBLE PRECISION muL,muR,mcL,mcR,mtL,mtR      
      DOUBLE PRECISION mdL,mdR,msL,msR,mbL,mbR
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR

      DOUBLE PRECISION result,g,e,eu,ed,el,alphhat
      DOUBLE PRECISION p,q,mHpm
      DOUBLE PRECISION thw,costhw,cos2thw,sinthw,cossqthw

      DOUBLE PRECISION thetatz,thetabz,thetatauz
      DOUBLE PRECISION thetacz,thetasz,thetamuz
      DOUBLE PRECISION thetauz,thetadz,thetaez
      
      double precision OCL(2,2),OCR(2,2),mchar(2)

      double precision smterm,higgsterm,charginoterm,sfermionterm
      data smterm/ 0.d0/, higgsterm/ 0.d0/, charginoterm/ 0.d0/,
     $     sfermionterm/0.d0/

      double precision b22tMWMW,b22tmtmt,b22tmcmc,b22tmumu,b22tmbmb
      double precision b22tmsms,b22tmdmd,b22tmtaumtau,b22tmmummu
      double precision b22tmeme,b0MWMW,b0mtmt,b0mcmc,b0mumu,b0mbmb
      double precision b0msms,b0mdmd,b0mtaumtau,b0mmummu,b0meme

      double precision b22tmhpmhp,b22tmcharchar(2),b0mcharchar(2)
      data b22tmhpmhp/ 0.d0/,b22tmcharchar/ 2* 0.d0/,
     $     b0mcharchar/ 2* 0.d0/ 

      double precision csqtht,ssqtht,csqthb,ssqthb,csqthtau,ssqthtau
      double precision csqthc,ssqthc,csqths,ssqths,csqthmu,ssqthmu
      double precision csqthu,ssqthu,csqthd,ssqthd,csqthe,ssqthe

      DOUBLE PRECISION b22tmuLmuL,b22tmuRmuR,b22tmcLmcL,b22tmcRmcR
      DOUBLE PRECISION b22tmtLmtL,b22tmtRmtR

      DOUBLE PRECISION b22tmdLmdL,b22tmdRmdR,b22tmsLmsL,b22tmsRmsR
      DOUBLE PRECISION b22tmbLmbL,b22tmbRmbR

      DOUBLE PRECISION b22tmeLmeL,b22tmeRmeR,b22tmmuLmmuL,b22tmmuRmmuR
      DOUBLE PRECISION b22tmtauLmtauL,b22tmtauRmtauR
      double precision b0mcharmchar(2),b22tmcharmchar(2), sinsqthw_mz

      double precision pi,MWpole,MW
      
!      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole


      common/sfmixing_mz/thetatz,thetabz,thetatauz,thetacz,thetasz,
     $     thetamuz,thetauz,thetadz,thetaez

      common/sinsq_mz/sinsqthw_mz

      include "stdinputs.h"

!--------------------

      MW = MWpole

      sinthw   = dsqrt(sinsqthw_mz)
      thw      = dasin(sinthw)
      costhw   = dcos(thw)
      cos2thw  = dcos(2.d0*thw)
      cossqthw = (dcos(thw))**2.d0
      pi = 4.d0 * datan(1.d0)
      e = dsqrt(alphhat * 4.d0 * pi)
      eu = (2.d0/3.d0)
      ed = - (1.d0/3.d0)
      el = - 1.d0

      gnuL = 0.5d0
      guL = 0.5d0 - 2.d0*sinsqthw_mz/ 3.d0 
      gdL = -0.5d0 + sinsqthw_mz/3.d0
      geL = -0.5d0 + sinsqthw_mz
      guR =  2.d0*sinsqthw_mz/3.d0 
      gdR = -sinsqthw_mz/3.d0
      geR = -sinsqthw_mz 
      yuL = 1.d0/3.d0 
      yuR = -4.d0/3.d0 
      ydL = 1.d0/3.d0 
      ydR = 2.d0/3.d0 
      yeL = -1.d0 
      yeR = 2.d0 
      ynuL = -1.d0 

     
      result = 0.d0

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

      mchar(1) = Ceg(1)
      mchar(2) = Ceg(2)

      mHpm = dsqrt((mHpmsq))

!------------------------------------------------------------------------

      call b22t(p,Mw,Mw,q,b22tMWMW)
      call b0(p,Mw,Mw,q,b0MWMW)

      call b22t(p,mt,mt,q,b22tmtmt)
      call b22t(p,mc,mc,q,b22tmcmc)
      call b22t(p,muq,muq,q,b22tmumu)

      call b22t(p,mb,mb,q,b22tmbmb)
      call b22t(p,ms,ms,q,b22tmsms)
      call b22t(p,md,md,q,b22tmdmd)

      call b22t(p,mtau,mtau,q,b22tmtaumtau)
      call b22t(p,mmu,mmu,q,b22tmmummu)
      call b22t(p,me,me,q,b22tmeme)

      call b0(p,mt,mt,q,b0mtmt)
      call b0(p,mc,mc,q,b0mcmc)
      call b0(p,muq,muq,q,b0mumu)

      call b0(p,mb,mb,q,b0mbmb)
      call b0(p,ms,ms,q,b0msms)
      call b0(p,md,md,q,b0mdmd)

      call b0(p,mtau,mtau,q,b0mtaumtau)
      call b0(p,mmu,mmu,q,b0mmummu)
      call b0(p,me,me,q,b0meme)


      smterm = 0.d0 + 
     $     (12.d0 * sinsqthw_mz - 10.d0) * b22tMWMW - 2.d0 * (MW*MW + 
     $     2.d0 * (1.d0 - sinsqthw_mz) * p*p) * b0MWMW +
     $     ((3.d0 * ((eu * (guL - guR) * ((4.d0 * (b22tmtmt + b22tmcmc + 
     $     b22tmumu)) + (p*p * (b0mtmt + b0mcmc + b0mumu)))) +
     $     (ed * (gdL - gdR) * ((4.d0 * (b22tmdmd + b22tmsms + 
     $     b22tmdmd)) + (p*p * (b0mbmb + b0msms + b0mdmd)))))) +
     $     (el * (geL - geR) * ((4.d0 * (b22tmtaumtau + b22tmmummu + 
     $     b22tmeme)) + p*p * (b0mtaumtau + b0mmummu + b0meme))))

!-------------------------
      
      call b22t(p,mHpm,mHpm,q,b22tmhpmhp)

      higgsterm =  2.d0 * cos2thw * b22tmhpmhp

!------------------------
      charginoterm = 0.d0

      loopchar: do i = 1, 2

      call b22t(p,dabs(mchar(i)),dabs(mchar(i)),q,b22tmcharmchar(i))
      call b0(p,dabs(mchar(i)),dabs(mchar(i)),q,b0mcharmchar(i))

      charginoterm = charginoterm + 0.5d0 * ((OCR(i,1)**2.d0 + 
     $     OCL(i,1)**2.d0 + 2.d0 * cos2thw) * (4.d0 * 
     $     b22tmcharmchar(i) + p*p * b0mcharmchar(i)))

      enddo loopchar

!------------------------

      call b22t(p,muL,muL,q,b22tmuLmuL)
      call b22t(p,muR,muR,q,b22tmuRmuR)
      call b22t(p,mcL,mcL,q,b22tmcLmcL)
      call b22t(p,mcR,mcR,q,b22tmcRmcR)
      call b22t(p,mtL,mtL,q,b22tmtLmtL)
      call b22t(p,mtR,mtR,q,b22tmtRmtR)

      call b22t(p,mdL,mdL,q,b22tmdLmdL)
      call b22t(p,mdR,mdR,q,b22tmdRmdR)
      call b22t(p,msL,msL,q,b22tmsLmsL)
      call b22t(p,msR,msR,q,b22tmsRmsR)
      call b22t(p,mbL,mbL,q,b22tmbLmbL)
      call b22t(p,mbR,mbR,q,b22tmbRmbR)
      
      call b22t(p,meL,meL,q,b22tmeLmeL)
      call b22t(p,meR,meR,q,b22tmeRmeR)
      call b22t(p,mmuL,mmuL,q,b22tmmuLmmuL)
      call b22t(p,mmuR,mmuR,q,b22tmmuRmmuR)
      call b22t(p,mtauL,mtauL,q,b22tmtauLmtauL)
      call b22t(p,mtauR,mtauR,q,b22tmtauRmtauR)

      csqtht = (dcos(thetatz))**2.d0
      ssqtht = (dsin(thetatz))**2.d0
      csqthb = (dcos(thetabz))**2.d0
      ssqthb = (dsin(thetabz))**2.d0
      csqthtau = (dcos(thetatauz))**2.d0
      ssqthtau = (dsin(thetatauz))**2.d0

      csqthc = (dcos(thetacz))**2.d0
      ssqthc = (dsin(thetacz))**2.d0
      csqths = (dcos(thetasz))**2.d0
      ssqths = (dsin(thetasz))**2.d0
      csqthmu = (dcos(thetamuz))**2.d0
      ssqthmu = (dsin(thetamuz))**2.d0

      csqthu = (dcos(thetauz))**2.d0
      ssqthu = (dsin(thetauz))**2.d0
      csqthd = (dcos(thetadz))**2.d0
      ssqthd = (dsin(thetadz))**2.d0
      csqthe = (dcos(thetaez))**2.d0
      ssqthe = (dsin(thetaez))**2.d0


      sfermionterm = 0.d0 +  4.d0 * 
     $     ((3.d0 * ((eu * 
     $     (((guL * csqtht - ssqtht * guR) * b22tmtLmtL + 
     $     (guL * ssqtht - csqtht * guR) * b22tmtRmtR) + 
     $     ((guL * csqthc - ssqthc * guR) * b22tmcLmcL + 
     $     (guL * ssqthc - csqthc * guR) * b22tmcRmcR) + 
     $     ((guL * csqthu - ssqthu * guR) * b22tmuLmuL + 
     $     (guL * ssqthu - csqthu * guR) * b22tmuRmuR))) + 
     $     (ed * 
     $     (((gdL * csqthb - ssqthb * gdR) * b22tmbLmbL + 
     $     (gdL * ssqthb - csqthb * gdR) * b22tmbRmbR) + 
     $     ((gdL * csqths - ssqths * gdR) * b22tmsLmsL + 
     $     (gdL * ssqths - csqths * gdR) * b22tmsRmsR) + 
     $     ((gdL * csqthd - ssqthd * gdR) * b22tmdLmdL + 
     $     (gdL * ssqthd - csqthd * gdR) * b22tmdRmdR))))) +
     $     (1.d0 * (el * 
     $     (((geL * csqthtau - ssqthtau * geR) * b22tmtauLmtauL + 
     $     (geL * ssqthtau - csqthtau * geR) * b22tmtauRmtauR) + 
     $     ((geL * csqthmu - ssqthmu * geR) * b22tmmuLmmuL + 
     $     (geL * ssqthmu - csqthmu * geR) * b22tmmuRmmuR) + 
     $     ((geL * csqthe - ssqthe * geR) * b22tmeLmeL + 
     $     (geL * ssqthe - csqthe * geR) * b22tmeRmeR)))))

!------------------------

      result = (((smterm - higgsterm + charginoterm - 
     $     sfermionterm)*e*g)/(16.d0 * pi*pi*costhw))

      return
      end subroutine pizgamma

!================================================================================

