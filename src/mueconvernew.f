****f*SuSeFLAV/mueconvernew.f 
*  NAME
*     Subroutine muegamma
*  SYNOPSIS
*    . Calculates decay rates for rare lfv processes
*  FUNCTION
*     THIS PROGRAM CALCULATES THE DECAY RATES OF l_j -> l_i, Gamma;
*     l_j -> l_i,l_i,l_i  and Mu-E Conversion in Nuclei. 
*     THE AMPLITUDES AND DECAY RATES ARE CLOSELY FOLLOWED FROM THE
*     PAPERS : J. Hisano and Y. Nomura, PR D59 (1999) 116005 and 
*     J. Hisano et. al, Phys.Rev.D53: (1996) 2442 and our NOTES.
*  
*  INPUTS
*
*     Neg = Neutralino Mass Eigenvalues (4) (of dim 1)  ; 
*     ON = Diagonalisation Matrix of Neutralinos (4 \times 4) such that
*     ON*MNeut*Transpose[ON] = DiagonalMatrix[Neg]
*
*     Ceg = Chargino Mass Eigenvalues (2) (of dim 1) ; 
*     OCL = Left diagonalisation Matrix of Charginos (2 \times 2) ;
*     OCR = Right diagonalisation Matrix of Charginos (2 \times 2) ; 
*     suchthat OCR*MChar*Transpose[OCL] = DiagonalMatrix[Ceg] 
*
*     SNeg = Sneutrino Mass Eigenvalues (3) (of  dim 2) ; 
*     USN = diagonalisation Matrix of Sneutrinos  (3 \times 3); such that
*     USN*MSN*Transpose[USN] = DiagonalMatrix[SNeg]
*
*     SLeg = Slepton Mass Eigenvalues (6) (of dim 2) .  
*     USL = diagonalisation Matrix of Sleptons  (6 \times 6); such that
*     USL*MSL*Transpose[USL] = Diag[MSL]
*   
*  RESULT
*
*     ALN = Mu-E-Gamma left couplings amplitude with Neutralinos (LL + LR)
*     ARN = Mu-E-Gamma right couplings amplitude with Neutralinos (RL + RR)
*     ALC = Mu-E-Gamma left couplings amplitude with Charginos (LL + RL)
*     ARC = Mu-E-Gamma right couplings amplitude with Charginos (LR + LL)
*     megrate = Mu-E-Gamma decay rate
*     Bmeg = Branching ratio Mu-E-Gamma
*     TauALN = Tau-Mu-Gamma left couplings amplitude with Neutralinos (LL + LR)
*     TauARN = Tau-Mu-Gamma right couplings amplitude with Neutralinos (RL + RR)
*     TauALC = Tau-Mu-Gamma left couplings amplitude with Charginos (LL + RL)
*     TauARC = Tau-Mu-Gamma right couplings amplitude with Charginos (LR + LL)
*     Tmugrate = Tau-Mu-Gamma decay rate
*     Btmg = Branching ratio Tau-Mu-Gamma
*     TEALN = Tau-E-Gamma left couplings amplitude with Neutralinos (LL + LR)
*     TEARN = Tau-E-Gamma right couplings amplitude with Neutralinos (RL + RR)
*     TEALC = Tau-E-Gamma left couplings amplitude with Charginos (LL + RL)
*     TEARC = Tau-E-Gamma right couplings amplitude with Charginos (LR + LL)
*     Tegrate = Tau-E-Gamma decay rate
*     Bteg = Branching ratio Tau-E-Gamma
*
*  EXAMPLE
*     subroutine muegamma(tanbeta,mTau,alph2,ON,Neg,Ceg,OCL,OCR,
*     $     SLeg,USL,SNeg,USN,SUeg,USU,SDeg,USD,megrate,Bmeg,tmugrate,
*     $     Btmug,tegrate,Bteg,mu3erate,Brmu3E,tau3murate,brtau3mu,
*     $     tau3erate,brtau3e,mueconverrate,brmueconver,gm2)
*
*  NOTES
*
*	Suggestions : f1,f2,f3,f4 can be possibly redefined interms of 
*	              double precision functions  1/04/03. 	 
*       ALL Masses and Parameters are in GeV
*  BUGS
*    ---
*  SEE ALSO
*    -----
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C  
!=====================================================================

      subroutine muegamma(tanbeta,mTau,alph2,ON,Neg,Ceg,OCL,OCR,
     $     SLeg,USL,SNeg,USN,SUeg,USU,SDeg,USD,megrate,Bmeg,tmugrate,
     $     Btmug,tegrate,Bteg,mu3erate,Brmu3E,tau3murate,brtau3mu,
     $     tau3erate,brtau3e,mueconverrate,brmueconver,gm2)
      
      implicit none 

      integer a,x,b,t,j 

      double precision mTau
      double precision stw,ctw,beta,SUeg(6),tw,sw2
      double precision tanbeta,SLeg(6),USL(6,6)
      double precision SDeg(6),alph2
      double precision USN(3,3),SNeg(3), Neg(4), ON(4,4)
      double precision OCR(2,2), OCL(2,2), Ceg(2),g2,tegrate
      double precision CRlE(2,3),CRlMU(2,3),CRlTAU(2,3),Bteg
      double precision CLlE(2,3),CLlMU(2,3),CLlTAU(2,3),cdenc
      double precision g2t, NRlE(4,6),NRlMU(4,6),NRlTAU(4,6)
      double precision NLlE(4,6),NLlMU(4,6),NLlTAU(4,6),y(4,6)
      double precision z(2,3),f1(4,6),f2(4,6)
      double precision AmpARN(4,6),ARN,f3(2,3),f4(2,3),AmpLC(2,3)
      double precision AmpRC(2,3),ALC,ARC,Bmeg,megrate,AmpALN(4,6)
      double precision AmpTauALN(4,6),AmpTauARN(4,6),TauALN, TauARN
      double precision AmpTauALC(2,3),TauALC,AmpTauARC(2,3),TauARC
      double precision tmugrate,Btmug,cdenn,piconst,ALN
      double precision AmpTEALN(4,6),TEALN,AmpTEARN(4,6),TEARN
      double precision AmpTEALC(2,3),TEALC,AmpTEARC(2,3),TEARC
      double precision fzn(6,4,4),gzn(6,4,4), fzc(3,2,2),gzc(3,2,2)
      double precision ZpengTauMuMuMuC, ZpengMuEEEC,ZpengTauEEEC
      double precision ZpengTauMueeC,fa1(4,6),fa2(2,3),VAN1LMUE
      double precision VAN1LTAUMU,VAN1LTAUE,VAC1LMUE,VAC1LTAUMU
      double precision VAC1LTAUE,ZpengMUEEENL,ZpengTauMuEENL 
      double precision ZpengTauMUMUMUNL,ZpengTauEEENL,ZpengMuEEENR
      double precision ZpengTauMuEENR,ZpengTauMuMuMuNR,ZpengTauEEENR
      double precision VAN1RMUE, VAN1RTAUMU,VAN1RTAUE,VAC1RMUE
      double precision VAC1RTAUMU,VAC1RTAUE,J4n(4,4,6,6),J4c(2,2,3,3)
      double precision I4n(4,4,6,6),I4c(2,2,3,3),I3n(4,4,6),I3c(2,2,3)
      double precision B1nLMu3E, B2nLMu3E,B3nLMu3E,B4nLMu3E,B1nRMu3E
      double precision B2nRMu3E,B3nRMu3E,B4nRMu3E,B1cLMu3E,B2cLMu3E
      double precision B3cLMu3E,B4cLMu3E,B1cRMu3E,B2cRMu3E,B3cRMu3E
      double precision B4cRMu3E,B1nLTau3E,B2nLTau3E,B3nLTau3E,B4nLTau3E
      double precision B2nRTau3E,B3nRTau3E,B4nRTau3E,B1cLTau3E,B2cLTau3E
      double precision B3cLTau3E,B4cLTau3E,B1cRTau3E,B2cRTau3E,B3cRTau3E
      double precision B4cRTau3E,B1nRTau3E,B1nLTau3Mu,B2nLTau3Mu
      double precision B3nLTau3Mu,B4nLTau3Mu, B2nRTau3Mu,B3nRTau3Mu
      double precision B4nRTau3Mu,B1cLTau3Mu,B2cLTau3Mu,B3cLTau3Mu
      double precision B4cLTau3Mu,B1cRTau3Mu,B2cRTau3Mu,B3cRTau3Mu
      double precision B4cRTau3Mu,B1nRTau3Mu,FLLMu3E,FRRMu3E,FLRMu3E
      double precision FRLMu3E,B1LMu3E,B1RMu3E,B2LMu3E,B2RMu3E,B3LMu3E
      double precision B3RMu3E,B4LMu3E,B4RMu3E,ALMuE,ARMuE,VLMuE
      double precision VRMuE,Mu3ERate,Brmu3E,b1ltau3mu,b2ltau3mu
      double precision b3ltau3mu,b4ltau3mu,b1rtau3mu,b2rtau3mu,b3rtau3mu
      double precision b4rtau3mu,b1ltau3e,b2ltau3e,b3ltau3e,b4ltau3e
      double precision b1rtau3e,b2rtau3e,b3rtau3e,b4rtau3e,flltau3e
      double precision flrtau3e,frltau3e,frrtau3e,flltau3mu,flrtau3mu
      double precision frltau3mu,frrtau3mu,vltaumu,vrtaumu,altaumu
      double precision artaumu,vltaue,vrtaue,altaue,artaue,tau3murate
      double precision tau3erate,brtau3mu,brtau3e,Zti,Nti,Zeffti,Fqsqti 
      double precision Brmueconver,mueconverrate,DbarLdMuE,DbarLuMuE
      double precision DcLuMuE, DcLdMuE,DnRdMuE,DnRuMuE,DnLdMuE,DnLuMuE
      double precision USU(6,6),USD(6,6),CRdD(2,6),CLdD(2,6),CRuU(2,6)
      double precision CLuU(2,6),NRdD(4,6),NLdD(4,6),NRuU(4,6),NLuU(4,6)
      double precision ScNeg(4),ScCeg(2),ScSLeg(6),ScSNeg(3),I3f2d
      double precision ScI3n(4,4,6),ScI3c(2,2,3),ScI4n(4,4,6,6)
      double precision scI4c(2,2,3,3),diff1,diff2,I4f3,I4f3d1,I4f2d2d
      double precision qqI4cd(2,2,3,6),qqI4nd(4,4,6,6)
      double precision qqJ4nd(4,4,6,6),qqJ4cd(2,2,3,6)
      double precision qqI4cu(2,2,3,6),qqJ4cu(2,2,3,6)
      double precision ScSUeg(6),ScSDeg(6),qqI4nu(4,4,6,6)
      double precision qqJ4nu(4,4,6,6)
      double precision ScqqI4nu(4,4,6,6)
      double precision ScqqI4nd(4,4,6,6),ScqqI4cu(2,2,3,6),DbarRdMuE
      double precision ScqqI4cd(2,2,3,6),DcRdMuE,DcRuMuE,DbarRuMuE
      double precision Ampg2LN(4,6),g2ALN,Ampg2RN(4,6),g2ARN
      double precision Ampg2LC(2,3),g2ALC,Ampg2RC(2,3),g2ARC,gminus2
      double precision gmC2, gmN2, gm2, gm2C(2,3), gm2N(4,6)

      double precision negarray(4),o1,o2,o3,o4,ONtemp(4,4),ONew(4,4)
      double precision alph,Gf,alphas,MW,MZ,pi,sinsqtw
      double precision mbpole, mtaupole, Mtpole, MZpole,MWpole

      integer i,o1pos,o2pos,o3pos,o4pos
      double precision m0, m12, m10, m20, sgnmu, a0
      common/mssminputs/ m0, m12, m10, m20, sgnmu,  a0
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole
      common/gauge/alph,Gf,alphas
 
!----------------------------------------------------------------------------

      include 'stdinputs.h'
!--------    

      MW = MWpole
      MZ = MZpole

      pi = 4.d0 * datan(1.d0)

      sw2 = (1.d0 - (MW / MZ)**2.d0)
      stw = dsqrt(sw2)
      ctw = dsqrt(1.d0 - sw2)
      tw = dasin(dsqrt(sw2))
      sinsqtw = sw2
      
      beta = datan(tanbeta)


      cdenc = dsqrt(2.d0)*MW*cos(beta)
      cdenn = MW*cos(beta)
      g2 = dsqrt(16.d0 * pi*pi * alph2)

      
C     Chargino-Fermion-sFermion Couplings 
C     ***********************************

      loopi01: do a = 1, 2 
      loopi02: do x = 1, 3 
      
      CRlE(a,x)   = 0.d0
      CRlMU(a,x)  = 0.d0
      CRlTAU(a,x) = 0.d0
      CLlE(a,x)   = 0.d0
      CLlMU(a,x)  = 0.d0
      CLlTAU(a,x) = 0.d0

      enddo loopi02
      enddo loopi01
      
C     --------------------------------------------------

      loop101: do a = 1, 2 
      loop102: do x = 1, 3 
C     -------------------------------------------------
C     LL-TYPE
C     -------------------------------------------------

      CRlE(a,x)   =  - (g2*OCR(a,1)*USN(x,1))

      CRlMU(a,x)  =  - (g2*OCR(a,1)*USN(x,2))

      CRlTAU(a,x) =  - (g2*OCR(a,1)*USN(x,3))

C     --------------------------------------------------
C     LR-TYPE
C     -------------------------------------------------

      CLlE(a,x)   = g2*mE/cdenc*OCL(a,2)*USN(x,1)
      
      CLlMU(a,x)  = g2*mMU/cdenc*OCL(a,2)*USN(x,2)
      
      CLlTAU(a,x) = g2*mTAU/cdenc*OCL(a,2)*USN(x,3)
      
      enddo loop102
      enddo loop101
      
C     Neutralino-Fermion-sFermion Couplings 
C     *************************************
      
      g2t = g2/dsqrt(2.d0)
      
C     --------------------------------------------------
      loopi03: do a = 1, 4
      loopi04: do x = 1, 6
      
      NRlE(a,x)   = 0.d0
      NRlMU(a,x)  = 0.d0
      NRlTAU(a,x) = 0.d0
      NLlE(a,x)   = 0.d0
      NLlMU(a,x)  = 0.d0
      NLlTAU(a,x) = 0.d0
      
      enddo loopi04
      enddo loopi03 
C     --------------------------------------------------
      
      loop103: do a = 1, 4 
      loop104: do x = 1, 6

C     --------------------------------------------------
C     LL + RL TYPE MI
C     --------------------------------------------------
      
      NRlE(a,x) = - (g2t)*((-ON(a,2) - ON(a,1)*dtan(tw))*USL(x,1) 
     $     + (mE/cdenn)*ON(a,3)*USL(x,4))
      
      NRlMU(a,x) = -(g2t)*((-ON(a,2) - ON(a,1)*dtan(tw))*USL(x,2) 
     $     + (mMU/cdenn)*ON(a,3)*USL(x,5))
      
      NRlTAU(a,x) = -(g2t)*((-ON(a,2) - ON(a,1)*dtan(tw))*USL(x,3) 
     $     + (mTAU/cdenn)*ON(a,3)*USL(x,6))
      
C     ---------------------------------------------------------
C     RL + RR TYPE MI 
C     ---------------------------------------------------------
      
      NLlE(a,x) =  -(g2t)*((mE/cdenn)*ON(a,3)*USL(x,1) 
     $     + 2.d0*ON(a,1)*dtan(tw)*USL(x,4))
      
      NLlMU(a,x) =  -(g2t)*((mMU/cdenn)*ON(a,3)*USL(x,2) 
     $     + 2.d0*ON(a,1)*dtan(tw)*USL(x,5))
      
      NLlTAU(a,x) =  -(g2t)*((mTAU/cdenn)*ON(a,3)*USL(x,3) 
     $     + 2.d0*ON(a,1)*dtan(tw)*USL(x,6))
      


      enddo loop104
      enddo loop103


C\\\\\\Checked @ 14:27 on 14/05/2010\\\\\\\\\\\\\\\\\\\\\\
      
C     ---------------------------------------------------------
C     Mu--->E-Gamma Amplitudes : 
C     ---------------------------------------------------------
      
C     ------------------------
C     Dipole Neutralino Contributions 
C     ------------------------
      
      piconst = 1.d0/(32.d0 *pi*pi) 
      
C     Defintion of y - remember it is a dimenionless quantity. 
C     ---------------------------------------------------------
      loopi05: do a = 1, 4
      loopi06: do x = 1, 6
      y(a,x) = 0.d0
      enddo loopi06
      enddo loopi05
C     ---------------------------------------------------------
      
      loop105: do a = 1, 4
      loop106: do x = 1, 6 
      
      y(a,x) = (Neg(a)**2.d0)/SLeg(x)
      
      enddo loop106 
      enddo loop105
      
C     f1(x) in the notes
C     ------------------
      
      loopi07: do a = 1, 4 
      loopi08: do x = 1, 6 
      f1(a,x) = 0.d0
      enddo loopi08
      enddo loopi07 

C     ----------------------------------------------------
      
      loop107: do a = 1, 4 
      loop108: do x = 1, 6 
      
      cond901: if(dabs(1.d0 - y(a,x)).gt.0.01d0)then
         
         f1(a,x) = ((1.d0 - 6.d0*y(a,x) + 3.d0*(y(a,x)**2.d0) + 
     $        2.d0*(y(a,x)**3.d0) - 6.d0*(y(a,x)**2.d0)*dlog(y(a,x)))/
     $        (6.d0 * (1.d0 - y(a,x))**4.d0))
         
      else cond901 
         
         f1(a,x) = 1.d0/12.d0 - (y(a,x) - 1.d0)/30.d0
         
      endif cond901 
      
      enddo loop108
      enddo loop107 

      
C     f2(x) in the notes
C     ------------------
      loopi09: do a = 1,4
      loopi10: do x = 1,6 
      f2(a,x)  = 0.d0
      enddo loopi10
      enddo loopi09 
C     ------------------
      
      loop109: do a = 1,4
      loop110: do x = 1,6 
      
      cond902:  if(dabs(1.d0 - y(a,x)).gt.0.01d0)then
         
         f2(a,x)  = ((1.d0 - y(a,x)**2.d0 +  2.d0 * y(a,x) *
     $        dLog(y(a,x)))/((1.d0-y(a,x))**3.d0))
         
      else cond902 
         
         f2(a,x) = 1.d0/3.d0 - (y(a,x) - 1.d0)/6.d0
         
      endif cond902 
      
      enddo loop110
      enddo loop109
      
C     Amplitude (L)
C     ----------------
      loopi11: do a = 1, 4 
      loopi12: do x = 1, 6
      
      AmpALN(a,x) = 0.d0
      
      enddo loopi12
      enddo loopi11
      
      loop111: do a = 1, 4 
      loop112: do x = 1, 6
      
      AmpALN(a,x) =  (NLlE(a,x) * NLlMU(a,x) * f1(a,x) + 
     $     NLlE(a,x) * NRlMU(a,x) * (Neg(a)/mMU) * f2(a,x)) * 
     $     (1.d0/SLeg(x))
      
      enddo loop112
      enddo loop111
      
      ALN = 0.d0 

C     Sum Neutralino (L)
C     ------------------
      
      loop113: do a = 1, 4 
      loop114: do x = 1, 6 

      ALN = ALN + AmpALN(a,x)

      enddo loop114
      enddo loop113
      
      ALN = piconst*ALN
      
C     Amplitude (R)
C     ------------------
      loopi15: do a = 1, 4
      loopi16: do x = 1, 6
      AmpARN(a,x) = 0.d0
      enddo loopi16
      enddo loopi15
      
      
      loop115: do a = 1, 4
      loop116: do x = 1, 6
      
      AmpARN(a,x) = (NRlE(a,x) * NRlMU(a,x) * f1(a,x) + 
     $     NRlE(a,x) * NLlMU(a,x) * (Neg(a)/mMU) * f2(a,x)) * 
     $     (1.d0/SLeg(x))
      
      enddo loop116
      enddo loop115
      
      
      ARN = 0.d0 

C     Sum Neutralino (R)
C     ------------------

      loop117: do a = 1, 4 
      loop118: do x = 1, 6 

      ARN = ARN + AmpARN(a,x)

      enddo loop118 
      enddo loop117
      
      ARN = piconst*ARN

C     ------------------------------------------------------
C     Dipole Chargino Contributions  
C     ------------------------------------------------------

      loopi19:do a = 1, 2
      loopi20:do x = 1, 3
      z(a,x) = 0.d0
      enddo loopi20
      enddo loopi19

      loop119: do a = 1, 2
      loop120: do x = 1, 3
      
      z(a,x) = (Ceg(a)**2.d0)/SNeg(x)
      
      enddo loop120
      enddo loop119

C     function f3 in notes
C     --------------------

      loopi21: do a = 1,2
      loopi22: do x = 1,3
      f3(a,x) = 0.d0
      enddo loopi22
      enddo loopi21

C     --------------------

      loop121: do a = 1,2
      loop122: do x = 1,3
      
      cond903: if(dabs(1.d0-z(a,x)).gt.0.01d0)then
         
         f3(a,x) = ((2.d0 + 3.d0*z(a,x) - 6.d0*(z(a,x)**2.d0) +  
     $        z(a,x)**3.d0 + 6.d0*z(a,x)*dlog(z(a,x)))/
     $        (6.d0*(1.d0 - z(a,x))**4.d0)) 
         
      else cond903

         f3(a,x) = 1.d0/12.d0 - (z(a,x) - 1.d0)/20.d0

      endif cond903

      enddo loop122
      enddo loop121
      
C     function f4 in notes
C     --------------------

      loopi23:do a = 1,2
      loopi24:do x = 1,3
      f4(a,x) = 0.d0
      enddo loopi24
      enddo loopi23

C     --------------------

      loop123: do a = 1,2
      loop124: do x = 1,3
      
      cond904: if(dabs(1.d0-z(a,x)).gt.0.01d0)then 
         
         f4(a,x) = ((-3.d0 + 4.d0*z(a,x) - z(a,x)**2.d0 - 
     $        2.d0*dLog(z(a,x)))/((1.d0 - z(a,x))**3.d0))
         
      else cond904
         
         f4(a,x) = 2.d0/3.d0 - (z(a,x) - 1.d0)/2.d0
         
      endif cond904
      
      enddo loop124
      enddo loop123
      
C     Chargino Contribution (L)
C     -------------------------

      loopi25: do a = 1,2
      loopi26: do x = 1, 3
      AmpLC(a,x) = 0.d0
      enddo loopi26
      enddo loopi25
C     -------------------------

      loop125: do a = 1, 2
      loop126: do x = 1, 3

      AmpLC(a,x) = -(piconst) * (1.d0/SNeg(x)) * (CLlE(a,x)*CLlMU(a,x) * 
     $     f3(a,x) + CLlE(a,x) * CRlMU(a,x) * (Ceg(a)/mMU) * f4(a,x))

      enddo loop126
      enddo loop125

      
      ALC = 0.d0
      
C     Sum Chargino (L)
C     ----------------

      loop127: do a = 1, 2
      loop128:   do x = 1, 3 
      
      ALC = ALC + AmpLC(a,x)
      
      enddo loop128
      enddo loop127


C     Chargino Contribution (R)
C     -------------------------
      loopi29:do a = 1,2
      loopi30:do x = 1, 3

      AmpRC(a,x) = 0.d0

      enddo loopi30
      enddo loopi29
C     -------------------------
      loop129: do a = 1,2
      loop130: do x = 1, 3

      AmpRC(a,x) = -(piconst) * (1.d0/SNeg(x)) * (CRlE(a,x)*CRlMU(a,x) * 
     $     f3(a,x) + CRlE(a,x) * CLlMU(a,x) * (Ceg(a)/mMU) * f4(a,x))

      enddo loop130
      enddo loop129

      ARC = 0.d0

C     Sum Chargino (R)
C     ----------------

      loop131: do a = 1, 2
      loop132: do x = 1, 3 
      
      ARC = ARC + AmpRC(a,x)
      
      enddo loop132
      enddo loop131
      
C     ------------------------------------------------
C     Mu-E-Gamma Decay Rate
C     ------------------------------------------------
      megrate = 0.d0
      Bmeg    = 0.d0

C     ------------------------------------------------
      
      megrate = (alph/4.d0) * (mMU**5.d0) * ((ALN + ALC )**(2.d0)
     $     + (ARN + ARC)**(2.d0))
      
      Bmeg = megrate/(3.d0 * (10.d0**(-19.d0)))


C     ================================================================
C     Tau-Mu-Gamma	        
C     ================================================================


C     Neutralino Amplitude (L)
C     ------------------------
      
      loopi33: do a = 1, 4
      loopi34: do x = 1, 6
      AmpTauALN(a,x) = 0.d0
      enddo loopi34
      enddo loopi33
      
C     -----------------------------------------------
      
      loop133: do a = 1, 4
      loop134: do x = 1, 6
      
      AmpTauALN(a,x) = piconst * (NLlMU(a,x) * NLlTAU(a,x) * f1(a,x) +
     $     NLlMU(a,x) * NRlTAU(a,x) * (Neg(a)/mTAU) * f2(a,x)) * 
     $     (1.d0/SLeg(x))
      
      enddo loop134
      enddo loop133
      
      
C     Sum Neutralino Amplitude (L)
C     ---------------------------

      TauALN = 0.d0
      
      loop135: do a = 1, 4
      loop136: do x = 1, 6

      TauALN = TauALN + AmpTauALN(a,x)
      
      enddo loop136
      enddo loop135

C     Neutralino Amplitude (R) 
C     ------------------------
      loopi37: do a = 1, 4
      loopi38:  do x = 1, 6
      AmpTauARN(a,x) = 0.d0
      enddo loopi38             
      enddo loopi37
C     ------------------------

      loop137: do a = 1, 4
      loop138: do x = 1, 6

      AmpTauARN(a,x) = piconst * (NRlMU(a,x) * NRlTAU(a,x) * f1(a,x) +
     $     NRlMU(a,x)*NLlTAU(a,x)*(Neg(a)/mTAU) * f2(a,x)) *
     $     (1.d0/SLeg(x))

      enddo loop138
      enddo loop137

C     Sum Neutralinos (R)
C     ---------------------

      TauARN = 0.d0
      
      loop139: do a = 1, 4
      loop140: do x = 1, 6

      TauARN = TauARN + AmpTauARN(a,x)
      
      enddo loop140
      enddo loop139


C     Chargino Amplitudes (L)
C     -------------------------------------

      loopi41: do a = 1,2
      loopi42: do x = 1, 3 
      AmpTauALC(a,x) = 0.d0
      enddo loopi42
      enddo loopi41

C     -------------------------------------

      loop141: do a = 1,2
      loop142: do x = 1, 3 

      AmpTauALC(a,x) = -piconst * (1/SNeg(x)) * (CLlMU(a,x)*CLlTAU(a,x)*
     $     f3(a,x) + CLlMU(a,x) * CRlTAU(a,x) * (Ceg(a)/mTAU) * f4(a,x))

      enddo loop142
      enddo loop141

      TauALC = 0.d0

C     Sum Chargino Amplitude (L)
C     --------------------------

      loop143: do a = 1,2
      loop144: do x = 1, 3
      
      TauALC = TauALC + AmpTauALC(a,x)

      enddo loop144
      enddo loop143


C     Chargino Amplitudes (R)
C     ---------------------------------------------
      
      loopi45: do a = 1,2
      loopi46: do x = 1,3 

      AmpTauARC(a,x) = 0.d0

      enddo loopi46
      enddo loopi45

C     ---------------------------------------------

      loop145: do a = 1,2
      loop146: do x = 1,3 

      AmpTauARC(a,x) = -piconst * (1/SNeg(x)) * (CRlMU(a,x)*CRlTAU(a,x)*
     $     f3(a,x) + CRlMU(a,x)*CLlTAU(a,x) * (Ceg(a)/mTAU) * f4(a,x))

      enddo loop146
      enddo loop145

      TauARC = 0.d0

C     Sum Chargino Amplitude (R)
C     ---------------------------

      loop147: do a = 1,2
      loop148: do x = 1, 3
      
      TauARC = TauARC + AmpTauARC(a,x)

      enddo loop148
      enddo loop147

      
C     Tau-Mu-Gamma Decay Rate 
C     ------------------------------------------------
      tmugrate = 0.d0
      Btmug = 0.d0
C     ------------------------

      tmugrate = (alph/4.d0) * (mTAU**5.d0) * ((TauALC + TauALN)**2.d0 + 
     $     (TauARC + TauARN)**2.d0)


      Btmug = tmugrate/(2.26d0*10.d0**(-12.d0))

C     ================================================================
C     Tau-E-Gamma	        
C     ================================================================
      
C     Neutralino Amplitude (L)
C     ------------------------

      loopi49: do a = 1, 4
      loopi50: do x = 1, 6

      AmpTEALN(a,x) = 0.d0

      enddo loopi50
      enddo loopi49

C     ------------------------

      loop149: do a = 1, 4
      loop150: do x = 1, 6

      AmpTEALN(a,x) = piconst * (NLlE(a,x)*NLlTAU(a,x) * f1(a,x) + 
     $     NLlE(a,x) * NRlTAU(a,x) * (Neg(a)/mTAU) * f2(a,x)) * 
     $     (1.d0/SLeg(x))

      enddo loop150
      enddo loop149

C     Sum Neutralino Amplitude (L)
C     ---------------------------

      TEALN = 0.d0
      
      loop151: do a = 1, 4
      loop152: do x = 1, 6

      TEALN = TEALN + AmpTEALN(a,x)
      
      enddo loop152
      enddo loop151

C     Neutralino Amplitude (R) 
C     ------------------------------------------------------

      loopi53: do a = 1, 4
      loopi54: do x = 1, 6

      AmpTEARN(a,x) = 0.d0

      enddo loopi54
      enddo loopi53

C     ------------------------------------------------------
      
      loop153: do a = 1, 4
      loop154: do x = 1, 6
      
      AmpTEARN(a,x) = piconst * (NRlE(a,x)*NRlTAU(a,x) * f1(a,x) +
     $     NRlE(a,x)*NLlTAU(a,x)*(Neg(a)/mTAU) * f2(a,x)) * 
     $     (1.d0/SLeg(x))
      
      enddo loop154
      enddo loop153 
      
      
C     Sum Neutralinos (R)
C     ---------------------

      TEARN = 0.d0
      
      loop155: do a = 1, 4
      loop156: do x = 1, 6

      TEARN = TEARN + AmpTEARN(a,x)
      
      enddo loop156
      enddo loop155


C     Chargino Amplitudes (L)
C     -------------------------------------

      loopi57: do a = 1,2
      loopi58: do x = 1, 3 

      AmpTEALC(a,x) = 0.d0

      enddo loopi58
      enddo loopi57
      
C     -------------------------------------
      
      loop157: do a = 1, 2
      loop158: do x = 1, 3 
      
      AmpTEALC(a,x)= -piconst * (1/SNeg(x)) * (CLlE(a,x)*CLlTAU(a,x) *
     $     f3(a,x) + CLlE(a,x)*CRlTAU(a,x) * (Ceg(a)/mTAU) * f4(a,x))
      
      enddo loop158
      enddo loop157
      
      TEALC = 0.d0
      
C     Sum Chargino Amplitude (L)
C     --------------------------

      loop159: do a = 1,2
      loop160: do x = 1,3
      
      TEALC = TEALC + AmpTEALC(a,x)
      enddo loop160
      enddo loop159
      


C     Chargino Amplitudes (R)
C     ---------------------------------------------

      loopi61: do a = 1, 2
      loopi62: do x = 1, 3 

      AmpTEARC(a,x) = 0.d0

      enddo loopi62
      enddo loopi61

C     ---------------------------------------------
      
      loop161: do a = 1, 2
      loop162: do x = 1, 3 


      AmpTEARC(a,x) = -piconst * (1/SNeg(x)) * (CRlE(a,x)*CRlTAU(a,x) *
     $     f3(a,x) + CRlE(a,x)*CLlTAU(a,x) *(Ceg(a)/mTAU) * f4(a,x))

      enddo loop162
      enddo loop161

      TEARC = 0.d0


C     Sum Chargino Amplitude (R)
C     ---------------------------

      loop163: do a = 1,2
      loop164: do x = 1, 3
      
      TEARC = TEARC + AmpTEARC(a,x)

      enddo loop164
      enddo loop163
      
C     ------------------------
C     Tau-E-Gamma Decay Rate 
C     ------------------------
      tegrate = 0.d0
      Bteg = 0.d0
C     ------------------------

      tegrate = (alph/4.d0) * (mTAU**5.d0) * ((TEALC + TEALN)**2.d0 + 
     $     (TEARC + TEARN)**2.d0)

      Bteg = tegrate/(2.26d0*10.d0**(-12.d0))

C\\\\\\\\\\\\Checked @ 11:30 on 16/05/2010\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


C     =====================  CHARGINO CONTRIBUTION  ======================
      
      loopgm2cx: do x = 1,3
      loopgm2ca: do a = 1,2

c$$$      gm2C(a,x) = (4.d0 * piconst * (mmu*mmu/SNeg(x)) * (CLlMU(a,x)*
c$$$     $     CLlMU(a,x)) * f3(a,x) + 2.d0 * piconst * (mmu*Ceg(a)/SNeg(x))
c$$$     $     * CLlMU(a,x)*CRlMU(a,x) * f4(a,x))

c$$$      gm2C(a,x) = (2.d0 * piconst * (mmu*mmu/SNeg(x)) * 
c$$$     $     (CLlMU(a,x)*CLlMU(a,x) + CRlMU(a,x)*CRlMU(a,x)) * f3(a,x) + 
c$$$     $     2.d0 * piconst * (mmu*Ceg(a)/SNeg(x)) * 
c$$$     $     (CLlMU(a,x)*CRlMU(a,x)) * f4(a,x))

      gm2C(a,x) = ((2.d0/3.d0) * piconst * (mmu*mmu/SNeg(x)) * 
     $     (CLlMU(a,x)*CLlMU(a,x) + CRlMU(a,x)*CRlMU(a,x)) * f3(a,x) + 
     $     2.d0 * piconst * (mmu*Ceg(a)/SNeg(x)) * 
     $     (CLlMU(a,x)*CRlMU(a,x)) * f4(a,x)) 
      
      enddo loopgm2ca 
      enddo loopgm2cx
      
      gmC2 = 0.d0
      
      loopcx: do x = 1, 3
      loopca: do a = 1, 2
      gmC2 = gmC2 + gm2C(a,x)
      enddo loopca
      enddo loopcx
      
C     ===================== NEUTRALINO CONTRIBUTION ====================
      
      loopgm2nx: do x = 1, 6
      loopgm2na: do a = 1, 4

c$$$      gm2N(a,x) = - (4.d0 * piconst * (mmu*mmu/SLeg(x)) * NLlMU(a,x)*
c$$$     $     NLlMU(a,x) * f1(a,x) + 2.d0 * piconst * (mmu*Neg(a)/SLeg(x))
c$$$     $     *NLlMU(a,x)*NRlMU(a,x) * f2(a,x)) 

c$$$      gm2N(a,x) = - (2.d0 * piconst * (mmu*mmu/SLeg(x)) * 
c$$$     $     (NLlMU(a,x)*NLlMU(a,x) + NRlMU(a,x)*NRlMU(a,x)) * f1(a,x) + 
c$$$     $     2.d0 * piconst * (mmu*Neg(a)/SLeg(x)) * 
c$$$     $     (NLlMU(a,x)*NRlMU(a,x)) * f2(a,x)) 


      gm2N(a,x) = - ((2.d0/3.d0) * piconst * (mmu*mmu/SLeg(x)) * 
     $     (NLlMU(a,x)*NLlMU(a,x) + NRlMU(a,x)*NRlMU(a,x)) * f1(a,x) + 
     $     2.d0 * piconst * (mmu*Neg(a)/SLeg(x)) * 
     $     (NLlMU(a,x)*NRlMU(a,x)) * f2(a,x)) 
      
      enddo loopgm2na
      enddo loopgm2nx
      
      gmN2 = 0.d0
      
      loopnx: do x = 1,6
      loopna: do a = 1,4
      gmN2 = gmN2 + gm2N(a,x)
      enddo loopna
      enddo loopnx

      gm2 = (gmN2 + gmC2)


C-----------------------------------------------------------------------

!      goto 456

C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\      
C     =============================================
C     Vector (Monopole)  Contributions (Neutralino) 
C     =============================================
C     (Left)
C     ==========

      
        loop165:Do a = 1,4
        loop166:   Do x = 1, 6 

	fa1(a,x) = 0.d0
	
        cond905: if(dabs(1.d0-y(a,x)).gt.0.01d0)then

            fa1(a,x) = (1.d0/SLeg(x))*(1.d0/(1.d0 - y(a,x))**4.d0)*
     . (2.d0 - 9.d0*y(a,x) + 18.d0*y(a,x)**2.0 - 11.d0*y(a,x)**3.d0
     . + 6.d0*y(a,x)**3.d0*dlog(y(a,x)))
            
            else cond905
               fa1(a,x) = (3.d0/2.d0)*(1.d0/SLeg(x))
            
               endif cond905

            enddo loop166
            enddo loop165

	VAN1LMUE = 0.d0
	VAN1LTAUMU = 0.d0
	VAN1LTAUE = 0.d0

         loop167:   do a = 1, 4
         loop168:   do x = 1, 6 

        VAN1LMUE = VAN1LMUE + (piconst/18.d0)*NRlE(a,x)*NRlMU(a,x)
     .  *fa1(a,x) 

        VAN1LTAUMU = VAN1LTAUMU +(piconst/18.d0)*NRlMU(a,x)*NRlTAU(a,x)
     .  *fa1(a,x) 

        VAN1LTAUE = VAN1LTAUE +(piconst/18.d0)*NRlE(a,x)*NRlTAU(a,x)
     .  *fa1(a,x) 

        enddo loop168
        enddo loop167

C	(Right)
C	========

	VAN1RMUE = 0.d0
	VAN1RTAUMU = 0.d0
	VAN1RTAUE = 0.d0

	loop169: Do a = 1,4
	loop170: Do x = 1, 6

        VAN1RMUE = VAN1RMUE + (piconst/18.d0)*NLlE(a,x)*NLlMU(a,x)
     .  *fa1(a,x) 

        VAN1RTAUMU = VAN1RTAUMU +(piconst/18.d0)*NLlMU(a,x)*NLlTAU(a,x)
     .  *fa1(a,x) 

        VAN1RTAUE = VAN1RTAUE +(piconst/18.d0)*NLlE(a,x)*NLlTAU(a,x)
     .  *fa1(a,x) 

        enddo loop170
        enddo loop169

C     ==================================
C     Vector (Monopole) Contributions (Chargino) 
C     ==================================

	
        loop171: Do a = 1,2
        loop172:   Do x = 1, 3 

	fa2(a,x) = 0.d0

       cond906: if((1.d0-z(a,x)).gt.0.01d0)then

         fa2(a,x) = (1.d0/SNeg(x))*(1.d0/(1.d0 - z(a,x))**4.d0)*
     . (16.d0 - 45.d0*z(a,x) + 36.d0*z(a,x)**2.0 - 7.d0*z(a,x)**3.d0
     . + 6.d0*(2.d0 - 3.d0*z(a,x))*dlog(z(a,x)))
            
         else cond906
            fa2(a,x) = (-9.d0/2.d0)*(1.d0/SNeg(x))

            endif cond906

            enddo loop172
            enddo loop171


C	 LEFT	
C	=======	

	VAC1LMUE = 0.d0
	VAC1LTAUMU = 0.d0
	VAC1LTAUE = 0.d0

          loop173:  do a = 1, 2
          loop174:  do x = 1, 3 

        VAC1LMUE = VAC1LMUE - (piconst/18.d0)*CRlE(a,x)*CRlMu(a,x)
     .  *fa2(a,x) 

        VAC1LTAUMU = VAC1LTAUMU - (piconst/18.d0)*CRlMU(a,x)*
     .	CRlTAU(a,x)*fa2(a,x) 

        VAC1LTAUE = VAC1LTAUE - (piconst/18.d0)*CRlE(a,x)*CRlTAU(a,x)
     .  *fa2(a,x) 

        enddo loop174
        enddo loop173

C	 RIGHT 
C	=======	

	VAC1RMUE = 0.d0
	VAC1RTAUMU = 0.d0
	VAC1RTAUE = 0.d0

         loop175: do a = 1, 2
         loop176: do x = 1, 3 

        VAC1RMUE = VAC1RMUE - (piconst/18.d0)*CLlE(a,x)*CLlMu(a,x)
     .  *fa2(a,x) 

        VAC1RTAUMU = VAC1RTAUMU - (piconst/18.d0)*CLlMU(a,x)*
     .	CLlTAU(a,x)*fa2(a,x) 

        VAC1RTAUE = VAC1RTAUE - (piconst/18.d0)*CLlE(a,x)*CLlTAU(a,x)
     .  *fa2(a,x) 

        enddo loop176
        enddo loop175


C     ===========================================================================
C     Loop functions for neutralino Z- Penguins 
C     ==========================================================================
C     Note: Here contributions proportional to the Yukawa couplings are set to 
C     zero. In particular this leads to Right Chargino Contributions to be zero. 
C     ===========================================================================

	loop177: Do x = 1,6 
	loop178: Do a = 1, 4	
	loop179: Do b = 1, 4

	fzn(x,a,b) = 0.d0
	gzn(x,a,b) = 0.d0
	
	cond907: If(a.ne.b.and.dabs(y(a,x)-y(b,x)).gt.0.01d0.and.
     .   dabs(1.d0-y(a,x)).gt.0.01d0.and.
     .   dabs(1.d0-y(b,x)).gt.0.01d0)then

	fzn(x,a,b) = dLog(y(a,x)) + 1.d0/(y(a,x) - y(b,x))*
     . 	( (y(a,x)**2.d0*dLog(y(a,x)))/(1.d0 - y(a,x)) -
     . 	 (y(b,x)**2.d0*dLog(y(b,x)))/(1.d0 - y(b,x)) )

	gzn(x,a,b) = (Neg(a)*Neg(b)/SLeg(x))*(1.d0/(y(a,x)- y(b,x)))*
     . 	( (y(a,x)*dLog(y(a,x)))/(1.d0 - y(a,x)) -
     . 	 (y(b,x)*dLog(y(b,x)))/(1.d0 - y(b,x)) )
	
         elseif(a.ne.b.and.dabs(y(a,x)-y(b,x)).gt.0.01d0.and.
     .   dabs(1.d0-y(a,x)).gt.0.01d0.and.dabs(1.d0-y(b,x)).le.
     .   0.01d0)then cond907

        fzn(x,a,b) = (-2.d0*dLog(y(a,x))*y(a,x) + dlog(y(a,x)) 
     .   + y(a,x) - 1.d0)/(y(a,x) -1.d0)**2.d0 

        gzn(x,a,b) = -(Neg(a)*Neg(b)/Sleg(x))*(dLog(y(a,x))*y(a,x)
     .   - y(a,x) + 1.d0)/((y(a,x) - 1.d0)**2.d0)

        elseif(a.ne.b.and.dabs(y(a,x)-y(b,x)).gt.0.01d0.and.
     .  dabs(1.d0-y(b,x)).gt.0.01d0.and.dabs(1.d0-y(a,x)).
     .  le.0.01d0)then cond907

         fzn(x,a,b) = (- dLog(y(b,x))*y(b,x)**2.d0 + y(b,x) - 1.d0)/
     .   (y(b,x) - 1.d0)**2.d0 
              
         gzn(x,a,b) = - (Neg(a)*Neg(b)/Sleg(x))*(dLog(y(b,x))*y(b,x)
     .   - y(b,x) + 1.d0)/((y(b,x) - 1.d0)**2.d0)

         elseif((a.eq.b.or.dabs(y(a,x)-y(b,x)).le.0.01d0).and.
     .   dabs(1.d0-y(a,x)).gt.0.01d0.and.
     .   dabs(1.d0-y(b,x)).gt.0.01d0)then cond907
 
         fzn(x,a,b) = (-y(a,x)**2.d0 + y(a,x) + dLog(y(a,x)))/
     .   (y(a,x) - 1.d0)**2.d0 

         gzn(x,a,b) = (Neg(b)**2.d0/Sleg(x))*((-y(a,x) + dLog(y(a,x)) 
     .   + 1.d0))/(y(a,x) - 1.d0)**2.d0 

            else cond907

               fzn(x,a,b) = - 3.d0/2.d0 - (y(a,x) - 1.d0)/3.d0

               gzn(x,a,b) = -Neg(a)*Neg(b)/(2.d0*SLeg(x))

            endif cond907
	

	enddo loop179
	enddo loop178
	enddo loop177

C     =======================================
C     Loop Functions for Chargino Z-Penguins.
C     =======================================

	
	loop180: Do x = 1,3 
        loop181: Do a = 1,2	
	loop182: Do b = 1,2

	fzc(x,a,b) = 0.d0
	gzc(x,a,b) = 0.d0
	
	cond908: If(a.ne.b.and.dabs(z(a,x)-z(b,x)).gt.0.01d0.and.
     .   dabs(1.d0-z(a,x)).gt.0.01d0.and.
     .   dabs(1.d0-z(b,x)).gt.0.01d0)then

	 fzc(x,a,b) = dLog(z(a,x)) + 1.d0/(z(a,x) - z(b,x))*
     . 	( (z(a,x)**2.d0*dLog(z(a,x)))/(1.d0 - z(a,x)) -
     . 	 (z(b,x)**2.d0*dLog(z(b,x)))/(1.d0 - z(b,x)) )

	 gzc(x,a,b) = (Ceg(a)*Ceg(b)/SNeg(x))*(1.d0/(z(a,x)- z(b,x)))*
     . 	( (z(a,x)*dLog(z(a,x)))/(1.d0 - z(a,x)) -
     . 	 (z(b,x)*dLog(z(b,x)))/(1.d0 - z(b,x)) )
	
         elseif(a.ne.b.and.dabs(z(a,x)-z(b,x)).gt.0.01d0.and.
     .   dabs(1.d0-z(a,x)).gt.0.01d0.and.dabs(1.d0-z(b,x)).le.
     .   0.01d0)then cond908

         fzc(x,a,b) = (-2.d0*dLog(z(a,x))*z(a,x) + z(a,x) + dlog(z(a,x))
     .    - 1.d0)/((z(a,x) -1.d0)**2.d0) 

          gzc(x,a,b) = - (Ceg(a)*Ceg(b)/SNeg(x))*(dLog(z(a,x))*z(a,x)
     .   - z(a,x) + 1.d0)/((z(a,x) - 1.d0)**2.d0)

         elseif(a.ne.b.and.dabs(z(a,x)-z(b,x)).gt.0.01d0.and.
     .   dabs(1.d0-z(b,x)).gt.0.01d0.and.dabs(1.d0-z(a,x)).le.
     .   0.01d0)then cond908 

         fzc(x,a,b) = (- dLog(z(b,x))*z(b,x)**2.d0 + z(b,x) - 
     .    1.d0)/(z(b,x) - 1.d0)**2.d0 
              
         gzc(x,a,b) = - (Ceg(a)*Ceg(b)/SNeg(x))*(dLog(z(b,x))*z(b,x)
     .   - z(b,x) + 1.d0)/((z(b,x) - 1.d0)**2.d0)

         elseif((a.eq.b.or.dabs(z(a,x)-z(b,x)).le.0.01d0).and.
     .   dabs(1.d0-z(a,x)).gt.0.01d0.and.
     .   dabs(1.d0-z(b,x)).gt.0.01d0)then cond908
 
         fzc(x,a,b) = (-z(a,x)**2.d0 + z(a,x) + dLog(z(a,x)))/
     .   (z(a,x) - 1.d0)**2.d0 

          gzc(x,a,b) = (Ceg(b)**2.d0/SNeg(x))*((-z(a,x) + dLog(z(a,x)) 
     .   + 1.d0))/(z(a,x) - 1.d0)**2.d0 

            else cond908

            fzc(x,a,b) = - 3.d0/2.d0 - (z(a,x) - 1.d0)/3.d0

            gzc(x,a,b) = -Ceg(a)*Ceg(b)/(2.d0*SNeg(x))

            endif cond908
	
	enddo loop182
	enddo loop181
	enddo loop180

C	-----------------------------------------------
C	Left  Contributions (Charginos) (LL MI)
C	------------------------------------------

        ZpengMueeeC = 0.d0
        ZpengTAUMueeC = 0.d0
        ZpengTauMUMUMuC = 0.d0
        ZpengTAUeeeC = 0.d0

	loop183: Do x = 1, 3
	loop184: Do a = 1, 2 
	loop185: Do b = 1, 2 

	ZpengMueeeC = ZpengMueeeC + CRlE(a,x)*CRlMU(b,x)*
     .  piconst*(OCL(a,2)*OCL(b,2)*gzc(x,a,b) - OCR(a,2)*OCR(b,2)*
     .  fzc(x,a,b)/2.d0)

C	ZpengTAUMueeC = ZpengTAUMueeC +  CRlMU(a,x)*CRlTAU(b,x)*
C    .  piconst*(OCL(a,2)*OCL(b,2)*gzc(x,a,b) - OCR(a,2)*OCR(b,2)*
C    .  fzc(x,a,b)/2.d0)

	ZpengTAUMuMuMuC = ZpengTAUMuMuMuC +  CRlMU(a,x)*CRlTAU(b,x)*
     .  piconst*(OCL(a,2)*OCL(b,2)*gzc(x,a,b) - OCR(a,2)*OCR(b,2)*
     .  fzc(x,a,b)/2.d0)

	ZpengTAUeeeC = ZpengTAUeeeC + CRlE(a,x)*CRlTAU(b,x)*piconst*
     .  (OCL(a,2)*OCL(b,2)*gzc(x,a,b) - OCR(a,2)*OCR(b,2)*
     .  fzc(x,a,b)/2.d0)
	
	enddo loop185
	enddo loop184
	enddo loop183


C	-----------------------------------------------
C	Left  Contributions (Neutralinos) (LL MI)
C	----------------------------------------

        ZpengMueeeNL = 0.d0 
        ZpengTauMueeNL = 0.d0
        ZpengTauMuMuMuNL = 0.d0
        ZpengTauEEENL = 0.d0


	loop186: Do x = 1,6
	loop187: Do a = 1,4
	loop188: Do b = 1,4 

	ZpengMueeeNL = ZpengMueeeNL + NRlE(a,x)*NRlMU(b,x)*piconst*
     .   (ON(a,3)*ON(b,3) - ON(a,4)*ON(b,4))*(fzn(x,a,b) + 
     .    2.d0*gzn(x,a,b))	

C	ZpengTauMueeNL = ZpengTauMueeNL + NRlMu(a,x)*NRlTAU(b,x)*
C    .   piconst*(ON(a,3)*ON(b,3) - ON(a,4)*ON(b,4))*(fzn(x,a,b) 
C    .   + 2.d0*gzn(x,a,b))	

	ZpengTauMuMuMuNL = ZpengTauMuMuMuNL + NRlMu(a,x)*NRlTAU(b,x)*
     .   piconst*(ON(a,3)*ON(b,3) - ON(a,4)*ON(b,4))*(fzn(x,a,b) 
     .   + 2.d0*gzn(x,a,b))	

	ZpengTauEEENL = ZpengTauEEENL + NRlE(a,x)*NRlTAU(b,x)*
     .   piconst*(ON(a,3)*ON(b,3)  - ON(a,4)*ON(b,4))*(fzn(x,a,b) 
     .   + 2.d0*gzn(x,a,b))	

	enddo loop188
	enddo loop187
	enddo loop186

C	-----------------------------------------------
C	Right  Contributions (Neutralinos Only) (RR MI)
C	-----------------------------------------------

        ZpengMueeeNR = 0.d0
        ZpengTauMueeNR = 0.d0
        ZpengTauMuMuMuNR = 0.d0
        ZpengTauEEENR = 0.d0

	loop189: Do x = 1,6
	loop190: Do a = 1,4
	loop191: Do b = 1,4 

	ZpengMueeeNR = ZpengMueeeNR - NLlE(a,x)*NLlMU(b,x)*piconst
     .   *(ON(a,3)*ON(b,3) - ON(a,4)*ON(b,4))*(fzn(x,a,b) 
     .   + 2.d0*gzn(x,a,b))	

C	ZpengTauMueeNR = ZpengTauMueeNR -  NLlMu(a,x)*NLlTAU(b,x)
C     .  *piconst*(ON(a,3)*ON(b,3) - ON(a,4)*ON(b,4))*(fzn(x,a,b) 
C     .	+ 2.d0*gzn(x,a,b))	

        ZpengTauMuMuMuNR = ZpengTauMuMuMuNR -  NLlMu(a,x)*NLlTAU(b,x)
     .  *piconst*(ON(a,3)*ON(b,3) - ON(a,4)*ON(b,4))*(fzn(x,a,b) 
     .	+ 2.d0*gzn(x,a,b))	

	ZpengTauEEENR = ZpengTauEEENR - NLlE(a,x)*NLlTAU(b,x)*piconst
     .  *(ON(a,3)*ON(b,3) - ON(a,4)*ON(b,4))*(fzn(x,a,b) + 
     .  2.d0*gzn(x,a,b))	

	enddo loop191
	enddo loop190
	enddo loop189

C     =====================================================
C     Loop Functions for the BOX Contributions 
C     =====================================================
C     Scaling of the masses so that it is numerically stable!!
C     ========================================================

        loopsc1 : do a = 1,4
        ScNeg(a) = Neg(a)/1000.d0
        enddo loopsc1

        loopsc2: do a = 1, 2
        ScCeg(a) = Ceg(a)/1000.d0
        enddo loopsc2

        loopsc3: do x = 1,6
        ScSLeg(x) = SLeg(x)/1000000.d0
        enddo loopsc3
        
        loopsc4: do x = 1, 3
        ScSNeg(x) = SNeg(x)/1000000.d0
        enddo loopsc4

        loopsc5 : do x = 1, 6
        ScSUeg(x) = SUeg(x)/1000000.d0
        enddo loopsc5 

        loopsc6: do x = 1, 6
        ScSDeg(x) = SDeg(x)/1000000.d0
        enddo loopsc6

C     ========================================================
C     I3-type functions for chargino/neutralino contributions.
C     ========================================================
C       I3n for neutralinos. 
C       --------------------

        diff1 = 0.000001d0

        loop192: do a = 1, 4
         loop193: do b = 1, 4
           loop194: do x = 1, 6

          ScI3n(a,b,x) = 0.d0 

        cond909:If(a.ne.b.and.dabs(ScNeg(a)**2.d0-ScNeg(b)**2.d0).
     .  gt.diff1.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff1
     .  .and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).gt.diff1)then

        ScI3n(a,b,x) = (-2.d0*piconst)*(ScNeg(b)*ScNeg(b)*
     . dLog(ScNeg(b)*ScNeg(b)/(ScNeg(a)*ScNeg(a)))/((ScNeg(b)*
     . ScNeg(b) - ScSleg(x))*(ScNeg(b)*ScNeg(b)-ScNeg(a)*
     . ScNeg(a)))+(ScSLeg(x)*dLog(ScSLeg(x)/(ScNeg(a)*ScNeg(a))))
     . /((ScSleg(x) - ScNeg(b)*ScNeg(b))*(ScSLeg(x) - ScNeg(a)*
     . ScNeg(a))))  

        elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a) - ScNeg(b)*
     .  ScNeg(b)).gt.diff1.and.dabs(ScNeg(a)*ScNeg(a)-
     .  ScSLeg(x)).gt.diff1.and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x))
     .  .le.diff1)then cond909

C        ScI3n(a,b,x) = (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*
C     .  dLog(ScNeg(a)*ScNeg(a)/(ScNeg(b)*ScNeg(b))))/((ScNeg(a)
C     .  *ScNeg(a) - ScNeg(b)*ScNeg(b))**2.d0) + 1.d0/(ScNeg(b)*
c     .  ScNeg(b)- ScNeg(a)*ScNeg(a)) )

          ScI3n(a,b,x) = I3f2d(ScNeg(a)*ScNeg(a),ScSLeg(x))

        elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a) - ScNeg(b)*
     .  ScNeg(b)).gt.diff1.and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x))
     .  .gt.diff1.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).le.
     .  diff1)then cond909

C        ScI3n(a,b,x) = (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*
c     .  dLog(ScNeg(b)*ScNeg(b)/(ScNeg(a)*ScNeg(a))))/((ScNeg(b)
c     .  *ScNeg(b)-ScNeg(a)*ScNeg(a))**2.d0) +1.d0/(ScNeg(a)*ScNeg(a)
C     .  -ScNeg(b)*ScNeg(b)) )

           ScI3n(a,b,x) = I3f2d(ScNeg(b)*ScNeg(b),ScSLeg(x))

       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .le.diff1).and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).gt.diff1
     . .and.dabs(ScNeg(a)*ScNeg(a)-ScSleg(x)).gt.diff1)then cond909

c        ScI3n(a,b,x) = (-2.d0*piconst)*((ScSLeg(x)*dLog(ScSLeg(x)/
C     . (ScNeg(a)*ScNeg(a))))/((ScSLeg(x)-ScSNeg(a)*ScSNeg(a))**2.d0)
C     . + 1.d0/(ScNeg(a)*ScNeg(a)- ScSLeg(x)) )

          ScI3n(a,b,x) = I3f2d(ScSLeg(x),ScNeg(b)*ScNeg(b))
         
       else cond909

         ScI3n(a,b,x) = -1.d0/(32.d0*ScNeg(b)*ScNeg(b)*pi*pi)

         endif cond909

         enddo loop194
         enddo loop193
         enddo loop192

C     =======================rescaling of the loop functions==========

         loopresclI3na: do a = 1, 4
         loopresclI3nb: do b = 1, 4
         loopresclI3nc: do x = 1, 6

         I3n(a,b,x) = ScI3n(a,b,x)/1000000.d0

         enddo loopresclI3nc
         enddo loopresclI3nb
         enddo loopresclI3na

C     ===================================================================
C       ==================
C       I3c for Charginos
C       ==================
        loop195:do a = 1, 2
        loop196:do b = 1, 2
        loop197:do x = 1, 3

        ScI3c(a,b,x) = 0.d0 

        cond910:If(a.ne.b.and.dabs(ScCeg(a)**2.d0-ScCeg(b)**2.d0).
     .  gt.diff1.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff1
     .  .and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).gt.diff1)then

        ScI3c(a,b,x) = (-2.d0*piconst)*(ScCeg(b)*ScCeg(b)*
     . dLog(ScCeg(b)*ScCeg(b)/(ScCeg(a)*ScCeg(a)))/((ScCeg(b)*
     . ScCeg(b) - ScSNeg(x))*(ScCeg(b)*ScCeg(b)-ScCeg(a)*
     . ScCeg(a)))+(ScSNeg(x)*dLog(ScSNeg(x)/(ScCeg(a)*ScCeg(a))))
     . /((ScSNeg(x) - ScCeg(b)*ScCeg(b))*(ScSNeg(x) - ScCeg(a)*
     . ScCeg(a))))  

        elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a) - ScCeg(b)*
     .  ScCeg(b)).gt.diff1.and.dabs(ScCeg(a)*ScCeg(a)-
     .  ScSNeg(x)).gt.diff1.and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x))
     .  .le.diff1)then cond910

          ScI3c(a,b,x) = I3f2d(ScCeg(a)*ScCeg(a),ScSNeg(x))

        elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a) - ScCeg(b)*
     .  ScCeg(b)).gt.diff1.and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x))
     .  .gt.diff1.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).le.
     .  diff1)then cond910

           ScI3c(a,b,x) = I3f2d(ScCeg(b)*ScCeg(b),ScSNeg(x))

       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .le.diff1).and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).gt.diff1
     . .and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff1)then cond910

          ScI3c(a,b,x) = I3f2d(ScSNeg(x),ScCeg(b)*ScCeg(b))
         
       else cond910

         ScI3c(a,b,x) = -1.d0/(32.d0*ScCeg(b)*ScCeg(b)*pi*pi)

         endif cond910

         enddo loop197
         enddo loop196
         enddo loop195

C     =======================rescaling of the loop functions==========

         loopresclI3ca: do a = 1, 2
         loopresclI3cb: do b = 1, 2
         loopresclI3cc: do x = 1, 3

         I3c(a,b,x) = ScI3c(a,b,x)/1000000.d0

         enddo loopresclI3cc
         enddo loopresclI3cb
         enddo loopresclI3ca

C     ===================================================================
C     =================================================
C     I4 type functions for neutralinos and Charginos. 
C     =================================================
C       For Neutralinos 
C       ---------------
         
         diff2 = 0.000001d0

       loop198: do a = 1,4
       loop199: do b = 1,4
       loop200: do x = 1,6
       loop201: do t = 1,6

         ScI4n(a,b,x,t)=0.d0
        
        cond911: If(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*
     . ScNeg(b)).gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(t)).gt.diff2
     . .and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(b)*ScNeg(b)-ScSLeg(t)).gt.diff2.and.x.ne.t.and.
     . dabs(ScSLeg(x) -ScSLeg(t)).gt.diff2)then

        ScI4n(a,b,x,t) = (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*
     . dLog(ScNeg(b)*ScNeg(b)/(ScNeg(a)*ScNeg(a))))/((ScNeg(b)*ScNeg(b)
     .  -ScNeg(a)*ScNeg(a))*(ScNeg(b)*ScNeg(b)-ScSleg(x))*(ScNeg(b)*
     .  ScNeg(b)-ScSLeg(t))) + (ScSLeg(x)*dLog(ScSLeg(x)/(ScNeg(a)*
     .  ScNeg(a))))/((ScSleg(x)-ScNeg(a)*ScNeg(a))*(ScSleg(x)-ScNeg(b)
     .  *ScNeg(b))*(ScSleg(x)-ScSLeg(t)) ) + (ScSLeg(t)*dLog(ScSLeg(t)
     . /(ScNeg(a)*ScNeg(a))))/((ScSleg(t)-ScNeg(a)*ScNeg(a))*(ScSleg(t)
     .  -ScNeg(b)*ScNeg(b))*(ScSleg(t)-ScSLeg(x)))) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSLeg(t)).gt.diff2.and.(x.eq.t.or.dabs(ScSLeg(x)-ScSLeg(t))
     . .le.diff2)) then cond911

       ScI4n(a,b,x,t) = I4f3(ScNeg(a)*ScNeg(a),ScNeg(b)*ScNeg(b),
     . ScSLeg(t)) 

C       (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*
C     . dLog(ScNeg(a)*ScNeg(a)/ScSLeg(t)))/((ScNeg(a)*ScNeg(a)-
C     . ScNeg(b)*ScNeg(b))*(ScNeg(a)*ScNeg(a)-ScSleg(t))**2.d0) + 
C     . (ScNeg(b)*ScNeg(b)*dLog(ScNeg(b)*ScNeg(b)/ScSLeg(t)))/((ScNeg(b)
C     . *ScNeg(b)-ScNeg(a)*ScNeg(a))*(ScNeg(b)*ScNeg(b)-ScSleg(t))
C     . **2.d0) + 1.d0/( (ScSleg(t)-ScNeg(a)*ScNeg(a))*(ScSleg(t)-
C     . ScNeg(b)*ScNeg(b))) ) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     .  dabs(ScNeg(a)*ScNeg(a)-ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)
     . *ScNeg(b)-ScSLeg(x)).gt.diff2.and.x.ne.t.and.dabs(ScSLeg(x)
     . -ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(t))
     . .le.diff2)then cond911

       ScI4n(a,b,x,t) = I4f3(ScNeg(a)*ScNeg(a),ScSLeg(x),ScSLeg(t))

C       ScI4n(a,b,x,t) = (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*dLog(
C     . ScNeg(a)*ScNeg(a)/ScSLeg(t)))/((ScNeg(a)*ScNeg(a)-ScSLeg(x))
C     . *(ScNeg(a)*ScNeg(a)-ScSleg(t))**2.d0)+(ScSLeg(x)*dLog(ScSLeg(x)
C     . /ScSLeg(t)))/((ScSLeg(x)-ScNeg(a)*ScNeg(a))*(ScSLeg(x)-
C     . ScSleg(t))**2.d0) + 1.d0/((ScSleg(t)-ScNeg(a)*ScNeg(a))*
C     . (ScSleg(t)-ScSLeg(x))) ) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSLeg(t)).gt.diff2.and.x.ne.t.and.dabs(ScSLeg(x)-
     . ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).le.
     . diff2)then cond911 

C       ScI4n(a,b,x,t) = (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*dLog(
C     . ScNeg(a)*ScNeg(a)/ScSLeg(x)))/((ScNeg(a)*ScNeg(a)-ScSLeg(t))*
C     . (ScNeg(a)*ScNeg(a)-ScSleg(x))**2.d0)+(ScSLeg(t)*dLog(ScSLeg(t)
C     . /ScSLeg(x)))/((ScSLeg(t)-ScNeg(a)*ScNeg(a))*(ScSLeg(t)-ScSleg(x))
C     . **2.d0) + 1.d0/((ScSleg(x)-ScNeg(a)*ScNeg(a))*(ScSleg(x)-
C     . ScSLeg(t))) ) 

       ScI4n(a,b,x,t) = I4f3(ScNeg(a)*ScNeg(a),ScSLeg(t),ScSLeg(x))

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(b)*ScNeg(b)-ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)
     . *ScNeg(b)-ScSLeg(x)).gt.diff2.and.x.ne.t.and.dabs(ScSLeg(x)-
     . ScSLeg(t)).gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(t)).
     . le.diff2)then cond911 

       ScI4n(a,b,x,t) = I4f3(ScNeg(b)*ScNeg(b),ScSLeg(x),ScSLeg(t))

C      (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*dLog(
C    . ScNeg(b)*ScNeg(b)/ScSLeg(t)))/((ScNeg(b)*ScNeg(b)-ScSLeg(x))
C    . *(ScNeg(b)*ScNeg(b)-ScSleg(t))**2.d0)+(ScSLeg(x)*dLog(ScSLeg(x)
C    . /ScSLeg(t)))/((ScSLeg(x)-ScNeg(b)*ScNeg(b))*(ScSLeg(x)-
C    . ScSLeg(t))**2.d0)+1.d0/( (ScSleg(t)-ScNeg(b)*ScNeg(b))*
C    . (ScSleg(t)-ScSLeg(x))) ) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).gt.
     . diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(t)).gt.diff2.and.dabs(
     . ScNeg(b)*ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)
     . -ScSLeg(t)).gt.diff2.and.x.ne.t.and.dabs(ScSLeg(x)-ScSLeg(t))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).le.diff2)
     . then cond911 

       ScI4n(a,b,x,t) = I4f3(ScNeg(b)*ScNeg(b),ScSLeg(t),ScSLeg(x)) 

C       (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*dLog(
C     . ScNeg(b)*ScNeg(b)/ScSLeg(x)))/((ScNeg(b)*ScNeg(b)-ScSLeg(t))
C     . *(ScNeg(b)*ScNeg(b)-ScSleg(x))**2.d0)+(ScSLeg(t)*dLog(ScSLeg(t)
C     . /ScSLeg(x)))/((ScSLeg(t)-ScNeg(b)*ScNeg(b))*(ScSLeg(t)-ScSleg(x)
C     . )**2.d0) + 1.d0/((ScSleg(x)-ScNeg(b)*ScNeg(b))*(ScSleg(x)- 
C     .  ScSLeg(t))) ) 

       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .le.diff2).and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSLeg(x)).gt.diff2.and.x.ne.t.and.dabs(ScSLeg(x)-ScSLeg(t))
     . .gt.diff2)then cond911

       ScI4n(a,b,x,t) = I4f3(ScSLeg(x),ScSLeg(t),ScNeg(b)*ScNeg(b))

C       (-2.d0*piconst)*((ScSLeg(x)*dLog(ScSLeg(x)/
C     . (ScNeg(b)*ScNeg(b))))/(((ScSleg(x)-ScNeg(b)*ScNeg(b))**2.d0)*
C     . (ScSleg(x)-ScSLeg(t))) + (ScSLeg(t)*dLog(ScSLeg(t)/(ScNeg(b)*
C     . ScNeg(b))))/(((ScSleg(t)-ScNeg(b)*ScNeg(b))**2.d0)*(ScSleg(t)- 
C     . ScSLeg(x))) + 1.d0/( (ScNeg(b)*ScNeg(b)-ScSLeg(x))*(ScNeg(b)*
C     . ScNeg(b)-ScSLeg(t))) )

       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).le.
     . diff2).and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).le.diff2.and.
     . x.ne.t.and.dabs(ScSLeg(x)-ScSLeg(t)).gt.diff2)then cond911

C       ScI4n(a,b,x,t)=(-2.d0*piconst)*((ScSLeg(t)*dLog(ScSLeg(t)/
C     . ScSleg(x)))/((ScSLeg(t)-ScSLeg(x))**3.d0)-(ScSLeg(t) + 
C     . ScSLeg(x))/(2.d0*ScSLeg(x)*(ScSLeg(x)-ScSLeg(t))**2.d0 ))

         ScI4n(a,b,x,t) = I4f3d1(ScSLeg(t),ScSLeg(x))

       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).le.
     . diff2).and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(t)).le.diff2.and.
     . x.ne.t.and.dabs(ScSLeg(x)-ScSLeg(t)).gt.diff2)then cond911

       ScI4n(a,b,x,t) = I4f3d1(ScSLeg(x),ScSLeg(t))

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).gt.
     . diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).le.diff2.and.
     . (x.eq.t.or.dabs(ScSLeg(x)-ScSLeg(t)).le.diff2))then cond911

       ScI4n(a,b,x,t) =  I4f3d1(ScNeg(b)*ScNeg(b),ScSLeg(t))

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).gt.
     . diff2.and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).le.diff2.and.
     . (x.eq.t.or.dabs(ScSLeg(x)-ScSLeg(t)).le.diff2))then cond911

       ScI4n(a,b,x,t) = I4f3d1(ScNeg(a)*ScNeg(a),ScSLeg(t))
 
      elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .le.diff2).and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSLeg(x)).gt.diff2.and.(x.eq.t.or.dabs(ScSLeg(x)-ScSLeg(t))
     . .le.diff2))then cond911
       
         ScI4n(a,b,x,t) = I4f2d2d(ScNeg(b)*ScNeg(b),ScSLeg(t))
         
      elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).le.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSLeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSLeg(x)).le.diff2.and.x.ne.t.and.dabs(ScSLeg(x)-ScSLeg(t))
     . .gt.diff2)then cond911
  
         ScI4n(a,b,x,t) = I4f2d2d(ScSLeg(t),ScSLeg(x))
  
      elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(t)).le.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSLeg(t)).le.diff2.and.x.ne.t.and.dabs(ScSLeg(x)-ScSLeg(t))
     . .gt.diff2)then cond911

       ScI4n(a,b,x,t) = I4f2d2d(ScSLeg(x),ScSLeg(t))

       else cond911
        
       ScI4n(a,b,x,t) = 1.d0/(96.d0*ScSLeg(x)*ScSLeg(x)*pi*pi)

        endif cond911 

       enddo loop201
       enddo loop200
       enddo loop199
       enddo loop198


c     ===============================================================

       loopresclI4na: do a = 1, 4
         loopresclI4nb: do b = 1, 4
         loopresclI4nc: do x = 1, 6
         loopresclI4nd: do t = 1,6

         I4n(a,b,x,t) = ScI4n(a,b,x,t)*(10.d0**(-12.d0))
         enddo loopresclI4nd
         enddo loopresclI4nc
         enddo loopresclI4nb
         enddo loopresclI4na

C     ==================================================================

C       --------------
C       For Charginos
C       --------------

       loop202:do a = 1,2
       loop203:do b = 1,2
       loop204:do x = 1,3
       loop205:do t = 1,3

         ScI4c(a,b,x,t)=0.d0
        
        cond912:If(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     .  .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     .  dabs(ScCeg(a)*ScCeg(a)-ScSNeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     .  ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     .  ScSNeg(t)).gt.diff2.and.dabs(ScSNeg(x)-ScSNeg(t)).gt.diff2)then

       ScI4c(a,b,x,t) = (-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*dLog(
     . ScCeg(b)*ScCeg(b)/(ScCeg(a)*ScCeg(a))))/((ScCeg(b)*ScCeg(b) -
     . ScCeg(a)*ScCeg(a))*(ScCeg(b)*ScCeg(b)-ScSNeg(x))*(ScCeg(b)*
     . ScCeg(b)-ScSNeg(t))) + (ScSNeg(x)*dLog(ScSNeg(x)/(ScCeg(a)*
     . ScCeg(a))))/((ScSNeg(x)-ScCeg(a)*ScCeg(a))*(ScSNeg(x)-ScCeg(b)
     . *ScCeg(b))*(ScSNeg(x)-ScSNeg(t)) ) + (ScSNeg(t)*dLog(ScSNeg(t)/
     . (ScCeg(a)*ScCeg(a))))/((ScSNeg(t)-ScCeg(a)*ScCeg(a))*(ScSNeg(t)- 
     . ScCeg(b)*ScCeg(b))*(ScSNeg(t)-ScSNeg(x)) ) ) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).gt.
     . diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.dabs(
     . ScCeg(a)*ScCeg(a)-ScSNeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)
     . -SNeg(x)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(t)).gt.
     . diff2.and.dabs(ScSNeg(x)-ScSNeg(t)).le.diff2)then cond912

      ScI4c(a,b,x,t) = I4f3(ScCeg(a)*ScCeg(a),ScCeg(b)*ScCeg(b),
     . ScSNeg(t))

C      ScI4c(a,b,x,t) = (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*dLog(
C     . ScCeg(a)*ScCeg(a)/ScSNeg(t)))/((ScCeg(a)*ScCeg(a)-ScCeg(b)*
C     . ScCeg(b))*(ScCeg(a)*ScCeg(a)-ScSNeg(t))**2.d0) + (ScCeg(b)*
C     . ScCeg(b)*dLog(ScCeg(b)*ScCeg(b)/ScSNeg(t)))/((ScCeg(b)*ScCeg(b)
c     . -ScCeg(a)*ScCeg(a))*(ScCeg(b)*ScCeg(b)-ScSNeg(t))**2.d0) + 
c     . 1.d0/( (ScSNeg(t)-ScCeg(a)*ScCeg(a))*(ScSNeg(t)-ScCeg(b)*
C    .  ScCeg(b))) ) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSNeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScSNeg(x)-ScSNeg(t)).gt.
     . diff2.and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(t)).le.diff2)then cond912

         ScI4c(a,b,x,t) = I4f3(ScCeg(a)*ScCeg(a),ScSNeg(x),ScSNeg(t))

c       ScI4c(a,b,x,t) = (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*dLog(
c     . ScCeg(a)*ScCeg(a)/ScSNeg(t)))/((ScCeg(a)*ScCeg(a)-ScSNeg(x))
C     . *(ScCeg(a)*ScCeg(a)-ScSNeg(t))**2.d0)+(ScSNeg(x)*dLog(ScSNeg(x)
C     . /ScSNeg(t)))/((ScSNeg(x)-ScCeg(a)*ScCeg(a))*(ScSNeg(x)-
C     . ScSNeg(t))**2.d0)+1.d0/((ScSNeg(t)-ScCeg(a)*ScCeg(a))*(ScSNeg(t)
C     . - SNeg(x))) ) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).gt.
     . diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.dabs(
     . ScCeg(a)*ScCeg(a)-ScSNeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)
     . -ScSNeg(t)).gt.diff2.and.dabs(ScSNeg(x)-ScSNeg(t)).gt.diff2.and.
     . dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).le.diff2)then cond912

         ScI4c(a,b,x,t) = I4f3(ScCeg(a)*ScCeg(a),ScSNeg(t),ScSNeg(x))

C       ScI4c(a,b,x,t) = (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*dLog(
C     . ScCeg(a)*ScCeg(a)/ScSNeg(x)))/((ScCeg(a)*ScCeg(a)-ScSNeg(t))*
C     . (ScCeg(a)*ScCeg(a)-ScSNeg(x))**2.d0)+(ScSNeg(t)*dLog(ScSNeg(t)
C     . /ScSNeg(x)))/((ScSNeg(t)-ScCeg(a)*ScCeg(a))*(ScSNeg(t)-
C     . ScSNeg(x))**2.d0) + 1.d0/( (ScSNeg(x)- ScCeg(a)*ScCeg(a))*
C     . (ScSNeg(x)- ScSNeg(t))) ) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).gt.
     . diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.dabs(
     . ScCeg(b)*ScCeg(b)-ScSNeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)
     . -ScSNeg(x)).gt.diff2.and.dabs(ScSNeg(x)-ScSNeg(t)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSNeg(t)).le.diff2)then cond912

         ScI4c(a,b,x,t) = I4f3(ScCeg(b)*ScCeg(b),ScSNeg(x),ScSNeg(t))

C       ScI4c(a,b,x,t) = (-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*dLog(
c     . ScCeg(b)*ScCeg(b)/ScSNeg(t)))/((ScCeg(b)*ScCeg(b)-ScSNeg(x))
C     . *(ScCeg(b)*ScCeg(b)-ScSNeg(t))**2.d0) + (ScSNeg(x)*dLog(
C     . ScSNeg(x))/ScSNeg(t))/((ScSNeg(x)-ScCeg(b)*ScCeg(b))*(ScSNeg(x)
C     . -ScSNeg(t))**2.d0) + 1.d0/((ScSNeg(t)-ScCeg(b)*ScCeg(b))*
C     . (ScSNeg(t)-ScSNeg(x))) )

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).gt.
     . diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(t)).gt.diff2.and.dabs(
     . ScCeg(b)*ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)
     . -ScSNeg(t)).gt.diff2.and.dabs(ScSNeg(x)-ScSNeg(t)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).le.diff2)then cond912

         ScI4n(a,b,x,t) = I4f3(ScCeg(b)*ScCeg(b),ScSNeg(t),ScSNeg(x))

C       ScI4c(a,b,x,t) = (-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*dLog(
C     . ScCeg(b)*ScCeg(b)/ScSNeg(x)))/((ScCeg(b)*ScCeg(b)-ScSNeg(t))
C     . *(ScCeg(b)*ScCeg(b)-ScSNeg(x))**2.d0) + (ScSNeg(t)*dLog(
c     . ScSNeg(t)/ScSNeg(x)))/((ScSNeg(t)-ScCeg(b)*ScCeg(b))*(ScSNeg(t)
c     . -ScSNeg(x))**2.d0) + 1.d0/( (ScSNeg(x)-ScCeg(b)*ScCeg(b))*
C     . (ScSNeg(x)- ScSNeg(t))) ) 

       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).le.
     . diff2).and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSNeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSNeg(x)).gt.diff2.and.dabs(ScSNeg(x)-ScSNeg(t)).gt.diff2)
     . then cond912

       ScI4c(a,b,x,t) = I4f3(ScSNeg(x),ScSNeg(t),ScCeg(b)*ScCeg(b))

C       ScI4c(a,b,x,t) = (-2.d0*piconst)*((ScSNeg(x)*dLog(ScSNeg(x)/(
C     . ScCeg(b)*ScCeg(b))))/(((ScSNeg(x)-ScCeg(b)*ScCeg(b))**2.d0)*(
C     . ScSNeg(x)-ScSNeg(t))) + (ScSNeg(t)*dLog(ScSNeg(t)/(ScCeg(b)*
c     . ScCeg(b))))/(((ScSNeg(t)-ScCeg(b)*ScCeg(b))**2.d0)*(ScSNeg(t)- 
C     . ScSNeg(x))) + 1.d0/( (ScCeg(b)*ScCeg(b)-ScSNeg(x))*(ScCeg(b)
c     . *ScCeg(b)-ScSNeg(t))) )

       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).le.
     . diff2).and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSNeg(t)).gt.diff2)then cond912

          ScI4c(a,b,x,t) = I4f3d1(ScSNeg(t),ScSNeg(x))

C      ScI4c(a,b,x,t) = (-2.d0*piconst)*((ScSNeg(t)*dLog(ScSNeg(t)/
c    . ScSNeg(x)))/((ScSNeg(t)-ScSNeg(x))**3.d0)-(ScSNeg(t) + 
C    . ScSNeg(x))/(2.d0*ScSNeg(x)*(ScSNeg(x) - ScSNeg(t))**2.d0 ))

       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).le.
     . diff2).and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(t)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSNeg(t)).gt.diff2)then cond912

          ScI4c(a,b,x,t) = I4f3d1(ScSNeg(x),ScSNeg(t))

C       ScI4c(a,b,x,t)=(-2.d0*piconst)*((ScSNeg(x)*dLog(ScSNeg(x)/
c     . ScSNeg(t)))/((ScSNeg(x)-ScSNeg(t))**3.d0)-(ScSNeg(x)+ScSNeg(t))
c     . /(2.d0*ScSNeg(t)*(ScSNeg(t)-ScSNeg(x))**2.d0))

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSNeg(t)).le.diff2)then cond912

          ScI4c(a,b,x,t) =  I4f3d1(ScCeg(b)*ScCeg(b),ScSNeg(t))

C       ScI4c(a,b,x,t)=(-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*dLog(ScCeg(b)
C     . *ScCeg(b)/(ScCeg(a)*ScCeg(a))))/((ScCeg(b)*ScCeg(b)-ScCeg(a)*
C     . ScCeg(a))**3.d0)-(ScCeg(a)*ScCeg(a)+ScCeg(b)*ScCeg(b))/(2.d0*
C     . ScCeg(a)*ScCeg(a)*(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))**2.d0))

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).gt.
     . diff2.and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSNeg(t)).le.diff2)then cond912

          ScI4c(a,b,x,t) = I4f3d1(ScCeg(a)*ScCeg(a),ScSNeg(t))

C       ScI4c(a,b,x,t) = (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*dLog(
C     . ScCeg(a)*ScCeg(a)/(ScCeg(b)*ScCeg(b))))/((ScCeg(a)*ScCeg(a)- 
C     . ScCeg(b)*ScCeg(b))**3.d0)-(ScCeg(b)*ScCeg(b) +ScCeg(a)*ScCeg(a))
c     . /(2.d0*ScCeg(a)*ScCeg(a)*(ScCeg(b)*ScCeg(b) - ScCeg(a)*
C     . ScCeg(a))**2.d0))

      elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .le.diff2).and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSNeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSNeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSNeg(x)).gt.diff2.and.(x.eq.t.or.dabs(ScSNeg(x)-ScSNeg(t))
     . .le.diff2))then cond912
       
         ScI4c(a,b,x,t) = I4f2d2d(ScCeg(b)*ScCeg(b),ScSNeg(t))
         
      elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).le.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSNeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSNeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSNeg(x)).le.diff2.and.x.ne.t.and.dabs(ScSNeg(x)-ScSNeg(t))
     . .gt.diff2)then cond912
  
         ScI4c(a,b,x,t) = I4f2d2d(ScSNeg(t),ScSNeg(x))
  
      elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(t)).le.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSNeg(t)).le.diff2.and.x.ne.t.and.dabs(ScSNeg(x)-ScSNeg(t))
     . .gt.diff2)then cond912

       ScI4c(a,b,x,t) = I4f2d2d(ScSNeg(x),ScSNeg(t))
       
       else cond912
        
       ScI4c(a,b,x,t) = 1.d0/(96.d0*ScSNeg(x)*ScSNeg(x)*pi*pi)

        endif cond912

       enddo loop205
       enddo loop204
       enddo loop203
       enddo loop202
       
c     ===============================================================

         loopresclI4ca: do a = 1, 2
         loopresclI4cb: do b = 1, 2
         loopresclI4cc: do x = 1, 3
         loopresclI4cd: do t = 1,3

         I4c(a,b,x,t) = ScI4c(a,b,x,t)*(10.d0**(-12.d0))

         enddo loopresclI4cd
         enddo loopresclI4cc
         enddo loopresclI4cb
         enddo loopresclI4ca

C     ==================================================================
C     ==============================
C       J4n and J4c loop functions. 
C     ==============================


      loopi106: do a = 1,4 
      loopi107: do b = 1,4
      loopi108: do x = 1, 6
      loopi109: do t = 1, 6

       J4n(a,b,x,t) = 0.d0

        J4n(a,b,x,t) = SLeg(t)*I4n(a,b,x,t) + I3n(a,b,x)

          enddo loopi109
          enddo loopi108
          enddo loopi107
          enddo loopi106


         loopi110: do a = 1,2 
         loopi111: do b = 1,2
         loopi112: do x = 1,3
         loopi113: do t = 1,3

       J4c(a,b,x,t) = 0.d0

       J4c(a,b,x,t) = SNeg(t)*I4c(a,b,x,t) + I3c(a,b,x)

          enddo loopi113
          enddo loopi112
          enddo loopi111
          enddo loopi110


C     =========================================================
C     BOX CONTRIBUTIONS
C     =========================================================
C     Neutralino Contributions. 
C     ========================

        B1nLMU3E = 0.d0
        B1nLTAu3Mu = 0.d0
        B1nLTau3E = 0.d0

        B2nLMu3E = 0.d0
        B2nLTau3Mu = 0.d0
        B2nLTau3E = 0.d0

        B3nLMu3E = 0.d0
        B3nLTau3Mu = 0.d0
        B3nLTau3E = 0.d0

        B4nLMu3E = 0.d0
        B4nLTau3Mu = 0.d0
        B4nLTau3E = 0.d0

        loop214:  do a = 1, 4
        loop215:  do b = 1, 4
        loop216:  do x = 1, 6
        loop217:  do t = 1, 6

C       -------------------------------------------------------------

        B1nLMu3E = B1nLMu3E + (1.d0/(4.d0*pi*alph))*(0.5d0*NRlMU(a,x)*
     .  NRlE(a,t)*NRlE(b,t)*NRlE(b,x)*J4n(a,b,x,t) + Neg(a)*Neg(b)*
     .  NRlMu(a,x)*NRlE(a,t)*NRlE(b,t)*NRlE(b,x)*I4n(a,b,x,t))

        B1nLTAU3MU = B1nLTAU3MU + (1.d0/(4.d0*pi*alph))*(0.5d0*
     .  NRlTAU(a,x)*NRlMU(a,t)*NRlMU(b,t)*NRlMU(b,x)*J4n(a,b,x,t) 
     .  + Neg(a)*Neg(b)*NRlTAu(a,x)*NRlMU(a,t)*NRlMU(b,t)*NRlMU(b,x)
     .  *I4n(a,b,x,t))

        B1nLTAU3E = B1nLTAU3E + (1.d0/(4.d0*pi*alph))*(0.5d0*
     .  NRlTAU(a,x)*NRlE(a,t)*NRlE(b,t)*NRlE(b,x)*J4n(a,b,x,t) + 
     .  Neg(a)*Neg(b)*NRlTAU(a,x)*NRlE(a,t)*NRlE(b,t)*NRlE(b,x)*
     .  I4n(a,b,x,t))

C       -------------------------------------------------------

        B2nLMu3E = B2nLMu3E+(1.d0/(4.d0*pi*alph))*(0.25d0*(NRlMU(a,x)*
     .  NRlE(a,t)*NLlE(b,t)*NLlE(b,x) + NRlMU(a,x)*NLlE(a,t)*NRlE(b,t)
     .  *NLlE(b,x) - NRlMu(a,x)*NLlE(a,t)*NLlE(b,t)*NRlE(b,x))*
     .  J4n(a,b,x,t) - 0.5d0*Neg(a)*Neg(b)*NRlMu(a,x)*NLlE(a,t)*
     .   NLlE(b,t)*NRlE(b,x)*I4n(a,b,x,t))

        B2nLTAu3MU = B2nLTau3Mu+(1.d0/(4.d0*pi*alph))*(0.25d0*(
     .  NRlTAU(a,x)*NRlMU(a,t)*NLlMU(b,t)*NLlMU(b,x) + NRlTAU(a,x)
     .  *NLlMU(a,t)*NRlMU(b,t)*NLlMU(b,x) - NRlTAu(a,x)*NLlMU(a,t)*
     .  NLlMU(b,t)*NRlMU(b,x))*J4n(a,b,x,t) - 0.5d0*Neg(a)*Neg(b)*
     .  NRlTAU(a,x)*NLlMU(a,t)*NLlMU(b,t)*NRlMU(b,x)*I4n(a,b,x,t))

         B2nLTAu3E = B2nLTau3E+(1.d0/(4.d0*pi*alph))*(0.25d0*(
     .   NRlTAU(a,x)*NRlE(a,t)*NLlE(b,t)*NLlE(b,x) + NRlTAU(a,x)*
     .   NLlE(a,t)*NRlE(b,t)*NLlE(b,x) - NRlTAu(a,x)*NLlE(a,t)*
     .   NLlE(b,t)*NRlE(b,x))*J4n(a,b,x,t) - 0.5d0*Neg(a)*Neg(b)*
     .   NRlTAU(a,x)*NLlE(a,t)*NLlE(b,t)*NRlE(b,x)*I4n(a,b,x,t))

C       -------------------------------------------------------------

        B3nLMu3E = B3nLMu3E + (1.d0/(4.d0*pi*alph))*(Neg(a)*Neg(b))*
     .  (NRlMu(a,x)*NLlE(a,t)*NRlE(b,t)*NLlE(b,x) + 0.5d0*NRlMU(a,x)
     .  *NRlE(a,t)*NLlE(b,t)*NLlE(b,x))*I4n(a,b,x,t)

        B3nLTau3MU = B3nLTAu3MU+(1.d0/(4.d0*pi*alph))*(Neg(a)*Neg(b))*
     .  (NRlTAu(a,x)*NLlMU(a,t)*NRlMU(b,t)*NLlMU(b,x)+0.5d0*
     .  NRlTAU(a,x)*NRlMU(a,t)*NLlMU(b,t)*NLlMU(b,x))*I4n(a,b,x,t)
        
        B3nLTau3E = B3nLTAu3E+(1.d0/(4.d0*pi*alph))*(Neg(a)*Neg(b))*
     .  (NRlTAu(a,x)*NLlE(a,t)*NRlE(b,t)*NLlE(b,x)+0.5d0*NRlTAU(a,x)
     .  *NRlE(a,t)*NLlE(b,t)*NLlE(b,x))*I4n(a,b,x,t)
        
C       -------------------------------------------------------------

        B4nLMu3E = B4nLMu3E + (1.d0/(32.d0*pi*alph))*(Neg(a)*Neg(b))
     .  *(NRlMu(a,x)*NRlE(a,t)*NLlE(b,t)*NLlE(b,x))*I4n(a,b,x,t)

        B4nLTAU3MU = B4nLTAU3MU + (1.d0/(32.d0*pi*alph))*(Neg(a)*Neg(b))
     .  *(NRlTAU(a,x)*NRlMU(a,t)*NLlMU(b,t)*NLlMU(b,x))*I4n(a,b,x,t)

        B4nLTAU3E = B4nLTAU3E + (1.d0/(32.d0*pi*alph))*(Neg(a)*Neg(b))
     .  *(NRlTAU(a,x)*NRlE(a,t)*NLlE(b,t)*NLlE(b,x))*I4n(a,b,x,t)

C       -------------------------------------------------------------
      
        enddo loop217
        enddo loop216
        enddo loop215
        enddo loop214

C     ===================
C     Right Contributions
C     ===================

        B1nRMu3E = 0.d0
        B1nRTau3E = 0.d0
        B1nRTau3Mu = 0.d0

        B2nRMu3E = 0.d0
        B2nRTau3E = 0.d0
        B2nRTau3Mu = 0.d0

        B3nRMu3E = 0.d0
        B3nRTau3E = 0.d0
        B3nRTau3Mu = 0.d0

        B4nRMu3E = 0.d0
        B4nRTau3E = 0.d0
        B4nRTau3Mu = 0.d0

         loop218: do a = 1, 4
         loop219: do b = 1, 4
         loop220: do x = 1, 6
         loop221: do t = 1, 6

C       --------------------------------------------------------

        B1nRMu3E = B1nRMu3E + (1.d0/(4.d0*pi*alph))*(0.5d0*NLlMU(a,x)*
     .  NLlE(a,t)*NLlE(b,t)*NLlE(b,x)*J4n(a,b,x,t) + Neg(a)*Neg(b)*
     .  NLlMu(a,x)*NLlE(a,t)*NLlE(b,t)*NLlE(b,x)*I4n(a,b,x,t))

        B1nRTAu3MU = B1nRTAu3MU + (1.d0/(4.d0*pi*alph))*
     .  (0.5d0*NLlTaU(a,x)*NLlMU(a,t)*NLlMU(b,t)*NLlMU(b,x)*J4n(a,b,x,t) 
     .  + Neg(a)*Neg(b)*NLlTAu(a,x)*NLlMU(a,t)*NLlMU(b,t)*NLlMU(b,x)
     .  *I4n(a,b,x,t))

        B1nRTAu3E = B1nRTAu3E + (1.d0/(4.d0*pi*alph))*(0.5d0*
     .  NLlTAU(a,x)*NLlE(a,t)*NLlE(b,t)*NLlE(b,x)*J4n(a,b,x,t) + 
     .  Neg(a)*Neg(b)*NLlTAu(a,x)*NLlE(a,t)*NLlE(b,t)*NLlE(b,x)*
     .  I4n(a,b,x,t))

C       --------------------------------------------------------

        B2nRMu3E = B2nRMu3E + (1.d0/(4.d0*pi*alph))*(0.25d0*(
     .  NLlMU(a,x)*NLlE(a,t)*NRlE(b,t)*NRlE(b,x) + NLlMU(a,x)*
     .  NRlE(a,t)*NLlE(b,t)*NRlE(b,x) - NLlMu(a,x)*NRlE(a,t)*
     .  NRlE(b,t)*NLlE(b,x))*J4n(a,b,x,t) - 0.5d0*Neg(a)*Neg(b)*
     .  NLlMu(a,x)*NRlE(a,t)*NRlE(b,t)*NLlE(b,x)*I4n(a,b,x,t))

        B2nRTAu3MU = B2nRTAu3MU + (1.d0/(4.d0*pi*alph))*(0.25d0*(
     .  NLlTAU(a,x)*NLlMU(a,t)*NRlMU(b,t)*NRlMU(b,x) + NLlTAU(a,x)*
     .  NRlMU(a,t)*NLlMU(b,t)*NRlMU(b,x) - NLlTAu(a,x)*NRlMU(a,t)*
     .  NRlMU(b,t)*NLlMU(b,x))*J4n(a,b,x,t) - 0.5d0*Neg(a)*Neg(b)*
     .  NLlTAu(a,x)*NRlMU(a,t)*NRlMU(b,t)*NLlMU(b,x)*I4n(a,b,x,t))

        B2nRTAu3E = B2nRTAu3E + (1.d0/(4.d0*pi*alph))*(0.25d0*(
     .  NLlTAU(a,x)*NLlE(a,t)*NRlE(b,t)*NRlE(b,x) + NLlTAU(a,x)*
     .  NRlE(a,t)*NLlE(b,t)*NRlE(b,x) - NLlTAu(a,x)*NRlE(a,t)*
     .  NRlE(b,t)*NLlE(b,x))*J4n(a,b,x,t) - 0.5d0*Neg(a)*Neg(b)*
     .  NLlTAu(a,x)*NRlE(a,t)*NRlE(b,t)*NLlE(b,x)*I4n(a,b,x,t))
  
C       --------------------------------------------------------

        B3nRMu3E = B3nRMu3E + (1.d0/(4.d0*pi*alph))*(Neg(a)*Neg(b))
     .  *(NLlMu(a,x)*NRlE(a,t)*NLlE(b,t)*NRlE(b,x) + 0.5d0*NLlMU(a,x)
     .  *NLlE(a,t)*NRlE(b,t)*NRlE(b,x))*I4n(a,b,x,t)

        B3nRTAu3MU = B3nRTAu3MU + (1.d0/(4.d0*pi*alph))*(Neg(a)*Neg(b))
     .  *(NLlTau(a,x)*NRlMU(a,t)*NLlMu(b,t)*NRlMU(b,x) + 
     .   0.5d0*NLlTAU(a,x)*NLlMu(a,t)*NRlMu(b,t)*NRlMu(b,x))
     .   *I4n(a,b,x,t)

         B3nRTAu3E = B3nRTAu3E + (1.d0/(4.d0*pi*alph))*(Neg(a)*Neg(b))
     .  *(NLlTAu(a,x)*NRlE(a,t)*NLlE(b,t)*NRlE(b,x) + 0.5d0*NLlTAU(a,x)
     .  *NLlE(a,t)*NRlE(b,t)*NRlE(b,x))*I4n(a,b,x,t)
        
C       --------------------------------------------------------

        B4nRMu3E = B4nRMu3E + (1.d0/(32.d0*pi*alph))*(Neg(a)*Neg(b))
     .  *(NLlMu(a,x)*NLlE(a,t)*NRlE(b,t)*NRlE(b,x))*I4n(a,b,x,t)

        B4nRTAu3Mu = B4nRTau3Mu +(1.d0/(32.d0*pi*alph))*(Neg(a)*Neg(b))
     .  *(NLlTAu(a,x)*NLlMU(a,t)*NRlMU(b,t)*NRlMU(b,x))*I4n(a,b,x,t)

        B4nRTAu3E = B4nRTAu3E + (1.d0/(32.d0*pi*alph))*(Neg(a)*Neg(b))
     .  *(NLlTAu(a,x)*NLlE(a,t)*NRlE(b,t)*NRlE(b,x))*I4n(a,b,x,t)

C       --------------------------------------------------------
        enddo loop221
        enddo loop220
        enddo loop219
        enddo loop218

C     ======================
C     Chargino Contributions
C     =======================
C     Left Amplitudes
C     ===============


        B1cLMu3E = 0.d0
        B1cLTau3E = 0.d0
        B1cLTau3Mu = 0.d0

        B2CLMu3E = 0.d0
        B2CLTau3E = 0.d0
        B2CLTau3Mu = 0.d0

        B3CLMu3E = 0.d0
        B3CLTau3E = 0.d0
        B3cLTau3Mu = 0.d0

        loop222: do a = 1, 2
        loop223: do b = 1, 2
        loop224: do x = 1, 3
        loop225: do t = 1, 3

C       --------------------------------------------------
        B1cLMu3E = B1CLMu3E + (1.d0/(8.d0*pi*alph))*(CRlMU(a,x)
     .  *CRlE(a,t)*CRlE(b,t)*CRlE(b,x)*J4c(a,b,x,t))

        B1cLTAu3MU= B1CLTAu3MU + (1.d0/(8.d0*pi*alph))*CRlTAU(a,x)
     .  *CRlMU(a,t)*CRlMU(b,t)*CRlMU(b,x)*J4c(a,b,x,t)

         B1cLTAu3E = B1CLTau3E + (1.d0/(8.d0*pi*alph))*CRlTAU(a,x)
     .   *CRlE(a,t)*CRlE(b,t)*CRlE(b,x)*J4c(a,b,x,t)

C       --------------------------------------------------
        
        B2cLMu3E = B2CLMu3E + (1.d0/(4.d0*pi*alph))*(0.25d0*
     .  CRlMU(a,x)*CRlE(a,t)*CLlE(b,t)*CLlE(b,x)*J4c(a,b,x,t)
     .  - 0.5d0*Ceg(a)*Ceg(b)*CRlMu(a,x)*CLlE(a,t)*CLlE(b,t)*
     .  CRlE(b,x)*I4c(a,b,x,t))

        B2cLTAu3MU = B2CLTAu3MU + (1.d0/(4.d0*pi*alph))*(0.25d0*
     .  CRlTAU(a,x)*CRlMU(a,t)*CLlMU(b,t)*CLlMu(b,x)*J4c(a,b,x,t)
     .  - 0.5d0*Ceg(a)*Ceg(b)*CRlTAu(a,x)*CLlMu(a,t)*CLlMu(b,t)*
     .  CRlMu(b,x)*I4c(a,b,x,t))

        B2cLTAu3E = B2CLTau3E + (1.d0/(4.d0*pi*alph))*(0.25d0*
     .  CRlTAU(a,x)*CRlE(a,t)*CLlE(b,t)*CLlE(b,x)*J4c(a,b,x,t)
     .  - 0.5d0*Ceg(a)*Ceg(b)*CRlTAu(a,x)*CLlE(a,t)*CLlE(b,t)*
     .  CRlE(b,x)*I4c(a,b,x,t))

C       --------------------------------------------------

        B3CLMu3E = B3CLMu3E + (1.d0/(4.d0*pi*alph))*Ceg(a)*Ceg(b)
     .   *CRlMu(a,x)*CLlE(a,t)*CRlE(b,t)*CLlE(b,x)*I4c(a,b,x,t)

       B3CLTAu3MU = B3CLTAu3MU + (1.d0/(4.d0*pi*alph))*Ceg(a)*Ceg(b)
     .   *CRlTAu(a,x)*CLlMU(a,t)*CRlMU(b,t)*CLlMU(b,x)*I4c(a,b,x,t)

        B3CLTAu3E = B3CLTAu3E + (1.d0/(4.d0*pi*alph))*Ceg(a)*Ceg(b)
     .   *CRlTAu(a,x)*CLlE(a,t)*CRlE(b,t)*CLlE(b,x)*I4c(a,b,x,t)

C       --------------------------------------------------

        enddo loop225
        enddo loop224
        enddo loop223
        enddo loop222

        B4CLMu3E = 0.d0
        B4CLTau3E = 0.d0
        B4CLTau3Mu = 0.d0
        
C     =================
C     Right Amplitudes
C     ================ 

        B1cRMu3E = 0.d0
        B1cRTau3E = 0.d0
        B1cRTau3Mu = 0.d0

        B2CRMu3E = 0.d0
        B2CRTau3Mu = 0.d0
        B2CRTau3e = 0.d0

        B3CRMu3E = 0.d0
        B3cRTau3Mu = 0.d0
        B3cRTau3E = 0.d0

        loop226: do a = 1, 2
        loop227: do b = 1, 2
        loop228: do x = 1, 3
        loop229: do t = 1, 3

C       -----------------------------------------------------

        B1cRMu3E = B1CRMu3E + (1.d0/(8.d0*pi*alph))*CLlMU(a,x)*
     .  CLlE(a,t)*CLlE(b,t)*CLlE(b,x)*J4c(a,b,x,t)

        B1cRTAu3Mu = B1CRTau3Mu + (1.d0/(8.d0*pi*alph))*CLlTaU(a,x)*
     .  CLlMu(a,t)*CLlMu(b,t)*CLlMu(b,x)*J4c(a,b,x,t)

        B1cRTau3E = B1CRTau3E + (1.d0/(8.d0*pi*alph))*CLlTaU(a,x)*
     .  CLlE(a,t)*CLlE(b,t)*CLlE(b,x)*J4c(a,b,x,t)

C       -----------------------------------------------------

        B2cRMu3E = B2CRMu3E + (1.d0/(4.d0*pi*alph))*(0.25d0*
     .  CLlMU(a,x)*CLlE(a,t)*CRlE(b,t)*CRlE(b,x)*J4c(a,b,x,t) 
     .  - 0.5d0* Ceg(a)*Ceg(b)*CLlMu(a,x)*CRlE(a,t)*CRlE(b,t)*
     .  CLlE(b,x)*I4c(a,b,x,t))

        B2cRTAu3Mu = B2CRTau3Mu + (1.d0/(4.d0*pi*alph))*(0.25d0*
     .  CLlTaU(a,x)*CLlMu(a,t)*CRlMu(b,t)*CRlMu(b,x)*J4c(a,b,x,t) 
     .  - 0.5d0* Ceg(a)*Ceg(b)*CLlTau(a,x)*CRlMu(a,t)*CRlMu(b,t)*
     .  CLlMu(b,x)*I4c(a,b,x,t))

        B2cRTau3E = B2CRTau3E + (1.d0/(4.d0*pi*alph))*(0.25d0*
     .  CLlTaU(a,x)*CLlE(a,t)*CRlE(b,t)*CRlE(b,x)*J4c(a,b,x,t) 
     .  - 0.5d0* Ceg(a)*Ceg(b)*CLlTau(a,x)*CRlE(a,t)*CRlE(b,t)*
     .  CLlE(b,x)*I4c(a,b,x,t))

C       -----------------------------------------------------

        B3CRMu3E = B3CRMu3E + (1.d0/(4.d0*pi*alph))*Ceg(a)*Ceg(b)*
     .  CLlMu(a,x)*CRlE(a,t)*CLlE(b,t)*CRlE(b,x)*I4c(a,b,x,t)

       B3CRTau3Mu = B3CRTau3Mu + (1.d0/(4.d0*pi*alph))*Ceg(a)*Ceg(b)*
     .  CLlTau(a,x)*CRlMu(a,t)*CLlMu(b,t)*CRlMu(b,x)*I4c(a,b,x,t)

       B3CRTau3E = B3CRTau3E + (1.d0/(4.d0*pi*alph))*Ceg(a)*Ceg(b)*
     .  CLlTau(a,x)*CRlE(a,t)*CLlE(b,t)*CRlE(b,x)*I4c(a,b,x,t)

C       -----------------------------------------------------

        enddo loop229
        enddo loop228
        enddo loop227
        enddo loop226
        
        B4cRMu3E = 0.d0
        B4cRTau3E = 0.d0
        B4cRTau3Mu = 0.d0


C     ============================
C     Amplitudes and Rates for MU
C     ============================

C     ==========
C     =: Mu3E :=
C     ==========

        B1LMu3E = B1nLMu3E + B1cLMu3E
        B2LMu3E = B2nLMu3E + B2cLMu3E
        B3LMu3E = B3nLMu3E + B3cLMu3E
        B4LMu3E = B4nLMu3E + B4cLMu3E
        
        B1RMu3E = B1nRMu3E + B1cRMu3E
        B2RMu3E = B2nRMu3E + B2cRMu3E
        B3RMu3E = B3nRMu3E + B3cRMu3E
        B4RMu3E = B4nRMu3E + B4cRMu3E

C     ============
C     =: Tau3Mu :=
C     ============

        B1LTau3Mu = B1nLTau3Mu + B1cLTAu3Mu
        B2LTau3Mu = B2nLTau3Mu + B2cLTAu3Mu
        B3LTau3Mu = B3nLTau3Mu + B3cLTAu3Mu
        B4LTau3Mu = B4nLTau3Mu + B4cLTAu3Mu
        
        B1RTau3Mu = B1nRTau3Mu + B1cRTAu3Mu
        B2RTau3Mu = B2nRTau3Mu + B2cRTAu3Mu
        B3RTau3Mu = B3nRTau3Mu + B3cRTAu3Mu
        B4RTau3Mu = B4nRTau3Mu + B4cRTAu3Mu

C     ============
C     =: Tau3E :=
C     ============

        B1LTau3E = B1nLTau3E + B1cLTAu3E 
        B2LTau3E = B2nLTau3E + B2cLTAu3E
        B3LTau3E = B3nLTau3E + B3cLTAu3E
        B4LTau3E = B4nLTau3E + B4cLTAu3E

        B1RTau3E = B1nRTau3E + B1cRTAu3E
        B2RTau3E = B2nRTau3E + B2cRTAu3E
        B3RTau3E = B3nRTau3E + B3cRTAu3E
        B4RTau3E = B4nRTau3E + B4cRTAu3E
        
C	======================================================
C	Check this more properly...only preliminary!!!!!
C	======================================================

	FLLMu3E = (ZpengMuEEENL + ZpengMuEEEc)*(-0.5d0 + sinsqtw)/
     .	(mZ*mZ*stw*stw*ctw*ctw)
	FRRMu3E = (ZpengMuEEENR)*(sinsqtw)/(mZ*mZ*stw*stw*ctw*ctw)
	FLRMu3E = (ZpengMuEEENL + ZpengMuEEEc)*( sinsqtw)/
     .	(mZ*mZ*stw*stw*ctw*ctw)
	FRLMu3E = (ZpengMuEEENR )*(-0.5d0 + sinsqtw)/(mZ*mZ*stw*
     .	stw*ctw*ctw)

        FLLTau3E = (ZpengTauEEENL + ZpengTauEEEc)*(-0.5d0 + sinsqtw)/
     .	(mZ*mZ*stw*stw*ctw*ctw)
	FRRTau3E = (ZpengTauEEENR)*(sinsqtw)/(mZ*mZ*stw*stw*ctw*ctw)
	FLRTau3E = (ZpengTauEEENL + ZpengTauEEEc)*( sinsqtw)/
     .	(mZ*mZ*stw*stw*ctw*ctw)
	FRLTau3E = (ZpengTauEEENR )*(-0.5d0 + sinsqtw)/(mZ*mZ*stw*
     .	stw*ctw*ctw)

        FLLTau3Mu = (ZpengTauMuMuMuNL + ZpengTauMuMuMuc)*
     .  (-0.5d0 + sinsqtw)/(mZ*mZ*stw*stw*ctw*ctw)
	FRRTau3Mu = (ZpengTauMuMuMuNR)*(sinsqtw)/
     .   (mZ*mZ*stw*stw*ctw*ctw)
	FLRTau3Mu = (ZpengTauMuMuMuNL + ZpengTauMuMuMuc)*
     .   ( sinsqtw)/(mZ*mZ*stw*stw*ctw*ctw)
	FRLTau3Mu = (ZpengTauMuMuMuNR )*(-0.5d0 + sinsqtw)/
     .  (mZ*mZ*stw*stw*ctw*ctw)

        	
C	=========================================================

	VLMuE = VAN1LMUE + VAC1LMUE
	VRMuE = VAN1RMUE + VAC1RMUE
	ALMuE = ALN + ALC 
	ARMuE = ARN + ARC 

C      ============================================================
        
        VLTAUMU = VAN1LTAUMU + VAC1LTAUMU
        VRTAUMU = VAN1RTaUMU + VAC1RTAUMU
        ALTAUMU = TAUALN + TAUALC
        ARTAUMU = TAUARN + TAUARC

C     ===============================================================

        VLTAUE = VAN1LTAUE + VAC1LTAUE
        VRTAUE = VAN1RTaUE + VAC1RTAUE
        ALTAUE = TEALN + TEALC
        ARTAUE = TEARN + TEARC

C	======================================================
	Mu3Erate = 0.d0

	Mu3Erate = (alph**2.d0/(32.d0*pi))*(mMu**5.d0)*(VLMuE**2.d0 
     .	+ VRMuE**2.d0 - 4.d0*(VLMuE*ARMuE + VRMuE*ALMuE) + 
     . 	(ALMuE**2.d0 + ARMuE**2.d0)*((16.d0/3.d0)*dlog(mMU/(2.d0*mE))- 
     . 14.d0/9.d0) + (1.d0/6.d0)*(B1LMu3E**2.d0 + B1RMu3E**2.d0) +
     . (1.d0/3.d0)*(B2LMu3E**2.d0 + B2RMu3E**2.d0) + (1.d0/24.d0)* 
     . (B3LMu3E**2.d0 + B3RMu3E**2.d0) + 6.d0*(B4LMu3E**2.d0 + 
     . B4RMu3E**2.d0) - (B3LMu3E*B4LMu3E + B3RMu3E*B4RMu3E) + 
     . (2.d0/3.d0)*(VLMuE*B1LMu3E + VRMuE*B1RMu3E + VLMuE*B2LMu3E 
     . + VRMuE*B2RMu3E) - (4.d0/3.d0)*(ARMuE*B1LMu3E + ALMuE*B1RMu3E 
     . + ARMuE*B2LMu3E + ALMuE*B2RMu3E) + (1.d0/3.d0)*(2.d0*(
     . FLLMu3E**2.d0 + FRRMu3E**2.d0) + FLRMu3E**2.d0 + FRLMu3E**2.d0
     . + 2.d0*(B1LMu3E*FLLMu3E + B1RMu3E*FRRMu3E + B2LMu3E*FLRMu3E 
     . + B2RMu3E*FRLMu3E) + 4.d0*(VLMuE*FLLMu3E + VRMuE*FRRMu3E)
     . + 2.d0*(VLMuE*FLRMu3E + VRMuE*FRLMu3E) - 8.d0*(ARMuE*FLLMu3E 
     . + ALMuE*FRRMu3E) - 4.d0*(ALMuE*FRLMu3E + ARMuE*FLRMu3E)))

	Brmu3E = Mu3Erate/(3.d0*(10.d0**(-19.d0)))

C     ==================================================================
        TAu3MUrate = 0.d0

	TAu3MUrate = (alph**2.d0/(32.d0*pi))*(mTau**5.d0)*(VLTauMU**2.d0 
     .	+ VRTauMu**2.d0 - 4.d0*(VLTauMu*ARTauMu + VRTauMu*ALTauMu) + 
     . 	(ALTauMu**2.d0 + ARTauMu**2.d0)*((16.d0/3.d0)*
     .  dlog(mTau/(2.d0*mmu))- 14.d0/9.d0) + (1.d0/6.d0)*(B1LTau3Mu
     .  **2.d0 + B1RTau3Mu**2.d0) + (1.d0/3.d0)*(B2LTau3Mu**2.d0 + 
     .  B2RTau3Mu**2.d0) + (1.d0/24.d0)* (B3LTau3Mu**2.d0 +B3RTau3Mu
     .  **2.d0) + 6.d0*(B4LTau3Mu**2.d0 + B4RTau3Mu**2.d0) - 
     .  (B3LTau3Mu*B4LTau3Mu + B3RTau3Mu*B4RTau3Mu) + (2.d0/3.d0)*
     .   (VLTauMu*B1LTau3Mu + VRTauMu*B1RTau3Mu + VLTauMu*B2LTau3Mu 
     . + VRTauMu*B2RTau3Mu) - (4.d0/3.d0)*(ARTauMu*B1LTau3Mu + ALTauMu*
     .   B1RTau3Mu + ARTauMu*B2LTau3Mu + ALTauMu*B2RTau3Mu) + 
     .   (1.d0/3.d0)*(2.d0*(FLLTau3Mu**2.d0 + FRRTau3Mu**2.d0) + 
     .   FLRTau3Mu**2.d0 + FRLTau3Mu**2.d0 + 2.d0*(B1LTau3Mu*FLLTau3Mu + 
     .   B1RTau3Mu*FRRTau3Mu + B2LTau3Mu*FLRTau3Mu +B2RTau3Mu*FRLTau3Mu) 
     .   + 4.d0*(VLTauMu*FLLTau3Mu + VRTauMu*FRRTau3Mu) + 2.d0*(VLTauMu*
     .   FLRTau3Mu + VRTauMu*FRLTau3Mu) - 8.d0*(ARTauMu*FLLTau3Mu 
     . + ALTauMu*FRRTau3Mu) - 4.d0*(ALTauMu*FRLTau3Mu + 
     .   ARTauMu*FLRTau3Mu)))


	BrTau3Mu = Tau3Murate/(2.26d0*(10.d0**(-12.d0)))
C     =======================================================================

      TAu3Erate = 0.d0

	TAu3Erate = (alph**2.d0/(32.d0*pi))*(mTau**5.d0)*(VLTauE**2.d0 
     .	+ VRTauE**2.d0 - 4.d0*(VLTauE*ARTauE + VRTauE*ALTauE) + 
     . 	(ALTauE**2.d0 + ARTauE**2.d0)*((16.d0/3.d0)*
     .  dlog(mTau/(2.d0*mE))- 14.d0/9.d0) + (1.d0/6.d0)*(B1LTau3E
     .  **2.d0 + B1RTau3E**2.d0) + (1.d0/3.d0)*(B2LTau3E**2.d0 + 
     .  B2RTau3E**2.d0) + (1.d0/24.d0)* (B3LTau3E**2.d0 +B3RTau3E
     .  **2.d0) + 6.d0*(B4LTau3E**2.d0 + B4RTau3E**2.d0) - 
     .  (B3LTau3E*B4LTau3E + B3RTau3E*B4RTau3E) + (2.d0/3.d0)*
     .   (VLTauE*B1LTau3E + VRTauE*B1RTau3E + VLTauE*B2LTau3E 
     . + VRTauE*B2RTau3E) - (4.d0/3.d0)*(ARTauE*B1LTau3E + ALTauE*
     .   B1RTau3E + ARTauE*B2LTau3E + ALTauE*B2RTau3E) + 
     .   (1.d0/3.d0)*(2.d0*(FLLTau3E**2.d0 + FRRTau3E**2.d0) + 
     .   FLRTau3E**2.d0 + FRLTau3E**2.d0 + 2.d0*(B1LTau3E*FLLTau3E + 
     .   B1RTau3E*FRRTau3E + B2LTau3E*FLRTau3E + B2RTau3E*FRLTau3E) 
     .   + 4.d0*(VLTauE*FLLTau3E + VRTauE*FRRTau3E) + 2.d0*(VLTauE*
     .   FLRTau3E + VRTauE*FRLTau3E) - 8.d0*(ARTauE*FLLTau3E 
     . + ALTauE*FRRTau3E) - 4.d0*(ALTauE*FRLTau3E + 
     .   ARTauE*FLRTau3E)))


	BrTau3E = Tau3Erate/(2.26d0*(10.d0**(-12.d0)))

C     ========================================================================
C       Mu-E Conversion in Nuclei
C     ========================================================================
C       ===================================
C       I3n for Neutralinos/UP-type SQUARKS
C       ===================================

C        loopnu192: do a = 1, 4
C         loopnu193: do b = 1, 4
C           loopnu194: do x = 1, 6
C
C          ScqqI3nu(a,b,x) = 0.d0 
C
C        condnu909:If(a.ne.b.and.dabs(ScNeg(a)**2.d0-ScNeg(b)**2.d0).
C     .  gt.diff1.and.dabs(ScNeg(a)*ScNeg(a)-ScSUeg(x)).gt.diff1
C     .  .and.dabs(ScNeg(b)*ScNeg(b)-ScSUeg(x)).gt.diff1)then
C
C        ScqqI3nu(a,b,x) = (-2.d0*piconst)*(ScNeg(b)*ScNeg(b)*
C     . dLog(ScNeg(b)*ScNeg(b)/(ScNeg(a)*ScNeg(a)))/((ScNeg(b)*
C     . ScNeg(b) - ScSleg(x))*(ScNeg(b)*ScNeg(b)-ScNeg(a)*
C     . ScNeg(a)))+(ScSUeg(x)*dLog(ScSUeg(x)/(ScNeg(a)*ScNeg(a))))
C     . /((ScSUeg(x) - ScNeg(b)*ScNeg(b))*(ScSUeg(x) - ScNeg(a)*
C     . ScNeg(a))))  
C
C        elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*
C     .  ScNeg(b)).gt.diff1.and.dabs(ScNeg(a)*ScNeg(a)-
C     .  ScSUeg(x)).gt.diff1.and.dabs(ScNeg(b)*ScNeg(b)-ScSUeg(x))
C     .  .le.diff1)then condnu909
C
CC        ScqqI3nu(a,b,x) = (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*
CC     .  dLog(ScNeg(a)*ScNeg(a)/(ScNeg(b)*ScNeg(b))))/((ScNeg(a)
CC     .  *ScNeg(a) - ScNeg(b)*ScNeg(b))**2.d0) + 1.d0/(ScNeg(b)*
Cc     .  ScNeg(b)- ScNeg(a)*ScNeg(a)) )
C
C          ScqqI3nu(a,b,x) = I3f2d(ScNeg(a)*ScNeg(a),ScSUeg(x))
C
C        elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a) - ScNeg(b)*
C     .  ScNeg(b)).gt.diff1.and.dabs(ScNeg(b)*ScNeg(b)-ScSUeg(x))
C     .  .gt.diff1.and.dabs(ScNeg(a)*ScNeg(a)-ScSUeg(x)).le.
C     .  diff1)then condnu909
C
CC        ScqqI3nu(a,b,x) = (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*
Cc     .  dLog(ScNeg(b)*ScNeg(b)/(ScNeg(a)*ScNeg(a))))/((ScNeg(b)
Cc     .  *ScNeg(b)-ScNeg(a)*ScNeg(a))**2.d0) +1.d0/(ScNeg(a)*ScNeg(a)
CC     .  -ScNeg(b)*ScNeg(b)) )
C
C           ScqqI3nu(a,b,x) = I3f2d(ScNeg(b)*ScNeg(b),ScSUeg(x))
C
C       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
C     . .le.diff1).and.dabs(ScNeg(b)*ScNeg(b)-ScSUeg(x)).gt.diff1
C     . .and.dabs(ScNeg(a)*ScNeg(a)-ScSUeg(x)).gt.diff1)then condnu909
C
Cc        ScqqI3nu(a,b,x) = (-2.d0*piconst)*((ScSUeg(x)*dLog(ScSUeg(x)/
CC     . (ScNeg(a)*ScNeg(a))))/((ScSUeg(x)-ScSNeg(a)*ScSNeg(a))**2.d0)
CC     . + 1.d0/(ScNeg(a)*ScNeg(a)- ScSUeg(x)) )
C
C          ScqqI3nu(a,b,x) = I3f2d(ScSUeg(x),ScNeg(b)*ScNeg(b))
C         
C       else condnu909
C
C         ScqqI3nu(a,b,x) = -1.d0/(32.d0*ScNeg(b)*ScNeg(b)*pi*pi)
C
C         endif condnu909
C
C         enddo loopnu194
C         enddo loopnu193
C         enddo loopnu192
C
CC     =======================rescaling of the loop functions==========
C
C         lpresclI3nua: do a = 1, 4
C         lpresclI3nub: do b = 1, 4
C         lpresclI3nuc: do x = 1, 6
C
C         qqI3nu(a,b,x) = ScqqI3nu(a,b,x)/1000000.d0
C
C         enddo lpresclI3nuc
C         enddo lpresclI3nub
C         enddo lpresclI3nua
C
CC   =================================================================
CC     I3 Loop Functions for the Down-Type SQUARKS and Neutralinos 
CC     ==========================================================
C
C        loopnd192: do a = 1, 4
C         loopnd193: do b = 1, 4
C           loopnd194: do x = 1, 6
C
C          ScqqI3nd(a,b,x) = 0.d0 
C
C        condnd909:If(a.ne.b.and.dabs(ScNeg(a)**2.d0-ScNeg(b)**2.d0).
C     .  gt.diff1.and.dabs(ScNeg(a)*ScNeg(a)-ScSDeg(x)).gt.diff1
C     .  .and.dabs(ScNeg(b)*ScNeg(b)-ScSDeg(x)).gt.diff1)then
C
C        ScqqI3nd(a,b,x) = (-2.d0*piconst)*(ScNeg(b)*ScNeg(b)*
C     . dLog(ScNeg(b)*ScNeg(b)/(ScNeg(a)*ScNeg(a)))/((ScNeg(b)*
C     . ScNeg(b) - ScSDeg(x))*(ScNeg(b)*ScNeg(b)-ScNeg(a)*
C     . ScNeg(a)))+(ScSDeg(x)*dLog(ScSDeg(x)/(ScNeg(a)*ScNeg(a))))
C     . /((ScSDeg(x) - ScNeg(b)*ScNeg(b))*(ScSDeg(x) - ScNeg(a)*
C     . ScNeg(a))))  
C
C        elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*
C     .  ScNeg(b)).gt.diff1.and.dabs(ScNeg(a)*ScNeg(a)-
C     .  ScSDeg(x)).gt.diff1.and.dabs(ScNeg(b)*ScNeg(b)-ScSDeg(x))
C     .  .le.diff1)then condnd909
C
CC        ScqqI3nd(a,b,x) = (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*
CC     .  dLog(ScNeg(a)*ScNeg(a)/(ScNeg(b)*ScNeg(b))))/((ScNeg(a)
CC     .  *ScNeg(a) - ScNeg(b)*ScNeg(b))**2.d0) + 1.d0/(ScNeg(b)*
Cc     .  ScNeg(b)- ScNeg(a)*ScNeg(a)) )
C
C          ScqqI3nd(a,b,x) = I3f2d(ScNeg(a)*ScNeg(a),ScSDeg(x))
C
C        elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a) - ScNeg(b)*
C     .  ScNeg(b)).gt.diff1.and.dabs(ScNeg(b)*ScNeg(b)-ScSDeg(x))
C     .  .gt.diff1.and.dabs(ScNeg(a)*ScNeg(a)-ScSDeg(x)).le.
C     .  diff1)then condnd909
C
CC        ScqqI3nd(a,b,x) = (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*
Cc     .  dLog(ScNeg(b)*ScNeg(b)/(ScNeg(a)*ScNeg(a))))/((ScNeg(b)
Cc     .  *ScNeg(b)-ScNeg(a)*ScNeg(a))**2.d0) +1.d0/(ScNeg(a)*ScNeg(a)
CC     .  -ScNeg(b)*ScNeg(b)) )
C
C           ScqqI3nd(a,b,x) = I3f2d(ScNeg(b)*ScNeg(b),ScSDeg(x))
C
C       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
C     . .le.diff1).and.dabs(ScNeg(b)*ScNeg(b)-ScSDeg(x)).gt.diff1
C     . .and.dabs(ScNeg(a)*ScNeg(a)-ScSDeg(x)).gt.diff1)then condnd909
C
Cc        ScqqI3nd(a,b,x) = (-2.d0*piconst)*((ScSDeg(x)*dLog(ScSDeg(x)/
CC     . (ScNeg(a)*ScNeg(a))))/((ScSDeg(x)-ScSNeg(a)*ScSNeg(a))**2.d0)
CC     . + 1.d0/(ScNeg(a)*ScNeg(a)- ScSDeg(x)) )
C
C          ScqqI3nd(a,b,x) = I3f2d(ScSDeg(x),ScNeg(b)*ScNeg(b))
C         
C       else condnd909
C
C         ScqqI3nd(a,b,x) = -1.d0/(32.d0*ScNeg(b)*ScNeg(b)*pi*pi)
C
C         endif condnd909
C
C         enddo loopnd194
C         enddo loopnd193
C         enddo loopnd192
C
CC     =======================rescaling of the loop functions==========
C
C         lpresclI3nda: do a = 1, 4
C         lpresclI3ndb: do b = 1, 4
C         lpresclI3ndc: do x = 1, 6
C
C         qqI3nd(a,b,x) = ScqqI3nd(a,b,x)/1000000.d0
C
C         enddo lpresclI3ndc
C         enddo lpresclI3ndb
C         enddo lpresclI3nda
C
CC     ===================================================================
CC	I3c functions for UP-SQUARKS and Charginos
CC     ===================================================================
C
C        loopcu192: do a = 1, 2
C         loopcu193: do b = 1, 2
C           loopcu194: do x = 1, 6
C
C          ScqqI3cu(a,b,x) = 0.d0 
C
C        condcu909:If(a.ne.b.and.dabs(ScCeg(a)**2.d0-ScCeg(b)**2.d0).
C     .  gt.diff1.and.dabs(ScCeg(a)*ScCeg(a)-ScSUeg(x)).gt.diff1
C     .  .and.dabs(ScCeg(b)*ScCeg(b)-ScSUeg(x)).gt.diff1)then
C
C        ScqqI3cu(a,b,x) = (-2.d0*piconst)*(ScCeg(b)*ScCeg(b)*
C     . dLog(ScCeg(b)*ScCeg(b)/(ScCeg(a)*ScCeg(a)))/((ScCeg(b)*
C     . ScCeg(b) - ScSUeg(x))*(ScCeg(b)*ScCeg(b)-ScCeg(a)*
C     . ScCeg(a)))+(ScSUeg(x)*dLog(ScSUeg(x)/(ScCeg(a)*ScCeg(a))))
C     . /((ScSUeg(x) - ScCeg(b)*ScCeg(b))*(ScSUeg(x) - ScCeg(a)*
C     . ScCeg(a))))  
C
C        elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*
C     .  ScCeg(b)).gt.diff1.and.dabs(ScCeg(a)*ScCeg(a)-
C     .  ScSUeg(x)).gt.diff1.and.dabs(ScCeg(b)*ScCeg(b)-ScSUeg(x))
C     .  .le.diff1)then condcu909
C
CC        ScqqI3cu(a,b,x) = (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*
CC     .  dLog(ScCeg(a)*ScCeg(a)/(ScCeg(b)*ScCeg(b))))/((ScCeg(a)
CC     .  *ScCeg(a) - ScCeg(b)*ScCeg(b))**2.d0) + 1.d0/(ScCeg(b)*
Cc     .  ScCeg(b)- ScCeg(a)*ScCeg(a)) )
C
C          ScqqI3cu(a,b,x) = I3f2d(ScCeg(a)*ScCeg(a),ScSUeg(x))
C
C        elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a) - ScCeg(b)*
C     .  ScCeg(b)).gt.diff1.and.dabs(ScCeg(b)*ScCeg(b)-ScSUeg(x))
C     .  .gt.diff1.and.dabs(ScCeg(a)*ScCeg(a)-ScSUeg(x)).le.
C     .  diff1)then condcu909
C
CC        ScqqI3cu(a,b,x) = (-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*
Cc     .  dLog(ScCeg(b)*ScCeg(b)/(ScCeg(a)*ScCeg(a))))/((ScCeg(b)
Cc     .  *ScCeg(b)-ScCeg(a)*ScCeg(a))**2.d0) +1.d0/(ScCeg(a)*ScCeg(a)
CC     .  -ScCeg(b)*ScCeg(b)) )
C
C           ScqqI3cu(a,b,x) = I3f2d(ScCeg(b)*ScCeg(b),ScSUeg(x))
C
C       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
C     . .le.diff1).and.dabs(ScCeg(b)*ScCeg(b)-ScSUeg(x)).gt.diff1
C     . .and.dabs(ScCeg(a)*ScCeg(a)-ScSUeg(x)).gt.diff1)then condcu909
C
Cc        ScqqI3cu(a,b,x) = (-2.d0*piconst)*((ScSUeg(x)*dLog(ScSUeg(x)/
CC     . (ScCeg(a)*ScCeg(a))))/((ScSUeg(x)-ScSNeg(a)*ScSNeg(a))**2.d0)
CC     . + 1.d0/(ScCeg(a)*ScCeg(a)- ScSUeg(x)) )
C
C          ScqqI3cu(a,b,x) = I3f2d(ScSUeg(x),ScCeg(b)*ScCeg(b))
C         
C       else condcu909
C
C         ScqqI3cu(a,b,x) = -1.d0/(32.d0*ScCeg(b)*ScCeg(b)*pi*pi)
C
C         endif condcu909
C
C         enddo loopcu194
C         enddo loopcu193
C         enddo loopcu192
C
CC     =======================rescaling of the loop functions==========
C
C         lpresclI3cua: do a = 1, 2
CC         lpresclI3cub: do b = 1, 2
CC         lpresclI3cuc: do x = 1, 6
CC
C         qqI3cu(a,b,x) = ScqqI3cu(a,b,x)/1000000.d0
CC
C         enddo lpresclI3cuc
C         enddo lpresclI3cub
C         enddo lpresclI3cua
CC     ===================================================================
CC	I3c functions for DOWN-SQUARKS and Charginos
CC     ===================================================================
C
C        loopcd192: do a = 1, 2
C         loopcd193: do b = 1, 2
C           loopcd194: do x = 1, 6
C
C          ScqqI3cd(a,b,x) = 0.d0 
C
C        condcd909:If(a.ne.b.and.dabs(ScCeg(a)**2.d0-ScCeg(b)**2.d0).
C     .  gt.diff1.and.dabs(ScCeg(a)*ScCeg(a)-ScSDeg(x)).gt.diff1
C     .  .and.dabs(ScCeg(b)*ScCeg(b)-ScSDeg(x)).gt.diff1)then
C
C        ScqqI3cd(a,b,x) = (-2.d0*piconst)*(ScCeg(b)*ScCeg(b)*
CC     . dLog(ScCeg(b)*ScCeg(b)/(ScCeg(a)*ScCeg(a)))/((ScCeg(b)*
CC     . ScCeg(b) - ScSDeg(x))*(ScCeg(b)*ScCeg(b)-ScCeg(a)*
C     . ScCeg(a)))+(ScSDeg(x)*dLog(ScSDeg(x)/(ScCeg(a)*ScCeg(a))))
C     . /((ScSDeg(x) - ScCeg(b)*ScCeg(b))*(ScSDeg(x)-ScCeg(a)*
C     . ScCeg(a))))  
C
C        elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*
C     .  ScCeg(b)).gt.diff1.and.dabs(ScCeg(a)*ScCeg(a)-
C     .  ScSDeg(x)).gt.diff1.and.dabs(ScCeg(b)*ScCeg(b)-ScSDeg(x))
C     .  .le.diff1)then condcd909
C
CC        ScqqI3cd(a,b,x) = (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*
CC     .  dLog(ScCeg(a)*ScCeg(a)/(ScCeg(b)*ScCeg(b))))/((ScCeg(a)
C     .  *ScCeg(a) - ScCeg(b)*ScCeg(b))**2.d0) + 1.d0/(ScCeg(b)*
Cc     .  ScCeg(b)- ScCeg(a)*ScCeg(a)) )
CC
C          ScqqI3cd(a,b,x) = I3f2d(ScCeg(a)*ScCeg(a),ScSDeg(x))
C
C        elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a) - ScCeg(b)*
C     .  ScCeg(b)).gt.diff1.and.dabs(ScCeg(b)*ScCeg(b)-ScSDeg(x))
C     .  .gt.diff1.and.dabs(ScCeg(a)*ScCeg(a)-ScSDeg(x)).le.
C     .  diff1)then condcd909
C
CC        ScqqI3cd(a,b,x) = (-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*
Cc     .  dLog(ScCeg(b)*ScCeg(b)/(ScCeg(a)*ScCeg(a))))/((ScCeg(b)
Cc     .  *ScCeg(b)-ScCeg(a)*ScCeg(a))**2.d0) +1.d0/(ScCeg(a)*ScCeg(a)
CC     .  -ScCeg(b)*ScCeg(b)) )
C
CC           ScqqI3cd(a,b,x) = I3f2d(ScCeg(b)*ScCeg(b),ScSDeg(x))
CC
C       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
C     . .le.diff1).and.dabs(ScCeg(b)*ScCeg(b)-ScSDeg(x)).gt.diff1
C     . .and.dabs(ScCeg(a)*ScCeg(a)-ScSDeg(x)).gt.diff1)then condcd909
C
Cc        ScqqI3cd(a,b,x) = (-2.d0*piconst)*((ScSDeg(x)*dLog(ScSDeg(x)/
CC     . (ScCeg(a)*ScCeg(a))))/((ScSDeg(x)-ScSNeg(a)*ScSNeg(a))**2.d0)
CC     . + 1.d0/(ScCeg(a)*ScCeg(a)- ScSDeg(x)) )
C
C          ScqqI3cd(a,b,x) = I3f2d(ScSDeg(x),ScCeg(b)*ScCeg(b))
C         
CC       else condcd909
C
C         ScqqI3cd(a,b,x) = -1.d0/(32.d0*ScCeg(b)*ScCeg(b)*pi*pi)
C
C         endif condcd909
C
C         enddo loopcd194
C         enddo loopcd193
C         enddo loopcd192
C
CC     =======================rescaling of the loop functions==========
CC
C         lpresclI3cda: do a = 1, 2
C         lpresclI3cdb: do b = 1, 2
C         lpresclI3cdc: do x = 1, 6
C
C         qqI3cd(a,b,x) = ScqqI3cd(a,b,x)/1000000.d0
C
C         enddo lpresclI3cdc
C         enddo lpresclI3cdb
C         enddo lpresclI3cda
CC     ======================================================================
CC     ============================================================
C     I4 type functions for neutralinos and Sleptons and UP-SQUARKS
C     =============================================================
       lpnu198: do a = 1,4
       lpnu199: do b = 1,4
       lpnu200: do x = 1,6
       lpnu201: do t = 1,6

         ScqqI4nu(a,b,x,t)=0.d0
        
        cndi4nu911: If(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*
     . ScNeg(b)).gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSUeg(t)).gt.diff2
     . .and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(b)*ScNeg(b)-ScSUeg(t)).gt.diff2.and.
     . dabs(ScSLeg(x) -ScSUeg(t)).gt.diff2)then

        ScqqI4nu(a,b,x,t) = (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*
     . dLog(ScNeg(b)*ScNeg(b)/(ScNeg(a)*ScNeg(a))))/((ScNeg(b)*ScNeg(b)
     .  -ScNeg(a)*ScNeg(a))*(ScNeg(b)*ScNeg(b)-ScSLeg(x))*(ScNeg(b)*
     .  ScNeg(b)-ScSUeg(t))) + (ScSLeg(x)*dLog(ScSLeg(x)/(ScNeg(a)*
     .  ScNeg(a))))/((ScSLeg(x)-ScNeg(a)*ScNeg(a))*(ScSLeg(x)-ScNeg(b)
     .  *ScNeg(b))*(ScSLeg(x)-ScSUeg(t)) ) + (ScSUeg(t)*dLog(ScSUeg(t)
     . /(ScNeg(a)*ScNeg(a))))/((ScSUeg(t)-ScNeg(a)*ScNeg(a))*(ScSUeg(t)
     .  -ScNeg(b)*ScNeg(b))*(ScSUeg(t)-ScSLeg(x)))) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSUeg(t)).gt.diff2.and.dabs(ScSLeg(x)-ScSUeg(t))
     . .le.diff2) then cndi4nu911

       ScqqI4nu(a,b,x,t) = I4f3(ScNeg(a)*ScNeg(a),ScNeg(b)*ScNeg(b),
     . ScSUeg(t)) 

C       (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*
C     . dLog(ScNeg(a)*ScNeg(a)/ScSUeg(t)))/((ScNeg(a)*ScNeg(a)-
C     . ScNeg(b)*ScNeg(b))*(ScNeg(a)*ScNeg(a)-ScSUeg(t))**2.d0) + 
C     . (ScNeg(b)*ScNeg(b)*dLog(ScNeg(b)*ScNeg(b)/ScSUeg(t)))/((ScNeg(b)
C     . *ScNeg(b)-ScNeg(a)*ScNeg(a))*(ScNeg(b)*ScNeg(b)-ScSUeg(t))
C     . **2.d0) + 1.d0/( (ScSUeg(t)-ScNeg(a)*ScNeg(a))*(ScSUeg(t)-
C     . ScNeg(b)*ScNeg(b))) ) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     .  dabs(ScNeg(a)*ScNeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)
     . *ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScSLeg(x)
     . -ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-ScSUeg(t))
     . .le.diff2)then cndi4nu911

       ScqqI4nu(a,b,x,t) =I4f3(ScNeg(a)*ScNeg(a),ScSLeg(x),ScSUeg(t))

C       ScqqI4nu(a,b,x,t) = (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*dLog(
C     . ScNeg(a)*ScNeg(a)/ScSUeg(t)))/((ScNeg(a)*ScNeg(a)-ScSLeg(x))
C     . *(ScNeg(a)*ScNeg(a)-ScSUeg(t))**2.d0)+(ScSLeg(x)*dLog(ScSLeg(x)
C     . /ScSUeg(t)))/((ScSLeg(x)-ScNeg(a)*ScNeg(a))*(ScSLeg(x)-
C     . ScSUeg(t))**2.d0) + 1.d0/((ScSUeg(t)-ScNeg(a)*ScNeg(a))*
C     . (ScSUeg(t)-ScSLeg(x))) ) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSUeg(t)).gt.diff2.and.dabs(ScSLeg(x)-
     . ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).le.
     . diff2)then cndi4nu911 

C       ScqqI4nu(a,b,x,t) = (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*dLog(
C     . ScNeg(a)*ScNeg(a)/ScSLeg(x)))/((ScNeg(a)*ScNeg(a)-ScSUeg(t))*
C     . (ScNeg(a)*ScNeg(a)-ScSLeg(x))**2.d0)+(ScSUeg(t)*dLog(ScSUeg(t)
C     . /ScSLeg(x)))/((ScSUeg(t)-ScNeg(a)*ScNeg(a))*(ScSUeg(t)-ScSLeg(x))
C     . **2.d0) + 1.d0/((ScSLeg(x)-ScNeg(a)*ScNeg(a))*(ScSLeg(x)-
C     . ScSUeg(t))) ) 

       ScqqI4nu(a,b,x,t) = I4f3(ScNeg(a)*ScNeg(a),ScSUeg(t),ScSLeg(x))

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(b)*ScNeg(b)-ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)
     . *ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScSLeg(x)-
     . ScSUeg(t)).gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSUeg(t)).
     . le.diff2)then cndi4nu911 

       ScqqI4nu(a,b,x,t) = I4f3(ScNeg(b)*ScNeg(b),ScSLeg(x),ScSUeg(t))

C      (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*dLog(
C    . ScNeg(b)*ScNeg(b)/ScSUeg(t)))/((ScNeg(b)*ScNeg(b)-ScSLeg(x))
C    . *(ScNeg(b)*ScNeg(b)-ScSUeg(t))**2.d0)+(ScSLeg(x)*dLog(ScSLeg(x)
C    . /ScSUeg(t)))/((ScSLeg(x)-ScNeg(b)*ScNeg(b))*(ScSLeg(x)-
C    . ScSUeg(t))**2.d0)+1.d0/( (ScSUeg(t)-ScNeg(b)*ScNeg(b))*
C    . (ScSUeg(t)-ScSLeg(x))) ) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).gt.
     . diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSUeg(t)).gt.diff2.and.dabs(
     . ScNeg(b)*ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)
     . -ScSUeg(t)).gt.diff2.and.dabs(ScSLeg(x)-ScSUeg(t))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).le.diff2)
     . then cndi4nu911 

       ScqqI4nu(a,b,x,t) = I4f3(ScNeg(b)*ScNeg(b),ScSUeg(t),ScSLeg(x)) 

C       (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*dLog(
C     . ScNeg(b)*ScNeg(b)/ScSLeg(x)))/((ScNeg(b)*ScNeg(b)-ScSUeg(t))
C     . *(ScNeg(b)*ScNeg(b)-ScSLeg(x))**2.d0)+(ScSUeg(t)*dLog(ScSUeg(t)
C     . /ScSLeg(x)))/((ScSUeg(t)-ScNeg(b)*ScNeg(b))*(ScSUeg(t)-ScSLeg(x)
C     . )**2.d0) + 1.d0/((ScSLeg(x)-ScNeg(b)*ScNeg(b))*(ScSLeg(x)- 
C     .  ScSUeg(t))) ) 

       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .le.diff2).and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSLeg(x)).gt.diff2.and.dabs(ScSLeg(x)-ScSUeg(t))
     . .gt.diff2)then cndi4nu911

       ScqqI4nu(a,b,x,t) = I4f3(ScSLeg(x),ScSUeg(t),ScNeg(b)*ScNeg(b))

C       (-2.d0*piconst)*((ScSLeg(x)*dLog(ScSLeg(x)/
C     . (ScNeg(b)*ScNeg(b))))/(((ScSLeg(x)-ScNeg(b)*ScNeg(b))**2.d0)*
C     . (ScSLeg(x)-ScSUeg(t))) + (ScSUeg(t)*dLog(ScSUeg(t)/(ScNeg(b)*
C     . ScNeg(b))))/(((ScSUeg(t)-ScNeg(b)*ScNeg(b))**2.d0)*(ScSUeg(t)- 
C     . ScSLeg(x))) + 1.d0/( (ScNeg(b)*ScNeg(b)-ScSLeg(x))*(ScNeg(b)*
C     . ScNeg(b)-ScSUeg(t))) )

       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).le.
     . diff2).and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).le.diff2.and.
     . dabs(ScSLeg(x)-ScSUeg(t)).gt.diff2)then cndi4nu911

C       ScqqI4nu(a,b,x,t)=(-2.d0*piconst)*((ScSUeg(t)*dLog(ScSUeg(t)/
C     . ScSLeg(x)))/((ScSUeg(t)-ScSLeg(x))**3.d0)-(ScSUeg(t) + 
C     . ScSLeg(x))/(2.d0*ScSLeg(x)*(ScSLeg(x)-ScSUeg(t))**2.d0 ))

         ScqqI4nu(a,b,x,t) = I4f3d1(ScSUeg(t),ScSLeg(x))

       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).le.
     . diff2).and.dabs(ScNeg(b)*ScNeg(b)-ScSUeg(t)).le.diff2.and.
     . dabs(ScSLeg(x)-ScSUeg(t)).gt.diff2)then cndi4nu911

       ScqqI4nu(a,b,x,t) = I4f3d1(ScSLeg(x),ScSUeg(t))

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).gt.
     . diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).le.diff2.and.
     . dabs(ScSLeg(x)-ScSUeg(t)).le.diff2)then cndi4nu911

       ScqqI4nu(a,b,x,t) =  I4f3d1(ScNeg(b)*ScNeg(b),ScSUeg(t))

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).gt.
     . diff2.and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).le.diff2.and.
     . dabs(ScSLeg(x)-ScSUeg(t)).le.diff2)then cndi4nu911

       ScqqI4nu(a,b,x,t) = I4f3d1(ScNeg(a)*ScNeg(a),ScSUeg(t))
 
      elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .le.diff2).and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSLeg(x)).gt.diff2.and.dabs(ScSLeg(x)-ScSUeg(t))
     . .le.diff2)then cndi4nu911
       
         ScqqI4nu(a,b,x,t) = I4f2d2d(ScNeg(b)*ScNeg(b),ScSUeg(t))
         
      elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).le.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSUeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSLeg(x)).le.diff2.and.dabs(ScSLeg(x)-ScSUeg(t))
     . .gt.diff2)then cndi4nu911
  
         ScqqI4nu(a,b,x,t) = I4f2d2d(ScSUeg(t),ScSLeg(x))
  
      elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSUeg(t)).le.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSUeg(t)).le.diff2.and.dabs(ScSLeg(x)-ScSUeg(t))
     . .gt.diff2)then cndi4nu911

       ScqqI4nu(a,b,x,t) = I4f2d2d(ScSLeg(x),ScSUeg(t))

       else cndi4nu911
        
       ScqqI4nu(a,b,x,t) = 1.d0/(96.d0*ScSLeg(x)*ScSLeg(x)*pi*pi)

        endif cndi4nu911 

       enddo lpnu201
       enddo lpnu200
       enddo lpnu199
       enddo lpnu198


c     ===============================================================

        lpnua: do a = 1, 4
         lpnub: do b = 1, 4
         lpnuc: do x = 1, 6
         lpnud: do t = 1,6

         qqI4nu(a,b,x,t) = ScqqI4nu(a,b,x,t)*(10.d0**(-12.d0))

         enddo lpnud
         enddo lpnuc
         enddo lpnub
         enddo lpnua

C     ==================================================================
C     =================================================
C     I4 type functions for neutralinos ,sleptons and DOWN-SQUARKS
C     =================================================
       lpnd198: do a = 1,4
       lpnd199: do b = 1,4
       lpnd200: do x = 1,6
       lpnd201: do t = 1,6

         ScqqI4nd(a,b,x,t)=0.d0
        
        cndi4nd911: If(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*
     . ScNeg(b)).gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSDeg(t)).gt.diff2
     . .and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(b)*ScNeg(b)-ScSDeg(t)).gt.diff2.and.
     . dabs(ScSLeg(x) -ScSDeg(t)).gt.diff2)then

        ScqqI4nd(a,b,x,t) = (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*
     . dLog(ScNeg(b)*ScNeg(b)/(ScNeg(a)*ScNeg(a))))/((ScNeg(b)*ScNeg(b)
     .  -ScNeg(a)*ScNeg(a))*(ScNeg(b)*ScNeg(b)-ScSLeg(x))*(ScNeg(b)*
     .  ScNeg(b)-ScSDeg(t))) + (ScSLeg(x)*dLog(ScSLeg(x)/(ScNeg(a)*
     .  ScNeg(a))))/((ScSLeg(x)-ScNeg(a)*ScNeg(a))*(ScSLeg(x)-ScNeg(b)
     .  *ScNeg(b))*(ScSLeg(x)-ScSDeg(t)) ) + (ScSDeg(t)*dLog(ScSDeg(t)
     . /(ScNeg(a)*ScNeg(a))))/((ScSDeg(t)-ScNeg(a)*ScNeg(a))*(ScSDeg(t)
     .  -ScNeg(b)*ScNeg(b))*(ScSDeg(t)-ScSLeg(x)))) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSDeg(t)).gt.diff2.and.dabs(ScSLeg(x)-ScSDeg(t))
     . .le.diff2) then cndi4nd911

       ScqqI4nd(a,b,x,t) = I4f3(ScNeg(a)*ScNeg(a),ScNeg(b)*ScNeg(b),
     . ScSDeg(t)) 

C       (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*
C     . dLog(ScNeg(a)*ScNeg(a)/ScSDeg(t)))/((ScNeg(a)*ScNeg(a)-
C     . ScNeg(b)*ScNeg(b))*(ScNeg(a)*ScNeg(a)-ScSDeg(t))**2.d0) + 
C     . (ScNeg(b)*ScNeg(b)*dLog(ScNeg(b)*ScNeg(b)/ScSDeg(t)))/((ScNeg(b)
C     . *ScNeg(b)-ScNeg(a)*ScNeg(a))*(ScNeg(b)*ScNeg(b)-ScSDeg(t))
C     . **2.d0) + 1.d0/( (ScSDeg(t)-ScNeg(a)*ScNeg(a))*(ScSDeg(t)-
C     . ScNeg(b)*ScNeg(b))) ) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     .  dabs(ScNeg(a)*ScNeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)
     . *ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScSLeg(x)
     . -ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-ScSDeg(t))
     . .le.diff2)then cndi4nd911

       ScqqI4nd(a,b,x,t) =I4f3(ScNeg(a)*ScNeg(a),ScSLeg(x),ScSDeg(t))

C       ScqqI4nd(a,b,x,t) = (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*dLog(
C     . ScNeg(a)*ScNeg(a)/ScSDeg(t)))/((ScNeg(a)*ScNeg(a)-ScSLeg(x))
C     . *(ScNeg(a)*ScNeg(a)-ScSDeg(t))**2.d0)+(ScSLeg(x)*dLog(ScSLeg(x)
C     . /ScSDeg(t)))/((ScSLeg(x)-ScNeg(a)*ScNeg(a))*(ScSLeg(x)-
C     . ScSDeg(t))**2.d0) + 1.d0/((ScSDeg(t)-ScNeg(a)*ScNeg(a))*
C     . (ScSDeg(t)-ScSLeg(x))) ) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSDeg(t)).gt.diff2.and.dabs(ScSLeg(x)-
     . ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).le.
     . diff2)then cndi4nd911 

C       ScqqI4nd(a,b,x,t) = (-2.d0*piconst)*((ScNeg(a)*ScNeg(a)*dLog(
C     . ScNeg(a)*ScNeg(a)/ScSLeg(x)))/((ScNeg(a)*ScNeg(a)-ScSDeg(t))*
C     . (ScNeg(a)*ScNeg(a)-ScSLeg(x))**2.d0)+(ScSDeg(t)*dLog(ScSDeg(t)
C     . /ScSLeg(x)))/((ScSDeg(t)-ScNeg(a)*ScNeg(a))*(ScSDeg(t)-ScSLeg(x))
C     . **2.d0) + 1.d0/((ScSLeg(x)-ScNeg(a)*ScNeg(a))*(ScSLeg(x)-
C     . ScSDeg(t))) ) 

       ScqqI4nd(a,b,x,t) = I4f3(ScNeg(a)*ScNeg(a),ScSDeg(t),ScSLeg(x))

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(b)*ScNeg(b)-ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)
     . *ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScSLeg(x)-
     . ScSDeg(t)).gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSDeg(t)).
     . le.diff2)then cndi4nd911 

       ScqqI4nd(a,b,x,t) = I4f3(ScNeg(b)*ScNeg(b),ScSLeg(x),ScSDeg(t))

C      (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*dLog(
C    . ScNeg(b)*ScNeg(b)/ScSDeg(t)))/((ScNeg(b)*ScNeg(b)-ScSLeg(x))
C    . *(ScNeg(b)*ScNeg(b)-ScSDeg(t))**2.d0)+(ScSLeg(x)*dLog(ScSLeg(x)
C    . /ScSDeg(t)))/((ScSLeg(x)-ScNeg(b)*ScNeg(b))*(ScSLeg(x)-
C    . ScSDeg(t))**2.d0)+1.d0/( (ScSDeg(t)-ScNeg(b)*ScNeg(b))*
C    . (ScSDeg(t)-ScSLeg(x))) ) 

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).gt.
     . diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSDeg(t)).gt.diff2.and.dabs(
     . ScNeg(b)*ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)
     . -ScSDeg(t)).gt.diff2.and.dabs(ScSLeg(x)-ScSDeg(t))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).le.diff2)
     . then cndi4nd911 

       ScqqI4nd(a,b,x,t) = I4f3(ScNeg(b)*ScNeg(b),ScSDeg(t),ScSLeg(x)) 

C       (-2.d0*piconst)*((ScNeg(b)*ScNeg(b)*dLog(
C     . ScNeg(b)*ScNeg(b)/ScSLeg(x)))/((ScNeg(b)*ScNeg(b)-ScSDeg(t))
C     . *(ScNeg(b)*ScNeg(b)-ScSLeg(x))**2.d0)+(ScSDeg(t)*dLog(ScSDeg(t)
C     . /ScSLeg(x)))/((ScSDeg(t)-ScNeg(b)*ScNeg(b))*(ScSDeg(t)-ScSLeg(x)
C     . )**2.d0) + 1.d0/((ScSLeg(x)-ScNeg(b)*ScNeg(b))*(ScSLeg(x)- 
C     .  ScSDeg(t))) ) 

       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .le.diff2).and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSLeg(x)).gt.diff2.and.dabs(ScSLeg(x)-ScSDeg(t))
     . .gt.diff2)then cndi4nd911

       ScqqI4nd(a,b,x,t) = I4f3(ScSLeg(x),ScSDeg(t),ScNeg(b)*ScNeg(b))

C       (-2.d0*piconst)*((ScSLeg(x)*dLog(ScSLeg(x)/
C     . (ScNeg(b)*ScNeg(b))))/(((ScSLeg(x)-ScNeg(b)*ScNeg(b))**2.d0)*
C     . (ScSLeg(x)-ScSDeg(t))) + (ScSDeg(t)*dLog(ScSDeg(t)/(ScNeg(b)*
C     . ScNeg(b))))/(((ScSDeg(t)-ScNeg(b)*ScNeg(b))**2.d0)*(ScSDeg(t)- 
C     . ScSLeg(x))) + 1.d0/( (ScNeg(b)*ScNeg(b)-ScSLeg(x))*(ScNeg(b)*
C     . ScNeg(b)-ScSDeg(t))) )

       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).le.
     . diff2).and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).le.diff2.and.
     . dabs(ScSLeg(x)-ScSDeg(t)).gt.diff2)then cndi4nd911

C       ScqqI4nd(a,b,x,t)=(-2.d0*piconst)*((ScSDeg(t)*dLog(ScSDeg(t)/
C     . ScSLeg(x)))/((ScSDeg(t)-ScSLeg(x))**3.d0)-(ScSDeg(t) + 
C     . ScSLeg(x))/(2.d0*ScSLeg(x)*(ScSLeg(x)-ScSDeg(t))**2.d0 ))

         ScqqI4nd(a,b,x,t) = I4f3d1(ScSDeg(t),ScSLeg(x))

       elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).le.
     . diff2).and.dabs(ScNeg(b)*ScNeg(b)-ScSDeg(t)).le.diff2.and.
     . dabs(ScSLeg(x)-ScSDeg(t)).gt.diff2)then cndi4nd911

       ScqqI4nd(a,b,x,t) = I4f3d1(ScSLeg(x),ScSDeg(t))

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).gt.
     . diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).le.diff2.and.
     . dabs(ScSLeg(x)-ScSDeg(t)).le.diff2)then cndi4nd911

       ScqqI4nd(a,b,x,t) =  I4f3d1(ScNeg(b)*ScNeg(b),ScSDeg(t))

       elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b)).gt.
     . diff2.and.dabs(ScNeg(b)*ScNeg(b)-ScSLeg(x)).le.diff2.and.
     . dabs(ScSLeg(x)-ScSDeg(t)).le.diff2)then cndi4nd911

       ScqqI4nd(a,b,x,t) = I4f3d1(ScNeg(a)*ScNeg(a),ScSDeg(t))
 
      elseif((a.eq.b.or.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .le.diff2).and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSLeg(x)).gt.diff2.and.dabs(ScSLeg(x)-ScSDeg(t))
     . .le.diff2)then cndi4nd911
       
         ScqqI4nd(a,b,x,t) = I4f2d2d(ScNeg(b)*ScNeg(b),ScSDeg(t))
         
      elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).le.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSDeg(t)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSLeg(x)).le.diff2.and.dabs(ScSLeg(x)-ScSDeg(t))
     . .gt.diff2)then cndi4nd911
  
         ScqqI4nd(a,b,x,t) = I4f2d2d(ScSDeg(t),ScSLeg(x))
  
      elseif(a.ne.b.and.dabs(ScNeg(a)*ScNeg(a)-ScNeg(b)*ScNeg(b))
     . .gt.diff2.and.dabs(ScNeg(a)*ScNeg(a)-ScSDeg(t)).le.diff2.and.
     . dabs(ScNeg(a)*ScNeg(a)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*
     . ScNeg(b)-ScSLeg(x)).gt.diff2.and.dabs(ScNeg(b)*ScNeg(b)-
     . ScSDeg(t)).le.diff2.and.dabs(ScSLeg(x)-ScSDeg(t))
     . .gt.diff2)then cndi4nd911

       ScqqI4nd(a,b,x,t) = I4f2d2d(ScSLeg(x),ScSDeg(t))

       else cndi4nd911
        
       ScqqI4nd(a,b,x,t) = 1.d0/(96.d0*ScSLeg(x)*ScSLeg(x)*pi*pi)

        endif cndi4nd911 

       enddo lpnd201
       enddo lpnd200
       enddo lpnd199
       enddo lpnd198


c     ===============================================================

        lpnda: do a = 1, 4
         lpndb: do b = 1, 4
         lpndc: do x = 1, 6
         lpndd: do t = 1,6

         qqI4nd(a,b,x,t) = ScqqI4nd(a,b,x,t)*(10.d0**(-12.d0))

         enddo lpndd
         enddo lpndc
         enddo lpndb
         enddo lpnda

C     ==================================================================
C     =================================================
C     I4 type functions for charginos, Sneutrinos and UP-SQUARKS
C     =================================================
       lpcu198: do a = 1,2
       lpcu199: do b = 1,2
       lpcu200: do x = 1,3
       lpcu201: do t = 1,6

         ScqqI4cu(a,b,x,t)=0.d0
        
        cndi4cu911: If(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*
     . ScCeg(b)).gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSUeg(t)).gt.diff2
     . .and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(b)*ScCeg(b)-ScSUeg(t)).gt.diff2.and.
     . dabs(ScSNeg(x) -ScSUeg(t)).gt.diff2)then

        ScqqI4cu(a,b,x,t) = (-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*
     . dLog(ScCeg(b)*ScCeg(b)/(ScCeg(a)*ScCeg(a))))/((ScCeg(b)*ScCeg(b)
     .  -ScCeg(a)*ScCeg(a))*(ScCeg(b)*ScCeg(b)-ScSNeg(x))*(ScCeg(b)*
     .  ScCeg(b)-ScSUeg(t))) + (ScSNeg(x)*dLog(ScSNeg(x)/(ScCeg(a)*
     .  ScCeg(a))))/((ScSNeg(x)-ScCeg(a)*ScCeg(a))*(ScSNeg(x)-ScCeg(b)
     .  *ScCeg(b))*(ScSNeg(x)-ScSUeg(t)) ) + (ScSUeg(t)*dLog(ScSUeg(t)
     . /(ScCeg(a)*ScCeg(a))))/((ScSUeg(t)-ScCeg(a)*ScCeg(a))*(ScSUeg(t)
     .  -ScCeg(b)*ScCeg(b))*(ScSUeg(t)-ScSNeg(x)))) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSUeg(t)).gt.diff2.and.dabs(ScSNeg(x)-ScSUeg(t))
     . .le.diff2) then cndi4cu911

       ScqqI4cu(a,b,x,t) = I4f3(ScCeg(a)*ScCeg(a),ScCeg(b)*ScCeg(b),
     . ScSUeg(t)) 

C       (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*
C     . dLog(ScCeg(a)*ScCeg(a)/ScSUeg(t)))/((ScCeg(a)*ScCeg(a)-
C     . ScCeg(b)*ScCeg(b))*(ScCeg(a)*ScCeg(a)-ScSUeg(t))**2.d0) + 
C     . (ScCeg(b)*ScCeg(b)*dLog(ScCeg(b)*ScCeg(b)/ScSUeg(t)))/((ScCeg(b)
C     . *ScCeg(b)-ScCeg(a)*ScCeg(a))*(ScCeg(b)*ScCeg(b)-ScSUeg(t))
C     . **2.d0) + 1.d0/( (ScSUeg(t)-ScCeg(a)*ScCeg(a))*(ScSUeg(t)-
C     . ScCeg(b)*ScCeg(b))) ) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     .  dabs(ScCeg(a)*ScCeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)
     . *ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScSNeg(x)
     . -ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-ScSUeg(t))
     . .le.diff2)then cndi4cu911

       ScqqI4cu(a,b,x,t) =I4f3(ScCeg(a)*ScCeg(a),ScSNeg(x),ScSUeg(t))

C       ScqqI4cu(a,b,x,t) = (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*dLog(
C     . ScCeg(a)*ScCeg(a)/ScSUeg(t)))/((ScCeg(a)*ScCeg(a)-ScSNeg(x))
C     . *(ScCeg(a)*ScCeg(a)-ScSUeg(t))**2.d0)+(ScSNeg(x)*dLog(ScSNeg(x)
C     . /ScSUeg(t)))/((ScSNeg(x)-ScCeg(a)*ScCeg(a))*(ScSNeg(x)-
C     . ScSUeg(t))**2.d0) + 1.d0/((ScSUeg(t)-ScCeg(a)*ScCeg(a))*
C     . (ScSUeg(t)-ScSNeg(x))) ) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSUeg(t)).gt.diff2.and.dabs(ScSNeg(x)-
     . ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).le.
     . diff2)then cndi4cu911 

C       ScqqI4cu(a,b,x,t) = (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*dLog(
C     . ScCeg(a)*ScCeg(a)/ScSNeg(x)))/((ScCeg(a)*ScCeg(a)-ScSUeg(t))*
C     . (ScCeg(a)*ScCeg(a)-ScSNeg(x))**2.d0)+(ScSUeg(t)*dLog(ScSUeg(t)
C     . /ScSNeg(x)))/((ScSUeg(t)-ScCeg(a)*ScCeg(a))*(ScSUeg(t)-ScSNeg(x))
C     . **2.d0) + 1.d0/((ScSNeg(x)-ScCeg(a)*ScCeg(a))*(ScSNeg(x)-
C     . ScSUeg(t))) ) 

       ScqqI4cu(a,b,x,t) = I4f3(ScCeg(a)*ScCeg(a),ScSUeg(t),ScSNeg(x))

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(b)*ScCeg(b)-ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)
     . *ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(SLeg(x)-
     . SLeg(t)).gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSUeg(t)).
     . le.diff2)then cndi4cu911 

       ScqqI4cu(a,b,x,t) = I4f3(ScCeg(b)*ScCeg(b),ScSNeg(x),ScSUeg(t))

C      (-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*dLog(
C    . ScCeg(b)*ScCeg(b)/ScSUeg(t)))/((ScCeg(b)*ScCeg(b)-ScSNeg(x))
C    . *(ScCeg(b)*ScCeg(b)-ScSUeg(t))**2.d0)+(ScSNeg(x)*dLog(ScSNeg(x)
C    . /ScSUeg(t)))/((ScSNeg(x)-ScCeg(b)*ScCeg(b))*(ScSNeg(x)-
C    . ScSUeg(t))**2.d0)+1.d0/( (ScSUeg(t)-ScCeg(b)*ScCeg(b))*
C    . (ScSUeg(t)-ScSNeg(x))) ) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).gt.
     . diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSUeg(t)).gt.diff2.and.dabs(
     . ScCeg(b)*ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)
     . -ScSUeg(t)).gt.diff2.and.dabs(ScSNeg(x)-ScSUeg(t))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).le.diff2)
     . then cndi4cu911 

       ScqqI4cu(a,b,x,t) = I4f3(ScCeg(b)*ScCeg(b),ScSUeg(t),ScSNeg(x)) 

C       (-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*dLog(
C     . ScCeg(b)*ScCeg(b)/ScSNeg(x)))/((ScCeg(b)*ScCeg(b)-ScSUeg(t))
C     . *(ScCeg(b)*ScCeg(b)-ScSNeg(x))**2.d0)+(ScSUeg(t)*dLog(ScSUeg(t)
C     . /ScSNeg(x)))/((ScSUeg(t)-ScCeg(b)*ScCeg(b))*(ScSUeg(t)-ScSNeg(x)
C     . )**2.d0) + 1.d0/((ScSNeg(x)-ScCeg(b)*ScCeg(b))*(ScSNeg(x)- 
C     .  ScSUeg(t))) ) 

       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .le.diff2).and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSNeg(x)).gt.diff2.and.dabs(ScSNeg(x)-ScSUeg(t))
     . .gt.diff2)then cndi4cu911

       ScqqI4cu(a,b,x,t) = I4f3(ScSNeg(x),ScSUeg(t),ScCeg(b)*ScCeg(b))

C       (-2.d0*piconst)*((ScSNeg(x)*dLog(ScSNeg(x)/
C     . (ScCeg(b)*ScCeg(b))))/(((ScSNeg(x)-ScCeg(b)*ScCeg(b))**2.d0)*
C     . (ScSNeg(x)-ScSUeg(t))) + (ScSUeg(t)*dLog(ScSUeg(t)/(ScCeg(b)*
C     . ScCeg(b))))/(((ScSUeg(t)-ScCeg(b)*ScCeg(b))**2.d0)*(ScSUeg(t)- 
C     . ScSNeg(x))) + 1.d0/( (ScCeg(b)*ScCeg(b)-ScSNeg(x))*(ScCeg(b)*
C     . ScCeg(b)-ScSUeg(t))) )

       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).le.
     . diff2).and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSUeg(t)).gt.diff2)then cndi4cu911

C       ScqqI4cu(a,b,x,t)=(-2.d0*piconst)*((ScSUeg(t)*dLog(ScSUeg(t)/
C     . ScSNeg(x)))/((ScSUeg(t)-ScSNeg(x))**3.d0)-(ScSUeg(t) + 
C     . ScSNeg(x))/(2.d0*ScSNeg(x)*(ScSNeg(x)-ScSUeg(t))**2.d0 ))

         ScqqI4cu(a,b,x,t) = I4f3d1(ScSUeg(t),ScSNeg(x))

       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).le.
     . diff2).and.dabs(ScCeg(b)*ScCeg(b)-ScSUeg(t)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSUeg(t)).gt.diff2)then cndi4cu911

       ScqqI4cu(a,b,x,t) = I4f3d1(ScSNeg(x),ScSUeg(t))

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).gt.
     . diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSUeg(t)).le.diff2)then cndi4cu911

       ScqqI4cu(a,b,x,t) =  I4f3d1(ScCeg(b)*ScCeg(b),ScSUeg(t))

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).gt.
     . diff2.and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSUeg(t)).le.diff2)then cndi4cu911

       ScqqI4cu(a,b,x,t) = I4f3d1(ScCeg(a)*ScCeg(a),ScSUeg(t))
 
      elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .le.diff2).and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSNeg(x)).gt.diff2.and.dabs(ScSNeg(x)-ScSUeg(t))
     . .le.diff2)then cndi4cu911
       
         ScqqI4cu(a,b,x,t) = I4f2d2d(ScCeg(b)*ScCeg(b),ScSUeg(t))
         
      elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).le.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSUeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSNeg(x)).le.diff2.and.dabs(ScSNeg(x)-ScSUeg(t))
     . .gt.diff2)then cndi4cu911
  
         ScqqI4cu(a,b,x,t) = I4f2d2d(ScSUeg(t),ScSNeg(x))
  
      elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSUeg(t)).le.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSUeg(t)).le.diff2.and.dabs(ScSNeg(x)-ScSUeg(t))
     . .gt.diff2)then cndi4cu911

       ScqqI4cu(a,b,x,t) = I4f2d2d(ScSNeg(x),ScSUeg(t))

       else cndi4cu911
        
       ScqqI4cu(a,b,x,t) = 1.d0/(96.d0*ScSNeg(x)*ScSNeg(x)*pi*pi)

        endif cndi4cu911 

       enddo lpcu201
       enddo lpcu200
       enddo lpcu199
       enddo lpcu198
c     ===============================================================

        lpcua: do a = 1,2
         lpcub: do b = 1,2
         lpcuc: do x = 1,3
         lpcud: do t = 1,6

         qqI4cu(a,b,x,t) = ScqqI4cu(a,b,x,t)*(10.d0**(-12.d0))

         enddo lpcud
         enddo lpcuc
         enddo lpcub
         enddo lpcua

C	======================================================================
C     =================================================
C     I4 type functions for charginos, sneutrinos and DOWN-SQUARKS
C     =================================================
       lpcd198: do a = 1,2
       lpcd199: do b = 1,2
       lpcd200: do x = 1,3
       lpcd201: do t = 1,6

         ScqqI4cd(a,b,x,t)=0.d0
        
        cndi4cd911: If(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*
     . ScCeg(b)).gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSDeg(t)).gt.diff2
     . .and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(b)*ScCeg(b)-ScSDeg(t)).gt.diff2.and.
     . dabs(ScSNeg(x) -ScSDeg(t)).gt.diff2)then

        ScqqI4cd(a,b,x,t) = (-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*
     . dLog(ScCeg(b)*ScCeg(b)/(ScCeg(a)*ScCeg(a))))/((ScCeg(b)*ScCeg(b)
     .  -ScCeg(a)*ScCeg(a))*(ScCeg(b)*ScCeg(b)-ScSNeg(x))*(ScCeg(b)*
     .  ScCeg(b)-ScSDeg(t))) + (ScSNeg(x)*dLog(ScSNeg(x)/(ScCeg(a)*
     .  ScCeg(a))))/((ScSNeg(x)-ScCeg(a)*ScCeg(a))*(ScSNeg(x)-ScCeg(b)
     .  *ScCeg(b))*(ScSNeg(x)-ScSDeg(t)) ) + (ScSDeg(t)*dLog(ScSDeg(t)
     . /(ScCeg(a)*ScCeg(a))))/((ScSDeg(t)-ScCeg(a)*ScCeg(a))*(ScSDeg(t)
     .  -ScCeg(b)*ScCeg(b))*(ScSDeg(t)-ScSNeg(x)))) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSDeg(t)).gt.diff2.and.dabs(ScSNeg(x)-ScSDeg(t))
     . .le.diff2) then cndi4cd911

       ScqqI4cd(a,b,x,t) = I4f3(ScCeg(a)*ScCeg(a),ScCeg(b)*ScCeg(b),
     . ScSDeg(t)) 

C       (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*
C     . dLog(ScCeg(a)*ScCeg(a)/ScSDeg(t)))/((ScCeg(a)*ScCeg(a)-
C     . ScCeg(b)*ScCeg(b))*(ScCeg(a)*ScCeg(a)-ScSDeg(t))**2.d0) + 
C     . (ScCeg(b)*ScCeg(b)*dLog(ScCeg(b)*ScCeg(b)/ScSDeg(t)))/((ScCeg(b)
C     . *ScCeg(b)-ScCeg(a)*ScCeg(a))*(ScCeg(b)*ScCeg(b)-ScSDeg(t))
C     . **2.d0) + 1.d0/( (ScSDeg(t)-ScCeg(a)*ScCeg(a))*(ScSDeg(t)-
C     . ScCeg(b)*ScCeg(b))) ) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     .  dabs(ScCeg(a)*ScCeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)
     . *ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScSNeg(x)
     . -ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-ScSDeg(t))
     . .le.diff2)then cndi4cd911

       ScqqI4cd(a,b,x,t) =I4f3(ScCeg(a)*ScCeg(a),ScSNeg(x),ScSDeg(t))

C       ScqqI4cd(a,b,x,t) = (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*dLog(
C     . ScCeg(a)*ScCeg(a)/ScSDeg(t)))/((ScCeg(a)*ScCeg(a)-ScSNeg(x))
C     . *(ScCeg(a)*ScCeg(a)-ScSDeg(t))**2.d0)+(ScSNeg(x)*dLog(ScSNeg(x)
C     . /ScSDeg(t)))/((ScSNeg(x)-ScCeg(a)*ScCeg(a))*(ScSNeg(x)-
C     . ScSDeg(t))**2.d0) + 1.d0/((ScSDeg(t)-ScCeg(a)*ScCeg(a))*
C     . (ScSDeg(t)-ScSNeg(x))) ) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSDeg(t)).gt.diff2.and.dabs(ScSNeg(x)-
     . ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).le.
     . diff2)then cndi4cd911 

C       ScqqI4cd(a,b,x,t) = (-2.d0*piconst)*((ScCeg(a)*ScCeg(a)*dLog(
C     . ScCeg(a)*ScCeg(a)/ScSNeg(x)))/((ScCeg(a)*ScCeg(a)-ScSDeg(t))*
C     . (ScCeg(a)*ScCeg(a)-ScSNeg(x))**2.d0)+(ScSDeg(t)*dLog(ScSDeg(t)
C     . /ScSNeg(x)))/((ScSDeg(t)-ScCeg(a)*ScCeg(a))*(ScSDeg(t)-ScSNeg(x))
C     . **2.d0) + 1.d0/((ScSNeg(x)-ScCeg(a)*ScCeg(a))*(ScSNeg(x)-
C     . ScSDeg(t))) ) 

       ScqqI4cd(a,b,x,t) = I4f3(ScCeg(a)*ScCeg(a),ScSDeg(t),ScSNeg(x))

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(b)*ScCeg(b)-ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)
     . *ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(SLeg(x)-
     . SLeg(t)).gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSDeg(t)).
     . le.diff2)then cndi4cd911 

       ScqqI4cd(a,b,x,t) = I4f3(ScCeg(b)*ScCeg(b),ScSNeg(x),ScSDeg(t))

C      (-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*dLog(
C    . ScCeg(b)*ScCeg(b)/ScSDeg(t)))/((ScCeg(b)*ScCeg(b)-ScSNeg(x))
C    . *(ScCeg(b)*ScCeg(b)-ScSDeg(t))**2.d0)+(ScSNeg(x)*dLog(ScSNeg(x)
C    . /ScSDeg(t)))/((ScSNeg(x)-ScCeg(b)*ScCeg(b))*(ScSNeg(x)-
C    . ScSDeg(t))**2.d0)+1.d0/( (ScSDeg(t)-ScCeg(b)*ScCeg(b))*
C    . (ScSDeg(t)-ScSNeg(x))) ) 

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).gt.
     . diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSDeg(t)).gt.diff2.and.dabs(
     . ScCeg(b)*ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)
     . -ScSDeg(t)).gt.diff2.and.dabs(ScSNeg(x)-ScSDeg(t))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).le.diff2)
     . then cndi4cd911 

       ScqqI4cd(a,b,x,t) = I4f3(ScCeg(b)*ScCeg(b),ScSDeg(t),ScSNeg(x)) 

C       (-2.d0*piconst)*((ScCeg(b)*ScCeg(b)*dLog(
C     . ScCeg(b)*ScCeg(b)/ScSNeg(x)))/((ScCeg(b)*ScCeg(b)-ScSDeg(t))
C     . *(ScCeg(b)*ScCeg(b)-ScSNeg(x))**2.d0)+(ScSDeg(t)*dLog(ScSDeg(t)
C     . /ScSNeg(x)))/((ScSDeg(t)-ScCeg(b)*ScCeg(b))*(ScSDeg(t)-ScSNeg(x)
C     . )**2.d0) + 1.d0/((ScSNeg(x)-ScCeg(b)*ScCeg(b))*(ScSNeg(x)- 
C     .  ScSDeg(t))) ) 

       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .le.diff2).and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSNeg(x)).gt.diff2.and.dabs(ScSNeg(x)-ScSDeg(t))
     . .gt.diff2)then cndi4cd911

       ScqqI4cd(a,b,x,t) = I4f3(ScSNeg(x),ScSDeg(t),ScCeg(b)*ScCeg(b))

C       (-2.d0*piconst)*((ScSNeg(x)*dLog(ScSNeg(x)/
C     . (ScCeg(b)*ScCeg(b))))/(((ScSNeg(x)-ScCeg(b)*ScCeg(b))**2.d0)*
C     . (ScSNeg(x)-ScSDeg(t))) + (ScSDeg(t)*dLog(ScSDeg(t)/(ScCeg(b)*
C     . ScCeg(b))))/(((ScSDeg(t)-ScCeg(b)*ScCeg(b))**2.d0)*(ScSDeg(t)- 
C     . ScSNeg(x))) + 1.d0/( (ScCeg(b)*ScCeg(b)-ScSNeg(x))*(ScCeg(b)*
C     . ScCeg(b)-ScSDeg(t))) )

       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).le.
     . diff2).and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSDeg(t)).gt.diff2)then cndi4cd911

C       ScqqI4cd(a,b,x,t)=(-2.d0*piconst)*((ScSDeg(t)*dLog(ScSDeg(t)/
C     . ScSNeg(x)))/((ScSDeg(t)-ScSNeg(x))**3.d0)-(ScSDeg(t) + 
C     . ScSNeg(x))/(2.d0*ScSNeg(x)*(ScSNeg(x)-ScSDeg(t))**2.d0 ))

         ScqqI4cd(a,b,x,t) = I4f3d1(ScSDeg(t),ScSNeg(x))

       elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).le.
     . diff2).and.dabs(ScCeg(b)*ScCeg(b)-ScSDeg(t)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSDeg(t)).gt.diff2)then cndi4cd911

       ScqqI4cd(a,b,x,t) = I4f3d1(ScSNeg(x),ScSDeg(t))

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).gt.
     . diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSDeg(t)).le.diff2)then cndi4cd911

       ScqqI4cd(a,b,x,t) =  I4f3d1(ScCeg(b)*ScCeg(b),ScSDeg(t))

       elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b)).gt.
     . diff2.and.dabs(ScCeg(b)*ScCeg(b)-ScSNeg(x)).le.diff2.and.
     . dabs(ScSNeg(x)-ScSDeg(t)).le.diff2)then cndi4cd911

       ScqqI4cd(a,b,x,t) = I4f3d1(ScCeg(a)*ScCeg(a),ScSDeg(t))
 
      elseif((a.eq.b.or.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .le.diff2).and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSNeg(x)).gt.diff2.and.dabs(ScSNeg(x)-ScSDeg(t))
     . .le.diff2)then cndi4cd911
       
         ScqqI4cd(a,b,x,t) = I4f2d2d(ScCeg(b)*ScCeg(b),ScSDeg(t))
         
      elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).le.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSDeg(t)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSNeg(x)).le.diff2.and.dabs(ScSNeg(x)-ScSDeg(t))
     . .gt.diff2)then cndi4cd911
  
         ScqqI4cd(a,b,x,t) = I4f2d2d(ScSDeg(t),ScSNeg(x))
  
      elseif(a.ne.b.and.dabs(ScCeg(a)*ScCeg(a)-ScCeg(b)*ScCeg(b))
     . .gt.diff2.and.dabs(ScCeg(a)*ScCeg(a)-ScSDeg(t)).le.diff2.and.
     . dabs(ScCeg(a)*ScCeg(a)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*
     . ScCeg(b)-ScSNeg(x)).gt.diff2.and.dabs(ScCeg(b)*ScCeg(b)-
     . ScSDeg(t)).le.diff2.and.dabs(ScSNeg(x)-ScSDeg(t))
     . .gt.diff2)then cndi4cd911

       ScqqI4cd(a,b,x,t) = I4f2d2d(ScSNeg(x),ScSDeg(t))

       else cndi4cd911
        
       ScqqI4cd(a,b,x,t) = 1.d0/(96.d0*ScSNeg(x)*ScSNeg(x)*pi*pi)

        endif cndi4cd911 

       enddo lpcd201
       enddo lpcd200
       enddo lpcd199
       enddo lpcd198
c     ===============================================================

        lpcda: do a = 1,2
         lpcdb: do b = 1,2
         lpcdc: do x = 1,3
         lpcdd: do t = 1,6

         qqI4cd(a,b,x,t) = ScqqI4cd(a,b,x,t)*(10.d0**(-12.d0))

         enddo lpcdd
         enddo lpcdc
         enddo lpcdb
         enddo lpcda
C    ================================================================
C     J4-type functions for up-type squarks and neutralinos 
C     ============================================================

      lpj4nu6: do a = 1,4 
      lpj4nu7: do b = 1,4
      lpj4nu8: do x = 1, 6
      lpj4nu9: do t = 1, 6

       qqJ4nu(a,b,x,t) = 0.d0

        qqJ4nu(a,b,x,t) = SUeg(t)*qqI4nu(a,b,x,t) + I3n(a,b,x)

          enddo lpj4nu9
          enddo lpj4nu8
          enddo lpj4nu7
          enddo lpj4nu6

C     ===========================================================
C     J4-type functions for Down-type squarks and neutralinos 
C     ===========================================================

      lpj4nd6: do a = 1,4 
      lpj4nd7: do b = 1,4
      lpj4nd8: do x = 1, 6
      lpj4nd9: do t = 1, 6

       qqJ4nd(a,b,x,t) = 0.d0

        qqJ4nd(a,b,x,t) = SDeg(t)*qqI4nd(a,b,x,t) + I3n(a,b,x)

          enddo lpj4nd9
          enddo lpj4nd8
          enddo lpj4nd7
          enddo lpj4nd6
C     ============================================================
C     J4-type functions for Up-type Squarks and Charginos
C     =========================================================

      lpj4cu6: do a = 1,2
      lpj4cu7: do b = 1,2
      lpj4cu8: do x = 1, 3
      lpj4cu9: do t = 1, 6

       qqJ4cu(a,b,x,t) = 0.d0

        qqJ4cu(a,b,x,t) = SUeg(t)*qqI4cu(a,b,x,t) + I3c(a,b,x)

          enddo lpj4cu9
          enddo lpj4cu8
          enddo lpj4cu7
          enddo lpj4cu6
C     =========================================================
C     J4-type functions for Down-type Squarks and Charginos
C     =========================================================

      lpj4cd6: do a = 1,2
      lpj4cd7: do b = 1,2
      lpj4cd8: do x = 1, 3
      lpj4cd9: do t = 1, 6

       qqJ4cd(a,b,x,t) = 0.d0

        qqJ4cd(a,b,x,t) = SDeg(t)*qqI4cd(a,b,x,t) + I3c(a,b,x)

          enddo lpj4cd9
          enddo lpj4cd8
          enddo lpj4cd7
          enddo lpj4cd6

C     ===================================================================
C       Up and Down Quark Couplings. 
C       ---------------------------

C       With Charginos. 
C       -------------- 

        intlpx1: Do a = 1,2 
        intlpx2: do x = 1,6

        CRdD(a,x) = 0.d0
        CLdD(a,x) = 0.d0
        CRuU(a,x) = 0.d0
        CLuU(a,x) = 0.d0

        enddo intlpx2
        enddo intlpx1


        lpx1: Do a = 1,2 
        lpx2: do x = 1,6
        
        CRdD(a,x) = g2*(-OCR(a,1)*USU(x,1) + mUQ*
     .       (1.d0/(dtan(beta)*cdenc))*OCR(a,2)*USU(x,4) ) 
        
        CLdD(a,x) = g2*md/cdenc*OCL(a,2)*USU(x,1) 
        
        CRuU(a,x) = g2*(- OCL(a,1)*USD(x,1) + 
     .       md/cdenc*OCL(a,2)*USD(x,4))
        
        CLuU(a,x) = g2*muq*(1.d0/(dtan(beta)*cdenc))*
     .       OCR(a,2)*USD(x,1)
        
        enddo lpx2 
        enddo lpx1

C       --------------------------------------------       
C       (* Neutralino-QUARK-sQUARK Couplings *)
C       --------------------------------------------       

        intlpnx1: Do a = 1,4 
        intlpnx2: do x = 1,6

        NRdD(a,x) = 0.d0
        NLdD(a,x) = 0.d0
        NRuU(a,x) = 0.d0
        NLuU(a,x) = 0.d0

        enddo intlpnx2
        enddo intlpnx1



        lpnx1: do a = 1, 4
        lpnx2: do x = 1, 6

        NRdD(a,x) = -g2t*((- ON(a,2) + (1.d0/3.d0)*ON(a,1)*
     .  dTan(tw))*USD(x,1) + md/cdenn*ON(a,3)*USD(x,4) )

        NLdD(a,x) = -g2t*(md/cdenn*ON(a,3)*USD(x,1) 
     .	 + (2.d0/3.d0)*dTan(tw)*ON(a,1)*USD(x,4))

        NRuU(a,x) =  -g2t*((ON(a,2) + (1.d0/3.d0)*ON(a,1)*
     .  dTan(tw))*USU(x,1) + muq*(1.d0/(dtan(beta)*cdenn))*
     .  ON(a,4)*USU(x,4))

        NLuU(a,x) = -g2t*(muq*(1.d0/(dtan(beta)*cdenn))*ON(a,4)*
     .  USU(x,1) - (4.d0/3.d0)*dTan(tw)*ON(a,1)*USU(x,4))

         enddo lpnx2 
         enddo lpnx1 

C    =================================================================
C     Box Contributions from Neutralinos
C     ================================================================

         DnLuMuE = 0.d0
         DnLdMuE = 0.d0
         DnRuMuE = 0.d0
         DnRdMuE = 0.d0

         ddlp1: do a = 1,4
            ddlp2: do b = 1, 4
               ddlp3: do x = 1, 6
                  ddlp4: do t = 1, 6

        DnLuMuE = DnLuMuE + (1.d0/(16.d0*pi*alph))*(0.5d0*(NRlMu(a,x)*
     .  NRlE(b,x)*NRuU(a,t)*NRuU(b,t) - NRlMu(a,x)*NRlE(b,x)*NLuU(a,t)
     .  *NLuU(b,t))*qqJ4nu(a,b,x,t) - Neg(a)*Neg(b)*qqI4nu(a,b,x,t)*
     .  (NRlMu(a,x)*NRlE(b,x)*NLuU(a,t)*NLuU(b,t) - NRlMu(a,x)*
     .  NRlE(b,x)*NRuU(a,t)*NRuU(b,t)))

        DnLdMuE = DnLdMuE + (1.d0/(16.d0*pi*alph))*(0.5d0*(NRlMu(a,x)*
     .  NRlE(b,x)*NRdD(a,t)*NRdD(b,t) - NRlMu(a,x)*NRlE(b,x)*NLdD(a,t)
     .  *NLdD(b,t))*qqJ4nd(a,b,x,t) - Neg(a)*Neg(b)*qqI4nd(a,b,x,t)*
     .  (NRlMu(a,x)*NRlE(b,x)*NLdD(a,t)*NLdD(b,t) - NRlMu(a,x)*
     .  NRlE(b,x)*NRdD(a,t)*NRdD(b,t)))
 
         DnRuMuE = DnRuMuE + (1.d0/(16.d0*pi*alph))*(0.5d0*(NLlMu(a,x)*
     .   NLlE(b,x)*NLuU(a,t)*NLuU(b,t) - NLlMu(a,x)*NLlE(b,x)*NRuU(a,t)
     .   *NRuU(b,t))*qqJ4nu(a,b,x,t) - Neg(a)*Neg(b)*qqI4nu(a,b,x,t)*
     .   (NLlMu(a,x)*NLlE(b,x)*NRuU(a,t)*NRuU(b,t) - NLlMu(a,x)*
     .    NLlE(b,x)*NLuU(a,t)*NLuU(b,t)))

        DnRdMuE = DnRdMuE + (1.d0/(16.d0*pi*alph))*(0.5d0*(NLlMu(a,x)*
     .  NLlE(b,x)*NLdD(a,t)*NLdD(b,t) - NLlMu(a,x)*NLlE(b,x)*NRdD(a,t)
     .  *NRdD(b,t))*qqJ4nd(a,b,x,t) - Neg(a)*Neg(b)*qqI4nd(a,b,x,t)*
     .  (NLlMu(a,x)*NLlE(b,x)*NRdD(a,t)*NRdD(b,t) - NLlMu(a,x)*
     .   NLlE(b,x)*NLdD(a,t)*NLdD(b,t))) 


          enddo ddlp4
          enddo ddlp3
          enddo ddlp2
          enddo ddlp1


C     =================================================================
C     Contributions From Charginos
C     ================================================================

          DcLdMuE = 0.d0
          DcLuMuE = 0.d0
	  DcRdMuE = 0.d0
	  DcRuMuE = 0.d0
                    
         dclp1:  do a = 1, 2
           dclp2:  do b = 1, 2
            dclp3:    do x = 1, 3
              dclp4:   do t = 1, 6

                      
C         DcLdMuE = DcLdMuE + (1.d0/(16.d0*pi*alph))*(0.5d0*CRlMu(a,x)*
C     .  CRlE(b,x)*CRdD(a,t)*CRdD(b,t)*qqJ4cd(a,b,x,t) - Ceg(a)*Ceg(b)*
C     .  qqI4cd(a,b,x,t)*CRlMu(a,x)*CRlE(b,x)*CLdD(a,t)*CLdD(b,t))

         DcLdMuE = DcLdMuE + (1.d0/(16.d0*pi*alph))*(0.5d0*CRlMu(a,x)*
     .  CRlE(b,x)*CRdD(a,t)*CRdD(b,t)*qqJ4cu(a,b,x,t) - Ceg(a)*Ceg(b)*
     .  qqI4cu(a,b,x,t)*CRlMu(a,x)*CRlE(b,x)*CLdD(a,t)*CLdD(b,t))

C         DcLuMuE = DcLuMuE - (1.d0/(16.d0*pi*alph))*(0.5d0*CRlMu(a,x)*
C     .  CRlE(b,x)*CLuU(a,t)*CLuU(b,t)*qqJ4cu(a,b,x,t) - Ceg(a)*Ceg(b)* 
C     .  qqI4cu(a,b,x,t)*CRlMu(a,x)*CRlE(b,x)*CRuU(a,t)*CRuU(b,t))

         DcLuMuE = DcLuMuE - (1.d0/(16.d0*pi*alph))*(0.5d0*CRlMu(a,x)*
     .  CRlE(b,x)*CLuU(a,t)*CLuU(b,t)*qqJ4cd(a,b,x,t) - Ceg(a)*Ceg(b)* 
     .  qqI4cd(a,b,x,t)*CRlMu(a,x)*CRlE(b,x)*CRuU(a,t)*CRuU(b,t))

C         DcRdMuE = DcRdMuE + (1.d0/(16.d0*pi*alph))*(0.5d0*CLlMu(a,x)*
C     .  CLlE(b,x)*CLdD(a,t)*CLdD(b,t)*qqJ4cd(a,b,x,t) - Ceg(a)*Ceg(b)*
C     .  qqI4cd(a,b,x,t)*CLlMu(a,x)*CLlE(b,x)*CRdD(a,t)*CRdD(b,t))

         DcRdMuE = DcRdMuE + (1.d0/(16.d0*pi*alph))*(0.5d0*CLlMu(a,x)*
     .  CLlE(b,x)*CLdD(a,t)*CLdD(b,t)*qqJ4cu(a,b,x,t) - Ceg(a)*Ceg(b)*
     .  qqI4cu(a,b,x,t)*CLlMu(a,x)*CLlE(b,x)*CRdD(a,t)*CRdD(b,t))

C         DcRuMuE = DcRuMuE - (1.d0/(16.d0*pi*alph))*(0.5d0*CLlMu(a,x)*
C     .  CLlE(b,x)*CRuU(a,t)*CRuU(b,t)*qqJ4cu(a,b,x,t) - Ceg(a)*Ceg(b)* 
C     .  qqI4cu(a,b,x,t)*CLlMu(a,x)*CLlE(b,x)*CLuU(a,t)*CLuU(b,t))

         DcRuMuE = DcRuMuE - (1.d0/(16.d0*pi*alph))*(0.5d0*CLlMu(a,x)*
     .  CLlE(b,x)*CRuU(a,t)*CRuU(b,t)*qqJ4cd(a,b,x,t) - Ceg(a)*Ceg(b)* 
     .  qqI4cd(a,b,x,t)*CLlMu(a,x)*CLlE(b,x)*CLuU(a,t)*CLuU(b,t))

          enddo dclp4
          enddo dclp3
          enddo dclp2
          enddo dclp1
         
C     ========================================================================

          DbarLuMuE =0.d0
          DbarLdMuE =0.d0
          DbarRuMuE =0.d0
          DbarRdMuE =0.d0

          DbarLuMuE = DcLuMuE + DnLuMuE + ((0.5d0 - (4.d0/3.d0)*
     .    sinsqtw)/2.d0)*((ZpengMuEEENL + ZpengMuEEEc)/(mZ*mZ*stw*stw*
     .    ctw*ctw))

          DbarLdMuE = DcLdMuE + DnLdMuE + ((-0.5d0 + (2.d0/3.d0)*
     .    sinsqtw)/2.d0)*((ZpengMuEEENL + ZpengMuEEEc)/(mZ*mZ*stw*stw*
     .    ctw*ctw))

          DbarRuMuE = DcRuMuE + DnRuMuE + ((0.5d0 - (4.d0/3.d0)*
     .    sinsqtw)/2.d0)*((ZpengMuEEENR)/(mZ*mZ*stw*stw*
     .    ctw*ctw))

          DbarRdMuE = DcRdMuE + DnRdMuE + ((-0.5d0 + (2.d0/3.d0)*
     .    sinsqtw)/2.d0)*((ZpengMuEEENR)/(mZ*mZ*stw*stw*
     .    ctw*ctw))

C     ===========================================================
C     Some numbers for Titanium as per Bernabeau, Nardi et.al, check also
C     Okadas numbers
C     ============================================================

        mueconverrate = 0.d0
        Brmueconver = 0.d0

          Zti = 22.d0 ! should be integers ? 
          Nti = 26.d0 ! should be integers ?
          Zeffti = 17.6d0
          Fqsqti = 0.54d0


          mueconverrate =  4.d0*(alph**5.d0)*((Zeffti**4.d0)/Zti)*
     .  Fqsqti*(mmu**5.d0)*( (Zti*(VLMuE - ARMuE) - (2.d0*Zti +
     .  Nti)*DbarLuMuE - (Zti + 2.d0*Nti)*DbarLdMuE)**2.d0 + 
     .  (Zti*(VRMuE - ALMuE) - (2.d0*Zti + Nti)*DbarRuMuE - 
     .  (Zti + 2.d0*Nti)*DbarRdMuE)**2.d0)
       
          Brmueconver = mueconverrate/(3.d0*(10.d0**(-19.d0)))



C     ====================================================================================
C	!!!!!!!!!!!!!! G-2 !!!!	
C     ====================================================================================
C       Neutralino Amplitude (L)
C       ------------------------

	loopgi33: do a = 1, 4
	 loopgi34:  do x = 1, 6
	Ampg2LN(a,x) = 0.d0
                enddo loopgi34
                enddo loopgi33

C       -----------------------------------------------

	loopg133: do a = 1, 4
	   loopg134: do x = 1, 6

	Ampg2LN(a,x) = -(mMU**2.d0)*(4.d0*piconst*(NLlMU(a,x)*NLlMU(a,x)
     .	*f1(a,x)+ 2.d0*NLlMU(a,x)*NRlMU(a,x)*(Neg(a)/mMU)*f2(a,x))*
     .	(1.d0/SLeg(x)))

	enddo loopg134
	enddo loopg133

	
C       Sum Neutralino Amplitude (L)
C       ---------------------------

	g2ALN = 0.d0
	
	loopg135: do a = 1, 4
	loopg136: do x = 1, 6

	g2ALN = g2ALN + Ampg2LN(a,x)
	
	enddo loopg136
	enddo loopg135

C       Neutralino Amplitude (R) 
C       ------------------------
	loopgi37: do a = 1, 4
	 loopgi38:  do x = 1, 6
	Ampg2RN(a,x) = 0.d0
                enddo loopgi38             
                enddo loopgi37
C       ------------------------

	loopg137: do a = 1, 4
	loopg138: do x = 1, 6

	Ampg2RN(a,x)= -(mMU**2.d0)*piconst*(4.d0*NRlMU(a,x)*NRlMU(a,x)
     .	*f1(a,x)+ 2.d0*NRlMU(a,x)*NLlMU(a,x)*(Neg(a)/mMU)*f2(a,x))*
     .  (1.d0/SLeg(x))

	enddo loopg138
	enddo loopg137

C       Sum Neutralinos (R)
C       ---------------------

	g2ARN = 0.d0
	
	loopg139: do a = 1, 4
	loopg140: do x = 1, 6

	g2ARN = g2ARN + Ampg2RN(a,x)
	
	enddo loopg140
	enddo loopg139


C       Chargino Amplitudes (L)
C       -------------------------------------

        loopgi41:do a = 1,2
	   loopgi42: do x = 1, 3 
	Ampg2LC(a,x) = 0.d0
                enddo loopgi42
                enddo loopgi41

C       -------------------------------------

	loopg141: do a = 1,2
	 loopg142:   do x = 1, 3 

        Ampg2LC(a,x) = mMU*mMU*piconst*(1/SNeg(x))*(4.d0*CLlMU(a,x)
     .	*CLlMU(a,x)*f3(a,x) + 2.d0*CLlMU(a,x)*CRlMU(a,x)*(Ceg(a)/mMU)*
     .  f4(a,x))

	enddo loopg142
	enddo loopg141

	g2ALC = 0.d0

C	Sum Chargino Amplitude (L)
C       --------------------------

	loopg143: do a = 1,2
	 loopg144:  do x = 1, 3
	      
	      g2ALC = g2ALC + Ampg2LC(a,x)

	      enddo loopg144
	      enddo loopg143


C       Chargino Amplitudes (R)
C       ---------------------------------------------
	      
	loopgi45: do a = 1,2
	loopgi46: do x = 1,3 

	Ampg2RC(a,x) = 0.d0

              enddo loopgi46
              enddo loopgi45

C       ---------------------------------------------

	loopg145: do a = 1,2
	loopg146: do x = 1,3 

	Ampg2RC(a,x) = mMU*mMU*piconst*(1/SNeg(x))*(4.d0*CRlMU(a,x)
     .	*CRlMU(a,x)*f3(a,x) + 2.d0*CRlMU(a,x)*CLlMU(a,x)*(Ceg(a)/mMU)
     .       *f4(a,x))
        
      enddo loopg146
      enddo loopg145
        
	g2ARC = 0.d0
        
C	Sum Chargino Amplitude (R)
C     ---------------------------
        
	loopg147: do a = 1,2
	loopg148: do x = 1, 3
        
        g2ARC = g2ARC + Ampg2RC(a,x)
        
      enddo loopg148
      enddo loopg147
      
              
C     ------------------------------------------------
C     G-2  is!!
C     ------------------------------------------------
      gminus2 = g2ALC + g2ALN  + g2ARC + g2ARN
      
C     ======================================================================================
      
      
! 456  return
      return
      end subroutine muegamma
      
C     ==============================================================
C     Limit function for I4(m1,m2,m3,m4) called I4f3(x1,x2,x3).
C     Notation x1,x2 are the two non-degenerate mass SQUAREDs where as
C     x3 denotes the degenerate mass.
C     ==============================================================

        double precision function I4f3(x1,x2,x3)
        implicit none
        double precision x1,x2,x3,pi

        pi = 3.14159265358979d0

        I4f3 = -( 1.d0/(16.d0*pi*pi) )*((x1*dlog(x1/x3))/
     .  ((x1-x2)*(x1-x3)**2.d0) + (x2*dlog(x2/x3))/((x2-x1)*
     .  ((x2-x3)**2.d0)) + 1.d0/((x3-x1)*(x3-x2)))
        
        end function  
 
C     ==============================================================
C     Limit function for I4(m1,m2,m3,m4) called I4f3d1(x1,x2).
C     Notation x1,x2 is such that x1 represents the non-degenerate
C     mass SQUARED where as x2 denotes the degenerate mass SQUARED 
C     of three particles
C     ==============================================================

        double precision function I4f3d1(x1,x2)
        implicit none
        double precision x1,x2,pi

        pi = 3.14159265358979d0


        I4f3d1 = -( 1.d0/(16.d0*pi*pi) )*((x1*dlog(x1/x2))/
     .  ((x1-x2)**3.d0) - (x1 + x2)/(2.d0*x2*(x2-x1)**2.d0))
        
        end function  

C     ==============================================================
C     Limit function for I4(m1,m2,m3,m4) called I4f2d2d(x1,x2).
C     Notation x1,x2 is such that x1 represents one pari of degenerate
C     mass SQUARED where as x2 denotes the other pair of degenerate mass 
C     SQUARED 
C     ==============================================================

        double precision function I4f2d2d(x1,x2)
        implicit none
        double precision x1,x2,pi

        pi = 3.14159265358979d0


        I4f2d2d = -( 1.d0/(16.d0*pi*pi) )*(((x1 + x2)*dlog(x1/x2))/
     .  ((x2-x1)**3.d0) + 2.d0/((x2-x1)**2.d0))
        
        end function  
C     =========================================================
C     ==============================================================
C     Limit function for I3(m1,m2,m3) called I3f2d(x1,x2).
C     Notation x1,x2 is such that x1 represents one pair of degenerate
C     mass SQUARED where as x2 denotes the other pair of degenerate mass 
C     SQUARED 
C     ==============================================================

        double precision function I3f2d(x1,x2)
        implicit none
        double precision x1,x2,pi

        pi = 3.14159265358979d0

        I3f2d = -( 1.d0/(16.d0*pi*pi) )*((x1*dlog(x1/x2))/
     .  ((x1-x2)**2.d0) + 1.d0/(x2-x1))
        
        end function  
C     =========================================================

