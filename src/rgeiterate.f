****f* SuSeFLAV/rgeiterate.f 
* Program : SuSeFLAV: Supersymmetry seesaw and flavor violation calculator
* Version : 1.0
* Website : http://cts.iisc.ernet.in/Suseflav/main.html
* Authors : Debtosh Choudhury  debtosh@cts.iisc.ernet.in 
*           Raghuveer Garani   rgarani@cts.iisc.ernet.in
*           Sudhir Vempati     vempati@cts.iisc.ernet.in
*  NAME
*    RECURSIVE Subroutine rgeiterate
*  SYNOPSIS
*    Run mssm and sm  RGEs, computes tree level physical masses 
*    and the corresponding one loop correction 
*    for MSSM parameters. 
*
*  FUNCTION
*    This suborutine integrates MSSM RGEs and computes low 
*    energy spectrum at msusy.  
*
*  INPUTS
*     MW         - Mass of W boson
*     MZ         - Mass of Z boson 
*     MT         - top quark mass in \overline{DR} scheme
*     mb         - bottom quark  mass in \overline{DR} scheme 
*     mtau       - \tau lepton  mass in \overline{DR} scheme 
*     msusyold   - Initial guess value of msusy = sqrt(m0^2 + 4 m12^2)
*     vevsc      - scaled vev. vev/root2
*     vevin      - vev at MZ 
*     yuin       - (3x3) up type yukawa matrix 
*     ydin       - (3x3) down type yukawa matrix 
*     yein       - (3x3) matrix yukawa for leptons  
*     alphaDR    - alpha_{\overline{DR}}
*     alph1in    -  g1^2/(16 pi^2) at MZ
*     alph2in    -  g2^2/(16 pi^2) at MZ
*     alph3in    -  g3^2/(16 pi^2) at MZ
*     mur        -  \mu at mz 
*     bmur       -  \bmu at mz 
*     msusy      - susy breaking scale
*     prnstat    - Print control
*     itcount    - Rge iteration count
*
*  RESULT
*     MTc_mz     -  one loop correction to top quark 
*     mBc_mz     -  one loop correction to bottom  quark
*     mtauc_mz   -  one loop correction to tau lepton
*     delalphem  -  one loop correction to em coupling
*     delalphas  -  one loop correction to strong coupling
*     murge      -  RGE output: \mu at M_{susy} scale
*     bmurge     -  RGE output: b_\mu at M_{susy} scale
*     newtbeta   -  Ratio of vev at msusy from rge running.
*     msusynew   -  geometric mean of stop1 and stop2
*     flags      -  flags problem with rge running, if any.
*     itcount    - iteration count
*     stopratu 
*     stopratd   - variables used to check for global convergence of \mu
*     sinsqtheff - one loop corrected effective sinsqthw
*
*    Calculated tree level and 1-loop masses are stored in 
*    common blocks        
*
*  EXAMPLE
*
*      RECURSIVE SUBROUTINE rgeit(MW,MZ,MT,MTc_mz,mB,mBc_mz,mTau,
*     $     mtauc_mz,msusyold,vevsc,vevin,vev1in,vev2in,yuin,
*     $     ydin,yein,alphaDR,alph1in,alph2in,delalphem,
*     $     alph3in,delalphas,mur,bmur,murge,bmurge,prnstat,check,
*     $     newtbeta, MTatMZ,msusynew,mursq,try,newmh0sq,cheg,
*     $     flags,runum,itcount,stopratu,stopratd,sinsqtheff)
*
*  NOTES
*     Common Blocks used:
*      common/mascorr/ MT_qcd_corr
*      common/sminputs/ mbpole, mtaupole, Mtpole
*      common/loops/ lopt,rhn
*      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
*      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
*      common/charinputs/case, model
*      common/rgeinput_high/ MX, M1X,M2X,M3X,mh10,mh20,
*     $     mQ0,mU0,mD0,mE0,mL0,mNU0
*      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
*     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut
*      common/softout_mat_mz/mSQRGz,mSURGz,mSDRGz,AURGz,ADRGz,
*     $     AERGz,ANURGz,mSLRGz,mSERGz, mSNURGz,ONz,OCLz,OCRz,
*     $     MCharz, MNeutz
*      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
*     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG
*      common/rgeoutput_MZ/ mh1mzz,mh2mzz,M1tmz,M2tmz,M3tmz,
*     $     alph1MZ,alph2MZ,alph3MZ,vev1mz,vev2mz,yumz,ydmz,yemz
*      common/sparticles_susy/SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
*     $     mh0sq,mhu0sq,mhpmsq,mA0sq
*      common/sparticles_MZ/SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
*     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz
*      common/soft_mat_susy/ mSL1, SLeg, USL, SNeg, USN, SUeg, USU,
*     $     SDeg, USD
*      common/soft_mat_mz/ mSL1z, SLegz, USLz, SNegz, USNz, SUegz, USUz,
*     $     SDegz, USDz
*      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
*     $     thetau,thetad,thetae    
*      common/mu_mz/ murgemz
*      common/rgeopt_mz/mt_mz, mb_mz, mtau_mz
*      common/runningew/MWsq_mz,MZsq_mz
*      common/hzV/ delta1,delta2
*      common/mu_rge/mufrge
*      common/runningmass_susy/MT_susy, MB_susy, Mtau_susy
*      common/mAflag/flag_bmu
* 
*     External Routines used: 
*      EXTERNAL completerun, runtomz,iterate,rewsbcor, coratmz
*
*  BUGS
*    ---
*  SEE ALSO
*    -----
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!================================================================================
      RECURSIVE SUBROUTINE rgeit(MW,MZ,MT,MTc_mz,mB,mBc_mz,mTau,
     $     mtauc_mz,msusyold,vevsc,vevin,vev1in,vev2in,yuin,
     $     ydin,yein,alphaDR,alph1in,alph2in,delalphem,
     $     alph3in,delalphas,mur,bmur,murge,bmurge,prnstat,check,
     $     newtbeta, MTatMZ,msusynew,mursq,try,M3t,
     $     flags,runum,itcount,stopratu,stopratd,sinsqtheff,
     $     exitcalc)

      IMPLICIT NONE

      integer num, try,rhn, try1,i,j,k
      integer prnstat,lopt,AOK,check
      integer  runum, itcount,muflag

      character*1 exitcalc
      character*3 case
      character*4 model
      CHARACTER*10 flag_bmu
      CHARACTER*100 flags
      
      real sinsqtheffr
      double precision alphaDR
      double precision mur0,bmur0,bmur,mur,bmurge,murge,mAm3
      double precision M1X,M2X,M3X,MX,mh1mz,mh2mz,mh10,mh20,mursq
      double precision M1tz,M2tz,M3tz,m10,m20, sgnmu

      double precision MT_susy,MB_susy,mTau_susy

      double precision m0,m12,tanbeta,a0,ue3
      double precision mh0sq,mhu0sq,mhpmsq,mA0sq
      double precision mq0(3,3),mu0(3,3),MR1,MR2,MR3
      double precision md0(3,3),ml0(3,3),me0(3,3),mnu0(3,3) !,b1,b2,b3
      double precision msusy    !,alph10,alph20,alph30,tZ,tq0
      double precision ANURG(3,3)
      double precision yuRG(3,3),ydRG(3,3),yeRG(3,3),Ynui(3,3)

      double precision SUeg(6),SDeg(6),SLeg(6),SNeg(3),mSL1(6,6)
      double precision mSQRG(3,3),mSURG(3,3),mSDRG(3,3),mSLRG(3,3),
     $     mSERG(3,3),mSNURG(3,3)
      double precision ADRG(3,3),AURG(3,3),AERG(3,3)
      double precision SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      double precision USD(6,6),USL(6,6),USU(6,6),USN(3,3) !,mHcharsq, mHchar,Bbsg


      double precision mTatMz
      
      DOUBLE PRECISION alph1,alph2,alph3,newtbeta
      DOUBLE PRECISION Mtaupole,MTpole,Mbpole
      

      double precision mT, mB, mTau
      double precision vev1,vev2
      double precision vevin, vev1in, vev2in
      double precision vevsc
      double precision alph1in,alph2in,alph3in
      double precision yuin(3,3),ydin(3,3),yein(3,3)

      double precision MChar(2,2)
      double precision Neuevi(4),ONL(4,4),ONR(4,4),MNeut(4,4),
     $     Neg(4),ON(4,4)
      double precision OCR(2,2),OCL(2,2),Ceg(2)

C--------------------------------------------------------------------------
      double precision mtn,mbn,mtaun,mtno,mbno,mtauno
      DOUBLE PRECISION alph3no, alph2no,alph1no,mbdrbar
      DOUBLE PRECISION ratio,msnew,ratio1
      data mtno/ 1 * 0.d0/, mbno/ 1 * 0.d0/, mtauno/ 1 * 0.d0/
!-------------------------------------------
      DOUBLE PRECISION  yumz(3,3), ydmz(3,3), yemz(3,3)
      DOUBLE PRECISION alph3MZ, alph2MZ, alph1MZ
      DOUBLE PRECISION alphas1, alphaem
      DOUBLE PRECISION delta1,delta2
      DOUBLE PRECISION newvev
      DOUBLE PRECISION M3t

      DOUBLE PRECISION msusyold, msusynew 
      DOUBLE PRECISION MWc_mz,MZc_mz,MTc_mz,mBc_mz,mTauc_mz
      DOUBLE PRECISION AURGz(3,3),ADRGz(3,3),AERGz(3,3)
      DOUBLE PRECISION mSQRGz(3,3), mSURGz(3,3),mSDRGz(3,3)
      DOUBLE PRECISION mSLRGz(3,3),mSNURGz(3,3), mSERGz(3,3),bmur_conv
      DOUBLE PRECISION bmurgemz,M1tmz,M2tmz,M3tmz, mh1mzz,mh2mzz,mu_conv
      DOUBLE PRECISION murgemz, newtbetamz,vev1mz,vev2mz, ANURGz(3,3)

      double precision mh0sqz,mhu0sqz,mhpmsqz,mA0sqz


      double precision SUegz(6),SDegz(6),SLegz(6),SNegz(3)
      double precision SUeggz(6),SDeggz(6),SLeggz(6),SNeggz(3)
      double precision USDz(6,6),USLz(6,6),USUz(6,6),USNz(3,3) 
      double precision MCharz(2,2)
      double precision MNeutz(4,4),
     $     Negz(4),ONz(4,4)
      double precision OCRz(2,2),OCLz(2,2),Cegz(2)
      DOUBLE PRECISION mSL1z(6,6)
      DOUBLE PRECISION delalphas, delalphem,alph2n,alph3n,alph1n


      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae

      DOUBLE PRECISION newmA0sq,newmhpmsq,newmh0sq,newmHu0sq
      DOUBLE PRECISION STeg(2),SBeg(2)
      DOUBLE PRECISION SCeg(2),SUqeg(2),STaueg(2),SEeg(2),SMueg(2)
      DOUBLE PRECISION SSTeg(2),SDneg(2),tsnu,musnu,elsnu
      DOUBLE PRECISION stopratu(55),stopratd(55), scale,vr

      double precision ratiou, ratiod, sinsqtheff,spectol
      
      data STeg/ 2 * 0.d0/, SBeg/ 2 * 0.d0/, SCeg/ 2 *0.d0/
      data SUqeg/ 2 * 0.d0/, SSTeg/ 2 * 0.d0/, STaueg/ 2 *0.d0/
      data SEeg/ 2 * 0.d0/, SMueg/ 2 * 0.d0/, SDneg/ 2 *0.d0/

      data tsnu/ 1 * 0.d0/, musnu/ 1 * 0.d0/, elsnu/ 1 *0.d0/

      double precision MT_qcd_corr, mt_mz, mb_mz, mtau_mz,mufrge
      DOUBLE PRECISION Cheg(2),neuteg(4),mtaupoledrbar,mTauMZmsbar,
     $     mTauMZdrbar
      DOUBLE PRECISION MWsq_mz,mZsq_mz,e1,yukgut(126)
      double precision VCKM(3,3), MWpole, MZpole,MW,MZ,pi

      integer tacsup,tacsdn,tacslp,tacsnu,tachiggs,tact
      integer tacsupz,tacsdnz,tacslpz,tacsnuz,tachiggsz
!---------------------------------------------------------------------
!--define by priyanka    

      integer nok,nbad,i0
      double precision tq0, tZ,yy_sm(31),n0,x1,x2,eps,hmin,h1
      double precision yuMSbar(3,3),ydMSbar(3,3),yeMSbar(3,3),vevMSbar
      double precision yudiag(3),yddiag(3),UQuMS(3,3),UTUMS(3,3)
      double precision UQdTMS(3,3),UQudr(3,3),UTUdr(3,3),UQddr(3,3)
      double precision UTDdr(3,3),UQddrT(3,3),VCKMdrbar(3,3)
      double precision yudrdig(3),yddrdiag(3),g3MS,GF,p,q,g3
      double precision UQdMS(3,3),UTDMS(3,3),UQdMST(3,3),VCKMMSbar(3,3)
      double precision alph3MSbar,alph2MSbar,alph1MSbar,alpha2dr
      double precision gpMS,g2MS,alphemMS,alphemdr,sinthwMS,alpha3dr
      double precision delrMS,sincorrec,sinthwsqdr,vevdr,deldowndr
      double precision delupdr,delelecdr,Mtdr,Mcdr,Mudr,Mbdr,Msdr,Mddr
      double precision Mtaudr,Mmudr,Medr,vevscdr,beta,alpha1dr
      double precision alpha1MZ,alpha2MZ,alpha3MZ,mbmzdrbar,mbMZmsbar

      DOUBLE PRECISION alph3drbar,alph2drbar,alph1drbar,alpha1MS
      DOUBLE PRECISION alpha2MS,alpha3MS,delupdrdown,deldowndrdown
      DOUBLE PRECISION delelecdrdown,MtMS,McMS,MuMS,MbMS,MsMS,MdMS
      DOUBLE PRECISION MtauMS,MmuMS,MeMS,vevMS,alpha1susy
      DOUBLE PRECISION alpha2susy,alpha3susy,sinsqthw_susy,vevscMS

      DOUBLE PRECISION gpdr,gdr,sinsqthwdr_susy,b1vecint(18),tautil
      DOUBLE PRECISION b3vecint(18),ttil,ctil,util,btil,stil,dtil,Mhpm
      DOUBLE PRECISION mutil,etil,taunutil,munutil,enutil,deltasusy3int
      DOUBLE PRECISION Massvecint(18),deltasusy1int,deltasusy2int,mA
      DOUBLE PRECISION correction_susy,MTc_msusy,mtsusy,mBc_msusy
      DOUBLE PRECISION mbcor_susy,mbsusy,mTauc_msusy,mtaucor_susy
      DOUBLE PRECISION mtaususy,pizzT_susy,tqsusy,muvalue,b2vecint(18)
      DOUBLE PRECISION MZ_susyup,MZ_susyd,pizzT_susyd
      DOUBLE PRECISION delalph1,delalph2,delalph3,g,gp
      DOUBLE PRECISION sinsqthw_Mz,mbcor_Mz,correction_Mz,mtaucor_Mz
      DOUBLE PRECISION pizzTSM
!---------------------------------------------
      common/qcd_cor/mbmzdrbar,mbMZmsbar
      common/yukawa_MZ/yuMZ,ydMZ,yeMZ
      common/alpha_MZ/alpha1MZ,alpha2MZ,alpha3MZ
!---------------------------------------------------      

      common/VCKMparam/ VCKM
      common/mascorr/ MT_qcd_corr
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole
      common/mtaurmass/ mtaupoledrbar,mTauMZmsbar,mTauMZdrbar
      common/loops/ lopt,rhn
      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
      common/charinputs/case, model
      common/rgeinput_high/ MX, M1X,M2X,M3X,mh10,mh20,
     $     mQ0,mU0,mD0,mE0,mL0,mNU0

      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut

      common/softout_mat_mz/mSQRGz,mSURGz,mSDRGz,AURGz,ADRGz,
     $     AERGz,ANURGz,mSLRGz,mSERGz, mSNURGz,ONz,OCLz,OCRz,
     $     MCharz, MNeutz

      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG

!      common/rgeoutput_MZ/ mh1mzz,mh2mzz,M1tmz,M2tmz,M3tmz,
!     $     alph1MZ,alph2MZ,alph3MZ,vev1mz,vev2mz,yumz,ydmz,yemz

      common/sparticles_susy/SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq

      common/sparticles_MZ/SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz

      common/soft_mat_susy/ mSL1, SLeg, USL, SNeg, USN, SUeg, USU,
     $     SDeg, USD

      common/soft_mat_mz/ mSL1z, SLegz, USLz, SNegz, USNz, SUegz, USUz,
     $     SDegz, USDz

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae
      
      common/ma/mAm3

      common/mu_mz/ murgemz
      common/rgeopt_mz/mt_mz, mb_mz, mtau_mz
      common/runningew/MWsq_mz,MZsq_mz

      common/hzV/ delta1,delta2
      common/mu_rge/mufrge
      common/runningmass_susy/MT_susy, MB_susy, Mtau_susy
      common/mAflag/flag_bmu
      common/unif/ e1,yukgut
      common/ptolerance/spectol

      common/tac_susy/ tacsup,tacsdn,tacslp,tacsnu,tachiggs
      common/tac_mz/ tacsupz,tacsdnz,tacslpz,tacsnuz,tachiggsz



!-------------------------
      EXTERNAL completerun, runtomz,iterate,rewsbcor, coratmz
      external RK4ROUTINE,QMSRK4,smrge,smrgemt,SVD,dag,matmulti
      external topcor,bottomcor,taucor
!----------------------------------

      include 'stdinputs.h'
!------------------


      M3t = 0.d0

      tacsup = 0
      tacsdn = 0
      tacslp = 0
      tacsnu = 0
      tachiggs = 0

      tacsupz = 0
      tacsdnz = 0
      tacslpz = 0
      tacsnuz = 0
      tachiggsz = 0

      MW = MWpole
      MZ = MZpole
      pi = 4.d0 * datan(1.d0)
      GF = 1.16639*10.d0**(-5)
      beta = datan(tanbeta)
      
      if(exitcalc.eq.'F') flags=' AOK'
!!!!!!!!!SM threshold correction of alpha's by priyanka
      If(itcount.eq.1)then

      delalph1=0.0
      delalph2=0.0
      delalph3=0.0
!      print*,"test output"
      else 
    
      delalph1= (1.0/(2.0*pi))*((17.d0/30.d0)*log(Mtpole/Mz)
     $          +(1.d0/10.d0)*log(125.35/Mz))
      delalph2= (1.0/(2.0*pi))*(log(Mtpole/Mz)
     $          +(1.d0/6.d0)*log(125.35/Mz))
      delalph3= (1.0/(2.0*pi))*((2.d0/3.d0)*log(Mtpole/Mz))

      endif
      alph1no =(1.0/(4.0*pi))*(alph1in/(1.0+alph1in*delalph1))
      alph2no =(1.0/(4.0*pi))*(alph2in/(1.0+alph2in*delalph2))
      alph3no =(1.0/(4.0*pi))*(alph3in/(1.0+alph3in*delalph3))
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      print*,"rgeiterate apha1,2,3", alph1no,alph2no,alph3no
      print*,"in rgeieterate.f",yuin(3,3),ydin(3,3),yein(3,3), 
     $        alph1no,alph2no,alph3no
!!--------------------------------------------------------------
!---------------------------------------------------------------
!modify by priyanka in june 2021

      
C     ===================================
C     Running of the Yukawas of the SM !!
C     ===================================
      MX = 5.d0*(10.d0**19.d0)
      tq0   =  dLog(MX**2.d0/msusyold**2.d0)
      tZ    =  dLog(MX**2.d0/MZ**2.d0)
C-----------------------------
C	top Yukawa
C-----------------------------

      i0 = 0

      do i = 1, 3
      yy_sm(i0 + i)   = yuin(1,i)
      yy_sm(i0+3 + i) = yuin(2,i)
      yy_sm(i0+6 + i) = yuin(3,i)
      enddo

C-------------------------------
C	bottom Yukawa
C-------------------------------

      i0 = 9

      do i = 1, 3
      yy_sm(i0 + i)   = ydin(1,i)
      yy_sm(i0+3 + i) = ydin(2,i)
      yy_sm(i0+6 + i) = ydin(3,i)
      enddo 


C-----------------------------
C     tau Yukawa
C-----------------------------

      i0 = 18

      do i = 1, 3
      yy_sm(i0 + i)   = yein(1,i)
      yy_sm(i0+3 + i) = yein(2,i)
      yy_sm(i0+6 + i) = yein(3,i)
      enddo 

      yy_sm(28) = alph3no
      yy_sm(29) = alph2no
      yy_sm(30) = alph1no

      yy_sm(31) = vevin
      print*,"before SM rge runing"
!--------------------------------------------------------------------------
!     RUNNING UP THE SM YUKAWAS AND GAUGE COUPLINGS
!--------------------------------------------------------------------------
      n0 = 31
!      check=0
      x2 = tq0
      x1 = tz
      h1   =  -1.d-5
      hmin =  2.d-8
      eps  =  1.d-6
      print*,"before",yy_sm(31),yy_sm(9)
!      call RK4ROUTINE(yy_sm,n0,x1,x2,eps,h1,hmin,nok,nbad,smrge,
!     $        QMSRK4,check)
         if(check.eq.100)then
            flags = 'variable underflow '
            return
         endif
      print*,"after",n0,x2,x1,tz,tq0,MZ,MX,msusyold,h1,hmin,eps,check,
     $       yy_sm(31),yy_sm(9)
      i0 = 0

      do i = 1, 3
      yuMSbar(1,i) = yy_sm(i0 + i)*(4.d0 * pi) 
      yuMSbar(2,i) = yy_sm(i0+3 + i)*(4.d0 * pi) 
      yuMSbar(3,i) = yy_sm(i0+6 + i)*(4.d0 * pi) 
      enddo

C-------------------------------
C	bottom Yukawa
C-------------------------------

      i0 = 9

      do i = 1, 3
      ydMSbar(1,i) = yy_sm(i0 + i)*(4.d0 * pi)    
      ydMSbar(2,i) = yy_sm(i0+3 + i)*(4.d0 * pi)  
      ydMSbar(3,i) = yy_sm(i0+6 + i)*(4.d0 * pi)  
      enddo 


C-----------------------------
C     tau Yukawa
C-----------------------------

      i0 = 18

      do i = 1, 3
      yeMSbar(1,i) = yy_sm(i0 + i)*(4.d0 * pi)    
      yeMSbar(2,i) = yy_sm(i0+3 + i)*(4.d0 * pi)  
      yeMSbar(3,i) = yy_sm(i0+6 + i)*(4.d0 * pi)  
      enddo 

      alph3MSbar = yy_sm(28)*(4.d0 * pi)  
      alph2MSbar = yy_sm(29)*(4.d0 * pi)  
      alph1MSbar = yy_sm(30)*(4.d0 * pi)  

      vevMSbar = yy_sm(31) 
      
      
      call SVD(3, 3, yuMSbar,3, yudiag, UQuMS,3, UTUMS,3, 0)
      call SVD(3, 3, ydMSbar,3, yddiag, UQdMS,3, UTDMS,3, 0)
      

      call dag(UQdMS,UQdTMS)
      call matmult(UQuMS,UQdTMS,VCKMMSbar)
      
!------------------------------------------------------------------------------
!-----DR conversion of gauge couplings , sinsqthw and vev
      If(itcount.eq.1)then
      
      alpha1dr=alph1MSbar
      alpha2dr= alph2MSbar/(1.d0-alph2MSbar/(6.d0*pi))
      alpha3dr= alph3MSbar/(1.d0-alph3MSbar/(4.d0*pi))       
      
      gpdr=Sqrt(3.d0/5.d0)*Sqrt(4.d0*pi*alpha1dr)
      gdr=Sqrt(4.d0*pi*alpha2dr)
!      g1MS=Sqrt(4.d0*pi*alph1MSbar)
!      g2MS=Sqrt(4.d0*pi*alph2MSbar)
      sinsqthwdr_susy= gpdr**2.d0/((gpdr**2.d0+gdr**2.d0))
      alphemdr=sinsqthwdr_susy*alpha2dr
!      alphemMS= ((g1MS*g2MS)**2.d0)/(4.d0*pi*(g1MS**2.d0+g2MS**2.d0))
!      alphemdr=alphemMS/(1.d0-alphemMS/(6.d0*pi))
!      sinthwMS=g1MS/sqrt(g1MS**2.d0+g2MS**2.d0)
!      delrMS=1.d0-pi*alphemMS/(sqrt(2.d0)*GF*MZ*MZ*sinthwMS**2.d0*
!     $ (1.d0-sinthwMS**2.d0))
!      sincorrec=pi*alphemdr/(sqrt(2.d0)*MZ*MZ*GF*(1.d0-delrMS))
!      sinthwsqdr=0.5 -sqrt(0.25 - sincorrec)
!      vevdr= sqrt(MZ*MZ*sinthwsqdr*(1.d0-sinthwsqdr)/(pi*alphemdr))
       vevdr= vevMSbar                                                 !(*Note: we are neglecting the difference between MS bar and DR for vev*)
!      MZ_susyup=(1.0/2.0)*dsqrt((gpdr**2.d0+gdr**2.d0) *
!     $           (vevMSbar**2.d0))                                          !(*We need to check...MZ is defined as a mixer of MS bar and DR bar*)
!      vevdr=sqrt(4.d0*(MZ_susyup**2.d0)/
!     $      (gpdr**2.d0+gdr**2.d0))
!       vevdr=sqrt(4.d0*MZ**2.d0/(gpdr**2.d0+gdr**2.d0))
      print*,"vevdr",vevdr    
!-----------DR conversion of fermion masses-------------------------------------
      delupdr = 1.d0-alpha3dr/(3.d0*pi)-4.d0*alphemdr/(36.d0*pi)-
     $          43.d0*alpha3dr**2.d0/(144.d0*pi*pi)      
      deldowndr = 1.d0-alpha3dr/(3.d0*pi)-1.d0*alphemdr/(36.d0*pi)-
     $          43.d0*alpha3dr**2.d0/(144.d0*pi*pi)  
      delelecdr = 1.d0-alphemdr/(4.d0*pi)
      
      Mtdr = (yuMSbar(3,3)*vevMSbar/Sqrt(2.d0))*delupdr
      Mcdr = (yuMSbar(2,2)*vevMSbar/Sqrt(2.d0))*delupdr
      Mudr = (yuMSbar(1,1)*vevMSbar/Sqrt(2.d0))*delupdr
      
      Mbdr = (yddiag(3)*vevMSbar/Sqrt(2.d0))*deldowndr
      Msdr = (yddiag(2)*vevMSbar/Sqrt(2.d0))*deldowndr
      Mddr = (yddiag(1)*vevMSbar/Sqrt(2.d0))*deldowndr
      
      Mtaudr = (yeMSbar(3,3)*vevMSbar/Sqrt(2.d0))*delelecdr
      Mmudr = (yeMSbar(2,2)*vevMSbar/Sqrt(2.d0))*delelecdr
      Medr = (yeMSbar(1,1)*vevMSbar/Sqrt(2.d0))*delelecdr
!-------------------------------------------------------------------------     
      else 
      
      b1vecint = (/17.0/60.0,17.0/60.0,17.0/60.0,1.0/12.0,
     $            1.0/12.0,1.0/12.0,1.0/4.0,1.0/4.0,1.0/4.0,1.0/20.0,
     $            1.0/20.0,1.0/20.0,0.0,0.0,0.0,2.0/5.0,1.0/20.0,
     $            1.0/10.0/)
      b2vecint = (/1.0/2.0,1.0/2.0,1.0/2.0,1.0/2.0,
     $            1.0/2.0,1.0/2.0,1.0/6.0,1.0/6.0,1.0/6.0,1.0/6.0,
     $            1.0/6.0,1.0/6.0,0.0,4.0/3.0,0.0,2.0/3.0,1.0/6.0,
     $            1.0/3.0/)
      b3vecint = (/1.0/3.0,1.0/3.0,1.0/3.0,1.0/3.0,
     $            1.0/3.0,1.0/3.0,0.0,0.0,0.0,0.0,
     $            0.0,0.0,0.0,0.0,2.0,0.0,0.0,0.0/)         

      ttil = Sqrt((SUegg(1)+SUegg(2))/2.d0)
      ctil = Sqrt((SUegg(3)+SUegg(4))/2.d0)
      util = Sqrt((SUegg(5)+SUegg(6))/2.d0)
      btil = Sqrt((SDegg(1)+SDegg(2))/2.d0)
      stil = Sqrt((SDegg(3)+SDegg(4))/2.d0)
      dtil = Sqrt((SDegg(5)+SDegg(6))/2.d0)
      tautil = Sqrt((SLegg(1)+SLegg(2))/2.d0)
      mutil = Sqrt((SLegg(3)+SLegg(4))/2.d0)
      etil = Sqrt((SLegg(5)+SLegg(6))/2.d0)
      taunutil = Sqrt(SNegg(1))
      munutil = Sqrt(SNegg(2))
      enutil = Sqrt(SNegg(3))
      mA = sqrt(mA0sq)
      Mhpm = sqrt(mhpmsq)
      muvalue=abs(murge)
    
      Massvecint = (/ttil,ctil,util,btil,stil,dtil,tautil,mutil,etil,
     $               taunutil,munutil,enutil,M1tz,M2tz,M3t,muvalue,
     $               mA,mhpm/)

      deltasusy1int = 1.0d0/(2.0d0*Pi)*DOT_product(b1vecint,
     $                Log(Massvecint/msusyold))
      deltasusy2int = 1.0d0/(2.0d0*Pi)*DOT_product(b2vecint,
     $                Log(Massvecint/msusyold))
      deltasusy3int = 1.0d0/(2.0d0*Pi)*DOT_product(b3vecint,
     $                Log(Massvecint/msusyold))


      alpha1susy=(alph1MSbar/(1.0d0+alph1MSbar*deltasusy1int))
      alpha2susy=(alph2MSbar/(1.0d0+alph2MSbar*deltasusy2int))
      alpha3susy=(alph3MSbar/(1.0d0+alph3MSbar*deltasusy3int))

      alpha1dr=alpha1susy
      alpha2dr= alpha2susy/(1.d0-alpha2susy/(6.d0*pi))
      alpha3dr= alpha3susy/(1.d0-alpha3susy/(4.d0*pi))       
      
  
!      g3 = Sqrt(4.d0*pi*alpha3dr)
      gpMS=Sqrt(3.d0/5.d0)*Sqrt(4.d0*pi*alpha1susy)
      g2MS=Sqrt(4.d0*pi*alpha2susy)
      g3MS=Sqrt(4.d0*pi*alph3MSbar)
      sinsqthw_susy= gpMS**2.d0/((gpMS**2.d0+g2MS**2.d0))
      
      mt_susy=(yuMSbar(3,3)*vevMSbar/Sqrt(2.d0))
      mb_susy=(ydMSbar(3,3)*vevMSbar/Sqrt(2.d0))
      mtau_susy = (yeMSbar(3,3)*vevMSbar/Sqrt(2.d0))
      
!----------------------------------------------------------------------------
!       self energy corrections of fermion   
!---------------------------------------------------------------------------      
      p = mtpole

      q = Msusyold

!      call topcor(p,q,g,gp,g3,M3tmz,mt_mz,mb_mz,tanbeta,
!     $     yuMZ,ydMZ,SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
!     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz,ONz,OCLz,
!     $     OCRz,sinsqthw_mz,correction)

      call topcor(p,q,g2MS,gpMS,g3MS,M3t,mt_susy,mb_susy,tanbeta,
     $     yuMSbar,ydMSbar,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,ON,OCL,
     $     OCR,sinsqthw_susy,correction_susy)
!----------------------------------------------------------------------------


      p = mbmzMSbar

      q = Msusyold

      call bottomcor(p,q,g2MS,gpMS,g3MS,M3t,mt_susy,mb_susy,tanbeta,
     $     yuMSbar,ydMSbar,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,ON,OCL,
     $     OCR,sinsqthw_susy,mbcor_susy,mbdrbar)
      

!----------------------------------------------------------------------------

      p = mtaupole

      q = Msusyold

      call taucor(p,q,g2MS,gpMS,M3t,mtau_susy,tanbeta,yeMSbar,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,
     $     mhu0sq,mhpmsq,mA0sq,ON,OCL,
     $     OCR,sinsqthw_susy,mtaucor_susy)


      MTc_msusy =  correction_susy

      mBc_msusy =  mbcor_susy

      mTauc_msusy = mtaucor_susy


      mtsusy = mt_susy* (1.d0 + MTc_msusy)

      mbsusy = mb_susy / (1.d0 + mBc_msusy)

      mtaususy = mtau_susy * (1.d0 + mtauc_msusy)
      
!-----------DR conversion of fermion masses-------------------------------------
      gpdr=Sqrt(3.d0/5.d0)*Sqrt(4.d0*pi*alpha1dr)
      gdr=Sqrt(4.d0*pi*alpha2dr)
      sinsqthwdr_susy= gpdr**2.d0/((gpdr**2.d0+gdr**2.d0))
      alphemdr=sinsqthwdr_susy*alpha2dr
      
      delupdr = 1.d0-alpha3dr/(3.d0*pi)-4.d0*alphemdr/(36.d0*pi)-
     $          43.d0*alpha3dr**2.d0/(144.d0*pi*pi)      
      deldowndr = 1.d0-alpha3dr/(3.d0*pi)-1.d0*alphemdr/(36.d0*pi)-
     $          43.d0*alpha3dr**2.d0/(144.d0*pi*pi)  
      delelecdr = 1.d0-alphemdr/(4.d0*pi)
      
      Mtdr = mtsusy*delupdr
      Mcdr = (yuMSbar(2,2)*vevMSbar/Sqrt(2.d0))*delupdr
      Mudr = (yuMSbar(1,1)*vevMSbar/Sqrt(2.d0))*delupdr
      
      Mbdr = mbsusy*deldowndr
      Msdr = (yddiag(2)*vevMSbar/Sqrt(2.d0))*deldowndr
      Mddr = (yddiag(1)*vevMSbar/Sqrt(2.d0))*deldowndr
      
      Mtaudr = mtaususy*delelecdr
      Mmudr = (yeMSbar(2,2)*vevMSbar/Sqrt(2.d0))*delelecdr
      Medr = (yeMSbar(1,1)*vevMSbar/Sqrt(2.d0))*delelecdr
!-----------------------------------------------------------------------------------
!    vev SUSY correction
!-----------------------------------------------------------------------------------           
      p = MZpole
 
      q = msusyold

      call pizz(p,q,gdr,mtdr,mbdr,mtaudr,newtbeta,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     ON,OCL,OCR,sinsqthwdr_susy,pizzT_susy) 
      
!      call pizz(p,q,g2MS,mtsusy,mbsusy,mtaususy,newtbeta,SUegg,SDegg,
!     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
!     $     ON,OCL,OCR,sinsqthw_susy,pizzT_susy)      
!------here I took MZ pole mass but I think it should mass at $M_{susy}$ scale, which using vev value 
!        
      MZ_susyup=(1.0/2.0)*dsqrt((gpdr**2.d0+gdr**2.d0) *
     $           (vevMSbar**2.d0))
      vevdr=sqrt(4.d0*(MZ_susyup**2.d0+pizzT_susy)/
     $      (gpdr**2.d0+gdr**2.d0))                              !(*Note: we are neglecting the difference between MS bar and DR for vev*)
           
      endif
!------define yukawa and gauge couplings at msusy---------------

      vevscdr=vevdr/sqrt(2.d0)
 
      ydin(1,1) = (VCKMMSbar(1,1)*mDdr)/(4.d0*pi*vevscdr)
      ydin(1,2) = (VCKMMSbar(1,2)*mSdr)/(4.d0*pi*vevscdr)
      ydin(1,3) = (VCKMMSbar(1,3)*mBdr)/(4.d0*pi*vevscdr)
      ydin(2,1) = (VCKMMSbar(2,1)*mDdr)/(4.d0*pi*vevscdr)
      ydin(2,2) = (VCKMMSbar(2,2)*mSdr)/(4.d0*pi*vevscdr)
      ydin(2,3) = (VCKMMSbar(2,3)*mBdr)/(4.d0*pi*vevscdr)
      ydin(3,1) = (VCKMMSbar(3,1)*mDdr)/(4.d0*pi*vevscdr)
      ydin(3,2) = (VCKMMSbar(3,2)*mSdr)/(4.d0*pi*vevscdr)
      ydin(3,3) = (VCKMMSbar(3,3)*mBdr)/(4.d0*pi*vevscdr) 


      yein(1,1) = mEdr/(4.d0*pi*vevscdr)
      yein(1,2) = 0.d0
      yein(1,3) = 0.d0
      yein(2,1) = 0.d0
      yein(2,2) = mMUdr/(4.d0*pi*vevscdr)
      yein(2,3) = 0.d0
      yein(3,1) = 0.d0
      yein(3,2) = 0.d0
      yein(3,3) = mTaudr/(4.d0*pi*vevscdr)

      yuin(1,1) = mUdr/(4.d0*pi*vevscdr)
      yuin(1,2) = 0.d0
      yuin(1,3) = 0.d0
      yuin(2,1) = 0.d0
      yuin(2,2) = mCdr/(4.d0*pi*vevscdr)
      yuin(2,3) = 0.d0
      yuin(3,1) = 0.d0
      yuin(3,2) = 0.d0
      yuin(3,3) = MTdr/(4.d0*pi*vevscdr)      
      
      vevin = vevdr
      alph1n = alpha1dr/(4.d0 * pi)
      alph2n = alpha2dr/(4.d0 * pi)
      alph3n = alpha3dr/(4.d0 * pi)
      
 !---------------------end priyanka correction-------------------------------------------
      
!===============================================================================
!	                 STANDARD MODEL RUNNING ENDS
!===============================================================================
!--------------------------------------------------------------
!      alph1n = (alphaDR/(4.d0 * pi * (1.d0 - delalphem) * 
!     $     (1.d0 - sinsqtheff))) * (5.d0/3.d0) 
!      alph2n = alphaDR/(4.d0 * pi * (1.d0 - delalphem) * sinsqtheff)
!      alph3n = alph3in/(1.d0 - delalphas) 

!            write(*,*)"sinsqtheff",sinsqtheff

      if(itcount.eq.1)then
         mtn = MT
         mbn = mb
         mtaun = mtau
      else

         mtn = MT * (1.d0 + MTc_mz)

         mbn = mb/(1 + mbc_mz)

         mtaun =  mtau * (1.d0 + mtauc_mz)

!         mtaun =  mtauMZdrbar * (1.d0 + mtauc_mz)

      endif

!--------------------
      try = itcount

      if(itcount.eq.1)then

         mur0  = 0.d0
         bmur0 = 0.d0
         mur   = mur0              
         bmur  = bmur0             

      else

         msusyold = msusynew
         
      endif

!      call completerun(msusyold,vevin, yuin,ydin,yein,alph1in,
!     $     alph2in,alph3in,mur,bmur,murge,bmurge,prnstat,check,
!     $     newtbeta, msusynew, mursq,try,
!     $     flags,runum,itcount)
!#      do i=1,3
!#        do j=1,3
      print*,"complete run input msusy", yuin(3,3),ydin(3,3),yein(3,3),
     $        alph1n,alph2n,alph3n,msusyold,vevin,mur,bmur,itcount
!#      enddo
!#      enddo
      print*,"before calling spectral.l"
      call completerun(msusyold,vevin, yuin,ydin,yein,alph1n,
     $     alph2n,alph3n,mur,bmur,murge,bmurge,prnstat,check,
     $     newtbeta, msusynew, mursq,try,
     $     flags,runum,itcount)
      
      if(flags.eq.'variable underflow '.or.flags.eq.'NPERTYUK')then
         return
      endif

!--------------------------------



      num = 1
      if((Suegg(1).lt.0.d0).OR.(Suegg(2).lt.0.d0).OR.
     $     (Suegg(3).lt.0.d0).OR.(Suegg(4).lt.0.d0).OR.
     $     (Suegg(5).lt.0.d0).OR.(Suegg(6).lt.0.d0))then

         Suegg(1) = dabs(Suegg(1))
         Suegg(2) = dabs(Suegg(2))
         Suegg(3) = dabs(Suegg(3))
         Suegg(4) = dabs(Suegg(4))
         Suegg(5) = dabs(Suegg(5))
         Suegg(6) = dabs(Suegg(6))

         tacsup = 1


      else

         continue
      endif
!-----
!     Down sector

      if((SDegg(1).lt.0.d0).OR.(SDegg(2).lt.0.d0).OR.
     $     (SDegg(3).lt.0.d0).OR.(SDegg(4).lt.0.d0).OR.
     $     (SDegg(5).lt.0.d0).OR.(SDegg(6).lt.0.d0))then

         SDegg(1) = dabs(SDegg(1))
         SDegg(2) = dabs(SDegg(2))
         SDegg(3) = dabs(SDegg(3))
         SDegg(4) = dabs(SDegg(4))
         SDegg(5) = dabs(SDegg(5))
         SDegg(6) = dabs(SDegg(6))

         tacsdn = 1

      else

         continue
      endif

!-------
!     lepton sector

      if((SLegg(1).lt.0.d0).OR.(SLegg(2).lt.0.d0).OR.
     $     (SLegg(3).lt.0.d0).OR.(SLegg(4).lt.0.d0).OR.
     $     (SLegg(5).lt.0.d0).OR.(SLegg(6).lt.0.d0))then

         SLegg(1) = dabs(SLegg(1))
         SLegg(2) = dabs(SLegg(2))
         SLegg(3) = dabs(SLegg(3))
         SLegg(4) = dabs(SLegg(4))
         SLegg(5) = dabs(SLegg(5))
         SLegg(6) = dabs(SLegg(6))

         tacslp = 1

      else

         continue
      endif
!-------------------------
!     sneutrinos

      if((SNegg(1).lt.0.d0).OR.(SNegg(2).lt.0.d0).OR.
     $     (SNegg(3).lt.0.d0))then

         SNegg(1) = dabs(SNegg(1))
         SNegg(2) = dabs(SNegg(2))
         SNegg(3) = dabs(SNegg(3))

         tacsnu = 1

      else

         continue
      endif

!----------------------------
C     higgs FLAG

      if(mh0sq.lt.0.d0.or.mhu0sq.lt.0.d0.or.mhpmsq.lt.0.d0.or.
     $     mA0sq.lt.0.d0)then

         mh0sq = dabs(mh0sq)
         mhu0sq = dabs(mhu0sq)
         mhpmsq = dabs(mhpmsq)
         mA0sq  = dabs(mA0sq)
         
         tachiggs = 1
         
      else
         continue
      ENDIF

!-----------------------------

      iteratemu: if(itcount.eq.1)then

         mu_conv = murge
         bmur_conv = bmurge

      else iteratemu
         
         try1  = 0

         call iterate(msusyold,mtn,murge,bmur,newtbeta,
     $        msusyold,muflag,try1,mur,exitcalc,flags,mursq)

         mu_conv = mur
         bmur_conv = bmur
         
      endif iteratemu

      if(itcount.gt.15)then
         mu_conv =  dsqrt(murge*mur)
      endif

!      print*,"erwrerr",mu_conv,bmur_conv

      if(isnan(bmur_conv).or.isnan(mu_conv))then 

        do i = 1, 6
            SUegg(i) = 0.d0
            SDegg(i) = 0.d0
            SLegg(i) = 0.d0
         enddo

         SNegg(1) = 0.d0
         SNegg(2) = 0.d0
         SNegg(3) = 0.d0

         mA0sq = 0.d0
         mh0sq = 0.d0
         mHu0sq = 0.d0
         mHpmsq = 0.d0    

         Neg(1) = 0.d0
         Neg(2) = 0.d0
         Neg(3) = 0.d0
         Neg(4) = 0.d0
         
         Ceg(1) = 0.d0
         Ceg(2) = 0.d0

         mur = 0.d0

         flags = "MUNOC"
!         print*,"NAN"
         return
      endif

!      if(isnan(bmur_conv)) return
!-------------

      exitp: if(exitcalc.eq.'T')then

         tact = 0

         tact = tacsup + tacsdn + tacslp + tacsnu + tachiggs
         
!         print*,"tact = ", tact
         
         tac:if(tact.eq.0)then

!            print*,"exitcalc is True"
            
c$$$  print*,"vev1, vev2 = ", vev1, vev2
c$$$  print*,"vev1in, vev2in = ", vev1in, vev2in

            MW = dsqrt((alph2*4.d0*pi*pi)*(vev1**2.d0 + vev2**2.d0))
            
            MZ = dsqrt(((alph2*4.d0*pi*pi) + 
     $           (alph1*4.d0*pi*pi) * (3.d0/5.d0)) *
     $           (vev1**2.d0 + vev2**2.d0))
            
c$$$            print*,"convergered mu at last run = ", mu_conv
c$$$            print*,"convergered bmu at last run = ", bmur_conv
c$$$            print*,"treelevel mu at last run = ", murge
c$$$            print*,"treelevel bmu at last run = ", bmurge

            call softspectrum(msusyold,sgnmu,newtbeta,mSQRG,mSDRG,mSURG,
     $           AURG,ADRG,vev1,vev2,MW*MW,MZ*MZ,mSLRG,mSERG,AERG,yuRG,
     $           yeRG,ydRG,M1tz,M2tz,mu_conv,bmur_conv,SUegg,USU,SUeg,
     $           SDegg,USD,SDeg,SLegg,USL,SLeg,SNegg,USN,SNeg,MNeut,ON,
     $           Neg,MChar,mSL1,OCR,OCL,Ceg,AOK,MT_susy,MB_susy,
     $           Mtau_susy,mh0sq,mhu0sq,mhpmsq,mA0sq,Neuevi,ONL,ONR)


!---------
!     Up sector

            if((Suegg(1).lt.0.d0).OR.(Suegg(2).lt.0.d0).OR.
     $           (Suegg(3).lt.0.d0).OR.(Suegg(4).lt.0.d0).OR.
     $           (Suegg(5).lt.0.d0).OR.(Suegg(6).lt.0.d0))then

               Suegg(1) = dabs(Suegg(1))
               Suegg(2) = dabs(Suegg(2))
               Suegg(3) = dabs(Suegg(3))
               Suegg(4) = dabs(Suegg(4))
               Suegg(5) = dabs(Suegg(5))
               Suegg(6) = dabs(Suegg(6))

               tacsup = 1


            else

               continue
            endif
!-----
!     Down sector

            if((SDegg(1).lt.0.d0).OR.(SDegg(2).lt.0.d0).OR.
     $           (SDegg(3).lt.0.d0).OR.(SDegg(4).lt.0.d0).OR.
     $           (SDegg(5).lt.0.d0).OR.(SDegg(6).lt.0.d0))then

               SDegg(1) = dabs(SDegg(1))
               SDegg(2) = dabs(SDegg(2))
               SDegg(3) = dabs(SDegg(3))
               SDegg(4) = dabs(SDegg(4))
               SDegg(5) = dabs(SDegg(5))
               SDegg(6) = dabs(SDegg(6))

               tacsdn = 1

            else

               continue
            endif

!-------
!     lepton sector

            if((SLegg(1).lt.0.d0).OR.(SLegg(2).lt.0.d0).OR.
     $           (SLegg(3).lt.0.d0).OR.(SLegg(4).lt.0.d0).OR.
     $           (SLegg(5).lt.0.d0).OR.(SLegg(6).lt.0.d0))then

               SLegg(1) = dabs(SLegg(1))
               SLegg(2) = dabs(SLegg(2))
               SLegg(3) = dabs(SLegg(3))
               SLegg(4) = dabs(SLegg(4))
               SLegg(5) = dabs(SLegg(5))
               SLegg(6) = dabs(SLegg(6))

               tacslp = 1

            else

               continue
            endif
!-------------------------
!     sneutrinos

            if((SNegg(1).lt.0.d0).OR.(SNegg(2).lt.0.d0).OR.
     $           (SNegg(3).lt.0.d0))then

               SNegg(1) = dabs(SNegg(1))
               SNegg(2) = dabs(SNegg(2))
               SNegg(3) = dabs(SNegg(3))

               tacsnu = 1

            else

               continue
            endif

!----------------------------
C     higgs FLAG

            if(mh0sq.lt.0.d0.or.mhu0sq.lt.0.d0.or.mhpmsq.lt.0.d0.or.
     $           mA0sq.lt.0.d0)then

               mh0sq = dabs(mh0sq)
               mhu0sq = dabs(mhu0sq)
               mhpmsq = dabs(mhpmsq)
               mA0sq  = dabs(mA0sq)
               
               tachiggs = 1
               
            else
               continue
            ENDIF


            tact = 0
            
            tact = tacsup + tacsdn + tacslp + tacsnu + tachiggs
            
            tac2: if(tact.eq.0.and.flag_bmu.eq.' AOK')then
               
               
               call REWSBCOR(sgnmu,mu_conv,bmur_conv,newtbeta,MTn,
     $              msusyold,msnew,M3t,mAm3,tanbeta,STeg,SCeg,SUqeg,
     $              SBeg,SSTeg,SDneg,STaueg,SMUeg,SEeg,tsnu,musnu,elsnu, 
     $              newmA0sq,newmh0sq,newmhpmsq,newmHu0sq,Cheg,neuteg) 
               
            else tac2
               flags = flag_bmu
!               print*,"mstop1 = ", SUegg(1), flags
               return
            endif tac2
            
         else tac

            flags = "TACSPEC"

!            print*,"Tachyonic Spectrum at MSUSY"


            do i = 1, 6
               SUegg(i) = 0.d0
               SDegg(i) = 0.d0
               SLegg(i) = 0.d0
            enddo

            SNegg(1) = 0.d0
            SNegg(2) = 0.d0
            SNegg(3) = 0.d0

            mA0sq = 0.d0
            mh0sq = 0.d0
            mHu0sq = 0.d0
            mHpmsq = 0.d0    

            Neg(1) = 0.d0
            Neg(2) = 0.d0
            Neg(3) = 0.d0
            Neg(4) = 0.d0
            
            Ceg(1) = 0.d0
            Ceg(2) = 0.d0

            mur = 0.d0


            return
         endif tac

      endif exitp



!-----------------------------------------
      msusy = msusyold 

      
!      call runtomz(MX,msusy,mu_conv,bmur_conv,
!     $     murgemz,bmurgemz,newtbetamz,flags)


!      if(flags.eq.'variable underflow ')then
!         return
!      endif


!      num = 1
!      if((Sueggz(1).lt.0.d0).OR.(Sueggz(2).lt.0.d0).OR.
!     $     (Sueggz(3).lt.0.d0).OR.(Sueggz(4).lt.0.d0).OR.
!     $     (Sueggz(5).lt.0.d0).OR.(Sueggz(6).lt.0.d0))then

!         Sueggz(1) = dabs(Sueggz(1))
!         Sueggz(2) = dabs(Sueggz(2))
!         Sueggz(3) = dabs(Sueggz(3))
!         Sueggz(4) = dabs(Sueggz(4))
!         Sueggz(5) = dabs(Sueggz(5))
!         Sueggz(6) = dabs(Sueggz(6))

!         tacsupz = 1

!      else

!         continue
!      endif

!!-----
!!     Down sector

!      if((SDeggz(1).lt.0.d0).OR.(SDeggz(2).lt.0.d0).OR.
!     $     (SDeggz(3).lt.0.d0).OR.(SDeggz(4).lt.0.d0).OR.
!     $     (SDeggz(5).lt.0.d0).OR.(SDeggz(6).lt.0.d0))then

!         SDeggz(1) = dabs(SDeggz(1))
!         SDeggz(2) = dabs(SDeggz(2))
!         SDeggz(3) = dabs(SDeggz(3))
!         SDeggz(4) = dabs(SDeggz(4))
!         SDeggz(5) = dabs(SDeggz(5))
!         SDeggz(6) = dabs(SDeggz(6))
!         
!         tacsdnz = 1
!         
!      else

!         continue
!      endif

!!-------
!!     lepton sector

!      if((SLeggz(1).lt.0.d0).OR.(SLeggz(2).lt.0.d0).OR.
!     $     (SLeggz(3).lt.0.d0).OR.(SLeggz(4).lt.0.d0).OR.
!     $     (SLeggz(5).lt.0.d0).OR.(SLeggz(6).lt.0.d0))then

!         SLeggz(1) = dabs(SLeggz(1))
!         SLeggz(2) = dabs(SLeggz(2))
!         SLeggz(3) = dabs(SLeggz(3))
!         SLeggz(4) = dabs(SLeggz(4))
!         SLeggz(5) = dabs(SLeggz(5))
!         SLeggz(6) = dabs(SLeggz(6))

!         tacslpz = 1
!         
!      else

!         continue
!      endif
!!-------------------------
!!     sneutrinos

!      if((SNeggz(1).lt.0.d0).OR.(SNeggz(2).lt.0.d0).OR.
!     $     (SNeggz(3).lt.0.d0))then

!         SNeggz(1) = dabs(SNeggz(1))
!         SNeggz(2) = dabs(SNeggz(2))
!         SNeggz(3) = dabs(SNeggz(3))
!         
!         tacsnuz = 1
!         
!      else

!         continue
!      endif

!!------------------------
!C     higgs FLAG-



!      if(mh0sqz.lt.0.d0.or.mhu0sqz.lt.0.d0.or.mhpmsqz.lt.0.d0.or.
!     $     mA0sqz.lt.0.d0)then

!         mh0sqz = dabs(mh0sqz)
!         mhu0sqz = dabs(mhu0sqz)
!         mhpmsqz = dabs(mhpmsqz)
!         mA0sqz  = dabs(mA0sqz)
!         
!         tachiggsz = 1

!      else
!         continue
!      ENDIF


!!      print*,"at MZ",tacsupz,tacsdnz,tacslpz,tacsnuz,tachiggsz
!---------------------------------------------------------------
!      alpha1dr=alph1MSbar
!      alpha2dr= alph2MSbar/(1.d0-alph2MSbar/(6.d0*pi))
!      alpha3dr= alph3MSbar/(1.d0-alph3MSbar/(4.d0*pi))       
      
!      gpdr=Sqrt(3.d0/5.d0)*Sqrt(4.d0*pi*alpha1dr)
!      gdr=Sqrt(4.d0*pi*alpha2dr)
!      g1MS=Sqrt(4.d0*pi*alph1MSbar)
!      g2MS=Sqrt(4.d0*pi*alph2MSbar)
!      sinsqthwdr_susy= gpdr**2.d0/((gpdr**2.d0+gdr**2.d0))
!      alphemdr=sinsqthwdr_susy*alpha2dr
!      alphemMS= ((g1MS*g2MS)**2.d0)/(4.d0*pi*(g1MS**2.d0+g2MS**2.d0))
!      alphemdr=alphemMS/(1.d0-alphemMS/(6.d0*pi))
!      sinthwMS=g1MS/sqrt(g1MS**2.d0+g2MS**2.d0)
!      delrMS=1.d0-pi*alphemMS/(sqrt(2.d0)*GF*MZ*MZ*sinthwMS**2.d0*
!     $ (1.d0-sinthwMS**2.d0))
!      sincorrec=pi*alphemdr/(sqrt(2.d0)*MZ*MZ*GF*(1.d0-delrMS))
!      sinthwsqdr=0.5 -sqrt(0.25 - sincorrec)
!      vevdr= sqrt(MZ*MZ*sinthwsqdr*(1.d0-sinthwsqdr)/(pi*alphemdr))
!       vevdr=sqrt(4.d0*MZ**2.d0/(gpdr**2.d0+gdr**2.d0))
            
      
      call SVD(3, 3, yuRG,3, yudrdig, UQudr,3, UTUdr,3, 0)
      call SVD(3, 3, ydRG,3, yddrdiag, UQddr,3, UTDdr,3, 0)

      call dag(UQddr,UQddrT)
      call matmult(UQudr,UQddrT,VCKMdrbar)
!!!!change by priyanka
!---------------------------------------------------------------------------
! subtract the self energy corrections and threshould correction
!---------------------------------------------------------------------------

      alph3drbar = alph3*(4.d0 * pi)  
      alph2drbar = alph2*(4.d0 * pi)  
      alph1drbar = alph1*(4.d0 * pi)  
      
      gpdr= Sqrt(3.d0/5.d0)*Sqrt(4.d0*pi*alph1drbar)
      gdr = Sqrt(4.d0*pi*alph2drbar)

!-----DR to MS bar conversion of gauge couplings , sinsqthw and vev
      alpha1MS = alph1drbar
      alpha2MS= alph2drbar/(1.d0+alph2drbar/(6.d0*pi))
      alpha3MS= alph3drbar/(1.d0+alph3drbar/(4.d0*pi))                    
      sinsqthwdr_susy= gpdr**2.d0/((gpdr**2.d0+gdr**2.d0))
      alphemdr=sinsqthwdr_susy*alpha2dr
!-----------DR to MS conversion of fermion masses-------------------------------------
      delupdrdown =1.d0+(alph3drbar/(3.d0*pi)+4.d0*alphemdr/(36.d0*pi)+
     $          43.d0*alph3drbar**2.d0/(144.d0*pi*pi))      
      deldowndrdown=1.d0+alph3drbar/(3.d0*pi)+1.d0*alphemdr/(36.d0*pi)+
     $          43.d0*alph3drbar**2.d0/(144.d0*pi*pi)  
      delelecdrdown = 1.d0+alphemdr/(4.d0*pi)
      
      MtMS = (yuRG(3,3)*vev2/Sqrt(2.d0))*delupdrdown
      McMS = (yuRG(2,2)*vev2/Sqrt(2.d0))*delupdrdown
      MuMS = (yuRG(1,1)*vev2/Sqrt(2.d0))*delupdrdown
      
      MbMS = (yddrdiag(3)*vev1/Sqrt(2.d0))*deldowndrdown
      MsMS = (yddrdiag(2)*vev1/Sqrt(2.d0))*deldowndrdown
      MdMS = (yddrdiag(1)*vev1/Sqrt(2.d0))*deldowndrdown
      
      MtauMS = (yeRG(3,3)*vev1/Sqrt(2.d0))*delelecdrdown
      MmuMS = (yeRG(2,2)*vev1/Sqrt(2.d0))*delelecdrdown
      MeMS = (yeRG(1,1)*vev1/Sqrt(2.d0))*delelecdrdown

      If(itcount.eq.1)then      
!      g1MS=Sqrt(3.d0/5.d0)*Sqrt(4.d0*pi*alpha1MS)
!      g2MS=Sqrt(4.d0*pi*alpha2MS)      
!      vevMS=sqrt(4.d0*MZ**2.d0/(g1MS**2.d0+g2MS**2.d0)) 
       vevMS=sqrt(vev1**2.d0+vev2**2.d0)     
!-------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
      else
      
      alpha1susy=(alpha1MS/(1.0d0-alpha1MS*deltasusy1int))
      alpha2susy=(alpha2MS/(1.0d0-alpha2MS*deltasusy2int))
      alpha3susy=(alpha3MS/(1.0d0-alpha3MS*deltasusy3int))
    
      
!      gpdr=Sqrt(3.d0/5.d0)*Sqrt(4.d0*pi*alpha1dr)
!      gdr=Sqrt(4.d0*pi*alpha2dr)  
!      g3 = Sqrt(4.d0*pi*alpha3dr)
      gpMS=Sqrt(3.d0/5.d0)*Sqrt(4.d0*pi*alpha1susy)
      g2MS=Sqrt(4.d0*pi*alpha2susy)
      g3MS=Sqrt(4.d0*pi*alpha3susy)
      sinsqthw_susy= gpMS**2.d0/((gpMS**2.d0+g2MS**2.d0))
      
!      mt_susy=(yuMSbar(3,3)*vevMSbar/Sqrt(2.d0))
!      mb_susy=(ydMSbar(3,3)*vevMSbar/Sqrt(2.d0))
!      mtau_susy = (yeMSbar(3,3)*vevMSbar/Sqrt(2.d0))
     
      
 
      MTc_msusy =  correction_susy

      mBc_msusy =  mbcor_susy

      mTauc_msusy = mtaucor_susy


      MtMS = MtMS* (1.d0 - MTc_msusy)

      MbMS = MbMS / (1.d0 - mBc_msusy)

      MtauMS = MtauMS * (1.d0 - mtauc_msusy)

      MZ_susyd =(0.5)* dsqrt((gpMS**2.d0+g2MS**2.d0) *
     $           (vev1**2.d0 + vev2**2.d0))

      p = MZpole
 
      q = msusyold

      call pizz(p,q,g2MS,mtMS,mbMS,mtauMS,newtbeta,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     ON,OCL,OCR,sinsqthw_susy,pizzT_susyd)   
      
      vevMS=sqrt(4.d0*(MZ_susyd**2.d0-pizzT_susyd)/
     $      (gpMS**2.d0+g2MS**2.d0)) 
      endif      
              
!------define yukawa and gauge couplings at msusy---------------
      vevscMS=vevMS/sqrt(2.d0)
 
      ydin(1,1) = (VCKMdrbar(1,1)*mDMS)/(4.d0*pi*vevscMS)
      ydin(1,2) = (VCKMdrbar(1,2)*mSMS)/(4.d0*pi*vevscMS)
      ydin(1,3) = (VCKMdrbar(1,3)*mBMS)/(4.d0*pi*vevscMS)
      ydin(2,1) = (VCKMdrbar(2,1)*mDMS)/(4.d0*pi*vevscMS)
      ydin(2,2) = (VCKMdrbar(2,2)*mSMS)/(4.d0*pi*vevscMS)
      ydin(2,3) = (VCKMdrbar(2,3)*mBMS)/(4.d0*pi*vevscMS)
      ydin(3,1) = (VCKMdrbar(3,1)*mDMS)/(4.d0*pi*vevscMS)
      ydin(3,2) = (VCKMdrbar(3,2)*mSMS)/(4.d0*pi*vevscMS)
      ydin(3,3) = (VCKMdrbar(3,3)*mBMS)/(4.d0*pi*vevscMS) 


      yein(1,1) = mEMS/(4.d0*pi*vevscMS)
      yein(1,2) = 0.d0
      yein(1,3) = 0.d0
      yein(2,1) = 0.d0
      yein(2,2) = mMUMS/(4.d0*pi*vevscMS)
      yein(2,3) = 0.d0
      yein(3,1) = 0.d0
      yein(3,2) = 0.d0
      yein(3,3) = mTauMS/(4.d0*pi*vevscMS)

      yuin(1,1) = mUMS/(4.d0*pi*vevscMS)
      yuin(1,2) = 0.d0
      yuin(1,3) = 0.d0
      yuin(2,1) = 0.d0
      yuin(2,2) = mCMS/(4.d0*pi*vevscMS)
      yuin(2,3) = 0.d0
      yuin(3,1) = 0.d0
      yuin(3,2) = 0.d0
      yuin(3,3) = MTMS/(4.d0*pi*vevscMS)      
      

      tqsusy   =  dLog(MX**2.d0/msusy**2.d0)
C-----------------------------
C	top Yukawa
C-----------------------------

      i0 = 0

      do i = 1, 3
      yy_sm(i0 + i)   = yuin(1,i)
      yy_sm(i0+3 + i) = yuin(2,i)
      yy_sm(i0+6 + i) = yuin(3,i)
      enddo

C-------------------------------
C	bottom Yukawa
C-------------------------------

      i0 = 9

      do i = 1, 3
      yy_sm(i0 + i)   = ydin(1,i)
      yy_sm(i0+3 + i) = ydin(2,i)
      yy_sm(i0+6 + i) = ydin(3,i)
      enddo 


C-----------------------------
C     tau Yukawa
C-----------------------------

      i0 = 18

      do i = 1, 3
      yy_sm(i0 + i)   = yein(1,i)
      yy_sm(i0+3 + i) = yein(2,i)
      yy_sm(i0+6 + i) = yein(3,i)
      enddo 

      yy_sm(28) = alpha3susy/(4.d0 * pi)
      yy_sm(29) = alpha2susy/(4.d0 * pi)
      yy_sm(30) = alpha1susy/(4.d0 * pi)
      yy_sm(31) = vevMS

!--------------------------------------------------------------------------
!     RUNNING UP THE SM YUKAWAS AND GAUGE COUPLINGS
!--------------------------------------------------------------------------
      n0 = 31
      check=0
      x2 = tz
      x1 = tqsusy
      h1   =  -1.d-5
      hmin =  2.d-8
      eps  =  1.d-6
         call RK4ROUTINE(yy_sm,n0,x1,x2,eps,h1,hmin,nok,nbad,smrgemt,
     $        QMSRK4,check)
         if(check.eq.100)then
            flags = 'variable underflow '
            return
         endif

      mt_mz = yy_sm(9)*(4.d0*pi)*yy_sm(31)/dsqrt(2.d0)
      mb_mz = yy_sm(18)*(4.d0*pi)*yy_sm(31)/dsqrt(2.d0)
      mtau_mz = yy_sm(27)*(4.d0*pi)*yy_sm(31)/dsqrt(2.d0)
      alpha1MZ = (4.d0 * pi)*yy_sm(30)
      alpha2MZ =  (4.d0 * pi)*yy_sm(29)
      alpha3MZ =  (4.d0 * pi)*yy_sm(28)
      
      i0 = 0

      do i = 1,3

      yuMZ(1,i) = (4.d0*pi)*yy_sm(i0 + i)
      j = 3 + i
      yuMZ(2,i) = (4.d0*pi)*yy_sm(i0 + j)
      k = 6 + i
      yuMZ(3,i) = (4.d0*pi)*yy_sm(i0 + k)

      enddo 


C     Bottom Yukawa !!
C     ----------------------------------

      i0 = 9

      do i = 1,3

      ydMZ(1,i) = (4.d0*pi)*yy_sm(i0 + i)
      j = 3 + i
      ydMZ(2,i) = (4.d0*pi)*yy_sm(i0 + j)
      k = 6 + i
      ydMZ(3,i) = (4.d0*pi)*yy_sm(i0 + k)

      enddo 



C     Tau Yukawa !!!
C     ----------------------------------

      i0 = 18

      do i = 1,3

      yeMZ(1,i) = (4.d0*pi)*yy_sm(i0 + i)
      j = 3 + i
      yeMZ(2,i) = (4.d0*pi)*yy_sm(i0 + j)
      k = 6 + i
      yeMZ(3,i) = (4.d0*pi)*yy_sm(i0 + k)

      enddo  

!---------------------------------------------------------------
!-------------------------------------------------
C     reinitializing one loop correction to SM inputs to zero

      MTc_mz = 0.d0
      mBc_mz = 0.d0
      mTauc_mz = 0.d0
      delalphas = 0.d0
      delalphem = 0.d0
!      newvev = vevsc*dsqrt(2.d0)
!      newvev = yy_sm(31)
!---------------------------------------------------------------
!---modified by priyanka--------added SM threshold correction to gauge couplings

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      call coratMZ(sgnmu,MW,MZ,newtbetamz,alphaDR,MWc_mz,MZc_mz,
!     $     MTc_mz,mBc_mz,mTauc_mz,alphas1,alphaem,delalphas,
!     $     delalphem,sinsqtheff,newvev,mbdrbar,flags)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      if((tacsupz+tacsdnz+tacslpz+tacsnuz+tachiggsz).eq.0)then

!      call coratMZ(sgnmu,MW,MZ,newtbetamz,alphaDR,MWc_mz,MZc_mz,
!     $     MTc_mz,mBc_mz,mTauc_mz,alphas1,alphaem,delalphas,
!     $     delalphem,sinsqtheff,newvev,mbdrbar,flags)

!      else


!         flags = "TACSPECMZ"
!         
!!         print*,"Tachyonic Spectrum at MZ"
!         
!         MTc_mz = 0.d0
!         mBc_mz = 0.d0
!         mTauc_mz = 0.d0
!         delalphas = 0.d0
!         delalphem = 0.d0

!      endif

!------------------------------------
C     itcount  => exit procedure

      ex: if(exitcalc.eq.'T')then

         Suegg(1) = STeg(1)
         Suegg(2) = STeg(2)
         Suegg(3) = SCeg(1)
         Suegg(4) = SCeg(2)
         Suegg(5) = SUqeg(1)
         Suegg(6) = SUqeg(2)

         SDegg(1) = SBeg(1)
         SDegg(2) = SBeg(2)
         SDegg(3) = SSTeg(1)
         SDegg(4) = SSTeg(2)
         SDegg(5) = SDneg(1)
         SDegg(6) = SDneg(2)

         SLegg(1) = STaueg(1)
         SLegg(2) = STaueg(2)
         SLegg(3) = SMueg(1)
         SLegg(4) = SMUeg(2)
         SLegg(5) = SEeg(1)
         SLegg(6) = SEeg(2)

         SNegg(1) = tsnu*tsnu
         SNegg(2) = musnu*musnu
         SNegg(3) = elsnu*elsnu

         Ceg(1) = Cheg(1)
         Ceg(2) = Cheg(2)

         Neg(1) = neuteg(1)
         Neg(2) = neuteg(2)
         Neg(3) = neuteg(3)
         Neg(4) = neuteg(4)

         mh0sq = newmh0sq
         mA0sq = newmA0sq
         mhpmsq = newmhpmsq
         mHu0sq = newmHu0sq


         return
      else ex
         continue
      endif ex
      


!-------------------------

!      alph1no = (alphaDR/(4.d0 * pi * (1.d0 - delalphem) * 
!     $     (1.d0 - sinsqtheff))) * (5.d0/3.d0) 
!      alph2no = alphaDR/(4.d0 * pi * (1.d0 - delalphem) * sinsqtheff)
!      alph3no = alph3in/(1.d0 - delalphas) 

!!!!!!!!!!!!fermion Sm one-loop corrections by priyanka------
      g  = dsqrt(alpha2MZ*16.d0*pi*pi)
      gp = dsqrt(alpha1MZ*16.d0*pi*pi) * dsqrt(3.d0/5.d0) 
      g3 = dsqrt(alpha3MZ*16.d0*pi*pi)
      sinsqthw_Mz= gp**2.d0/(gp**2.d0+g**2.d0)
      p = mtpole

      q = MZpole 

!      call topcor(p,q,g,gp,g3,M3tmz,mt_mz,mb_mz,tanbeta,
!     $     yuMZ,ydMZ,SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
!     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz,ONz,OCLz,
!     $     OCRz,sinsqthw_mz,correction)

      call topcorSM(p,q,g,gp,g3,M3tz,mt_mz,mb_mz,tanbeta,
     $     yuMZ,ydMZ,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,ON,OCL,
     $     OCR,sinsqthw_mz,correction_Mz)
!----------------------------------------------------------------------------


      p = mbmzMSbar

      q = MZpole

      call bottomcorSM(p,q,g,gp,g3,M3tz,mt_mz,mb_mz,tanbeta,
     $     yuMZ,ydMZ,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,ON,OCL,
     $     OCR,sinsqthw_mz,mbcor_Mz,mbdrbar)
      

!----------------------------------------------------------------------------

      p = mtaupole

      q = MZpole

      call taucorSM(p,q,g,gp,M3tz,mtau_mz,tanbeta,yeMZ,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,
     $     mhu0sq,mhpmsq,mA0sq,ON,OCL,
     $     OCR,sinsqthw_mz,mtaucor_Mz)

!----------------------------------------------------------------------------
      MTc_mz = correction_Mz
      mBc_mz = mbcor_Mz 
      mtauc_mz = mtaucor_Mz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      mtno = MT * (1.d0 + MTc_mz)
      mbno = mb / (1.d0 + mBc_mz)
      mtauno = mtau * (1.d0 + mtauc_mz)

      p = MZpole

      q = MZpole

      call pizzSM(p,q,g,mtpole,mb_mz,mtau_mz,tanbeta,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,ON,OCL,OCR,
     $     sinsqthw_mz,pizzTSM)

!      print*,"pizzT = ", pizzT

      MZc_mz = dsqrt(MZpole*MZpole + pizzTSM)
      newvev = sqrt(4.d0*MZc_mz**2.d0/(gp**2.d0+g**2.d0)) 

!!!!!!!Priyanka correction ends

c$$$      print*,"mt corr = ", mtno, mtc_mz, mtpole
c$$$      print*,"mb corr = ", mbno, mbc_mz, mb, mbc_mz
c$$$      print*,"mtau corr = ", mtauno, mtauc_mz, mtaupole

      ratio = 1 - min(mtn,mtno)/max(mtn,mtno) 

      ratio1 = 1 - min(mbn,mbno)/max(mbn,mbno) 


!      print*,"ratiomt, ratiomb = ", ratio, ratio1

      stopratd(itcount) = mh2mz

      stopratu(itcount) = mu_conv 

!----------------

      c1: if(itcount.ne.1)then

         ratiou = dabs(1 - min(stopratu(itcount),stopratu(itcount-1))/
     $        max(stopratu(itcount),stopratu(itcount-1))) 

         ratiod = 1 - min(stopratd(itcount),stopratd(itcount-1))/
     $        max(stopratd(itcount),stopratd(itcount-1)) 

!         print*,"ratiomu = ", ratiou, mu_conv 

         c2: if(itcount.eq.30.and.(muflag.eq.2.or.
     $        (ratiou.gt.spectol)))then

 97         format(1pe11.4,2x,1pe11.4,2x,1pe11.4,2x,1pe11.4)


            c3: if(flags.eq.' MUNOC'.or.flags.eq.' REWSB')then
               
               continue

            else c3

               flags = " FSNC"

!               write(15,97) m0,m12,tanbeta,a0


            endif c3

            exitcalc = 'T'

         else c2

            c4: if(ratiou.le.spectol.and.muflag.ne.2)then

               
               c5: if(flag_bmu.ne.' AOK')then
                  flags = flag_bmu                 
               endif c5

               
               exitcalc = 'T'
               
            else c4

               runum = runum + 1
               itcount = itcount + 1
               
            endif c4


         endif c2

      else c1

         itcount = itcount + 1
         
      endif c1

!-----------------------------------

!      yuin(1,1) = ((VCKM(1,1)*mUQ)/(4.d0*pi*newvev/dsqrt(2.d0)))
!      yuin(1,2) = (mUQ*VCKM(1,2))/(4.d0*pi*newvev/dsqrt(2.d0))
!      yuin(1,3) = (VCKM(1,3)*mUQ)/(4.d0*pi*newvev/dsqrt(2.d0))
!      yuin(2,1) = (mC*VCKM(2,1))/(4.d0*pi*newvev/dsqrt(2.d0))
!      yuin(2,2) = (mC*VCKM(2,2))/(4.d0*pi*newvev/dsqrt(2.d0))
!      yuin(2,3) = (mC*VCKM(2,3))/(4.d0*pi*newvev/dsqrt(2.d0))    
!      yuin(3,1) = (MTno*VCKM(3,1))/(4.d0*pi*newvev/dsqrt(2.d0))
!      yuin(3,2) = (MTno*VCKM(3,2))/(4.d0*pi*newvev/dsqrt(2.d0))
!      yuin(3,3) = (MTno*VCKM(3,3))/(4.d0*pi*newvev/dsqrt(2.d0))


!      ydin(1,1) = mD/(4.d0*pi*newvev/dsqrt(2.d0))
!      ydin(1,2) = 0.d0
!      ydin(1,3) = 0.d0
!      ydin(2,1) = 0.d0
!      ydin(2,2) = mS/(4.d0*pi*newvev/dsqrt(2.d0))
!      ydin(2,3) = 0.d0
!      ydin(3,1) = 0.d0
!      ydin(3,2) = 0.d0
!      ydin(3,3) = mBno/(4.d0*pi*newvev/dsqrt(2.d0)) 


!      yein(1,1) = mE/(4.d0*pi*newvev/dsqrt(2.d0))
!      yein(1,2) = 0.d0
!      yein(1,3) = 0.d0
!      yein(2,1) = 0.d0
!      yein(2,2) = mMU/(4.d0*pi*newvev/dsqrt(2.d0))
!      yein(2,3) = 0.d0
!      yein(3,1) = 0.d0
!      yein(3,2) = 0.d0
!      yein(3,3) = mtauno/(4.d0*pi*newvev/dsqrt(2.d0))

      yuin(1,1) = mUQ/(4.d0*pi*newvev/dsqrt(2.d0))
      yuin(1,2) = 0.d0
      yuin(1,3) = 0.d0
      yuin(2,1) = 0.d0
      yuin(2,2) = mC/(4.d0*pi*newvev/dsqrt(2.d0))
      yuin(2,3) = 0.d0    
      yuin(3,1) = 0.d0
      yuin(3,2) = 0.d0
      yuin(3,3) = MTno/(4.d0*pi*newvev/dsqrt(2.d0))


      ydin(1,1) = (VCKM(1,1)*mD)/(4.d0*pi*newvev/dsqrt(2.d0))
      ydin(1,2) = (VCKM(1,2)*mS)/(4.d0*pi*newvev/dsqrt(2.d0)) 
      ydin(1,3) = (VCKM(1,3)*mBno)/(4.d0*pi*newvev/dsqrt(2.d0)) 
      ydin(2,1) = (VCKM(2,1)*mD)/(4.d0*pi*newvev/dsqrt(2.d0)) 
      ydin(2,2) = (VCKM(2,2)*mS)/(4.d0*pi*newvev/dsqrt(2.d0))
      ydin(2,3) = (VCKM(2,3)*mBno)/(4.d0*pi*newvev/dsqrt(2.d0)) 
      ydin(3,1) = (VCKM(3,1)*mD)/(4.d0*pi*newvev/dsqrt(2.d0)) 
      ydin(3,2) = (VCKM(3,2)*mS)/(4.d0*pi*newvev/dsqrt(2.d0)) 
      ydin(3,3) = (VCKM(3,3)*mBno)/(4.d0*pi*newvev/dsqrt(2.d0)) 

      yein(1,1) = mE/(4.d0*pi*newvev/dsqrt(2.d0))
      yein(1,2) = 0.d0
      yein(1,3) = 0.d0
      yein(2,1) = 0.d0
      yein(2,2) = mMU/(4.d0*pi*newvev/dsqrt(2.d0))
      yein(2,3) = 0.d0
      yein(3,1) = 0.d0
      yein(3,2) = 0.d0
      yein(3,3) = mtauno/(4.d0*pi*newvev/dsqrt(2.d0))

c$$$      print*,"ytin = ", yuin(3,3), newvev
c$$$      print*,"ybin = ", ydin(3,3)
c$$$      print*,"ytauin = ", yein(3,3)
      

      runum = runum + 1
      msnew = dsqrt(dsqrt(STeg(1))*dsqrt(STeg(2)))
      msusyold = msusynew  
      murgemz= mu_conv
      bmurgemz = bmur_conv     
!      mur = murgemz
!      bmur = bmurgemz
      mur = mu_conv
      bmur = bmur_conv

      vevin = newvev
      
!      print*,"mur, bmur = ", mur, bmur

      call rgeit(MW,MZ,MT,MTc_mz,mB,mBc_mz,mTau,
     $     mtauc_mz,msusyold,vevsc,vevin,vev1in,vev2in,yuin,
     $     ydin,yein,alphaDR,alph1in,alph2in,delalphem,
     $     alph3in,delalphas,mur,bmur,murge,bmurge,prnstat,check,
     $     newtbeta,MTatMZ, msusynew,mursq,try,M3t,
     $     flags,runum,itcount,stopratu,stopratd,sinsqtheff,
     $     exitcalc)


      END SUBROUTINE rgeit
!--------------------------------------------------------------------------
