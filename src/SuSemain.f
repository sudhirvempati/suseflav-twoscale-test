C==============================================================================
****f* SuSeFLAV/SuSemain.f 
* Program : SuSeFLAV: Supersymmetry seesaw and flavor violation calculator
* Version : 1.0
* Website : http://cts.iisc.ernet.in/Suseflav/main.html
* Authors : Debtosh Choudhury  debtosh@cts.iisc.ernet.in 
*           Raghuveer Garani   rgarani@cts.iisc.ernet.in
*           Sudhir Vempati     vempati@cts.iisc.ernet.in
*  NAME
*    subroutine SuSemain
*  SYNOPSIS
*    Main routine for running SuSeFLAV. 
*
*  FUNCTION
*     Computes 1-loop corrected Supersymmetric particle spectrum for a given  
*     set of mSUGRA/NUHM/CNUM inputs. Also, the routine computes branching 
*     ratios and decay rates for rare lfv processes. 
*
*  INPUTS
*
*     prnstat    -    Print Control. 1= print statements 
*     Mg1        -    high energy input for bino
*     Mg2        -    high energy input for  wino
*     Mg3        -     high energy input for gluino
* 
*    Other relevant input parameters are stored in common block
*      common/sminputs/ mbpole, mtaupole, Mtpole
*      common/loops/ lopt,rhn
*      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
*      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
*      common/charinputs/case, model
*
*  RESULT
*     errge      - if any error is encountered errge =1 else =0.
*    Output is written in slha.out
*  EXAMPLE
*
*    SUBROUTINE SuSeFLAV(prnstat,mq20, mq30, mu20, mu30, md20, md30, 
*     $     ml20,ml30, me20,me30, mnu20,mnu30,Mg1,Mg2,Mg3,errge)
*
*  NOTES
*     Common blocks used:
*      common/sminputs/ mbpole, mtaupole, Mtpole
*      common/loops/ lopt,rhn
*      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
*      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
*      common/charinputs/case, model
*      common/rgeinput_high/ MX, M1X,M2X,M3X,mh10,mh20,
*     $     mQ0,mU0,mD0,mE0,mL0,mNU0
*      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
*     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut
*      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
*     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG
*      common/rgeoutput_MZ/ mh1mzz,mh2mzz,M1tmz,M2tmz,M3tmz,
*     $     alph1MZ,alph2MZ,alph3MZ,vev1mz,vev2mz,yumz,ydmz,yemz
*      common/sparticles_susy/SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
*     $     mh0sq,mhu0sq,mhpmsq,mA0sq
*      common/soft_mat_susy/ mSL1, SLeg, USL, SNeg, USN, SUeg, USU,
*     $     SDeg, USD
*      common/qcd_cor/mbmzdrbar,mbMZmsbar
*      common/sparticles_MZ/SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
*     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz
*      common/soft_mat_mz/ mSL1z, SLegz, USLz, SNegz, USNz, SUegz, USUz,
*     $     SDegz, USDz
*      common/softout_mat_mz/mSQRGz,mSURGz,mSDRGz,AURGz,ADRGz,
*     $     AERGz,ANURGz,mSLRGz,mSERGz, mSNURGz,ONz,OCLz,OCRz,
*     $     MCharz, MNeutz
*      common/rgeopt_mz/mt_mz, mb_mz, mtau_mz
*      common/mutaup/rrstaumu,USLsrt,USLTsrt
*      common/gauge/alph,Gf,alphas
*
*     External routines: 
*      external  rgeit
*
* Date of last modification and authors
***
** 15/12/2011: 1900 Hrs. Modified by Sudhir. (removed a write statmement
** in the Higgs flag) 
** 
***
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!==================================================================================================     

       SUBROUTINE SuSeFLAV(prnstat,mq11,mq12,mq13,mq21,mq22,mq23,
     $     mq31,mq32,mq33,mu11,mu12,mu13,mu21,mu22,mu23,mu31,mu32,
     $     mu33,md11,md12,md13,md21,md22,md23,md31,md32,md33,
     $     ml11,ml12,ml13,ml21,ml22,ml23,ml31,ml32,ml33, 
     $     me11,me12,me13,me21,me22,me23,me31,me32,me33 ,
     $     mnu11,mnu12,mnu13,mnu21,mnu22,mnu23,mnu31,mnu32,mnu33,
     $     Mg1,Mg2,Mg3,errge)


      implicit none 
      
      integer i,j,num, try, rhn, errge,errget
      integer prnstat,lopt,check
      integer runum, itcount
      INTEGER lsppos,nlsppos,iunit
      character*1 exitcalc
      character*3 case
      CHARACTER*4 model
      CHARACTER*100 flags,flagt

      real*8 rt,x1,x2,acc,t,t1,lqcd,alsmb

      double precision DeltaTZ,alph,Gf,alphas
      double precision gmsbsusyb, gmsbmess,nhat,gr
      
      double precision tansqtw,alph1tz,beta
      double precision ctw,stw, sgnmu
      double precision bmur,mur,bmurge,murge
      double precision M1X,M2X,M3X,MX,mh1mz,mh2mz,mh10,mh20,mursq
      double precision M1tz,M2tz,M3tz,m10,m20

      double precision m0,m12,tanbeta,a0,ue3
      double precision mh0sq,mhu0sq,mhpmsq,mA0sq
      double precision mq0(3,3),mu0(3,3),MR1,MR2,MR3
      double precision md0(3,3),ml0(3,3),me0(3,3),mnu0(3,3) 
      double precision ANURG(3,3)
      double precision yuRG(3,3),ydRG(3,3),yeRG(3,3),Ynui(3,3)

      double precision SUeg(6),SDeg(6),SLeg(6),SNeg(3),mSL1(6,6)
      double precision mSQRG(3,3),mSURG(3,3),mSDRG(3,3),mSLRG(3,3),
     $     mSERG(3,3),mSNURG(3,3)
      double precision ADRG(3,3),AURG(3,3),AERG(3,3)
      double precision SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      double precision USD(6,6),USL(6,6),USU(6,6),USN(3,3),Bbsg 

      double precision alphasmt,mTatMz, alphaDR
      
      DOUBLE PRECISION alph1,alph2,alph3,newtbeta

      double precision alphsMZdrbar,fbmz,fbmb,alphsmbMSbar,mbMZmsbar
      double precision mbMZdrbar,mTauMZdrbar,alphaemMZdrbar
      DOUBLE PRECISION Mtaupole,Mtpole,Mbpole


      double precision mT, mB, mTau
      double precision vev1,vev2
      double precision vevin, vev1in, vev2in
      double precision vevsc,vev1sc,vev2sc
      double precision alph1in,alph2in,alph3in
      double precision yuin(3,3),ydin(3,3),yein(3,3)

      double precision MChar(2,2)
      double precision MNeut(4,4),
     $     Neg(4),Negt(4),ON(4,4)
      double precision OCR(2,2),OCL(2,2),Ceg(2)

C--------------------------------------------------------------------------
      double precision Bmeg,Bteg,tmugrate,Btmug,scbmeg
      double precision scbtmug,megrate,scbteg
      double precision tegrate,gm2

      double precision mu3erate,Brmu3E,tau3murate,brtau3mu,tau3erate
      double precision brtau3e,mueconverrate,brmueconver
      double precision scbmu3e,scbmueconver,scbtau3mu,scbtau3e,scgminus2
      DOUBLE PRECISION minspec(31), lsp, nlsp 
!-------------------------------------------
      DOUBLE PRECISION  yumz(3,3), ydmz(3,3), yemz(3,3),mAm3
      DOUBLE PRECISION alph3MZ, alph2MZ, alph1MZ, alpha
      DOUBLE PRECISION piwwT,pizzT
      DOUBLE PRECISION  M3t
      DOUBLE PRECISION MTc_mz,mBc_mz,mTauc_mz

C--------------------------------------------------------------------------- 
      DOUBLE PRECISION mq11,mq12,mq13,mq21,mq22,mq23,mq31,mq32,mq33 
      DOUBLE PRECISION mu11,mu12,mu13,mu21,mu22,mu23,mu31,mu32,mu33 
      DOUBLE PRECISION md11,md12,md13,md21,md22,md23,md31,md32,md33 
      DOUBLE PRECISION ml11,ml12,ml13,ml21,ml22,ml23,ml31,ml32,ml33 
      DOUBLE PRECISION me11,me12,me13,me21,me22,me23,me31,me32,me33 
      DOUBLE PRECISION mnu11,mnu12,mnu13,mnu21,mnu22,mnu23,mnu31,mnu32,
     $     mnu33
      DOUBLE PRECISION Mg1,Mg2,Mg3
C-------------------------------------  


      double precision SUegz(6),SDegz(6),SLegz(6),SNegz(3)
      double precision SUeggz(6),SDeggz(6),SLeggz(6),SNeggz(3)
      double precision USDz(6,6),USLz(6,6),USUz(6,6),USNz(3,3) 
      double precision MCharz(2,2)
      double precision MNeutz(4,4),
     $     Negz(4),ONz(4,4)
      double precision OCRz(2,2),OCLz(2,2),Cegz(2)
      DOUBLE PRECISION mSL1z(6,6)

      DOUBLE PRECISION AURGz(3,3),ADRGz(3,3),AERGz(3,3)
      DOUBLE PRECISION mSQRGz(3,3), mSURGz(3,3),mSDRGz(3,3)
      DOUBLE PRECISION mSLRGz(3,3),mSNURGz(3,3), mSERGz(3,3)

      DOUBLE PRECISION ANURGz(3,3),mt_mz, mb_mz, mtau_mz
      double precision mh0sqz,mhu0sqz,mhpmsqz,mA0sqz

!------------------------------------
      DOUBLE PRECISION delalphas, delalphem,mt_qcd_corr,sinsqtheff
      double precision msusyold,msusynew, stopratu(55),
     $     stopratd(55)

      double precision M1tmz,M2tmz,M3tmz, mh1mzz,mh2mzz, vev1mz,vev2mz
!      double precision rrstaumu(2), USLsrt(6,6),USLTsrt(6,6)

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae
      DOUBLE PRECISION delta1z,delta2z, piwwT0,pizzT0
      DOUBLE PRECISION murgemz,ftmz,ftmt,drho
      DOUBLE PRECISION yugut(3,3),ydgut(3,3),yegut(3,3),ynugut(3,3)
      double precision SLegvd(6),SUegvd(6),SDegvd(6)
      double precision USLvd(6,6),USUvd(6,6),USDvd(6,6)
      double precision SNegvd(3),USNvd(3,3),MWsqpole_MZ
      double precision MWpole, MZpole,sinsqtw,MW,MZ,pi

      double precision mTaupoledrbar,mTauMZmsbar

      double precision mgrav

      real MHL,MHH,MHA,MG,MB1,MB2,AAB,MT1,MT2,AAT,MU,M2,
     $     ALPHAH,MSL,MCL,BRBS,BRBD,thetati,thetabi,mtisa,tbisa


!---VCM parameters
      
      integer qmix

      double precision VCKM(3,3),Vud,Vus,Vub,Vcd,Vcs,Vcb 
      double precision Vtd,Vts,Vtb,th12,th23,th13

      complex*8 VCKMW(3,3),Imag
      double precision lckm,Ackm,rhockm,etackm,rhobckm,etabckm
        
      common/VCKMparam/ VCKM
        
!---R parametrizastion

      Double precision DM(3,3),R(3,3),Dk(3,3),UPMNST(3,3),
     $     drd(3,3),UPMNS(3,3)

      double precision thl12,thl23,thl13,UPMNSR(3,3)


      integer tacsup,tasdn,tacslp,tacsnu,tachiggs
      integer tacsupz,tacsdnz,tacslpz,tacsnuz,tachiggsz

      
      common/rpar/ DM,R,Dk


!---
      common/lfvmixing/USLvd,USUvd,USDvd,USNvd
      common/lfveigval/SLegvd,SUegvd,SDegvd,SNegvd

!-----------------------------------------
      common/mascorr/ MT_qcd_corr

      integer pcount
      common/counter_tag/pcount
      common/counter/itcount
!------------
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mtaurmass/ mtaupoledrbar, mTauMZmsbar, mTauMZdrbar
      common/mwpole/ MWpole
      common/loops/ lopt,rhn
      common/quarkmix/ qmix
      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
      common/charinputs/case, model
      common/rgeinput_high/ MX, M1X,M2X,M3X,mh10,mh20,
     $     mQ0,mU0,mD0,mE0,mL0,mNU0

      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut

      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG

      common/rgeoutput_MZ/ mh1mzz,mh2mzz,M1tmz,M2tmz,M3tmz,
     $     alph1MZ,alph2MZ,alph3MZ,vev1mz,vev2mz,yumz,ydmz,yemz

      common/sparticles_susy/SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq

      common/soft_mat_susy/ mSL1, SLeg, USL, SNeg, USN, SUeg, USU,
     $     SDeg, USD

      common/qcd_cor/mbmzdrbar,mbMZmsbar

      common/sparticles_MZ/SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz

      common/soft_mat_mz/ mSL1z, SLegz, USLz, SNegz, USNz, SUegz, USUz,
     $     SDegz, USDz

      common/softout_mat_mz/mSQRGz,mSURGz,mSDRGz,AURGz,ADRGz,
     $     AERGz,ANURGz,mSLRGz,mSERGz, mSNURGz,ONz,OCLz,OCRz,
     $     MCharz, MNeutz

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

      double precision OCRm(2,2),OCLTm(2,2),ONm(4,4),alphahiggs

      common/opt_mixing/OCRm, OCLTm,ONm, alphahiggs

      common/ma/mAm3
      common/yukawagut/yugut,ydgut,yegut,ynugut

      common/rgeopt_mz/mt_mz, mb_mz, mtau_mz
      common/mu_mz/ murgemz

      common/tac_susy/ tacsup,tasdn,tacslp,tacsnu,tachiggs
      common/tac_mz/ tacsupz,tacsdnz,tacslpz,tacsnuz,tachiggsz


!      common/mutaup/rrstaumu,USLsrt,USLTsrt

      common/gauge/alph,Gf,alphas

      common/finetuning/delta1z,delta2z
      common/deltarho/ piwwT0,pizzT0,pizzT,piwwT,MWsqpole_MZ

      common/lowepar/drho,gm2,Bbsg,bmeg,brmu3e,btmug,brtau3mu,
     $     bteg,brtau3e

      common/gmsbinputs/ gmsbsusyb, gmsbmess,gr,nhat

      double precision mh0sqcor,mhu0sqcor,mA0sqcor,mHpmsqcor

      common/m32/mgrav

!--------

      common/higgs2loop/ mh0sqcor,mhu0sqcor,mA0sqcor,mHpmsqcor


!--------------------
      real tstart(2)            
      DOUBLE PRECISION total              ! For receiving total time
      
!-----------------------------------------------
      external  rgeit,gmsb,printslha,dag,matmult,mat3prod
!---------------------------------------------

C     WARNING: stdinputs MUST be AFTER all other declarations, 
C     but BEFORE any other statement

      include 'stdinputs.h'
!--------------------------------

      iunit = 0

      verb01: if(prnstat.eq.1)then

 25      format(2x,A17,2x,ES11.4)
 27      format(2x,A17,2x,A4)
 26      format(2x,A17,2x,I2)
 29      format(2x,A17,2x,A5)
 28      format(2x,1es11.4,2x,1es11.4,2x,1es11.4)  
 30      format(2x,A10)
         
c$$$  write(10,*) ''
c$$$  write(*,*) ''
         write(10,*) '******** Begin Program SuSeFLAV ************'
         write(*,*) '******** Begin Program SuSeFLAV ************'
         write(10,29) ' model = ', model
         write(*,29) ' model = ', model
         write(10,26) ' loop = ', lopt
         write(*,26) ' loop = ', lopt
         write(10,25) ' tanbeta = ',tanbeta
         write(*,25) ' tanbeta = ',tanbeta

         if(model.eq.'GMSB')then
            Write(10,25) ' lambda  = ', gmsbsusyb
            Write(*,25) ' lambda  = ', gmsbsusyb
            Write(10,25) ' Messenger scale  = ', gmsbmess
            Write(*,25) ' Messenger scale  = ', gmsbmess
            Write(10,25) ' number of messenger multiplets  = ', nhat
            Write(*,25) ' number of messenger multiplets  = ', nhat
            write(10,25) ' sign mu = ', sgnmu
            write(*,25) ' sign mu = ', sgnmu
            write(10,26) ' qmix = ', qmix 
            write(*,26) ' qmix = ', qmix
            write(10,25) ' top pole mass = ',Mtpole
            write(*,25) ' top pole mass = ',Mtpole
         endif

         if(model.eq.'mSUG'.or.model.eq.'NUHM')then
            Write(10,25) ' m0  = ', m0
            Write(*,25) ' m0  = ', m0
            Write(10,25) ' a0  = ', a0
            Write(*,25) ' a0  = ', a0
            Write(10,25) ' M12 = ', M12
            Write(*,25) ' M12 = ', M12
            Write(10,25) ' m10 = ', m10
            Write(*,25) ' m10 = ', m10
            Write(10,25) ' m20 = ', m20
            Write(*,25) ' m20 = ', m20
            write(10,25) ' sign mu = ', sgnmu
            write(*,25) ' sign mu = ', sgnmu
            write(10,26) ' qmix = ', qmix 
            write(*,26) ' qmix = ', qmix
            write(10,25) ' top pole mass = ',Mtpole
            write(*,25) ' top pole mass = ',Mtpole
         endif


         if(model.eq.'NUGM')then
            Write(10,25) ' m0  = ', m0
            Write(*,25) ' m0  = ', m0
            Write(10,25) ' a0  = ', a0
            Write(*,25) ' a0  = ', a0
            Write(10,25) ' M1 = ', Mg1
            Write(*,25) ' M1 = ', Mg1
            Write(10,25) ' M2 = ', Mg2
            Write(*,25) ' M2 = ', Mg2
            Write(10,25) ' M3 = ', Mg3
            Write(*,25) ' M3 = ', Mg3
            Write(10,25) ' m10 = ', m10
            Write(*,25) ' m10 = ', m10
            Write(10,25) ' m20 = ', m20
            Write(*,25) ' m20 = ', m20
            write(10,25) ' sign mu = ', sgnmu
            write(*,25) ' sign mu = ', sgnmu
            write(10,26) ' qmix = ', qmix 
            write(*,26) ' qmix = ', qmix
            write(10,25) ' top pole mass = ',Mtpole
            write(*,25) ' top pole mass = ',Mtpole
         endif

         
         if(model.eq.'CNUM')then
            Write(10,30) ' mQ  = '
            Write(10,28)  mq11,mq12,mq13
            Write(10,28)  mq21,mq22,mq23
            Write(10,28)  mq31,mq32,mq33
            Write(10,30) ' mU  = '
            Write(10,28)  mu11,mu12,mu13
            Write(10,28)  mu21,mu22,mu23
            Write(10,28)  mu31,mu32,mu33
            Write(10,30) ' mD  = '
            Write(10,28)  md11,md12,md13
            Write(10,28)  md21,md22,md23
            Write(10,28)  md31,md32,md33
            Write(10,30) ' mL  = '
            Write(10,28)  ml11,ml12,ml13
            Write(10,28)  ml21,ml22,ml23
            Write(10,28)  ml31,ml32,ml33
            Write(10,30) ' mE  = '
            Write(10,28)  me11,me12,me13
            Write(10,28)  me21,me22,me23
            Write(10,28)  me31,me32,me33
            Write(10,30) ' mNu  = '
            Write(10,28)  mnu11,mnu12,mnu13
            Write(10,28)  mnu21,mnu22,mnu23
            Write(10,28)  mnu31,mnu32,mnu33
            Write(*,30) ' mQ  = '
            Write(*,28)  mq11,mq12,mq13
            Write(*,28)  mq21,mq22,mq23
            Write(*,28)  mq31,mq32,mq33
            Write(*,30) ' mU  = '
            Write(*,28)  mu11,mu12,mu13
            Write(*,28)  mu21,mu22,mu23
            Write(*,28)  mu31,mu32,mu33
            Write(*,30) ' mD  = '
            Write(*,28)  md11,md12,md13
            Write(*,28)  md21,md22,md23
            Write(*,28)  md31,md32,md33
            Write(*,30) ' mL  = '
            Write(*,28)  ml11,ml12,ml13
            Write(*,28)  ml21,ml22,ml23
            Write(*,28)  ml31,ml32,ml33
            Write(*,30) ' mE  = '
            Write(*,28)  me11,me12,me13
            Write(*,28)  me21,me22,me23
            Write(*,28)  me31,me32,me33
            Write(*,30) ' mNu  = '
            Write(*,28)  mnu11,mnu12,mnu13
            Write(*,28)  mnu21,mnu22,mnu23
            Write(*,28)  mnu31,mnu32,mnu33
            Write(10,25) ' m10 = ', m10
            Write(*,25) ' m10 = ', m10
            Write(10,25) ' m20 = ', m20
            Write(*,25) ' m20 = ', m20
            Write(10,25) ' a0  = ', a0
            Write(*,25) ' a0  = ', a0
            Write(10,25) ' M1 = ', Mg1
            Write(*,25) ' M1 = ', Mg1
            Write(10,25) ' M2 = ', Mg2
            Write(*,25) ' M2 = ', Mg2
            Write(10,25) ' M3 = ', Mg3
            Write(*,25) ' M3 = ', Mg3
            write(10,25) ' sign mu = ', sgnmu
            write(*,25) ' sign mu = ', sgnmu
            write(10,26) ' qmix = ', qmix 
            write(*,26) ' qmix = ', qmix
            write(10,25) ' top pole mass = ',Mtpole
            write(*,25) ' top pole mass = ',Mtpole
         endif

         if(rhn.eq.1)then
            Write(10,26) ' rhn = ', rhn 
            Write(*,26) ' rhn = ', rhn 
!            Write(10,25) ' Ue3 = ', ue3 
!            Write(*,25) ' Ue3 = ', ue3 
            write(10,27) ' case = ', case
            write(*,27) ' case = ', case
            Write(10,25) ' MR1 = ', MR1
            Write(*,25) ' MR1 = ', MR1
            Write(10,25) ' MR2 = ', MR2
            Write(*,25) ' MR2 = ', MR2
            Write(10,25) ' MR3 = ', MR3
            Write(*,25) ' MR3 = ', MR3
            if(case.eq.'USD')then
               write(10,*) ' Dirac Neutrino matrix at high energy:'
               write(*,*) ' Dirac Neutrino matrix at high energy:'
               opu22: do i = 1, 3
 223           format(1x,3(2x,1es11.4))
               write(10,223) (Ynui(i,j),j = 1, 3)
               write(*,223) (Ynui(i,j),j = 1, 3)
            enddo opu22  
         endif
      else
         continue
      endif

      write(10,*) '********************************************'
      write(*,*) '********************************************'
      write(10,*) ''
      write(*,*) ''
      else if(prnstat.eq.0)then
         continue
      endif verb01
      

!-------------------------------------------------------------------------
!     SM gauge couplings -  MSbar to DRbar 
!-------------------------------------------------------------------------
      
!      print*,"mzpole = ", mzpole

      pi = 4.d0*datan(1.d0)
      sinsqtw = 0.2221d0
      stw = dsqrt(sinsqtw)
      ctw = dsqrt(1.d0 - sinsqtw)
      MWpole = MZpole * ctw
      tansqtw = (stw/ctw)**2.d0

      MW = MWpole
      MZ = MZpole


      alphsMZdrbar = (pi/((pi/alphas) - 0.25d0))

      alpha =  1.d0/((1.d0/alph) - (1.d0/(6.d0*pi)))
      
!      print*,"alpha conv = ", alpha

      alphaemMZdrbar =  alpha/sinsqtw

      alph1tz = (alpha/(ctw**2.d0))*(5.d0/3.d0)

      alphaDR = alpha
      
c$$$      print*,"alpha1tz = ",  (alpha/(ctw**2.d0)),alph1tz
c$$$      print*,"alpha2 = ", alphaemMZdrbar

!----------------------------------------------------------------
!     Top,tau masses  -  Pole mass to DRbar 
!----------------------------------------------------------------

      DeltaTZ = 2.d0 * dLog(Mtpole/MZ)

      alphasmt = alphsMZdrbar /(1.d0+3.d0*alphsMZdrbar*
     $     DeltaTZ/(4.d0* Pi))

      MT = Mtpole * (1.d0 -                                            !<---check expression.
     $     (alphasmt * (5.d0 - 3.d0 * DeltaTZ)/(3.d0*Pi)) -
     $     (alphasmt**2.d0) * (0.538d0 - 43.d0 * DeltaTZ/(24.d0 * pi*pi)
     $     + (3.d0 * (DeltaTZ**2.d0)/(8.d0 * pi * pi))))


!------------------------------

      mTaupoledrbar = mTaupole * (1.d0 - (3.d0/32.d0) * (alph1tz - 
     $     alphaemMZdrbar/1.d0))

      mTauMZmsbar = mTaupole * ( 1.d0 - (alph/pi) * (1.d0 + 
     $     (3.d0/4.d0) * dlog(MZ*MZ/mTaupole**2.d0))) 

c$$$      print*,"mTaupoledrbar,mTauMZmsbar = ", mTaupoledrbar,mTauMZmsbar
c$$$      print*,"mtauMZdrbar = ",  mTauMZmsbar * (1.d0 - (3.d0/32.d0) * 
c$$$     $     (alph1tz -  alphaemMZdrbar/1.d0))


      mTau = Mtaupole

      
!      mTauMZdrbar = mTau * (1 - (3.d0/8.d0) * (alph1tz - 
!     $     alphaemMZdrbar/4.d0))

      mTauMZdrbar = mTauMZmsbar * (1.d0 - (3.d0/32.d0) * (alph1tz - 
     $     alphaemMZdrbar/1.d0))

!-------------------------------------------------------------------
!     Bottom mass  -  MSbar to DRbar       
!     From hep-ph/0207126     (Baer et al)
!-------------------------------------------------------------------

      Mb = Mbpole      
     
      fbmz = ((23.d0 * alphas)/(6.d0 * pi))**(12.d0/23.d0) * (1.d0 + 
     $     (3731.d0 * alphas/(3174.d0 * pi)) + 1.500706d0 * 
     $     (alphas/pi)**2.d0)


!      print*,"root before = ", x1,alphas
      
      x1 = 1.d0
      x2 = 25.d0
      acc = 1.d-6

      t = rt(alphas,x1,x2,acc)

      if(t==0.d0)then
         print*,"Unable to Find LambdaQCD"
         errge = 1
         return
      endif

      lqcd = dsqrt(MZ*MZ/dexp(t)) 

!      print*,"root = ", t, MZ, lqcd !rt(alphas,x1,x2,acc)

      t1 = 2.d0 * dlog(Mbpole/lqcd)

      alphsmbMSbar = (12.d0 * pi * ( 1.d0 + 
     $     (121104.d0 *
     $     (-0.32233865107676046d0 + 
     $     (-0.5d0 + dLog(t1))**2.d0))/
     $     (279841.d0 * t1**2.d0) - 
     $     (348.d0 * dLog(t1))/(529.d0 * t1)))/(23.d0 * t1)

!      print*,"alsmb = ", alphsmbMSbar


      fbmb = ((23.d0 * alphsmbMSbar)/(6.d0 * pi))**(12.d0/23.d0) * 
     $     (1.d0 + (3731.d0 * alphsmbMSbar/(3174.d0 * pi)) + 
     $     1.500706d0 * (alphsmbMSbar/pi)**2.d0)

      mbMZmsbar = mb * (fbmz/fbmb)

      mbMZdrbar = mbMZmsbar * (1.d0 - (alphsMZdrbar/(3.d0 * pi)) - 
     $     (29.d0 * (alphsMZdrbar**2.d0)/(72.d0 * pi * pi)) + 
     $     (3.d0 * alphaemMZdrbar/(32.d0 * pi)) + 
     $     (13.d0 * (alpha/(ctw**2.d0))/(288.d0*pi)))  
!     $     (13.d0 * alph1tz/(288.d0*pi)))

!      print*,"alphas = ", alphsMZdrbar
      

!---------------------------------------------------------------------------

c$$$      print*,"mb = ", mbMZdrbar,mbpole,mbMZmsbar
c$$$      print*,"mt = ", mT,mtpole
c$$$      print*,"mtau = ",  mTauMZdrbar,mtaupole
      
      mTatMZ = MT               ! storing the value for later use
      mb = mbMZdrbar
      mTau = mTauMZdrbar

!---------------------------------------------------------------------------
!     ANGLES AND CKM MATRICES
!---------------------------------------------------------------------------

c$$$      th12 = dASin(0.2272d0)
c$$$      th23 = dASin(0.04221d0)
c$$$      th13 = dASin(0.003959d0)

      th12 = dASin(0.2246d0)
      th23 = dASin(0.0420d0)
      th13 = dASin(0.0035d0)

      Vud = dCos(th12)*dCos(th13)
      Vus = dSin(th12)*dCos(th13)
      Vub = dSin(th13) 

      Vcd = - dSin(th12)*dCos(th23) - dCos(th12)*dSin(th23)
     $     *dSin(th13)
      Vcs = dCos(th12)*dCos(th23) - dSin(th12)*dSin(th23)*dSin(th13)
      Vcb = dCos(th13)*dSin(th23)

      Vtd = dSin(th12)*dSin(th23)- dCos(th12)*dCos(th23)*dSin(th13)
      Vts = -dCos(th12)*dSin(th23)-dSin(th12)*dCos(th23)*dSin(th13) 
      Vtb = dCos(th23)*dCos(th13)


      VCKM(1,1) = Vud
      VCKM(1,2) = Vus
      VCKM(1,3) = Vub

      VCKM(2,1) = Vcd
      VCKM(2,2) = Vcs
      VCKM(2,3) = Vcb

      VCKM(3,1) = Vtd
      VCKM(3,2) = Vts
      VCKM(3,3) = Vtb

c$$$      print*,"CKM Matrix in Standard Parametrization"
c$$$
c$$$      stnckm: do i = 1, 3
c$$$!     223  format(1x,3(2x,1es11.4))
c$$$      write(*,223) (VCKM(i,j),j = 1, 3)
c$$$      enddo stnckm  



!--------------------------------------------------------------------------
!     Wolfenstein Parametrization
!--------------------------------------------------------------------------


      Imag = (0.d0, 1.d0)
      
      Ackm = 0.812d0
      lckm = 0.22543d0
      rhobckm = 0.145d0
      etabckm = 0.343d0

      rhockm = rhobckm/(1.d0 - lckm*lckm/2.d0)
      etackm = etabckm/(1.d0 - lckm*lckm/2.d0)

c$$$      VCKMW(1,1) = 1.d0 - lckm*lckm/2.d0
c$$$      VCKMW(1,2) = lckm
c$$$      VCKMW(1,3) = Ackm*lckm**3.d0*(rhockm - Imag * etackm)
c$$$      
c$$$      VCKMW(2,1) = - lckm
c$$$      VCKMW(2,2) = 1.d0 - lckm*lckm/2.d0
c$$$      VCKMW(2,3) = Ackm*lckm**2.d0
c$$$
c$$$      VCKMW(3,1) = Ackm*lckm**3.d0*(1.d0 - rhockm - Imag * etackm)
c$$$      VCKMW(3,2) = - Ackm*lckm**2.d0
c$$$      VCKMW(3,3) = 1.d0

!      qmix = 1

!      print*,"qmix = ", qmix

      if(qmix>0)then

         VCKM(1,1) = 1.d0 - lckm*lckm/2.d0
         VCKM(1,2) = lckm
         VCKM(1,3) = Ackm*lckm**3.d0*
     $        dsqrt(rhockm**2.d0+etackm**2.d0)
         
         VCKM(2,1) = - lckm
         VCKM(2,2) = 1.d0 - lckm*lckm/2.d0
         VCKM(2,3) = Ackm*lckm**2.d0

         VCKM(3,1) = Ackm*lckm**3.d0*
     $        dsqrt((1.d0 - rhockm)**2.d0+etackm**2.d0)
         VCKM(3,2) = - Ackm*lckm**2.d0
         VCKM(3,3) = 1.d0

      else

         VCKM(1,1) = 1.d0
         VCKM(1,2) = 0.d0
         VCKM(1,3) = 0.d0
         
         VCKM(2,1) = 0.d0
         VCKM(2,2) = 1.d0
         VCKM(2,3) = 0.d0

         VCKM(3,1) = 0.d0
         VCKM(3,2) = 0.d0
         VCKM(3,3) = 1.d0

      endif

c$$$      print*,"CKM Matrix in Wolfenstein Parametrization"
c$$$
c$$$      wolckm: do i = 1, 3
c$$$ 2024 format(1x,6(2x,1es11.4))
c$$$      write(*,2024) (VCKM(i,j),j = 1, 3)
c$$$      enddo wolckm  


!---------------------------------------------------------------------------
!     Yukawa RGE Inputs
!---------------------------------------------------------------------------
      
      beta = datan(tanbeta)

c$$$      vevin =  246.d0 

      vevin =  dsqrt(MZ*MZ/((alphaemMZdrbar + 
     $     (alph1tz * (3.d0/5.d0)))*pi))

c$$$      print*,'vev from g', dsqrt(MZ*MZ/((alphaemMZdrbar + 
c$$$     $     (alph1tz * (3.d0/5.d0)))*pi))

      vev1in = vevin * dcos(beta)
      vev2in = vevin * dsin(beta)

      vevsc = vevin/dsqrt(2.d0)
      vev1sc = vev1in/dsqrt(2.d0)
      vev2sc = vev2in/dsqrt(2.d0)

      ydin(1,1) = mD/(4.d0*pi*vevsc)
      ydin(1,2) = 0.d0
      ydin(1,3) = 0.d0
      ydin(2,1) = 0.d0
      ydin(2,2) = mS/(4.d0*pi*vevsc)
      ydin(2,3) = 0.d0
      ydin(3,1) = 0.d0
      ydin(3,2) = 0.d0
      ydin(3,3) = mB/(4.d0*pi*vevsc) 

      yein(1,1) = mE/(4.d0*pi*vevsc)
      yein(1,2) = 0.d0
      yein(1,3) = 0.d0
      yein(2,1) = 0.d0
      yein(2,2) = mMU/(4.d0*pi*vevsc)
      yein(2,3) = 0.d0
      yein(3,1) = 0.d0
      yein(3,2) = 0.d0
      yein(3,3) = mTau/(4.d0*pi*vevsc)

      yuin(1,1) = (mUQ*VCKM(1,1))/(4.d0*pi*vevsc)
      yuin(1,2) = (mUQ*VCKM(1,2))/(4.d0*pi*vevsc)
      yuin(1,3) = (mUQ*VCKM(1,3))/(4.d0*pi*vevsc)
      yuin(2,1) = (mC*VCKM(2,1))/(4.d0*pi*vevsc)
      yuin(2,2) = (mC*VCKM(2,2))/(4.d0*pi*vevsc)
      yuin(2,3) = (mC*VCKM(2,3))/(4.d0*pi*vevsc)

c$$$      yuin(3,1) = (MTpole*VCKM(3,1))/(4.d0*pi*vevsc)
c$$$      yuin(3,2) = (MTpole*VCKM(3,2))/(4.d0*pi*vevsc)
c$$$      yuin(3,3) = (MTpole*VCKM(3,3))/(4.d0*pi*vevsc)

      yuin(3,1) = (MT*VCKM(3,1))/(4.d0*pi*vevsc)
      yuin(3,2) = (MT*VCKM(3,2))/(4.d0*pi*vevsc)
      yuin(3,3) = (MT*VCKM(3,3))/(4.d0*pi*vevsc)

      alph3in = alphsMZdrbar/(4.d0 * pi)   !alph3tz
      alph2in = alphaemMZdrbar/(4.d0 * pi) !alph2tz
      alph1in = alph1tz/(4.d0 * pi)

c$$$      print*,"yuin(3,3) = ", yuin(3,3)*4.d0*pi/(dsin(beta)),mt/(vev2sc),
c$$$     $     yuin(3,3),mt
c$$$      print*,"ydin(3,3) = ", ydin(3,3)*4.d0*pi/(dcos(beta))
c$$$      print*,"yein(3,3) = ", yein(3,3)*4.d0*pi/(dcos(beta))
c$$$
c$$$      print*,"vev = ", vevsc
c$$$
c$$$      print*,"MWin = ", MZ,  MZ*dcos(dasin(dsqrt(sinsqtw))), MZpole, 
c$$$     $     MWpole               !dsqrt(alph2in*16.d0*pi*pi*vevin**2.d0/4.d0)

!--------------------------
!     Inputs at high energy
!--------------------------
            
      Rparametrization: if(case.eq.'Rpr')then

         thl12 = dASin(dsqrt(0.320d0))
         thl23 = dASin(dsqrt(0.49d0))
         thl13 = dASin(dsqrt(0.026d0))
         
         UPMNS(1,1) = dCos(thl12)*dCos(thl13)
         UPMNS(1,2) = dSin(thl12)*dCos(thl13)
         UPMNS(1,3) = dSin(thl13) 

         UPMNS(2,1) = - dSin(thl12)*dCos(thl23) - 
     $        dCos(thl12)*dSin(thl23)*dSin(thl13)
         UPMNS(2,2) = dCos(thl12)*dCos(thl23) -
     $        dSin(thl12)*dSin(thl23)*dSin(thl13)
         UPMNS(2,3) = dCos(thl13)*dSin(thl23)

         UPMNS(3,1) = dSin(thl12)*dSin(thl23) - 
     $        dCos(thl12)*dCos(thl23)*dSin(thl13)
         UPMNS(3,2) = - dCos(thl12)*dSin(thl23) - 
     $        dSin(thl12)*dCos(thl23)*dSin(thl13) 
         UPMNS(3,3) = dCos(thl23)*dCos(thl13)


         call dag(UPMNS,UPMNST)
         call mat3prod(DM,R,Dk/vev2sc,drd)
         call matmult(drd/(4.d0*pi),UPMNST,Ynui)


      endif Rparametrization

c$$$      print*,"----recent"
c$$$      recmns: do i = 1, 3
c$$$!     223  format(1x,3(2x,1es11.4))
c$$$      write(*,223) (UPMNSR(i,j),j = 1, 3)
c$$$      enddo recmns  
c$$$
c$$$      print*,"------ stan"
c$$$      stnmns: do i = 1, 3
c$$$!     223  format(1x,3(2x,1es11.4))
c$$$      write(*,223) (UPMNS(i,j),j = 1, 3)
c$$$      enddo stnmns  

!------------------

      MX = 5.d0*(10.d0**19.d0)


      if(model.eq.'mSUG'.or.model.eq.'NUHM'.or.
     $     model.eq.'NUGM'.or.model.eq.'CNUM')then
         

         M1X = Mg1
         M2X = Mg2
         M3X = Mg3

         mh10 = sign(1.d0,m10)*m10*m10
         mh20 = sign(1.d0,m20)*m20*m20


!--   mq                               !<------------- NOTE: -ve inputs to higgs m10,m20 implies -m10^2 
         mQ0(1,1) = sign(1.d0,mq11)*mq11**2.d0
         mQ0(1,2) = sign(1.d0,mq12)*mq12**2.d0
         mQ0(1,3) = sign(1.d0,mq13)*mq13**2.d0
         mQ0(2,1) = sign(1.d0,mq21)*mq21**2.d0
         mQ0(2,2) = sign(1.d0,mq22)*mq22**2.d0
         mQ0(2,3) = sign(1.d0,mq23)*mq23**2.d0
         mQ0(3,1) = sign(1.d0,mq31)*mq31**2.d0
         mQ0(3,2) = sign(1.d0,mq32)*mq32**2.d0
         mQ0(3,3) = sign(1.d0,mq33)*mq33**2.d0

!---  mu

         mU0(1,1) = sign(1.d0,mu11)*mu11**2.d0
         mU0(1,2) = sign(1.d0,mu12)*mu12**2.d0
         mU0(1,3) = sign(1.d0,mu13)*mu13**2.d0
         mU0(2,1) = sign(1.d0,mu21)*mu21**2.d0
         mU0(2,2) = sign(1.d0,mu22)*mu22**2.d0
         mU0(2,3) = sign(1.d0,mu23)*mu23**2.d0
         mU0(3,1) = sign(1.d0,mu31)*mu31**2.d0
         mU0(3,2) = sign(1.d0,mu32)*mu32**2.d0
         mU0(3,3) = sign(1.d0,mu33)*mu33**2.d0


!---  md

         mD0(1,1) = sign(1.d0,md11)*md11**2.d0
         mD0(1,2) = sign(1.d0,md12)*md12**2.d0
         mD0(1,3) = sign(1.d0,md13)*md13**2.d0
         mD0(2,1) = sign(1.d0,md21)*md21**2.d0
         mD0(2,2) = sign(1.d0,md22)*md22**2.d0
         mD0(2,3) = sign(1.d0,md23)*md23**2.d0
         mD0(3,1) = sign(1.d0,md31)*md31**2.d0
         mD0(3,2) = sign(1.d0,md32)*md32**2.d0
         mD0(3,3) = sign(1.d0,md33)*md33**2.d0

!---  ml
         mL0(1,1) = sign(1.d0,ml11)*ml11**2.d0
         mL0(1,2) = sign(1.d0,ml12)*ml12**2.d0
         mL0(1,3) = sign(1.d0,ml13)*ml13**2.d0
         mL0(2,1) = sign(1.d0,ml21)*ml21**2.d0
         mL0(2,2) = sign(1.d0,ml22)*ml22**2.d0
         mL0(2,3) = sign(1.d0,ml23)*ml23**2.d0
         mL0(3,1) = sign(1.d0,ml31)*ml31**2.d0
         mL0(3,2) = sign(1.d0,ml32)*ml32**2.d0
         mL0(3,3) = sign(1.d0,ml33)*ml33**2.d0

!---  me
         mE0(1,1) = sign(1.d0,me11)*me11**2.d0
         mE0(1,2) = sign(1.d0,me12)*me12**2.d0
         mE0(1,3) = sign(1.d0,me13)*me13**2.d0
         mE0(2,1) = sign(1.d0,me21)*me21**2.d0
         mE0(2,2) = sign(1.d0,me22)*me22**2.d0
         mE0(2,3) = sign(1.d0,me23)*me23**2.d0  
         mE0(3,1) = sign(1.d0,me31)*me31**2.d0
         mE0(3,2) = sign(1.d0,me32)*me32**2.d0
         mE0(3,3) = sign(1.d0,me33)*me33**2.d0


!--   mnu

         mNU0(1,1) = sign(1.d0,mnu11)*mnu11**2.d0
         mNU0(1,2) = sign(1.d0,mnu12)*mnu12**2.d0
         mNU0(1,3) = sign(1.d0,mnu13)*mnu13**2.d0
         mNU0(2,1) = sign(1.d0,mnu21)*mnu21**2.d0
         mNU0(2,2) = sign(1.d0,mnu22)*mnu22**2.d0
         mNU0(2,3) = sign(1.d0,mnu23)*mnu23**2.d0
         mNU0(3,1) = sign(1.d0,mnu31)*mnu31**2.d0
         mNU0(3,2) = sign(1.d0,mnu32)*mnu32**2.d0
         mNU0(3,3) = sign(1.d0,mnu33)*mnu33**2.d0
         
      endif

      

      mur   = 0.d0 
      bmur  = 0.d0 
      

C     ====================================
C     Calling the RGE running routine!!
C     ====================================
      
C     initialize to 0 the checking parameter
      
      check = 0

      msusyold = 1000.d0        ! initial guess for msusy


      itcount = 1
      runum = 1
      MTc_mz = 0.d0
      mBc_mz = 0.d0
      mtauc_mz = 0.d0

      delalphas = 0.d0
      delalphem = 0.d0
      sinsqtheff = sinsqtw
      flags = ' AOK'
      exitcalc = 'F'

      call rgeit(MW,MZ,MT,MTc_mz,mB,mBc_mz,mTau,
     $     mtauc_mz, msusyold,vevsc,vevin,vev1in,vev2in,yuin,
     $     ydin,yein,alphaDR,alph1in,alph2in,delalphem,alph3in,
     $     delalphas,mur,bmur, murge,bmurge,prnstat,check,
     $     newtbeta, MTatMZ, msusynew,mursq,try,
     $     M3t,flags,runum,itcount,stopratu,
     $     stopratd,sinsqtheff,exitcalc)

 19   format(1es11.4,2x,1es11.4,2x,1es11.4,2x,1es11.4)   

!-----------------------------------------------

      brbd = 0.d0
      brbs = 0.d0

      if(isnan(AURG(3,3)).or.isnan(yuRG(3,3)))then

         AURG(3,3) = 0.d0
         yuRG(3,3) = 0.d0

      endif

      if(isnan(mh0sqcor))then
         
         flags = 'TACMh'

         mh0sqcor = 0.d0
         
         flagt = flags

      endif

      if(isnan(mHu0sqcor))then
         
         flags = 'TACMH'

         mHu0sqcor = 0.d0
         
         flagt = flags

      endif

      if(isnan(mA0sqcor))then
         
         flags = 'TACMA'

         mA0sqcor = 0.d0
         
         flagt = flags

      endif

      if(isnan(mHpmsqcor))then
         
         flags = 'TACMHPM'

         mHpmsqcor = 0.d0
         
         flagt = flags

      endif

      if(isnan(mh0sq))then
         
         flags = 'TACMh'

         mh0sq = 0.d0
         
         flagt = flags

      endif

      if(isnan(mHu0sq))then
         
         flags = 'TACMH'

         mHu0sq = 0.d0
         
         flagt = flags

      endif

      if(isnan(mA0sq))then
         
         flags = 'TACMA'

         mA0sq = 0.d0
         
         flagt = flags

      endif

      if(isnan(mHpmsq))then
         
         flags = 'TACMHPM'

         mHpmsq = 0.d0
         
         flagt = flags

      endif


      if(flags.eq.'variable underflow ')then
         flags = 'VARUNDER'
         
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

!         write(93,19) m0,m12,a0,tanbeta
      endif

!      print*,tacsup,tasdn,tacslp,tacsnu,tachiggs
      
      flagt= 'AOK'
      
      if(flags.eq.' AOK')then
         num = 1
         if((Suegg(1).lt.0.d0).OR.(Suegg(2).lt.0.d0).OR.
     $        (Suegg(3).lt.0.d0).OR.(Suegg(4).lt.0.d0).OR.
     $        (Suegg(5).lt.0.d0).OR.(Suegg(6).lt.0.d0))then
            
            Suegg(1) = dabs(Suegg(1))
            Suegg(2) = dabs(Suegg(2))
            Suegg(3) = dabs(Suegg(3))
            Suegg(4) = dabs(Suegg(4))
            Suegg(5) = dabs(Suegg(5))
            Suegg(6) = dabs(Suegg(6))
            
            flags = "TACSUP"
            flagt = flags
         else
            
            continue
         endif
!-----
!     Down sector

         if((SDegg(1).lt.0.d0).OR.(SDegg(2).lt.0.d0).OR.
     $        (SDegg(3).lt.0.d0).OR.(SDegg(4).lt.0.d0).OR.
     $        (SDegg(5).lt.0.d0).OR.(SDegg(6).lt.0.d0))then

            SDegg(1) = dabs(SDegg(1))
            SDegg(2) = dabs(SDegg(2))
            SDegg(3) = dabs(SDegg(3))
            SDegg(4) = dabs(SDegg(4))
            SDegg(5) = dabs(SDegg(5))
            SDegg(6) = dabs(SDegg(6))
            
            flags = 'TACSDN'
            flagt = flags

         else

            continue
         endif

!-------
!     lepton sector

         if((SLegg(1).lt.0.d0).OR.(SLegg(2).lt.0.d0).OR.
     $        (SLegg(3).lt.0.d0).OR.(SLegg(4).lt.0.d0).OR.
     $        (SLegg(5).lt.0.d0).OR.(SLegg(6).lt.0.d0))then


            SLegg(1) = dabs(SLegg(1))
            SLegg(2) = dabs(SLegg(2))
            SLegg(3) = dabs(SLegg(3))
            SLegg(4) = dabs(SLegg(4))
            SLegg(5) = dabs(SLegg(5))
            SLegg(6) = dabs(SLegg(6))
            
            flags = 'TACSLP'
            flagt = flags

         else
            
            continue
         endif
!-------------------------
!     sneutrinos

         if((SNegg(1).lt.0.d0).OR.(SNegg(2).lt.0.d0).OR.
     $        (SNegg(3).lt.0.d0))then

            SNegg(1) = dabs(SNegg(1))
            SNegg(2) = dabs(SNegg(2))
            SNegg(3) = dabs(SNegg(3))

            flags = 'TACSNU'
            flagt = flags

         else

            continue
         endif

         if((mhpmsq.lt.0.d0))then

            mhpmsq = dabs(mhpmsq)

            flags = 'TACHPM'
            flagt = flags

         endif


         if((mA0sq.lt.0.d0))then

            mA0sq = dabs(mA0sq)

            flags = 'TACMA'

            flagt = flags
         endif

         if((mh0sq.lt.0.d0))then

            mh0sq = dabs(mh0sq)

            flags = 'TACMh'

            flagt = flags
         endif

         if((mHu0sq.lt.0.d0))then

            mHu0sq = dabs(mHu0sq)

            flags = 'TACMH'

            flagt = flags

         endif

         if(isnan(mh0sq))then
            
            flags = 'TACMh'

            mh0sq = 0.d0
         
            flagt = flags

         endif

         if(isnan(mHu0sq))then
            
            flags = 'TACMH'

            mHu0sq = 0.d0
         
            flagt = flags

         endif

         if(isnan(mA0sq))then
            
            flags = 'TACMA'

            mA0sq = 0.d0
         
            flagt = flags

         endif

         if(isnan(mHpmsq))then
            
            flags = 'TACMHPM'

            mHpmsq = 0.d0
         
            flagt = flags

         endif

!------------------------
C     higgs FLAG-

         if((mh0sq.lt.(114.5d0)**2.d0).and.mh0sq.gt.0.d0)then
            flags = 'LEPH'

            if(prnstat.eq.1)then
               write(10,*) 'Invalid point ', flags
            else if(prnstat.eq.0)then
               write(*,*) 'Invalid point ', flags
            endif
         else
            continue
         ENDIF

!--------

C     Chargino flag
         
         if((MIN(dabs(Ceg(1)), dabs(Ceg(2)))).le.103.5d0)then       
            num = 1
            flags = 'LEPC'

!            write(97,19) m0,m12,a0,tanbeta

            if(prnstat.eq.1)then
               write(10,*) 'Invalid point ', flags
            else if(prnstat.eq.0)then
               write(*,*) 'Invalid point ', flags
            endif
         else
            continue
         ENDIF

!--------------------------
C     finding lsp

         minspec(1) = SUegg(1)  !stops
         minspec(2) = SUegg(2)
         minspec(3) = SUegg(3)  !scharms
         minspec(4) = SUegg(4) 
         minspec(5) = SUegg(5)  !sups
         minspec(6) = SUegg(6) 
         minspec(7) = SDegg(1)  !sbottoms
         minspec(8) = SDegg(2) 
         minspec(9) = SDegg(3)  !sstranges
         minspec(10) = SDegg(4)
         minspec(11) = SDegg(5) ! sdowns
         minspec(12) = SDegg(6) 
         minspec(13) = SLegg(1) !staus
         minspec(14) = SLegg(2)
         minspec(15) = SLegg(3) !smu
         minspec(16) = SLegg(4)  
         minspec(17) = SLegg(5) !sel
         minspec(18) = SLegg(6)
         minspec(19) = SNegg(1) !sneutrinos
         minspec(20) = SNegg(2) 
         minspec(21) = SNegg(3)
         minspec(22) = Neg(1)*Neg(1) ! neutralino
         minspec(23) = Neg(2)*Neg(2)
         minspec(24) = Neg(3)*Neg(3)
         minspec(25) = Neg(4)*Neg(4)
         minspec(26) = Ceg(1)*Ceg(1) ! chargino
         minspec(27) = Ceg(2)*Ceg(2)
!     minspec(28) = newmh0sq       ! lightest higgs
         minspec(28) = mhu0sq   ! cp even higgs
         minspec(29) = mhpmsq   ! charged higgs
         minspec(30) = mA0sq    ! cp odd higgs

         if(model.eq.'GMSB')then
            minspec(31) = mgrav
            lsppos = MINLOC(minspec(1:31), DIM = 1) 
            lsp =  MINVAL(minspec(1:31))
            nlsp =  MINVAL(minspec(1:31), DIM =1, 
     $           MASK = minspec(1:31).gt.lsp)
            nlsppos = MINLOC(minspec(1:31), DIM = 1, 
     $           MASK = minspec(1:31).gt.lsp)
         else
            lsppos = MINLOC(minspec(1:30), DIM = 1) 
            lsp =  MINVAL(minspec(1:30))
            nlsp =  MINVAL(minspec(1:30), DIM =1, 
     $           MASK = minspec(1:30).gt.lsp)
            nlsppos = MINLOC(minspec(1:30), DIM = 1,
     $           MASK = minspec(1:30).gt.lsp)
         endif

         if(lsppos.eq.13.OR.lsppos.eq.14)then
            
            flags = "LSPSTAU"
            if(prnstat.eq.1)then
               write(10,*) 'Invalid point ', flags
            else if(prnstat.eq.0)then
               write(*,*) 'Invalid point ', flags
            endif
            
!            write(99,19) m0,m12,a0,tanbeta
            
         ENDIF
         

         
         
      else
         continue
      endif
      
      
!      if(flagt.ne.'AOK')then
!         write(94,19) m0,m12,a0,tanbeta
!      endif
      
      
      if(flags.eq.' AOK'.or.flags.eq.'LEPC'.or.flags.eq.'LEPH'.or.
     $     flags.eq."LSPSTAU")then

         errge = 0
         

         if(flagt.eq.'AOK')then
            errget = 0
         else
            errget = 1
            flags = flagt
         endif

      else

         errge = 1
         
         flagt = flags(1:2)//flags(3:4)//flags(5:6)//flags(7:8)//
     $        flags(9:9)

         print*,"flagt = ", flagt

         flags = flagt

         if(prnstat.eq.1)then
            write(10,*) 'Invalid point ', flags
            write(*,*) 'Invalid point ', flags
         elseif(prnstat.eq.0)then
            write(*,*) 'Invalid point ', flags
         endif
         
      endif
      
      
      mueconverrate = 0.d0
      brmueconver = 0.d0
      

!=======

c$$$          if(flagt.ne.'AOK')then
c$$$             write(*,*) m0,m12,a0, flagt
c$$$          else
c$$$             write(*,*) m0,m12,a0, flags
c$$$          endif

      specok: if(errge.eq.0.and.errget.eq.0)then

!---------initialization

         ftmz = 0.d0
         ftmt = 0.d0
         drho = 0.d0

         verbspec: if(prnstat.eq.1)then
 50      format(2x,A22,2x,1pe11.4)

         write(10,*) ' up-type yukawa at high energy:'
         write(*,*) ' up-type yukawa at high energy:'
         opu23: do i = 1, 3
 224     format(1x,3(2x,1pe12.4))
         write(10,224) (yugut(i,j),j = 1, 3)
         write(*,224) (yugut(i,j),j = 1, 3)
         enddo opu23  
         write(10,*) ''
         write(*,*) ''

         write(10,*) ' down-type yukawa at high energy:'
         write(*,*) ' down-type yukawa at high energy:'
         opu24: do i = 1, 3
         write(10,224) (ydgut(i,j),j = 1, 3)
         write(*,224) (ydgut(i,j),j = 1, 3)
         enddo opu24  
         write(10,*) ''
         write(*,*) ''

         write(10,*) ' lepton-type yukawa at high energy:'
         write(*,*) ' lepton-type yukawa at high energy:'
         opu25: do i = 1, 3
         write(10,224) (yegut(i,j),j = 1, 3)
         write(*,224) (yegut(i,j),j = 1, 3)
         enddo opu25  
         write(10,*) ''
         write(*,*) ''


         write(10,*) ' nuetrino yukawa at high energy:'
         write(*,*) ' neutrino yukawa at high energy:'
         opu25a: do i = 1, 3
         write(10,224) (ynugut(i,j),j = 1, 3)
         write(*,224) (ynugut(i,j),j = 1, 3)
         enddo opu25a  
         write(10,*) ''
         write(*,*) ''


         write(10,*) ' up-type yukawa at msusy :'
         write(*,*) ' up-type yukawa at msusy :'
         opu26: do i = 1, 3
         write(10,224) (yuRG(i,j),j = 1, 3)
         write(*,224) (yuRG(i,j),j = 1, 3)
         enddo opu26  
         write(10,*) ''
         write(*,*) ''


         write(10,*) ' down-type yukawa at msusy :'
         write(*,*) ' down-type yukawa at msusy :'
         opu27: do i = 1, 3
         write(10,224) (ydRG(i,j),j = 1, 3)
         write(*,224) (ydRG(i,j),j = 1, 3)
         enddo opu27  
         write(10,*) ''
         write(*,*) ''


         write(10,*) ' lepton-type yukawa at msusy :'
         write(*,*) ' lepton-type yukawa at msusy :'
         opu28: do i = 1, 3
         write(10,224) (yeRG(i,j),j = 1, 3)
         write(*,224) (yeRG(i,j),j = 1, 3)
         enddo opu28  
         write(10,*) ''
         write(*,*) ''
         
         Write(10,50) 'Spectrum at msusy, q =  ', msusyold
         Write(*,50) 'Spectrum at msusy, q =  ', msusyold
         Write(10,50) 'alpha_1  =', dsqrt(alph1*(16.d0*pi*pi)*0.6)
         Write(*,50) 'alpha_1  =', dsqrt(alph1*(16.d0*pi*pi)*0.6)
         Write(10,50) 'alpha_2  =', dsqrt(alph2*(16.d0*pi*pi))
         Write(*,50) 'alpha_2  =', dsqrt(alph2*(16.d0*pi*pi))
         Write(10,50) 'alpha_3  =', dsqrt(alph3*(16.d0*pi*pi))
         Write(*,50) 'alpha_3  =', dsqrt(alph3*(16.d0*pi*pi))
         write(10,50)'vev1    =',vev1
         write(*,50)'vev1    =',vev1
         write(10,50)'vev2    =',vev2           
         write(*,50)'vev2    =',vev2           
         write(10,50) 'newtbeta    =', (vev2/vev1)
         write(*,50) 'newtbeta    =', (vev2/vev1)
         write(10,50) '\mu  =', sgnmu*mur
         write(*,50) '\mu  =', sgnmu*mur
         write(10,50) '~gluino  =', M3t
         write(*,50) '~gluino  =', M3t
         write(10,50) '~stop_1 =', dsqrt(SUegg(1))
         write(*,50) '~stop_1 =', dsqrt(SUegg(1))
         write(10,50) '~stop_2 =', dsqrt(SUegg(2))
         write(*,50) '~stop_2 =', dsqrt(SUegg(2))
         write(10,50) '~scharm_R =', dsqrt(SUegg(3))
         write(*,50) '~scharm_R =', dsqrt(SUegg(3))
         write(10,50) '~scharm_L =', dsqrt(SUegg(4))
         write(*,50) '~scharm_L =', dsqrt(SUegg(4))
         write(10,50) '~sup_R =', dsqrt(SUegg(5))
         write(*,50) '~sup_R =', dsqrt(SUegg(5))
         write(10,50) '~sup_L =', dsqrt(SUegg(6))
         write(*,50) '~sup_L =', dsqrt(SUegg(6))
         write(10,50) '~sbottom_1  =', dsqrt(SDegg(1))
         write(*,50) '~sbottom_1  =', dsqrt(SDegg(1))
         write(10,50) '~sbottom_2  =', dsqrt(SDegg(2))
         write(*,50) '~sbottom_2  =', dsqrt(SDegg(2))
         write(10,50) '~sstrange_R =', dsqrt(SDegg(3))
         write(*,50) '~sstrange_R =', dsqrt(SDegg(3))
         write(10,50) '~sstrange_L =', dsqrt(SDegg(4))
         write(*,50) '~sstrange_L =', dsqrt(SDegg(4))
         write(10,50) '~sdown_R  =', dsqrt(SDegg(5))
         write(*,50) '~sdown_R  =', dsqrt(SDegg(5))
         write(10,50) '~sdown_L  =', dsqrt(SDegg(6))
         write(*,50) '~sdown_L  =', dsqrt(SDegg(6))
         write(10,50) '~stau_1 =', dsqrt(SLegg(1))
         write(*,50) '~stau_1 =', dsqrt(SLegg(1))
         write(10,50) '~stau_2 =', dsqrt(SLegg(2))
         write(*,50) '~stau_2 =', dsqrt(SLegg(2))
         write(10,50) '~smu_R  =', dsqrt(SLegg(3))
         write(*,50) '~smu_R  =', dsqrt(SLegg(3))
         write(10,50) '~smu_L  =', dsqrt(SLegg(4))
         write(*,50) '~smu_L  =', dsqrt(SLegg(4))
         write(10,50) '~sel_R   =', dsqrt(SLegg(5))
         write(*,50) '~sel_R   =', dsqrt(SLegg(5))
         write(10,50) '~sel_L   =', dsqrt(SLegg(6))
         write(*,50) '~sel_L   =', dsqrt(SLegg(6))
         write(10,50) '~tausnu  =', dsqrt(SNegg(1))
         write(*,50) '~tausnu  =', dsqrt(SNegg(1))
         write(10,50) '~musnu =', dsqrt(SNegg(2))
         write(*,50) '~musnu =', dsqrt(SNegg(2))
         write(10,50) '~elsnu =', dsqrt(SNegg(3))           
         write(*,50) '~elsnu =', dsqrt(SNegg(3))           
         write(10,*) ''
         write(*,*) ''
         
         write(10,*) 'Higgs Spectrum'
         write(*,*) 'Higgs Spectrum'
         write(10,50) 'mA0   =',dsqrt(mA0sq)
         write(*,50) 'mA0   =',dsqrt(mA0sq)
         write(10,50) 'mh_charged =',dsqrt(mhpmsq)
         write(*,50) 'mh_charged =',dsqrt(mhpmsq)
         write(10,50) 'mh0   =',dsqrt(dabs(mh0sq))
         write(*,50) 'mh0   =',dsqrt(dabs(mh0sq))
         write(10,50) 'mH   =',dsqrt(mHu0sq)           
         write(*,50) 'mH   =',dsqrt(mHu0sq)           
         write(10,50) '\alpha_h   =',alphahiggs
         write(*,50) '\alpha_h   =',alphahiggs
         write(10,*) ''
         write(*,*) ''

         Write(10,*) 'Neutralino spectrum'
         Write(*,*) 'Neutralino spectrum'
         write(10,50) 'N1   =', Neg(1)
         write(*,50) 'N1   =', Neg(1)
         write(10,50) 'N2   =', Neg(2)
         write(*,50) 'N2   =', Neg(2)
         write(10,50) 'N3   =', Neg(3)
         write(*,50) 'N3   =', Neg(3)
         write(10,50) 'N4   =', Neg(4)
         write(*,50) 'N4   =', Neg(4)
         write(10,*) ''
         write(*,*) ''

         Write(10,*) 'Chargino spectrum'
         Write(*,*) 'Chargino spectrum'
         write(10,50) 'C1   =', Ceg(1)
         write(*,50) 'C1   =', Ceg(1)
         write(10,50) 'C2   =', Ceg(2)
         write(*,50) 'C2   =', Ceg(2)
         write(10,*) ''
         write(*,*) ''

      else if(prnstat.eq.0)then
         continue
      endif verbspec

!------------------------------------------------------

       call bsg(tanbeta,MT,mB,Ceg,SUegvd,USUvd,OCL,OCR,
     $        dsqrt(mhpmsq),Bbsg)

       if(isnan(Bbsg)) Bbsg = 0.d0


      call muegamma(newtbeta,mTau,alph2,ONm,Neg,Ceg,OCL,OCR,
     $      SLegvd,USLvd,SNegvd,USNvd,SUegvd,USUvd,SDegvd,USDvd,
     $      megrate,Bmeg,tmugrate,Btmug,tegrate,Bteg,mu3erate,
     $      Brmu3E,tau3murate,brtau3mu,tau3erate,brtau3e,
     $      mueconverrate,brmueconver,gm2)

C     ------------------------------------------
C     Scaling the Branching Ratios for plotting 
C     ------------------------------------------

      
      scbmeg = (10.d0**(11.d0))*bmeg 
      scbmu3e = (10.d0**(12.d0))*brmu3e 
      scbmueconver = (10.d0**(12.d0))*brmueconver 
      
      
      scbtmug = (10.d0**(8.d0))*btmug
      scbtau3mu = (10.d0**(7.d0))*brtau3mu  
      
      scbteg = (10.d0**(7.d0))*bteg
      scbtau3e = (10.d0**(7.d0))*brtau3e
      
      scgminus2 = (10.d0**(10.d0))*gm2

      
 99   format(1pe11.4,2x,1pe11.4,2x,1pe11.4,2x,1pe11.4,2x,
     $     1pe11.4,2x,1pe11.4,2x,1pe11.4,2x,1pe11.4)
      
!      write(95,99) m0,m12,a0,bmeg,bteg,btmug,gm2,Bbsg
      

c      print*,"----------(mu -> e, gamma ends here)----------------"     

C       Fine-tuning stuff-----------------------------------
         

      call finetune(MZ,tanbeta,murgemz,mh1mzz,mh2mzz,delta1z,delta2z,
     $        ftmz,ftmt)

C      delrho


      call deltarhoew(MW,MZ,piwwT0,pizzT0,piwwT,pizzT,drho)


C--------------------------------------------------------------  
C     if prnstat writes

 51   format(A30,2x,1pe11.4,1x,A10)
 52   format(A30,2x,1pe11.4,1x)

      verb10: if(prnstat.eq.1) then
         Write(10,*)' Low Energy Observables '
         Write(*,*)' Low Energy Observables '
         Write(10,52)'Fine tuning, Cmz^2mu^2   = ', ftmz
         Write(*,52)'Fine tuning, Cmz^2mu^2   = ', ftmz
         Write(10,52)'Fine tuning, Cmt^2mu^2   = ', ftmt
         Write(*,52)'Fine tuning, Cmt^2mu^2   = ', ftmt
         Write(10,52)'Br(B  => s,gamma)   = ',Bbsg
         Write(*,52)'Br(B  => s,gamma)   = ',Bbsg
         Write(10,51)'Br(mu  => e,gamma)   = ',scbmeg,'X 10^(-11)'
         Write(*,51)'Br(mu  => e,gamma)   = ',scbmeg,'X 10^(-11)'
         Write(10,51)'Br(tau => mu,gamma)  = ',scbtmug,'X 10^(-08)'
         Write(*,51)'Br(tau => mu,gamma)  = ',scbtmug,'X 10^(-08)'
         Write(10,51)'Br(tau => e,gamma)   = ',scbteg,'X 10^(-07)'
         Write(*,51)'Br(tau => e,gamma)   = ',scbteg,'X 10^(-07)'
         Write(10,51)'Br(tau => e,e,e)   = ',scbtau3e,'X 10^(-07)'
         Write(*,51)'Br(tau => e,e,e)   = ',scbtau3e,'X 10^(-07)'
         Write(10,51)'Br(tau => mu,mu,mu)   = ',scbtau3mu,'X 10^(-07)'
         Write(*,51)'Br(tau => mu,mu,mu)   = ',scbtau3mu,'X 10^(-07)'
         Write(10,51)'Br(mu => e,e,e)   = ',scbmu3e,'X 10^(-12)'
         Write(*,51)'Br(mu => e,e,e)   = ',scbmu3e,'X 10^(-12)'
         Write(10,51)'Br(mu => e in Ti)   = ',scbmueconver,'X 10^(-12)'
         Write(*,51)'Br(mu => e in Ti)   = ',scbmueconver,'X 10^(-12)'
         write(10,51)'      (g_mu - 2)      = ',scgminus2,'X 10^(-10)'
         write(*,51)'       (g_mu - 2)      = ',scgminus2,'X 10^(-10)'
         Write(10,51)' -----------------------------------------------'
         Write(*,51)' -----------------------------------------------'
      endif verb10



!--------writing out in slha format


      call printslha(errge,flags,msusyold,mur,mAm3,
     $     MT,mb,mtau,MW,MZ,mh1mz,mh2mz,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,ONm,OCLTm,OCRm,alph3,alph2,alph1,
     $     M1tz,M2tz,M3tz,M3t,thetat,thetab,thetatau,
     $     mh0sq,mHu0sq,mA0sq,mhpmsq,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,mSQRG,mSURG,mSDRG,mSLRG,mSERG,
     $     USL,ftmz,ftmt,MWsqpole_MZ)
      
      else specok

      optsp: if(model.eq.'mSUG')then

         WRITE(999,82) tanbeta,m0, m12,a0,sgnmu,dsqrt(mh0sq),
     $        gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags
         
      else if(model.eq.'GMSB')then
         
         WRITE(999,83) tanbeta,gmsbsusyb, gmsbmess,nhat,sgnmu,
     $        dsqrt(mh0sq),gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,
     $        brtau3e,flags
         
      else if(model.eq.'NUHM')then
         
         WRITE(999,84) tanbeta,m0, m12,a0,m10,m20,mh0sq,
     $        gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags

      endif optsp
      
      if(model.eq.'mSUG') then 
C         write(iunit,97) tanbeta,m0,m12,a0,flags
         write(iunit,87) "flag for the point is",flags
         if(prnstat.eq.1)then
C            write(10,97) tanbeta,m0,m12,a0,flags
         write(10,87) "flag for the point is",flags
         endif
      else if(model.eq.'NUHM') then
C         write(iunit,972) tanbeta,m0,m12,a0,m10,m20,flags
         write(iunit,87) "flag for the point is",flags
         if(prnstat.eq.1)then
C            write(10,972) tanbeta,m0,m12,a0,m10,m20,flags
         write(10,87) "flag for the point is",flags
         endif
      else if(model.eq.'GMSB') then
C         write(iunit,97) tanbeta,gmsbsusyb,gmsbmess,nhat,flags
         write(iunit,87) "flag for the point is",flags
         if(prnstat.eq.1)then
C            write(10,97) tanbeta,gmsbsusyb,gmsbmess,nhat,flags 
         write(10,87) "flag for the point is",flags
         endif
      else
C         write(iunit,97) tanbeta,m0,m12,a0,flags
         write(iunit,87) "flag for the point is",flags
         if(prnstat.eq.1)then
C            write(10,97) tanbeta,m0,m12,a0,flags
         write(10,87) "flag for the point is",flags
         endif
      endif
      
      RETURN
      
      endif specok
      
 97   format(1x,1pe11.4,2x,1pe11.4,2x,1pe11.4,2x,1pe11.4,2x,A)   
 87   format(1x,A,1x,A)   
      
!----console output

!----console output

      opt: if(model.eq.'mSUG')then

 82   FORMAT(1pE11.4,2x,1pE11.4, 2x,1pE11.4, 2x,1pE11.4,2x,
     $        1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,
     $        1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,
     $        1pE11.4,2x,1pE11.4,2x,A)


      WRITE(999,82) tanbeta,m0, m12,a0,sgnmu,dsqrt(mh0sq),
     $     gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags

      if(prnstat.eq.1)then
      WRITE(*,*) "-------------------------------------------- " 
      WRITE(*,*) "Copy of the output in `suseflav.out' " 
      WRITE(*,*) "-------------------------------------------- " 
      else if(prnstat.eq.0)then
         continue
      endif
         
      WRITE(*,'(1x,2A,/,2A,/,2A)') "Observables written in",
     $     " tmp/output.txt in format: tanbeta m0 m12 a0",
     $     " sgnmu lightHiggs g-2 Brbsgamma Brmuegamma",
     $     " Brtaumugamma Brtauegamma", " Brmu3e Brtau3mu Brtau3e",
     $     " flags"

      WRITE(*,*) "--------------------------------------------------" 

C      write(iunit,97) tanbeta,m0,m12,a0,flags
         write(iunit,87) "flag for the point is",flags
      if(prnstat.eq.1)then
C         write(10,97) tanbeta,m0,m12,a0,flags
         write(10,87) "flag for the point is",flags
      endif

      else if(model.eq.'GMSB')then

 83      FORMAT(1pE11.4,2x,1pE11.4, 2x,1pE11.4, 2x,1pE11.4,2x,
     $        1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,
     $        1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,
     $        1pE11.4,2x,1pE11.4,2x,A)

         WRITE(999,83) tanbeta,gmsbsusyb, gmsbmess,nhat,sgnmu,
     $        dsqrt(mh0sq),gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,
     $        brtau3e,flags
         

         if(prnstat.eq.1)then
            WRITE(*,*) "-------------------------------------------- " 
            WRITE(*,*) "Copy of the output in `suseflav.out' " 
            WRITE(*,*) "-------------------------------------------- " 
         else if(prnstat.eq.0)then
            continue
         endif

         WRITE(*,'(1x,2A,/,2A,/,2A)') "Observables written in",
     $        " tmp/output.txt in format: tanbeta Lambda Mmess",
     $        " nmess sgnmu lightHiggs g-2 Brbsgamma Brmuegamma",
     $        " Brtaumugamma Brtauegamma", " Brmu3e Brtau3mu Brtau3e",
     $        " flags"

      WRITE(*,*) "--------------------------------------------------" 
         
C         write(iunit,97) tanbeta,gmsbsusyb,gmsbmess,nhat,flags
         write(iunit,87) "flag for the point is",flags
         if(prnstat.eq.1)then
C            write(10,97) tanbeta,gmsbsusyb,gmsbmess,nhat,flags 
         write(10,87) "flag for the point is",flags
         endif

      else if(model.eq.'NUHM')then
         
 84      format(15(1pE15.5,2x),2x,G8.0)
         
         WRITE(999,84) tanbeta,m0, m12,a0,m10,m20,mh0sq,
     $        gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags

         
 972     FORMAT(1x,1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,
     $        1pE11.4,2x,1pE11.4,2x,A)
         
         
         if(prnstat.eq.1)then
            WRITE(*,*) "-------------------------------------------- " 
            WRITE(*,*) "Copy of the output in `suseflav.out' " 
            WRITE(*,*) "-------------------------------------------- " 
         else if(prnstat.eq.0)then
            continue
         endif
         
      WRITE(*,'(1x,2A,/,2A,/,2A)') "Observables written in",
     $     " tmp/output.txt in format: tanbeta m0 m12 a0",
     $     " sgnmu mh10 mh20 lightHiggs g-2 Brbsgamma Brmuegamma",
     $     " Brtaumugamma", " Brtauegamma Brmu3e Brtau3mu Brtau3e",
     $     " flags"
      
      WRITE(*,*) "--------------------------------------------------" 
      
C     write(iunit,972) tanbeta,m0,m12,a0,m10,m20,flags
      write(iunit,87) "flag for the point is",flags
      if(prnstat.eq.1)then
C     write(10,972) tanbeta,m0,m12,a0,m10,m20,flags
         write(10,87) "flag for the point is",flags
      endif
      
      else
         
C     write(iunit,972) tanbeta,m0,m12,a0,m10,m20,flags
         write(iunit,87) "flag for the point is",flags
         if(prnstat.eq.1)then
C     write(10,972) tanbeta,m0,m12,a0,m10,m20,flags
            write(10,87) "flag for the point is",flags
         endif
         
         
      endif opt
      
C     -------------------------------------------------------------
C     The End
C     -------------------------------------------------------------      
      
      RETURN
      END subroutine SuSeFLAV
      
!====================================================================================
!
!     Flag to detect nan
!-----------------------------------------------------------------------------------

      Logical Function tachflag(num)

      Double precision  num

      tachflag = .NOT. (num.ge.0d0)

      end
!=========================================================================================================
C    GMSB boundary conditions adapted from hep-ph/9703211


      SUBROUTINE gmsb(gr,nhat,gmsbsusyb, gmsbmess,alph1,alph2,
     $     alph3,mQ0,mU0,mD0,mL0,mE0,mNU0,Mg1,Mg2,Mg3,m10,m20)
      
     
      DOUBLE PRECISION gmsbsusyb, gmsbmess,nhat,alph1,alph2,alph3
      DOUBLE PRECISION x,gx,fx
      DOUBLE PRECISION mq20, mq30, mu20, mu30, md20, md30, ml20,ml30
      DOUBLE PRECISION me20,me30, mnu20,mnu30, Mg1,Mg2,Mg3,m10,m20

      double precision mQ0(3,3),mU0(3,3),gravitino,gr
      double precision mD0(3,3),mL0(3,3),mE0(3,3),mNU0(3,3) !,b1,b2,b3

      DOUBLE PRECISION DDILOG

      EXTERNAL rungauge
      common/m32/gravitino
!-----------------------------
      x = gmsbsusyb/gmsbmess


      if(x.eq.1.d0) then
      gx = 1.38629
      else

      gx = 1.d0/(x*x) * ((1.d0 + x)*dlog(1.d0 + x) + 
     $        (1.d0 - x)*dlog(1.d0 - x))

      endif

      if(x.eq.1.d0)then
         fx = 0.702266
         else
      fx = (1.d0 + x)/(x*x) * (dlog(1.d0 + x) - 
     $           2.d0*DDILOG(x/(1.d0 +x)) + 
     $     0.5d0 * DDILOG(2.d0*x/(1.d0+x))) +
     $     (1.d0-x)/(x*x) * (dlog(1.d0-x) - 2.d0*DDILOG(-x/(1.d0-x)) + 
     $     0.5d0 * DDILOG(-2.d0*x/(1.d0-x)))

      endif
!-----------

         gravitino =  2.37d-19 * gmsbsusyb * gmsbmess * gr
 
c  gauginos 
c     

       Mg1 = nhat * alph1 * gmsbsusyb * gx
       Mg2 = nhat * alph2 * gmsbsusyb * gx
       Mg3 = nhat * alph3 * gmsbsusyb * gx

c Scalars



      mu30 = 2.d0 * gmsbsusyb**2.d0  * nhat * fx *
     $      (4.d0/3.d0 * alph3**2.d0 + 4.d0/15.d0 * alph1**2.d0)

      mu20 = mu30
!----------


      mq30 =  2.d0 * gmsbsusyb**2.d0 * nhat * fx *
     $     (4.d0/3.d0 * alph3**2.d0 + 3.d0/4.d0 * alph2**2.d0 +
     $     1.d0/60.d0 * alph1**2.d0)

      mq20 = mq30
!--------


      md30  = 2.d0 * gmsbsusyb**2.d0 * nhat * fx * 
     $     (4.d0/3.d0 * alph3**2 + 2.d0/30.d0 * alph1**2.d0)
      
      md20 = md30
!---------------



      ml30 = 2.d0 * gmsbsusyb**2.d0 * nhat * fx *
     $     (3.d0/4.d0 * alph2**2 + 3.d0/20.d0 * alph1**2.d0) 

      ml20 = ml30
!----------------


      me30 = 2.d0 * gmsbsusyb**2.d0 * nhat * fx *
     $     (3.d0/5.d0 * alph1**2.d0)

      me20 = me30

!------------
      mnu20 = me20
      mnu30 = me30
!----------------


C     Higgs 

      m10 = ml30
      m20 = ml30

      mQ0(1,1) = mq20
      mQ0(1,2) = 0.d0
      mQ0(1,3) = 0.d0
      mQ0(2,1) = 0.d0
      mQ0(2,2) = mq20
      mQ0(2,3) = 0.d0
      mQ0(3,1) = 0.d0
      mQ0(3,2) = 0.d0
      mQ0(3,3) = mq30

      mU0(1,1) = mu20
      mU0(1,2) = 0.d0
      mU0(1,3) = 0.d0
      mU0(2,1) = 0.d0
      mU0(2,2) = mu20
      mU0(2,3) = 0.d0
      mU0(3,1) = 0.d0
      mU0(3,2) = 0.d0
      mU0(3,3) = mu30

      mD0(1,1) = md20
      mD0(1,2) = 0.d0
      mD0(1,3) = 0.d0
      mD0(2,1) = 0.d0
      mD0(2,2) = md20
      mD0(2,3) = 0.d0
      mD0(3,1) = 0.d0
      mD0(3,2) = 0.d0
      mD0(3,3) = md30

      mL0(1,1) = ml20
      mL0(1,2) = 0.d0
      mL0(1,3) = 0.d0
      mL0(2,1) = 0.d0
      mL0(2,2) = ml20
      mL0(2,3) = 0.d0
      mL0(3,1) = 0.d0
      mL0(3,2) = 0.d0
      mL0(3,3) = ml30

      mE0(1,1) = me20
      mE0(1,2) = 0.d0
      mE0(1,3) = 0.d0
      mE0(2,1) = 0.d0
      mE0(2,2) = me20
      mE0(2,3) = 0.d0    
      mE0(3,1) = 0.d0
      mE0(3,2) = 0.d0    
      mE0(3,3) = me30
      
      mNU0(1,1) = mnu20
      mNU0(1,2) = 0.d0
      mNU0(1,3) = 0.d0
      mNU0(2,1) = 0.d0
      mNU0(2,2) = mnu20
      mNU0(2,3) = 0.d0
      mNU0(3,1) = 0.d0
      mNU0(3,2) = 0.d0
      mNU0(3,3) = mnu30


      END SUBROUTINE      
!================================================================================
      subroutine finetune(MZ,tanbeta,mur,mh1mz,mh2mz,delta1,delta2,
     $     ftmz,ftmt)
      
      implicit none
      double precision mur,MZ,mh1c,mh2c,tanbeta,mh2mz
      double precision ftmz,delta1,delta2,ftmt,mh1mz,tbsqr

      mh1c = mh1mz - delta1
      mh2c = mh2mz - delta2
      tbsqr = (tanbeta**2.d0 + 1.d0)/(tanbeta**2.d0 - 1.d0)
      
      ftmz = (2.d0 * mur**2.d0/(MZ*MZ)) * (1.d0 + 
     $     ((4.d0 * tanbeta**2.d0 * (mh1c - mh2c) * tbsqr)/
     $     (((mh1c - mh2c)*tbsqr - MZ*MZ) * 
     $     (tanbeta**2.d0 - 1.d0)**2.d0)))
     
      ftmt = ftmz*0.5d0 + (2.d0 * mur*mur)/((mh1c + mh2c) * 
     $     (tanbeta**2.d0 - 1.d0))

      return
      end subroutine finetune

!---------------------------------------
        
      subroutine deltarhoew(MW,MZ,piwwt0,pizzt0,piwwT,pizzT,drho)

      implicit none
      double precision piwwt0,pizzt0,MW,MZ,drho,piwwT,pizzT

      drho = ((pizzt0/(MZ*MZ - pizzT)) - (piwwt0/(MW*MW - piwwT)))  
      
      return
      end subroutine deltarhoew

!================================================
