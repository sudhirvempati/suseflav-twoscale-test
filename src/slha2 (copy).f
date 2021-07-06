!----------------------------------------------------
C     15.08.2010
!-----------------------------------------------------
* Program : SuSeFLAV: Supersymmetry seesaw and flavor violation calculator
* Version : 1.1.3
* Website : http://cts.iisc.ernet.in/Suseflav/main.html
* Authors : Debtosh Choudhury  debtosh@cts.iisc.ernet.in 
*           Raghuveer Garani   rgarani@cts.iisc.ernet.in
*           Sudhir Vempati     vempati@cts.iisc.ernet.in

      SUBROUTINE printslha2(errge,flags,msusyold,murge,mAm3,
     $     MT,mb,mtau,MW,MZ,mh1mz,mh2mz,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,ON,OCL,OCR,alph3,alph2,alph1,M1tz,
     $     M2tz,M3tz,M3t,thetat,thetab,thetatau,mh0sq,mHu0sq,
     $     mA0sq,mhpmsq,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mSQRG,mSURG,mSDRG,mSLRG,mSERG,USL,ftmz,ftmt,
     $     MWsqpole_MZ)

      
      implicit double precision (a-h,m,o-z)
      character*100 flags
      INTEGER lhaout,o1pos,o2pos,o3pos,o4pos,checkpi2,errge
      INTEGER lopt, rhn,mrhnum,qmix
      double precision extpar(0:60)
      character charprn(1:20)*20,softpar(0:60)*20
      Character(len=8) :: date
      Character(len=10) :: time


      CHARACTER*4 model
      CHARACTER*15 mrhn
      CHARACTER*3 case
      character*6 modeln
      DOUBLE PRECISION m0, m12, a0, tanbeta, sgnmu
      DOUBLE PRECISION msusyold, murge, mAm3, MT, mb, mtau
      DOUBLE PRECISION  mh1mz, mh2mz, vev1n,vev2n,MWsqpole_MZ
      DOUBLE PRECISION alph3, alph2, alph1, M1tz, M2tz, M3tz
      DOUBLE PRECISION thetat, thetab, thetatau, M3t
      DOUBLE PRECISION ON(4,4),OCL(2,2),OCR(2,2)
      DOUBLE PRECISION AURG(3,3),ADRG(3,3),AERG(3,3)
      DOUBLE PRECISION TURG(3,3),TDRG(3,3),TERG(3,3)
      DOUBLE PRECISION yuRG(3,3),ydRG(3,3),yeRG(3,3)
      DOUBLE PRECISION mh0sq, mA0sq, mhpmsq, mHu0sq
      DOUBLE PRECISION Ceg(2),Neg(4) !, Mneut(4), Mchar(2)
      DOUBLE PRECISION SUegg(6),SDegg(6),SLegg(6),SNegg(3),mSERG(3,3)
      DOUBLE PRECISION mSQRG(3,3),mSURG(3,3),mSDRG(3,3),mSLRG(3,3)
      DOUBLE PRECISION OCRm(2,2),OCLTm(2,2),ONm(4,4)
      DOUBLE PRECISION o1,o2,o3,o4,ONtemp(4,4),negarray(4),USL(6,6)
      double precision scale
        
      DOUBLE PRECISION MSQU3(2,2),MSQD3(2,2),MSQE3(2,2),
     $     th1,th2,th3,th4

!      double precision rrstaumu(2),USLsrt(6,6),USLTsrt(6,6)

      double precision sw2, gnuL,gul,gur,gdl,gdr,gel,ger,beta
      integer num,modnum,smnum(10),minparnum(10),counter_sminp
      integer counter_minpar,counter_algopar,algonum(20)
      double precision sminp(10), minpar(10), algopar(20)

      DOUBLE PRECISION alph,Gf,alphas,alphahiggs
      DOUBLE PRECISION gmsbsusyb, gmsbmess,gr,nhat
      DOUBLE PRECISION s1,s2,s3,s4
      DOUBLE PRECISION drho,gm2,Bbsg,bmeg,brmu3e,btmug,brtau3mu,
     $     bteg,brtau3e, ftmz,ftmt,gravitino
      DOUBLE PRECISION o1s,o2s,o3s,o4s
      DOUBLE PRECISION USLvd(6,6),USUvd(6,6),USDvd(6,6)
      DOUBLE PRECISION SLegvd(6),SUegvd(6),SDegvd(6)
      double precision SNegvd(3),USNvd(3,3)
      double precision gutscale,yukgut(126)
      double precision Mg1,Mg2,Mg3

      DOUBLE PRECISION mq11,mq12,mq13,mq21,mq22,mq23,mq31,mq32,mq33 
      DOUBLE PRECISION mu11,mu12,mu13,mu21,mu22,mu23,mu31,mu32,mu33 
      DOUBLE PRECISION md11,md12,md13,md21,md22,md23,md31,md32,md33 
      DOUBLE PRECISION ml11,ml12,ml13,ml21,ml22,ml23,ml31,ml32,ml33 
      DOUBLE PRECISION me11,me12,me13,me21,me22,me23,me31,me32,me33 
      DOUBLE PRECISION mnu11,mnu12,mnu13,mnu21,mnu22,mnu23,mnu31,
     $     mnu32,mnu33 
      DOUBLE PRECISION a0u11,a0u12,a0u13,a0u21,a0u22,a0u23,a0u31,
     $     a0u32,a0u33
      DOUBLE PRECISION a0d11,a0d12,a0d13,a0d21,a0d22,a0d23,a0d31,
     $     a0d32,a0d33
      DOUBLE PRECISION a0e11,a0e12,a0e13,a0e21,a0e22,a0e23,a0e31,
     $     a0e32,a0e33
      DOUBLE PRECISION a0nu11,a0nu12,a0nu13,a0nu21,a0nu22,a0nu23,
     $     a0nu31,a0nu32,a0nu33

      double precision mbpole, mtaupole, Mpole, MZpole, MWpole
!----------------------------------
      common/sminputs/ mbpole, mtaupole, Mpole, MZpole
      common/mwpole/MWpole
      common/unif/ gutscale,yukgut

      common/opt_mixing/OCRm,OCLTm,ONm,alphahiggs
      
!     common/opt_mixing/OCRTm, OCLm, ONm, alphahiggs
      common/vev_ewsb/ vev1n,vev2n
!     common/mutaup/rrstaumu,USLsrt,USLTsrt
      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
      common/gauge/alph,Gf,alphas
      common/charinputs/case, model
      common/gmsbinputs/ gmsbsusyb, gmsbmess,gr,nhat

      common/lowepar/drho,gm2,Bbsg,bmeg,brmu3e,btmug,brtau3mu,
     $     bteg,brtau3e

      common/lfvmixing/USLvd,USUvd,USDvd,USNvd
      common/lfveigval/SLegvd,SUegvd,SDegvd,SNegvd

      common/gauginos/Mg1,Mg2,Mg3

      common/nuaterms/a0u11,a0u12,a0u13,a0u21,a0u22,a0u23,a0u31,a0u32,
     $     a0u33,a0d11,a0d12,a0d13,a0d21,a0d22,a0d23,a0d31,a0d32,
     $     a0d33,a0e11,a0e12,a0e13,a0e21,a0e22,a0e23,a0e31,a0e32,
     $     a0e33,a0nu11,a0nu12,a0nu13,a0nu21,a0nu22,a0nu23,
     $     a0nu31,a0nu32,a0nu33

      common/nusoterms/mq11,mq12,mq13,mq21,mq22,mq23,
     $     mq31,mq32,mq33,mu11,mu12,mu13,mu21,mu22,mu23,mu31,
     $     mu32,mu33,md11,md12,md13,md21,md22,md23,md31,md32,
     $     md33,ml11,ml12,ml13,ml21,ml22,ml23,ml31,ml32,ml33, 
     $     me11,me12,me13,me21,me22,me23,me31,me32,me33,mnu11,
     $     mnu12,mnu13,mnu21,mnu22,mnu23,mnu31,mnu32,mnu33


!---SLHA common blocks---

      common/loops/ lopt,rhn
      common/slha_modsel/num, modnum
      common/slha_counters/counter_sminp,counter_minpar,counter_algopar
      common/slha_num/smnum,minparnum,algonum
      common/slha_variables/sminp,minpar,algopar
      common/m32/gravitino
      common/quarkmix/ qmix

!--------------------------------------------------------------------------------
      include 'stdinputs.h'
 
      lhaout = 222

      OPEN(lhaout, FILE = 'susy_flavor.in',ACCESS = 'APPEND',
     $     STATUS = 'replace')

!----------------------------------------------
C     

      pi = 4.d0*datan(1.d0)
      sw2 = 1.d0 - (MW / MZ)**2.d0
      gnuL = 0.5d0
      guL = 0.5d0 - 2.d0*sw2/ 3.d0 
      gdL = -0.5d0 + sw2/3.d0
      geL = -0.5d0 + sw2
      guR =  2.d0*sw2/3.d0 
      gdR = -sw2/3.d0
      geR = -sw2 

      beta = datan(tanbeta)

      do i = 1, 3
         do j = 1, 3
            TURG(i,j) = AURG(i,j)*YURG(i,j)*4.d0*pi
            TDRG(i,j) = ADRG(i,j)*YDRG(i,j)*4.d0*pi
            TERG(i,j) = AERG(i,j)*YERG(i,j)*4.d0*pi
         enddo
      enddo
      
C     for stop

      MSQU3(1,1) = mSQRG(3,3) + MT**2.d0 + 
     $           guL*MZ*MZ*dcos(2.d0*beta)

      MSQU3(1,2) = MT*((AURG(3,3) - murge/tanbeta))

      MSQU3(2,1) = MSQU3(1,2)

      MSQU3(2,2) = mSURG(3,3) + MT**2.d0 + 
     $     guR*MZ*MZ*dcos(2.d0*beta)



C     for sbottom

      MSQD3(1,1) = mSQRG(3,3) + MB**2.d0 + 
     $     gdL*MZ*MZ*dcos(2.d0*beta)

      MSQD3(1,2) = MB*((ADRG(3,3) - murge*tanbeta))
      MSQD3(2,1) = MSQD3(1,2)

      MSQD3(2,2) = mSDRG(3,3) + MB**2.d0 + 
     $     gdR*MZ*MZ*dcos(2.d0*beta)



C     for stau

      MSQE3(1,1) = mSLRG(3,3) + MTAU**2.d0 + 
     $     geL*MZ*MZ*dcos(2.d0*beta)
      MSQE3(1,2) = MTAU*((AERG(3,3) - murge*dtan(beta)))
      MSQE3(2,1) = MSQE3(1,2)
      MSQE3(2,2) = mSERG(3,3) + MTAU**2.d0 + 
     $     geR*MZ*MZ*dcos(2.d0*beta)



      thetat   = 0.d0 + 
     $     datan((2.d0*MSQU3(1,2))/(MSQU3(1,1)-MSQU3(2,2)))*0.5d0

      thetab   = 0.d0 + 
     $     datan((2.d0*MSQD3(1,2))/(MSQD3(1,1)-MSQD3(2,2)))*0.5d0

      thetatau = 0.d0 + 
     $     datan((2.d0*MSQE3(1,2))/(MSQE3(1,1)-MSQE3(2,2)))*0.5d0



      if(MSQE3(1,1).gt.MSQE3(2,2))then
         thetatau = thetatau + pi/2.d0
         checkpi2 = 1

      endif

      if(MSQU3(1,1).gt.MSQU3(2,2))then
         thetat = thetat + pi/2.d0
      endif

      if(MSQD3(1,1).gt.MSQD3(2,2))then
         thetab = thetab + pi/2.d0
      endif

!-----------------------------------------------


c PDG values:
      id =1
      idb=-1
      iu =2
      iub=-2
      is =3
      isb=-3
      ic =4
      icb=-4
      ib =5
      ibb=-5
      it =6
      itb=-6

      ie   =11
      ine  =12
      imu  =13
      inmu =14
      itau =15
      intau=16

      ihl=25
      ihh=35
      iha=36
      ihc=37
      igl=21
      iga=22
      iz =23
      iwc=24

      isdl=1000001
      isdr=2000001
      isul=1000002
      isur=2000002
      issl=1000003
      issr=2000003
      iscl=1000004
      iscr=2000004
      isb1=1000005
      isb2=2000005
      ist1=1000006
      ist2=2000006

      iglo=1000021
      in1 =1000022
      in2 =1000023
      in3 =1000025
      in4 =1000035
      ic1 =1000024
      ic2 =1000037

      intau1=1000016 
      intau2=2000016 
      inel  =1000012
      iner  =2000012
      inmul =1000014
      inmur =2000014
      
      isell =1000011
      iselr =2000011
      ismul =1000013
      ismur =2000013
      istau1=1000015
      istau2=2000015

      igrav =1000039
!---------------------------------      
      charprn(3) = ' tanbeta(mz)'
      charprn(4) = ' sign(mu)'

c input for msugra models
      if(model.eq.'mSUG')then
      charprn(1) = ' m0'
      charprn(2) = ' m_1/2'
      charprn(5) = ' A0'
      else if(model.eq.'GMSB')then
c    input for GMSB models:
      charprn(1) = ' Lambda_susy'
      charprn(2) = ' Lambda_mess'
      charprn(5) = ' N_mess'  
      charprn(6) = ' c_grav'
c    input for NUHM models:
      else if(model.eq.'NUHM')then
      charprn(1) = ' m0'
      charprn(2) = ' m_1/2'
      charprn(5) = ' A0'
      charprn(6) = ' m10'
      charprn(7) = ' m20'
c    input for NUGM models:
      else if(model.eq.'NUGM')then
      charprn(1) = ' m0'
!      charprn(2) = ' m_1/2'
      charprn(5) = ' A0'
      charprn(6) = ' m10'
      charprn(7) = ' m20'

      endif
            
c  SOFTPAR :
        softpar(0) = ' EWSB scale'          
        softpar(10) = ' GUT scale'
	softpar(23) = ' mu(EWSB)'
	softpar(24) = ' m^2_A_run(EWSB)'
        softpar(25) = ' tanbeta(in)'
	softpar(26) = ' MA_pole'
        softpar(1) = ' M_1'
        softpar(2) = ' M_2'
        softpar(3) = ' M_3'
        softpar(21) = ' M^2_Hd'
        softpar(22) = ' M^2_Hu'
        softpar(31) = ' M_eL'
        softpar(32) = ' M_muL'
        softpar(33) = ' M_tauL'
        softpar(34) = ' M_eR'
        softpar(35) = ' M_muR'
        softpar(36) = ' M_tauR'
        softpar(41) = ' M_q1L'
        softpar(42) = ' M_q2L'
        softpar(43) = ' M_q3L'
        softpar(44) = ' M_uR'
        softpar(45) = ' M_cR'
        softpar(46) = ' M_tR'
        softpar(47) = ' M_dR'
        softpar(48) = ' M_sR'
        softpar(49) = ' M_bR'
        softpar(11) = ' A_t'
        softpar(12) = ' A_b'
        softpar(13) = ' A_tau'
        softpar(14) = ' A_u'
        softpar(15) = ' A_d'
        softpar(16) = ' A_e'

        scale = msusyold
!--------------------------------------------------------


      if(model.eq.'mSUG')then
         modeln = 'mSUGRA'
         modnum = 1
         mrhn = 'mSUGRA + RHN'
         mrhnum = 4
      endif

      if(model.eq.'GMSB')then
         modeln = 'GMSB'
         modnum = 2
         mrhn = 'GMSB + RHN'
         mrhnum = 5
      endif

      if(model.eq.'NUHM')then
         modeln = 'NUHM'
         modnum = 1
         mrhn = 'NUHM + RHN'
         mrhnum = 6
      endif

      if(model.eq.'NUGM')then
         modeln = 'NUGM'
         modnum = 1
         mrhn = 'NUGM + RHN'
         mrhnum = 7
      endif

      if(model.eq.'CNUM')then
         modeln = 'CNUM'
         modnum = 1
         mrhn = 'CNUM + RHN'
         mrhnum = 8
      endif

      Call Date_and_time(date,time)

!---------------writing begins

!      Write(lhaout,106) '# Spectrum Output in SUSY Les Houches Accord 2'
!      Write(lhaout,106) '# SuSeFLAV v1.2.0 '
!      Write(lhaout,107) "# D. Chowdhury, R. Garani and S. K. Vempati,",
!     $     " hep-ph/1109.3551 "
!      Write(lhaout,107) '# For bug reports or any other queries please', 
!     $     ' send email to suseflav@cts.iisc.ernet.in '
!      Write(lhaout,106) '# Created on '//date(7:8)//'.'//date(5:6)//'.'
!     $     //date(1:4)// ' at '// time(1:2)//':'//time(3:4)//' Hrs'
!      write(lhaout,105)
!      write(lhaout,51) 'SPINFO',' Program information'
!      write(lhaout,61) 1,'SuSeFLAV     # Spectrum calculator'
!      write(lhaout,61) 2,'1.2.0        # Version number'
      if(errge.eq.1)then
         write(lhaout,61) 3,flags
      endif

      write(lhaout,105)
!     if(rhn.eq.0)then
      write(lhaout,51) 'MODSEL',' Select model'
      write(lhaout,611) 1, 0 , 'General MSSM'
      write(lhaout,611) 3, 0 , 'MSSM particle content'
      write(lhaout,611) 4, 0 , 'R-parity conserving MSSM'
      write(lhaout,611) 5, 2 , 'CP violated'
      write(lhaout,611) 6, 3 , 'Lepton and quark flavor violated'

      write(lhaout,51)'SOFTINP', 'Choose convention for the soft terms'
      write(lhaout,611) 1, 1, 'iconv (conventions for SLHA2 )'
      write(lhaout,611) 2, 2, 'input_type(dimension of soft mass )'
      write(lhaout,611) 3, 2, 'ilev (level of chiral resummation)'

      write(lhaout,51) 'SMINPUTS',' Standard Model inputs'
      write(lhaout,52) 1, 127.934d0, 'alpha^(-1) SM MSbar(MZ)'
      write(lhaout,52) 3, 0.1181d0, ' alpha_s(MZ) SM MSbar'
      write(lhaout,52) 4, 91.1876d0, 'MZ(pole)'
      write(lhaout,52) 5, 4.18d0, 'mb(mb) SM MSbar'
      write(lhaout,52) 6, 173.1d0, 'mtop(pole)' 
      write(lhaout,52) 7, 1.77684d0, 'mtau(pole)'
      write(lhaout,52) 11, 0.000511d0, 'me(pole)'
      write(lhaout,52) 13, 0.105658d0, 'mmu(pole)'
      write(lhaout,52) 21, 0.0047d0, 'md(2 GeV) MSbar'
      write(lhaout,52) 22, 0.0021d0, 'mu(2 GeV) MSbar'
      write(lhaout,52) 23, 0.0934d0, 'ms(2 GeV) MSbar'
      write(lhaout,52) 24, 1.279d0, 'mc(mc) MSbar'
      write(lhaout,52) 30, 80.398d0, 'MW (pole), '
      write(lhaout,52) 31, 0.23116d0, 's_W^2 (MSbar)'

      write(lhaout,51) 'VCKMIN',' CKM matrix'
      write(lhaout,52) 1, 0.2258d0, 'lambda'
      write(lhaout,52) 2, 0.808d0, 'A'
      write(lhaout,52) 3, 0.177d0, 'rho bar'
      write(lhaout,52) 4, 0.36d0, 'eta bar'
      errspec: if(errge.eq.1)then
!         write(lhaout,*),"Invalid point ", flags
         
      else errspec
      write(lhaout,51) 'EXTPAR',' input at susy scale, real part'
      write(lhaout,52) 0, -1.d0, 'input scale' 
      write(lhaout,52) 1,M1tz,'Re(m1), U(1) gaugino mass'
      write(lhaout,52) 2,M2tz,'Re(m2), SU(2) gaugino mass'
      write(lhaout,52) 3,M3tz,'m3, SU(3) gaugino mass'
      write(lhaout,55) 23,sgnmu*murge,'Re(mu)'
      write(lhaout,55) 25,vev2n/vev1n,'tan(beta)'
      write(lhaout,55) 26, sqrt(mAm3) ,'MA'

      write(lhaout,51) 'IMEXTPAR', 'imaginary part'
      write(lhaout,52) 1, 0.0, 'Im(m1), U(1) gaugino mass'
      write(lhaout,52) 2, 0.0, 'Im(m2), SU(2) gaugino mass'
      write(lhaout,52) 23, 0.0, 'Im(mu)'

      write(lhaout,51)'MSL2IN','Left softslepton mass matrix real part'
         write(lhaout,53) 1, 1, mSLRG(1,1), 'mSLRG_11'
         write(lhaout,53) 2, 2, mSLRG(2,2), 'mSLRG_22'
         write(lhaout,53) 3, 3, mSLRG(3,3), 'mSLRG_33'
         write(lhaout,53) 1, 2, mSLRG(1,2), 'mSLRG_12'
         write(lhaout,53) 2, 3, mSLRG(2,3), 'mSLRG_23'
         write(lhaout,53) 1, 3, mSLRG(1,3), 'mSLRG_13'

      write(lhaout,51)'IMMSL2IN',' imaginary part'
      write(lhaout,53) 1, 2, 0.0, 'mSLRG_12'
      write(lhaout,53) 2, 3, 0.0, 'mSLRG_23'
      write(lhaout,53) 1, 3, 0.0, 'mSLRG_13'

      write(lhaout,51)'MSE2IN','Right slepton mass matrix real part'
         write(lhaout,53) 1, 1, mSERG(1,1), 'mSERG_11'
         write(lhaout,53) 2, 2, mSERG(2,2), 'mSERG_22'
         write(lhaout,53) 3, 3, mSERG(3,3), 'mSERG_33'
         write(lhaout,53) 1, 2, mSERG(1,2), 'mSERG_12'
         write(lhaout,53) 2, 3, mSERG(2,3), 'mSERG_23'
         write(lhaout,53) 1, 3, mSERG(1,3), 'mSERG_13'

      write(lhaout,51)'IMMSE2IN',' imaginary part'
      write(lhaout,53) 1, 2,  0.0, 'mSERG_12'
      write(lhaout,53) 2, 3,  0.0, 'mSERG_23'
      write(lhaout,53) 1, 3,  0.0, 'mSERG_32'

      write(lhaout,51)'MSQ2IN','Left squark mass matrix real part'
         write(lhaout,53) 1, 1, mSQRG(1,1), 'mSQRG_11'
         write(lhaout,53) 2, 2, mSQRG(2,2), 'mSQRG_22'
         write(lhaout,53) 3, 3, mSQRG(3,3), 'mSQRG_33'
         write(lhaout,53) 1, 2, mSQRG(1,2), 'mSQRG_12'
         write(lhaout,53) 2, 3, mSQRG(2,3), 'mSQRG_23'
         write(lhaout,53) 1, 3, mSQRG(1,3), 'mSQRG_13'

      write(lhaout,51)'IMMSQ2IN',' imaginary part'
      write(lhaout,53) 1, 2,  0.0, 'mSQRG_12'
      write(lhaout,53) 2, 3,  0.0, 'mSQRG_23'
      write(lhaout,53) 1, 3,  0.0, 'mSQRG_13'

      write(lhaout,51)'MSU2IN','Right up squark mass matrix real part'
         write(lhaout,53) 1, 1, mSURG(1,1), 'mSURG_11'
         write(lhaout,53) 2, 2, mSURG(2,2), 'mSURG_22'
         write(lhaout,53) 3, 3, mSURG(3,3), 'mSURG_33'
         write(lhaout,53) 1, 2, mSURG(1,2), 'mSURG_12'
         write(lhaout,53) 2, 3, mSURG(2,3), 'mSURG_23'
         write(lhaout,53) 1, 3, mSURG(1,3), 'mSURG_13'

      write(lhaout,51)'IMMSU2IN',' imaginary part'
      write(lhaout,53) 1, 2,  0.0, 'mSURG_12'
      write(lhaout,53) 2, 3,  0.0, 'mSURG_23'
      write(lhaout,53) 1, 3,  0.0, 'mSURG_13'

      write(lhaout,51)'MSD2IN','Right downsquark mass matrix real part'
         write(lhaout,53) 1, 1, mSDRG(1,1), 'mSDRG_11'
         write(lhaout,53) 2, 2, mSDRG(2,2), 'mSDRG_22'
         write(lhaout,53) 3, 3, mSDRG(3,3), 'mSDRG_33'
         write(lhaout,53) 1, 2, mSDRG(1,2), 'mSDRG_12'
         write(lhaout,53) 2, 3, mSDRG(2,3), 'mSDRG_23'
         write(lhaout,53) 1, 3, mSDRG(1,3), 'mSDRG_13'

      write(lhaout,51)'IMMSD2IN',' imaginary part'
      write(lhaout,53) 1, 2,  0.0, 'mSDRG_12'
      write(lhaout,53) 2, 3,  0.0, 'mSDRG_23'
      write(lhaout,53) 1, 3,  0.0, 'mSDRG_13'

      write(lhaout,51)'TEIN ','slepton trilinear mixing, real part'
         write(lhaout,53) 1, 1, AERG(1,1), 'AERG_11'
         write(lhaout,53) 2, 2, AERG(2,2), 'AERG_22'
         write(lhaout,53) 3, 3, AERG(3,3), 'AERG_33'
         write(lhaout,53) 1, 2, AERG(1,2), 'AERG_12'
         write(lhaout,53) 2, 1, AERG(2,1), 'AERG_21'
         write(lhaout,53) 2, 3, AERG(2,3), 'AERG_23'
         write(lhaout,53) 3, 2, AERG(3,2), 'AERG_32'
         write(lhaout,53) 1, 3, AERG(1,3), 'AERG_13'
         write(lhaout,53) 3, 1, AERG(3,1), 'AERG_31'

      write(lhaout,51)'IMTEIN','slepton trilinear mixing,imaginarypart'
         write(lhaout,53) 1, 1, 0.0, 'AERG_11'
         write(lhaout,53) 2, 2, 0.0, 'AERG_22'
         write(lhaout,53) 3, 3, 0.0, 'AERG_33'
         write(lhaout,53) 1, 2, 0.0, 'AERG_12'
         write(lhaout,53) 2, 1, 0.0, 'AERG_21'
         write(lhaout,53) 2, 3, 0.0, 'AERG_23'
         write(lhaout,53) 3, 2, 0.0, 'AERG_32'
         write(lhaout,53) 1, 3, 0.0, 'AERG_13'
         write(lhaout,53) 3, 1, 0.0, 'AERG_31'

      write(lhaout,51)'TUIN ','up trilinear mixing, real part'
         write(lhaout,53) 1, 1, AURG(1,1), 'AURG_11'
         write(lhaout,53) 2, 2, AURG(2,2), 'AURG_22'
         write(lhaout,53) 3, 3, AURG(3,3), 'AURG_33'
         write(lhaout,53) 1, 2, AURG(1,2), 'AURG_12'
         write(lhaout,53) 2, 1, AURG(2,1), 'AURG_21'
         write(lhaout,53) 2, 3, AURG(2,3), 'AURG_23'
         write(lhaout,53) 3, 2, AURG(3,2), 'AURG_32'
         write(lhaout,53) 1, 3, AURG(1,3), 'AURG_13'
         write(lhaout,53) 3, 1, AURG(3,1), 'AURG_31'

      write(lhaout,51)'IMTUIN','up trilinear mixing,imaginarypart'
         write(lhaout,53) 1, 1, 0.0, 'AURG_11'
         write(lhaout,53) 2, 2, 0.0, 'AURG_22'
         write(lhaout,53) 3, 3, 0.0, 'AURG_33'
         write(lhaout,53) 1, 2, 0.0, 'AURG_12'
         write(lhaout,53) 2, 1, 0.0, 'AURG_21'
         write(lhaout,53) 2, 3, 0.0, 'AURG_23'
         write(lhaout,53) 3, 2, 0.0, 'AURG_32'
         write(lhaout,53) 1, 3, 0.0, 'AURG_13'
         write(lhaout,53) 3, 1, 0.0, 'AURG_31'

      write(lhaout,51)'TDIN ','down trilinear mixing, real part'
         write(lhaout,53) 1, 1, ADRG(1,1), 'ADRG_11'
         write(lhaout,53) 2, 2, ADRG(2,2), 'ADRG_22'
         write(lhaout,53) 3, 3, ADRG(3,3), 'ADRG_33'
         write(lhaout,53) 1, 2, ADRG(1,2), 'ADRG_12'
         write(lhaout,53) 2, 1, ADRG(2,1), 'ADRG_21'
         write(lhaout,53) 2, 3, ADRG(2,3), 'ADRG_23'
         write(lhaout,53) 3, 2, ADRG(3,2), 'ADRG_32'
         write(lhaout,53) 1, 3, ADRG(1,3), 'ADRG_13'
         write(lhaout,53) 3, 1, ADRG(3,1), 'ADRG_31'

      write(lhaout,51)'IMTDIN','down trilinear mixing,imaginarypart'
         write(lhaout,53) 1, 1, 0.0, 'ADRG_11'
         write(lhaout,53) 2, 2, 0.0, 'ADRG_22'
         write(lhaout,53) 3, 3, 0.0, 'ADRG_33'
         write(lhaout,53) 1, 2, 0.0, 'ADRG_12'
         write(lhaout,53) 2, 1, 0.0, 'ADRG_21'
         write(lhaout,53) 2, 3, 0.0, 'ADRG_23'
         write(lhaout,53) 3, 2, 0.0, 'ADRG_32'
         write(lhaout,53) 1, 3, 0.0, 'ADRG_13'
         write(lhaout,53) 3, 1, 0.0, 'ADRG_31'

      write(lhaout,511) 'SFLAV_HADRON'
      write(lhaout,52) 1, 0.156d0, 'f_K'
      write(lhaout,52) 2, 0.2d0,   'f_D'
      write(lhaout,52) 3, 0.193d0, 'f_B_d'
      write(lhaout,52) 4, 0.232d0, 'f_B_s'
      write(lhaout,52) 5, 0.724d0, 'B_K for SM contribution to KKbar'
      write(lhaout,52) 6, 1.87d0,  'eta_cc in KK mixing (SM)'
      write(lhaout,52) 7, 0.496d0, 'eta_ct in KK mixing (SM)'
      write(lhaout,52) 8, 0.577d0, 'eta_ct in KK mixing (SM)'
      write(lhaout,52) 9, 2.d0,    'scale for B_K (non-SM)' 
      write(lhaout,52) 10, 0.61d0, 'B_K for VLL (non-SM)'
      write(lhaout,52) 11, 0.76d0, 'B_K for SLL1'
      write(lhaout,52) 12, 0.51d0, 'B_K for SLL2'
      write(lhaout,52) 13, 0.96d0, 'B_K for LR1'
      write(lhaout,52) 14, 1.30d0, 'B_K for LR2'
      write(lhaout,52) 15, 1.d0,   'B_D for SM contribution'
      write(lhaout,52) 16, 2.d0,   'scale for B_D (non-SM)'
      write(lhaout,52) 17, 1.d0,   'B_D for VLL'
      write(lhaout,52) 18, 1.d0,   'B_D for SLL1'
      write(lhaout,52) 19, 1.d0,   'B_D for SLL2' 
      write(lhaout,52) 20, 1.d0,   'B_D for LR1' 
      write(lhaout,52) 21, 1.d0,   'B_D for LR2'
      write(lhaout,52) 22, 1.22d0, 'B_Bd for SM contribution' 
      write(lhaout,52) 23, 4.6d0,  'scaleforB_B(non-SM,both Bd and Bs)' 
      write(lhaout,52) 24, 0.87d0, 'B_Bd for VLL (non-SM)'
      write(lhaout,52) 25, 0.8d0,  'B_Bd for SLL1'
      write(lhaout,52) 26, 0.71d0, 'B_Bd for SLL2'
      write(lhaout,52) 27, 1.71d0, 'B_Bd for LR1'
      write(lhaout,52) 28, 1.16d0, 'B_Bd for LR2'
      write(lhaout,52) 29, 1.22d0, 'B_Bs for SM contribution'
      write(lhaout,52) 30, 0.55d0, 'eta_b for BsBs (SM)'
      write(lhaout,52) 31, 0.87d0, 'B_Bs for VLL (non-SM)'
      write(lhaout,52) 32, 0.8d0,  'B_Bs for SLL1'
      write(lhaout,52) 33, 0.71d0, 'B_Bs for SLL2'
      write(lhaout,52) 34, 1.71d0, 'B_Bs for LR1'
      write(lhaout,52) 35, 1.16d0, 'B_Bs for LR2' 
      write(lhaout,52) 36, 1.519d-12, 'Bd lifetime (experimental)'
      write(lhaout,52) 37, 1.512d-12, 'Bs lifetime (experimental)'
      write(lhaout,52) 38, 5.27958d0, 'Bd mass (experimental)'
      write(lhaout,52) 39, 5.36677d0, 'Bs mass (experimental)'
      write(lhaout,52) 40, 3.337d-13, 'Delta Bd (experimental)'
      write(lhaout,52) 41, 1.17d-11, 'Delta Bs (experimental)'
      write(lhaout,52) 42, 0.497614d0, 'K0 mass (experimental)'
      write(lhaout,52) 43, 3.483d-15,  'Delta mK (experimental)'
      write(lhaout,52) 44, 2.229d-3,   'eps_K (experimental)'
      write(lhaout,52) 45, 1.8645d0,   'D0 mass (experimental)'
      write(lhaout,52) 46, 1.56d-14,   'Delta mD (experimental)'
      write(lhaout,52) 47, 2.231d-10,  'parameter kappa in K^0->pi^0vv'
      write(lhaout,52) 48, 5.173d-11,  'parameter kappa in K^+->pi^+vv'
      write(lhaout,52) 49, 0.41d0,     'parameter P_c in K->pivv'
      write(lhaout,52) 50, 0.013d-10,  'error of ak0'
      write(lhaout,52) 51, 0.024d-11,  'error of akp'
      write(lhaout,52) 52, 0.03d0,     'error of pc' 
      write(lhaout,52) 53, 0.79d0,     'neutron EDM_d QCD coefficient'
      write(lhaout,52) 54, -0.2d0,     'neutron EDM_u QCD coefficient'
      write(lhaout,52) 55, 0.59d0,     'neutron CDM_d QCD coefficient'
      write(lhaout,52) 56, 0.3d0,      'neutron CDM_u QCD coefficient'
      write(lhaout,52) 57, 3.4d0,      'neutron CDM_g QCD coefficient'
      write(lhaout,52) 58, 1.18d0,     'neutron EDM chiral symmetry '
      write(lhaout,52) 59, 1.5d0,      'pole c quark mass'
      write(lhaout,52) 60, 0.1872d0,   'Br(tau->evv)'  
      write(lhaout,52) 61, 5.27917d0,  'M_B+'
      write(lhaout,52) 62,0.297d0,'Br(B->D tau nu)/Br(B->D l nu)in SM'
      write(lhaout,52) 63,0.017d0,'errorof Br(B->D taunu)/Br(B->Dlnu)'
      write(lhaout,52) 64,0.252d0,'Br(B->D* tau nu)/Br(B->D* l nu) '
      write(lhaout,52) 65,0.003d0,'errorofBr(B->D*taunu)/Br(B->D*lnu)'
      endif errspec

 50   format('#',1x,A)
 51   format('BLOCK',1x,A,2x,'#',1x,A)
 52   format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
 522  format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
 53   format(1x,I2,1x,I2,3x,1P,E16.8,0P,3x,'#',1x,A)
 54   format('BLOCK',1x,A,1P,E16.8,2x,'#',1x,A)
 55   format(1x,I5,3x,1P,E16.8,0P,3x,'#',1x,A)
 60   format(9x,1P,E16.8,0P,3x,'#',1x,A)
 61   format(1x,I5,3x,A)
 611  format(1x,I5,2x,I5,3x,'#',1x,A)
 62   format(1x,I5,1x,I5,3x,'#',A)
 72   format(1x,I9,3x,1P,E16.8,0P,3x,'#',1x,A)
 105  format('#') 
 106  format(A) 
 107  format(A,A) 
 108  format('BLOCK',1x,A,1P,E16.8,2x,'#',1x,A,A)
 109  format(9x,1P,E16.8,0P,3x,'#',1x,A,A)
 511   format('BLOCK',1x,A)

      close(lhaout)
c -------------------------------------------------------------------------
      end subroutine
!=====================================================================================================
