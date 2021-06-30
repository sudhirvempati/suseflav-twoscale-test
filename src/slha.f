!----------------------------------------------------
C     15.08.2010
!-----------------------------------------------------
* Program : SuSeFLAV: Supersymmetry seesaw and flavor violation calculator
* Version : 1.1.3
* Website : http://cts.iisc.ernet.in/Suseflav/main.html
* Authors : Debtosh Choudhury  debtosh@cts.iisc.ernet.in 
*           Raghuveer Garani   rgarani@cts.iisc.ernet.in
*           Sudhir Vempati     vempati@cts.iisc.ernet.in

      SUBROUTINE printslha(errge,flags,msusyold,murge,mAm3,
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
 
      lhaout = 2

      OPEN(lhaout, FILE = 'slha.out',ACCESS = 'APPEND',
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
      Write(lhaout,106) '# Spectrum Output in SUSY Les Houches Accord 2'
      Write(lhaout,106) '# SuSeFLAV v1.2.0 '
      Write(lhaout,107) "# D. Chowdhury, R. Garani and S. K. Vempati,",
     $     " hep-ph/1109.3551 "
      Write(lhaout,107) '# For bug reports or any other queries please', 
     $     ' send email to suseflav@cts.iisc.ernet.in '
      Write(lhaout,106) '# Created on '//date(7:8)//'.'//date(5:6)//'.'
     $     //date(1:4)// ' at '// time(1:2)//':'//time(3:4)//' Hrs'
      write(lhaout,105)
      write(lhaout,51) 'SPINFO',' Program information'
      write(lhaout,61) 1,'SuSeFLAV     # Spectrum calculator'
      write(lhaout,61) 2,'1.2.0        # Version number'
      if(errge.eq.1)then
         write(lhaout,61) 3,flags
      endif

      write(lhaout,105)
!     if(rhn.eq.0)then
      write(lhaout,51) 'MODSEL',' MODEL NAME'
      write(lhaout,611) 1, modnum , modeln

      if(qmix.eq.1.and.rhn.eq.0)then
         write(lhaout,611) 6, 1, 'Quark flavor is violated'
      endif

      if(rhn.eq.1)then
!         write(lhaout,51) 'MODSEL',' MODEL NAME'
!         write(lhaout,611) 1, modnum , modeln
         write(lhaout,611) 3, mrhnum , mrhn
         write(lhaout,611) 6, 3 , 'Lepton and quark flavor is violated'
      endif
      
      if(model.eq.'GMSB')then

      write(lhaout,105)
      write(lhaout,51) 'MINPAR','Input parameters'
      write(lhaout,52) 1, gmsbsusyb ,charprn(1)
      write(lhaout,52) 2, gmsbmess ,charprn(2)
      write(lhaout,52) 3, tanbeta ,charprn(3)
      write(lhaout,52) 4, sgnmu ,charprn(4)
      write(lhaout,52) 5, nhat ,charprn(5)      
      write(lhaout,52) 6, gr ,charprn(6)

      else if(model.eq.'mSUG')then
 
      write(lhaout,105)
      write(lhaout,51) 'MINPAR','Input parameters'
      write(lhaout,52) 1, m0 ,charprn(1)
      write(lhaout,52) 2, m12 ,charprn(2)
      write(lhaout,52) 3, tanbeta ,charprn(3)
      write(lhaout,52) 4, sgnmu ,charprn(4)
      write(lhaout,52) 5, a0 ,charprn(5)

      else if(model.eq.'NUHM')then

      write(lhaout,105)
      write(lhaout,51) 'MINPAR','Input parameters'
      write(lhaout,52) 1, m0 ,charprn(1)
      write(lhaout,52) 2, m12 ,charprn(2)
      write(lhaout,52) 3, tanbeta ,charprn(3)
      write(lhaout,52) 4, sgnmu ,charprn(4)
      write(lhaout,52) 5, a0 ,charprn(5)
c$$$      write(lhaout,52) 6, m10 ,charprn(6)
c$$$      write(lhaout,52) 7, m20 ,charprn(7)

      else if(model.eq.'NUGM')then

      write(lhaout,105)
      write(lhaout,51) 'MINPAR','Input parameters'
      write(lhaout,52) 1, m0 ,charprn(1)
!      write(lhaout,52) 2, m12 ,charprn(2)
      write(lhaout,52) 3, tanbeta ,charprn(3)
      write(lhaout,52) 4, sgnmu ,charprn(4)
      write(lhaout,52) 5, a0 ,charprn(5)

      else if(model.eq.'CNUM')then

      write(lhaout,105)

      endif

c ----------------------- c
c The SM input parameters c
c ----------------------- c

      write(lhaout,105)
      write(lhaout,51) 'SMINPUTS','Standard Model inputs'
      write(lhaout,52) 1, 1/alph  ,'alpha_em (M_Z)^MSbar'
      write(lhaout,52) 2, Gf    ,'G_F [GeV^-2]'
      write(lhaout,52) 3, alphas,'alpha_S(M_Z)^MSbar'
      write(lhaout,52) 4, MZ    ,'M_Z pole mass'
      write(lhaout,52) 5, mbpole    ,'mb(mb)^MSbar'
      write(lhaout,52) 6, Mpole    ,'mt pole mass'
      write(lhaout,52) 7, mtaupole  ,'mtau pole mass'

      errspec: if(errge.eq.1)then
!         write(lhaout,*),"Invalid point ", flags
         
      else errspec

c$$$  extpar(0) = scale
c$$$  softpar(0) ='EWSB scale' 

      write(lhaout,105)

      if(model.ne.'GMSB')then
         write(lhaout,51) 'EXTPAR','Extra Input parameters'
         write(lhaout,72) 0,gutscale,'Unification Scale'
      endif

      if(model.eq.'NUHM')then
         write(lhaout,72) 21,sign(1.d0,m10)*m10*m10,'m^2_Hd'
         write(lhaout,72) 22,sign(1.d0,m20)*m20*m20,'m^2_Hu'
      endif

      if(model.eq.'NUGM')then
         write(lhaout,72) 1,Mg1,'M_1'
         write(lhaout,72) 2,Mg2,'M_2'
         write(lhaout,72) 3,Mg3,'M_3'
      endif

      if(model.eq.'CNUM')then

         write(lhaout,72) 1,Mg1,'M_1'
         write(lhaout,72) 2,Mg2,'M_2'
         write(lhaout,72) 3,Mg3,'M_3'

         write(lhaout,72) 11,a0u33,'A_t'
         write(lhaout,72) 12,a0d33,'A_b'
         write(lhaout,72) 13,a0e33,'A_tau'
         
         write(lhaout,72) 21,sign(1.d0,m10)*m10*m10,'m^2_Hd'
         write(lhaout,72) 22,sign(1.d0,m20)*m20*m20,'m^2_Hu'

         write(lhaout,72) 31,ml11,'m_eL'
         write(lhaout,72) 32,ml22,'m_muL'
         write(lhaout,72) 33,ml33,'m_tauL'
         write(lhaout,72) 34,me11,'m_eR'
         write(lhaout,72) 35,me22,'m_muR'
         write(lhaout,72) 36,me33,'m_tauR'
         write(lhaout,72) 41,mq11,'m_q1L'
         write(lhaout,72) 42,mq22,'m_q2L'
         write(lhaout,72) 43,mq33,'m_q3L'
         write(lhaout,72) 44,mu11,'m_uR'
         write(lhaout,72) 45,mu22,'m_cR'
         write(lhaout,72) 46,mu33,'m_tR'
         write(lhaout,72) 47,md11,'m_dR'
         write(lhaout,72) 48,md22,'m_sR'
         write(lhaout,72) 49,md33,'m_bR'
        

      endif

c ----------------- c
c The mass spectrum c
c ----------------- c
      
      write(lhaout,105)
      write(lhaout,51) 'MASS','Mass Spectrum'
      write(lhaout,50) 'PDG code           mass       particle'
      write(lhaout,52) iwc,dsqrt(MWsqpole_MZ),'W+'
      write(lhaout,52) ihl,dsqrt(mh0sq),'h'
      write(lhaout,52) ihh,dsqrt(mHu0sq),'H'
      write(lhaout,52) iha,dsqrt(mA0sq),'A'
      write(lhaout,52) ihc,dsqrt(mhpmsq),'H+'

c$$$      write(lhaout,52) ib,mbpole,'b pole mass calculated 
c$$$     $from mb(mb)_MSbar'

      if(rhn.eq.0)then

         write(lhaout,52) isdl,dsqrt(SDegg(6)),'~d_L'
         write(lhaout,52) isdr,dsqrt(SDegg(5)),'~d_R'
         write(lhaout,52) isul,dsqrt(SUegg(6)),'~u_L'
         write(lhaout,52) isur,dsqrt(SUegg(5)),'~u_R'
         write(lhaout,52) issl,dsqrt(SDegg(4)),'~s_L'
         write(lhaout,52) issr,dsqrt(SDegg(3)),'~s_R'
         write(lhaout,52) iscl,dsqrt(SUegg(4)),'~c_L'
         write(lhaout,52) iscr,dsqrt(SUegg(3)),'~c_R'
         write(lhaout,52) isb1,dsqrt(SDegg(1)),
     $        '~b_1'
         write(lhaout,52) isb2,dsqrt(SDegg(2)),
     $        '~b_2'
         write(lhaout,52) ist1,dsqrt(SUegg(1)),
     $        '~t_1'
         write(lhaout,52) ist2,dsqrt(SUegg(2)),
     $        '~t_2'
         write(lhaout,52) isell,dsqrt(SLegg(6)),'~e_L'
         write(lhaout,52) iselr,dsqrt(SLegg(5)),'~e_R'
         write(lhaout,52) inel, dsqrt(SNegg(3)),'~nu_eL'
         write(lhaout,52) ismul,dsqrt(SLegg(4)),
     $        '~mu_L'
         write(lhaout,52) ismur,dsqrt(SLegg(3)),
     $        '~mu_R'
         write(lhaout,52) inmul,dsqrt(SNegg(2)),'~nu_muL'
         write(lhaout,52) istau1,dsqrt(SLegg(1)),
     $        '~tau_1'
         write(lhaout,52) istau2,dsqrt(SLegg(2)),
     $        '~tau_2'
         write(lhaout,52) intau1,dsqrt(SNegg(1)),'~nu_tauL'

      else

         write(lhaout,52) isdl,dsqrt(SDegvd(6)),'~d_1'
         write(lhaout,52) issl,dsqrt(SDegvd(5)),'~d_2'
         write(lhaout,52) isb1,dsqrt(SDegvd(4)),'~d_3'
         write(lhaout,52) isdr,dsqrt(SDegvd(3)),'~d_4'
         write(lhaout,52) issr,dsqrt(SDegvd(2)),
     $        '~d_5'
         write(lhaout,52) isb2,dsqrt(SDegvd(1)),
     $        '~d_6'

         write(lhaout,52) isul,dsqrt(SUegvd(6)),'~u_1'
         write(lhaout,52) iscl,dsqrt(SUegvd(5)),'~u_2'
         write(lhaout,52) ist1,dsqrt(SUegvd(4)),'~u_3'
         write(lhaout,52) isur,dsqrt(SUegvd(3)),'~u_4'
         write(lhaout,52) iscr,dsqrt(SUegvd(2)),
     $        '~u_5'
         write(lhaout,52) ist2,dsqrt(SUegvd(1)),
     $        '~u_6'

         write(lhaout,52) isell,dsqrt(SLegvd(6)),'~l_1'
         write(lhaout,52) ismul,dsqrt(SLegvd(5)),'~l_2'
         write(lhaout,52) istau1,dsqrt(SLegvd(4)),
     $        '~l_3'
         write(lhaout,52) iselr,dsqrt(SLegvd(3)),
     $        '~l_4'
         write(lhaout,52) ismur,dsqrt(SLegvd(2)),
     $        '~l_5'
         write(lhaout,52) istau2,dsqrt(SLegvd(1)),
     $        '~l_6'
         
         write(lhaout,52) inel, dsqrt(SNegvd(3)),'~nu_1'
         write(lhaout,52) inmul,dsqrt(SNegvd(2)),'~nu_2'
         write(lhaout,52) intau1,dsqrt(SNegvd(1)),'~nu_3'

      endif

      write(lhaout,52) iglo,M3t,'~g'
      write(lhaout,522) in1,Neg(1),'~chi_10'
      write(lhaout,522) in2,Neg(2),'~chi_20'
      write(lhaout,522) in3,Neg(3),'~chi_30'
      write(lhaout,522) in4,Neg(4),'~chi_40'
      write(lhaout,52) ic1,ceg(1),'~chi_1+'
      write(lhaout,52) ic2,ceg(2),'~chi_2+'
      if(model.eq.'GMSB')then
      write(lhaout,52) igrav,gravitino,'~gravitino'
      endif




c ------------------------------------------------------------------- c
c The neutralino mixing matrix N and the chargino mixing matrices U,V c
c ------------------------------------------------------------------- c
      
      write(lhaout,105)
      write(lhaout,51) 'NMIX','Neutralino Mixing Matrix'
      write(lhaout,53) 1,1,ONM(1,1),'N_11'
      write(lhaout,53) 1,2,ONM(1,2),'N_12'
      write(lhaout,53) 1,3,ONM(1,3),'N_13'
      write(lhaout,53) 1,4,ONM(1,4),'N_14'
      write(lhaout,53) 2,1,ONM(2,1),'N_21'
      write(lhaout,53) 2,2,ONM(2,2),'N_22'
      write(lhaout,53) 2,3,ONM(2,3),'N_23'
      write(lhaout,53) 2,4,ONM(2,4),'N_24'
      write(lhaout,53) 3,1,ONM(3,1),'N_31'
      write(lhaout,53) 3,2,ONM(3,2),'N_32'
      write(lhaout,53) 3,3,ONM(3,3),'N_33'
      write(lhaout,53) 3,4,ONM(3,4),'N_34'
      write(lhaout,53) 4,1,ONM(4,1),'N_41'
      write(lhaout,53) 4,2,ONM(4,2),'N_42'
      write(lhaout,53) 4,3,ONM(4,3),'N_43'
      write(lhaout,53) 4,4,ONM(4,4),'N_44'


      write(lhaout,105)
      write(lhaout,51) 'UMIX','Chargino Mixing Matrix U'
      write(lhaout,53) 1,1,OCLTm(1,1),'U_11'
      write(lhaout,53) 1,2,OCLTm(1,2),'U_12'
      write(lhaout,53) 2,1,OCLTm(2,1),'U_21'
      write(lhaout,53) 2,2,OCLTm(2,2),'U_22'

      write(lhaout,105)
      write(lhaout,51) 'VMIX','Chargino Mixing Matrix V'
      write(lhaout,53) 1,1,OCRm(1,1),'V_11'
      write(lhaout,53) 1,2,OCRm(1,2),'V_12'
      write(lhaout,53) 2,1,OCRm(2,1),'V_21'
      write(lhaout,53) 2,2,OCRm(2,2),'V_22'

c ------------------------------------------ c
c The stop, sbottom and stau mixing matrices c
c ------------------------------------------ c

      if(rhn.eq.0)then

         write(lhaout,105)
         write(lhaout,51) 'STOPMIX','Stop Mixing Matrix'
         write(lhaout,53) 1,1,dcos(thetat),'STOPMIX_11'
         write(lhaout,53) 1,2,dsin(thetat),'STOPMIX_12'
         write(lhaout,53) 2,1,-dsin(thetat),'STOPMIX_21'
         write(lhaout,53) 2,2,dcos(thetat),'STOPMIX_22'

         write(lhaout,105)
         write(lhaout,51) 'SBOTMIX','Sbottom Mixing Matrix'
         write(lhaout,53) 1,1,dcos(thetab),'SBOTMIX_11'
         write(lhaout,53) 1,2,dsin(thetab),'SBOTMIX_12'
         write(lhaout,53) 2,1,-dsin(thetab),'SBOTMIX_21'
         write(lhaout,53) 2,2,dcos(thetab),'SBOTMIX_22'


         write(lhaout,105)
         write(lhaout,51) 'STAUMIX','Stau Mixing Matrix'
         write(lhaout,53) 1,1,dcos(thetatau),'STAUMIX_11'
         write(lhaout,53) 1,2,dsin(thetatau) ,'STAUMIX_12'
         write(lhaout,53) 2,1,-dsin(thetatau),'STAUMIX_21'
         write(lhaout,53) 2,2,dcos(thetatau),'STAUMIX_22'

      else 

         write(lhaout,105)
         write(lhaout,51) 'USQMIX','squark Mixing Matrix'
         write(lhaout,53) 1,1,USUvd(6,1),'USQMIX_11'
         write(lhaout,53) 1,2,USUvd(6,2),'USQMIX_12'
         write(lhaout,53) 1,3,USUvd(6,3),'USQMIX_13'
         write(lhaout,53) 1,4,USUvd(6,4),'USQMIX_14'
         write(lhaout,53) 1,5,USUvd(6,5),'USQMIX_15'
         write(lhaout,53) 1,6,USUvd(6,6),'USQMIX_16'
         write(lhaout,53) 2,1,USUvd(5,1),'USQMIX_21'
         write(lhaout,53) 2,2,USUvd(5,2),'USQMIX_22'
         write(lhaout,53) 2,3,USUvd(5,3),'USQMIX_23'
         write(lhaout,53) 2,4,USUvd(5,4),'USQMIX_24'
         write(lhaout,53) 2,5,USUvd(5,5),'USQMIX_25'
         write(lhaout,53) 2,6,USUvd(5,6),'USQMIX_26'
         write(lhaout,53) 3,1,USUvd(4,1),'USQMIX_31'
         write(lhaout,53) 3,2,USUvd(4,2),'USQMIX_32'
         write(lhaout,53) 3,3,USUvd(4,3),'USQMIX_33'
         write(lhaout,53) 3,4,USUvd(4,4),'USQMIX_34'
         write(lhaout,53) 3,5,USUvd(4,5),'USQMIX_35'
         write(lhaout,53) 3,6,USUvd(4,6),'USQMIX_36'
         write(lhaout,53) 4,1,USUvd(3,1),'USQMIX_41'
         write(lhaout,53) 4,2,USUvd(3,2),'USQMIX_42'
         write(lhaout,53) 4,3,USUvd(3,3),'USQMIX_43'
         write(lhaout,53) 4,4,USUvd(3,4),'USQMIX_44'
         write(lhaout,53) 4,5,USUvd(3,5),'USQMIX_45'
         write(lhaout,53) 4,6,USUvd(3,6),'USQMIX_46'
         write(lhaout,53) 5,1,USUvd(2,1),'USQMIX_51'
         write(lhaout,53) 5,2,USUvd(2,2),'USQMIX_52'
         write(lhaout,53) 5,3,USUvd(2,3),'USQMIX_53'
         write(lhaout,53) 5,4,USUvd(2,4),'USQMIX_54'
         write(lhaout,53) 5,5,USUvd(2,5),'USQMIX_55'
         write(lhaout,53) 5,6,USUvd(2,6),'USQMIX_56'
         write(lhaout,53) 6,1,USUvd(1,1),'USQMIX_61'
         write(lhaout,53) 6,2,USUvd(1,2),'USQMIX_62'
         write(lhaout,53) 6,3,USUvd(1,3),'USQMIX_63'
         write(lhaout,53) 6,4,USUvd(1,4),'USQMIX_64'
         write(lhaout,53) 6,5,USUvd(1,5),'USQMIX_65'
         write(lhaout,53) 6,6,USUvd(1,6),'USQMIX_66'

         write(lhaout,105)
         write(lhaout,51) 'DSQMIX','squark Mixing Matrix'
         write(lhaout,53) 1,1,USDvd(6,1),'DSQMIX_11'
         write(lhaout,53) 1,2,USDvd(6,2),'DSQMIX_12'
         write(lhaout,53) 1,3,USDvd(6,3),'DSQMIX_13'
         write(lhaout,53) 1,4,USDvd(6,4),'DSQMIX_14'
         write(lhaout,53) 1,5,USDvd(6,5),'DSQMIX_15'
         write(lhaout,53) 1,6,USDvd(6,6),'DSQMIX_16'
         write(lhaout,53) 2,1,USDvd(5,1),'DSQMIX_21'
         write(lhaout,53) 2,2,USDvd(5,2),'DSQMIX_22'
         write(lhaout,53) 2,3,USDvd(5,3),'DSQMIX_23'
         write(lhaout,53) 2,4,USDvd(5,4),'DSQMIX_24'
         write(lhaout,53) 2,5,USDvd(5,5),'DSQMIX_25'
         write(lhaout,53) 2,6,USDvd(5,6),'DSQMIX_26'
         write(lhaout,53) 3,1,USDvd(4,1),'DSQMIX_31'
         write(lhaout,53) 3,2,USDvd(4,2),'DSQMIX_32'
         write(lhaout,53) 3,3,USDvd(4,3),'DSQMIX_33'
         write(lhaout,53) 3,4,USDvd(4,4),'DSQMIX_34'
         write(lhaout,53) 3,5,USDvd(4,5),'DSQMIX_35'
         write(lhaout,53) 3,6,USDvd(4,6),'DSQMIX_36'
         write(lhaout,53) 4,1,USDvd(3,1),'DSQMIX_41'
         write(lhaout,53) 4,2,USDvd(3,2),'DSQMIX_42'
         write(lhaout,53) 4,3,USDvd(3,3),'DSQMIX_43'
         write(lhaout,53) 4,4,USDvd(3,4),'DSQMIX_44'
         write(lhaout,53) 4,5,USDvd(3,5),'DSQMIX_45'
         write(lhaout,53) 4,6,USDvd(3,6),'DSQMIX_46'
         write(lhaout,53) 5,1,USDvd(2,1),'DSQMIX_51'
         write(lhaout,53) 5,2,USDvd(2,2),'DSQMIX_52'
         write(lhaout,53) 5,3,USDvd(2,3),'DSQMIX_53'
         write(lhaout,53) 5,4,USDvd(2,4),'DSQMIX_54'
         write(lhaout,53) 5,5,USDvd(2,5),'DSQMIX_55'
         write(lhaout,53) 5,6,USDvd(2,6),'DSQMIX_56'
         write(lhaout,53) 6,1,USDvd(1,1),'DSQMIX_61'
         write(lhaout,53) 6,2,USDvd(1,2),'DSQMIX_62'
         write(lhaout,53) 6,3,USDvd(1,3),'DSQMIX_63'
         write(lhaout,53) 6,4,USDvd(1,4),'DSQMIX_64'
         write(lhaout,53) 6,5,USDvd(1,5),'DSQMIX_65'
         write(lhaout,53) 6,6,USDvd(1,6),'DSQMIX_66'


         write(lhaout,105)
         write(lhaout,51) 'SELMIX','Slepton Mixing Matrix'
         write(lhaout,53) 1,1,USLvd(6,1),'SELMIX_11'
         write(lhaout,53) 1,2,USLvd(6,2),'SELMIX_12'
         write(lhaout,53) 1,3,USLvd(6,3),'SELMIX_13'
         write(lhaout,53) 1,4,USLvd(6,4),'SELMIX_14'
         write(lhaout,53) 1,5,USLvd(6,5),'SELMIX_15'
         write(lhaout,53) 1,6,USLvd(6,6),'SELMIX_16'
         write(lhaout,53) 2,1,USLvd(5,1),'SELMIX_21'
         write(lhaout,53) 2,2,USLvd(5,2),'SELMIX_22'
         write(lhaout,53) 2,3,USLvd(5,3),'SELMIX_23'
         write(lhaout,53) 2,4,USLvd(5,4),'SELMIX_24'
         write(lhaout,53) 2,5,USLvd(5,5),'SELMIX_25'
         write(lhaout,53) 2,6,USLvd(5,6),'SELMIX_26'
         write(lhaout,53) 3,1,USLvd(4,1),'SELMIX_31'
         write(lhaout,53) 3,2,USLvd(4,2),'SELMIX_32'
         write(lhaout,53) 3,3,USLvd(4,3),'SELMIX_33'
         write(lhaout,53) 3,4,USLvd(4,4),'SELMIX_34'
         write(lhaout,53) 3,5,USLvd(4,5),'SELMIX_35'
         write(lhaout,53) 3,6,USLvd(4,6),'SELMIX_36'
         write(lhaout,53) 4,1,USLvd(3,1),'SELMIX_41'
         write(lhaout,53) 4,2,USLvd(3,2),'SELMIX_42'
         write(lhaout,53) 4,3,USLvd(3,3),'SELMIX_43'
         write(lhaout,53) 4,4,USLvd(3,4),'SELMIX_44'
         write(lhaout,53) 4,5,USLvd(3,5),'SELMIX_45'
         write(lhaout,53) 4,6,USLvd(3,6),'SELMIX_46'
         write(lhaout,53) 5,1,USLvd(2,1),'SELMIX_51'
         write(lhaout,53) 5,2,USLvd(2,2),'SELMIX_52'
         write(lhaout,53) 5,3,USLvd(2,3),'SELMIX_53'
         write(lhaout,53) 5,4,USLvd(2,4),'SELMIX_54'
         write(lhaout,53) 5,5,USLvd(2,5),'SELMIX_55'
         write(lhaout,53) 5,6,USLvd(2,6),'SELMIX_56'
         write(lhaout,53) 6,1,USLvd(1,1),'SELMIX_61'
         write(lhaout,53) 6,2,USLvd(1,2),'SELMIX_62'
         write(lhaout,53) 6,3,USLvd(1,3),'SELMIX_63'
         write(lhaout,53) 6,4,USLvd(1,4),'SELMIX_64'
         write(lhaout,53) 6,5,USLvd(1,5),'SELMIX_65'
         write(lhaout,53) 6,6,USLvd(1,6),'SELMIX_66'
         

         write(lhaout,105)
         write(lhaout,51) 'SNUMIX','Sneutrino Mixing Matrix'
         write(lhaout,53) 1,1,USNvd(3,1),'SNUMIX_11'
         write(lhaout,53) 1,2,USNvd(3,2),'SNUMIX_12'
         write(lhaout,53) 1,3,USNvd(3,3),'SNUMIX_13'
         write(lhaout,53) 2,1,USNvd(2,1),'SNUMIX_21'
         write(lhaout,53) 2,2,USNvd(2,2),'SNUMIX_22'
         write(lhaout,53) 2,3,USNvd(2,3),'SNUMIX_23'
         write(lhaout,53) 3,1,USNvd(1,1),'SNUMIX_31'
         write(lhaout,53) 3,2,USNvd(1,2),'SNUMIX_32'
         write(lhaout,53) 3,3,USNvd(1,3),'SNUMIX_33'

      endif

c ------------------------------------------------------------------- c
c The angle alpha in the Higgs sector and the Higgs mixing parameters c
c ------------------------------------------------------------------- c

      write(lhaout,105)
      write(lhaout,51) 'ALPHA','Higgs mixing'
      write(lhaout,109) alphahiggs,'Mixing angle in the neutral', 
     $     ' Higgs boson sector'

      write(lhaout,105)
      write(lhaout,54) 'HMIX Q=',scale,'DRbar Higgs Parameters'
      write(lhaout,55) 1,sgnmu*murge,'mu(Q)'
      write(lhaout,55) 2,vev2n/vev1n,'tanbeta(Q)'
      write(lhaout,55) 3,dsqrt(vev2n**2.d0 + vev1n**2.d0),'vev(Q)'
      write(lhaout,55) 4, mAm3 ,'MA^2(Q)'


c ------------------- c
c The gauge couplings c
c ------------------- c

      write(lhaout,105)      
      write(lhaout,54) 'GAUGE Q=',scale,'The gauge couplings'
      write(lhaout,55) 1,dsqrt(alph1*16.d0*pi*pi)*dsqrt(3.d0/5.d0),
     $     'gprime(Q) DRbar'
      write(lhaout,55) 2,dsqrt(alph2*16.d0*pi*pi),'g(Q) DRbar'
      write(lhaout,55) 3,dsqrt(alph3*16.d0*pi*pi),'g_3(Q) DRbar'

c ------------------------------------- c
c The trilinear couplings Au, Ad and Ae c
c ------------------------------------- c
      scalesave=scale

      if(rhn.eq.0)then
         
         write(lhaout,105)
         write(lhaout,54) 'AU Q=',scale,'The trilinear couplings'
         write(lhaout,53) 1,1,AURG(1,1), 'A_u(Q) DRbar'
         write(lhaout,53) 2,2,AURG(2,2), 'A_c(Q) DRbar'
         write(lhaout,53) 3,3,AURG(3,3),'A_t(Q) DRbar'

         write(lhaout,105)
         write(lhaout,54) 'AD Q=',scale,'The trilinear couplings'
         write(lhaout,53) 1,1,ADRG(1,1),'A_d(Q) DRbar'
         write(lhaout,53) 2,2,ADRG(2,2),'A_s(Q) DRbar'
         write(lhaout,53) 3,3,ADRG(3,3),'A_b(Q) DRbar'
 
         write(lhaout,105)
         write(lhaout,54) 'AE Q=',scale,'The trilinear couplings'
         write(lhaout,53) 1,1,AERG(1,1),'A_e(Q) DRbar'
         write(lhaout,53) 2,2,AERG(2,2),'A_mu(Q) DRbar'
         write(lhaout,53) 3,3,AERG(3,3),'A_tau(Q) DRbar'

      else
         
         write(lhaout,105)
         write(lhaout,54) 'TU Q=',scale,'The trilinear couplings'
         write(lhaout,53) 1,1,TURG(1,1),'TURG_11'
         write(lhaout,53) 1,2,TURG(1,2),'TURG_12'
         write(lhaout,53) 1,3,TURG(1,3),'TURG_13'
         write(lhaout,53) 2,1,TURG(2,1),'TURG_21'
         write(lhaout,53) 2,2,TURG(2,2),'TURG_22'
         write(lhaout,53) 2,3,TURG(2,3),'TURG_23'
         write(lhaout,53) 3,1,TURG(3,1),'TURG_31'
         write(lhaout,53) 3,2,TURG(3,2),'TURG_32'
         write(lhaout,53) 3,3,TURG(3,3),'TURG_33'


         write(lhaout,105)
         write(lhaout,54) 'TD Q=',scale,'The trilinear couplings'
         write(lhaout,53) 1,1,TDRG(1,1),'TDRG_11'
         write(lhaout,53) 1,2,TDRG(1,2),'TDRG_12'
         write(lhaout,53) 1,3,TDRG(1,3),'TDRG_13'
         write(lhaout,53) 2,1,TDRG(2,1),'TDRG_21'
         write(lhaout,53) 2,2,TDRG(2,2),'TDRG_22'
         write(lhaout,53) 2,3,TDRG(2,3),'TDRG_23'
         write(lhaout,53) 3,1,TDRG(3,1),'TDRG_31'
         write(lhaout,53) 3,2,TDRG(3,2),'TDRG_32'
         write(lhaout,53) 3,3,TDRG(3,3),'TDRG_33'

         write(lhaout,105)
         write(lhaout,54) 'TE Q=',scale,'The trilinear couplings'
         write(lhaout,53) 1,1,TERG(1,1),'TERG_11'
         write(lhaout,53) 1,2,TERG(1,2),'TERG_12'
         write(lhaout,53) 1,3,TERG(1,3),'TERG_13'
         write(lhaout,53) 2,1,TERG(2,1),'TERG_21'
         write(lhaout,53) 2,2,TERG(2,2),'TERG_22'
         write(lhaout,53) 2,3,TERG(2,3),'TERG_23'
         write(lhaout,53) 3,1,TERG(3,1),'TERG_31'
         write(lhaout,53) 3,2,TERG(3,2),'TERG_32'
         write(lhaout,53) 3,3,TERG(3,3),'TERG_33'

      endif
c ---------------------------------- c
c The Yukawa couplings Yu, Yd and Ye c
c ---------------------------------- c

      if(rhn.eq.0.d0)then
         
         write(lhaout,105)
         write(lhaout,54) 'YU Q=',scalesave,'The top Yukawa coupling'
         write(lhaout,53) 3,3, yuRG(3,3),'y_top(Q) DRbar'
c     
         write(lhaout,105)
         write(lhaout,54) 'YD Q=',scalesave,'The down Yukawa coupling'
         write(lhaout,53) 3,3, ydRG(3,3),'y_b(Q) DRbar'
c     
         write(lhaout,105)
         write(lhaout,54) 'YE Q=',scalesave,'The tau Yukawa coupling'
         write(lhaout,53) 3,3, yeRG(3,3),'y_tau(Q) DRbar'
         
      else

         write(lhaout,105)
         write(lhaout,54) 'YU Q=',scalesave,'The top Yukawa coupling'
         write(lhaout,53) 1,1, yuRG(1,1),'YU_11'
         write(lhaout,53) 1,2, yuRG(1,2),'YU_12'
         write(lhaout,53) 1,3, yuRG(1,3),'YU_13'
         write(lhaout,53) 2,1, yuRG(2,1),'YU_21'
         write(lhaout,53) 2,2, yuRG(2,2),'YU_22'
         write(lhaout,53) 2,3, yuRG(2,3),'YU_23'
         write(lhaout,53) 3,1, yuRG(3,1),'YU_31'
         write(lhaout,53) 3,2, yuRG(3,2),'YU_32'
         write(lhaout,53) 3,3, yuRG(3,3),'YU_33'
c     
         write(lhaout,105)
         write(lhaout,54) 'YD Q=',scalesave,'The down Yukawa coupling'
         write(lhaout,53) 1,1, ydRG(1,1),'YD_11'
         write(lhaout,53) 1,2, ydRG(1,2),'YD_12'
         write(lhaout,53) 1,3, ydRG(1,3),'YD_13'
         write(lhaout,53) 2,1, ydRG(2,1),'YD_21'
         write(lhaout,53) 2,2, ydRG(2,2),'YD_22'
         write(lhaout,53) 2,3, ydRG(2,3),'YD_23'
         write(lhaout,53) 3,1, ydRG(3,1),'YD_31'
         write(lhaout,53) 3,2, ydRG(3,2),'YD_32'
         write(lhaout,53) 3,3, ydRG(3,3),'YD_33'

c     
         write(lhaout,105)
         write(lhaout,54) 'YE Q=',scalesave,'The tau Yukawa coupling'
         write(lhaout,53) 1,1, yeRG(1,1),'YE_11'
         write(lhaout,53) 1,2, yeRG(1,2),'YE_12'
         write(lhaout,53) 1,3, yeRG(1,3),'YE_13'
         write(lhaout,53) 2,1, yeRG(2,1),'YE_21'
         write(lhaout,53) 2,2, yeRG(2,2),'YE_22'
         write(lhaout,53) 2,3, yeRG(2,3),'YE_23'
         write(lhaout,53) 3,1, yeRG(3,1),'YE_31'
         write(lhaout,53) 3,2, yeRG(3,2),'YE_32'
         write(lhaout,53) 3,3, yeRG(3,3),'YE_33'
         
      endif

 
c ----------------------------- c
c The soft SUSY breaking masses c
c ----------------------------- c

      write(lhaout,105)
      write(lhaout,108) 'MSOFT Q=',scale,'soft SUSY breaking masses at', 
     $     ' scale Q'
      write(lhaout,52) 1,M1tz,'M_1'
      write(lhaout,52) 2,M2tz,'M_2'
      write(lhaout,52) 3,M3tz,'M_3'
      write(lhaout,52) 21,mh2mz,'M^2_Hd'
      write(lhaout,52) 22,mh1mz,'M^2_Hu'
      write(lhaout,52) 31,dsqrt(mSLRG(1,1)),'M_eL'
      write(lhaout,52) 32,dsqrt(mSLRG(2,2)),'M_muL'
      write(lhaout,52) 33,dsqrt(mSLRG(3,3)),'M_tauL'
      write(lhaout,52) 34,dsqrt(mSERG(1,1)),'M_eR'
      write(lhaout,52) 35,dsqrt(mSERG(2,2)),'M_muR'
      write(lhaout,52) 36,dsqrt(mSERG(3,3)),'M_tauR'
      write(lhaout,52) 41,dsqrt(mSQRG(1,1)),'M_q1L'
      write(lhaout,52) 42,dsqrt(mSQRG(2,2)),'M_q2L'
      write(lhaout,52) 43,dsqrt(mSQRG(3,3)),'M_q3L'
      write(lhaout,52) 44,dsqrt(mSURG(1,1)),'M_uR'
      write(lhaout,52) 45,dsqrt(mSURG(2,2)),'M_cR'
      write(lhaout,52) 46,dsqrt(mSURG(3,3)),'M_tR'
      write(lhaout,52) 47,dsqrt(mSDRG(1,1)),'M_dR'
      write(lhaout,52) 48,dsqrt(mSDRG(2,2)),'M_sR'
      write(lhaout,52) 49,dsqrt(mSDRG(3,3)),'M_bR'

      if(rhn.eq.1)then
         write(lhaout,105)
         write(lhaout,108) 'MSQ2 Q=',scale,'M^2_Q soft SUSY breaking', 
     $        ' masses'
         write(lhaout,53) 1,1,mSQRG(1,1),'mSQRG_11'
         write(lhaout,53) 1,2,mSQRG(1,2),'mSQRG_12'
         write(lhaout,53) 1,3,mSQRG(1,3),'mSQRG_13'
         write(lhaout,53) 2,1,mSQRG(2,1),'mSQRG_21'
         write(lhaout,53) 2,2,mSQRG(2,2),'mSQRG_22'
         write(lhaout,53) 2,3,mSQRG(2,3),'mSQRG_23'
         write(lhaout,53) 3,1,mSQRG(3,1),'mSQRG_31'
         write(lhaout,53) 3,2,mSQRG(3,2),'mSQRG_32'
         write(lhaout,53) 3,3,mSQRG(3,3),'mSQRG_33'

         write(lhaout,105)
         write(lhaout,108) 'MSU2 Q=',scale,'M^2_U soft SUSY breaking', 
     $        ' masses'
         write(lhaout,53) 1,1,mSURG(1,1),'mSURG_11'
         write(lhaout,53) 1,2,mSURG(1,2),'mSURG_12'
         write(lhaout,53) 1,3,mSURG(1,3),'mSURG_13'
         write(lhaout,53) 2,1,mSURG(2,1),'mSURG_21'
         write(lhaout,53) 2,2,mSURG(2,2),'mSURG_22'
         write(lhaout,53) 2,3,mSURG(2,3),'mSURG_23'
         write(lhaout,53) 3,1,mSURG(3,1),'mSURG_31'
         write(lhaout,53) 3,2,mSURG(3,2),'mSURG_32'
         write(lhaout,53) 3,3,mSURG(3,3),'mSURG_33'

         write(lhaout,105)
         write(lhaout,108) 'MSD2 Q=',scale,'M^2_D soft SUSY breaking', 
     $        ' masses'
         write(lhaout,53) 1,1,mSDRG(1,1),'mSDRG_11'
         write(lhaout,53) 1,2,mSDRG(1,2),'mSDRG_12'
         write(lhaout,53) 1,3,mSDRG(1,3),'mSDRG_13'
         write(lhaout,53) 2,1,mSDRG(2,1),'mSDRG_21'
         write(lhaout,53) 2,2,mSDRG(2,2),'mSDRG_22'
         write(lhaout,53) 2,3,mSDRG(2,3),'mSDRG_23'
         write(lhaout,53) 3,1,mSDRG(3,1),'mSDRG_31'
         write(lhaout,53) 3,2,mSDRG(3,2),'mSDRG_32'
         write(lhaout,53) 3,3,mSDRG(3,3),'mSDRG_33'

         write(lhaout,105)
         write(lhaout,108) 'MSL2 Q=',scale,'M^2_L soft SUSY breaking',
     $        ' masses'
         write(lhaout,53) 1,1,mSLRG(1,1),'mSLRG_11'
         write(lhaout,53) 1,2,mSLRG(1,2),'mSLRG_12'
         write(lhaout,53) 1,3,mSLRG(1,3),'mSLRG_13'
         write(lhaout,53) 2,1,mSLRG(2,1),'mSLRG_21'
         write(lhaout,53) 2,2,mSLRG(2,2),'mSLRG_22'
         write(lhaout,53) 2,3,mSLRG(2,3),'mSLRG_23'
         write(lhaout,53) 3,1,mSLRG(3,1),'mSLRG_31'
         write(lhaout,53) 3,2,mSLRG(3,2),'mSLRG_32'
         write(lhaout,53) 3,3,mSLRG(3,3),'mSLRG_33'

         write(lhaout,105)
         write(lhaout,108) 'MSE2 Q=',scale,'M^2_E soft SUSY breaking',
     $        ' masses'
         write(lhaout,53) 1,1,mSERG(1,1),'mSERG_11'
         write(lhaout,53) 1,2,mSERG(1,2),'mSERG_12'
         write(lhaout,53) 1,3,mSERG(1,3),'mSERG_13'
         write(lhaout,53) 2,1,mSERG(2,1),'mSERG_21'
         write(lhaout,53) 2,2,mSERG(2,2),'mSERG_22'
         write(lhaout,53) 2,3,mSERG(2,3),'mSERG_23'
         write(lhaout,53) 3,1,mSERG(3,1),'mSERG_31'
         write(lhaout,53) 3,2,mSERG(3,2),'mSERG_32'
         write(lhaout,53) 3,3,mSERG(3,3),'mSERG_33'

      endif

c  low energy parameters:
c
      write(lhaout,105)
      write(lhaout,51) 'SuSeFLAVLOWENERGY ',' PARAMETERS '
      write(lhaout,52) 1,Drho,'Delta rho parameter'
      write(lhaout,52) 2,gm2,'g_mu - 2'
      write(lhaout,52) 3,Bbsg,'Br(b -> s gamma)'
      write(lhaout,52) 4,Btmug,'Br(tau -> mu gamma)'
      write(lhaout,52) 5,Bteg,'Br(tau -> e gamma)'
      write(lhaout,52) 6,Bmeg,'Br(mu -> e gamma)'
      write(lhaout,52) 7,brtau3mu,'Br(tau -> mu mu mu)'
      write(lhaout,52) 8,brtau3e,'Br(tau -> e e e)'
      write(lhaout,52) 9,Brmu3e,'Br(mu -> e e e)'


c  fine-tuning parameters :
c
      write(lhaout,105)
      write(lhaout,51) 'FINETUNE','Fine Tuning'
      write(lhaout,52) 1,ftmz,'delta mZ^2/mZ^2 (mu^2)'
      write(lhaout,52) 2,ftmt,'delta mt/mt (mu^2)'


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

      close(lhaout)
c -------------------------------------------------------------------------
      end subroutine
!=====================================================================================================

      subroutine readslha(input)

      character line1*6,line2*60                          ! the line buffer
      integer, parameter:: linemax=50              ! max on a line
      real::val(linemax)
      integer i,k,num,modnum,smnum(10),minparnum(10),counter_sminp
      integer counter_minpar,counter_algopar,algonum(20)
      double precision sminp(10), minpar(10), algopar(20)

      character (len=80) :: line3
      integer i_test,i_mod,qmix

      common/slha_modsel/num, modnum
      common/slha_counters/counter_sminp,counter_minpar,counter_algopar
      common/slha_num/smnum,minparnum,algonum
      common/slha_variables/sminp,minpar,algopar
      common/quarkmix/ qmix

!--------------------------
      do k=1,linemax,1
         val(k)= 0.d0
      end do

      rewind(input)

       domodsel: do
                                                 ! read a line into memory
      read (input, "(a6,a60)", end=999) line1,line2
!      print*, "line number and value ", line1,line2 

!      print*,line11,line12

!      print*,line11,line12

!---------Block MODSEL------------

      if(line1(1:5).eq.'Block'.OR.line1(1:5).eq.'BLOCK')then

         if(line2(1:6).eq.'MODSEL')then

!            read(input,*) num, modnum

            do 
               read(input,*) line3
               
!               print*,"line3 = ", line3
               
               if (line3(1:1).eq."#") cycle

               backspace(input) ! rewind

               if ((line3(1:1).eq."B").or.
     $              (line3(1:1).eq."b")) exit

               read(input,*) inum,imodnum 

!               print*,"itest, imod =", inum,imodnum,trim(line3)
               
               if(inum.eq.1)then

                  num = inum
                  modnum = imodnum
                  
               elseif(inum.eq.6)then

                  if(imodnum.eq.0)then
                     qmix = 0
                  elseif(imodnum.eq.1)then
                     qmix = 1
                  else
                     qmix = 5
                  endif

               endif

            enddo

         endif

      endif

      enddo domodsel

 999  continue

!-------------Block sminputs---

      counter_sminp = 0

      rewind(input)

      dosminputs:  do
                                                  ! read a line into memory
       read (input, '(a6,a60)', end=1000,err=1000) line1,line2


      if(line1(1:5).eq.'Block'.OR.line1(1:5).eq.'BLOCK')then

         if(line2(1:8).eq.'SMINPUTS')then

            do i = 1,10

            read(input,*,end=110,err=110) smnum(i), sminp(i)
 !           print*,'smnum and sminp ', smnum(i), sminp(i)
            counter_sminp = counter_sminp + 1
            
            enddo
            endif

      endif

 110  continue

      enddo dosminputs
 1000 continue


!------------Block minpar

      counter_minpar = 0

      rewind(input)
        dominpar: do
                                                 ! read a line into memory
      read (input, '(a6,a60)', end=1010,err=1010) line1,line2


      if(line1(1:5).eq.'Block'.OR.line1(1:5).eq.'BLOCK')then

         if(line2(1:6).eq.'MINPAR')then

            do i = 1,10

            read(input,*,end=120,err=120) minparnum(i), minpar(i)
 !           print*,'minnum and minpar ', minparnum(i), minpar(i)
            counter_minpar = counter_minpar + 1
            
            enddo
            endif

      endif

 120  continue

      enddo  dominpar
 1010  continue

!---------Block SUSEFLAV

      counter_algopar = 0

      rewind(input)
       doalgo: do
                                                 ! read a line into memory
      read (input, '(a6,a60)', end=1020,err=1020) line1,line2


      if(line1(1:5).eq.'Block'.OR.line1(1:5).eq.'BLOCK')then

         if(line2(1:8).eq.'SUSEFLAV')then

            do i = 1,20

            read(input,*,end=130,err=130) algonum(i), algopar(i)
 !           print*,'algonum and algopar ', algonum(i), algopar(i)
            counter_algopar = counter_algopar + 1
            
            enddo
            endif

      endif

 130  continue

      enddo doalgo
 1020  continue


      end subroutine


!=================================================================================================
      subroutine uppercase(name)
      implicit none
      character(len=80), intent(inout) :: name
      integer :: len=80, i1

      do i1=1,len
         if (name(i1:i1).eq."a") name(i1:i1) = "A"
         if (name(i1:i1).eq."b") name(i1:i1) = "B"
         if (name(i1:i1).eq."c") name(i1:i1) = "C"
         if (name(i1:i1).eq."d") name(i1:i1) = "C"
         if (name(i1:i1).eq."e") name(i1:i1) = "E"
         if (name(i1:i1).eq."f") name(i1:i1) = "F"
         if (name(i1:i1).eq."g") name(i1:i1) = "G"
         if (name(i1:i1).eq."h") name(i1:i1) = "H"
         if (name(i1:i1).eq."i") name(i1:i1) = "I"
         if (name(i1:i1).eq."j") name(i1:i1) = "J"
         if (name(i1:i1).eq."k") name(i1:i1) = "K"
         if (name(i1:i1).eq."l") name(i1:i1) = "L"
         if (name(i1:i1).eq."m") name(i1:i1) = "M"
         if (name(i1:i1).eq."n") name(i1:i1) = "N"
         if (name(i1:i1).eq."o") name(i1:i1) = "O"
         if (name(i1:i1).eq."p") name(i1:i1) = "P"
         if (name(i1:i1).eq."q") name(i1:i1) = "Q"
         if (name(i1:i1).eq."r") name(i1:i1) = "R"
         if (name(i1:i1).eq."s") name(i1:i1) = "S"
         if (name(i1:i1).eq."t") name(i1:i1) = "T"
         if (name(i1:i1).eq."u") name(i1:i1) = "U"
         if (name(i1:i1).eq."v") name(i1:i1) = "V"
         if (name(i1:i1).eq."w") name(i1:i1) = "W"
         if (name(i1:i1).eq."x") name(i1:i1) = "X"
         if (name(i1:i1).eq."y") name(i1:i1) = "Y"
         if (name(i1:i1).eq."z") name(i1:i1) = "Z"
      end do

      end subroutine uppercase
!==============================================================================
