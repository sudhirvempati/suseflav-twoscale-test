! program for Goto and Nihei paper(9808255)
! Doing patch on 31/01/2020  
! 
      subroutine Proton_goto(tanbeta,msusy,e1,yukgut,mSQRG,mSURG,
     $           yuRG,AURG,mSDRG,ydRG,ADRG,mSLRG,mSERG,yeRG,AERG,Mchar,
     $           MNeut,M1tz,M2tz,mur,vev1,vev2,M3t,check,halfppil,
     $           halfpetal,halfpkl,halfpipnu,halfpkpnu)

      implicit none

      integer i,j,k,l,check,n0,nok,nbad,A,Abar,I0,M,N,m1
      integer i1,j1,n1,p

      CHARACTER*100 flags

      DOUBLE PRECISION Hfunc

      DOUBLE PRECISION mufrge,Mp,Mpion0,Mpionp,MEta0,Mk0,Mkp
      DOUBLE PRECISION PiRL0,PiLL0,PiRLp,PiLLp,EtaRL,EtaLL,KRL0,KLL0
      DOUBLE PRECISION KsRLp,KsLLp,KdRLp,KdLLp,KdsRLp,KdsLLp,beta
      DOUBLE PRECISION Mpl,tanbeta,MZ,MX,MZpole,msusy,QCDscale,MHC
      DOUBLE PRECISION Mtpole,Mbpole,tpl,tq0,tZ,mtscale,mbscale,Mgut
      DOUBLE PRECISION fuGUT(3,3),fdGUT(3,3),VCKMGUT(3,3),VCKMGUTT(3,3) 
      DOUBLE PRECISION yugut(3,3),ydgut(3,3),yegut(3,3),delta(3,3)
      DOUBLE PRECISION C5L(3,3,3,3),C5R(3,3,3,3),yy_sm(31),pyye(31)            
      DOUBLE PRECISION pyyb(31),pyyg(518),pyy2(31),alph3gut,alph2gut
      DOUBLE PRECISION alph1gut,vev1gut,vev2gut,x2,x1,h1,hmin,eps,pi
      DOUBLE PRECISION VNLd(6,4,3),VNRd(6,4,3),VNLu(6,4,3),VNRu(6,4,3)
      DOUBLE PRECISION VNLl(6,4,3),VNRl(6,4,3), VNLnu(3,4,3)
      DOUBLE PRECISION VGLD(6,3),VGRd(6,3),VGLu(6,3),VGRu(6,3)
      DOUBLE PRECISION VCLd(6,2,3),VCRd(6,2,3),VCLu(6,2,3),VCRu(6,2,3)
      DOUBLE PRECISION VCLl(3,2,3),VCRl(3,2,3),VCLnu(6,2,3)

      DOUBLE PRECISION UD(6,6),UDT(6,6),UU(6,6),UUT(6,6),UN(4,4)
      DOUBLE PRECISION UNT(4,4),Uplus(2,2),UplusT(2,2)
      DOUBLE PRECISION mu(3),MW,sb,md(3),cb,Uminus(2,2),UminusT(2,2)
      DOUBLE PRECISION UNeuT(3,3),UNeu(3,3),UL(6,6),ULT(6,6),ml(3),tw
      DOUBLE PRECISION CudulL(6,6,3,3),CuudlL(6,6,3,3),CdluuL(6,6,3,3)
      DOUBLE PRECISION CuudlR(6,6,3,3),CuddnuL(6,6,3,3)
      DOUBLE PRECISION CddunuL(6,6,3,3),CuludL(6,6,3,3)
      DOUBLE PRECISION CudulR(6,6,3,3),CuludR(6,6,3,3),CdluuR(6,6,3,3)
      DOUBLE PRECISION CdnuudL(6,3,3,3),CunuddL(6,3,3,3) 
      DOUBLE PRECISION xG(6),uG(6),xC(6,2),uC(6,2),wC(6,2),nuN(6,4)
      DOUBLE PRECISION yN(6,4),zN(6,4),zC(3,2),wN(3,4)
      DOUBLE PRECISION msdown(6),MG,msup(6),mslepton(6)
      DOUBLE PRECISION MCh(2),MN(4),msnu(3),msusyold
      Double Precision :: Vfd(3,3),fuV(3,3),VTfd(3,3)
 

      DOUBLE PRECISION CLLudul(3,3),CLLudulG(3,3),CLLudulchipm(3,3)
      DOUBLE PRECISION CLLudulchi0(3,3),CRLudul(3,3),CRLudulG(3,3)
      DOUBLE PRECISION CRLudulchipm(3,3),CRLudulchi0(3,3),CLRudul(3,3)
      DOUBLE PRECISION CLRudulG(3,3),CLRudulchipm(3,3)
      DOUBLE PRECISION CRRudul(3,3),CRRudulG(3,3),CRRudulchipm(3,3)
      DOUBLE PRECISION CRRudulchi0(3,3),CLRudulchi0(3,3)
      DOUBLE PRECISION CLLuddnu(3,3,3),CLLuddnuG(3,3,3)
      DOUBLE PRECISION CLLuddnuchipm(3,3,3),CLLuddnuchi0(3,3,3)
      DOUBLE PRECISION CRLuddnu(3,3,3),CRLuddnuG(3,3,3)
      DOUBLE PRECISION CRLuddnuchipm(3,3,3),CRLuddnuchi0(3,3,3)
      DOUBLE PRECISION CRLddunu(3,3,3),CRLddunuG(3,3,3)
      DOUBLE PRECISION CRLddunuchi0(3,3,3),g3,g2,g1
      DOUBLE PRECISION mneu,Fvalue,Dvalue, Fpi,VCKMST(3,3)
      DOUBLE PRECISION alphasmsusy,alphasmt,alphasmb,alphas2GeV
      DOUBLE PRECISION ALmusuymt0,ALmusuymt,ALmtmb0,ALmtmb,ALmb2GeV0
      DOUBLE PRECISION ALmb2GeV,ALLLdelta0,ALRLdelta,alpha2msusy
      DOUBLE PRECISION CLLudul2(3,3),CRLudul2(3,3),CLRudul2(3,3)
      DOUBLE PRECISION CRRudul2(3,3),CLLuddnu2(3,3,3),CRLuddnu2(3,3,3)
      DOUBLE PRECISION CRLddunu2(3,3,3),tgut,VTYUV(3,3),VTfu(3,3)

      DOUBLE PRECISION alp,bep,ALppil(3),ARppil(3),decayppil(3)
      DOUBLE PRECISION ALpetal(3),ARpetal(3),decaypetal(3),MB
      DOUBLE PRECISION ALpkl(3),ARpkl(3),decaypkl(3),ALppipnu(3)
      DOUBLE PRECISION decayppipnu(3),ALpkpnu(3),decaypkpnu(3)
      double precision halfppil(3),halfpetal(3),halfpkl(3),halfpipnu(3)
      double precision halfpkpnu(3),alphamsbar3,alphamsbar2
!--------------------------------------------------------------------
!-----Input coming from SUSeFLAV variable
!--------------------------------------------------------------------
      DOUBLE PRECISION e1,yukgut(126),fudiag(3),fddiag(3)
      DOUBLE PRECISION UQu(3,3),UTU(3,3),UQd(3,3),UTD(3,3),UQdT(3,3)
      DOUBLE PRECISION cos2beta,sinsqtw,stw,ctw,sgnmu,M3t,LOT(6,6)
      DOUBLE PRECISION mSQRG(3,3),AURG(3,3),vev2,vev1,mur,mSURG(3,3)
      DOUBLE PRECISION ADRG(3,3),mSDRG(3,3),mSLRG(3,3),AERG(3,3)
      DOUBLE PRECISION mSERG(3,3),M1tz,M2tz,MChar(2,2),MNeut(4,4)
      DOUBLE PRECISION mSN(3,3),mSL(6,6),mSD(6,6),mSU(6,6)
      DOUBLE PRECISION fu(3),fd(3),fe(3),UQus(3,3),UTUS(3,3),UQds(3,3)
      DOUBLE PRECISION UTDS(3,3),ULS(3,3),UTES(3,3),UQdTS(3,3)
      DOUBLE PRECISION VCKMS(3,3),mSUT(6,6),mSDT(6,6),msupsq(6)        
      DOUBLE PRECISION msdownsq(6),msleptonsq(6),Lo(6,6),msnusq(3)
      DOUBLE PRECISION usnu(3,3),Neuto(4,4),MCeg(2),OLC(2,2),ORC(2,2)
      DOUBLE PRECISION yuRG(3,3),ydRG(3,3),yeRG(3,3),tautotal,etotal

      DOUBLE PRECISION CRLddunu123,CRLddunu121,CRLddunuG123,
     $  CRLddunuG121,CRLuddnu123,CRLuddnu121,CRLuddnuG123,CRLuddnuG121,
     $  CRLuddnuchipm123,CRLuddnuchipm121,CLLuddnu123,CLLuddnu121,
     $  CLLuddnuG123,CLLuddnuG121,CLLuddnuchipm123,CLLuddnuchipm121,
     $ CRLddunuchi0123,CRLddunuchi0121,CRLuddnuchi0123,CRLuddnuchi0121,
     $  CLLuddnuchi0123,CLLuddnuchi0121,CRLuddnu213,CRLuddnu211,
     $  CRLuddnuG213,CRLuddnuG211,CRLuddnuchipm213,CRLuddnuchipm211,
     $  CLLuddnu213,CLLuddnu211,CLLuddnuG213,CLLuddnuG211,
     $  CLLuddnuchipm213,CLLuddnuchipm211,CRLuddnuchi0211,
     $  CRLuddnuchi0213,CLLuddnuchi0213,CLLuddnuchi0211

      DOUBLE PRECISION ALgut,ALRLdelta6d,CLRudul26d(3,3),
     $            CRLudul26d(3,3),CLRudul6d(3,3),CRLudul6d(3,3),
     $          CRLuddnu26d(3,3,3),CRLuddnu6d(3,3,3),ALgutmsusy
      DOUBLE PRECISION ALppil6d(3),ARppil6d(3),decayppil6d(3)
      DOUBLE PRECISION ALpetal6d(3),ARpetal6d(3),decaypetal6d(3)
      DOUBLE PRECISION ALpkl6d(3),ARpkl6d(3),decaypkpnu6d(3)
      DOUBLE PRECISION ALpkpnu6d(3),halfpipnu6d(3),halfpkpnu6d(3)
      DOUBLE PRECISION decaypkl6d(3),ALppipnu6d(3),decayppipnu6d(3)
      DOUBLE PRECISION halfppil6d(3),halfpetal6d(3),halfpkl6d(3)
      DOUBLE PRECISION halfppiepure6d,decayppiepure6d,halfppimupure6d
      DOUBLE PRECISION decayppimupure6d,halfppipnuepure6d,
     $  decayppipnuepure6d,halfknuepure6d,decaypknuepure6d,
     $  halfpknumupure6d,decaypknumupure6d,halfpknupure6d
      DOUBLE PRECISION decaypknupure6d,halfpipnufull
      DOUBLE PRECISION halfpk0epure6d,decaypk0epure6d,decaypk0mupure6d,
     $   halfpk0mupure6d,halfpipnu6dfull,halfpkpnu6dfull,halfpkpnufull
      DOUBLE PRECISION CLLddunuchi0(3,3,3),CLLddunu(3,3,3)
      DOUBLE PRECISION CLLddunuG(3,3,3),CLLddunu2(3,3,3)


      double precision ALppil1,ARppil1,CRLudul11,CRLudulG11,
     $  CRLudulchipm11,CRLudulchi011,CLLudul11,CLLudulG11,
     $  CLLudulchipm11,CLLudulchi011,CRRudul11,CRRudulG11,
     $  CRRudulchipm11,CRRudulchi011,CLRudul11,CLRudulG11,
     $  CLRudulchipm11,CLRudulchi011,alpha2mz

      common/pgoespi0e_plus/ALppil1,ARppil1,CRLudul11,CRLudulG11,
     $  CRLudulchipm11,CRLudulchi011,CLLudul11,CLLudulG11,
     $  CLLudulchipm11,CLLudulchi011,CRRudul11,CRRudulG11,
     $  CRRudulchipm11,CRRudulchi011,CLRudul11,CLRudulG11,
     $  CLRudulchipm11,CLRudulchi011

      common/gluino_chargino/CRLddunu123,CRLddunu121,CRLddunuG123,
     $  CRLddunuG121,CRLuddnu123,CRLuddnu121,CRLuddnuG123,CRLuddnuG121,
     $  CRLuddnuchipm123,CRLuddnuchipm121,CLLuddnu123,CLLuddnu121,
     $  CLLuddnuG123,CLLuddnuG121,CLLuddnuchipm123,CLLuddnuchipm121,
     $ CRLddunuchi0123,CRLddunuchi0121,CRLuddnuchi0123,CRLuddnuchi0121,
     $  CLLuddnuchi0123,CLLuddnuchi0121,CRLuddnu213,CRLuddnu211,
     $  CRLuddnuG213,CRLuddnuG211,CRLuddnuchipm213,CRLuddnuchipm211,
     $  CLLuddnu213,CLLuddnu211,CLLuddnuG213,CLLuddnuG211,
     $  CLLuddnuchipm213,CLLuddnuchipm211,CRLuddnuchi0211,
     $  CRLuddnuchi0213,CLLuddnuchi0213,CLLuddnuchi0211,tautotal,etotal

      common/protondecay_6d/halfppil6d,halfpetal6d,halfpkl6d,
     $   halfpipnu6d,halfpkpnu6d,
     $  halfpipnu6dfull,halfpkpnu6dfull,halfppiepure6d,halfppimupure6d,
     $  halfppipnuepure6d,halfknuepure6d,halfpknumupure6d,
     $  halfpknupure6d,halfpk0epure6d,halfpk0mupure6d,halfpipnufull,
     $  halfpkpnufull

      external RK4ROUTINE,QMSRK4,smrge,matmult,dag2
      external smrgemt,protonrge,SVD,dag,dag4,dag6, mat3prod
      external CEigensystem
!----------------------------------------------------------------------
!hadronic matrix elements value by lattice simulation (from arXiv:1705.01338)
! and hadron masses (from pdg(check))
! 
!----------------------------------------------------------------------
      pi= 4.d0*datan(1.d0)
!      tanbeta = 5.d0

!      Write(*,*)"proton",mufrge
!      Mp = 0.9382d0
!      Mpion0 = 0.1349d0      
!      Mpionp = 0.1395d0
!      MEta0 = 0.5478d0
!      Mk0 = 0.4976d0
!      Mkp = 0.4936d0
!      mneu = 0.939d0
!      Fvalue =0.48
!      Dvalue = 0.76
!      Fpi = 0.130d0
!--(In GeV**2)
      PiRL0 = -0.131
      PiLL0 =  0.134
      PiRLp = -0.186
      PiLLp =  0.189
      EtaRL =  0.006
      EtaLL =  0.113
      KRL0  =  0.103
      KLL0  =  0.057
 !--Kaon +
      KsRLp = -0.049
      KsLLp =  0.041
      KdRLp = -0.134 
      kdLLp =  0.139
      KdsRLp = -0.054
      KdsLLp = -0.098
!----------------------------------------------------------------------      pi = 4.d0*datan(1.d0)
      beta = datan(tanbeta)
      Mpl = 2.4d0*10**(18.d0)
      MZ = 91.2d0
      MX  = 5.d0*(10.d0**19.d0) 
      Mgut = e1
!      Mgut = 2.50d0*10**(17.d0)
!      msusy = 9.0*10**(4.d0)
      Mtpole=173.1d0
      Mbpole=4.78d0
      MW = 80.38d0
      cb = dCos(beta)
      sb = dSin(beta)
      tw = 0.546535d0                   !--------check it

      cos2beta = cb*cb - sb*sb
      sinsqtw = (1.d0 - (MW**2.d0/MZ**2.d0))      
      stw = dsqrt(sinsqtw)
      ctw = dsqrt(1.d0 - sinsqtw)
      sgnmu = 1.d0

C       Some definitions
C       ================

      tpl   =  dLog(MX**2.d0/Mpl**2.d0)
      tq0   =  dLog(MX**2.d0/msusy**2.d0)
      tgut   =  dLog(MX**2.d0/Mgut**2.d0)
      tZ    =  dLog(MX**2.d0/MZ**2.d0)
      mtscale = dlog(MX**2.d0/Mtpole**2.d0)
      mbscale = dlog(MX**2.d0/mbpole**2.d0)
      QCDscale = dlog(MX**2.d0/2.d0**2.d0)

C    ================================================================
C  fuGUT, fdGUT are diagonal yukawa at GUT scale,yugut, yegut and 
C  ydgut is yukawa coupling at Gut scale
C   ================================================================

!      MHC=2.50d0*10**(17.d0)
      MHC= e1

!      open(1001,FILE='gluinochargino.out',Access='Append')
       i0 = 0

       do i = 1,3

         yugut(1,i) = (4.d0*pi)*yukgut(i0 + i)
         j = 3 + i
         yugut(2,i) = (4.d0*pi)*yukgut(i0 + j)
         k = 6 + i
         yugut(3,i) = (4.d0*pi)*yukgut(i0 + k)
!      Write(*,*)'yugut1i,2i,3i',yugut(1,i),yugut(2,i),yugut(3,i),i
      enddo 

c$$$C     Bottom Yukawa !!!
c$$$C      ----------------------------------
         i0 = 9
       do i = 1,3

         ydgut(1,i) = (4.d0*pi)*yukgut(i0 + i)
         j = 3 + i
         ydgut(2,i) = (4.d0*pi)*yukgut(i0 + j)
         k = 6 + i
         ydgut(3,i) = (4.d0*pi)*yukgut(i0 + k)
!      Write(*,*)'ydgut1i,2i,3i',ydgut(1,i),ydgut(2,i),ydgut(3,i),i
      enddo

c$$$C     Tau Yukawa !!!
c$$$C      ----------------------------------

         i0 = 18

          do i = 1,3

         yegut(1,i) = (4.d0*pi)*yukgut(i0 + i)
         j = 3 + i
         yegut(2,i) = (4.d0*pi)*yukgut(i0 + j)
         k = 6 + i
         yegut(3,i) = (4.d0*pi)*yukgut(i0 + k)
!      Write(*,*)'yegut1i,2i,3i',yegut(1,i),yegut(2,i),yegut(3,i),i
      enddo 


      call SVD(3, 3, yugut,3, fudiag, UQu,3, UTU,3, 0)
      call SVD(3, 3, ydgut,3, fddiag, UQd,3, UTD,3, 0)

      call dag(UQd,UQdT)
      call matmult(UQu,UQdT,VCKMGUT)

      fuGUT(1,1)=fudiag(1)
      fuGUT(2,2)=fudiag(2)
      fuGUT(3,3)=fudiag(3)
      fdGUT(1,1)=fddiag(1)
      fdGUT(2,2)=fddiag(2)
      fdGUT(3,3)=fddiag(3)
!      write(*,*)"fudig",fuGUT(1,1),fuGUT(2,2),fuGUT(3,3)
!      write(*,*)"fddig",fdGUT(1,1),fdGUT(2,2),fdGUT(3,3)
!      write(*,*)"VCKMGUT",VCKMGUT(1,1),VCKMGUT(1,2),VCKMGUT(1,3),
!     $          VCKMGUT(2,1),VCKMGUT(2,2),VCKMGUT(2,3),VCKMGUT(3,1),
!     $           VCKMGUT(3,2),VCKMGUT(3,3)
      fuGUT(1,2)=0.0
      fuGUT(1,3)=0.0
      fuGUT(2,1)=0.0
      fuGUT(2,3)=0.0
      fuGUT(3,1)=0.0
      fuGUT(3,2)=0.0
      fdGUT(1,2)=0.0
      fdGUT(1,3)=0.0
      fdGUT(2,1)=0.0
      fdGUT(2,3)=0.0
      fdGUT(3,1)=0.0
      fdGUT(3,2)=0.0

         alph1gut = yukgut(121)*4.d0*pi
         alph2gut = yukgut(120)*4.d0*pi
         alph3gut= yukgut(119)*4.d0*pi
         vev1gut = yukgut(125)
         vev2gut = yukgut(126)
!      write(*,*)alph1gut,alph2gut,alph3gut,vev1gut,vev2gut     
!      VCKMGUT(1,1)=0.97431
!      VCKMGUT(1,2)=0.22516
!      VCKMGUT(1,3)=3.55131d0*10**(-3.d0)
!      VCKMGUT(2,1)=-0.22511
!      VCKMGUT(2,2)=0.9734d0
!      VCKMGUT(2,3)=4.12291d0*10**(-2.d0)
!      VCKMGUT(3,1)=5.8263d0*10**(-3.d0)
!      VCKMGUT(3,2)=-4.0969d0*10**(-2.d0)
!      VCKMGUT(3,3)=0.9991d0



!      yugut(1,1) = 2.8539d0*10**(-6.d0)
!      yugut(1,2) = -3.0427d0*10**(-17.d0)
!      yugut(1,3) = 1.1093d0*10**(-14.d0)
!      yugut(2,1) =-3.0427d0*10**(-14.d0)
!      yugut(2,2) = 2.8540d0*10**(-3.d0)
!      yugut(2,3) = -2.0092d0*10**(-11.d0)
!      yugut(3,1) = 1.41147d0*10**(-9.d0)
!      yugut(3,2) = -2.5564d0*10**(-9.d0)
!      yugut(3,3) = 0.5699d0

!      ydgut(1,1) = 3.7656d0*10**(-5.d0)
!      ydgut(1,2) = 1.7420d0*10**(-4.d0)
!      ydgut(1,3) = 8.9825d0*10**(-5.d0)
!      ydgut(2,1) = -8.7102d0*10**(-6.d0)
!      ydgut(2,2) = 7.5314d0*10**(-4.d0)
!      ydgut(2,3) = 1.0428d0*10**(-3.d0)
!      ydgut(3,1) = 3.3125d0*10**(-7.d0)
!      ydgut(3,2) =-3.1888d0*10**(-5.d0)
!      ydgut(3,3) = 2.5271d0*10**(-2.d0)

!      yegut(1,1) = 9.3223d0*10**(-6.d0)
!      yegut(1,2) =0.0
!      yegut(1,3) =0.0
!      yegut(2,1) =0.0
!      yegut(2,2) =1.9473d0*10**(-3.d0)
!      yegut(2,3) =0.0
!      yegut(3,1) =0.0
!      yegut(3,2) =0.0
!      yegut(3,3) =3.2356d0*10**(-2.d0)  

!      vev1gut = 51.64
!      vev2gut = 169.274

!       alph3gut= 0.0459d0
!       alph2gut= 0.0461d0
!       alph1gut= 0.04014d0

!      yugut(1,1) = 2.88d0*10**(-6.d0)
!      yugut(1,2) = 6.67d0*10**(-7.d0)
!      yugut(1,3) = 1.4344d0*10**(-8.d0)
!      yugut(2,1) =-6.67d0*10**(-4.d0)
!      yugut(2,2) = 2.885d0*10**(-3.d0)
!      yugut(2,3) = 1.21d0*10**(-4.d0)
!      yugut(3,1) = 5.06d0*10**(-3.d0)
!      yugut(3,2) = -2.43d0*10**(-2.d0)
!      yugut(3,3) = 0.5909d0

!      ydgut(1,1) = 3.927d0*10**(-5.d0)
!      ydgut(1,2) = -2.015d0*10**(-9.d0)
!      ydgut(1,3) = 4.78d0*10**(-8.d0)
!      ydgut(2,1) = -4.030d0*10**(-8.d0)
!      ydgut(2,2) = 7.856d0*10**(-4.d0)
!      ydgut(2,3) = -4.605d0*10**(-6.d0)
!      ydgut(3,1) = 2.733d0*10**(-5.d0)
!      ydgut(3,2) =-1.3166d0*10**(-4.d0)
!      ydgut(3,3) = 2.5675d0*10**(-2.d0)

!      yegut(1,1) = 9.413d0*10**(-6.d0)
!      yegut(1,2) =-6.94d0*10**(-11.d0)
!      yegut(1,3) =1.65d0*10**(-9.d0)
!      yegut(2,1) =-1.45d0*10**(-8.d0)
!      yegut(2,2) =1.96d0*10**(-3.d0)
!      yegut(2,3) =-1.66d0*10**(-6.d0)
!      yegut(3,1) =5.75d0*10**(-6.d0)
!      yegut(3,2) =-2.77d0*10**(-5.d0)
!      yegut(3,3) =3.3344d0*10**(-2.d0)

!      vev1gut = 51.56
!      vev2gut = 169.098

!       alph3gut= 0.0405d0
!       alph2gut= 0.0458d0
!       alph1gut= 0.0446d0

      call dag(VCKMGUT,VCKMGUTT)

!---------write value of all rotation matrix and VCKM at susy scale ,MW,sb,cb and mass of fermion at msusy
      do i=1,3
       do j=1,3
       if(i .eq.j)then 
       delta(i,j) = 1.d0
       else
       delta(i,j) = 0.d0

       endif
!       Write(*,*)"VCKMGUT",VCKMGUT(i,j)
       end do
      end do
!------------------------initializing Arrays to zero

      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
      C5L(i,j,k,l)=0.0
      C5R(i,j,k,l)=0.0
         end do
        end do
       end do
      end do
      call matmult(VCKMGUT,fdGUT,Vfd)
      call matmult(fuGUT,VCKMGUT,fuV)
      call matmult(VCKMGUTT,fuGUT,VTfu)
      call matmult(VCKMGUTT,fdGUT,VTfd)
      call mat3prod(VCKMGUTT,fuGUT,VCKMGUT,VTYUV)
!---------------------------------------------------
!-------check normalization factor of C's RG's for yukawa and alpha's
!      do i=1,3
!       do j=1,3
!        do k=1,3
!         do l=1,3
!!          do n=1,3
!!      C5L(i,j,k,l)=(1.d0/MHC)*(fuGUT(i)*delta(i,j)*VCKMGUTT(k,l)*
!!     $             fdGUT(l))
!!      C5R(i,j,k,l)=(1.d0/MHC)*(fuGUT(i)*VCKMGUT(i,j)*VCKMGUTT(k,l)*
!!    $              fdGUT(j)) 

!!      C5L(i,j,k,l)=(2.d0/MHC)*(fuGUT(i)*delta(i,j)*VCKMGUTT(k,l)*
!!     $             fdGUT(l))
!!      C5R(i,j,k,l)=(2.d0/MHC)*(VCKMGUTT(i,j)*fdGUT(j)*fuGUT(k)*
!!     $              VCKMGUT(k,l))             !<------------------------should check?
!!             do m=1,3
! !     C5L(i,j,k,l)=(1.d0/MHC)*fdGUT(i,j)*VTYUV(k,l)
!!             enddo 
!!      C5R(i,j,k,l)=(1.d0/MHC)*Vfd(i,j)*VTfu(k,l)
!!      C5L(i,j,k,l)=-(1.d0/MHC)*fuGUT(i,j)*VTfd(k,l)
!!      C5R(i,j,k,l)=(1.d0/MHC)*(fuV(i,j)*VTfd(k,l)) 

!       C5L(i,j,k,l)=(1.d0/MHC)*Vfd(i,j)*fuGUT(k,l)
!       C5R(i,j,k,l)=(1.d0/MHC)*Vfd(i,j)*fuV(l,k)

!!      C5L(i,j,k,l)=(1.d0/MHC)*fuGUT(i,j)*Vfd(k,l)
!!      C5R(i,j,k,l)=(1.d0/MHC)*(fuV(i,j)*Vfd(k,l)) 
!!       C5L(i,j,k,l)=C5L(i,j,k,l)+(1.d0/MHC)*Vfd(i,j)
!!     $               *fuGUT(n,l)*VCKMGUT(n,k)
!!          enddo
!!       Write(*,*)C5L(i,j,k,l),C5R(i,j,k,l)
!         enddo
!        end do
!       end do
!      end do
!      call mat3prod(VCKMGUTT,fuGUT,VCKMGUT,VTYUV)
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
       C5L(i,j,k,l)=(1.d0/MHC)*Vfd(i,j)*fuGUT(k,l)
       C5R(i,j,k,l)=(1.d0/MHC)*Vfd(i,j)*fuV(l,k)

!       C5L(i,j,k,l)=(1.d0/MHC)*fdGUT(i,j)*VTYUV(k,l)
!       C5R(i,j,k,l)=(1.d0/MHC)*Vfd(i,j)*VTfu(k,l)

         enddo
        end do
       end do
      end do

!       Write(*,*)C5L(1,1,1,1),C5R(1,1,1,1)
!------------------------------------------------------------------
!-----------------------------------------------------
      do i = 1,31
        yy_sm(i) = 0.d0
         pyye(i) = 0.d0
         pyyb(i) = 0.d0
         pyy2(i) = 0.d0
      enddo

      do i = 1,518
         pyyg(i) = 0.d0
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     ===================================
C     Running of the Yukawas of the SM !!
C     ===================================
C-----------------------------
C	top Yukawa
C-----------------------------
 
      i0 = 0

      inp1: do i = 1, 3
      pyyg(i0 + i)   = yugut(1,i)/(4*pi)
      pyyg(i0+3 + i) = yugut(2,i)/(4*pi)
      pyyg(i0+6 + i) = yugut(3,i)/(4*pi)
      enddo inp1
C-------------------------------
C	bottom Yukawa
C-------------------------------

      i0 = 9

      inp2: do i = 1, 3
      pyyg(i0 + i)   = ydgut(1,i)/(4*pi)
      pyyg(i0+3 + i) = ydgut(2,i)/(4*pi)
      pyyg(i0+6 + i) = ydgut(3,i)/(4*pi)
      enddo inp2
C-----------------------------
C     tau Yukawa
C-----------------------------

      i0 = 18

      inp3: do i = 1, 3
      pyyg(i0 + i)   = yegut(1,i)/(4*pi)
      pyyg(i0+3 + i) = yegut(2,i)/(4*pi)
      pyyg(i0+6 + i) = yegut(3,i)/(4*pi)
      enddo inp3
C-----------------------------
C     gauge coupling
C-----------------------------  
      pyyg(28) = alph3gut/(4*pi)
      pyyg(29) = alph2gut/(4*pi)
      pyyg(30) = alph1gut/(4*pi)

C-----------------------------
C     C5L,C5R
C-----------------------------
       a=1
       do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
       pyyg(30+a)=C5L(i,j,k,l)
       pyyg(30+81+a)=C5R(i,j,k,l)
       a=a+1
         end do
        end do
       end do
      end do 
      
      matchdo10: do i = 1,324
      pyyg(i+192) = 0.d0
      enddo matchdo10

      pyyg(517)=vev1gut
      pyyg(518)=vev2gut 
!---here we have to check we need modified mgut or above defined mgut       !<---check
!      tpl   =  dLog(MX**2.d0/e1**2.d0
!-----------------------Integrating Mgut->Msusy
!----------------------------------
!      write(*,*)"tq0,tgut",tq0,tgut,msusy,mgut
      x2=tq0
      x1=tgut
      n0   = 518
      h1   =  1.d-4
      hmin =  1.d-10
      eps  =  1.d-6
    
      call RK4ROUTINE(pyyg,n0,x1,x2,eps,h1,hmin,nok,nbad,protonrge,
     $     QMSRK4,check)

      if(check.eq.100)then
         flags = 'variable underflow '
         return
      endif

      if(maxval(pyyg(1:27)).gt.dsqrt(1.d0/(4.d0*pi)))then
         flags = "NPERTYUK"
         return
      endif
!      Write(*,*)"What the hell is this?"
!      do i =1,192
!        Write(*,*)pyyg(i)
!       enddo
!----------------------------------------------------------
!	matching  mssm at msusy scale
!---------------------------------------------------------- 
!-----------------------------------------------------------------------------------
!  some function define here for Matching at Msusy
! Mixing factors at each vertex are written in terms of mass-diagonalizing matrices
!-----------------------------------------------------------------------------------
!-------define array zero-----------------------
!----------
!  Gluino
!---------
      do I0 = 1,6
      do j = 1,3
      VGLd(I0,j) = 0.d0
      VGRd(I0,j) = 0.d0
      VGLu(I0,j) = 0.d0
      VGRu(I0,j) = 0.d0
      enddo
      enddo

!----------
!  Chargino
!---------
      do I0 = 1,6
      do A = 1,2
      do j = 1,3
      VCLd(I0,A,j) = 0.d0
      VCRd(I0,A,j) = 0.d0
      VCLu(I0,A,j) = 0.d0
      VCRu(I0,A,j) = 0.d0
      VCLnu(I0,A,j) = 0.d0
      enddo
      enddo
      enddo
      do i = 1,3
      do A = 1,2
      do j = 1,3
      VCLl(i,A,j) = 0.d0
      VCRl(i,A,j) = 0.d0
      enddo
      enddo
      enddo
!-------------
!  Neutralino
!-------------
      do I0 = 1,6
      do Abar = 1,4
      do j = 1,3
      VNLd(I0,Abar,j) = 0.d0
      VNRd(I0,Abar,j) = 0.d0
      VNLu(I0,Abar,j) = 0.d0
      VNRu(I0,Abar,j) = 0.d0
      VNLl(I0,Abar,j) = 0.d0
      VNRl(I0,Abar,j) = 0.d0
      enddo
      enddo
      enddo
      do i = 1,3
      do Abar = 1,4
      do j = 1,3
      VNLnu(i,Abar,j) = 0.d0
      enddo
      enddo
      enddo
!------------end of array---------
!----define value of vertex------ 
!-------------------------------define rotation matrix array =0

      do I=1,6
      do J=1,6
      UD(I,J)=0.0
      UDT(I,J)=0.0
      UU(I,J)=0.0
      UUT(I,J)=0.0
      UL(I,J)=0.0
      ULT(I,J)=0.0
      enddo
      enddo  

       do i=1,3
       do j=1,3
       UNeu(i,j) = 0.0
       UNeuT(i,j) = 0.0
       VCKMS(i,j) = 0.0
       VCKMST(i,j) = 0.0
       enddo
       enddo
          
       do i=1,4
       do j=1,4
       UN(i,j) = 0.0
       UNT(i,j) = 0.0
       enddo
       enddo
 
          
       do i=1,2
       do j=1,2
       Uplus(i,j) = 0.0
       UplusT(i,j) = 0.0
       Uminus(i,j) = 0.0
       UminusT(i,j) = 0.0
       enddo
       enddo

C-----------------------------
C     C5L,C5R
C-----------------------------
       a=1
       do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3
      C5L(i,j,k,l)= pyyg(30+a)
      C5R(i,j,k,l)= pyyg(30+81+a)
!      write(*,*)C5L(i,j,k,l),C5R(i,j,k,l)
       a=a+1
         end do
        end do
       end do
      end do 
 
      g3 = Sqrt(16.d0*pi**2.d0*pyyg(28))
      g2 = Sqrt(16.d0*pi**2.d0*pyyg(29))
      g1 = Sqrt(16.d0*pi**2.d0*pyyg(30))

!-----------------value of rotation matrixs

!      UU(1,1)=0.999d0
!      UU(1,2)=2.036971d0*10**(-6.d0)
!      UU(1,3)=-2.4960d0*10**(-8.d0)
!      UU(1,4)=0.0
!      UU(1,5)=0.0
!      UU(1,6)=0.0
!      UU(2,1)=-2.03697d0*10**(-6.d0)
!      UU(2,2)=0.9999
!      UU(2,3)=4.52098d0*10**(-8.d0)
!      UU(2,4)=0.0
!      UU(2,5)=0.0
!      UU(2,6)=0.0
!      UU(3,1)=2.49603d0*10**(-8.d0)
!      UU(3,2)=-4.52097d0*10**(-8.d0)
!      UU(3,3)=0.9999d0
!      UU(3,4)=0.0
!      UU(3,5)=0.0
!      UU(3,6)=0.0
!      UU(4,1)=0.0
!      UU(4,2)=0.0
!      UU(4,3)=0.0
!      UU(4,4)=1.00d0
!      UU(4,5)=8.82958d0*10**(-15.d0)
!      UU(4,6)=-1.40262d0*10**(-14.d0)
!      UU(5,1)=0.0
!      UU(5,2)=0.0
!      UU(5,3)=0.0
!      UU(5,4)=-8.8295d0*10**(-15.d0)
!      UU(5,5)=1.0d0
!      UU(5,6)=2.54053d0*10**(-11.d0)
!      UU(6,1)=0.0
!      UU(6,2)=0.0
!      UU(6,3)=0.0
!      UU(6,4)=1.40262d0*10**(-14.d0)
!      UU(6,5)=-2.54053d0*10**(-11.d0)
!      UU(6,6)=1.000d0

!      UD(1,1)=0.999d0
!      UD(1,2)=2.036971d0*10**(-6.d0)
!      UD(1,3)=-2.4960d0*10**(-8.d0)
!      UD(1,4)=0.0
!      UD(1,5)=0.0
!      UD(1,6)=0.0
!      UD(2,1)=-2.03697d0*10**(-6.d0)
!      UD(2,2)=0.9999d0
!      UD(2,3)=4.52098d0*10**(-8.d0)
!      UD(2,4)=0.0
!      UD(2,5)=0.0
!      UD(2,6)=0.0
!      UD(3,1)=2.49603d0*10**(-8.d0)
!      UD(3,2)=-4.52097d0*10**(-8.d0)
!      UD(3,3)=0.9999d0
!      UD(3,4)=0.0
!      UD(3,5)=0.0
!      UD(3,6)=0.0
!      UD(4,1)=0.0
!      UD(4,2)=0.0
!      UD(4,3)=0.0
!      UD(4,4)=0.9743d0
!      UD(4,5)=-0.2251d0
!      UD(4,6)=5.8263d0*10**(-3.d0)
!      UD(5,1)=0.0
!      UD(5,2)=0.0
!      UD(5,3)=0.0
!      UD(5,4)= 0.2251d0
!      UD(5,5)=0.9734d0
!      UD(5,6)=-4.0969d0*10**(-2.d0)
!      UD(6,1)=0.0
!      UD(6,2)=0.0
!      UD(6,3)=0.0
!      UD(6,4)=3.551d0*10**(-3.d0)
!      UD(6,5)=4.1197d0*10**(-2.d0)
!      UD(6,6)=0.9991d0

!      UL(1,1)=1.0d0
!      UL(1,2)=-3.5645d0*10**(-12.d0)
!      UL(1,3)=-2.0720d0*10**(-12.d0)
!      UL(1,4)=0.0
!      UL(1,5)=0.0
!      UL(1,6)=0.0
!      UL(2,1)=3.5647*10**(-12.d0)
!      UL(2,2)=0.9999d0
!      UL(2,3)=-5.0516d0*10**(-9.d0)
!      UL(2,4)=0.0
!      UL(2,5)=0.0
!      UL(2,6)=0.0
!      UL(3,1)=2.0983d0*10**(-12.d0)
!      UL(3,2)=5.1153d0*10**(-9.d0)
!      UL(3,3)=0.9999d0
!      UL(3,4)=0.0
!      UL(3,5)=0.0
!      UL(3,6)=0.0
!      UL(4,1)=0.0
!      UL(4,2)=0.0
!      UL(4,3)=0.0
!      UL(4,4)=0.9999d0
!      UL(4,5)=-3.8196d0*10**(-4.d0)
!      UL(4,6)=-1.3053d0*10**(-5.d0)
!      UL(5,1)=0.0
!      UL(5,2)=0.0
!      UL(5,3)=0.0
!      UL(5,4)=3.8196d0*10**(-4.d0)
!      UL(5,5)=0.9999
!      UL(5,6)=-1.5235d0*10**(-4.d0)
!      UL(6,1)=0.0
!      UL(6,2)=0.0
!      UL(6,3)=0.0
!      UL(6,4)=1.3111d0*10**(-5.d0)
!      UL(6,5)=1.5234d0*10**(-4.d0)
!      UL(6,6)=0.9999d0

!      UNeu(1,1)=1.0d0
!      UNeu(1,2)=-3.4766d0*10**(-12.d0)
!      UNeu(1,3)=-2.023d0*10**(-12.d0)
!      UNeu(2,1)=3.476d0*10**(-12.d0)
!      UNeu(2,2)=1.0d0
!      UNeu(2,3)=-4.93d0*10**(-9.d0)
!      UNeu(3,1)=2.023d0*10**(-12.d0)
!      UNeu(3,2)=4.932d0*10**(-9.d0)
!      UNeu(3,3)=1.0

!      UN(1,1)=0.9999d0
!      UN(1,2)=-2.7979d0*10**(-4.d0)
!      UN(1,3)=-7.6602d0*10**(-4.d0)
!      UN(1,4)=1.18673d0*10**(-3.d0)
!      UN(2,1)=2.8352d0*10**(-4.d0)
!      UN(2,2)=0.9999
!      UN(2,3)=1.4409d0*10**(-3.d0)
!      UN(2,4)=-2.2092d0*10**(-3.d0)
!      UN(3,1)=1.3800d0*10**(-3.d0)
!      UN(3,2)=-2.5815d0*10**(-3.d0)
!      UN(3,3)=0.7071d0
!      UN(3,4)=-0.7071d0
!      UN(4,1)=-2.9733d0*10**(-4.d0)
!      UN(4,2)=5.43346d0*10**(-4.d0)
!      UN(4,3)=0.7071d0
!      UN(4,4)=0.7071d0

!      Uplus(1,1)=0.9999d0
!      Uplus(1,2)=7.6829d0*10**(-4.d0)
!      Uplus(2,1)=-7.6829d0*10**(-4.d0)
!      Uplus(2,2)=0.9999d0

!      Uminus(1,1)=-0.9999d0
!      Uminus(1,2)=-3.6502d0*10**(-3.d0)
!      Uminus(2,1)=3.6502d0*10**(-3.d0)
!      Uminus(2,2)=-0.9999

!      msup(1)=97520.20d0
!      msup(2)=97519.33d0
!      msup(3)=78698.59d0
!      msup(4)=98032.43d0
!      msup(5)=98030.89d0
!      msup(6)= 54494.30d0

!      msdown(1)=97520.20d0
!      msdown(2)=97519.33d0
!      msdown(3)=78698.51d0
!      msdown(4)=98112.75d0
!      msdown(5)=98112.55d0
!      msdown(6)=97936.05d0

!      mslepton(1)=99339.28d0
!      mslepton(2)=99338.957d0
!      mslepton(3)=99248.226d0
!      mslepton(4)=99830.42d0
!      mslepton(5)=99829.758d0
!      mslepton(6)=99647.58d0

!      msnu(1)=99339.26d0
!      msnu(2)=99338.93d0
!      msnu(3)=99248.203d0

!      MN(1)=494.477d0
!      MN(2)=334.909d0
!      MN(3)=-30624.68d0
!      MN(4)=30624.794d0

!      MCh(1)=334.90d0
!      MCh(2)=30624.8d0



!      VCKMS(1,1)=0.9743d0
!      VCKMS(1,2)=0.2251d0
!      VCKMS(1,3)=3.5513d0*10**(-3.d0)
!      VCKMS(2,1)=-0.2251d0
!      VCKMS(2,2)=0.9734d0
!      VCKMS(2,3)= 4.1229d0*10**(-2.d0)
!      VCKMS(3,1)= 5.8263d0*10**(-3.d0)
!      VCKMS(3,2)= -4.0969d0*10**(-2.d0)
!      VCKMS(3,3)= 0.9991d0

!      ml(1)=4.7387d0*10**(-4.d0)
!      ml(2)=9.8985d0*10**(-2.d0)
!      ml(3)=1.6432d0

!      mu(1)=9.0852d0*10**(-4.d0)
!      mu(2)=0.9085d0
!      mu(3)=135.44d0
 
!      md(1)=3.6384d0*10**(-3.d0)
!      md(2)= 7.2828d0*10**(-2.d0)
!      md(3)=2.15957d0

      call SVD(3, 3, yuRG,3, fu, UQuS,3, UTUS,3, 0)
      call SVD(3, 3, ydRg,3, fd, UQdS,3, UTDS,3, 0)
      call SVD(3, 3, YeRG,3, fe, ULS,3, UTES,3, 0)

      do i=1,3
      mu(i)=fu(i)*vev2/Sqrt(2.d0)
      md(i)=fd(i)*vev1/Sqrt(2.d0)
      ml(i)=fe(i)*vev1/Sqrt(2.d0)
      enddo
      Print*,"mu_susy=",mu(1),mu(2),mu (3),"md_susy=",
     $        md(1),md(2),md(3),"me_susy=",ml(1),
     $        ml(2),ml(3)

      call dag(UQdS,UQdTS)
      call matmult(UQuS,UQdTS,VCKMS)

      loopmuli:  do i = 1, 3
      loopmulj:  do j = 1, 3
!      write(*,*)fu(i),fd(i),fe(i),VCKMS(i,j),AURG(i,j),
!     $    ADRG(i,j),AERG(i,j)
      AURG(i,j) = AURG(i,j) * yuRG(i,j)
      
      ADRG(i,j) = ADRG(i,j) * ydRG(i,j)

      AERG(i,j) = AERG(i,j) * yeRG(i,j)
!      write(*,*)fu(i),fd(i),fe(i),VCKMS(i,j),AURG(i,j),
!     $    ADRG(i,j),AERG(i,j)
      enddo loopmulj
      enddo loopmuli

      
      mSU(1,1) = mSQRG(1,1) + mu(1)**2.d0 + MZ*MZ*cos2beta*(0.5d0 - 
     .	         (2.d0/3.d0)*sinsqtw)
      mSU(1,2) = mSQRG(1,2) 
      mSU(1,3) = mSQRG(1,3)
      mSU(1,4) =  AURG(1,1)*vev2/dsqrt(2.d0)-
     $     mu(1)*sgnmu*mur*(1.d0/tanbeta) 
      mSU(1,5) =  AURG(1,2)*vev2/dsqrt(2.d0) 
      mSU(1,6) =  AURG(1,3)*vev2/dsqrt(2.d0)
      mSU(2,1) = mSQRG(2,1)
      mSU(2,2) = mSQRG(2,2) + mu(2)**2.d0 + MZ*MZ*cos2beta*(0.5d0 -
     .     (2.d0/3.d0)*sinsqtw)
      mSU(2,3) = mSQRG(2,3)
      mSU(2,4) =  AURG(2,1)*vev2/dsqrt(2.d0) 
      mSU(2,5) =  AURG(2,2)*vev2/dsqrt(2.d0)-
     $     mu(2)*sgnmu*mur*(1.d0/tanbeta)
      mSU(2,6) =  AURG(2,3)*vev2/dsqrt(2.d0) 
      mSU(3,1) = mSQRG(3,1)
      mSU(3,2) = mSQRG(3,2)
      mSU(3,3) = mSQRG(3,3) + mu(3)**2.d0 + MZ*MZ*cos2beta*(0.5d0 
     $     - (2.d0/3.d0)*sinsqtw)
      mSU(3,4) =  AURG(3,1)*vev2/dsqrt(2.d0) 
      mSU(3,5) =  AURG(3,2)*vev2/dsqrt(2.d0) 
      mSU(3,6) =  AURG(3,3)*vev2/dsqrt(2.d0) - 
     $     mu(3)*sgnmu*mur*(1.d0/tanbeta)
      mSU(4,1) =  AURG(1,1)*vev2/dsqrt(2.d0) - 
     $     mu(1)*sgnmu*mur*(1.d0/tanbeta)
      mSU(4,2) =  AURG(2,1)*vev2/dsqrt(2.d0)
      mSU(4,3) =  AURG(3,1)*vev2/dsqrt(2.d0)
      mSU(4,4) = mSURG(1,1)+mu(1)**2.d0+(2.d0/3.d0)*MZ*MZ*cos2beta*
     $           sinsqtw
      mSU(4,5) = mSURG(1,2)
      mSU(4,6) = mSURG(1,3)
      mSU(5,1) =  AURG(1,2)*vev2/dsqrt(2.d0)
      mSU(5,2) =  AURG(2,2)*vev2/dsqrt(2.d0) -
     $     mu(2)*sgnmu*mur*(1.d0/tanbeta)
      mSU(5,3) =  AURG(3,2)*vev2/dsqrt(2.d0)
      mSU(5,4) = mSURG(2,1)
      mSU(5,5) = mSURG(2,2)+mu(2)**2.d0+(2.d0/3.d0)*MZ*MZ*cos2beta*
     $           sinsqtw
      mSU(5,6) = mSURG(2,3)
      mSU(6,1) =  AURG(1,3)*vev2/dsqrt(2.d0)
      mSU(6,2) =  AURG(2,3)*vev2/dsqrt(2.d0)
      mSU(6,3) =  AURG(3,3)*vev2/dsqrt(2.d0) - 
     $            mu(2)*sgnmu*mur*(1.d0/tanbeta)
      mSU(6,4) = mSURG(3,1)
      mSU(6,5) = mSURG(3,2)
      mSU(6,6) = mSURG(3,3) + mu(3)**2.d0 + (2.d0/3.d0)*MZ*MZ*cos2beta*
     $           sinsqtw
      call dag6(mSU,mSUT)
      call CEigensystem(6, mSUT,6, msupsq, UU,6, 0)
      do i=1,6
      msup(i)=dSqrt(abs(msupsq(i)))
      enddo
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

      mSD(1,1) = mSQRG(1,1) + md(1)**2.d0 + MZ*MZ*cos2beta*(- 0.5d0 + 
     .	          (1.d0/3.d0)*sinsqtw) 
      mSD(1,2) = mSQRG(1,2) 
      mSD(1,3) = mSQRG(1,3)
      mSD(1,4) = ADRG(1,1)*vev1/dsqrt(2.d0) - md(1)*sgnmu*mur*tanbeta 
      mSD(1,5) = ADRG(1,2)*vev1/dsqrt(2.d0) 
      mSD(1,6) = ADRG(1,3)*vev1/dsqrt(2.d0)
      mSD(2,1) = mSQRG(2,1)
      mSD(2,2) = mSQRG(2,2) + md(2)**2.d0 + MZ*MZ*cos2beta*(- 0.5d0 + 
     . 	          (1.d0/3.d0)*sinsqtw)
      mSD(2,3) = mSQRG(2,3)
      mSD(2,4) = ADRG(2,1)*vev1/dsqrt(2.d0) 
      mSD(2,5) = ADRG(2,2)*vev1/dsqrt(2.d0) - md(2)*sgnmu*mur*tanbeta
      mSD(2,6) = ADRG(2,3)*vev1/dsqrt(2.d0) 
      mSD(3,1) = mSQRG(3,1)
      mSD(3,2) = mSQRG(3,2)
      mSD(3,3) = mSQRG(3,3) + md(3)**2.d0 + MZ*MZ*cos2beta*(- 0.5d0 + 
     .            (1.d0/3.d0)*sinsqtw)
      mSD(3,4) = ADRG(3,1)*vev1/dsqrt(2.d0) 
      mSD(3,5) = ADRG(3,2)*vev1/dsqrt(2.d0) 
      mSD(3,6) = ADRG(3,3)*vev1/dsqrt(2.d0) - md(3)*sgnmu*mur*tanbeta
      mSD(4,1) = ADRG(1,1)*vev1/dsqrt(2.d0) - md(1)*sgnmu*mur*tanbeta
      mSD(4,2) = ADRG(2,1)*vev1/dsqrt(2.d0)
      mSD(4,3) = ADRG(3,1)*vev1/dsqrt(2.d0)
      mSD(4,4) = mSDRG(1,1) + md(1)**2.d0 - (1.d0/3.d0)*MZ*MZ*
     .                  cos2beta*sinsqtw
      mSD(4,5) = mSDRG(1,2)
      mSD(4,6) = mSDRG(1,3)
      mSD(5,1) = ADRG(1,2)*vev1/dsqrt(2.d0)
      mSD(5,2) = ADRG(2,2)*vev1/dsqrt(2.d0) - md(2)*sgnmu*mur*tanbeta
      mSD(5,3) = ADRG(3,2)*vev1/dsqrt(2.d0)
      mSD(5,4) = mSDRG(2,1)
      mSD(5,5) = mSDRG(2,2) + md(2)**2.d0 - (1.d0/3.d0)*MZ*MZ*
     .            cos2beta*sinsqtw
      mSD(5,6) = mSDRG(2,3)
      mSD(6,1) = ADRG(1,3)*vev1/dsqrt(2.d0)
      mSD(6,2) = ADRG(2,3)*vev1/dsqrt(2.d0)
      mSD(6,3) = ADRG(3,3)*vev1/dsqrt(2.d0) - md(3)*sgnmu*mur*tanbeta
      mSD(6,4) = mSDRG(3,1)
      mSD(6,5) = mSDRG(3,2)
      mSD(6,6) = mSDRG(3,3) + md(3)**2.d0 - (1.d0/3.d0)*MZ*MZ*
     .            cos2beta*sinsqtw
      call dag6(mSD,mSDT)
      call CEigensystem(6, mSDT,6, msdownsq, UD,6, 0)
      do i=1,6
      msdown(i)=dSqrt(abs(msdownsq(i)))
      enddo
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------  

      mSL(1,1) = mSLRG(1,1)+ml(1)**2.d0+MZ*MZ*cos2beta*
     $           (- 0.5d0 + sinsqtw) 
      mSL(1,2) = mSLRG(1,2) 
      mSL(1,3) = mSLRG(1,3)
      mSL(1,4) = AERG(1,1)*vev1/dsqrt(2.d0) - ml(1)*sgnmu*mur*tanbeta 
      mSL(1,5) = AERG(1,2)*vev1/dsqrt(2.d0) 
      mSL(1,6) = AERG(1,3)*vev1/dsqrt(2.d0)
      mSL(2,1) = mSLRG(2,1)
      mSL(2,2) = mSLRG(2,2)+ml(2)**2.d0+MZ*MZ*cos2beta*(-0.5d0+sinsqtw)
      mSL(2,3) = mSLRG(2,3)
      mSL(2,4) = AERG(2,1)*vev1/dsqrt(2.d0) 
      mSL(2,5) = AERG(2,2)*vev1/dsqrt(2.d0) - ml(2)*sgnmu*mur*tanbeta
      mSL(2,6) = AERG(2,3)*vev1/dsqrt(2.d0) 
      mSL(3,1) = mSLRG(3,1)
      mSL(3,2) = mSLRG(3,2)
      mSL(3,3) = mSLRG(3,3)+ml(3)**2.d0+MZ*MZ*cos2beta*(-0.5d0+sinsqtw)
      mSL(3,4) = AERG(3,1)*vev1/dsqrt(2.d0) 
      mSL(3,5) = AERG(3,2)*vev1/dsqrt(2.d0) 
      mSL(3,6) = AERG(3,3)*vev1/dsqrt(2.d0) - ml(3)*sgnmu*mur*tanbeta
      mSL(4,1) = AERG(1,1)*vev1/dsqrt(2.d0) - ml(1)*sgnmu*mur*tanbeta
      mSL(4,2) = AERG(2,1)*vev1/dsqrt(2.d0)
      mSL(4,3) = AERG(3,1)*vev1/dsqrt(2.d0)
      mSL(4,4) = mSERG(1,1) + ml(1)**2.d0 - MZ*MZ*cos2beta*sinsqtw
      mSL(4,5) = mSERG(1,2)
      mSL(4,6) = mSERG(1,3)
      mSL(5,1) = AERG(1,2)*vev1/dsqrt(2.d0)
      mSL(5,2) = AERG(2,2)*vev1/dsqrt(2.d0) - ml(2)*sgnmu*mur*tanbeta
      mSL(5,3) = AERG(3,2)*vev1/dsqrt(2.d0)
      mSL(5,4) = mSERG(2,1)
      mSL(5,5) = mSERG(2,2) + ml(2)**2.d0 - MZ*MZ*cos2beta*sinsqtw
      mSL(5,6) = mSERG(2,3) 
      mSL(6,1) = AERG(1,3)*vev1/dsqrt(2.d0)
      mSL(6,2) = AERG(2,3)*vev1/dsqrt(2.d0)
      mSL(6,3) = AERG(3,3)*vev1/dsqrt(2.d0) - ml(3)*sgnmu*mur*tanbeta
      mSL(6,4) = mSERG(3,1)
      mSL(6,5) = mSERG(3,2) 
      mSL(6,6) = mSERG(3,3) + ml(3)**2.d0 - MZ*MZ*cos2beta*sinsqtw

!      do i=1,6

!      Write(*,*)"msl",mSL(1,i),mSL(2,i),mSL(3,i),mSL(4,i),mSL(5,i),
!     $          mSL(6,i)

!      enddo
      call CEigensystem(6, mSL,6, msleptonsq, Lo,6, 0)
!         call SVD(6, 6, mSL,6, msleptonsq, Lo,6, LoT,6, -1)
      call dag6(Lo,UL)
      do i=1,6
      mslepton(i)=dSqrt(abs(msleptonsq(i)))

!      write(*,*)mslepton(i),msleptonsq(i),i,vev1,cos2beta,sinsqtw,
!     $           mur,sgnmu,tanbeta,MZ
      enddo
!      write(*,*)"nsl",mSLRG(1,1),AERG(3,3),mSLRG(2,2),mSLRG(3,3)
!----------------------------------------------------------------------
!------------------------------------------------------------------------

      mSN(1,1) = mSLRG(1,1) + 0.5d0*MZ*MZ*cos2beta
      mSN(1,2) = mSLRG(1,2)
      mSN(1,3) = mSLRG(1,3)
      mSN(2,1) = mSLRG(2,1)
      mSN(2,2) = mSLRG(2,2) + 0.5d0*MZ*MZ*cos2beta
      mSN(2,3) = mSLRG(2,3)
      mSN(3,1) = mSLRG(3,1)
      mSN(3,2) = mSLRG(3,2)
      mSN(3,3) = mSLRG(3,3) + 0.5d0*MZ*MZ*cos2beta

      Call CEigensystem(3,mSN,3,msnusq,usnu,3,0)
      call dag(usnu,UNeu)
      do i=1,3
      msnu(i)=dsqrt(abs(msnusq(i)))
      enddo
     
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

      MNeut(1,1) = M1tz
      MNeut(1,2) = 0.d0
      MNeut(1,3) = - (MZ*stw*cb)
      MNeut(1,4) = MZ*stw*sb
      MNeut(2,1) = 0.d0
      MNeut(2,2) = M2tz
      MNeut(2,3) = MZ*ctw*cb
      MNeut(2,4) = -(MZ*ctw*sb)
      MNeut(3,1) = -(MZ*stw*cb)
      MNeut(3,2) =  MZ*ctw*cb
      MNeut(3,3) = 0.d0
      MNeut(3,4) =  -sgnmu*mur
      MNeut(4,1) = MZ*stw*sb
      MNeut(4,2) = - (MZ*ctw*sb)
      MNeut(4,3) = - sgnmu*mur 
      MNeut(4,4) = 0.d0 
 

      call CEigensystem(4,MNeut,4,MN,Neuto,4, 0)
      call dag4(Neuto,UN)
!------------------------------------------------------------------------
!------------------------------------------------------------------------
      MChar(1,1) = M2tz  
      MChar(1,2) = sqrt(2.d0)*MW*sb
      MChar(2,1) = sqrt(2.d0)*MW*cb
      MChar(2,2) = sgnmu*mur
      call SVD(2, 2, MChar,2, MCeg, OLC,2, ORC,2, 0)
      call dag2(OLC,Uminus) 
      call dag2(ORC,Uplus) 
      MCh(1)=Mceg(1)
      MCh(2)=MCeg(2)


      MG=M3t
!      write(*,*)Sqrt(mSU(1,1)),Sqrt(mSU(2,2)),Sqrt(mSU(3,3)),
!     $     Sqrt(mSU(4,4)),Sqrt(mSU(5,5)),Sqrt(mSU(6,6))
!      write(*,*)Sqrt(mSD(1,1)),Sqrt(mSD(2,2)),Sqrt(mSD(3,3)),
!     $     Sqrt(mSD(4,4)),Sqrt(mSD(5,5)),Sqrt(mSD(6,6))
!      write(*,*)Sqrt(mSL(1,1)),Sqrt(mSL(2,2)),Sqrt(mSL(3,3)),
!     $     Sqrt(mSL(4,4)),Sqrt(mSL(5,5)),Sqrt(mSL(6,6))
!      write(*,*)Sqrt(mSN(1,1)),Sqrt(mSN(2,2)),Sqrt(mSN(3,3)),
!     $      M1tz,M2tz,M3t,mur,g1,g2,g3
!      write(*,*)mSQRG(1,3)/Sqrt(mSU(1,1)*mSU(3,3)),
!     $          mSQRG(1,3)/Sqrt(mSD(1,1)*mSD(3,3))
!----------
!  Gluino
!---------

      do I0 = 1,6
      do j = 1,3
      loopk1: do k = 1,3
      VGLd(I0,j) = VGLd(I0,j)+UD(I0,k)*VCKMS(k,j)  !-------define rotation matrix array is equal zero and take value for this from suseflav
!      write(*,*)VGLd(I0,j)
              enddo loopk1                         !----VCKMS is VCKM at susy scale ...call from outside

      VGRd(I0,j) = UD(I0,j+3)
      VGLu(I0,j) = UU(I0,j)
      VGRu(I0,j) = UU(I0,j+3)
!      write(*,*)VGLd(I0,j),VGLu(I0,j),VGRu(I0,j)
      enddo
      enddo

      call dag6(UD,UDT)
      call dag6(UU,UUT)
      call dag6(UL,ULT)
      call dag(UNeu,UNeuT)
      call dag4(UN,UNT)
      call dag2(Uplus,UplusT)
      call dag2(Uminus,UminusT)
      call dag(VCKMS,VCKMST)
!----------
!  Chargino
!---------
      do I = 1,6
      do A = 1,2
      do j = 1,3

       loopk2 : do k = 1,3
      VCLd(I,A,j) = VCLd(I,A,j) + (UU(I,k)*Uplus(1,A)+UU(I,k+3)*
     $              (mu(k)/(Sqrt(2.d0)*MW*sb))*Uplus(2,A))*VCKMS(k,j)
!      VCLd(I,A,j) = VCLd(I,A,j) + (UU(I,k+3)*
!     $              (mu(k)/(Sqrt(2.d0)*MW*sb))*Uplus(2,A))*VCKMS(k,j)
      VCRd(I,A,j) = VCRd(I,A,j) - (UU(I,k)*VCKMS(k,j)*
     $              (md(j)/(Sqrt(2.d0)*MW*cb))*Uminus(2,A))
      VCLu(I,A,j) = VCLu(I,A,j)-(UD(I,k+3)*(md(k)/(Sqrt(2.d0)*MW*cb))
     $              *VCKMST(k,j)*UminusT(A,2))
!      write(*,*)"VCLd,VCRd,VCLU",VCLd(I,A,j),VCRd(I,A,j),VCLu(I,A,j),I,
!     $            A,j
                enddo loopk2
      VCLu(I,A,j) = VCLu(I,A,j)+UD(I,j)*UminusT(A,1)
      VCRu(I,A,j) = UD(I,j)*(mu(j)/(Sqrt(2.d0)*MW*sb))*UplusT(A,2)
!      write(*,*)VCLu(I,A,j),VCRu(I,A,j)
      enddo
      enddo
      enddo

      do i = 1,3
      do A = 1,2
      do j = 1,3
      VCLl(i,A,j) = -UNeuT(i,j)*Uplus(1,A)
      VCRl(i,A,j) = (ml(j)/(Sqrt(2.d0)*MW*cb))*UNeuT(i,j)*Uminus(2,A)
!      write(*,*)VCLl(i,A,j),VCRl(i,A,j)
      enddo
      enddo
      enddo

      do I = 1,6
      do A = 1,2
      do j = 1,3
      VCLnu(I,A,j) = -ULT(I,j)*UminusT(A,1)+(ml(j)/(Sqrt(2.d0)*MW*cb))
     $               *ULT(I,j+3)*UminusT(A,2)

      enddo
      enddo
      enddo
!-------------
!  Neutralino
!-------------
      do I = 1,6
      do Abar = 1,4
      do j = 1,3
       loopk3: do k = 1,3
      VNLd(I,Abar,j) = VNLd(I,Abar,j)+dSqrt(2.d0)*(0.5d0*UN(2,Abar)-
     $                 (1.d0/6.d0)*tw*UN(1,Abar))*(UD(I,k)*VCKMS(k,j))
      VNRd(I,Abar,j) = VNRd(I,Abar,j)-(md(j)/(Sqrt(2.d0)*MW*cb))
     $                 *UNT(Abar,3)*UD(I,k)*VCKMS(k,j)
!      Write(*,*)VNLd(I,Abar,j),VNRd(I,Abar,j)
               enddo loopk3
      VNLd(I,Abar,j) = VNLd(I,Abar,j)-(md(j)/(Sqrt(2.d0)*MW*cb))*
     $                 UN(3,Abar)*UD(I,j+3)
      VNRd(I,Abar,j) = VNRd(I,Abar,j)+Sqrt(2.d0)*(-1.d0/3.d0)*tw*
     $                 UNT(Abar,1)*UD(I,j+3)
      VNLu(I,Abar,j) = Sqrt(2.d0)*(-0.5d0*UN(2,Abar)-(1.d0/6.d0)*tw*
     $                 UN(1,Abar))*UU(I,j) - (mu(j)/(Sqrt(2.d0)*MW*sb))
     $                 *UN(4,Abar)*UU(I,j+3)
      VNRu(I,Abar,j) = Sqrt(2.d0)*((2.d0/3.d0)*tw*UNT(Abar,2))*
     $                 UU(I,j+3)-(mu(j)/(Sqrt(2.d0)*MW*sb))*
     $                 UNT(Abar,4)*UU(I,j)
      VNLl(I,Abar,j) = Sqrt(2.d0)*(0.5d0*UN(2,Abar)+(1.d0/2.d0)*tw*
     $                 UN(1,Abar))*ULT(I,j)-(ml(j)/(Sqrt(2.d0)*MW*cb))*
     $                 UN(3,Abar)*ULT(I,j+3)
      VNRl(I,Abar,j) = -Sqrt(2.d0)*tw*UNT(Abar,1)*ULT(I,j+3)-
     $                 (ml(j)/(Sqrt(2.d0)*MW*cb))*UNT(Abar,3)*ULT(I,j)
!      Write(*,*)VNLd(I,Abar,j),VNRd(I,Abar,j),VNLu(I,Abar,j),
!     $          VNRu(I,Abar,j),VNLl(I,Abar,j),VNRl(I,Abar,j)
      enddo
      enddo
      enddo
      do i = 1,3
      do Abar = 1,4
      do j = 1,3
      VNLnu(i,Abar,j) = Sqrt(2.d0)*(-0.5d0*UN(2,Abar)+(1.d0/2.d0)*tw*
     $                 UN(1,Abar))*UNeuT(i,j)
!      Write(*,*)VNLnu(i,Abar,j)
      enddo
      enddo
      enddo

!----------------------------------------------------------------------------
!  The coefficients C's five dimension operator are written in terms of C5L,R
!  First define array equal to zero.first two are with tilde
!--------------------------------------------------------------------------

      do M = 1,6
      do N = 1,6 
      do i = 1,3
      do j = 1,3
      CudulL(M,N,i,j)=0.d0
      CuudlL(M,N,i,j)=0.d0
      CudulR(M,N,i,j)=0.d0
      CuudlR(M,N,i,j)=0.d0
      CuddnuL(M,N,i,j)=0.d0
      CddunuL(M,N,i,j)=0.d0
      CuludL(M,N,i,j)=0.d0
      CdluuL(M,N,i,j)=0.d0
      CuludR(M,N,i,j)=0.d0
      CdluuR(M,N,i,j)=0.d0
      enddo
      enddo
      enddo
      enddo

      do M = 1,6
      do k = 1,3 
      do i = 1,3
      do j = 1,3
      CdnuudL(M,k,i,j)=0.d0
      CunuddL(M,k,i,j)=0.d0
      enddo
      enddo
      enddo
      enddo

!-----------------------------------------------------------------
!------value of Coefficients C------------------------------------

      do M = 1,6
      do N = 1,6 
      do i = 1,3
      do j = 1,3
       loopk4: do k = 1,3
       loopl1: do l = 1,3

      CudulL(M,N,i,j)=CudulL(M,N,i,j) + (C5L(i,j,k,l)-C5L(k,j,i,l))*
     $                UUT(k,M)*UDT(l,N)
      CudulR(M,N,i,j)=CudulR(M,N,i,j) + (C5R(k,l,j,i)-C5R(i,l,j,k))*
     $                UUT(k+3,M)*UDT(l+3,N)
      CuudlR(M,N,i,j)=CuudlR(M,N,i,j) + (C5R(l,i,j,k)-C5R(k,i,j,l))*
     $                UUT(k+3,M)*UUT(l+3,N)
      CddunuL(M,N,i,j)=CddunuL(M,N,i,j) + (C5L(l,j,i,k)-C5L(k,j,i,l))
     $                 *UDT(k,M)*UDT(l,N)
!       Write(*,*)CddunuL(M,N,i,j),CudulR(M,N,i,j),CuudlR(M,N,i,j),
!     $           CddunuL(M,N,i,j)
       loopm1: do m1 = 1,3

      CuudlL(M,N,i,j)=CuudlL(M,N,i,j) + (C5L(k,j,l,m1)-C5L(l,j,k,m1))
     $                *UUT(k,M)*UUT(l,N)*VCKMS(m1,i)
      CuddnuL(M,N,i,j)=CuddnuL(M,N,i,j) + (C5L(m1,j,k,l)-C5L(l,j,k,m1))
     $                 *UUT(k,M)*UDT(l,N)*VCKMS(m1,i)
!      write(*,*)CuudlL(M,N,i,j),CuddnuL(M,N,i,j)
               enddo loopm1
               enddo loopl1
               enddo loopk4
      enddo
      enddo
      enddo
      enddo
!       write(*,*)C5L(1,1,1,1),C5R(1,1,1,1)

      do I = 1,6
      do J = 1,6 
      do k = 1,3
      do l = 1,3
       loopi1: do i1 = 1,3
       loopj1: do j1 = 1,3


      CdluuL(I,J,k,l)=CdluuL(I,J,k,l) + (C5L(k,j1,l,i1)-C5L(l,j1,k,i1))
     $                *UDT(i1,I)*UL(j1,J)
      CuludR(I,J,k,l)=CuludR(I,J,k,l) + (C5R(k,l,j1,i1)-C5R(i1,l,j1,k))
     $                *UUT(i1+3,I)*UL(j1+3,J)
      CdluuR(I,J,k,l)=CdluuR(I,J,k,l) + (C5R(l,i1,j1,k)-C5R(k,i1,j1,l))
     $                *UDT(i1+3,I)*UL(j1+3,J)
!      Write(*,*)CdluuL(I,J,k,l),CuludR(I,J,k,l),CdluuR(I,J,k,l)
       loopm2: do m1 = 1,3

      CuludL(I,J,k,l)=CuludL(I,J,k,l) + (C5L(i1,j1,k,m1)-
     $                C5L(k,j1,i1,m1))*UUT(i1,I)*UL(j1,J)*VCKMS(m1,l)
!      write(*,*)CuludL(I,J,k,l)
               enddo loopm2
               enddo loopj1
               enddo loopi1
      enddo
      enddo
      enddo
      enddo


      do I = 1,6
      do j = 1,3 
      do k = 1,3
      do l = 1,3
           do n1 = 1,3
           do m1 = 1,3
           do i1 = 1,3
      CdnuudL(I,j,k,l)=CdnuudL(I,j,k,l)+(C5L(i1,n1,k,m1)-
     $                C5L(m1,n1,k,i1))*UDT(i1,I)*UNeu(n1,j)*VCKMS(m1,l)
!      write(*,*)CdnuudL(I,j,k,l)
           do p = 1,3
      CunuddL(I,j,k,l)=CunuddL(I,j,k,l)+(C5L(m1,n1,i1,p)-
     $     C5L(p,n1,i1,m1))*UUT(i1,I)*UNeu(n1,j)*VCKMS(p,k)*VCKMS(m1,l)
!      write(*,*)CunuddL(I,j,k,l) 
           enddo
           enddo
           enddo
           enddo
      enddo
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------------------
! after dressing effective four-fermion operator in term of five dimension operator
! at msusy scale
! first define array value is zero
!----------------------------------------------------------------------------------- 

      do M = 1,6
      xG(M) = 0.d0
      uG(M) = 0.d0
!      msdown(M)= 0.d0
!      msup(M)= 0.d0
!      mslepton(M)= 0.d0
       do A = 1,2
       xC(M,A) = 0.d0
       uC(M,A) = 0.d0
       wC(M,A) = 0.d0
       enddo
       do Abar = 1,4
       nuN(M,Abar) = 0.d0
       yN(M,Abar) = 0.d0
       zN(M,Abar) = 0.d0
       enddo
      enddo

      do m1 = 1,3
        do A = 1,2
        zC(m1,A) = 0.d0
        enddo
        do Abar = 1,4
        wN(m1,Abar) = 0.d0
        enddo
      enddo

!        do A = 1,2
!        MCh(A) = 0.d0
!        enddo
!        do Abar = 1,4
!        MN(Abar) = 0.d0
!        enddo
!        do m1 = 1,3
!        msnu(m1) = 0.d0
!        enddo



!----------egein value of susy particles
!---------------------------------------------------
! define value of these parameter
!---------------------------------------------------
      do M = 1,6
      xG(M) = msdown(M)**2.d0/MG**2.d0
      uG(M) = msup(M)**2.d0/MG**2.d0
!      write(*,*)xG(M),uG(M)
       do A = 1,2
       xC(M,A) = msup(M)**2.d0/MCh(A)**2.d0
       uC(M,A) = msdown(M)**2.d0/MCh(A)**2.d0
       wC(M,A) = mslepton(M)**2.d0/MCh(A)**2.d0
!      write(*,*)xC(M,A),uC(M,A),wC(M,A),mslepton(M),MCh(A),M,A
       enddo
       do Abar = 1,4
       nuN(M,Abar) = msup(M)**2.d0/MN(Abar)**2.d0
       yN(M,Abar) = msdown(M)**2.d0/MN(Abar)**2.d0
       zN(M,Abar) = mslepton(M)**2.d0/MN(Abar)**2.d0
!      write(*,*)nuN(M,Abar),yN(M,Abar),zN(M,Abar),M,Abar,mslepton(M),
!     $          MN(Abar)
       enddo
      enddo

      do m1 = 1,3
        do A = 1,2
        zC(m1,A) = msnu(m1)**2.d0/MCh(A)**2.d0
!      write(*,*)zC(m1,A)
        enddo
        do Abar = 1,4
        wN(m1,Abar) = msnu(m1)**2.d0/MN(Abar)**2.d0
!      write(*,*)wN(m1,Abar)
        enddo
      enddo

      alphasmsusy = pyyg(28)*4.d0*pi
      alpha2msusy= pyyg(29)*4.d0*pi
      ALgutmsusy=(alphasmsusy/alph2gut)**(4.d0/9.d0)*
     $           (alpha2msusy/alph2gut)**(-3.d0/2.d0)
!----------------------------------------------------
! define wilson cofficent array is zero
!----------------------------------------------------

      do i = 1,3
      do k = 1,3
      CLRudul6d(i,k)=0.0
      CLRudul6d(i,k)=0.0

      CLLudul(i,k) = 0.d0
      CLLudulG(i,k) = 0.d0
      CLLudulchipm(i,k) = 0.d0
      CLLudulchi0(i,k) = 0.d0

      CRLudul(i,k) = 0.d0
      CRLudulG(i,k) = 0.d0
      CRLudulchipm(i,k) = 0.d0
      CRLudulchi0(i,k) = 0.d0

      CLRudul(i,k) = 0.d0
      CLRudulG(i,k) = 0.d0
      CLRudulchipm(i,k) = 0.d0
      CLRudulchi0(i,k) = 0.d0

      CRRudul(i,k) = 0.d0
      CRRudulG(i,k) = 0.d0
      CRRudulchipm(i,k) = 0.d0
      CRRudulchi0(i,k) = 0.d0
!      write(*,*)CLLudul(i,k),CLRudul(i,k)
      enddo
      enddo

      do i = 1,3
      do j = 1,3
      do k = 1,3
      CRLuddnu6d(i,j,k)=0.0

      CLLuddnu(i,j,k) = 0.d0
      CLLuddnuG(i,j,k) = 0.d0
      CLLuddnuchipm(i,j,k) = 0.d0
      CLLuddnuchi0(i,j,k) = 0.d0

      CLLddunu(i,j,k) = 0.d0
      CLLddunuG(i,j,k) = 0.d0
      CLLddunuchi0(i,j,k) = 0.d0

      CRLuddnu(i,j,k) = 0.d0
      CRLuddnuG(i,j,k) = 0.d0
      CRLuddnuchipm(i,j,k) = 0.d0
      CRLuddnuchi0(i,j,k) = 0.d0

      CRLddunu(i,j,k) = 0.d0
      CRLddunuG(i,j,k) = 0.d0
      CRLddunuchi0(i,j,k) = 0.d0
      enddo
      enddo
      enddo

!---------------------------------------------------
! matching
!---------------------------------------------------

      do i = 1,3
      do k = 1,3
         do M = 1,6
         do N = 1,6
         CLLudulG(i,k) = CLLudulG(i,k) + (4.d0/3.d0)*(g3**2.d0/MG)*
     $          CudulL(M,N,1,k)*VGLu(M,1)*VGLd(N,i)*Hfunc(uG(M),xG(N))
         CRLudulG(i,k) = CRLudulG(i,k) + (4.d0/3.d0)*(g3**2.d0/MG)*
     $          CudulL(M,N,1,k)*VGRu(M,1)*VGRd(N,i)*Hfunc(uG(M),xG(N))
         CLRudulG(i,k) = CLRudulG(i,k) + (4.d0/3.d0)*(g3**2.d0/MG)*
     $          CudulR(M,N,1,k)*VGLu(M,1)*VGLd(N,i)*Hfunc(uG(M),xG(N))
         CRRudulG(i,k) = CRRudulG(i,k) + (4.d0/3.d0)*(g3**2.d0/MG)*
     $         CudulR(M,N,1,k)*VGRu(M,1)*VGRd(N,i)*Hfunc(uG(M),xG(N))
!       write(*,*)CLLudulG(i,k),g3,MG,CudulL(M,N,1,k),VGLu(M,1),
!     $           VGLd(N,i),Hfunc(uG(M),xG(N))
!      write(*,*)Hfunc(uG(M),xG(N)),Hfunc(uG(M),xG(N)),Hfunc(uG(M),xG(N)),
!     $          Hfunc(uG(M),xG(N))
!      Write(*,*)CudulL(M,N,1,k),g3,MG,VGLd(N,i),Hfunc(uG(M),xG(N))
!      Write(*,*)
           do A = 1,2
!      CLLudulchipm(i,k) = CLLudulchipm(i,k) + 
      CRLudulchipm(i,k) = CRLudulchipm(i,k) - (g2**2.d0/MCh(A))*
     $                    CudulL(M,N,1,k)*VCRu(N,A,1)*VCRd(M,A,i)
     $                    *Hfunc(xC(M,A),uC(N,A))
!      CLRudulchipm(i,k) = CLRudulchipm(i,k) 
      CRRudulchipm(i,k) = CRRudulchipm(i,k) - (g2**2.d0/MCh(A))*
     $                    CudulR(M,N,1,k)*VCRu(N,A,1)*VCRd(M,A,i)
     $                    *Hfunc(xC(M,A),uC(N,A))
!      Write(*,*)Hfunc(xC(M,A),uC(N,A)),Hfunc(xC(M,A),uC(N,A))
          enddo
          do Abar = 1,4
      CLLudulchi0(i,k) = CLLudulchi0(i,k) + (g2**2.d0/MN(Abar))*
     $           (CudulL(M,N,1,k)*VNLu(M,Abar,1)*VNLd(N,Abar,i)*
     $            Hfunc(nuN(M,Abar),yN(N,Abar)) + CuludL(M,N,1,i)*
     $    VNLu(M,Abar,1)*VNLl(N,Abar,k)*Hfunc(nuN(M,Abar),zN(N,Abar)))
      CRLudulchi0(i,k) = CRLudulchi0(i,k) + (g2**2.d0/MN(Abar))*
     $           (CudulL(M,N,1,k)*VNRu(M,Abar,1)*VNRd(N,Abar,i)*
     $            Hfunc(nuN(M,Abar),yN(N,Abar)) + CuludR(M,N,1,i)*
     $    VNLu(M,Abar,1)*VNLl(N,Abar,k)*Hfunc(nuN(M,Abar),zN(N,Abar)))
      CLRudulchi0(i,k) = CLRudulchi0(i,k) + (g2**2.d0/MN(Abar))*
     $           (CudulR(M,N,1,k)*VNLu(M,Abar,1)*VNLd(N,Abar,i)*
     $            Hfunc(nuN(M,Abar),yN(N,Abar)) + CuludL(M,N,1,i)*
     $    VNRu(M,Abar,1)*VNRl(N,Abar,k)*Hfunc(nuN(M,Abar),zN(N,Abar)))
      CRRudulchi0(i,k) = CRRudulchi0(i,k) + (g2**2.d0/MN(Abar))*
     $           (CudulR(M,N,1,k)*VNRu(M,Abar,1)*VNRd(N,Abar,i)*
     $            Hfunc(nuN(M,Abar),yN(N,Abar)) + CuludR(M,N,1,i)*
     $     VNRu(M,Abar,1)*VNRl(N,Abar,k)*Hfunc(nuN(M,Abar),zN(N,Abar)))
!      write(*,*)Hfunc(nuN(M,Abar),yN(N,Abar)),Hfunc(nuN(M,Abar),zN(N,Abar)),
!     $          Hfunc(nuN(M,Abar),yN(N,Abar)),Hfunc(nuN(M,Abar),zN(N,Abar)),
!     $          Hfunc(nuN(M,Abar),yN(N,Abar)),Hfunc(nuN(M,Abar),zN(N,Abar)),
!     $          Hfunc(nuN(M,Abar),yN(N,Abar)),Hfunc(nuN(M,Abar),zN(N,Abar))
          enddo
         enddo
         enddo
      CRLudul6d(i,k)=(ALgutmsusy*16.0*pi**2.d0*alph1gut*4*pi/MHC**2.d0
     $    *delta(i,k))+CRLudulG(i,k)+CRLudulchipm(i,k)+CRLudulchi0(i,k)
!      CLLudul(i,k) = 
      CRLudul(i,k) = CRLudulG(i,k)+CRLudulchipm(i,k)+CRLudulchi0(i,k)
!      CLRudul(i,k) =
      CRRudul(i,k) = CRRudulG(i,k)+CRRudulchipm(i,k)+CRRudulchi0(i,k)
!      write(*,*)"CRL",CRLudul(i,k),CRLudulG(i,k),CRLudulchipm(i,k),
!     $            CRLudulchi0(i,k)
      enddo
      enddo
!      write(*,*)"Cfull","gluino","chargino","higgsino"
!      write(*,*)"CRL",CRLudul(1,1),CRLudulG(1,1),CRLudulchipm(1,1),
!     $            CRLudulchi0(1,1)
!      write(*,*)"CRR",CRRudul(1,1),CRRudulG(1,1),CRRudulchipm(1,1),
!     $            CRRudulchi0(1,1)
      do i = 1,3
      do k = 1,3
         do N = 1,6
         do A =1,2
         do M = 1,6
      CLLudulchipm(i,k) = CLLudulchipm(i,k) - (g2**2.d0/MCh(A))*
     $                    CudulL(M,N,1,k)*VCLu(N,A,1)*VCLd(M,A,i)
     $                    *Hfunc(xC(M,A),uC(N,A))

      CLRudulchipm(i,k) = CLRudulchipm(i,k) - (g2**2.d0/MCh(A))*
     $                    CudulR(M,N,1,k)*VCLu(N,A,1)*VCLd(M,A,i)
     $                    *Hfunc(xC(M,A),uC(N,A))
!      write(*,*)Hfunc(xC(M,A),uC(N,A)),Hfunc(xC(M,A),uC(N,A))
!       write(*,*)M,N,A,k,i,CudulR(M,N,1,k)*VCLu(N,A,1)*VCLd(M,A,i)
!     $                    *Hfunc(xC(M,A),uC(N,A))
         enddo
         do m1 = 1,3
      CLLudulchipm(i,k) = CLLudulchipm(i,k) + (g2**2.d0/MCh(A))*
     $                    CdnuudL(N,m1,1,i)*VCLu(N,A,1)*VCLl(m1,A,k)
     $                    *Hfunc(uC(N,A),zC(m1,A))

      CLRudulchipm(i,k) = CLRudulchipm(i,k) + (g2**2.d0/MCh(A))*
     $                    CdnuudL(N,m1,1,i)*VCRu(N,A,1)*VCRl(m1,A,k)
     $                    *Hfunc(uC(N,A),zC(m1,A))
!      write(*,*)Hfunc(uC(N,A),zC(m1,A)),Hfunc(uC(N,A),zC(m1,A))

!      write(*,*)N,m1,A,k,i,(g2**2.d0/MCh(A))*
!     $                    CdnuudL(N,m1,1,i)*VCRu(N,A,1)*VCRl(m1,A,k)
!     $                    *Hfunc(uC(N,A),zC(m1,A))
         enddo
         enddo
         enddo
!         enddo
      CLLudul(i,k) = CLLudulG(i,k)+CLLudulchipm(i,k)+CLLudulchi0(i,k)
      CLRudul(i,k) = CLRudulG(i,k)+CLRudulchipm(i,k)+CLRudulchi0(i,k)

      CLRudul6d(i,k) =(ALgutmsusy*16.0*pi**2.d0*alph1gut*4*pi/MHC**2.d0
     $               *(delta(i,k)+VCKMGUT(1,i)*VCKMGUTT(k,1)))+ 
     $                CLRudulG(i,k)+CLRudulchipm(i,k)+CLRudulchi0(i,k)
!      write(*,*)CLLudul(i,k),CLRudul(i,k)
      enddo
      enddo
!      write(*,*)"CLL",CLLudul(1,1),CLLudulG(1,1),CLLudulchipm(1,1),
!     $            CLLudulchi0(1,1)
!      write(*,*)"CLR",CLRudul(1,1),CLRudulG(1,1),CLRudulchipm(1,1),
!     $            CLRudulchi0(1,1)
!      write(*,*)CudulR(6,6,1,1),VCLu(6,2,1),VCLd(6,2,2),
!     $                    Hfunc(xC(6,2),uC(6,2))
!      write(*,*)"R",CRRudul(1,1),CLLudul(1,1),CLRudul(1,1),CRLudul(1,1)
!------------------------------------------------------------------
      do i = 1,3
      do j = 1,3
      do k = 1,3
         do M = 1,6
      loopN:   do N = 1,6

!------CLLddunu and CLLuddnu for W_0 notation......
      CLLddunuG(i,j,k)=CLLddunuG(i,j,k)+(4.d0/3.d0)*(g3**2.d0/MG)* 
     $    (CddunuL(M,N,1,k)*VGLd(M,j)*VGLd(N,i)*Hfunc(xG(M),xG(N)))
      CLLuddnuG(i,j,k) = CLLuddnuG(i,j,k) + (4.d0/3.d0)*(g3**2.d0/MG)*
     $      (CuddnuL(M,N,j,k)*VGLu(M,1)*VGLd(N,i)*Hfunc(uG(M),xG(N)))

!--------goto expression---
!      CLLuddnuG(i,j,k) = CLLuddnuG(i,j,k) + (4.d0/3.d0)*(g3**2.d0/MG)* 
!     $        (CuddnuL(M,N,j,k)*VGLu(M,1)*VGLd(N,i)*Hfunc(uG(M),xG(N))
!     $       +CddunuL(M,N,1,k)*VGLd(M,j)*VGLd(N,i)*Hfunc(xG(M),xG(N)))

      CRLuddnuG(i,j,k) = CRLuddnuG(i,j,k) + (4.d0/3.d0)*(g3**2.d0/MG)*
     $         CuddnuL(M,N,j,k)*VGRu(M,1)*VGRd(N,i)*Hfunc(uG(M),xG(N))
      CRLddunuG(i,j,k) = CRLddunuG(i,j,k) + (4.d0/3.d0)*(g3**2.d0/MG)*
     $         CddunuL(M,N,1,k)*VGRd(M,i)*VGRd(N,j)*Hfunc(xG(M),xG(N))
!      write(*,*)Hfunc(uG(M),xG(N)),Hfunc(xG(M),xG(N)),Hfunc(uG(M),xG(N)),
!     $          Hfunc(xG(M),xG(N))  
!      write(*,*)CRLddunuG(i,j,k),CddunuL(M,N,1,k),Hfunc(xG(M),xG(N))
         do A = 1,2
      CLLuddnuchipm(i,j,k) = CLLuddnuchipm(i,j,k) + (g2**2.d0/MCh(A))*
     $             (-CuddnuL(M,N,j,k)*VCLu(N,A,1)*VCLd(M,A,i)*
     $             Hfunc(xC(M,A),uC(N,A))+CuludL(M,N,1,i)*VCLd(M,A,j)*
     $             VCLnu(N,A,k)*Hfunc(xC(M,A),wC(N,A)))
      CRLuddnuchipm(i,j,k) = CRLuddnuchipm(i,j,k) + (g2**2.d0/MCh(A))*
     $             (-CuddnuL(M,N,j,k)*VCRu(N,A,1)*VCRd(M,A,i)*
     $             Hfunc(xC(M,A),uC(N,A))+CuludR(M,N,1,i)*VCLd(M,A,j)*
     $             VCLnu(N,A,k)*Hfunc(xC(M,A),wC(N,A)))
!      write(*,*)VCLu(N,A,1),VCLd(M,A,i),
!     $             VCLd(M,A,j),
!     $             M,N,A,i,j,k,
!     $              CuddnuL(M,N,j,k),
!     $              CuludL(M,N,1,i)
         enddo
         do Abar = 1,4
!------CLLddunu and CLLuddnu for W_0 notation......
      CLLddunuchi0(i,j,k)=CLLddunuchi0(i,j,k)+(g2**2.d0/MN(Abar))*
     $               CddunuL(M,N,1,k)* VNLd(M,Abar,j)*VNLd(N,Abar,i)*
     $               Hfunc(yN(M,Abar),yN(N,Abar))
      CLLuddnuchi0(i,j,k) = CLLuddnuchi0(i,j,k) + (g2**2.d0/MN(Abar))*
     $              (CuddnuL(M,N,j,k)*VNLu(M,Abar,1)*VNLd(N,Abar,i)*
     $              Hfunc(nuN(M,Abar),yN(N,Abar)))

!--------goto expression---
!      CLLuddnuchi0(i,j,k) = CLLuddnuchi0(i,j,k) + (g2**2.d0/MN(Abar))*
!     $              (CuddnuL(M,N,j,k)*VNLu(M,Abar,1)*VNLd(N,Abar,i)*
!     $              Hfunc(nuN(M,Abar),yN(N,Abar)) + CddunuL(M,N,1,k)*
!     $     VNLd(M,Abar,j)*VNLd(N,Abar,i)*Hfunc(yN(M,Abar),yN(N,Abar)))

      CRLuddnuchi0(i,j,k) = CRLuddnuchi0(i,j,k) + (g2**2.d0/MN(Abar))*
     $              CuddnuL(M,N,j,k)*VNRu(M,Abar,1)*VNRd(N,Abar,i)*
     $              Hfunc(nuN(M,Abar),yN(N,Abar))
      CRLddunuchi0(i,j,k) = CRLddunuchi0(i,j,k) + (g2**2.d0/MN(Abar))*
     $              CddunuL(M,N,1,k)*VNRd(M,Abar,i)*VNRd(N,Abar,j)*
     $              Hfunc(yN(M,Abar),yN(N,Abar))
!      Write(*,*)Hfunc(nuN(M,Abar),yN(N,Abar)),Hfunc(yN(M,Abar),yN(N,Abar)),
!     $          Hfunc(nuN(M,Abar),yN(N,Abar)),Hfunc(yN(M,Abar),yN(N,Abar))                     
!      write(*,*)CRLddunuchi0(i,j,k),VNRd(M,Abar,i),VNRd(N,Abar,j),
!     $            Hfunc(yN(M,Abar),yN(N,Abar))
         enddo
         enddo loopN
         do n1 = 1,3
         do Abar = 1,4
!------CLLddunu and CLLuddnu for W_0 notation......
      CLLddunuchi0(i,j,k)=CLLddunuchi0(i,j,k)+(g2**2.d0/MN(Abar))*
     $          (CunuddL(M,n1,j,i)*VNLu(M,Abar,1)
     $           *VNLnu(n1,Abar,k)*Hfunc(nuN(M,Abar),wN(n1,Abar)))

      CLLuddnuchi0(i,j,k) = CLLuddnuchi0(i,j,k) + (g2**2.d0/MN(Abar))*
     $         (CdnuudL(M,n1,1,i)*VNLd(M,Abar,j)*VNLnu(n1,Abar,k)*
     $         Hfunc(yN(M,Abar),wN(n1,Abar)))
!--------goto expression---
!      CLLuddnuchi0(i,j,k) = CLLuddnuchi0(i,j,k) + (g2**2.d0/MN(Abar))*
!     $          (CdnuudL(M,n1,1,i)*VNLd(M,Abar,j)*VNLnu(n1,Abar,k)*
!     $   Hfunc(yN(M,Abar),wN(n1,Abar))+CunuddL(M,n1,j,i)*VNLu(M,Abar,1)
!     $           *VNLnu(n1,Abar,k)*Hfunc(nuN(M,Abar),wN(n1,Abar)))
!         write(*,*)Hfunc(yN(M,Abar),wN(n1,Abar)),Hfunc(nuN(M,Abar),wN(n1,Abar))

!         write(*,*)Hfunc(yN(M,Abar),wN(n1,Abar)),Hfunc(nuN(M,Abar),wN(n1,Abar))
         enddo
         enddo
         enddo
      CLLddunu(i,j,k) = - (CLLddunuG(i,j,k)+CLLddunuchi0(i,j,k))
      CLLuddnu(i,j,k) = CLLuddnuG(i,j,k) + CLLuddnuchipm(i,j,k) +
     $                  CLLuddnuchi0(i,j,k)
      CRLuddnu(i,j,k) = CRLuddnuG(i,j,k) + CRLuddnuchipm(i,j,k) +
     $                  CRLuddnuchi0(i,j,k)
      CRLddunu(i,j,k) = CRLddunuG(i,j,k) + CRLddunuchi0(i,j,k) 
      CRLuddnu6d(i,j,k)= (ALgutmsusy*16.0*pi**2.d0*VCKMGUT(1,j)*
     $                 delta(i,k)*alph1gut*4*pi/MHC**2.d0)+
     $     CRLuddnuG(i,j,k) + CRLuddnuchipm(i,j,k) +CRLuddnuchi0(i,j,k)
!      Write(*,*)CRLuddnu(1,2,3),CRLuddnuG(1,2,3),CRLuddnuchi0(1,2,3),
!     $       CRLuddnuchipm(1,2,3)   
!      Write(*,*)CRLddunu(i,j,k),CRLuddnu(i,j,k),CLLuddnu(i,j,k)
      enddo
      enddo
      enddo
!      Write(*,*)"CRLddunu123Gchi0",CRLddunu(1,2,3)/158.0,
!     $      CRLddunuG(1,2,3)/158.0,CRLddunuchi0(1,2,3)/158.0
!      Write(*,*)"CRLuddnu123Gchi0Chip",CRLuddnu(1,2,3)/158.0,
!     $         CRLuddnuG(1,2,3)/158.0,
!     $       CRLuddnuchi0(1,2,3)/158.0,CRLuddnuchipm(1,2,3)/158.0
!      Write(*,*)"CRLuddnu213Gchi0Chip",CRLuddnu(2,1,3)/158.0,
!     $          CRLuddnuG(2,1,3)/158.0,
!     $       CRLuddnuchi0(2,1,3)/158.0,CRLuddnuchipm(2,1,3)/158.0
!      Write(*,*)"CLLuddnu123Gchi0Chip",CLLuddnu(1,2,3)/158.0,
!     $         CLLuddnuG(1,2,3)/158.0,
!     $       CLLuddnuchi0(1,2,3)/158.0,CLLuddnuchipm(1,2,3)/158.0
!      Write(*,*)"CLLuddnu213Gchi0Chip",CLLuddnu(2,1,3)/158.0,
!     $          CLLuddnuG(2,1,3)/158.0,
!     $       CLLuddnuchi0(2,1,3)/158.0,CLLuddnuchipm(2,1,3)/158.0
!      Write(*,*)"LL111",CLLuddnu(1,1,1)*1.d0/(4.d0*pi)**2.d0,
!     $          CLLuddnuG(1,1,1)*1.d0/(4.d0*pi)**2.d0,
!     $        CLLuddnuchipm(1,1,1)*1.d0/(4.d0*pi)**2.d0
!     $        ,CLLuddnuchi0(1,1,1)*1.d0/(4.d0*pi)**2.d0
!      Write(*,*)"LL112",CLLuddnu(1,1,2)*1.d0/(4.d0*pi)**2.d0
!     $          ,CLLuddnuG(1,1,2)*1.d0/(4.d0*pi)**2.d0,
!     $        CLLuddnuchipm(1,1,2)*1.d0/(4.d0*pi)**2.d0
!     $        ,CLLuddnuchi0(1,1,2)*1.d0/(4.d0*pi)**2.d0
!      Write(*,*)"LL113",CLLuddnu(1,1,3)*1.d0/(4.d0*pi)**2.d0
!     $          ,CLLuddnuG(1,1,3)*1.d0/(4.d0*pi)**2.d0,
!     $        CLLuddnuchipm(1,1,3)*1.d0/(4.d0*pi)**2.d0
!     $        ,CLLuddnuchi0(1,1,3)*1.d0/(4.d0*pi)**2.d0
!      Write(*,*)"RL111",CRLuddnu(1,1,1)*1.d0/(4.d0*pi)**2.d0,
!     $          CRLuddnuG(1,1,1)*1.d0/(4.d0*pi)**2.d0,
!     $        CRLuddnuchipm(1,1,1)*1.d0/(4.d0*pi)**2.d0
 !    $       ,CRLuddnuchi0(1,1,1)*1.d0/(4.d0*pi)**2.d0
!      Write(*,*)"RL112",CRLuddnu(1,1,2)*1.d0/(4.d0*pi)**2.d0,
!     $        CRLuddnuG(1,1,2)*1.d0/(4.d0*pi)**2.d0,
!     $        CRLuddnuchipm(1,1,2)*1.d0/(4.d0*pi)**2.d0,
!     $        CRLuddnuchi0(1,1,2)*1.d0/(4.d0*pi)**2.d0
!      Write(*,*)"RL113",CRLuddnu(1,1,3)*1.d0/(4.d0*pi)**2.d0,
!     $        CRLuddnuG(1,1,3)*1.d0/(4.d0*pi)**2.d0,
!     $        CRLuddnuchipm(1,1,3)*1.d0/(4.d0*pi)**2.d0,
!     $        CRLuddnuchi0(1,1,3)*1.d0/(4.d0*pi)**2.d0
!---------------------------------------------------------------------
!------RG runing from Msusy to 2GeV..Using analytical equation (B.4)
! ---from shirai paper. we need alphas at msusy,mt,mb,2GeV... 
! took from my other program
!--
!      alphasmsusy = pyyg(28)*4.d0*pi
!      write(*,*)"alphasmsusy=",alphasmsusy
       do i = 1, 9
      pyye(i) = pyyg(i)*sb
       enddo 

       do i = 1, 18
      pyye(i+9) = pyyg(i+9)*cb
       enddo 

      alphamsbar3=(pyyg(28)*4.d0*pi)/(1.d0+pyyg(28))
!      write(*,*)"alphasmsusydr=",alphamsbar3
      alphamsbar2=(pyyg(29)*4.d0*pi)/
     $            (1.d0+(pyyg(29)*4.d0*pi)/(6.d0*pi))

      pyye(28)=alphamsbar3/(4.d0*pi)
      pyye(29)=alphamsbar2/(4.d0*pi)
      pyye(30)=pyyg(30)
      pyye(31)=Sqrt(pyyg(517)**2.d0+pyyg(518)**2.d0)
      alphasmt=0.10836d0
      alphasmb = 0.2171d0
      alphas2GeV = 0.30367d0
      alpha2mz=0.0352d0
!      write(*,*)"alphas2GeV=",alphas2GeV
!--------runing from msusy to mtop
!      n0 = 31
!      x2 = mtscale
!      x1 = tq0
!      h1   =  -1.d-5
!      hmin =  2.d-8
!      eps  =  1.d-6

!      call RK4ROUTINE(pyye,n0,x1,x2,eps,h1,hmin,nok,nbad,smrge,
!     .        QMSRK4,check)

!       if(check.eq.100)then
!       flags = 'variable underflow '
!        return
!       endif

!      if(maxval(pyye(1:27)).gt.dsqrt(1.d0/(4.d0*pi)))then
!         flags = "NPERTYUK"
!         return
!      endif
!-------------------------------
!      alphasmt=pyye(28)*4.d0*pi
!      write(*,*)"alphasmt=",alphasmt
!      inps7: do i = 1, 31
!      pyyb(i) = pyye(i)
!      write(*,*) pyyb(i)     
!       enddo inps7
!------------------------------- runing from mt to mb
!      n0 = 31
!      x2 = mbscale
!      x1 = mtscale
!      h1   =  -1.d-5
!      hmin =  2.d-8
!      eps  =  1.d-6

!      call RK4ROUTINE(pyyb,n0,x1,x2,eps,h1,hmin,nok,nbad,smrgemt,
!     .        QMSRK4,check)

!       if(check.eq.120)then
!       flags = 'variable underflow '
!       Write(*,*)'variable underflow '
!        return
!       endif

!      if(maxval(pyyb(1:27)).gt.dsqrt(1.d0/(4.d0*pi)))then
!         flags = "NPERTYUK"
!       Write(*,*)"NPERTYUK"
!         return
!      endif

 !-----alpha's at bottom mass      
!      alphasmb = pyyb(28)*4.d0*pi
!      write(*,*)"alphasmb=",alphasmb
!       do i = 1, 31
!      pyy2(i) = pyyb(i)
!      write(*,*)pyyb(i)
!       enddo 
!----------------runing from mb to 2GeV
!      n0 = 31
!      x2 = QCDscale
!      x1 = mbscale
!      h1   =  -1.d-5
!      hmin =  2.d-8
!      eps  =  1.d-6

!      call RK4ROUTINE(pyy2,n0,x1,x2,eps,h1,hmin,nok,nbad,smrgemt,
!     .        QMSRK4,check)

!       if(check.eq.100)then
!       flags = 'variable underflow '
!        return
!       endif

!      if(maxval(pyy2(1:27)).gt.dsqrt(1.d0/(4.d0*pi)))then
!         flags = "NPERTYUK"
!         return
!      endif
!----alpha's at 2 GeV
!      alphas2GeV = pyy2(28)*4.d0*pi
!      write(*,*)"alphas2GeV=",alphas2GeV
!-------------RG runing of wilson coiff. from msusy to 2GeV
! RG runing factor from msusy to mt,mt to mb, mb to 2 GeV for LL and RR
!      ALmusuymt0=(alphasmt/alphasmsusy)**(2.d0/7.d0)*((-28.d0*pi-
!     $          26.d0*alphasmt)/(-28.d0*pi-26.d0*alphasmsusy))**
!     $          (-79.d0/546.d0)
!      ALmtmb0=(alphasmb/alphasmt)**(6.d0/23.d0)*((4.d0*23.d0*pi+
!     $          116.d0*alphasmb)/(4.d0*23.d0*pi+116.d0*alphasmt))**
!     $          (-1375.d0/8004.d0)
!      ALmb2GeV0=(alphas2GeV/alphasmb)**(6.d0/25.d0)*((4.d0*25.d0*pi+
!     $          154.d0*alphas2GeV)/(4.d0*25.d0*pi+154.d0*alphasmb))**
!     $          (-2047.d0/11550.d0)
!      ALLLdelta0=ALmusuymt0*ALmtmb0*ALmb2GeV0
      ALLLdelta0=(alphasmt/alphasmsusy)**(2.d0/7.d0)*
     $            (alphasmb/alphasmt)**(6.d0/23.d0)*
     $            (alphas2GeV/alphasmb)**(6.d0/25.d0)*
     $            (alpha2mz/alpha2msusy)**(27.d0/28.d0)
! RG runing factor from msusy to mt,mt to mb, mb to 2 GeV for LR and RL
!      ALmusuymt=(alphasmt/alphasmsusy)**(2.d0/7.d0)*((-28.d0*pi-
!     $          26.d0*alphasmt)/(-28.d0*pi-26.d0*alphasmsusy))**
!     $          (-19.d0/91.d0)
!      ALmtmb=(alphasmb/alphasmt)**(6.d0/23.d0)*((4.d0*23.d0*pi+
!     $          116.d0*alphasmb)/(4.d0*23.d0*pi+116.d0*alphasmt))**
!     $          (-430.d0/2001.d0)
!      ALmb2GeV=(alphas2GeV/alphasmb)**(6.d0/25.d0)*((4.d0*25.d0*pi+
!     $          154.d0*alphas2GeV)/(4.d0*25.d0*pi+154.d0*alphasmb))**
!     $          (-173.d0/825.d0)
!      ALRLdelta=ALmusuymt*ALmtmb*ALmb2GeV
      ALRLdelta=(alphasmt/alphasmsusy)**(2.d0/7.d0)*
     $            (alphasmb/alphasmt)**(6.d0/23.d0)*
     $            (alphas2GeV/alphasmb)**(6.d0/25.d0)*
     $            (alpha2mz/alpha2msusy)**(27.d0/28.d0)
      ALRLdelta6d=ALgutmsusy*ALRLdelta
       

!      Write(*,*)"ALLLdelta0",ALLLdelta0,ALRLdelta
!----------------------------------------------------------------------
! wilson coiff. at 2 GeV-----------
      do i=1,3
      do j=1,3
      do k=1,3
      CRLuddnu26d(i,j,k)=0.0
      CLLuddnu2(i,j,k)=0.0
      CLLddunu2(i,j,k)=0.0
      CRLuddnu2(i,j,k)=0.0
      CRLddunu2(i,j,k)=0.0
      enddo
      enddo
      enddo

      do i=1,3
      do k=1,3
       CLRudul26d(i,k)=0.0
       CRLudul26d(i,k)=0.0
       CLLudul2(i,k)=0.0
       CLRudul2(i,k)=0.0
       CRLudul2(i,k)=0.0
       CRRudul2(i,k)=0.0
      enddo
      enddo

      do i=1,3
      do j=1,3
      do k=1,3
      CRLuddnu26d(i,j,k)=ALRLdelta*CRLuddnu6d(i,j,k)*1.d0/
     $                  ((4.d0*pi)**2.d0)

      CLLddunu2(i,j,k)=ALLLdelta0*CLLddunu(i,j,k)*1.d0/
     $                 ((4.d0*pi)**2.d0)
      CLLuddnu2(i,j,k)=ALLLdelta0*CLLuddnu(i,j,k)*1.d0/
     $                 ((4.d0*pi)**2.d0)
      CRLuddnu2(i,j,k)=ALRLdelta*CRLuddnu(i,j,k)*1.d0/
     $                 ((4.d0*pi)**2.d0)
      CRLddunu2(i,j,k)=ALRLdelta*CRLddunu(i,j,k)*1.d0/
     $                 ((4.d0*pi)**2.d0)
!      Write(*,*)CRLuddnu2(i,j,k),ALRLdelta,CRLuddnu(i,j,k),
!     $          ALLLdelta0,CLLuddnu(i,j,k)        
!      Write(*,*)CRLddunu2(i,j,k),ALRLdelta,CRLddunu(i,j,k),
!     $          CLLuddnu2(i,j,k)
      enddo
      enddo
      enddo


!      write(1001,*)CRLddunu(1,2,3)/(158.0*total),CRLddunuG(1,2,3)/
!     $      (158.0*total),CRLuddnu(1,2,3)/(158.0*total),
!     $       CRLuddnuG(1,2,3)/(158.0*total),CRLuddnuchipm(1,2,3)/
!     $      (158.0*total),CRLuddnu(2,1,3)/(158.0*total),
!     $       CRLuddnuG(2,1,3)/(158.0*total),CRLuddnuchipm(2,1,3)/
!     $      (158.0*total),CLLuddnu(1,2,3)/(158.0*total),
!     $       CLLuddnuG(1,2,3)/(158.0*total),CLLuddnuchipm(1,2,3)/
!     $      (158.0*total),CLLuddnu(2,1,3)/(158.0*total),
!     $       CLLuddnuG(2,1,3)/(158.0*total),CLLuddnuchipm(2,1,3)/
!     $      (158.0*total)
      do i=1,3
      do k=1,3
       CLRudul26d(i,k)=(ALRLdelta*CLRudul6d(i,k))/((4.d0*pi)**2.d0)
       CRLudul26d(i,k)=(ALRLdelta*CRLudul6d(i,k))/((4.d0*pi)**2.d0)

       CLLudul2(i,k)=ALLLdelta0*CLLudul(i,k)*1.d0/((4.d0*pi)**2.d0)
       CLRudul2(i,k)=ALRLdelta*CLRudul(i,k)*1.d0/((4.d0*pi)**2.d0)
       CRLudul2(i,k)=ALRLdelta*CRLudul(i,k)*1.d0/((4.d0*pi)**2.d0)
       CRRudul2(i,k)=ALLLdelta0*CRRudul(i,k)*1.d0/((4.d0*pi)**2.d0)
!      Write(*,*)ALLLdelta0,ALRLdelta
      enddo
      enddo
!      write(*,*)"2 GeV",CLLudul2(1,1),CRLudul2(1,1),CLRudul2(1,1),
!     $          CRRudul2(1,1)
!----------Partial decay widths--------------------------------------
!first----Proton decay into charged kaon and neutrino-------
!------Matching condition at Electroweak Scale-----
!-----------------------------------------------------------------
      Mp = 0.9382d0
      Mpion0 = 0.1349d0      
      Mpionp = 0.1395d0
      MEta0 = 0.5478d0
      Mk0 = 0.4976d0
      Mkp = 0.4936d0
      mneu = 0.94d0
      MB = 1.15d0
      Fvalue =0.47d0
      Dvalue = 0.80d0
      Fpi = 0.131d0
      alp = -0.0144d0
      bep = 0.0144d0

      PiRL0 = -0.131
      PiLL0 =  0.134
      PiRLp = -0.186
      PiLLp =  0.189
      EtaRL =  0.006
      EtaLL =  0.113
      KRL0  =  0.103
      KLL0  =  0.057
 !--Kaon +
      KsRLp = -0.049
      KsLLp =  0.041
      KdRLp = -0.134 
      kdLLp =  0.139
      KdsRLp = -0.054
      KdsLLp = -0.098

!------------------------------------------------------------------
      decay: do k=1,3
!----------Partial decay width of p goes to Pion0 and charged lepton---
      ALppil(k)=PiRL0*CRLudul2(1,k) + PiLL0*CLLudul2(1,k)
      ARppil(k)=PiRL0*CLRudul2(1,k) + PiLL0*CRRudul2(1,k)
      decayppil(k)=(Mp/(32.d0*pi))*(1.d0-Mpion0**2.d0/Mp**2.d0)**2.d0*
     $            (ALppil(k)**2.d0+ARppil(k)**2.d0)

!----------Partial decay width of p goes to Eta0 and charged lepton----
      ALpetal(k)=CRLudul2(1,k)*EtaRL+CLLudul2(1,k)*EtaLL
      ARpetal(k)=CLRudul2(1,k)*EtaRL+CRRudul2(1,k)*EtaLL
      decaypetal(k)=(Mp/(32.d0*pi))*(1.d0-MEta0**2.d0/Mp**2.d0)**2.d0*
     $           (ALpetal(k)**2.d0+ARpetal(k)**2.d0)

!----------Partial decay width of p goes to K0 and charged lepton------

       ALpkl(k)=KRL0*CRLudul2(2,k) + KLL0*CLLudul2(2,k)
       ARpkl(k)=KRL0*CLRudul2(2,k) + KLL0*CRRudul2(2,k)

      decaypkl(k)=(Mp/(32.d0*pi))*(1.d0-Mk0**2.d0/Mp**2.d0)**2.d0*
     $           (ALpkl(k)**2.d0+ARpkl(k)**2.d0)
!----------Partial decay width of p goes to Pion+ and anti neutrino----
      ALppipnu(k)=CRLuddnu2(1,1,k)*PiRLp+ CLLuddnu2(1,1,k)*PiLLp

      decayppipnu(k)=(Mp/(32.d0*pi))*(1.d0-Mpionp**2.d0/Mp**2.d0)**2.d0
     $               *(ALppipnu(k)**2.d0)
!----------Partial decay width of p goes to K+ and anti neutrino-------
      ALpkpnu(k)= KdsRLp*CRLddunu2(1,2,k) + KdRLp*CRLuddnu2(1,2,k) +
     $     KsRLp*CRLuddnu2(2,1,k) + KdLLp*CLLuddnu2(1,2,k) + 
     $     KsLLp* CLLuddnu2(2,1,k)+CLLddunu2(1,2,k)*KdsLLp

      decaypkpnu(k)=(Mp/(32.d0*pi))*(1.d0-Mkp**2.d0/Mp**2.d0)**2.d0
     $               *(ALpkpnu(k)**2.d0)

       halfppil(k)=2.1126d0*10**(-32.d0)/decayppil(k)
       halfpetal(k)=2.1126d0*10**(-32.d0)/decaypetal(k)
       halfpkl(k)=2.1126d0*10**(-32.d0)/decaypkl(k)
       halfpipnu(k)=2.1126d0*10**(-32.d0)/decayppipnu(k)
       halfpkpnu(k)=2.1126d0*10**(-32.d0)/decaypkpnu(k)

       enddo decay

       halfpipnufull=2.1126d0*10**(-32.d0)/(decayppipnu(1)+
     $             decayppipnu(2)+decayppipnu(3))
       halfpkpnufull=2.1126d0*10**(-32.d0)/(decaypkpnu(1)+
     $             decaypkpnu(2)+decaypkpnu(3))

!=======================================================================================
!-------with 6-D operator---------------
!------------------------------------------------------------------
      decay6d: do k=1,3
!----------Partial decay width of p goes to Pion0 and charged lepton---
      ALppil6d(k)=PiRL0*CRLudul26d(1,k) + PiLL0*CLLudul2(1,k)
      ARppil6d(k)=PiRL0*CLRudul26d(1,k) + PiLL0*CRRudul2(1,k)
      decayppil6d(k)=(ALppil6d(k)**2.d0+ARppil6d(k)**2.d0)*
     $      (Mp/(32.d0*pi))*(1.d0-Mpion0**2.d0/Mp**2.d0)**2.d0
              

!----------Partial decay width of p goes to Eta0 and charged lepton----
      ALpetal6d(k)=CRLudul26d(1,k)*EtaRL+CLLudul2(1,k)*EtaLL
      ARpetal6d(k)=CLRudul26d(1,k)*EtaRL+CRRudul2(1,k)*EtaLL
      decaypetal6d(k)=(ALpetal6d(k)**2.d0+ARpetal6d(k)**2.d0)*
     $     (Mp/(32.d0*pi))*(1.d0-MEta0**2.d0/Mp**2.d0)**2.d0      

!----------Partial decay width of p goes to K0 and charged lepton------
      ALpkl6d(k)=KRL0*CRLudul26d(2,k)+KLL0*CLLudul2(2,k)
      ARpkl6d(k)=KRL0*CLRudul26d(2,k)+KLL0*CRRudul2(2,k)

      decaypkl6d(k)=(Mp/(32.d0*pi))*(1.d0-Mk0**2.d0/Mp**2.d0)**2.d0*
     $           (ALpkl6d(k)**2.d0+ARpkl6d(k)**2.d0)
!----------Partial decay width of p goes to Pion+ and anti neutrino----
      ALppipnu6d(k)=CRLuddnu26d(1,1,k)*PiRLp+ CLLuddnu2(1,1,k)*PiLLp

      decayppipnu6d(k)=(ALppipnu6d(k)**2.d0)*(Mp/(32.d0*pi))*
     $               (1.d0-Mpionp**2.d0/Mp**2.d0)**2.d0
!----------Partial decay width of p goes to K+ and anti neutrino-------
      ALpkpnu6d(k)= KdsRLp*CRLddunu2(1,2,k)+KdRLp*CRLuddnu26d(1,2,k)+
     $     KsRLp*CRLuddnu26d(2,1,k) + KdLLp*CLLuddnu2(1,2,k) + 
     $     KsLLp* CLLuddnu2(2,1,k)+CLLddunu2(1,2,k)*KdsLLp

      decaypkpnu6d(k)=(Mp/(32.d0*pi))*(1.d0-Mkp**2.d0/Mp**2.d0)**2.d0
     $               *(ALpkpnu6d(k)**2.d0)

       halfppil6d(k)=2.1126d0*10**(-32.d0)/decayppil6d(k)
       halfpetal6d(k)=2.1126d0*10**(-32.d0)/decaypetal6d(k)
       halfpkl6d(k)=2.1126d0*10**(-32.d0)/decaypkl6d(k)
       halfpipnu6d(k)=2.1126d0*10**(-32.d0)/decayppipnu6d(k)
       halfpkpnu6d(k)=2.1126d0*10**(-32.d0)/decaypkpnu6d(k)

       enddo decay6d

      halfpipnu6dfull=2.1126d0*10**(-32.d0)/(decayppipnu6d(1)+
     $                decayppipnu6d(2)+decayppipnu6d(3))
      halfpkpnu6dfull=2.1126d0*10**(-32.d0)/(decaypkpnu6d(1)+
     $                 decaypkpnu6d(2)+decaypkpnu6d(3))
!=======================================================================================

!--------------------only 6-D contribution----------------------------------------
      decayppiepure6d=pi*mp/2.d0*(1.d0-Mpion0**2.d0/Mp**2.d0)**2.d0*
     $                ALRLdelta6d**2.d0*alph1gut**2.d0/MHC**4.d0*
     $             (1.d0+(1.d0+VCKMGUT(1,1)**2.d0)**2.d0)*PiRL0**2.d0

      decayppimupure6d=pi*mp/2.d0*(1.d0-Mpion0**2.d0/Mp**2.d0)**2.d0*
     $                ALRLdelta6d**2.d0*alph1gut**2.d0/MHC**4.d0*
     $             (1.d0+VCKMGUT(1,1)*VCKMGUTT(2,1))**2.d0*PiRL0**2.d0

      decayppipnuepure6d=pi*mp/2.d0*(1.d0-Mpionp**2.d0/Mp**2.d0)**2.d0*
     $                ALRLdelta6d**2.d0*alph1gut**2.d0/MHC**4.d0*
     $             (VCKMGUT(1,1)*PiRLp)**2.d0  
 !     decayppipnuepure6d=pi*mp/2.d0*(1.d0-Mpionp**2.d0/Mp**2.d0)**2.d0*
!     $                ALRLdelta6d**2.d0*alph1gut**2.d0/MHC**4.d0*
!     $             (VCKMGUT(1,1)*PiRLp)**2.d0  
      decaypknuepure6d=pi*mp/2.d0*(1.d0-Mkp**2.d0/Mp**2.d0)**2.d0*
     $                ALRLdelta6d**2.d0*alph1gut**2.d0/MHC**4.d0*
     $             (VCKMGUT(1,2)*KdRLp)**2.d0 
      decaypknumupure6d=pi*mp/2.d0*(1.d0-Mkp**2.d0/Mp**2.d0)**2.d0*
     $                ALRLdelta6d**2.d0*alph1gut**2.d0/MHC**4.d0*
     $             (VCKMGUT(1,1)*KsRLp)**2.d0 
      decaypknupure6d=pi*mp/2.d0*(1.d0-Mkp**2.d0/Mp**2.d0)**2.d0*
     $                ALRLdelta6d**2.d0*alph1gut**2.d0/MHC**4.d0*
     $        ((VCKMGUT(1,1)*KsRLp)**2.d0+(VCKMGUT(1,2)*KdRLp)**2.d0)
      decaypk0epure6d=pi*mp/2.d0*(1.d0-Mk0**2.d0/Mp**2.d0)**2.d0*
     $                ALRLdelta6d**2.d0*alph1gut**2.d0/MHC**4.d0*
     $        (VCKMGUT(1,2)*VCKMGUTT(1,1)*KRL0)**2.d0
      decaypk0mupure6d=pi*mp/2.d0*(1.d0-Mk0**2.d0/Mp**2.d0)**2.d0*
     $                ALRLdelta6d**2.d0*alph1gut**2.d0/MHC**4.d0*
     $        (1.d0+(1.d0+VCKMGUT(1,2)*VCKMGUTT(2,1))**2.d0)*KRL0**2.d0

      halfppiepure6d=2.1126d0*10**(-32.d0)/decayppiepure6d
      halfppimupure6d=2.1126d0*10**(-32.d0)/decayppimupure6d
      halfppipnuepure6d=2.1126d0*10**(-32.d0)/decayppipnuepure6d
      halfknuepure6d=2.1126d0*10**(-32.d0)/decaypknuepure6d
      halfpknumupure6d=2.1126d0*10**(-32.d0)/decaypknumupure6d
      halfpknupure6d=2.1126d0*10**(-32.d0)/decaypknupure6d
      halfpk0epure6d=2.1126d0*10**(-32.d0)/decaypk0epure6d
      halfpk0mupure6d=2.1126d0*10**(-32.d0)/decaypk0mupure6d











   
!      decay: do k=1,3
!!----------Partial decay width of p goes to Pion0 and charged lepton---
!      ALppil(k)=(1.d0/sqrt(2.d0))*(1+Fvalue+Dvalue)*(alp*CRLudul2(1,k)+
!     $           bep*CLLudul2(1,k))
!      ARppil(k)=-(1.d0/sqrt(2.d0))*(1+Fvalue+Dvalue)*(alp*CLRudul2(1,k)
!     $           +bep*CRRudul2(1,k))

!      decayppil(k)=(Mp/(32.d0*pi))*(1.d0-Mpion0**2.d0/Mp**2.d0)**2.d0*
!     $             (1.d0/(Fpi**2.d0))*(ALppil(k)**2.d0+ARppil(k)**2.d0)
!!----------Partial decay width of p goes to Eta0 and charged lepton----
!      ALpetal(k)=(sqrt(3.d0/2.d0))*((-1.d0/3.d0+Fvalue-Dvalue/3.d0)*
!     $           alp*CRLudul2(1,k) + (1.d0+Fvalue-Dvalue*1.d0/3.d0)*
!     $           bep*CLLudul2(1,k))
!      ARpetal(k)=-(sqrt(3.d0/2.d0))*((-1.d0/3.d0+Fvalue-Dvalue/3.d0)*
!     $           alp*CLRudul2(1,k) + (1.d0+Fvalue-Dvalue*1.d0/3.d0)*
!     $           bep*CRRudul2(1,k))

!      decaypetal(k)=(Mp/(32.d0*pi))*(1.d0-MEta0**2.d0/Mp**2.d0)**2.d0*
!     $           (1.d0/(Fpi**2.d0))*(ALpetal(k)**2.d0+ARpetal(k)**2.d0)
!!----------Partial decay width of p goes to K0 and charged lepton------
!      ALpkl(k)=((-1.d0+(mneu/MB)*(Fvalue-Dvalue))*alp*CRLudul2(2,k)+
!     $         (1.d0+(mneu/MB)*(Fvalue-Dvalue))*bep*CLLudul2(2,k))
!      ARpkl(k)=-((-1.d0+(mneu/MB)*(Fvalue-Dvalue))*alp*CLRudul2(2,k)+
!     $         (1.d0+(mneu/MB)*(Fvalue-Dvalue))*bep*CRRudul2(2,k))

!      decaypkl(k)=(Mp/(32.d0*pi))*(1.d0-Mk0**2.d0/Mp**2.d0)**2.d0*
!     $           (1.d0/(Fpi**2.d0))*(ALpkl(k)**2.d0+ARpkl(k)**2.d0)

!!----------Partial decay width of p goes to Pion+ and anti neutrino----
!      ALppipnu(k)=(1+Fvalue+Dvalue)*(alp*CRLuddnu2(1,1,k)+bep*
!     $             CLLuddnu2(1,1,k))

!      decayppipnu(k)=(Mp/(32.d0*pi))*(1.d0-Mpionp**2.d0/Mp**2.d0)**2.d0
!     $               *(1.d0/(Fpi**2.d0))*(ALppipnu(k)**2.d0)
!!----------Partial decay width of p goes to K+ and anti neutrino-------
!      ALpkpnu(k)=(((1.d0-(mneu/MB)*(Fvalue-Dvalue/3.d0))*alp*
!     $     CRLddunu2(1,2,k)) + ((1.d0+(mneu/MB)*(Fvalue+Dvalue/3.d0))*
!     $     alp*CRLuddnu2(1,2,k)) + ((Dvalue*2.d0*mneu/(3.d0*MB))*alp*
!     $     CRLuddnu2(2,1,k)) + ((1.d0+(mneu/MB)*(Fvalue+Dvalue/3.d0))*
!     $     bep*CLLuddnu2(1,2,k)) + ((Dvalue*2.d0*mneu/(3.d0*MB))*
!     $     bep*CLLuddnu2(2,1,k)))

!      decaypkpnu(k)=(Mp/(32.d0*pi))*(1.d0-Mkp**2.d0/Mp**2.d0)**2.d0
!     $               *(1.d0/(Fpi**2.d0))*(ALpkpnu(k)**2.d0)
!       halfppil(k)=2.1126d0*10**(-32.d0)/decayppil(k)
!       halfpetal(k)=2.1126d0*10**(-32.d0)/decaypetal(k)
!       halfpkl(k)=2.1126d0*10**(-32.d0)/decaypkl(k)
!       halfpipnu(k)=2.1126d0*10**(-32.d0)/decayppipnu(k)
!       halfpkpnu(k)=2.1126d0*10**(-32.d0)/decaypkpnu(k)

!!       write(*,*)"halflife",halfppil(k),halfpetal(k),halfpkl(k),
!!     $           halfpipnu(k),halfpkpnu(k)
!!       write(*,*)decayppil(k),decaypetal(k),decaypkl(k),decayppipnu(k),
!!     $           decaypkpnu(k)

!       enddo decay

      tautotal=ALpkpnu(3)
      etotal=ALpkpnu(1)

      CRLddunu123=ALRLdelta*CRLddunu(1,2,3)/((4.d0*pi)**2.d0*tautotal)
      CRLddunuG123=ALRLdelta*CRLddunuG(1,2,3)/((4.d0*pi)**2.d0
     $             *tautotal)
      CRLddunuchi0123=ALRLdelta*CRLddunuchi0(1,2,3)/
     $               ((4.d0*pi)**2.d0*tautotal)
      CRLuddnu123=ALRLdelta*CRLuddnu(1,2,3)/((4.d0*pi)**2.d0*tautotal)
      CRLuddnuG123=ALRLdelta*CRLuddnuG(1,2,3)/((4.d0*pi)**2.d0
     $            *tautotal)
      CRLuddnuchipm123=ALRLdelta*CRLuddnuchipm(1,2,3)/((4.d0*pi)**2.d0
     $                *tautotal)
      CRLuddnuchi0123=ALRLdelta*CRLuddnuchi0(1,2,3)/((4.d0*pi)**2.d0
     $                *tautotal)
      CLLuddnu123=ALLLdelta0*CLLuddnu(1,2,3)/((4.d0*pi)**2.d0*tautotal)
      CLLuddnuG123=ALLLdelta0*CLLuddnuG(1,2,3)/((4.d0*pi)**2.d0
     $             *tautotal)
      CLLuddnuchipm123=ALLLdelta0*CLLuddnuchipm(1,2,3)/((4.d0*pi)**2.d0
     $                  *tautotal)
      CLLuddnuchi0123=ALLLdelta0*CLLuddnuchi0(1,2,3)/((4.d0*pi)**2.d0
     $                 *tautotal)

      CRLuddnu213=ALRLdelta*CRLuddnu(2,1,3)/((4.d0*pi)**2.d0*tautotal)
      CRLuddnuG213=ALRLdelta*CRLuddnuG(2,1,3)/((4.d0*pi)**2.d0
     $            *tautotal)
      CRLuddnuchipm213=ALRLdelta*CRLuddnuchipm(2,1,3)/((4.d0*pi)**2.d0
     $                  *tautotal)
      CRLuddnuchi0213=ALRLdelta*CRLuddnuchi0(2,1,3)/((4.d0*pi)**2.d0
     $                *tautotal)
      CLLuddnu213=ALLLdelta0*CLLuddnu(2,1,3)/((4.d0*pi)**2.d0*tautotal)
      CLLuddnuG213=ALLLdelta0*CLLuddnuG(2,1,3)/((4.d0*pi)**2.d0
     $             *tautotal)
      CLLuddnuchipm213=ALLLdelta0*CLLuddnuchipm(2,1,3)/((4.d0*pi)**2.d0
     $                *tautotal)
      CLLuddnuchi0213=ALLLdelta0*CLLuddnuchi0(2,1,3)/((4.d0*pi)**2.d0
     $               *tautotal)

      CRLddunu121=ALRLdelta*CRLddunu(1,2,1)/((4.d0*pi)**2.d0*etotal)
      CRLddunuG121=ALRLdelta*CRLddunuG(1,2,1)/((4.d0*pi)**2.d0*etotal)   
      CRLddunuchi0121=ALRLdelta*CRLddunuchi0(1,2,1)/((4.d0*pi)**2.d0
     $                *etotal)
      CRLuddnu121=ALRLdelta*CRLuddnu(1,2,1)/((4.d0*pi)**2.d0*etotal)
      CRLuddnuG121=ALRLdelta*CRLuddnuG(1,2,1)/((4.d0*pi)**2.d0*etotal)
      CRLuddnuchipm121=ALRLdelta*CRLuddnuchipm(1,2,1)/((4.d0*pi)**2.d0
     $                 *etotal)
      CRLuddnuchi0121=ALRLdelta*CRLuddnuchi0(1,2,1)/((4.d0*pi)**2.d0
     $                *etotal)
      CLLuddnu121=ALLLdelta0*CLLuddnu(1,2,1)/((4.d0*pi)**2.d0*etotal)
      CLLuddnuG121=ALLLdelta0*CLLuddnuG(1,2,1)/((4.d0*pi)**2.d0*etotal)
      CLLuddnuchipm121=ALLLdelta0*CLLuddnuchipm(1,2,1)/((4.d0*pi)**2.d0
     $               *etotal)
      CLLuddnuchi0121=ALLLdelta0*CLLuddnuchi0(1,2,1)/((4.d0*pi)**2.d0
     $               *etotal)

      CRLuddnu211=ALRLdelta*CRLuddnu(2,1,1)/((4.d0*pi)**2.d0*etotal)
      CRLuddnuG211=ALRLdelta*CRLuddnuG(2,1,1)/((4.d0*pi)**2.d0*etotal)
      CRLuddnuchipm211=ALRLdelta*CRLuddnuchipm(2,1,1)/((4.d0*pi)**2.d0
     $                 *etotal)
      CRLuddnuchi0211=ALRLdelta*CRLuddnuchi0(2,1,1)/((4.d0*pi)**2.d0
     $                *etotal)
      CLLuddnu211=ALLLdelta0*CLLuddnu(2,1,1)/((4.d0*pi)**2.d0*etotal)
      CLLuddnuG211=ALLLdelta0*CLLuddnuG(2,1,1)/((4.d0*pi)**2.d0*etotal)
      CLLuddnuchipm211=ALLLdelta0*CLLuddnuchipm(2,1,1)/((4.d0*pi)**2.d0
     $                 *etotal)
      CLLuddnuchi0211=ALLLdelta0*CLLuddnuchi0(2,1,1)/((4.d0*pi)**2.d0
     $                *etotal)


      ALppil1=ALppil(1)
      ARppil1=ARppil(1)

      CRLudul11=ALRLdelta*CRLudul(1,1)/((4.d0*pi)**2.d0*ALppil1)
      CRLudulG11=ALRLdelta*CRLudulG(1,1)/((4.d0*pi)**2.d0*ALppil1)
      CRLudulchipm11=ALRLdelta*CRLudulchipm(1,1)/((4.d0*pi)**2.d0
     $               *ALppil1)
      CRLudulchi011=ALRLdelta*CRLudulchi0(1,1)/((4.d0*pi)**2.d0
     $              *ALppil1)

      CLLudul11=ALLLdelta0*CLLudul(1,1)/((4.d0*pi)**2.d0*ALppil1)
      CLLudulG11=ALLLdelta0*CLLudulG(1,1)/((4.d0*pi)**2.d0*ALppil1)
      CLLudulchipm11=ALLLdelta0*CLLudulchipm(1,1)/((4.d0*pi)**2.d0
     $               *ALppil1)
      CLLudulchi011=ALLLdelta0*CLLudulchi0(1,1)/((4.d0*pi)**2.d0
     $              *ALppil1)

      CRRudul11=ALLLdelta0*CRRudul(1,1)/((4.d0*pi)**2.d0*ARppil1)
      CRRudulG11=ALLLdelta0*CRRudulG(1,1)/((4.d0*pi)**2.d0*ARppil1)
      CRRudulchipm11=ALLLdelta0*CRRudulchipm(1,1)/((4.d0*pi)**2.d0
     $               *ARppil1)
      CRRudulchi011=ALLLdelta0*CRRudulchi0(1,1)/((4.d0*pi)**2.d0
     $              *ARppil1)

      CLRudul11=ALRLdelta*CLRudul(1,1)/((4.d0*pi)**2.d0*ARppil1)
      CLRudulG11=ALRLdelta*CLRudulG(1,1)/((4.d0*pi)**2.d0*ARppil1)
      CLRudulchipm11=ALRLdelta*CLRudulchipm(1,1)/((4.d0*pi)**2.d0
     $               *ARppil1)
      CLRudulchi011=ALRLdelta*CLRudulchi0(1,1)/((4.d0*pi)**2.d0
     $              *ARppil1)
      
!       write(*,*)"pi0",halfppil(1),ALppil(1),ARppil(1),decayppil(1)
!      return 

!      Close(1001)
      end subroutine Proton_goto

!---------------------------------------------------
! Function Hfunc(x,y) in different limits-------------
!---------------------------------------------------
      double precision function Hfunc(x,y)
      implicit none
      double precision x,y
       If(abs(x-y)<10**(-10.d0))then
       Hfunc = ((y-Log(y)-1.d0)/((y-1.d0)**2.d0))
         If(abs(y-1.d0)<10**(-10.d0))then
         Hfunc = 0.5d0
         endif
       else if(abs(x-1.d0)<10**(-10.d0))then
        Hfunc = (-y+y*Log(y)+1.d0)/((y-1)**2.d0)
         If(abs(y-1.d0)<10**(-10.d0))then
         Hfunc = 0.5d0
         endif  
       else if(abs(y-1.d0)<10**(-10.d0))then
        Hfunc = (-x+x*Log(x)+1.d0)/((x-1)**2.d0)
         If(abs(x-1.d0)<10**(-10.d0))then
         Hfunc = 0.5d0
         endif  
       else
       Hfunc = (1.d0/(x-y))*((x*Log(x))/(x-1.d0)-(y*Log(y))/(y-1.d0))
      endif
!       write(*,*)Hfunc,x,y
        end function 

!------(x4=alphas at 2 GeV, x5=alphas at bottom mass x6=alphamz)
!-----------this funtion for LL, RR delta=0
      double precision function ALRR(x4,x5,x6)
      implicit none
      double precision a1,a2,a3,pi
      double precision x4,x5,x6

      pi = 4.d0*datan(1.d0) 
        
         a1 = ((x4/x5)**(6.d0/25.d0))*((x5/x6)**(6.d0/23.d0))
         a2 = ((x4+(50.d0*pi/77.d0))/(x5+(50.d0*pi/77.d0)))**
     $        (-2047.d0/11550.d0)
         a3 = ((x5+(23.d0*pi/29.d0))/(x6+(23.d0*pi/29.d0)))**
     $        (-1375.d0/8004.d0)
      
        ALRR = a1*a2*a3
      
        end function 
!------(x4=alphas at 2 GeV, x5=alphas at bottom mass x6=alphamz)
!-----------this funtion for LR delta=-10/3
      double precision function ALLR(x4,x5,x6)
      implicit none
      double precision a1,a2,a3,pi
      double precision x4,x5,x6

      pi = 4.d0*datan(1.d0)

         a1 = ((x4/x5)**(6.d0/25.d0))*((x5/x6)**(6.d0/23.d0))
         a2 = ((x4+(50.d0*pi/77.d0))/(x5+(50.d0*pi/77.d0)))**
     $        (-173.d0/825.d0)
         a3 = ((x5+(23.d0*pi/29.d0))/(x6+(23.d0*pi/29.d0)))**
     $        (-430.d0/2001.d0)
      
        ALLR = a1*a2*a3
      
        end function 
!============================================================================
!------------------------------------------------------------------------
!                              ProtonRGE BEGINS
!------------------------------------------------------------------------
!---------In tihs program we take care of C5L and C5R rges of shirai. we used 
!  indices of shirai that is different for goto(in lepton yukawa) and we also
!  take care of -1/2 factor of sudhir code defination of log(mu).
!  here I am taking 2/5 and 12/5 coeff. of g1 square but in goto paper
!  it is 2/3 and 12/3 for C5L,C5R respectively.   
****f* susyflav/proton decay rge.f 
*  NAME
*    mssmrge
*  SYNOPSIS
*     In this subroutine I write all the proton decay wilson coefficent renormalization group equations, 
*     from GUT to 2 GeV scale. 
*  INPUTS
*     pyy(518)    - Initial values for all RGEs
*     t          - energy scale 
*  RESULT
*     pyy(518)    - RGE output at a scale t
*  EXAMPLE
*     subroutine mssmrge(t,pyy,dydx)
*  NOTES
*     The notation I use follows Natsumi Nagata and Santoshi
*     Shirai paper(arXiv :1312.7854v2) 
*     Note that alpha1, alpha2, alpha3 in the following are 
*     further normalised by a "4 pi " factor. 
*     Remember all the yukawa matrices are also normalised by this "4 pi"
*     factor. thus : yu(i,j) = [ 1/(4 pi) ] (yu(i,j) ); where yu = yukawa
*     in the lagrangian and yu == yukawa in the program. 
*	
*	dydx(1)-dydx(3) : yu(1,1) - yu(1,3) 
*	dydx(4)-dydx(6) : yu(2,1) - yu(2,3) 
*	dydx(7)-dydx(9) : yu(3,1) - yu(3,3) 
*
*	dydx(10)-dydx(12) : yd(1,1) - yd(1,3) 
*	dydx(13)-dydx(15) : yd(2,1) - yd(2,3) 
*	dydx(16)-dydx(18) : yd(3,1) - yd(3,3) 
*
*	dydx(19)-dydx(21) : ye(1,1) - ye(1,3) 
*	dydx(22)-dydx(24) : ye(2,1) - ye(2,3) 
*	dydx(25)-dydx(27) : ye(3,1) - ye(3,3) 
*       dydx(28)-dydx(30) : aplh3-alph1 

*	dydx(31)-dydx(111) : C5L(1,1,1,1) - C5L(3,3,3,3) 
*	dydx(112)-dydx(192) : C5R(1,1,1,1) - C5R(3,3,3,3)
*
*	dydx(193)-dydx(273) : C1(1,1,1,1) - C1(3,3,3,3)
*	dydx(274)-dydx(354) : C2(1,1,1,1) - C2(3,3,3,3) 
*	dydx(355)-dydx(435) : C3(1,1,1,1) - C3(3,3,3,3) 
*	dydx(436)-dydx(516) : C4(1,1,1,1) - C4(3,3,3,3) 
*	 	 
*       dydx(517)-dydx(518) : vev1-vev2
*    
C       Modified date 27th may 2018
c================================================================================= 
      subroutine protonrge(t,pyy,dydx)
 
      implicit none

      integer i, j, k, l, i0, q, c, ip 

      DOUBLE PRECISION pi, b1, b2, b3, MX, yu(3,3), yd(3,3), ye(3,3)
      DOUBLE PRECISION C5L1(3,3,3,3), C5R1(3,3,3,3), C11(3,3,3,3)
      DOUBLE PRECISION C21(3,3,3,3), C31(3,3,3,3),C41(3,3,3,3),ynu(3,3)
      DOUBLE PRECISION gauge5L(3,3,3,3), gauge5R(3,3,3,3), pyy(518)

      DOUBLE PRECISION yukawai5L(3,3,3,3), yukawaj5L(3,3,3,3)
      DOUBLE PRECISION yukawak5L(3,3,3,3), yukawal5L(3,3,3,3)
      DOUBLE PRECISION yukawai5R(3,3,3,3), yukawaj5R(3,3,3,3)
      DOUBLE PRECISION yukawak5R(3,3,3,3), yukawal5R(3,3,3,3)
      DOUBLE PRECISION alph3, alph2, alph1, vev1, vev2

      DOUBLE PRECISION yudag(3,3), yuyudag(3,3), tryuyudag
      DOUBLE PRECISION yddag(3,3), ydyddag(3,3), trydyddag
      DOUBLE PRECISION yedag(3,3), yeyedag(3,3), tryeyedag

      DOUBLE PRECISION ynudag(3,3), ynuynudag(3,3), ynudagynu(3,3)
      DOUBLE PRECISION trynudagynu, trynuynudag, tryudagyu, tryddagyd
      DOUBLE PRECISION yudagyu(3,3), yddagyd(3,3), yedagye(3,3)
      DOUBLE PRECISION tryedagye, fuplusfd(3,3), yeyedagye(3,3)
      DOUBLE PRECISION yuyudagyu(3,3), ydyddagyd(3,3), ydyudagyu(3,3) 
      DOUBLE PRECISION yuyddagyd(3,3), ydyddagydyddag(3,3)
      DOUBLE PRECISION yuyddagydyudag(3,3), yeyedagyeyedag(3,3)
      DOUBLE PRECISION ynuyedagyeynudag(3,3), yddagydyddagyd(3,3)
      DOUBLE PRECISION yudagyuyudagyu(3,3), yudagyuyddagyd(3,3)
      DOUBLE PRECISION yuyudagyuyudag(3,3), ynuynudagynuynudag(3,3)
      DOUBLE PRECISION yddagydyudagyu(3,3), tryuyudagyuyudag
      DOUBLE PRECISION ynudagynuynudagynu(3,3), yedagyeyedagye(3,3) 
      DOUBLE PRECISION ynudagynuyedagye(3,3), trydyddagydydda 
      DOUBLE PRECISION tryeyedagyeyedag, ydyddagyuyudag(3,3)
      DOUBLE PRECISION trydyddagyuyudag, yuyudagydyddag(3,3)
      DOUBLE PRECISION tryuyudagydyddag, trydyddagydyddag

      DOUBLE PRECISION gut, gdt, get, ydb1(3,3), yub1(3,3), yeb1(3,3)
      DOUBLE PRECISION ydb2(3,3), yub2(3,3), yeb2(3,3), id(3,3)
      DOUBLE PRECISION beta1yu(9), beta2yu(9), beta1yd(9), beta2yd(9)
      DOUBLE PRECISION beta1ye(9), beta2ye(9), dydx(518)

      DOUBLE PRECISION betaC5L(3,3,3,3), betaC5R(3,3,3,3)
      DOUBLE PRECISION betaC1(3,3,3,3), betaC2(3,3,3,3) 
      DOUBLE PRECISION betaC3(3,3,3,3), betaC4(3,3,3,3)
      DOUBLE PRECISION t, a1, a2, a3, a1_unif, a2_unif, r
      DOUBLE PRECISION tol, tol1,rc, e, e1, e_next, yukgut(518)

!----------------------------------------------------------
      
      external trace,matmult,dag
      external add, mat3prod, mat4pr, rgeb1, rgeb2
!---------------------------------------------------
     
      pi = datan(1.d0) * 4.d0
      b1 = 33.d0/5.d0 
      b2 = 1.d0
      b3 = -3.d0
      MX  = 5.d0*(10.d0**19.d0)
!---------------------------------------------------
!------------------------initializing matrices to zero
      do i=1,3
       do j=1,3 

      yu(i,j) = 0.0  
      yd(i,j) = 0.0
      ye(i,j) = 0.0    
      ynu(i,j)= 0.0

         end do
        end do
!------------------------------------
C  wilson coeff.arry zero
!------------------------------------
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3

      C5L1(i,j,k,l)=0.0
      C5R1(i,j,k,l)=0.0
      C11(i,j,k,l)=0.0
      C21(i,j,k,l)=0.0
      C31(i,j,k,l)=0.0
      C41(i,j,k,l)=0.0

      gauge5L(i,j,k,l)=0.0
      gauge5R(i,j,k,l)=0.0
      yukawai5L(i,j,k,l)=0.0
      yukawaj5L(i,j,k,l)=0.0
      yukawak5L(i,j,k,l)=0.0
      yukawal5L(i,j,k,l)=0.0
      yukawai5R(i,j,k,l)=0.0
      yukawaj5R(i,j,k,l)=0.0
      yukawak5R(i,j,k,l)=0.0
      yukawal5R(i,j,k,l)=0.0
      betaC5R(i,j,k,l)=0.0
      betaC5L(i,j,k,l)=0.0
      betaC1(i,j,k,l)=0.0
      betaC2(i,j,k,l)=0.0
      betaC3(i,j,k,l)=0.0
      betaC4(i,j,k,l)=0.0 
         end do
        end do
       end do
      end do

C     Top Yukawa !!!
C     ----------------------------------

      i0 = 0 

      mssm01: do i = 1,3
      
      yu(1,i) = pyy(i0 + i)
      j = 3 + i
      yu(2,i) = pyy(i0 + j)
      k = 6 + i
      yu(3,i) = pyy(i0 + k)

      enddo mssm01

C     Bottom Yukawa !!!
C     ----------------------------------

      i0 = 9 
      mssm02: do i = 1,3
      
      yd(1,i) = pyy(i0 + i)
      j = 3 + i
      yd(2,i) = pyy(i0 + j)
      k = 6 + i
      yd(3,i) = pyy(i0 + k)
      
      enddo mssm02     
      
C     Tau Yukawa !!!
C     ----------------------------------

      i0 = 18

      mssm03: do i = 1,3
      
      ye(1,i) = pyy(i0 + i)
      j = 3 + i
      ye(2,i) = pyy(i0 + j)
      k = 6 + i
      ye(3,i) = pyy(i0 + k)

      enddo  mssm03
!-----------------------------------------------------
C------------------------------------------
C     Gauge Couplings
C-------------------------------------------
      
      alph3 = pyy(28)           !! defination, alphi = gi
      alph2 = pyy(29)
      alph1 = pyy(30)

C-----------------------------------------------

      vev1 =pyy(517)
      vev2 = pyy(518)

!------------------------------------
C  wilson coeff. intial value
!------------------------------------
      q = 31

      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3

      C5L1(i,j,k,l)=pyy(q)
      C5R1(i,j,k,l)=pyy(q+81)
      C11(i,j,k,l)=pyy(q+2*81)
      C21(i,j,k,l)=pyy(q+3*81)
      C31(i,j,k,l)=pyy(q+4*81)
      C41(i,j,k,l)=pyy(q+5*81)

       q=q+1

         end do
        end do
       end do
      end do
C-----------------------------------------------

C     Up-sector :
C     ----------
      call dag(yu,yudag)
      call matmult(yu,yudag,yuyudag)    
      call trace(yuyudag,tryuyudag)

C     Down-sector : 
C     ------------ 
      call dag(yd,yddag)
      call matmult(yd,yddag,ydyddag)  
      call trace(ydyddag,trydyddag)

C     Charged Lepton sector : 
C     ----------------------
      call dag(ye,yedag)
      call matmult(ye,yedag,yeyedag)    
      call trace(yeyedag,tryeyedag)

C     Neutrino-sector :
C     ----------
      call dag(ynu,ynudag)
      call matmult(ynu,ynudag,ynuynudag)
      call matmult(ynudag,ynu,ynudagynu)
      call trace(ynudagynu,trynudagynu)    
      call trace(ynuynudag,trynuynudag)

C     Up-sector :
C     ----------
!      call dag(yu,yudag)
      call matmult(yudag,yu,yudagyu)    
      call trace(yudagyu,tryudagyu)

C     Down-sector : 
C     ------------ 
!      call dag(yd,yddag)
      call matmult(yddag,yd,yddagyd)    
      call trace(yddagyd,tryddagyd)

C     Charged Lepton sector : 
C     ----------------------
!      call dag(ye,yedag)
      call matmult(yedag,ye,yedagye)    
      call trace(yedagye,tryedagye)

!------------------------------------
      call add(yuyudag,ydyddag,fuplusfd)
C------------------------------------------------------------------------------------------
C     Remaining matrix products(product of 3,4,5 (3x3)matrices) and their respective traces
C------------------------------------------------------------------------------------------

      call mat3prod(yu,yudag,yu,yuyudagyu)
      call mat3prod(yd,yddag,yd,ydyddagyd)
      call mat3prod(ye,yedag,ye,yeyedagye)
      call mat3prod(yd,yudag,yu,ydyudagyu)
      call mat3prod(yu,yddag,yd,yuyddagyd)

      call mat4pr(yd,yddag,yd,yddag,ydyddagydyddag)
      call mat4pr(yu,yddag,yd,yudag,yuyddagydyudag)
      call mat4pr(ye,yedag,ye,yedag,yeyedagyeyedag)
      call mat4pr(ynu,yedag,ye,ynudag,ynuyedagyeynudag)
      call mat4pr(yddag,yd,yddag,yd,yddagydyddagyd)
      call mat4pr(yudag,yu,yudag,yu,yudagyuyudagyu)
      call mat4pr(yudag,yu,yddag,yd,yudagyuyddagyd)

      call mat4pr(yu,yudag,yu,yudag,yuyudagyuyudag)
      call mat4pr(ynu,ynudag,ynu,ynudag,ynuynudagynuynudag)
      call mat4pr(yddag,yd,yudag,yu,yddagydyudagyu)
      call trace(yuyudagyuyudag,tryuyudagyuyudag)

      call mat4pr(ynudag,ynu,ynudag,ynu,ynudagynuynudagynu) 
      call mat4pr(yedag,ye,yedag,ye,yedagyeyedagye)

      call mat4pr(ynudag,ynu,yedag,ye,ynudagynuyedagye)
      call trace(ydyddagydyddag,trydyddagydyddag)
      call trace(yeyedagyeyedag,tryeyedagyeyedag)

      call mat4pr(yd,yddag,yu,yudag,ydyddagyuyudag)
      call trace(ydyddagyuyudag,trydyddagyuyudag)
      call mat4pr(yu,yudag,yd,yddag,yuyudagydyddag)
      call trace(yuyudagydyddag,tryuyudagydyddag)

!----------------------------------------------

C     -----------------------------------------------
C     Gauge Part of the Yukawa RGES (one loop terms)
C     ----------------------------------------------

      gut  = (8.d0/3.d0)*(alph3) + (1.5d0)*(alph2) + 
     $     (13.d0/30.d0)*(alph1)
      
      gdt  = (8.d0/3.d0)*(alph3) + (1.5d0)*(alph2) +
     $     (7.d0/30.d0)*(alph1)
      
      get  = (1.5d0)*(alph2) + (9.d0/10.d0)*(alph1) 
 
!---------------------------------------------------
!  C-----------------------------------------------------------setting yb1 to zero
       loopset01: do i=1,3
       loopset02:   do j=1,3
     
       ydb1(i,j)=0.d0
       yub1(i,j)=0.d0
       yeb1(i,j)=0.d0
       ydb2(i,j)=0.d0
       yub2(i,j)=0.d0
       yeb2(i,j)=0.d0

      enddo loopset02
      enddo loopset01
C-----------------------------------------------------------      
C identity matrix

      id(1,1)=1.d0          
      id(1,2)=0.d0
      id(1,3)=0.d0
      id(2,1)=0.d0
      id(2,2)=1.d0
      id(2,3)=0.d0
      id(3,1)=0.d0
      id(3,2)=0.d0
      id(3,3)=1.d0   

C--------------------------------------------------------------------
C     one loop- yukawas-matrix terms : y_b1 is the output matrix
C--------------------------------------------------------------------
      call rgeb1(yuyudag,ynuynudag,yudagyu,yddagyd,yu,yub1)
      call rgeb1(ydyddag,yeyedag,yedagye,ynudagynu,ye,yeb1)      
      call rgeb1(ydyddag,yeyedag,yddagyd,yudagyu,yd,ydb1)
      
C--------------------------------------------------------------------------
C     two loops- yukawas - matrix terms : y_b2 is the output matrix
C---------------------------------------------------------------------------

      call rgeb2(ydyddagydyddag,yuyddagydyudag,yeyedagyeyedag,
     $     ynuyedagyeynudag,yudagyu,yuyudag,ynuynudag,yddagyd,ydyddag,
     $    yeyedag,yddagydyddagyd,yudagyuyudagyu,yudagyuyddagyd,yd,ydb2)
      
      call rgeb2(yuyudagyuyudag,yuyddagydyudag,ynuynudagynuynudag,
     $     ynuyedagyeynudag,yddagyd,ydyddag,yeyedag,yudagyu,yuyudag,
     $     ynuynudag,yudagyuyudagyu,yddagydyddagyd,yddagydyudagyu,
     $     yu,yub2)
     
      call rgeb2(ydyddagydyddag,yuyddagydyudag,yeyedagyeyedag,
     $    ynuyedagyeynudag,ynudagynu,yuyudag,ynuynudag,yedagye,yddagyd,
     $     yeyedag,yedagyeyedagye,ynudagynuynudagynu,ynudagynuyedagye,
     $     ye,yeb2)

C     --------------------------------------------------------------
C     BETA1 $ BETA2 FOR YUKAWAS
C     --------------------------------------------------------------
       do i=1,9
       beta1yu(i) = 0.0
       beta2yu(i) = 0.0
       beta1yd(i) = 0.0
       beta2yd(i) = 0.0
       beta1ye(i) = 0.0
       beta2ye(i) = 0.0
       enddo

        c=1
        
        rgei: do i=1,3
        rgej: do j=1,3
        
C     yu terms--------------------------------------------------------
        beta1yu(c) = + (gut*yu(i,j)) - (yub1(i,j)) !<--------------- one loop
        
        beta2yu(c) = yub2(i,j)  !<--------------- two loop begins
     $       -(8.d0*(alph3)+((2.d0/5.d0)*(alph1)))*tryuyudag*yu(i,j)-   
     $       ((3.d0*(alph2))+((1.d0/5.d0)*(alph1)))*yuyudagyu(i,j) 
     $       -((1.d0/5.d0)*(alph1)*yuyddagyd(i,j)) +
     $       ((8.d0/9.d0)*(alph3**2)*yu(i,j)) -
     $       (4.d0*(alph3)*(alph2)*yu(i,j))- 
     $       (68.d0/45.d0)*(alph3)*(alph1)*yu(i,j) - 
     $       (15.d0/4.d0)*(alph2**2)*yu(i,j)-
     $       (1.d0/2.d0)*(alph2)*(alph1)*yu(i,j)-
     $       (2743.d0/900.d0)*(alph1**2)*yu(i,j)
       
C     yd terms--------------------------------------------------------

      beta1yd(c) = (gdt*yd(i,j)) - (ydb1(i,j))                !<--------------- one loop 
      
      beta2yd(c) = ydb2(i,j)                                    !<--------------- two loops begin  
     $ -((8.d0*(alph3))-((1.d0/5.d0)*(alph1)))*trydyddag*yd(i,j) -  
     $ ((3.d0*(alph2))+((2.d0/5.d0)*(alph1)))*ydyddagyd(i,j)  
     $  -((2.d0/5.d0)*(alph1)*ydyudagyu(i,j)) -                 
     $  (3.d0/5.d0)*(alph1)*tryeyedag*yd(i,j)+
     $  ((8.d0/9.d0)*(alph3**2)*yd(i,j)) -
     $  (4.d0*(alph3)*(alph2)*yd(i,j))- 
     $  (4.d0/9.d0)*(alph3)*(alph1)*yd(i,j) - 
     $  (15.d0/4.d0)*(alph2**2)*yd(i,j)-
     $  (0.5d0)*(alph2)*(alph1)*yd(i,j)-
     $  (287.d0/180.d0)*(alph1**2)*yd(i,j)

C ye terms--------------------------------------------------------------

      beta1ye(c) =  get*ye(i,j) - (yeb1(i,j))     !<--------------- one loop 

      beta2ye(c) = yeb2(i,j)                 !<--------------- two loops begin
     $  -((8.d0*(alph3))-((1.d0/5.d0)*(alph1)))*trydyddag*ye(i,j)-  
     $  (3.d0/5.d0)*(alph1)*tryeyedag*ye(i,j)-
     $  (3.d0)*(alph2)*yeyedagye(i,j) - 
     $  ((15.d0/4.d0)*(alph2**2)*ye(i,j)-
     $  (9.d0/10.d0)*((alph2*alph1))*ye(i,j)-
     $  (27.d0/4.d0)*(alph1**2)*ye(i,j))
      
        c = c+1

        enddo rgej
        enddo rgei
  
!=================================================================================================================
!                                       RGEs BEGIN
!==================================================================================================================

      do i = 1,518
         dydx(i) = 0.d0
      enddo

C-----------------------------------------------------
C     RGE FOR YUKAWAS
C-----------------------------------------------------
c      RGE FOR Yu
c---------------------------------------------------------
      k =1
     
      loopyuc:  do c = 1,9
     
       dydx(k) = beta1yu(c) + beta2yu(c)                 

      k = k+1
      enddo loopyuc
   
c------------------------------------------------------
c  RGE FOR Yd
c------------------------------------------------------ 
!      k = 10
      loopydc: do c=1,9

       dydx(k) = beta1yd(c) + beta2yd(c)
       
      k = k+1
      enddo loopydc           

c------------------------------------------------------
c RGE FOR Ye
c------------------------------------------------------ 
!      k =19          
      loopyec:     do c=1,9

      dydx(k) = beta1ye(c) + beta2ye(c)

      k = k+1
    
      enddo loopyec


      dydx(28)= -((alph3**(2.d0))*b3)
     $     - ((alph3**(2.d0)*
     $     (((11.d0/5.d0)*alph1) + 
     $     ((9.d0)*alph2)+
     $     ((14.d0)*alph3)-((4.d0)*tryudagyu)-((4.d0)*tryddagyd))))


      dydx(29)= -((alph2**(2.d0))*b2)
     $     - ((alph2**(2.d0)*
     $     (((9.d0/5.d0)*alph1) + 
     $     ((25.d0)*alph2) +
     $     ((24.d0)*alph3)
     $     -  ((6.d0)*tryudagyu)-((6.d0)*tryddagyd)-((2.d0)*tryedagye)
     $     - ((2.d0)*trynudagynu))))


      dydx(30)= -((alph1**(2.d0))*b1)
     $     -((alph1**(2.d0)*
     $     (((199.d0/25.d0)*alph1) + 
     $     ((27.d0/5.d0)*alph2) + ((88.d0/5.d0)*alph3)
     $     -((26.d0/5.d0)*tryudagyu)-((14.d0/5.d0)*tryddagyd)
     $     -((18.d0/5.d0)*tryedagye)-((6.d0/5.d0)*trynudagynu))))

!----------------------------------------
!-------------------------------------------------------------------
C     RGEs for vev
!-------------------------------------------------------------------
     
      dydx(517) = ((-3.d0/8.d0) * ((alph1/5.d0) + alph2) + 
     $     1.5d0 * trydyddag + 0.5d0 * tryeyedag) * vev1   +
     $     (( - (3.d0/8.d0) *                                         !<------ 2 loop starts here
     $     (3.d0 * trydyddagydyddag + 3.d0 * trydyddagyuyudag + 
     $     tryeyedagyeyedag)) + ((0.4d0 * alph1 + 4.5d0 * alph2 + 
     $     20.d0 * alph3) * trydyddag/2.d0) + ((1.8d0 * alph1 + 
     $     1.5d0 * alph2) * tryeyedag/2.d0) + ((5967.d0/3200.d0) * 
     $     alph1 * alph1) + ((1485.d0/128.d0) * alph2 * alph2) +
     $     ((27.d0/80.d0) * alph1 * alph2)) * vev1


      dydx(518) = ((-3.d0/8.d0) * ((alph1/5.d0) + alph2) + 
     $     1.5d0*tryuyudag) * vev2 +
     $     (( - (3.d0/8.d0) * (3.d0 * tryuyudagyuyudag + 3.d0 *       !<------- 2 loop starts here
     $     tryuyudagydyddag)) + ((1.9d0 * alph1 + 4.5d0 * alph2 + 
     $     20.d0 * alph3) * tryuyudag/2.d0) + ((5967.d0/3200.d0) * 
     $     alph1 * alph1) + ((1485.d0/128.d0) * alph2 * alph2) +
     $     ((27.d0/80.d0) * alph1 * alph2)) * vev2


!----------------------------------------------------------------------
!------Yukawa part
           
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3


      loopi: do ip=1,3

      yukawai5L(i,j,k,l)=yukawai5L(i,j,k,l)+C5L1(ip,j,k,l)
     $                   *fuplusfd(ip,i)
      yukawaj5L(i,j,k,l)=yukawaj5L(i,j,k,l)+C5L1(i,ip,k,l)
     $                   *yedagye(ip,j)
      yukawak5L(i,j,k,l)=yukawak5L(i,j,k,l)+C5L1(i,j,ip,l)
     $                   *fuplusfd(ip,k)
      yukawal5L(i,j,k,l)=yukawal5L(i,j,k,l)+C5L1(i,j,k,ip)
     $                  *fuplusfd(ip,l)


      yukawai5R(i,j,k,l)=yukawai5R(i,j,k,l)+C5R1(ip,j,k,l)
     $                   *(2.d0*yudagyu(ip,i))
      yukawaj5R(i,j,k,l)=yukawaj5R(i,j,k,l)+C5R1(i,ip,k,l)
     $                  *(2.d0*yddagyd(ip,j))
      yukawak5R(i,j,k,l)=yukawak5R(i,j,k,l)+C5R1(i,j,ip,l)
     $                   *(2.d0*yeyedag(ip,k))
      yukawal5R(i,j,k,l)=yukawal5R(i,j,k,l)+C5R1(i,j,k,ip)
     $                   *(2.d0*yudagyu(ip,l))

       enddo loopi

!---------gauge part of rge
!----------------------------------------------------------------------
!--here fududag=16*pi**2*yuyudag(sir notatio) and g^2=16*pi**2*alph
!------Yukawa part
!-----------------------------------------------
C g1**2=16*pi**2*alph1
!------------------------------------------------
!Note....Note.......here i am using 2/5 and 12/5 with g1 square.....need to be check??????? 
      gauge5L(i,j,k,l) = -((2.d0/5.d0)*alph1+6.d0*alph2+8.d0*alph3)
     $                  *16.d0*pi**2.d0*C5L1(i,j,k,l)

      gauge5R(i,j,k,l) = -((12.d0/5.d0)*alph1+8.d0*alph3)
     $                  *16.d0*pi**2.d0*C5R1(i,j,k,l) 


         end do
        end do
       end do
      end do
!-----------------------------------------------
C fu*fu(RGE)=16*pi**2*fu*fu(program)
!------------------------------------------------
!------------------------RGEs of the wilson Coeff. from GUT to Msusy  
             
      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3

      betaC5L(i,j,k,l)= (1.d0/(16.d0*pi**2.d0))*(gauge5L(i,j,k,l)+
     $                  (16.d0*pi**2.d0)*(yukawai5L(i,j,k,l)+
     $                  yukawaj5L(i,j,k,l)+yukawak5L(i,j,k,l)
     $                  +yukawal5L(i,j,k,l)))
         
      betaC5R(i,j,k,l)= (1.d0/(16.d0*pi**2.d0))*(gauge5R(i,j,k,l)+
     $                  (16.d0*pi**2.d0)*(yukawai5R(i,j,k,l)+
     $                  yukawaj5R(i,j,k,l)+yukawak5R(i,j,k,l)
     $                  +yukawal5R(i,j,k,l)))    
!------------------------RGEs of the wilson Coeff. from Msusy to Mz
!------------------------here yukawa coupling contribution is neglected
!------------------------alph=alpha/4*pi

      betaC1(i,j,k,l)=((alph1)*(-11.d0/10.d0)+
     $   (alph2)*(-9.d0/2.d0)+(alph3)*(-4.d0))
     $   *C11(i,j,k,l)

      betaC2(i,j,k,l)=((alph1)*(-23.d0/10.d0)+
     $   (alph2)*(-9.d0/2.d0)+(alph3)*(-4.d0))
     $   *C21(i,j,k,l)

      betaC3(i,j,k,l)=(((alph1)*(-1.d0/5.d0)+
     $   (alph2)*(-3.d0)+(alph3)*(-4.d0))
     $   *C31(i,j,k,l))+((alph2)*(-4.d0)*(C31(j,i,k,l)
     $   +C31(k,j,i,l)+C31(i,k,j,l)))

      betaC4(i,j,k,l) = ((-6.d0/5.d0)*alph1 + (-4.d0)*alph3)
     $                *C41(i,j,k,l) + (-4.d0)*alph1*C41(k,j,i,l)   

         end do
        end do
       end do
      end do
      
      q = 31

      do i=1,3
       do j=1,3
        do k=1,3
         do l=1,3

       dydx(q) = -0.5*betaC5L(i,j,k,l)
       dydx(q+81) = -0.5*betaC5R(i,j,k,l)
       dydx(q+2*81) = -0.5*betaC1(i,j,k,l)
       dydx(q+3*81) = -0.5*betaC2(i,j,k,l)
       dydx(q+4*81) = -0.5*betaC3(i,j,k,l)
       dydx(q+5*81) =-0.5* betaC4(i,j,k,l)
       q = q+1

         end do
        end do
       end do
      end do

!-----------------------------------------------------------------------------
C     RGES END
C-----------------------------------------------------------------------------      
      e = (MX/(dexp(t/2.d0)))
      
      a1 = 1.d0/(alph1*4.d0*pi)
      a2 = 1.d0/(alph2*4.d0*pi)
      a3 = 1.d0/(alph3*4.d0*pi)
     
      a1_unif = a1
      a2_unif = a2
            
      r = (1.d0 - (MIN(a1_unif,a2_unif)/
     $     MAX(a1_unif,a2_unif)))      
!------------------------------------------------------

      tol  = 0.001d0
      tol1 = 10.d0
      rc = 0
            e1 = 0.d0
            e1 = e             
            e_next = e 
!      fuscale = 0
      
      
      RETURN
 
      end subroutine protonrge
!===============================================================================================
C                                   Proton RGE ENDS
!==============================================================================================

!     subroutine dagger
!1------------------------------------------------------
      SUBROUTINE dag2(A, B)
      IMPLICIT NONE
      INTEGER  i,j              !Matrix Dimensions
      DOUBLE PRECISION A(2,2)             !Matrix A 
      DOUBLE PRECISION B(2,2)             !Matrix B
      loopdi: DO i = 1, 2
      loopdj: DO j = 1, 2
      B(i,j) = 0.d0
      B(i,j) = A(j,i)
      
      END DO loopdj
      ENDDO loopdi
      RETURN
      
      END SUBROUTINE dag2

!-------------------------------------------------------
!     subroutine dagger
!1------------------------------------------------------
      SUBROUTINE dag4(A, B)
      IMPLICIT NONE
      INTEGER  i,j              !Matrix Dimensions
      DOUBLE PRECISION A(4,4)             !Matrix A 
      DOUBLE PRECISION B(4,4)             !Matrix B
      loopdi: DO i = 1, 4
      loopdj: DO j = 1, 4
      B(i,j) = 0.d0
      B(i,j) = A(j,i)
      
      END DO loopdj
      ENDDO loopdi
      RETURN
      
      END SUBROUTINE dag4

!-------------------------------------------------------
!     subroutine dagger
!1------------------------------------------------------
      SUBROUTINE dag6(A, B)
      IMPLICIT NONE
      INTEGER  i,j              !Matrix Dimensions
      DOUBLE PRECISION A(6,6)             !Matrix A 
      DOUBLE PRECISION B(6,6)             !Matrix B
      loopdi: DO i = 1, 6
      loopdj: DO j = 1, 6
      B(i,j) = 0.d0
      B(i,j) = A(j,i)
      
      END DO loopdj
      ENDDO loopdi
      RETURN
      
      END SUBROUTINE dag6









 


