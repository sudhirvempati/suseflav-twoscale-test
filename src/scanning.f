****f*SuSeFLAV/scanning.f 
* Program : SuSeFLAV: Supersymmetry seesaw and flavor violation calculator
* Version : 1.0
* Website : http://cts.iisc.ernet.in/Suseflav/main.html
* Authors : Debtosh Choudhury  debtosh@cts.iisc.ernet.in 
*           Raghuveer Garani   rgarani@cts.iisc.ernet.in
*           Sudhir Vempati     vempati@cts.iisc.ernet.in
*  NAME
*    Program scanspace
*  SYNOPSIS
*    Main program to scan parameter space for SuSeFLAV. 
*  FUNCTION
*     Collects points with cuts using random number generator
*  INPUTS
*    Inputs read from sinput_scan.in
*    tanbeta  -- ratio of vevs, input at msusy scale 
*   
*  RESULT
*     Set of input parameters  (m0,m12,a0,m10,m20) are written in the following .txt files, 
*     each of these files check for constraints such as lep limits, rewsb
*
*    flaghiggs.txt    -  points with light higgs < 114.1 GeV (Lep limit) 
*    chargino.txt     -  points with chargino < 103.5 GeV 

*  EXAMPLE
*    ----
*  NOTES
*  
*  BUGS
*    ---
*  SEE ALSO
*    -----
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\      
      PROGRAM scanspace

      IMPLICIT NONE

      INTEGER pcount,rhn,input,qmix,k,k1
      DOUBLE PRECISION m0,m12,m10,m20,sgnmu

      double precision tanbeta,a0,ue3,Mtpole,MR1,MR2,MR3
      double precision Ynui(3,3),Mbpole,Mtaupole,tanbetanew
      double precision mh0sq,gm2,Bbsg,bmeg,btmug,bteg,brmu3e,
     $     brtau3mu,brtau3e,MZpole,mueconver

      integer lopt, prnstat,m,p

      character*3 case
      CHARACTER*4 model
      character*100 flags,flags7

      
      real knran

      double precision m0avg,m0diff,m0max,m0min
      double precision M12avg,M12min,M12max,M12diff
      DOUBLE PRECISION m10avg,m10max,m10min, m10diff
      DOUBLE PRECISION m20avg,m20max,m20min,m20diff
      DOUBLE PRECISION a0avg, a0diff, a0max,a0min
      double precision tbmin, tbmax, tbavg, tbdiff
      double precision alph,Gf,alphas, alphemi,spectol
      double precision  mg1,a0u11,C1,msusy,gutscale
      
      double precision brbd,brbs,SUegg(6),SDegg(2),
     $     SLegg(2),mA0sq,mHu0sq,mHpmsq,Neg1,Ceg1,delm,
     $     mh0sqcor,mhu0sqcor,mA0sqcor,mHpmsqcor,AURG33,yuRG33

      integer tacsupz,tacsdnz,tacslpz,tacsnuz,tachiggsz
      double precision Cmax,Cmin,Cavg,Cdiff,C,Mpl

      integer N,i,j,Niter,Nmax,idum1,idum2,idum3,idum4,idum9,idum10
      integer idum5,idum6,idum7,idum8,idum11,m1,b1,b2,b3,b4,Ido

!      double precision,allocatable,dimension(:,:,:,:):: dacaii,dacaij
!      double precision,allocatable,dimension(:,:,:,:):: aacaij
!      double precision,allocatable,dimension(:,:,:):: dapii,dapij,alpij
!      double precision,allocatable,dimension(:,:,:):: calpha,m0ii,m0ij
!      double precision,allocatable,dimension(:,:,:):: m0iisq,m0ijsq
!      double precision,allocatable,dimension(:,:,:):: saca3
!     double precision,allocatable,dimension(:,:,:):: A0ij,saca1,saca2
!      double precision,allocatable,dimension(:,:):: sap3,sap1,sap2
!      double precision,allocatable,dimension(:,:):: M121,M122,M123

      DOUBLE PRECISION mten0avg,mten0max,mten0min,mten0diff,mfive0avg
      DOUBLE PRECISION mfive0max,mfive0min,mfive0diff,n10avg,n10max
      DOUBLE PRECISION n10min,n10diff,n5avg,n5diff,n5max,n5min
      DOUBLE PRECISION del10avg,del10max,del10min,del10diff
      DOUBLE PRECISION del5avg,del5diff,del5max,del5min,mten0
      DOUBLE PRECISION mfive0,n10,n5,del10,del5

      double precision mtilde,g0,daltaij

!      INTEGER,allocatable,dimension(:,:) :: idum2,idum3
!      double precision,allocatable,dimension(:,:) :: x,z,dalpha,salpha
!      double precision,allocatable,dimension(:,:,:) :: w,dalpha
!      double precision,allocatable,dimension(:) :: m0gut
!      double precision,allocatable,dimension(:):: y
!      double precision,allocatable,dimension(:,:):: z,x


!-------------------------------------------------------------
      double precision mu11,mu22,mu33,mu12,mu23,mu13,m01
      double precision md11,md22,md33,md12,md23,md13,m1211
      double precision mq11,mq22,mq33,mq12,mq23,mq13,a011
      double precision ml11,ml22,ml33,ml12,ml23,ml13
      double precision me11,me22,me33,me12,me23,me13
      double precision au11,au22,au33,au12,au23,au13
      double precision ad11,ad22,ad33,ad12,ad23,ad13
      double precision ae11,ae22,ae33,ae12,ae23,ae13
      double precision M1tz,M2tz,M3tz,mur,newtan,mA
      double precision a7,mqtilsq,M3t,deltaLL12,deltaRR12,deltaLR12,
     $            deltaRL12,deltaELL12,deltaERR12,deltaELR12,
     $            deltaERL12
      double precision ONM1,ONM2,ONM3,ONM4,thetat,gluino
      double precision mh,mH0,mA0,mHpm,snutau,snumu,snuelec
      double precision chi01,chi02,chi03,chi04,chipm1,chipm2
      double precision stop1,stop2,scharm1,scharm2,sup1,sup2
      double precision sbot1,sbot2,sstrange1,sstrange2,sdown1,sdown2
      double precision stau1,stau2,smu1,smu2,selec1,selec2
      double precision halfppil(3),halfpetal(3),halfpkl(3),halfpipnu(3)
      double precision halfpkpnu(3),deltaLL13,deltaRR13,deltaERR13

      DOUBLE PRECISION CRLddunu123,CRLddunu121,CRLddunuG123,
     $  CRLddunuG121,CRLuddnu123,CRLuddnu121,CRLuddnuG123,CRLuddnuG121,
     $  CRLuddnuchipm123,CRLuddnuchipm121,CLLuddnu123,CLLuddnu121,
     $  CLLuddnuG123,CLLuddnuG121,CLLuddnuchipm123,CLLuddnuchipm121,
     $ CRLddunuchi0123,CRLddunuchi0121,CRLuddnuchi0123,CRLuddnuchi0121,
     $  CLLuddnuchi0123,CLLuddnuchi0121,CRLuddnu213,CRLuddnu211,
     $  CRLuddnuG213,CRLuddnuG211,CRLuddnuchipm213,CRLuddnuchipm211,
     $  CLLuddnu213,CLLuddnu211,CLLuddnuG213,CLLuddnuG211,
     $  CLLuddnuchipm213,CLLuddnuchipm211,CRLuddnuchi0211,
     $  CRLuddnuchi0213,CLLuddnuchi0213,CLLuddnuchi0211,tautotal,etotal

      double precision ALppil1,ARppil1,CRLudul11,CRLudulG11,
     $  CRLudulchipm11,CRLudulchi011,CLLudul11,CLLudulG11,
     $  CLLudulchipm11,CLLudulchi011,CRRudul11,CRRudulG11,
     $  CRRudulchipm11,CRRudulchi011,CLRudul11,CLRudulG11,
     $  CLRudulchipm11,CLRudulchi011

      DOUBLE PRECISION halfppil6d(3),halfpetal6d(3),halfpkl6d(3),
     $   halfpipnu6d(3),halfpkpnu6d(3),
     $  halfpipnu6dfull,halfpkpnu6dfull,halfppiepure6d,halfppimupure6d,
     $  halfppipnuepure6d,halfknuepure6d,halfpknumupure6d,
     $  halfpknupure6d,halfpk0epure6d,halfpk0mupure6d,halfpipnufull,
     $  halfpkpnufull
C------           
      real tstart(2),tstop(2)            
      double precision total                  ! For receiving total time

      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/loops/ lopt, rhn
      common/counter_tag/pcount
      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
      common/charinputs/case, model
      common/gauge/alph,Gf,alphas

      EXTERNAL SuSeFLAV
      
!------------------------------------------------------------------------------
C     -----------------------------------------------------------------
C     INPUTS -  read from input files
C     ----------------------------------------------------------------
      input = 666
      OPEN( input, FILE='sinputs_scan.in') !,FORM='FORMATTED')

      READ( input, *) prnstat
      
      write(*,*) 'one loop or two loops'
      READ( input, * ) lopt
      
      write(*,*) 'Input tanbeta'
      READ( input, * ) tanbetanew
      
      write(*,*) 'Model- NUHm or mSUG or CNUM'
      READ( input, *)  model
      
      write(*,*) 'Input sign of mu (either 1 or -1)'
      READ( input, *) sgnmu

C     Pole Masses

      write(*,*) 'Input Pole MT'
      READ( input, *) Mtpole

      write(*,*) 'Input Pole Mb'
      READ( input, *) mbpole
            
      write(*,*) 'Input Pole MTau'     
      READ( input, *) mtaupole
            
      write(*,*) 'Input gf'
      READ( input, *) gf
      
      write(*,*) 'Input alphaem^(-1)'
      READ( input, *) alphemi
      alph = 1/alphemi

      write(*,*) 'Input alphas'      
      READ( input, *) alphas

      write(*,*) 'MZ_pole '
      read(input,*) MZpole

      write(*,*) "Number of Iterations"    
      Write(*,*)"(make sure it is even)"
      READ( input,* )Niter
      print*,"Niter", Niter

C     ----------------------------------------------------------------
C     scanning INPUTS : 
C     ----------------------------------------------------------------
 
      Write(*,*)"SUSY breaking Inputs"

      write(*,*) 'Input mten0max'      
      READ( input, *) mten0max
      print*,"mten10max=",mten0max

      write(*,*) 'Input mten0min'      
      READ( input, *) mten0min
      print*,"mten0min=",mten0min

      write(*,*) 'Input mfive0max'      
      READ( input, *) mfive0max
      print*,"mfive0max=",mfive0max

      write(*,*) 'Input mfive0min'      
      READ( input, *) mfive0min
      print*,"mfive0min=",mfive0min

      write(*,*) 'Input M12max'      
      READ( input, *) M12max
      print*,"M12max=",M12max

      write(*,*) 'Input M12min'      
      READ( input, *) M12min
      print*,"M12min=",M12min

      write(*,*) 'Input m10max'      
      READ( input, *) m10max
      print*,"m10max=",m10max

      write(*,*) 'Input m10min'      
      READ( input, *) m10min
      print*,"m10min=",m10min

      write(*,*) 'Input m20max'      
      READ( input, *) m20max
      print*,"m20max=",m20max

      write(*,*) 'Input m20min'      
      READ( input, *) m20min
      print*,"m20min=",m20min

      write(*,*) 'Input a0max'      
      READ( input, *) a0max
      print*,"a0max=",a0max

      write(*,*) 'Input a0min'      
      READ( input, *) a0min
      print*,"a0min=",a0min

      write(*,*) 'Input n10max'      
      READ( input, *) n10max
      print*,"n10max=",n10max

      write(*,*) 'Input n10min'      
      READ( input, *) n10min
      print*,"n10min=",n10min

      write(*,*) 'Input n5max'      
      READ( input, *) n5max
      print*,"n5max=",n5max

      write(*,*) 'Input n5min'      
      READ( input, *) n5min
      print*,"n5min=",n5min

      write(*,*) 'Input del10max'      
      READ( input, *) del10max
      print*,"del10max=",del10max

      write(*,*) 'Input del10min'      
      READ( input, *) del10min
      print*,"del10min=",del10min

      write(*,*) 'Input del5max'      
      READ( input, *) del5max
      print*,"del5max=",del5max

      write(*,*) 'Input del5min'      
      READ( input, *) del5min
      print*,"del5min=",del5min


C     ----------------------------------------------------------------
C     scanning INPUTS : 
C     ----------------------------------------------------------------
!     write(*,*) 'quark mixing'

      read(input,*) qmix     
C     Right handed neutrinos

      write(*,*) 'Input 1 for rhn. 0 for no rhn'
      READ( input, *) rhn           
      
C     case is either CKM or MNS
      write(*,*) 'Input case (either CKM or MNS or USD)'
      READ( input, * ) case

      write(*,*) 'Spectrum tolerance'
      READ( input, * ) spectol

      write(*,*) 'MR1'
      READ( input, * ) MR1

      write(*,*) 'MR2'
      READ( input, * ) MR2

      write(*,*) 'MR3'
      READ( input, * ) MR3
      
      write(*,*) 'Dirac Neutrino matrix at high energy-'
      READ( input, * ) Ynui(1,1)
      READ( input, * ) Ynui(1,2)
      READ( input, * ) Ynui(1,3)
      READ( input, * ) Ynui(2,1)
      READ( input, * ) Ynui(2,2)
      READ( input, * ) Ynui(2,3)
      READ( input, * ) Ynui(3,1)
      READ( input, * ) Ynui(3,2)
      READ( input, * ) Ynui(3,3)

      CLOSE(input)
!-----------------------------------------------------------------------------
      wrmod:if(model.ne.'mSUG'.and.model.ne.'NUHM'.and.
     $         model.ne.'CNUM')then
!     $     .and.model.ne.'CNUM')then !.and.model.ne.'NUGM'
!     $     .and.model.ne.'CNUM')then
         
         write(*,*)'=========================='
         write(*,*) 'Wrong Input for model'
         write(*,*)'=========================='       
!     quit = 'T'
         call exit(1)
      endif wrmod

      wrinq2:if(qmix.ne.0.and.qmix.ne.1)then
         write(*,*)'=============================='
         write(*,*) 'Wrong Input for Quark Mixing'
         write(*,*)'=============================='       
!         quit = 'T'
         call exit(1)
      endif wrinq2

      wrinp41:if(rhn.eq.1.and.qmix.eq.0)then
         write(*,*)'============================================'
         write(*,'(1x,A,A,I2,A,I2)')'Invalid Combination for Right',
     $        " Handed Neutrino Mixing", rhn, " and Quark Mixing",qmix
         write(*,*) "Setting Quark Mixing 1"
         write(*,*)'============================================'       
         qmix = 1
      endif wrinp41

      wrinp:if(prnstat.ne.1.d0.and.prnstat.ne.0.d0)then
         write(*,*)'=========================='
         write(*,*) 'Wrong Input for Printstat'
         write(*,*)'=========================='       
!     quit = 'T'
         call exit(1)
      endif wrinp

      tbeta:if(tanbetanew>=1.d2.or.tanbetanew<=1.d0)then
         write(*,*)'=========================='
         write(*,*) 'Bad Values of TANBETA'
         write(*,*)'=========================='       
!     quit = 'T'
         call exit(1)
      endif tbeta

      musgn:if(sgnmu.ne.1.d0.and.sgnmu.ne.-1.d0)then
         write(*,*)'=========================='
         write(*,*) 'Wrong Sign of MU'
         write(*,*)'=========================='       
!     quit = 'T'
         call exit(1)
      endif musgn
      
      wrinp1:if(rhn.ne.1.and.rhn.ne.0)then
         write(*,*)'======================================='
         write(*,*) 'Wrong Input for Right Handed Neutrino '
         write(*,*)'======================================='       
!     quit = 'T'
         call exit(1)
      endif wrinp1
      
      wrinp3:if(rhn.eq.1.and.case.ne.'CKM'.and.case.ne.'MNS'
     $     .and.case.ne.'USD')then
         write(*,*)'============================================='
         write(*,*) 'Wrong Input for Right Handed Neutrino Mixing'
         write(*,*)'============================================='       
!     quit = 'T'
         call exit(1)
      endif wrinp3
           
      wrinp4:if(lopt.ne.1.and.lopt.ne.2)then
         write(*,*)'=========================='
         write(*,*) 'Wrong Input for RGE Loops'
         write(*,*)'=========================='       
!     quit = 'T'
         call exit(1)
      endif wrinp4     
!-----------------------------------------------------------------------------    
      open(170,FILE='output.out')
      open(171,FILE='protondecay.out')
      open(172,FILE='protondecayhiggs.out')
      open(173,FILE='gluino_chargino.out')
      open(174,FILE='pgoespieplus.out')
      open(179,FILE='M1M2M3.out')
      open(181,FILE='mixing.out')
      open(182,FILE='higgs.out')
      open(183,FILE='neutrilino.out')
      open(184,FILE='stop.out')
      open(185,FILE='sbot.out')
      open(186,FILE='slepton.out') 
      open(187,FILE='protondecayfull.out')
      open(188,FILE='protondecay5d.out')
      open(189,FILE='protondecay5d6d.out')
      open(190,FILE='protondecay5d6dwithouthiggs.out')
C     inputs done

      mten0avg     = (mten0max + mten0min)/2.d0
     
      mten0diff    = (mten0max - mten0min)/2.d0
c********
c********
      mfive0avg     = (mfive0max + mfive0min)/2.d0
      mfive0diff    = (mfive0max - mfive0min)/2.d0
c********
c********
      M12avg    = (M12max + M12min)/2.d0
      M12diff   = (M12max - M12min)/2.d0
      
      m10avg     = (m10max + m10min)/2.d0
      m10diff    = (m10max - m10min)/2.d0
      m20avg     = (m20max + m20min)/2.d0
      m20diff    = (m20max - m20min)/2.d0
      
      a0avg  = (a0max + a0min)/2.d0
      a0diff  = (a0max - a0min)/2.d0

      n10avg = (n10max + n10min)/2.d0
      n10diff= (n10max - n10min)/2.d0

      n5avg = (n5max + n5min)/2.d0
      n5diff= (n5max - n5min)/2.d0

      del10avg = (del10max + del10min)/2.d0
      del10diff= (del10max - del10min)/2.d0

      del5avg = (del5max + del5min)/2.d0
      del5diff= (del5max - del5min)/2.d0

      
 
!      tbmin = 10.d0
!      tbmax = 30.d0

!      tbavg = (tbmin + tbmax)/2.d0
!      tbdiff = (-tbmin + tbmax)/2.d0
      
      
      pcount = 0 
C     set the seeds for the Random Number Generator. 
      
      idum1 = -13 
      idum2 = -26
      idum3 = -9
      idum4 = -29
      idum5 = -3
      idum6 = -19
      idum7=-2
      idum8=-5
      idum9=-37
      idum10=-8
      idum11=-3

      pcount = 0

C     ----------------------------
C     Main Loop : Doing the SCAN. 
C     BEGIN SCAN !!!!	
C     ----------------------------
      
      open(999,FILE='scan.out',ACCESS='APPEND',
     $     status='replace')


      loop0: do Ido = 1, Niter
       Write(*,*)"Ido" , Ido
           
!      Nmid = Niter/2  

      pcount = pcount + 1

C     --------------------------------------------------------------
C     Choose a point in the SUSY breaking parameter space 
C     --------------------------------------------------------------
      
      mten0 = mten0avg  + (2.d0*knran(idum1) - 1.d0)*mten0diff 
      mfive0=mfive0avg  + (2.d0*knran(idum2) - 1.d0)*mfive0diff !m0
     
      M12 =M12avg + (2.d0*knran(idum3) - 1.d0)*M12diff
      a0  =a0avg  + (2.d0*knran(idum4) - 1.d0)*a0diff
      M10 = M10avg + (2.d0*knran(idum5) - 1.d0)*M10diff
      M20 = M20avg  + (2.d0*knran(idum6) - 1.d0)*M20diff

      n10=n10avg+ (2.d0*knran(idum7) - 1.d0)*n10diff
      mfive0=mten0
      n5=n10
!      n10=1
!      del10=0.0
!      del5=0.0
!      n5=n5avg+ (2.d0*knran(idum8) - 1.d0)*n5diff

      del10=del10avg+ (2.d0*knran(idum9) - 1.d0)*del10diff
      del5=del5avg+ (2.d0*knran(idum10) - 1.d0)*del5diff
!      mfive0=20000.d0
!      mten0 =20000.d0
!      n5=15
!      n10=15
!!      m12=1000.d0
!      a0=4000.d0

!       M10= -4921.040d0
!       M20= 3884.94d0
!       M12=2887.690d0
!       a0= 4000.0d0
!       n10= 15.0d0
!       n5=15.0d0
!       del10= 0.8769d0
!       del5 =0.8304d0
!       mten0= 20000.0d0
!       mfive0=20000.0d0 

!       M10= -2697.59d0
!       M20= 1530.16d0
!       M12=2207.040d0
!       a0= 4000.0d0
!       n10= 15.0d0
!       n5=15.0d0
!       del10= 0.858d0
!       del5 =0.3492d0
!       mten0= 20000.0d0
!       mfive0=20000.0d0

!       M10= 1266.819d0
!       M20= 995.0689d0
!       M12=1960.02d0
!       a0= 4000.0d0
!       n10= 15.0d0
!       n5=15.0d0
!       del10= 0.8936d0
!       del5 =0.337d0
!       mten0= 20000.0d0
!       mfive0=20000.0d0
           
!       M10= 2210.559d0
!       M20= 6910.300d0
!       M12=2079.60d0
!       a0= -7710.69d0
!       n10= 2.9215d0
!       n5=7.61129d0
!       del10= 0.3432d0
!       del5 =0.7933d0
!       mten0= 7277.632d0
!       mfive0=12092.830d0

!       M10=-6531.819d0
!       M20=-1717.190d0
!       M12=2932.050d0
!       a0= -6102.086d0
!       n10= 2.92887d0
!       n5=5.0008d0
!       del10= 0.7979d0
!       del5 =0.18996d0
!       mten0= 25484.799d0
!       mfive0=28321.126d0

!       M10=-7898.2500d0
!       M20=-6268.930d0
!       M12=1980.132d0
!       a0= -7737.8451d0
!       n10= 6.2085d0
!       n5=8.8346d0
!       del10= 0.14992d0
!       del5 =0.24760d0
!       mten0= 13491.917d0
!       mfive0=2562.820d0

!       M10=-4358.390d0
!       M20=7866.930d0
!       M12=2809.74d0
!       a0= -5946.60d0
!       n10= 6.6091d0
!       n5=9.5844d0
!       del10= 0.6535d0
!       del5 =0.2402d0
!       mten0= 22354.647d0
!       mfive0=17891.91d0
!      a0max=sqrt(3.0d0)*m0
!      a0min=-sqrt(3.0d0)*m0

!       M10= -13555.6d0
!       M20= -9239.7d0
!       M12=11821.26d0
!       a0= -49017.7d0
!       n10= 6.76d0
!       n5=6.76d0
!       del10= -0.459d0
!       del5 =0.0d0
!       mten0= 70343.35d0
!       mfive0=70343.35d0

!      M10=571.764
!      M20=2556.01
!      mten0=87968.04
!      n10=1.902
!      mfive0=87968.044
!      n5=1.902355
!      a0=-72629.543
!      M12=10058.539
!      del10=0.202
!      del5=-0.0855

      M10=16353.90
      M20=22911.79
      mten0=71299.52
      n10=2.53739
      mfive0=71299.52
      n5=2.5373
      a0=-40563.55
      M12=4925.59
      del10=0.0047
      del5=-0.7127

      pcount = pcount + 1  
C     --------------------------------------------------------------
C     Choose a point in the SUSY breaking parameter space 
C     --------------------------------------------------------------    


!----------------------------------------------------------------------
!----------------------------------------------------------------------


 11   FORMAT(I2,14x,'#',A)
 12   FORMAT(1pE13.5,4x,'#',A)
 13   FORMAT(1pE13.5,4x,'#',A)
 14   FORMAT(1pE13.5,4x,'#',A)
 15   FORMAT(1pE13.5,4x,'#',A)
 16   FORMAT(A3,13x,'#',A)
 17   FORMAT(1pE13.5,4x,'#',A)
 18   FORMAT(A4,12x,'#',A)
 29   FORMAT(1pE13.5,4x,1pE13.5,4x,1pE13.5,4x,'#',A)

      if(model.eq.'mSUG')then
      OPEN(692, FILE = 'scan-msug.in',status='replace')
      else if(model.eq.'NUHM')then
      OPEN(692, FILE = 'scan-nuhm.in',status='replace')
      else if(model.eq.'CNUM')then
      OPEN(692, FILE = 'scan-cnum.in',status='replace')
      endif

      write(692,18) model , ' MODEL name'
      write(692,11) prnstat , ' prinstat: 1 = yes, 0 = no'
      write(692,12) tanbetanew, 'tanbeta'
      
      if(model.eq.'NUHM' .and. model.eq.'mSUG')then
      write(692,12) m0, ' m0'
      write(692,12) a0, ' a0'
      write(692,12) m12, ' M1/2'
      write(692,12) sgnmu, ' sgn(mu)'    
      endif

      if(model.eq.'NUHM')then
         write(692,12) m10, ' m10 = m_hu'
         write(692,12) m20, ' m20 = m_hd' 
      endif

      if(model.eq.'CNUM')then

         write(692,12) M10, ' m10 = m_hu'
         write(692,12) M20, ' m20 = m_hd' 

      write(692,29) mten0,0.0,sign(1.d0,del10)*Sqrt(abs(del10)/n10)*
     $              mten0, 'mq11,mq12,mq13'
      write(692,29) 0.0,mten0,0.0, 'mq21,mq22,mq23'
      write(692,29) sign(1.d0,del10)*Sqrt(abs(del10)/n10)*mten0,0.0, 
     $             mten0/n10,  'mq31,mq32,mq33'

      write(692,29) mten0,0.0,sign(1.d0,del10)*Sqrt(abs(del10)/n10)*
     $              mten0, 'mu11,mu12,mu13'
      write(692,29) 0.0,mten0,0.0, 'mu21,mu22,mu23'
      write(692,29) sign(1.d0,del10)*Sqrt(abs(del10)/n10)*mten0,0.0,
     $            mten0/n10,    'mu31,mu32,mu33'

      write(692,29) mfive0,0.0,0.0, 'md11,md12,md13'
      write(692,29) 0.0,mfive0,0.0, 'md21,md22,md23'
      write(692,29) 0.0,0.0,mfive0/n5, 
     $             'md31,md32,md33'

      write(692,29) mfive0,0.0,0.0, 'ml11,ml12,ml13'
      write(692,29) 0.0,mfive0,0.0, 'ml21,ml22,ml23'
      write(692,29) 0.0,0.0,mfive0/n5,
     $               'ml31,ml32,ml33'

      write(692,29) mten0,0.0,sign(1.d0,del10)*Sqrt(abs(del10)/n10)*
     $              mten0, 'me11,me12,me13'
      write(692,29) 0.0,mten0,0.0, 'me21,me22,me23'
      write(692,29) sign(1.d0,del10)*Sqrt(abs(del10)/n10)*mten0,0.0, 
     $             mten0/n10,    'me31,me32,me33'

      write(692,29) 0.0,0.0,0.0, 'mnu11,mnu12,mnu13'
      write(692,29) 0.0,0.0,0.0, 'mnu21,mnu22,mnu23'
      write(692,29) 0.0,0.0,0.0, 'mnu31,mnu32,mnu33'

      write(692,29) a0,0.0,0.0, 'au11,au12,au13'
      write(692,29) 0.0,a0,0.0, 'au21,au22,au23'
      write(692,29) 0.0,0.0,a0, 'au31,au32,au33'

      write(692,29) a0,0.0,0.0, 'ad11,ad12,ad13'
      write(692,29) 0.0,a0,0.0, 'ad21,ad22,ad23'
      write(692,29) 0.0,0.0,a0, 'ad31,ad32,ad33'

      write(692,29) a0,0.0,0.0, 'ae11,ae12,ae13'
      write(692,29) 0.0,a0,0.0, 'ae21,ae22,ae23'
      write(692,29) 0.0,0.0,a0, 'ae31,ae32,ae33'

      write(692,29) 0.0,0.0,0.0, 'anu11,anu12,anu13'
      write(692,29) 0.0,0.0,0.0, 'anu21,anu22,anu23'
      write(692,29) 0.0,0.0,0.0, 'anu31,anu32,anu33'

      write(692,12) M12, 'M1_1/2'
      write(692,12) M12, 'M2_1/2' 
      write(692,12) M12, 'M3_1/2'

      write(692,12) sgnmu, ' sgn(mu)'
      endif

      write(692,14) gf, ' G_fermi'
      write(692,15) alphemi, ' alpha_em^(-1)'    
      write(692,13) alphas, ' alpha_strong'    
      write(692,13) mzpole, ' mzpole'    
      write(692,11) 2, ' one loop or two loops'
      write(692,11) qmix, ' quark mixing'
      write(692,11) rhn, ' 1= rhn on; 0 = rhn off'
      write(692,16) case, ' case: CKM/MNS/USD'

      write(692,13) spectol, ' spectrum tolerance'
      write(692,12) MTpole , ' Mtpole'
      write(692,13) MBpole , ' Mbpole'
      write(692,13) mTaupole , ' Mtaupole'
      write(692,17) MR1,' MR1'
      write(692,17) MR2,' MR2'
      write(692,17) MR3,' MR3'
      write(692,17) Ynui(1,1), ' Ynu(1,1)'
      write(692,17) Ynui(1,2), ' Ynu(1,2)'
      write(692,17) Ynui(1,3), ' Ynu(1,3)'
      write(692,17) Ynui(2,1), ' Ynu(2,1)'
      write(692,17) Ynui(2,2), ' Ynu(2,2)'
      write(692,17) Ynui(2,3), ' Ynu(2,3)'
      write(692,17) Ynui(3,1), ' Ynu(3,1)'
      write(692,17) Ynui(3,2), ' Ynu(3,2)'
      write(692,17) Ynui(3,3), ' Ynu(3,3)'


      if(model.eq.'mSUG')then

      CALL SYSTEM('./suseflav <scan-msug.in')
      open(998,FILE='../tmp/output.txt')

      read(998,*) tanbeta,m0, m12,a0,sgnmu,mh0sq,
     $     gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags


 82   FORMAT(1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,
     $     1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,
     $     1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,1pE11.4,2x,1pE11.4,
     $     2x,1pE11.4,2x,A)

      WRITE(999,*) tanbeta,m0, m12,a0,sgnmu,mh0sq,
     $     gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags

      else if(model.eq.'NUHM')then


         CALL SYSTEM('./suseflav <scan-nuhm.in')
      
         open(998,FILE='../tmp/output.txt')
         
         
         read(998,84) tanbeta,m0,m12,a0,m10,m20,mh0sq,
     $        gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags

 84      format(15(1pE15.5,2x),A)

         WRITE(999,84) tanbeta,m0,m12,a0,m10,m20,mh0sq,
     $        gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags

      else if(model.eq.'CNUM')then


         CALL SYSTEM('./suseflav <scan-cnum.in')
      
         open(998,FILE='../tmp/output.txt')
         
         
         read(998,*) tanbeta,m01, m1211,a011,m10,m20,mh0sq,
     $        gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags

 85      format(15(1pE15.5,2x),I6,2x,G8.0)

      WRITE(999,*)M10,M20,mten0,n10,mfive0,n5,a0,M12,del10,del5,
     $            Sqrt(mh0sq),flags
       If(flags .eq. 'AOK')then
      open(141,FILE= '../tmp/check1.txt')
      open(142,FILE= '../tmp/check2.txt')
      open(143,FILE= '../tmp/check3.txt')
      open(144,FILE= '../tmp/check4.txt')
!      open(145,FILE= '../tmp/check5.txt')
!      open(146,FILE= '../tmp/check6.txt')
!      open(147,FILE= '../tmp/check7.txt')
!      open(148,FILE= '../tmp/check8.txt')
      open(149,FILE= '../tmp/check9.txt')
!      open(150,FILE= '../tmp/check10.txt')
      open(151,FILE= '../tmp/check11.txt')
      open(152,file= '../tmp/check12.txt')
      open(153,FILE= '../tmp/check13.txt')
      open(154,file= '../tmp/check14.txt')
      open(155,FILE= '../tmp/check15.txt')
      open(156,file= '../tmp/check16.txt')

!-----------------------------------------------------------
! here mq11 is scalar mass ,mA is mass of CP-odd higgs


      read(141,*)halfppil(1),halfpetal(1),halfpkl(1),
     $          halfpipnu(1),halfpkpnu(1),halfppil(2),halfpetal(2),
     $          halfpkl(2),halfpipnu(2),halfpkpnu(2),halfppil(3),
     $          halfpetal(3),halfpkl(3),halfpipnu(3),halfpkpnu(3),
     $          msusy,gutscale,deltaLL13,deltaRR13,deltaERR13
      read(142,*)CRLddunu123,CRLddunu121,CRLddunuG123,
     $  CRLddunuG121,CRLuddnu123,CRLuddnu121,CRLuddnuG123,CRLuddnuG121,
     $  CRLuddnuchipm123,CRLuddnuchipm121,CLLuddnu123,CLLuddnu121,
     $  CLLuddnuG123,CLLuddnuG121,CLLuddnuchipm123,CLLuddnuchipm121,
     $ CRLddunuchi0123,CRLddunuchi0121,CRLuddnuchi0123,CRLuddnuchi0121,
     $  CLLuddnuchi0123,CLLuddnuchi0121,CRLuddnu213,CRLuddnu211,
     $  CRLuddnuG213,CRLuddnuG211,CRLuddnuchipm213,CRLuddnuchipm211,
     $  CLLuddnu213,CLLuddnu211,CLLuddnuG213,CLLuddnuG211,
     $  CLLuddnuchipm213,CLLuddnuchipm211,CRLuddnuchi0211,
     $  CRLuddnuchi0213,CLLuddnuchi0213,CLLuddnuchi0211,tautotal,etotal

      read(143,*) ALppil1,ARppil1,CRLudul11,CRLudulG11,
     $  CRLudulchipm11,CRLudulchi011,CLLudul11,CLLudulG11,
     $  CLLudulchipm11,CLLudulchi011,CRRudul11,CRRudulG11,
     $  CRRudulchipm11,CRRudulchi011,CLRudul11,CLRudulG11,
     $  CLRudulchipm11,CLRudulchi011

      read(144,*)halfppil6d(1),halfppil6d(2),halfppil6d(3),
     $   halfpetal6d(1),halfpetal6d(2),halfpetal6d(3),halfpkl6d(1),
     $   halfpkl6d(2),halfpkl6d(3),halfpipnu6d(1),halfpipnu6d(2),
     $   halfpipnu6d(3),halfpkpnu6d(1),halfpkpnu6d(2),halfpkpnu6d(3),
     $  halfpipnu6dfull,halfpkpnu6dfull,halfppiepure6d,halfppimupure6d,
     $  halfppipnuepure6d,halfknuepure6d,halfpknumupure6d,
     $  halfpknupure6d,halfpk0epure6d,halfpk0mupure6d,halfpipnufull,
     $  halfpkpnufull,au33,ad33,ae33
!      read(143,*)mq11,mq22,mq33,mq12,mq23,mq13
!      read(144,*)ml11,ml22,ml33,ml12,ml23,ml13
!      read(145,*)me11,me22,me33,me12,me23,me13
!      read(146,*)au11,au22,au33,au12,au23,au13
!      read(147,*)ad11,ad22,ad33,ad12,ad23,ad13
!      read(148,*)ae11,ae22,ae33,ae12,ae23,ae13
      read(149,*)M1tz,M2tz,M3tz,mur,newtan,mA
!      read(150,*)a7,mqtilsq,M3t,deltaLL12,deltaRR12,deltaLR12,
!     $            deltaRL12,deltaELL12,deltaERR12,deltaELR12,
!     $            deltaERL12
      read(151,*)ONM1,ONM2,ONM3,ONM4,thetat,gluino
      read(152,*)mh,mH0,mA0,mHpm,snutau,snumu,snuelec
      read(153,*)chi01,chi02,chi03,chi04,chipm1,chipm2
      read(154,*)stop1,stop2,scharm1,scharm2,sup1,sup2
      read(155,*)sbot1,sbot2,sstrange1,sstrange2,sdown1,sdown2
      read(156,*)stau1,stau2,smu1,smu2,selec1,selec2
!      write(160,*)stop1,stop2,scharm1,scharm2,sup1,sup2,flags7
      write(171,*)halfppil(1),halfpetal(1),halfpkl(1),
     $          halfpipnu(1),halfpkpnu(1),halfppil(2),halfpetal(2),
     $          halfpkl(2),halfpipnu(2),halfpkpnu(2),halfppil(3),
     $          halfpetal(3),halfpkl(3),halfpipnu(3),halfpkpnu(3),n10,
     $          mten0,sup1/stop1,sdown1/sbot1,gluino,chi01,msusy,
     $          gutscale,del10,del5,n5,stop1,stop2,scharm1,scharm2,
     $          sup1,sup2,M1tz,M2tz,M3tz,mur,newtan,mA,sbot1,sbot2,
     $          sdown1,sdown2,stau1,stau2,selec1,selec2,a0,M12,M10,M20,
     $          del10/n10,deltaLL13,deltaRR13,deltaERR13,bteg,
     $          brtau3e,gm2
      Write(190,*)n10,mten0,del10,halfppil6d(1),halfppil6d(2),
     $        halfpkl6d(1),halfpkl6d(2),halfpkpnu6d(1),halfpkpnu6d(3),
     $        halfpipnu6dfull,halfpkpnu6dfull,gutscale,
     $  deltaLL13,deltaRR13,deltaERR13,bteg,brtau3e,gm2,au33,
     $  snutau,snumu,snuelec,stop1,stop2,scharm1,scharm2,sup1,sup2,
     $  sbot1,sbot2,sstrange1,sstrange2,sdown1,sdown2,
     $  stau1,stau2,smu1,smu2,selec1,selec2 
      if(mh .ge. 123.0d0 .and. mh .le. 128.0d0)then
      if(gluino .ge. 2200.0 .and. stop1 .ge. 1200.0)then
      write(170,*)M10,M20,mten0,n10,mfive0,n5,a0,M12,del10,del5
      write(172,*)halfppil(1),halfpetal(1),halfpkl(1),
     $          halfpipnu(1),halfpkpnu(1),halfppil(2),halfpetal(2),
     $          halfpkl(2),halfpipnu(2),halfpkpnu(2),halfppil(3),
     $          halfpetal(3),halfpkl(3),halfpipnu(3),halfpkpnu(3),n10,
     $          mten0,sup1/stop1,sdown1/sbot1,gluino,chi01,msusy,
     $          gutscale,del10,del5,n5,stop1,stop2,scharm1,scharm2,
     $          sup1,sup2,M1tz,M2tz,M3tz,mur,newtan,mA,sbot1,sbot2,
     $          sdown1,sdown2,stau1,stau2,selec1,selec2,a0,M12,M10,M20,
     $          del10/n10,deltaLL13,deltaRR13,deltaERR13,bteg,
     $          brtau3e,gm2,au33,ad33,ae33


      Write(187,*)halfppil6d(1),halfppil6d(2),halfppil6d(3),
     $   halfpetal6d(1),halfpetal6d(2),halfpetal6d(3),halfpkl6d(1),
     $   halfpkl6d(2),halfpkl6d(3),halfpipnu6d(1),halfpipnu6d(2),
     $   halfpipnu6d(3),halfpkpnu6d(1),halfpkpnu6d(2),halfpkpnu6d(3),
     $  halfpipnu6dfull,halfpkpnu6dfull,halfppiepure6d,halfppimupure6d,
     $  halfppipnuepure6d,halfknuepure6d,halfpknumupure6d,
     $  halfpknupure6d,halfpk0epure6d,halfpk0mupure6d,halfpipnufull,
     $  halfpkpnufull,n10,mten0,sup1/stop1,sdown1/sbot1,gluino,chi01,
     $  msusy, gutscale,del10,del5,n5,stop1,stop2,scharm1,scharm2,
     $  sup1,sup2,M1tz,M2tz,M3tz,mur,newtan,mA,sbot1,sbot2,
     $  sdown1,sdown2,stau1,stau2,selec1,selec2,a0,M12,M10,M20,
     $  deltaLL13,deltaRR13,deltaERR13,bteg,brtau3e,gm2,au33,ad33,ae33 
      Write(188,*)n10,mten0,del10,halfppil(1),halfppil(2),halfpkl(1),
     $           halfpkl(2),halfpkpnu(1),halfpkpnu(3),halfpipnufull,
     $           halfpkpnufull,halfppiepure6d,halfppimupure6d,
     $ halfpk0epure6d,halfpk0mupure6d,halfknuepure6d,halfpknumupure6d,
     $  halfpknupure6d,halfppipnuepure6d,
     $  deltaLL13,deltaRR13,deltaERR13,bteg,brtau3e,gm2,au33 ,
     $  snutau,snumu,snuelec,stop1,stop2,scharm1,scharm2,sup1,sup2,
     $  sbot1,sbot2,sstrange1,sstrange2,sdown1,sdown2,
     $  stau1,stau2,smu1,smu2,selec1,selec2 
      Write(189,*)n10,mten0,del10,halfppil6d(1),halfppil6d(2),
     $        halfpkl6d(1),halfpkl6d(2),halfpkpnu6d(1),halfpkpnu6d(3),
     $        halfpipnu6dfull,halfpkpnu6dfull,gutscale,
     $  deltaLL13,deltaRR13,deltaERR13,bteg,brtau3e,gm2,au33,
     $  snutau,snumu,snuelec,stop1,stop2,scharm1,scharm2,sup1,sup2,
     $  sbot1,sbot2,sstrange1,sstrange2,sdown1,sdown2,
     $  stau1,stau2,smu1,smu2,selec1,selec2 

      write(173,*)CRLddunu123,CRLddunu121,CRLddunuG123,
     $  CRLddunuG121,CRLuddnu123,CRLuddnu121,CRLuddnuG123,CRLuddnuG121,
     $  CRLuddnuchipm123,CRLuddnuchipm121,CLLuddnu123,CLLuddnu121,
     $  CLLuddnuG123,CLLuddnuG121,CLLuddnuchipm123,CLLuddnuchipm121,
     $ CRLddunuchi0123,CRLddunuchi0121,CRLuddnuchi0123,CRLuddnuchi0121,
     $  CLLuddnuchi0123,CLLuddnuchi0121,CRLuddnu213,CRLuddnu211,
     $  CRLuddnuG213,CRLuddnuG211,CRLuddnuchipm213,CRLuddnuchipm211,
     $  CLLuddnu213,CLLuddnu211,CLLuddnuG213,CLLuddnuG211,
     $  CLLuddnuchipm213,CLLuddnuchipm211,CRLuddnuchi0211,
     $  CRLuddnuchi0213,CLLuddnuchi0213,CLLuddnuchi0211,halfpkpnu(1),
     $  halfpkpnu(3),n10,n5,mten0,sup1/stop1,sdown1/sbot1,gluino,chi01,
     $  msusy,gutscale,del10,del5,n5,stop1,stop2,scharm1,scharm2,
     $  sup1,sup2,M1tz,M2tz,M3tz,mur,newtan,mA,sbot1,sbot2,
     $  sdown1,sdown2,stau1,stau2,selec1,selec2,a0,M12,M10,M20,bteg,
     $          brtau3e,gm2,tautotal,etotal

      write(174,*)ALppil1,ARppil1,CRLudul11,CRLudulG11,
     $  CRLudulchipm11,CRLudulchi011,CLLudul11,CLLudulG11,
     $  CLLudulchipm11,CLLudulchi011,CRRudul11,CRRudulG11,
     $  CRRudulchipm11,CRRudulchi011,CLRudul11,CLRudulG11,
     $  CLRudulchipm11,CLRudulchi011,halfppil(1),
     $  n10,n5,mten0,sup1/stop1,sdown1/sbot1,gluino,chi01,
     $  msusy,gutscale,del10,del5,n5,stop1,stop2,scharm1,scharm2,
     $  sup1,sup2,M1tz,M2tz,M3tz,mur,newtan,mA,sbot1,sbot2,
     $  sdown1,sdown2,stau1,stau2,selec1,selec2,a0,M12,M10,M20,bteg,
     $  brtau3e,gm2
!      write(172,*)md11,md22,md33,md12,md23,md13,100
!      write(173,*)mq11,mq22,mq33,mq12,mq23,mq13,100
!      write(174,*)ml11,ml22,ml33,ml12,ml23,ml13,100
!      write(175,*)me11,me22,me33,me12,me23,me13,100
!      write(176,*)au11,au22,au33,au12,au23,au13,100
!      write(177,*)ad11,ad22,ad33,ad12,ad23,ad13,100
!      write(178,*)ae11,ae22,ae33,ae12,ae23,ae13,100
      write(179,*)M1tz,M2tz,M3tz,mur,newtan,mA
!      write(180,*)a7,mqtilsq,M3t,deltaLL12,deltaRR12,deltaLR12,
!     $            deltaRL12,100
      write(181,*)ONM1,ONM2,ONM3,ONM4,thetat,gluino
      write(182,*)mh,mH0,mA0,mHpm,snutau,snumu,snuelec
      write(183,*)chi01,chi02,chi03,chi04,chipm1,chipm2
      write(184,*)stop1,stop2,scharm1,scharm2,sup1,sup2
      write(185,*)sbot1,sbot2,sstrange1,sstrange2,sdown1,sdown2
      write(186,*)stau1,stau2,smu1,smu2,selec1,selec2,I

!      CALL SYSTEM('./sfile susy_flavor.in')
!      write(161,*)gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e
!      write(162,*)deltaELL12,deltaERR12,deltaELR12,deltaERL12,bmeg
!      If(chi01 .le. 2000.0 .and. chi01 .lt. chipm1)then
!      CALL SYSTEM('./main slha.out')
!      write(163,*)deltaELL12,deltaERR12,deltaELR12,deltaERL12,bmeg,mur,
!     $            M3t,stop1
!      write(164,*)gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e
!      write(165,*)a7,mqtilsq,M3t,deltaLL12,deltaRR12,deltaLR12,
!     $            deltaRL12,100
      endif
      endif
      endif
      endif


      close(692)
      close(998)


      close(141)
      close(142)
      close(143)
      close(144)
!      close(145)
!      close(146)
!      close(147)
!      close(148)
      close(149)
!      close(150)
      close(151)
      close(152)
      close(153)
      close(154)
      close(155)
      close(156)

      cycle
!      endif vevcond

      enddo loop0

      close(101)
      close(999)
      close(160)
      close(161)
      close(162)
      close(163)
      close(164)
      close(165)

      close(170)
      close(171)
      close(172)
      close(173)
      close(174)
!      close(175)
!      close(176)
!      close(177)
!      close(178)
      close(179)
!      close(180)
      close(181)
      close(182)
      close(183)
      close(184)
      close(185)
      close(186)
      close(187)
      close(188)
      close(189)
      close(190)

      END PROGRAM scanspace
      
!===============================================================================

C     ====================================================================
C     -------------------------------------------------------------------
C     Portable Random Number Generator. Recommended by D. Knuth. 
C     iseed is an integer seed which can be set to any negative value. 
C     Use different seeds if you calling the same random number generator
C     at different points. 
C     -------------------------------------------------------------------
      
      FUNCTION knran(iseed)
      INTEGER iseed
      INTEGER RNBIG,RNSEED,RMZ
C     REAL RNBIG,RNSEED,RMZ
      REAL knran,RNFAC 
C      REAL RNFAC ! for modification
      PARAMETER (RNBIG=1000000000,RNSEED=161803398,RMZ=0,RNFAC=1./RNBIG)
C     PARAMETER (RNBIG=4000000.,RNSEED=1618033.,RMZ=0.,RNFAC=1./RNBIG)
      INTEGER i,iff,rnii,rinxt,rinxtp,k
      INTEGER rmj,rmk,rma(55)
C     REAL rmj,rmk,rma(55)
      SAVE iff,rinxt,rinxtp,rma
      DATA iff /0/
      if(iseed.lt.0.or.iff.eq.0)then
        iff=1
        rmj=RNSEED-iabs(iseed)
        rmj=mod(rmj,RNBIG)
        rma(55)=rmj 
        rmk=1
        do 11 i=1,54
          rnii=mod(21*i,55)
          rma(rnii)=rmk
          rmk=rmj-rmk
          if(rmk.lt.RMZ)rmk=rmk+RNBIG
          rmj=rma(rnii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            rma(i)=rma(i)-rma(1+mod(i+30,55))
            if(rma(i).lt.RMZ)rma(i)=rma(i)+RNBIG
12        continue
13      continue
        rinxt=0
        rinxtp=31
        iseed=1
      endif
      rinxt=rinxt+1
      if(rinxt.eq.56)rinxt=1
      rinxtp=rinxtp+1
      if(rinxtp.eq.56)rinxtp=1
      rmj=rma(rinxt)-rma(rinxtp)
      if(rmj.lt.RMZ)rmj=rmj+RNBIG
      rma(rinxt)=rmj
      knran=rmj*RNFAC
      return
      END

C=============================================================================
