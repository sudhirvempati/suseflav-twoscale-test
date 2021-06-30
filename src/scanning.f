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

      INTEGER pcount,rhn,input,qmix
      DOUBLE PRECISION m0,m12,m10,m20,sgnmu

      double precision tanbeta,a0,ue3,Mtpole,MR1,MR2,MR3
      double precision Ynui(3,3),Mbpole,Mtaupole,tanbetanew
      double precision mh0sq,gm2,Bbsg,bmeg,btmug,bteg,brmu3e,
     $     brtau3mu,brtau3e,MZpole,mueconver

      integer lopt, prnstat

      character*3 case
      CHARACTER*4 model
      character*100 flags
      integer Niter,ido,idum1,idum2,idum3,idum4,idum5
      
      real knran

      double precision m0avg,m0diff,m0max,m0min
      double precision M12avg,M12min,M12max,M12diff
      DOUBLE PRECISION m10avg,m10max,m10min, m10diff
      DOUBLE PRECISION m20avg,m20max,m20min,m20diff
      DOUBLE PRECISION a0avg, a0diff, a0max,a0min
      double precision tbmin, tbmax, tbavg, tbdiff
      double precision alph,Gf,alphas, alphemi,spectol
      
      double precision brbd,brbs,M3t,SUegg(6),SDegg(2),
     $     SLegg(2),mA0sq,mHu0sq,mHpmsq,mur,Neg1,Ceg1,delm,
     $     mh0sqcor,mhu0sqcor,mA0sqcor,mHpmsqcor,AURG33,yuRG33

      integer tacsupz,tacsdnz,tacslpz,tacsnuz,tachiggsz
      
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

!     write(*,*) 'MZ_pole '
      read(input,*) MZpole

C     ----------------------------------------------------------------
C     scanning INPUTS : 
C     ----------------------------------------------------------------
      
      Write(*,*)"SUSY breaking Inputs"

      Write(*,*)"Number of Iterations" 
      Write(*,*)"(make sure it is even)"
      READ( input,* )Niter
      print*,"Nitr", Niter

      Write(*,*)"Maximum value of m0"
      READ( input,*) m0max
      print*,"m0max", m0max

      Write(*,*)"Minimum value of m0"
      READ( input, * ) m0min
      print*,"Nitr", m0min

      Write(*,*)"Maximum value of M_(1/2)"
      READ( input, *) M12max 
      print*,"Nitr", m12max

      Write(*,*)"Minimum value of M_(1/2)"
      READ( input, * ) M12min 
      print*,"Nitr", m12min

      Write(*,*)"Maximum value of A0"
      READ( input, *) a0max 
      print*,"Nitr", a0max

      Write(*,*)"Minimum value of A0"
      READ( input, *) a0min 
      print*,"Nitr", a0min


      Write(*,*)"Minimum value of m10"
      READ( input, *) m10max 
      Write(*,*)"Minimum value of m10"
      READ( input, * ) m10min 
      Write(*,*)"Minimum value of m20"
      READ( input, *) m20max 
      Write(*,*)"Minimum value of m20"
      READ( input, *) m20min 


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

      wrmod:if(model.ne.'mSUG'.and.model.ne.'NUHM')then
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
      
C     inputs done
!-----------------------------

      m0avg     = (m0max + m0min)/2.d0
      m0diff    = (m0max - m0min)/2.d0
      M12avg    = (M12max + M12min)/2.d0
      M12diff   = (M12max - M12min)/2.d0
      a0avg     = (a0max + a0min)/2.d0
      a0diff    = (a0max - a0min)/2.d0
      m10avg     = (m10max + m10min)/2.d0
      m10diff    = (m10max - m10min)/2.d0
      m20avg     = (m20max + m20min)/2.d0
      m20diff    = (m20max - m20min)/2.d0
      

 
!      tbmin = 3.d0
!      tbmax = 45.d0

!      tbavg = (tbmin + tbmax)/2.d0
!      tbdiff = (-tbmin + tbmax)/2.d0

       
C     ==================================================================
C     WRITE STATEMENTS 
C     ==================================================================
      
      Write(*,*)"------------------------------------------------------"
      print*,"Program running for Scanning Parameter Space valid for", 
     $     " SUSY"
      Write(*,*)"------------------------------------------------------"
      Write(*,*)"INPUTS :"
      print '(1x,A,ES12.5)', "tanbeta =",tanbetanew 
      print '(1x,2(A,ES12.5,A,ES12.5,A))', "m0 -> {",m0min,",", m0max,
     $     " };"," M_1/2 -> {",M12min,",", M12max," }"
      
      Write(*,*)"**************************************************"
!      Write(*,*)"SAMPLE POINT"
      
C     ==================================================================
      
      
C     set the seeds for the Random Number Generator. 
      
      idum1 = -5 
      idum2 = -15 
      idum3 = -25 
      idum4 = -35 
      idum5 = -45
   
      pcount = 0

C     ----------------------------
C     Main Loop : Doing the SCAN. 
C     BEGIN SCAN !!!!	
C     ----------------------------

      open(999,FILE='scan.out',ACCESS='APPEND',
     $     status='replace')


      loop0: do Ido = 1, Niter
           
!      Nmid = Niter/2  

      pcount = pcount + 1

C     --------------------------------------------------------------
C     Choose a point in the SUSY breaking parameter space 
C     --------------------------------------------------------------
      
      m0  = m0avg  + (2.d0*knran(idum1) - 1.d0)*m0diff 
      M12 = M12avg + (2.d0*knran(idum2) - 1.d0)*M12diff
      a0  = a0avg  + (2.d0*knran(idum3) - 1.d0)*a0diff
      M10 = M10avg + (2.d0*knran(idum4) - 1.d0)*M10diff
      M20 = M20avg  + (2.d0*knran(idum5) - 1.d0)*M20diff

!      m10 = m0
!      m20 = m0
!      a0 = 0.d0

C     Uncomment the following if you want to scan in tanbeta parameter
C     space also. 
C     tanbetanew  = tbavg  + (2.d0*knran(idum4) - 1.d0)*tbdiff

 
!      print*,"m0,m12,A0 in scanning =", m0, m12, A0

 11   FORMAT(I2,14x,'#',A)
 12   FORMAT(1pE13.5,4x,'#',A)
 13   FORMAT(1pE13.5,4x,'#',A)
 14   FORMAT(1pE13.5,4x,'#',A)
 15   FORMAT(1pE13.5,4x,'#',A)
 16   FORMAT(A3,13x,'#',A)
 17   FORMAT(1pE13.5,4x,'#',A)
 18   FORMAT(A4,12x,'#',A)

      if(model.eq.'mSUG')then
      OPEN(692, FILE = 'scan-msug.in',status='replace')
      else if(model.eq.'NUHM')then
      OPEN(692, FILE = 'scan-nuhm.in',status='replace')
      endif

      write(692,18) model , ' MODEL name'
      write(692,11) prnstat , ' prinstat: 1 = yes, 0 = no'
      write(692,12) tanbetanew, 'tanbeta'
      write(692,12) m0, ' m0'
      write(692,12) a0, ' a0'
      write(692,12) m12, ' M1/2'
      write(692,12) sgnmu, ' sgn(mu)'    

      if(model.eq.'NUHM')then
         write(692,12) m10, ' m10 = m_hu'
         write(692,12) m20, ' m20 = m_hd' 
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

      WRITE(999,82) tanbeta,m0, m12,a0,sgnmu,mh0sq,
     $     gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags


      else if(model.eq.'NUHM')then


         CALL SYSTEM('./suseflav <scan-nuhm.in')
      
         open(998,FILE='../tmp/output.txt')
         
         
         read(998,84) tanbeta,m0,m12,a0,m10,m20,mh0sq,
     $        gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags

 84      format(15(1pE15.5,2x),A)

         WRITE(999,84) tanbeta,m0,m12,a0,m10,m20,mh0sq,
     $        gm2,Bbsg,bmeg,btmug,bteg,brmu3e,brtau3mu,brtau3e,flags

      endif

      close(692)
      close(998)

      cycle


      enddo loop0

      close(999)

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
