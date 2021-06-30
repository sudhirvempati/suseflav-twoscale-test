****f*SuSeFLAV/runslha.f 
* Program : SuSeFLAV: Supersymmetry seesaw and flavor violation calculator
* Version : 1.0
* Website : http://cts.iisc.ernet.in/Suseflav/main.html
* Authors : Debtosh Choudhury  debtosh@cts.iisc.ernet.in 
*           Raghuveer Garani   rgarani@cts.iisc.ernet.in
*           Sudhir Vempati     vempati@cts.iisc.ernet.in
*  NAME
*    Program runonce
*  SYNOPSIS
*    Main program to compute spectrum at a single point for SuSeFLAV. 
*  FUNCTION
*     
*  INPUTS
*    Inputs read from sinput.in
*    tanbeta  -- ratio of vevs, input at msusy scale 
*   
*  RESULT
*    
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
C===================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     1. This will give result for a single point in the parameter space.
C     2. Checked on 10/05/2010.
C     
C     
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
      

      Program runslha

      IMPLICIT NONE 

      INTEGER rhn,errge,qmix

      DOUBLE PRECISION m0,m12,m10,m20, sgnmu,prnstat_r
      double precision tanbeta,a0,ue3,Mtpole,MR1,MR2,MR3
      double precision Ynui(3,3),Mbpole,Mtaupole,MZpole

      integer lopt, prnstat

      character (len=255) :: cwd

      character*1 quit
      character*3 case 
      CHARACTER*4 model
      INTEGER input
      real tstart(2),tstop(2)            
      double precision total              ! For receiving total time

C--------------------------------------------------------------------------- 
      DOUBLE PRECISION Mg1,Mg2,Mg3
      DOUBLE PRECISION mq11,mq12,mq13,mq21,mq22,mq23,mq31,mq32,mq33 
      DOUBLE PRECISION mu11,mu12,mu13,mu21,mu22,mu23,mu31,mu32,mu33 
      DOUBLE PRECISION md11,md12,md13,md21,md22,md23,md31,md32,md33 
      DOUBLE PRECISION ml11,ml12,ml13,ml21,ml22,ml23,ml31,ml32,ml33 
      DOUBLE PRECISION me11,me12,me13,me21,me22,me23,me31,me32,me33 
      DOUBLE PRECISION mnu11,mnu12,mnu13,mnu21,mnu22,mnu23,mnu31,mnu32,
     $     mnu33 

      integer i,num,modnum,smnum(10),minparnum(10),counter_sminp
      integer counter_minpar,counter_algopar,algonum(20)
      double precision sminp(10), minpar(10), algopar(20)
      double precision alph,gf,alphas,spectol,rhn_r
      DOUBLE PRECISION runloops, rhnmix
      character*100 mr3er,filepath
C---GMSB

      DOUBLE PRECISION gmsbsusyb, gmsbmess,gr,nhat

!--------------------------------------      
      common/sminputs/ mbpole, mtaupole, Mtpole,MZpole
      common/loops/ lopt,rhn
      common/quarkmix/ qmix
      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
      common/charinputs/case, model
      common/gauge/alph,Gf,alphas
      common/gmsbinputs/ gmsbsusyb, gmsbmess,gr,nhat
C----------------SLHA common blocks-   

      common/slha_modsel/num, modnum
      common/slha_counters/counter_sminp,counter_minpar,counter_algopar
      common/slha_num/smnum,minparnum,algonum
      common/slha_variables/sminp,minpar,algopar
      common/ptolerance/spectol

C---------------------------------------
      EXTERNAL SuSeFLAV,readslha
C------------------------------------


C     ----------------------------------------------------------------
C     INPUTS -  read from init# files
C     ----------------------------------------------------------------
      quit = 'F'

      filepath = 'slha.in'  

      input=666


      OPEN(input, FILE=filepath) !,FORM='FORMATTED')

      call readslha(input)

!------inputs read from slha --> program variables
C     Modsel

      if(num.eq.1.and.(modnum.eq.0.or.modnum.eq.1))then
         model = 'mSUG'
      endif

      if(num.eq.1.and.modnum.eq.2)then
         model = 'GMSB'
      endif

      if(num.eq.1.and.modnum.eq.4)then
         model = 'NUHM'
      endif

!      print*,"qmix in runslha = ", qmix

C sm inputs

      sminpblock: do i=1,counter_sminp     
      
      if(smnum(i).eq.1) alph = 1/sminp(i)
      if(smnum(i).eq.2) Gf = sminp(i)
      if(smnum(i).eq.3) alphas = sminp(i)
      if(smnum(i).eq.4) mzpole = sminp(i)
      if(smnum(i).eq.5) mbpole = sminp(i)
      if(smnum(i).eq.6) mtpole = sminp(i)
      if(smnum(i).eq.7) mtaupole = sminp(i)

      enddo sminpblock
C  minpar

      BCmsugra: if(model.eq.'mSUG')then
      minparblock: do i = 1, counter_minpar

      if(minparnum(i).eq.1) m0 =  minpar(i)
      if(minparnum(i).eq.2) m12 =  minpar(i)
      if(minparnum(i).eq.3) tanbeta =  minpar(i)
      if(minparnum(i).eq.4) sgnmu =  minpar(i)
      if(minparnum(i).eq.5) a0 =  minpar(i)

          enddo minparblock
      endif BCmsugra
 
      BCmgmsb: if(model.eq.'GMSB')then
      minparblock1: do i = 1, counter_minpar

      if(minparnum(i).eq.1) gmsbsusyb =  minpar(i)
      if(minparnum(i).eq.2) gmsbmess =  minpar(i)
      if(minparnum(i).eq.3) tanbeta =  minpar(i)
      if(minparnum(i).eq.4) sgnmu =  minpar(i)
      if(minparnum(i).eq.5) nhat =  minpar(i)
      if(minparnum(i).eq.6) gr =  minpar(i)

          enddo minparblock1
      endif BCmgmsb

      BCnuhm: if(model.eq.'NUHM')then
      minparblock2: do i = 1, counter_minpar

      if(minparnum(i).eq.1) m0 =  minpar(i)
      if(minparnum(i).eq.2) m12 =  minpar(i)
      if(minparnum(i).eq.3) tanbeta =  minpar(i)
      if(minparnum(i).eq.4) sgnmu =  minpar(i)
      if(minparnum(i).eq.5) a0 =  minpar(i)
      if(minparnum(i).eq.6) m10 =  minpar(i)
      if(minparnum(i).eq.7) m20 =  minpar(i)

          enddo minparblock2
      endif BCnuhm

C algorithm parameters

      sflavalgo: do i = 1, counter_algopar

      if(algonum(i).eq.1) spectol = algopar(i)

      if(algonum(i).eq.2)then
         rhn_r = algopar(i)
         if(rhn_r.eq.1.d0) rhn = 1
         if(rhn_r.eq.0.d0) rhn = 0
      endif

      if(algonum(i).eq.3)then
         prnstat_r = algopar(i)
         if(prnstat_r.eq.1.d0) prnstat = 1
         if(prnstat_r.eq.0.d0) prnstat = 0
      endif
      
      rhnmix = 0.d0
      
      if(algonum(i).eq.4) rhnmix = algopar(i)

      if(rhnmix.eq.1.d0) case = 'CKM'
      if(rhnmix.eq.2.d0) case = 'MNS'
      if(rhnmix.eq.3.d0) case = 'USD'

      runloops = 0.d0
          
      if(algonum(i).eq.5) runloops = algopar(i)
      
      if(runloops.eq.1.d0) lopt = 1
      if(runloops.eq.2.d0) lopt = 2
      
!      ue3 = 0.d0
      
!      if(algonum(i).eq.6) ue3 = algopar(i)
      if(algonum(i).eq.7) MR1 = algopar(i)
      if(algonum(i).eq.8) MR2 = algopar(i)
      if(algonum(i).eq.9) MR3 = algopar(i)

      if(algonum(i).eq.10) Ynui(1,1)= algopar(i)
      if(algonum(i).eq.11) Ynui(1,2)= algopar(i)
      if(algonum(i).eq.12) Ynui(1,3)= algopar(i)
      if(algonum(i).eq.13) Ynui(2,1)= algopar(i)
      if(algonum(i).eq.14) Ynui(2,2)= algopar(i)
      if(algonum(i).eq.15) Ynui(2,3)= algopar(i)
      if(algonum(i).eq.16) Ynui(3,1)= algopar(i)
      if(algonum(i).eq.17) Ynui(3,2)= algopar(i)
      if(algonum(i).eq.18) Ynui(3,3)= algopar(i)

      enddo sflavalgo
!-----------------

      wrinq2:if(qmix.ne.0.and.qmix.ne.1)then
         write(*,*)'=============================='
         write(*,*) 'Wrong Input for Quark Mixing'
         write(*,*)'=============================='       
         quit = 'T'
         call exit(1)
      endif wrinq2

      wrmod:if(model.ne.'mSUG'.and.model.ne.'NUHM'
     $     .and.model.ne.'GMSB')then !.and.model.ne.'NUGM'
!     $     .and.model.ne.'CNUM')then
         
         write(*,*)'=========================='
         write(*,*) 'Wrong Input for model'
         write(*,*)'=========================='       
         quit = 'T'
         call exit(1)
      endif wrmod

      wrinp:if(prnstat_r.ne.1.d0.and.prnstat_r.ne.0.d0)then
         write(*,*)'=========================='
         write(*,*) 'Wrong Input for Printstat'
         write(*,*)'=========================='       
         quit = 'T'
         call exit(1)
      endif wrinp


      tbeta:if(tanbeta>=1.d2.or.tanbeta<=1.d0)then
         write(*,*)'=========================='
         write(*,*) 'Bad Values of TANBETA'
         write(*,*)'=========================='       
         quit = 'T'
         call exit(1)
      endif tbeta

      musgn:if(sgnmu.ne.1.d0.and.sgnmu.ne.-1.d0)then
         write(*,*)'=========================='
         write(*,*) 'Wrong Sign of MU'
         write(*,*)'=========================='       
         quit = 'T'
         call exit(1)
      endif musgn
      

      wrinp1:if(rhn_r.ne.1.d0.and.rhn_r.ne.0.d0)then
         write(*,*)'======================================='
         write(*,*) 'Wrong Input for Right Handed Neutrino '
         write(*,*)'======================================='       
         quit = 'T'
         call exit(1)
      endif wrinp1
      
      wrinp3:if(rhn_r.eq.1.d0.and.case.ne.'CKM'.and.case.ne.'MNS'
     $     .and.case.ne.'USD')then
         write(*,*)'============================================='
         write(*,*) 'Wrong Input for Right Handed Neutrino Mixing'
         write(*,*)'============================================='       
         quit = 'T'
         call exit(1)
      endif wrinp3
      
      
      wrinp4:if(lopt.ne.1.and.lopt.ne.2)then
         write(*,*)'=========================='
         write(*,*) 'Wrong Input for RGE Loops'
         write(*,*)'=========================='       
         quit = 'T'
         call exit(1)
      endif wrinp4

      wrinp41:if(rhn_r.eq.1.d0.and.qmix.eq.0)then
         write(*,*)'============================================'
         write(*,'(1x,A,A,I2,A,I2)')'Invalid Combination for Right',
     $        " Handed Neutrino Mixing", rhn, " and Quark Mixing",qmix
         write(*,*) "Setting Quark Mixing 1"
         write(*,*)'============================================'       
         qmix = 1
      endif wrinp41


!-------------------------

      if(model.eq.'GMSB'.and.rhn.eq.1)then

         if(gmsbsusyb/gmsbmess.gt.1.d0)then

            mr3er = " lambda/messenger should be less than 1"
            print*,"========================="
            print*, mr3er
            print*," Program aborted"
            print*,"=========================="
            quit = 'T'
            call exit(1)
         endif

         if(MR3.ge.gmsbmess)then
            mr3er = " Error: Messenger scale lesser than seesaw scale"
            print*,"========================="
            print*, mr3er
            print*," Program aborted"
            print*,"=========================="
            quit = 'T'   
            call exit(1)
         ENDIF

      endif
      
      CALL getcwd(cwd)
C      open(999,FILE=trim(cwd)//'/../'//'/output.txt',status='replace')
      open(999,FILE='../tmp/output.txt',status='replace')

      if(prnstat.eq.1)then
         open(10,FILE='suseflav.out',access='append',status='replace')
      endif

!-----------------------------------------------------------------------------
      
!-----------------------------------------------------------------------------



      if(model.eq.'mSUG')then

         m10 = m0
         m20 = m0

         mq11 = m0
         mq12 = 0.d0
         mq13 = 0.d0
         mq21 = 0.d0
         mq22 = m0
         mq23 = 0.d0
         mq31 = 0.d0
         mq32 = 0.d0
         mq33 = m0

         mu11 = m0
         mu12 = 0.d0
         mu13 = 0.d0
         mu21 = 0.d0
         mu22 = m0
         mu23 = 0.d0
         mu31 = 0.d0
         mu32 = 0.d0
         mu33 = m0

         md11 = m0
         md12 = 0.d0
         md13 = 0.d0
         md21 = 0.d0
         md22 = m0
         md23 = 0.d0
         md31 = 0.d0
         md32 = 0.d0
         md33 = m0

         ml11 = m0
         ml12 = 0.d0
         ml13 = 0.d0
         ml21 = 0.d0
         ml22 = m0
         ml23 = 0.d0
         ml31 = 0.d0
         ml32 = 0.d0
         ml33 = m0

         me11 = m0
         me12 = 0.d0
         me13 = 0.d0
         me21 = 0.d0
         me22 = m0
         me23 = 0.d0
         me31 = 0.d0
         me32 = 0.d0
         me33 = m0

         mnu11 = m0
         mnu12 = 0.d0
         mnu13 = 0.d0
         mnu21 = 0.d0
         mnu22 = m0
         mnu23 = 0.d0
         mnu31 = 0.d0
         mnu32 = 0.d0
         mnu33 = m0
         
         Mg1 = m12
         Mg2 = m12
         Mg3 = m12

      elseif(model.eq.'NUHM')then


         mq11 = m0
         mq12 = 0.d0
         mq13 = 0.d0
         mq21 = 0.d0
         mq22 = m0
         mq23 = 0.d0
         mq31 = 0.d0
         mq32 = 0.d0
         mq33 = m0

         mu11 = m0
         mu12 = 0.d0
         mu13 = 0.d0
         mu21 = 0.d0
         mu22 = m0
         mu23 = 0.d0
         mu31 = 0.d0
         mu32 = 0.d0
         mu33 = m0

         md11 = m0
         md12 = 0.d0
         md13 = 0.d0
         md21 = 0.d0
         md22 = m0
         md23 = 0.d0
         md31 = 0.d0
         md32 = 0.d0
         md33 = m0

         ml11 = m0
         ml12 = 0.d0
         ml13 = 0.d0
         ml21 = 0.d0
         ml22 = m0
         ml23 = 0.d0
         ml31 = 0.d0
         ml32 = 0.d0
         ml33 = m0

         me11 = m0
         me12 = 0.d0
         me13 = 0.d0
         me21 = 0.d0
         me22 = m0
         me23 = 0.d0
         me31 = 0.d0
         me32 = 0.d0
         me33 = m0

         mnu11 = m0
         mnu12 = 0.d0
         mnu13 = 0.d0
         mnu21 = 0.d0
         mnu22 = m0
         mnu23 = 0.d0
         mnu31 = 0.d0
         mnu32 = 0.d0
         mnu33 = m0
         
         Mg1 = m12
         Mg2 = m12
         Mg3 = m12
         
         
c$$$  elseif(model.eq.'NUGM')then
c$$$  
c$$$  m10 = m0
c$$$  m20 = m0
c$$$  
c$$$  mq11 = m0
c$$$  mq12 = 0.d0
c$$$  mq13 = 0.d0
c$$$  mq21 = 0.d0
c$$$  mq22 = m0
c$$$  mq31 = 0.d0
c$$$  mq32 = 0.d0
c$$$  mq33 = m0
c$$$  
c$$$  mu11 = m0
c$$$  mu12 = m0
c$$$  mu13 = m0
c$$$  mu21 = m0
c$$$  mu22 = m0
c$$$  mu23 = m0
c$$$  mu31 = m0
c$$$  mu32 = m0
c$$$  mu33 = m0
c$$$  
c$$$  md11 = m0
c$$$  md12 = m0
c$$$  md13 = m0
c$$$  md21 = m0
c$$$  md22 = m0
c$$$  md23 = m0
c$$$  md31 = m0
c$$$  md32 = m0
c$$$  md33 = m0
c$$$  
c$$$  ml11 = m0
c$$$  ml12 = m0
c$$$  ml13 = m0
c$$$  ml21 = m0
c$$$  ml22 = m0
c$$$  ml23 = m0
c$$$  ml31 = m0
c$$$  ml32 = m0
c$$$  ml33 = m0
c$$$  
c$$$  me11 = m0
c$$$  me12 = m0
c$$$  me13 = m0
c$$$  me21 = m0
c$$$  me22 = m0
c$$$  me23 = m0
c$$$  me31 = m0
c$$$  me32 = m0
c$$$  me33 = m0
c$$$  
c$$$  mnu11 = m0
c$$$  mnu12 = m0
c$$$  mnu13 = m0
c$$$  mnu21 = m0
c$$$  mnu22 = m0
c$$$  mnu23 = m0
c$$$  mnu31 = m0
c$$$  mnu32 = m0
c$$$  mnu33 = m0
c$$$  
c$$$  
c$$$  elseif(model.eq.'CNUM')then
c$$$  continue

      endif

      
      if(quit.eq.'F')then

         call SuSeFLAV(prnstat,mq11,mq12,mq13,mq21,mq22,mq23,
     $        mq31,mq32,mq33,mu11,mu12,mu13,mu21,mu22,mu23,mu31,mu32,
     $        mu33,md11,md12,md13,md21,md22,md23,md31,md32,md33,
     $        ml11,ml12,ml13,ml21,ml22,ml23,ml31,ml32,ml33, 
     $        me11,me12,me13,me21,me22,me23,me31,me32,me33 ,
     $        mnu11,mnu12,mnu13,mnu21,mnu22,mnu23,mnu31,mnu32,mnu33,
     $        Mg1,Mg2,Mg3,errge)


c$$$         print*,"--------------------------------------------------"
c$$$         print*,"FLAGS AND THEIR MEANINGS"
c$$$         print*,"AOK = Everything is fine."
c$$$         print*,"BMUNEG = B_mu is negative at Msusy."
c$$$         print*,"REWSB = |\mu|^2 < 0 at Msusy."
c$$$         print*,"MUNOC = Non-convergent |\mu| at Msusy."
c$$$         print*,"SW2NOC = Non-convergent Sin^2_thetaw at Mz."
c$$$         print*,"NPERTYUK = Non-perturbative yukawa."
c$$$         print*,"TACSPEC = Spectrum is tachyonic at Msusy."
c$$$         print*,"TACSPECMZ = Spectrum is tachyonic at Mz."
c$$$         print*,"FSNC = Final spectrum non-convergent."
c$$$         print*,"TACMh = Lightest CP-even neutral higgs tachyonic."
c$$$         print*,"TACMH = Heaviest CP-even neutral higgs tachyonic."
c$$$         print*,"TACMA = CP-odd neutral higgs tachyonic."
c$$$         print*,"TACMHCH = Charged higgs tachyonic."
c$$$         print*,"VARUNDER = Stepsize is zero while integrating the", 
c$$$     $        " RGEs."
c$$$         print*,"TACSUP = SUP sector tachyonic."
c$$$         print*,"TACSDN = SDOWN sector tachyonic."
c$$$         print*,"TACSLP = SLEPTON sector tachyonic."
c$$$         print*,"TACSNU = SNEUTRINO sector tachyonic."
c$$$         print*,"LEPH = Lightest higgs mass below LEP limit."
c$$$         print*,"LEPC = Lightest chargino mass < 103.5 GeV."
c$$$         print*,"LSPSTAU = Lightest stau is LSP."
c$$$         print*,"--------------------------------------------------"

      if(prnstat.eq.1)then
      write(10,*)"--------------------------------------------------"
      write(10,*)"FLAGS AND THEIR MEANINGS"
      write(10,*)"AOK = Everything is fine."
      write(10,*)"BMUNEG = B_mu is negative at Msusy."
      write(10,*)"REWSB = |\mu|^2 < 0 at Msusy."
      write(10,*)"MUNOC = Non-convergent |\mu| at Msusy."
      write(10,*)"SW2NOC = Non-convergent Sin^2_thetaw at Mz."
      write(10,*)"NPERTYUK = Non-perturbative yukawa."
      write(10,*)"TACSPEC = Spectrum is tachyonic at Msusy."
      write(10,*)"TACSPECMZ = Spectrum is tachyonic at Mz."
      write(10,*)"FSNC = Final spectrum non-convergent."
      write(10,*)"TACMh = Lightest CP-even neutral higgs tachyonic."
      write(10,*)"TACMH = Heaviest CP-even neutral higgs tachyonic."
      write(10,*)"TACMA = CP-odd neutral higgs tachyonic."
      write(10,*)"TACMHCH = Charged higgs tachyonic."
      write(10,*)"VARUNDER = Stepsize is zero while integrating the",
     $     " RGEs."
      write(10,*)"TACSUP = SUP sector tachyonic."
      write(10,*)"TACSDN = SDOWN sector tachyonic."
      write(10,*)"TACSLP = SLEPTON sector tachyonic."
      write(10,*)"TACSNU = SNEUTRINO sector tachyonic."
      write(10,*)"LEPH = Lightest higgs mass below LEP limit."
      write(10,*)"LEPC = Lightest chargino mass < 103.5 GeV."
      write(10,*)"LSPSTAU = Lightest stau is LSP."
      write(10,*)"--------------------------------------------------"
      endif

      elseif(quit.eq.'T')then
         print*,'WRONG INPUTS'
         call exit(1)
      endif


      close(999)

      if(prnstat.eq.1)then
         close(10)
      endif
      
      END 
      
!===============================================================================
