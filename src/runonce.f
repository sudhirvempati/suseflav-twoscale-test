****f*SuSeFLAV/runonce.f 
* Program : SuSeFLAV: Supersymmetry seesaw and flavor violation calculator
* Version : 1.0
* Website : http://cts.iisc.ernet.in/Suseflav/main.html
* Authors : Debtosh Choudhury  debtosh@cts.iisc.ernet.in 
**           Raghuveer Garani   rgarani@cts.iisc.ernet.in
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

      Program runonce

      IMPLICIT NONE 

      INTEGER i,j,rhn, errge,qmix

      DOUBLE PRECISION m0,m12,m10,m20, sgnmu
      double precision tanbeta,a0,ue3,Mtpole,MR1,MR2,MR3
      double precision Ynui(3,3),Mbpole,Mtaupole,MZpole
      DOUBLE PRECISION alph,Gf,alphas,alphemi
      character (len=255) :: cwd

      integer lopt, prnstat
      character*1 quit,orthcheck
      character*3 case
      CHARACTER*4 model
      character*100 mr3er
      real tstart(2),tstop(2)            
      double precision total,spectol ! For receiving total time

C--------------------------------------------------------------------------- 
      DOUBLE PRECISION Mg1,Mg2,Mg3
      DOUBLE PRECISION mq11,mq12,mq13,mq21,mq22,mq23,mq31,mq32,mq33 
      DOUBLE PRECISION mu11,mu12,mu13,mu21,mu22,mu23,mu31,mu32,mu33 
      DOUBLE PRECISION md11,md12,md13,md21,md22,md23,md31,md32,md33 
      DOUBLE PRECISION ml11,ml12,ml13,ml21,ml22,ml23,ml31,ml32,ml33 
      DOUBLE PRECISION me11,me12,me13,me21,me22,me23,me31,me32,me33 
      DOUBLE PRECISION mnu11,mnu12,mnu13,mnu21,mnu22,mnu23,mnu31,mnu32,
     $     mnu33 
      DOUBLE PRECISION a0u11,a0u12,a0u13,a0u21,a0u22,a0u23,a0u31,a0u32,
     $     a0u33

      DOUBLE PRECISION a0d11,a0d12,a0d13,a0d21,a0d22,a0d23,a0d31,a0d32,
     $     a0d33

      DOUBLE PRECISION a0e11,a0e12,a0e13,a0e21,a0e22,a0e23,a0e31,a0e32,
     $     a0e33

      DOUBLE PRECISION a0nu11,a0nu12,a0nu13,a0nu21,a0nu22,a0nu23,
     $     a0nu31,a0nu32,a0nu33

C---GMSB

      DOUBLE PRECISION gmsbsusyb, gmsbmess,gr,nhat
!-------R parametrization

      Double precision DM(3,3),R(3,3),Dk(3,3),nueig1,nueig2,nueig3      
      common/rpar/ DM,R,Dk
!--------------------------------------      
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/loops/ lopt,rhn
      common/quarkmix/ qmix
      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
      common/charinputs/case, model
      common/gauge/alph,Gf,alphas
      common/gmsbinputs/ gmsbsusyb, gmsbmess,gr,nhat
      common/ptolerance/spectol

!-------non universal A terms

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

      common/gauginos/Mg1,Mg2,Mg3

C----------------------------------------   

      EXTERNAL SuSeFLAV,chkorth

C     ----------------------------------------------------------------
C     INPUTS -  read from init# files
C     ----------------------------------------------------------------

!-------------------------------------------------------------
      quit = 'F'

      read(*,*) model

      wrinp1:if(model.ne.'mSUG'.and.model.ne.'NUHM'
     $     .and.model.ne.'GMSB'.and.model.ne.'NUGM'
     $     .and.model.ne.'CNUM')then

         write(*,*)'=========================='
         write(*,*) 'Wrong Input for Model'
         write(*,*)'=========================='       
         quit = 'T'
         call exit(1)
      endif wrinp1

      read(*,*) prnstat

      wrinp:if(prnstat.ne.1.and.prnstat.ne.0)then
         write(*,*)'=========================='
         write(*,*) 'Wrong Input for Printstat'
         write(*,*)'=========================='       
         quit = 'T'
         call exit(1)
      endif wrinp
      

!write(*,*) 'Input tanbeta'
      read(*,*) tanbeta

      tbeta:if(tanbeta>=1.d2.or.tanbeta<=1.d0)then
         write(*,*)'=========================='
         write(*,*) 'Bad Values of TANBETA'
         write(*,*)'=========================='       
         quit = 'T'
         call exit(1)
      endif tbeta


      BCmsugra: if(model.eq.'mSUG')then

!write(*,*) 'Input m0'
         read(*,*) m0


!write(*,*) 'Input a0'
         read(*,*)a0

!write(*,*) 'Input M_{1/2}'         
         read(*,*) M12 

!write(*,*) 'Input sign of mu (either POS or NEG)'
         read(*,*) sgnmu


      endif BCmsugra
!----------------

      BCnuhm: if(model.eq.'NUHM')then

!write(*,*) 'Input m0'
         read(*,*)m0


!write(*,*) 'Input a0'
         read(*,*)a0


!write(*,*) 'Input M_{1/2}'
         read(*,*) M12 


!write(*,*) 'Input sign of mu (either POS or NEG)'
         read(*,*) sgnmu

!write(*,*) 'Input mh10'
         read(*,*) m10


!write(*,*) 'Input mh20'
         read(*,*) m20

      endif BCnuhm
!-----------------

      BCgmsb: if(model.eq.'GMSB')then

!write(*,*) 'Input lambda'
         read(*,*) gmsbsusyb

!write(*,*) 'Input messenger scale'
         read(*,*) gmsbmess
         
!write(*,*) 'Input sign of mu (either 1 or -1)'
         read(*,*) sgnmu


!write(*,*) 'number of messengers in loop'
         read(*,*) nhat

!write(*,*) 'cgrav'
         read(*,*) gr

         endif BCgmsb

      BCnugauginos: if(model.eq.'NUGM')then

!write(*,*) 'Input m0'
         read(*,*) m0


!write(*,*) 'Input a0'
         read(*,*)a0

!write(*,*) 'Input M1_{1/2}'         
         read(*,*) Mg1 

!write(*,*) 'Input M2_{1/2}'         
         read(*,*) Mg2 

!write(*,*) 'Input M3_{1/2}'         
         read(*,*) Mg3 

!write(*,*) 'Input sign of mu (either POS or NEG)'
         read(*,*) sgnmu

      endif BCnugauginos


      BCcnum: if(model.eq.'CNUM')then

!write(*,*) 'Input mhu'
         read(*,*) m10

!write(*,*) 'Input mhd
         read(*,*) m20

!write(*,*) 'Input mq11'
         read(*,*) mq11, mq12, mq13

!write(*,*) 'Input mq22
         read(*,*) mq21,mq22, mq23

!write(*,*) 'Input mq33'
         read(*,*) mq31,mq32,mq33


!write(*,*) 'Input mu11'
         read(*,*) mu11,mu12,mu13

!write(*,*) 'Input mu22
         read(*,*) mu21,mu22,mu23

!write(*,*) 'Input mu33'
         read(*,*) mu31,mu32,mu33

!write(*,*) 'Input md11'
         read(*,*) md11,md12,md13

!write(*,*) 'Input md22
         read(*,*) md21,md22,md23

!write(*,*) 'Input md33'
         read(*,*) md31,md32,md33

!write(*,*) 'Input ml11'
         read(*,*) ml11,ml12,ml13

!write(*,*) 'Input ml22
         read(*,*) ml21,ml22,ml23

!write(*,*) 'Input ml33'
         read(*,*) ml31,ml32,ml33

!write(*,*) 'Input me11'
         read(*,*) me11,me12,me13

!write(*,*) 'Input me22
         read(*,*) me21,me22,me23

!write(*,*) 'Input me33'
         read(*,*) me31,me32,me33

!write(*,*) 'Input mnu11'
         read(*,*) mnu11,mnu12,mnu13

!write(*,*) 'Input mnu22
         read(*,*) mnu21,mnu22,mnu23

!write(*,*) 'Input mnu33'
         read(*,*) mnu31,mnu32,mnu33


!write(*,*) 'Input a011'
         read(*,*) a0u11, a0u12,a0u13

!write(*,*) 'Input a0u22'
         read(*,*) a0u21, a0u22,a0u23

!write(*,*) 'Input a0u33'
         read(*,*) a0u31, a0u32,a0u33

!write(*,*) 'Input a0d11'
         read(*,*) a0d11, a0d12,a0d13

!write(*,*) 'Input a0d22'
         read(*,*) a0d21, a0d22,a0d23

!write(*,*) 'Input a0d33'
         read(*,*) a0d31, a0d32,a0d33

!write(*,*) 'Input a0e11'
         read(*,*) a0e11, a0e12,a0e13

!write(*,*) 'Input a0e22'
         read(*,*) a0e21, a0e22,a0e23

!write(*,*) 'Input a0e33'
         read(*,*) a0e31, a0e32,a0e33

!write(*,*) 'Input a0nu11'
         read(*,*) a0nu11, a0nu12,a0nu13

!write(*,*) 'Input a0nu22'
         read(*,*) a0nu21, a0nu22,a0nu23

!write(*,*) 'Input a0nu33'
         read(*,*) a0nu31, a0nu32,a0nu33
         
         a0 = 0.d0

!write(*,*) 'Input M1_{1/2}'         
         read(*,*) Mg1 

!write(*,*) 'Input M2_{1/2}'         
         read(*,*) Mg2 

!write(*,*) 'Input M3_{1/2}'         
         read(*,*) Mg3 

!write(*,*) 'Input sign of mu (either POS or NEG)'
         read(*,*) sgnmu

      endif BCcnum

!------------------

      musgn:if(sgnmu.ne.1.d0.and.sgnmu.ne.-1.d0)then
         write(*,*)'=========================='
         write(*,*) 'Wrong Sign of MU'
         write(*,*)'=========================='       
         quit = 'T'
         call exit(1)
      endif musgn


!---------------

!write(*,*) 'Gf fermi constant'         
         read(*,*) Gf

!write(*,*) 'alpha_em '         
         read(*,*) alphemi
         alph = 1/alphemi


!write(*,*) 'alpha_strong '
         read(*,*) alphas

!write(*,*) 'MZ_pole '
         read(*,*) MZpole


!write(*,*) 'one loop or two loops'
         read(*,*) lopt

      wrinp2:if(lopt.ne.1.and.lopt.ne.2)then
         write(*,*)'=========================='
         write(*,*) 'Wrong Input for RGE Loops'
         write(*,*)'=========================='       
         quit = 'T'
         call exit(1)
      endif wrinp2


!     write(*,*) 'quark mixing'
      read(*,*) qmix
      
      wrinq2:if(qmix.ne.0.and.qmix.ne.1)then
         write(*,*)'=============================='
         write(*,*) 'Wrong Input for Quark Mixing'
         write(*,*)'=============================='       
         quit = 'T'
         call exit(1)
      endif wrinq2

!     write(*,*) 'Input 1 for rhn. 0 for no rhn'
      read(*,*) rhn

      wrinp3:if(rhn.ne.1.and.rhn.ne.0)then
         write(*,*)'======================================'
         write(*,*) 'Wrong Input for Right Handed Neutrino'
         write(*,*)'======================================'       
         quit = 'T'
         call exit(1)
      endif wrinp3

!write(*,*) 'Input case (either CKM or MNS or USD)'
         read(*,*) case

      wrinp4:if(rhn.eq.1.and.case.ne.'CKM'.and.case.ne.'MNS'
     $     .and.case.ne.'USD'.and.case.ne.'Rpr')then
         write(*,*)'============================================'
         write(*,*)'Wrong Input for Right Handed Neutrino Mixing'
         write(*,*)'============================================'       
         quit = 'T'
         call exit(1)
      endif wrinp4


      wrinp41:if(rhn.eq.1.and.qmix.eq.0)then
         write(*,*)'============================================'
         write(*,'(1x,A,A,I2,A,I2)')'Invalid Combination for Right',
     $        " Handed Neutrino Mixing", rhn, " and Quark Mixing",qmix
         write(*,*) "Setting Quark Mixing 1"
         write(*,*)'============================================'       
         qmix = 1
      endif wrinp41


! spectrum tolerance 

         read(*,*) spectol

!     write(*,*) 'Input Pole MT'
         read(*,*) Mtpole 

!     write(*,*) 'Input Pole Mb'
         read(*,*) Mbpole 

!     write(*,*) 'Input Pole MTau'
         read(*,*) Mtaupole 

         
!     write(*,*) 'MR1'
         read(*,*) MR1

!     write(*,*) 'MR2'
         read(*,*) MR2


!     write(*,*) 'MR3'
         read(*,*) MR3 

!write(*,*) 'Dirac Neutrino matrix at high energy-'
       Rpr:  if(case.ne.'Rpr')then
            ynuini:do i = 1,3
            ynuinj:   do j = 1,3
            read(*,*) Ynui(i,j) 
         enddo ynuinj
      enddo ynuini

      elseif(case.eq.'Rpr')then


! Read R matrix first

         read(*,*) R(1,1), R(1,2), R(1,3)

         read(*,*) R(2,1), R(2,2), R(2,3)

         read(*,*) R(3,1), R(3,2), R(3,3)

! Read light neutrino mass matrix

         call chkorth(R,orthcheck)

         ochk: if(orthcheck.eq.'F')then
            print*,"R should be orthogonal"
            quit = 'T'
            call exit(1)
         else
            quit = 'F'
         endif ochk


         do i = 1,3
            DM(1,i) = 0.d0
            DM(2,i) = 0.d0
            DM(3,i) = 0.d0

            Dk(1,i) = 0.d0
            Dk(2,i) = 0.d0
            Dk(3,i) = 0.d0
            
         enddo 

         DM(1,1) = dsqrt(MR1)
         DM(2,2) = dsqrt(MR2)
         DM(3,3) = dsqrt(MR3)
      

         read(*,*) nueig1 

         read(*,*) nueig2 

         read(*,*) nueig3 

         Dk(1,1) = dsqrt(nueig1*1.d-9)
         Dk(2,2) = dsqrt(nueig2*1.d-9)
         Dk(3,3) = dsqrt(nueig3*1.d-9)
      
      endif Rpr


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

!-----------------------------------------------------------------------------
      
!-----------------------------------------------------------------------------

      call getcwd(cwd)

C      open(999,FILE=trim(cwd)//'/../'//'/output.txt',status='replace')

      open(999,FILE= '../tmp/output.txt',status='replace')

      if(prnstat.eq.1)then
         open(10,FILE='suseflav.out',access='append',status='replace')
      endif


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
         
      elseif(model.eq.'NUGM')then

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
         

      elseif(model.eq.'CNUM')then
         continue
      endif
         
         
      if(quit.eq.'F')then
         
         call SuSeFLAV(prnstat,mq11,mq12,mq13,mq21,mq22,mq23,
     $        mq31,mq32,mq33,mu11,mu12,mu13,mu21,mu22,mu23,mu31,
     $        mu32,mu33,md11,md12,md13,md21,md22,md23,md31,md32,
     $        md33,ml11,ml12,ml13,ml21,ml22,ml23,ml31,ml32,ml33, 
     $        me11,me12,me13,me21,me22,me23,me31,me32,me33,mnu11,
     $        mnu12,mnu13,mnu21,mnu22,mnu23,mnu31,mnu32,mnu33,
     $        Mg1,Mg2,Mg3,errge)

C         print*,"--------------------------------------------------"
C         print*,"FLAGS AND THEIR MEANINGS"
C         print*,"AOK = Everything is fine."
C         print*,"BMUNEG = B_mu is negative at Msusy."
C         print*,"REWSB = |\mu|^2 < 0 at Msusy."
C         print*,"MUNOC = Non-convergent |\mu| at Msusy."
C         print*,"SW2NOC = Non-convergent Sin^2_thetaw at Mz."
C         print*,"TACSPEC = Spectrum is tachyonic at Msusy."
C         print*,"NPERTYUK = Non-perturbative yukawa."
C         print*,"TACSPECMZ = Spectrum is tachyonic at Mz."
C         print*,"FSNC = Final spectrum non-convergent."
C         print*,"TACMh = Lightest CP-even neutral higgs tachyonic."
C         print*,"TACMH = Heaviest CP-even neutral higgs tachyonic."
C         print*,"TACMA = CP-odd neutral higgs tachyonic."
C         print*,"TACMHCH = Charged higgs tachyonic."
C         print*,"VARUNDER = Stepsize is zero while integrating the", 
C     $        " RGEs."
C         print*,"TACSUP = SUP sector tachyonic."
C         print*,"TACSDN = SDOWN sector tachyonic."
C         print*,"TACSLP = SLEPTON sector tachyonic."
C         print*,"TACSNU = SNEUTRINO sector tachyonic."
C         print*,"LEPH = Lightest higgs mass below LEP limit."
C         print*,"LEPC = Lightest chargino mass < 103.5 GeV."
C         print*,"LSPSTAU = Lightest stau is LSP."
C         print*,"--------------------------------------------------"
         
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
