****f*SuSeFLAV/runrges.f/MSSMRUN
*  NAME
*    SUBROUTINE MSSMRUN
*  SYNOPSIS
*    Integrates mssm RGEs from MZ to M_{GUT} and M_{GUT} to msusy.
*  FUNCTION
*    Integrates mssm RGEs from MZ to M_{GUT} and M_{GUT} to msusy for
*    given inputs at ew scale.
*    Runs smrge and mssm rges at with complete 2-loops and 3x3 flavor
*    structure.
*  INPUTS
*     vevin      - vev at MZ
*     yuin       - (3x3) up type yukawa matrix
*     ydin       - (3x3) down type yukawa matrix
*     yein       - (3x3) matrix yukawa for leptons
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
*     murge    -  RGE output: \mu at M_{susy} scale
*     bmurge   -  RGE output: b_\mu at M_{susy} scale
*     newtbeta -  Ratio of vev at msusy from rge running.
*     flags    -  flags problem with rge running, if any.
*
*  EXAMPLE
*      subroutine MSSMRUN(vevin,yuin,ydin,yein,alph1in,alph2in,
*     $     alph3in,mur,bmur,murge,bmurge,prnstat,check,newtbeta,
*     $     msusy,runum,itcount,flags)
*
*  NOTES
*     Common blocks and external routines used:
*      COMMON/rgeyy/yy_d
*      common/unif/ e1, yukgut
*      common/sminputs/ mbpole, mtaupole, Mtpole
*      common/loops/ lopt,rhn
*      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
*      common/charinputs/case, model
*      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
*      common/rgeinput_high/ MX, M1X,M2X,M3X,mh10,mh20,
*     $     mQ0,mU0,mD0,mE0,mL0,mNU0
*      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
*     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut
*      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
*     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG
*      common/gauge/alph,Gf,alphas
*
*      External routines:
*       external RK4ROUTINE,QMSRK4,smrge,mssmrge,matmult,dag
*       external smrgemt
*
*  BUGS
*    ---
*  SEE ALSO
*    RK4ROUTINE,smrgemt,math.f
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C

C     MATRIX TERMS IN RGE IS COMPUTED USING USER DEFINED SUBROUTINES RGEB1 AND RGEB2
C
C     AUTHOR:-Debtosh Chowdhury,Raghuveer Garani & Sudhir K Vempati         MODIFIED :- 01st December 2009 7:30PM
!----------------------------------------------------------------------------------------------------



      subroutine MSSMRUN(vevin,yuin,ydin,yein,alph1in,alph2in,
     $     alph3in,mur,bmur,murge,bmurge,prnstat,check,newtbeta,
     $     msusy,runum,itcount,flags)

      implicit none
      integer i0,i,n0,j0,j,l,nok,nbad,k
      integer check,lopt,rhn,itcount

      integer LWORK,info,fuscale
      parameter(LWORK = 35)
      double precision WORK(LWORK)

      character*1 quitint
      CHARACTER*4 model
      CHARACTER*100 flags
      DOUBLE PRECISION e1
      double precision yy(126),yy_d(126),MX
      data yy/126 * 0.d0/, yy_d/126 * 0.d0/
      double precision mQ0(3,3),mU0(3,3),mD0(3,3),mL0(3,3),mE0(3,3)
      double precision mh10,mh20,mSQRG(3,3),mSURG(3,3)
      double precision mSDRG(3,3),AURG(3,3),ADRG(3,3),AERG(3,3)
      double precision mSLRG(3,3),mSERG(3,3),mh1mz,mh2mz,murge
      double precision tZ,yuin(3,3),ydin(3,3),yein(3,3),bmurge
      double precision beta,tanbeta
      double precision vev1,vev2
      double precision x1,x2,h1,hmin,eps,M10,M20,m0,m12
      double precision YERG(3,3),YURG(3,3),YDRG(3,3),tq0,msusy
      double precision M1X,M2X,M3X
      double precision ANURG(3,3),mNU0(3,3),mSNURG(3,3)
      double precision mnu(3,3),mnuev(3)
      double precision alph1,alph2,alph3,mnu_mz(3,3),perm(3,3)
      DOUBLE PRECISION ydgut(3,3),yegut(3,3),Ynui(3,3),ynugut(3,3)
      DOUBLE PRECISION M1tz,M2tz,M3tz,newtbeta
      double precision alph1in,alph2in,alph3in

!     Neutrino stuff
      character*3 case
      double precision Mpl,tpl,t2,t1

      double precision DM_atm, DM_sol, mnu_scale ! neutrinos Delta M
      character*3 HIE
      double precision thl12, thl23, thl13
      DOUBLE PRECISION ynuD1(3,3),ynuD2(3,3),ynuD3(3,3)

      double precision MR3,MR2,MR1,t3,mur,bmur
      double precision yugut(3,3)

!     Effective theory stuff
      DOUBLE PRECISION yy_sm(31)

      INTEGER prnstat,runum
      DOUBLE PRECISION a0

      double precision  ue3

      double precision vevin,mbpole, mtaupole, Mtpole
      double precision mtscale

      double precision VCKM(3,3),MZpole,MZ,pi

!------------------------------------------

      double precision sgnmu,yukgut(126), newgut
      double precision alph,Gf,alphas
      double precision on(4,4),OCL(2,2),OCR(2,2),MChar(2,2),Mneut(4,4)
!--non-universal a terms
      DOUBLE PRECISION a0u11,a0u12,a0u13,a0u21,a0u22,a0u23,a0u31,a0u32,
     $     a0u33

      DOUBLE PRECISION a0d11,a0d12,a0d13,a0d21,a0d22,a0d23,a0d31,a0d32,
     $     a0d33

      DOUBLE PRECISION a0e11,a0e12,a0e13,a0e21,a0e22,a0e23,a0e31,a0e32,
     $     a0e33

      DOUBLE PRECISION a0nu11,a0nu12,a0nu13,a0nu21,a0nu22,a0nu23,
     $     a0nu31,a0nu32,a0nu33

      double precision yuR3(3,3),yuDMR3(3),yuDMR3i(3),
     $                yuev(3,3),yuevT(3,3),UPMNS(3,3)

!--GMSB
      double precision gmsbsusyb, gmsbmess,nhat,gr

!----------------
      common/gmsbinputs/ gmsbsusyb, gmsbmess,gr,nhat
      COMMON/rgeyy/yy_d
      common/unif/ e1,yukgut
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole

      common/loops/ lopt,rhn
      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
      common/charinputs/case, model

      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0

      common/rgeinput_high/ MX, M1X,M2X,M3X,mh10,mh20,
     $     mQ0,mU0,mD0,mE0,mL0,mNU0

      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut

      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG

      common/gauge/alph,Gf,alphas

      common/yukawagut/yugut,ydgut,yegut,ynugut
      common/quitrge/quitint
      common/unifs/fuscale

      common/VCKMparam/ VCKM

!     common/rgeinput/M10,M20,M30,MX,Murka,M_UNI

!-------non universal A terms

      common/nuaterms/a0u11,a0u12,a0u13,a0u21,a0u22,a0u23,a0u31,a0u32,
     $     a0u33,a0d11,a0d12,a0d13,a0d21,a0d22,a0d23,a0d31,a0d32,
     $     a0d33,a0e11,a0e12,a0e13,a0e21,a0e22,a0e23,a0e31,a0e32,
     $     a0e33,a0nu11,a0nu12,a0nu13,a0nu21,a0nu22,a0nu23,
     $     a0nu31,a0nu32,a0nu33

!-------------------------------------------------------------

      external RK4ROUTINE,QMSRK4,smrge,mssmrge,matmult,dag
      external smrgemt,gmsb

!----------------------------------------------------------------
!     WARNING: no statements before these calls
      include 'stdinputs.h'

!----------------------------------------------------------------

!      fuscale = 0

      pi = 4.d0*datan(1.d0)

!      UPMNS(1,3) = ue3

      beta = datan(tanbeta)

      Mpl = 5.d18
      MZ = MZpole 

C       Some definitions
C       ================

      tpl   =  dLog(MX**2.d0/Mpl**2.d0)


      tZ    =  dLog(MX**2.d0/MZ**2.d0)
      mtscale = dlog(MX**2.d0/mtpole**2.d0)


      tq0   =  dLog(MX**2.d0/msusy**2.d0)


      t1 = -2.d0*dlog(MR1/MX)
      t2 = -2.d0*dlog(MR2/MX)
      t3 = -2.d0*dlog(MR3/MX)


      do i = 1,126
         yy(i) = 0.d0
         yy_d(i) = 0.d0
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
      yy_sm(i0 + i)   = yuin(1,i)
      yy_sm(i0+3 + i) = yuin(2,i)
      yy_sm(i0+6 + i) = yuin(3,i)
      enddo inp1

C-------------------------------
C	bottom Yukawa
C-------------------------------

      i0 = 9

      inp2: do i = 1, 3
      yy_sm(i0 + i)   = ydin(1,i)
      yy_sm(i0+3 + i) = ydin(2,i)
      yy_sm(i0+6 + i) = ydin(3,i)
      enddo inp2


C-----------------------------
C     tau Yukawa
C-----------------------------

      i0 = 18

      inp3: do i = 1, 3
      yy_sm(i0 + i)   = yein(1,i)
      yy_sm(i0+3 + i) = yein(2,i)
      yy_sm(i0+6 + i) = yein(3,i)
      enddo inp3


      yy_sm(28) = alph3in
      yy_sm(29) = alph2in
      yy_sm(30) = alph1in

      yy_sm(31) = vevin

!      print*,"vevin = ", vevin
!--------------------------------------------------------------------------
!     RUNNING UP THE SM YUKAWAS AND GAUGE COUPLINGS
!--------------------------------------------------------------------------
      n0 = 31

c$$$      if(runum.eq.1)then
c$$$
c$$$         x2 = mtscale
c$$$         x1 = tz
c$$$
c$$$         h1   =  -1.d-5
c$$$         hmin =  2.d-8
c$$$         eps  =  1.d-6
c$$$!     
c$$$         call RK4ROUTINE(yy_sm,n0,x1,x2,eps,h1,hmin,nok,nbad,smrgemt,
c$$$     .        QMSRK4,check)
c$$$         if(check.eq.100)then
c$$$            flags = 'variable underflow '
c$$$            return
c$$$         endif
c$$$
c$$$      endif


!===============================================================================
!	                 STANDARD MODEL RUNNING ENDS
!===============================================================================
!----------------------------------------
!     matching sm and mssm
!----------------------------------------

      matchdo1: do i = 1,9
      yy(i) = ((yy_sm(i))/dsin(beta))
      enddo matchdo1

      matchdo2: do i = 1,18
      yy(i+9) = ((yy_sm(i+9))/dcos(beta))
      enddo matchdo2


      matchdo3: do i = 1,9
      yy(i+27) = 0.d0
      enddo matchdo3


      matchdo4: do i = 1,72
      yy(i+36) = 0.d0
      enddo matchdo4

!      print*,"max yuk = ", maxval(yy(1:36)),dsqrt(1.d0/(4.d0*pi)),
!     $     maxloc(yy(1:36))

      yy(109) = 0.d0
      yy(110) = 0.d0
      yy(111) = mur
      yy(112) = bmur

      if(maxval(yy(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
         flags = "NPERTYUK"
         return
      endif

!      if(runum.eq.1) then
c$$$         print*,"in runrges yu at mz = ", yy(9)*4.d0*pi,dtan(beta)
c$$$         print*,"yd at mz = ", yy(18)*4.d0*pi
c$$$         print*,"ye at mz = ", yy(27)*4.d0*pi
c$$$         print*,"alpha3in = " , dsqrt(alph3in*16.d0*pi*pi)
c$$$         print*,"alpha2in = " , dsqrt(alph2in*16.d0*pi*pi)
c$$$         print*,"alpha1in = " , dsqrt(alph1in*16.d0*pi*pi)
!      endif
!===========================================================================
!     NOT CHANGED BY D & R 15/12/09 16:58
!     ---------------------------
!     neutino stuff at low energy
!     ---------------------------
!===========================================================================

c     neutrino's patterns from hep-ph/0405172 v5 30 jun 2006

      dm_sol = 7.62d0*10**(-5.d0-18.d0)
      dm_atm = 2.53d0*10.d0**(-3.d0-18.d0)

      mnu_scale = 0.001d0 * 10.d0**(-9.d0)

c     normal hierarchy

      hieif: if(hie.eq.'nor')then

         mnuev(1) = mnu_scale
         mnuev(2) = dsqrt(mnuev(1)**2.d0 + dm_sol)
         mnuev(3) = dsqrt(mnuev(2)**2.d0 + dm_atm)

!     mnuev(2) = mnuev(1) + dsqrt(dm_sol)
!     mnuev(3) = mnuev(2) + dsqrt(dm_atm)

      else hieif

c     inverted hierarchy

         mnuev(1) = dsqrt(mnuev(2)**2.d0 + dm_atm)
         mnuev(2) = dsqrt(mnuev(1)**2.d0 + dm_sol)
         mnuev(3) = mnu_scale

      endif hieif


!     ------------------------------------------------
!     redefinition of upmns so that it is orthogonal
!     ------------------------------------------------


      thl12 = dASin(dsqrt(0.320d0))
      thl23 = dASin(dsqrt(0.49d0))
      thl13 = dASin(dsqrt(0.026d0))
      
      UPMNS(1,1) = dCos(thl12)*dCos(thl13)
      UPMNS(1,2) = dSin(thl12)*dCos(thl13)
      UPMNS(1,3) = dSin(thl13) 

      UPMNS(2,1) = - dSin(thl12)*dCos(thl23) - 
     $     dCos(thl12)*dSin(thl23)*dSin(thl13)
      UPMNS(2,2) = dCos(thl12)*dCos(thl23) -
     $     dSin(thl12)*dSin(thl23)*dSin(thl13)
      UPMNS(2,3) = dCos(thl13)*dSin(thl23)

      UPMNS(3,1) = dSin(thl12)*dSin(thl23) - 
     $     dCos(thl12)*dCos(thl23)*dSin(thl13)
      UPMNS(3,2) = - dCos(thl12)*dSin(thl23) - 
     $     dSin(thl12)*dCos(thl23)*dSin(thl13) 
      UPMNS(3,3) = dCos(thl23)*dCos(thl13)


c$$$	upmns(1,2) = 0.547723d0
c$$$	upmns(2,3) = 0.707107d0
c$$$	upmns(1,3) = ue3
c$$$
c$$$	s13 = upmns(1,3)
c$$$	c13 = dsqrt(1.d0 - s13*s13)
c$$$
c$$$	s23 = upmns(2,3)/c13
c$$$	c23 = dsqrt(1.d0 - s23*s23)
c$$$
c$$$	s12 = upmns(1,2)/c13
c$$$	c12 = dsqrt(1.d0 - s12*s12)
c$$$
c$$$
c$$$	upmns(1,1) = c12*c13
c$$$
c$$$	upmns(2,1) = - (s12*c23 + c12*s23*s13)
c$$$	upmns(2,2) = c12*c23 - s12*s23*s13
c$$$
c$$$	upmns(3,1) = s12*s23 - c12*c23*s13
c$$$	upmns(3,2) = - (c12*s23 + s12*c23*s13)
c$$$	upmns(3,3) = c23*c13

      perm_def_i: do i=1,3
      perm_def_j: do j=1,3
      perm(i,j)=0.d0
      enddo perm_def_j
      enddo perm_def_i


      daunnomeallif: if(hie.eq.'nor')then
         perm(1,3)=1.d0
         perm(2,2)=1.d0
         perm(3,1)=1.d0
      else daunnomeallif
!     perm(1,2)=1.d0
!     perm(2,1)=1.d0
!     perm(3,3)=1.d0
         perm(1,3)=1.d0
         perm(2,1)=1.d0
         perm(3,2)=1.d0
      endif daunnomeallif

c=======================================
c     mnu defintion
c=======================================

      neutr1: do i=1,3
      neutr2: do j=1,3
      mnu(i,j) = 0
      enddo neutr2
      enddo neutr1

      neutr3: do i=1,3
      neutr4: do j=1,3
      neutr5: do l=1,3
      mnu(i,j) = mnu(i,j) + mnuev(l)*upmns(i,l)*upmns(j,l)
      mnu_mz(i,j) = mnu(i,j)
      enddo neutr5
      enddo neutr4
      enddo neutr3


c$$$  !	setting them for the running
c$$$  
      yy(113) = mnu_mz(1,1)
      yy(114) = mnu_mz(1,2)
      yy(115) = mnu_mz(1,3)
      yy(116) = mnu_mz(2,2)
      yy(117) = mnu_mz(2,3)
      yy(118) = mnu_mz(3,3)

      yy(119) = yy_sm(28)
      yy(120) = yy_sm(29)
      yy(121) = yy_sm(30)


      yy(122) = 0.d0
      yy(123) = 0.d0
      yy(124) = 0.d0

      yy(125) = yy_sm(31)*dcos(beta)
      yy(126) = yy_sm(31)*dsin(beta)

c$$$  print*,"vev2 = ", yy(126)
c$$$  print*,"vev1 = ", yy(125)

!============================================================================
!     mssm gauge coupling running from m_susy to MR1
!============================================================================

c------------------------------------------------------------------------
c	           calling the integrator RK4ROUTINE
c------------------------------------------------------------------------

      n0 = 126

      rhn0tz: if(rhn.eq.0)then

         if(model.eq.'GMSB')then
            x2 = -2.d0*dlog(gmsbmess/MX)
         else
            x2 = tpl

!     if(runum.eq.1)then
!     x2 = 2.d0 * dLog(MX/2.d16)
!     else
!     x2 = 2.d0 * dLog(MX/2.24d16)
!     endif

         endif

      else rhn0tz

         if(MR3.eq.0.d0.and.rhn.eq.1)then
            x2 = t2
         else
            x2 = t3
         endif
      endif rhn0tz

      if(runum.eq.1)then
         x1 = tz                !mtscale
      else
         x1 = tz
      endif


      h1   =  -1.d-4
      hmin =  1.d-10
      eps  =  1.d-6


      call RK4ROUTINE(yy,n0,x1,x2,eps,h1,hmin,nok,nbad,mssmrge,QMSRK4,
     $     check)

      if(check.eq.100)then
         flags = 'variable underflow '
         return
      endif

      if(maxval(yy(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
         flags = "NPERTYUK"
         return
      endif

!     print*,"*********Mz --> heaviest rhn scale  done***********"




!============================================================================
!     	mssm gauge coupling running from MR3 to Mgut
!============================================================================

      rhn0gut: if(rhn.eq.0)then
         continue
      else rhn0gut
c     
         if(case.eq.'USD'.or.case.eq.'Rpr')then
            continue
         else


            if(case.eq.'CKM')then
               

!     CKM case:
!     setting Y_nu = Y_up at MR3

!                 if(runum.eq.1)then
               
               Ynui(1,1) = 1.d0 * yy(1)
               Ynui(1,2) = 1.d0 * yy(2)
                    Ynui(1,3) = 1.d0 * yy(3)
                    Ynui(2,1) = 1.d0 * yy(4)
                    Ynui(2,2) = 1.d0 * yy(5)
                    Ynui(2,3) = 1.d0 * yy(6)
                    Ynui(3,1) = 1.d0 * yy(7)
                    Ynui(3,2) = 1.d0 * yy(8)
                    Ynui(3,3) = 1.d0 * yy(9)

!                 else

c$$$                    Ynui(1,1) = 0.d0 * yy(1)
c$$$                    Ynui(1,2) = 0.d0 * yy(2)
c$$$                    Ynui(1,3) = 0.d0 * yy(3)
c$$$                    Ynui(2,1) = 0.d0 * yy(4)
c$$$                    Ynui(2,2) = 0.d0 * yy(5)
c$$$                    Ynui(2,3) = 0.d0 * yy(6)
c$$$                    Ynui(3,1) = 0.d0 * yy(7)
c$$$                    Ynui(3,2) = 0.d0 * yy(8)
c$$$                    Ynui(3,3) = 0.d0 * yy(9)

!                 endif
                    
c$$$              Ynui(1,1) = coeff * yy(1)
c$$$              Ynui(1,2) = coeff * yy(2)
c$$$              Ynui(1,3) = coeff * yy(3)
c$$$              Ynui(2,1) = coeff * yy(4)
c$$$              Ynui(2,2) = coeff * yy(5)
c$$$              Ynui(2,3) = coeff * yy(6)
c$$$              Ynui(3,1) = coeff * yy(7)
c$$$              Ynui(3,2) = coeff * yy(8)
c$$$              Ynui(3,3) = coeff * yy(9)

c$$$              print*,"coeff = ", coeff

           else

              if(case.eq.'MNS')then

C     PMNS case:
!     Ynu = UPMNS*yu(diagonal)


                 i0 = 0

                 MR3u: do i = 1,3

                 yuR3(1,i) = yy(i0 + i)
                 j = 3 + i
                 yuR3(2,i) = yy(i0 + j)
                 k = 6 + i
                 yuR3(3,i) = yy(i0 + k)

              enddo MR3u


              Call CEigensystem(3,yuR3,3,yuDMR3,yuev,3,0)
c$$$
c$$$              call DGEEV( 'V', 'V', 3, yuR3, 3, yuDMR3, yuDMR3i,yuev,3,
c$$$     $             yuevT, 3, WORK, LWORK, INFO )

              Ynui(1,1) = UPMNS(1,1)*yuDMR3(1)
              Ynui(1,2) = UPMNS(1,2)*yuDMR3(1)
              Ynui(1,3) = UPMNS(1,3)*yuDMR3(1)
              Ynui(2,1) = UPMNS(2,1)*yuDMR3(2)
              Ynui(2,2) = UPMNS(2,2)*yuDMR3(2)
              Ynui(2,3) = UPMNS(2,3)*yuDMR3(2)
              Ynui(3,1) = UPMNS(3,1)*yuDMR3(3)
              Ynui(3,2) = UPMNS(3,2)*yuDMR3(3)
              Ynui(3,3) = UPMNS(3,3)*yuDMR3(3)

              ENDIF
           ENDIF
        ENDIF

!---------------------------------------------------

      if(rhn.eq.1)then


         ynuD3(1,1) = Ynui(1,1)
         ynuD3(1,2) = Ynui(1,2)
         ynuD3(1,3) = Ynui(1,3)
         ynuD3(2,1) = Ynui(2,1)
         ynuD3(2,2) = Ynui(2,2)
         ynuD3(2,3) = Ynui(2,3)
         ynuD3(3,1) = Ynui(3,1)
         ynuD3(3,2) = Ynui(3,2)
         ynuD3(3,3) = Ynui(3,3)

         if(MR3.eq.0.d0.and.rhn.eq.1)then

            ynuD3(3,1) = 0.d0
            ynuD3(3,2) = 0.d0
            ynuD3(3,3) = 0.d0

         elseif(MR2.eq.0.d0.and.rhn.eq.1)then

            ynuD3(2,1) = 0.d0
            ynuD3(2,2) = 0.d0
            ynuD3(2,3) = 0.d0

         elseif(MR1.eq.0.d0.and.rhn.eq.1)then

            ynuD3(1,1) = 0.d0
            ynuD3(1,2) = 0.d0
            ynuD3(1,3) = 0.d0

         endif



      else
!         print*," rhn--off"
         ynuD3(1,1) = 0.d0
         ynuD3(1,2) = 0.d0
         ynuD3(1,3) = 0.d0
         ynuD3(2,1) = 0.d0
         ynuD3(2,2) = 0.d0
         ynuD3(2,3) = 0.d0
         ynuD3(3,1) = 0.d0
         ynuD3(3,2) = 0.d0
         ynuD3(3,3) = 0.d0
      endif




       i0 = 27
       match3: do i = 1,3

         yy(i0 + i) = ynuD3(1,i)
         j = 3 + i
         yy(i0 + j) = ynuD3(2,i)
         k = 6 + i
         yy(i0 + k) = ynuD3(3,i)

      enddo match3


!	-------------------------------
!	Calling ODEINT from MR3 to MGUT
!	-------------------------------

 	n0 = 126

        if(MR3.eq.0.d0.and.rhn.eq.1)then
           x1 = t2              !!!!upper limit
        else
           x1 = t3
        endif


        if(model.eq.'GMSB')then
           x2 = -2.d0*dlog(gmsbmess/MX)
        else
           x2 = tpl
        endif

      if(maxval(yy(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
         flags = "NPERTYUK"
         return
      endif

      call RK4ROUTINE(yy,n0,x1,x2,eps,h1,hmin,nok,nbad,mssmrge,
     .     QMSRK4,check)

      if(check.eq.100)then
         flags = 'variable underflow '
         return
      endif

      if(maxval(yy(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
         flags = "NPERTYUK"
         return
      endif

      endif rhn0gut


c$$$
c$$$C     Top Yukawa !!!
c$$$C      ----------------------------------
c$$$
       i0 = 0

       sm01: do i = 1,3

         yugut(1,i) = (4.d0*pi)*yy(i0 + i)
         j = 3 + i
         yugut(2,i) = (4.d0*pi)*yy(i0 + j)
         k = 6 + i
         yugut(3,i) = (4.d0*pi)*yy(i0 + k)

      enddo sm01
c$$$
c$$$
c$$$C     Bottom Yukawa !!!
c$$$C      ----------------------------------
c$$$
         i0 = 9
      sm02: do i = 1,3

         ydgut(1,i) = (4.d0*pi)*yy(i0 + i)
         j = 3 + i
         ydgut(2,i) = (4.d0*pi)*yy(i0 + j)
         k = 6 + i
         ydgut(3,i) = (4.d0*pi)*yy(i0 + k)

      enddo sm02
c$$$
c$$$
c$$$C     Tau Yukawa !!!
c$$$C      ----------------------------------
c$$$
         i0 = 18

         sm03: do i = 1,3

         yegut(1,i) = (4.d0*pi)*yy(i0 + i)
         j = 3 + i
         yegut(2,i) = (4.d0*pi)*yy(i0 + j)
         k = 6 + i
         yegut(3,i) = (4.d0*pi)*yy(i0 + k)

      enddo  sm03


         i0 = 27

         sm04: do i = 1,3

         ynugut(1,i) = (4.d0*pi)*yy(i0 + i)
         j = 3 + i
         ynugut(2,i) = (4.d0*pi)*yy(i0 + j)
         k = 6 + i
         ynugut(3,i) = (4.d0*pi)*yy(i0 + k)

      enddo  sm04


!------------------------------------------------------------------------------
!     Picking unification scale : intersection of a1 and a2
!------------------------------------------------------------------------------

      gmsbhigh: if(model.ne.'GMSB')then

         loopreadp: do i = 1, 126

         yy(i) = yukgut(i)

      enddo loopreadp

      tpl   =  dLog(MX**2.d0/e1**2.d0)

 72   format(1x,A,1x,1pe15.6)

!      print 72, "unification scale = ", e1

      endif gmsbhigh
!----------------------------------------------------------
!	matching  mssm at gut scale:- inputs at high energy
!----------------------------------------------------------

	matchdo11: do i = 1,36
	yy_d(i) = yy(i)
      enddo matchdo11


c-------------

      if(case.eq.'USD'.or.case.eq.'Rpr')then
         continue
      else

         if(case.eq.'CKM')then
!     CKM case:
!     setting Y_nu = Y_up at MR3
!     print*,"CKM case"


c$$$            Ynui(1,1) = yy_d(1)
c$$$            Ynui(1,2) = yy_d(2)
c$$$            Ynui(1,3) = yy_d(3)
c$$$            Ynui(2,1) = yy_d(4)
c$$$            Ynui(2,2) = yy_d(5)
c$$$            Ynui(2,3) = yy_d(6)
c$$$            Ynui(3,1) = yy_d(7)
c$$$            Ynui(3,2) = yy_d(8)
c$$$            Ynui(3,3) = yy_d(9)

            Ynui(1,1) = yy_d(1)
            Ynui(1,2) = yy_d(2)
            Ynui(1,3) = yy_d(3)
            Ynui(2,1) = yy_d(4)
            Ynui(2,2) = yy_d(5)
            Ynui(2,3) = yy_d(6)
            Ynui(3,1) = yy_d(7)
            Ynui(3,2) = yy_d(8)
            Ynui(3,3) = yy_d(9)

         else

            if(case.eq.'MNS')then
C     PMNS case:
!     Ynu = UPMNS*yu(diagonal)
!               print*,"MNS"
                 i0 = 0

                 MR3gut: do i = 1,3

                 yuR3(1,i) = yy_d(i0 + i)
                 j = 3 + i
                 yuR3(2,i) = yy_d(i0 + j)
                 k = 6 + i
                 yuR3(3,i) = yy_d(i0 + k)

              enddo MR3gut

              Call CEigensystem(3,yuR3,3,yuDMR3,yuev,3,0)


c$$$              call DGEEV( 'V', 'V', 3, yuR3, 3, yuDMR3, yuDMR3i,yuev,3,
c$$$     $             yuevT, 3, WORK, LWORK, INFO )

              Ynui(1,1) = UPMNS(1,1)*yuDMR3(1)
              Ynui(1,2) = UPMNS(1,2)*yuDMR3(1)
              Ynui(1,3) = UPMNS(1,3)*yuDMR3(1)
              Ynui(2,1) = UPMNS(2,1)*yuDMR3(2)
              Ynui(2,2) = UPMNS(2,2)*yuDMR3(2)
              Ynui(2,3) = UPMNS(2,3)*yuDMR3(2)
              Ynui(3,1) = UPMNS(3,1)*yuDMR3(3)
              Ynui(3,2) = UPMNS(3,2)*yuDMR3(3)
              Ynui(3,3) = UPMNS(3,3)*yuDMR3(3)


            ENDIF
         ENDIF
      ENDIF

!---------------------------------------------------

      if(rhn.eq.1)then
!         print*," rhn--on"

         if(MR3.eq.0.d0)then

            ynui(3,1) = 0.d0
            ynui(3,2) = 0.d0
            ynui(3,3) = 0.d0

         elseif(MR2.eq.0.d0)then

            ynui(2,1) = 0.d0
            ynui(2,2) = 0.d0
            ynui(2,3) = 0.d0

         elseif(MR1.eq.0.d0)then

            ynui(1,1) = 0.d0
            ynui(1,2) = 0.d0
            ynui(1,3) = 0.d0

         endif

      yy_d(28) = Ynui(1,1)
      yy_d(29) = Ynui(1,2)
      yy_d(30) = Ynui(1,3)
      yy_d(31) = Ynui(2,1)
      yy_d(32) = Ynui(2,2)
      yy_d(33) = Ynui(2,3)
      yy_d(34) = Ynui(3,1)
      yy_d(35) = Ynui(3,2)
      yy_d(36) = Ynui(3,3)




      else
!         print*," rhn--off"

      yy_d(28) = 0.d0
      yy_d(29) = 0.d0
      yy_d(30) = 0.d0
      yy_d(31) = 0.d0
      yy_d(32) = 0.d0
      yy_d(33) = 0.d0
      yy_d(34) = 0.d0
      yy_d(35) = 0.d0
      yy_d(36) = 0.d0

      endif

!----------------------------------

      gmsbBC: if(model.eq.'GMSB')then


       call gmsb(gr,nhat,gmsbsusyb, gmsbmess,yy(121),yy(120),
     $     yy(119),mQ0,mU0,mD0,mL0,mE0,mNU0,M1X,M2X,M3X,m10,m20)

         a0 = 0.d0
         matchdoau1g: do i=1,9
         yy_d(36+i) = a0*yy_d(i)
      enddo matchdoau1g

      matchdoad1g: do i=1,9
      yy_d(45+i) = a0*yy_d(9+i)
      enddo matchdoad1g

      matchdoae1g: do i=1,9
      yy_d(54+i) = a0*yy_d(18+i)
      enddo matchdoae1g

      matchdoanu1g: do i=1,9
      yy_d(63+i) =a0*yy_d(27+i)
      enddo matchdoanu1g

      if(MR3.eq.0.d0.and.rhn.eq.1)then
         yy_d(63 + 7) = 0.d0
         yy_d(63 + 8) = 0.d0
         yy_d(63 + 9) = 0.d0
         endif

      if(MR2.eq.0.d0.and.rhn.eq.1)then
         yy_d(63 + 4) = 0.d0
         yy_d(63 + 5) = 0.d0
         yy_d(63 + 6) = 0.d0
         endif

      if(MR1.eq.0.d0.and.rhn.eq.1)then
         yy_d(63 + 1) = 0.d0
         yy_d(63 + 2) = 0.d0
         yy_d(63 + 3) = 0.d0
         endif

      yy_d(73) = mQ0(1,1)
      yy_d(74) = mQ0(1,2)
      yy_d(75) = mQ0(1,3)
      yy_d(76) = mQ0(2,2)
      yy_d(77) = mQ0(2,3)
      yy_d(78) = mQ0(3,3)

      yy_d(79) = mU0(1,1)
      yy_d(80) = mU0(1,2)
      yy_d(81) = mU0(1,3)
      yy_d(82) = mU0(2,2)
      yy_d(83) = mU0(2,3)
      yy_d(84) = mU0(3,3)

      yy_d(85) = mD0(1,1)
      yy_d(86) = mD0(1,2)
      yy_d(87) = mD0(1,3)
      yy_d(88) = mD0(2,2)
      yy_d(89) = mD0(2,3)
      yy_d(90) = mD0(3,3)

      yy_d(91) = mL0(1,1)
      yy_d(92) = mL0(1,2)
      yy_d(93) = mL0(1,3)
      yy_d(94) = mL0(2,2)
      yy_d(95) = mL0(2,3)
      yy_d(96) = mL0(3,3)

      yy_d(97)  = mE0(1,1)
      yy_d(98)  = mE0(1,2)
      yy_d(99)  = mE0(1,3)
      yy_d(100) = mE0(2,2)
      yy_d(101) = mE0(2,3)
      yy_d(102) = mE0(3,3)

      yy_d(103) = mNU0(1,1)
      yy_d(104) = mNU0(1,2)
      yy_d(105) = mNU0(1,3)
      yy_d(106) = mNU0(2,2)
      yy_d(107) = mNU0(2,3)
      yy_d(108) = mNU0(3,3)


!      mh10 = m10
!      mh20 = m20

!      print*,"m10, m20 = ", m10, m20

      yy_d(109) = mh10
      yy_d(110) = mh20

      yy_d(111) = 10.d0
      yy_d(112) = 10.d0

      yy_d(113) = yy(113)       !mnu_mz(1,1)
      yy_d(114) = yy(114)       !mnu_mz(1,2)
      yy_d(115) = yy(115)       !mnu_mz(1,3)
      yy_d(116) = yy(116)       !mnu_mz(2,2)
      yy_d(117) = yy(117)       !mnu_mz(2,3)
      yy_d(118) = yy(118)       !mnu_mz(3,3)

      yy_d(119) = yy(119)
      yy_d(120) = yy(120)
      yy_d(121) = yy(121)

      yy_d(122) = M1X
      yy_d(123) = M2X
      yy_d(124) = M3X

      yy_d(125) = yy(125)
      yy_d(126) = yy(126)



      else gmsbBC

         matchdoau1: do i=1,9
         yy_d(36+i) = a0*yy_d(i)
      enddo matchdoau1

      matchdoad1: do i=1,9
      yy_d(45+i) = a0*yy_d(9+i)
      enddo matchdoad1

      matchdoae1: do i=1,9
      yy_d(54+i) = a0*yy_d(18+i)
      enddo matchdoae1

      matchdoanu1: do i=1,9
      yy_d(63+i) =a0*yy_d(27+i)
      enddo matchdoanu1


         nuaterms: if(model.eq.'CNUM')then

            yy_d(37) = a0u11
            yy_d(38) = a0u12
            yy_d(39) = a0u13
            yy_d(40) = a0u21
            yy_d(41) = a0u22
            yy_d(42) = a0u23
            yy_d(43) = a0u31
            yy_d(44) = a0u32
            yy_d(45) = a0u33

            yy_d(46) = a0d11
            yy_d(47) = a0d12
            yy_d(48) = a0d13
            yy_d(49) = a0d21
            yy_d(50) = a0d22
            yy_d(51) = a0d23
            yy_d(52) = a0d31
            yy_d(53) = a0d32
            yy_d(54) = a0d33

            yy_d(55) = a0e11
            yy_d(56) = a0e12
            yy_d(57) = a0e13
            yy_d(58) = a0e21
            yy_d(59) = a0e22
            yy_d(60) = a0e23
            yy_d(61) = a0e31
            yy_d(62) = a0e32
            yy_d(63) = a0e33

            yy_d(64) = a0nu11
            yy_d(65) = a0nu12
            yy_d(66) = a0nu13
            yy_d(67) = a0nu21
            yy_d(68) = a0nu22
            yy_d(69) = a0nu23
            yy_d(70) = a0nu31
            yy_d(71) = a0nu32
            yy_d(72) = a0nu33

         endif nuaterms

      if(MR3.eq.0.d0.and.rhn.eq.1)then
         yy_d(63 + 7) = 0.d0
         yy_d(63 + 8) = 0.d0
         yy_d(63 + 9) = 0.d0
         endif

      if(MR2.eq.0.d0.and.rhn.eq.1)then
         yy_d(63 + 4) = 0.d0
         yy_d(63 + 5) = 0.d0
         yy_d(63 + 6) = 0.d0
         endif

      if(MR1.eq.0.d0.and.rhn.eq.1)then
         yy_d(63 + 1) = 0.d0
         yy_d(63 + 2) = 0.d0
         yy_d(63 + 3) = 0.d0
         endif


      yy_d(73) = mQ0(1,1)
      yy_d(74) = mQ0(1,2)
      yy_d(75) = mQ0(1,3)
      yy_d(76) = mQ0(2,2)
      yy_d(77) = mQ0(2,3)
      yy_d(78) = mQ0(3,3)

      yy_d(79) = mU0(1,1)
      yy_d(80) = mU0(1,2)
      yy_d(81) = mU0(1,3)
      yy_d(82) = mU0(2,2)
      yy_d(83) = mU0(2,3)
      yy_d(84) = mU0(3,3)

      yy_d(85) = mD0(1,1)
      yy_d(86) = mD0(1,2)
      yy_d(87) = mD0(1,3)
      yy_d(88) = mD0(2,2)
      yy_d(89) = mD0(2,3)
      yy_d(90) = mD0(3,3)

      yy_d(91) = mL0(1,1)
      yy_d(92) = mL0(1,2)
      yy_d(93) = mL0(1,3)
      yy_d(94) = mL0(2,2)
      yy_d(95) = mL0(2,3)
      yy_d(96) = mL0(3,3)

      yy_d(97)  = mE0(1,1)
      yy_d(98)  = mE0(1,2)
      yy_d(99)  = mE0(1,3)
      yy_d(100) = mE0(2,2)
      yy_d(101) = mE0(2,3)
      yy_d(102) = mE0(3,3)

      yy_d(103) = mNU0(1,1)
      yy_d(104) = mNU0(1,2)
      yy_d(105) = mNU0(1,3)
      yy_d(106) = mNU0(2,2)
      yy_d(107) = mNU0(2,3)
      yy_d(108) = mNU0(3,3)

!      print*,"mh10, mh20 = ", mh10, mh20

      yy_d(109) = mh10
      yy_d(110) = mh20

      if(itcount.eq.1)then
         yy_d(111) = sgnmu*1.1d0*m0
         yy_d(112) = sgnmu*2.d0*m0
         if(m0.eq.0.d0) yy_d(111) = sgnmu*1.1d0 * m12
         if(m0.eq.0.d0) yy_d(112) = sgnmu*2.d0 * m12

      else
         yy_d(111) = yy(111)
         yy_d(112) = yy(112)
      endif


      yy_d(113) = yy(113)       !mnu_mz(1,1)
      yy_d(114) = yy(114)       !mnu_mz(1,2)
      yy_d(115) = yy(115)       !mnu_mz(1,3)
      yy_d(116) = yy(116)       !mnu_mz(2,2)
      yy_d(117) = yy(117)       !mnu_mz(2,3)
      yy_d(118) = yy(118)       !mnu_mz(3,3)

      yy_d(119) = yy(119)
      yy_d(120) = yy(120)
      yy_d(121) = yy(121)

      yy_d(122) = M1X
      yy_d(123) = M2X
      yy_d(124) = M3X

      yy_d(125) = yy(125)
      yy_d(126) = yy(126)



      endif gmsbBC

!-----------------------Integrating Mgut->MR3

      n0 = 126

      rhn0gutd: if(rhn.eq.0)then
         x2 = tq0
      else rhn0gutd
         if(MR3.eq.0.d0.and.rhn.eq.1)then
            x2 = t2             !!!!upper limit
         else
            x2 = t3
         endif
!     ynu last col = 0
      endif rhn0gutd

      if(model.eq.'GMSB')then
         x1 = -2.d0*dlog(gmsbmess/MX)
      else
         x1 = tpl               !!!!lower limit
      endif


!      print*,"at GUT scale mu, bmu = ", yy_d(111),yy_d(112)

      h1   =  1.d-4
      hmin =  1.d-10
      eps  =  1.d-6

      if(maxval(yy_d(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
         flags = "NPERTYUK"
         return
      endif

      
      call RK4ROUTINE(yy_d,n0,x1,x2,eps,h1,hmin,nok,nbad,mssmrge,
     $     QMSRK4,check)
      if(check.eq.100)then
         flags = 'variable underflow '
         return
      endif

      if(maxval(yy_d(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
         flags = "NPERTYUK"
         return
      endif

!     print*,"***************Mgut --> heaviest rhn done ***********"


!============================================================================
!     mssm gauge coupling running from MR3 to MR2
!============================================================================
! setting up Ynu for running : two rows non zero
c-------------
        rhn0susy: if(rhn.eq.0)then
           continue
           else rhn0susy
        mr30: if(MR3.eq.0.d0.and.rhn.eq.1)then
            continue
            else mr30

        if(case.eq.'USD'.or.case.eq.'Rpr')then
           continue
        else

           if(case.eq.'CKM')then
!     CKM case:
!     setting Y_nu = Y_up at MR3

              Ynui(1,1) = yy_d(1)
              Ynui(1,2) = yy_d(2)
              Ynui(1,3) = yy_d(3)
              Ynui(2,1) = yy_d(4)
              Ynui(2,2) = yy_d(5)
              Ynui(2,3) = yy_d(6)
              Ynui(3,1) = yy_d(7)
              Ynui(3,2) = yy_d(8)
              Ynui(3,3) = yy_d(9)

c$$$              Ynui(1,1) = coeff * yy_d(1)
c$$$              Ynui(1,2) = coeff * yy_d(2)
c$$$              Ynui(1,3) = coeff * yy_d(3)
c$$$              Ynui(2,1) = coeff * yy_d(4)
c$$$              Ynui(2,2) = coeff * yy_d(5)
c$$$              Ynui(2,3) = coeff * yy_d(6)
c$$$              Ynui(3,1) = coeff * yy_d(7)
c$$$              Ynui(3,2) = coeff * yy_d(8)
c$$$              Ynui(3,3) = coeff * yy_d(9)

           else

              if(case.eq.'MNS')then
C     PMNS case:
!     Ynu = UPMNS*yu(diagonal)


                 i0 = 0

                 MR3d: do i = 1,3

                 yuR3(1,i) = yy_d(i0 + i)
                 j = 3 + i
                 yuR3(2,i) = yy_d(i0 + j)
                 k = 6 + i
                 yuR3(3,i) = yy_d(i0 + k)

              enddo MR3d


             Call CEigensystem(3,yuR3,3,yuDMR3,yuev,3,0)

c$$$ 
c$$$              call DGEEV( 'V', 'V', 3, yuR3, 3, yuDMR3, yuDMR3i,yuev,3,
c$$$     $             yuevT, 3, WORK, LWORK, INFO )

              Ynui(1,1) = UPMNS(1,1)*yuDMR3(1)
              Ynui(1,2) = UPMNS(1,2)*yuDMR3(1)
              Ynui(1,3) = UPMNS(1,3)*yuDMR3(1)
              Ynui(2,1) = UPMNS(2,1)*yuDMR3(2)
              Ynui(2,2) = UPMNS(2,2)*yuDMR3(2)
              Ynui(2,3) = UPMNS(2,3)*yuDMR3(2)
              Ynui(3,1) = UPMNS(3,1)*yuDMR3(3)
              Ynui(3,2) = UPMNS(3,2)*yuDMR3(3)
              Ynui(3,3) = UPMNS(3,3)*yuDMR3(3)

              ENDIF
           ENDIF
        ENDIF

!---------------------------------------------------

        if(rhn.eq.1)then
!         print*," rhn--on"


           ynuD2(1,1) = Ynui(1,1)
           ynuD2(1,2) = Ynui(1,2)
           ynuD2(1,3) = Ynui(1,3)
           ynuD2(2,1) = Ynui(2,1)
           ynuD2(2,2) = Ynui(2,2)
           ynuD2(2,3) = Ynui(2,3)
           ynuD2(3,1) = 0.d0
           ynuD2(3,2) = 0.d0
           ynuD2(3,3) = 0.d0

           if(MR2.eq.0.d0)then

           ynuD2(2,1) = 0.d0
           ynuD2(2,2) = 0.d0
           ynuD2(2,3) = 0.d0

           endif

           if(MR1.eq.0.d0)then

           ynuD2(1,1) = 0.d0
           ynuD2(1,2) = 0.d0
           ynuD2(1,3) = 0.d0

           endif

        else
!         print*," rhn--off"

           ynuD2(1,1) = 0.d0
           ynuD2(1,2) = 0.d0
           ynuD2(1,3) = 0.d0
           ynuD2(2,1) = 0.d0
           ynuD2(2,2) = 0.d0
           ynuD2(2,3) = 0.d0
           ynuD2(3,1) = 0.d0
           ynuD2(3,2) = 0.d0
           ynuD2(3,3) = 0.d0
        endif
!--------
!----------------------------------------------------
C decoupling of Anu, Ynu (i,3) elements at MR3

         i0 = 27
         matchd031d: do i = 1,3

         yy_d(i0 + i) = ynuD2(1,i)
         j = 3 + i
         yy_d(i0 + j) = ynuD2(2,i)
         k = 6 + i
         yy_d(i0 + k) = ynuD2(3,i)

         enddo matchd031d


         i0 = 63

           if(MR2.eq.0.d0.and.rhn.eq.1)then

            yy_d(i0 + 4)= 0.d0
            yy_d(i0 + 5)= 0.d0
            yy_d(i0 + 6)= 0.d0

           endif

           if(MR1.eq.0.d0.and.rhn.eq.1)then

            yy_d(i0 + 1)= 0.d0
            yy_d(i0 + 2)= 0.d0
            yy_d(i0 + 3)= 0.d0

           endif

         yy_d(i0 + 7) = 0.d0
         yy_d(i0 + 8) = 0.d0
         yy_d(i0 + 9) = 0.d0



C-------------------------------------------------------------------
C calling RK4ROUTINE- running down from MR3 to MR2
C--------------------------------------------------------------------

        n0 = 126

        if(MR2.eq.0.d0.and.rhn.eq.1)then
           x2 = t1
        else
           x2 = t2              !!!!upper limit
        endif

	x1 = t3                 !!!!lower limit


        if(maxval(yy_d(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
           flags = "NPERTYUK"
           return
        endif


        call RK4ROUTINE(yy_d,n0,x1,x2,eps,h1,hmin,nok,nbad,mssmrge,
     .     QMSRK4,check)
	if(check.eq.100)then
           flags = 'variable underflow '
           return
        endif

        if(maxval(yy_d(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
           flags = "NPERTYUK"
           return
        endif

!        print*,"***************MR3 --> MR2 done ***********"


        endif mr30

!============================================================================
!     	mssm  running from MR2 to MR1
!============================================================================
! setting up Ynu for running : one row non zero

        mr20: if(MR2.eq.0.d0.and.rhn.eq.1)then
           continue
           else mr20

        if(case.eq.'USD'.or.case.eq.'Rpr')then
           continue
        else

           if(case.eq.'CKM')then
!     CKM case:
!     setting Y_nu = Y_up at MR3

              Ynui(1,1) = yy_d(1)
              Ynui(1,2) = yy_d(2)
              Ynui(1,3) = yy_d(3)
              Ynui(2,1) = yy_d(4)
              Ynui(2,2) = yy_d(5)
              Ynui(2,3) = yy_d(6)
              Ynui(3,1) = yy_d(7)
              Ynui(3,2) = yy_d(8)
              Ynui(3,3) = yy_d(9)

c$$$              Ynui(1,1) = coeff * yy_d(1)
c$$$              Ynui(1,2) = coeff * yy_d(2)
c$$$              Ynui(1,3) = coeff * yy_d(3)
c$$$              Ynui(2,1) = coeff * yy_d(4)
c$$$              Ynui(2,2) = coeff * yy_d(5)
c$$$              Ynui(2,3) = coeff * yy_d(6)
c$$$              Ynui(3,1) = coeff * yy_d(7)
c$$$              Ynui(3,2) = coeff * yy_d(8)
c$$$              Ynui(3,3) = coeff * yy_d(9)

           else

C     PMNS case:

              if(case.eq.'MNS')then


                 i0 = 0

                 MR2d: do i = 1,3

                 yuR3(1,i) = yy_d(i0 + i)
                 j = 3 + i
                 yuR3(2,i) = yy_d(i0 + j)
                 k = 6 + i
                 yuR3(3,i) = yy_d(i0 + k)

              enddo MR2d

              Call CEigensystem(3,yuR3,3,yuDMR3,yuev,3,0)


c$$$              call DGEEV( 'V', 'V', 3, yuR3, 3, yuDMR3, yuDMR3i,yuev,3,
c$$$     $             yuevT, 3, WORK, LWORK, INFO )

              Ynui(1,1) = UPMNS(1,1)*yuDMR3(1)
              Ynui(1,2) = UPMNS(1,2)*yuDMR3(1)
              Ynui(1,3) = UPMNS(1,3)*yuDMR3(1)
              Ynui(2,1) = UPMNS(2,1)*yuDMR3(2)
              Ynui(2,2) = UPMNS(2,2)*yuDMR3(2)
              Ynui(2,3) = UPMNS(2,3)*yuDMR3(2)
              Ynui(3,1) = UPMNS(3,1)*yuDMR3(3)
              Ynui(3,2) = UPMNS(3,2)*yuDMR3(3)
              Ynui(3,3) = UPMNS(3,3)*yuDMR3(3)

              ENDIF
           ENDIF
        ENDIF

!---------------------------------------------------

        if(rhn.eq.1)then
!         print*," rhn--on"

           ynuD1(1,1) = Ynui(1,1)
           ynuD1(1,2) = Ynui(1,2)
           ynuD1(1,3) = Ynui(1,3)
           ynuD1(2,1) = 0.d0
           ynuD1(2,2) = 0.d0
           ynuD1(2,3) = 0.d0
           ynuD1(3,1) = 0.d0
           ynuD1(3,2) = 0.d0
           ynuD1(3,3) = 0.d0

        else
!           print*," rhn--off"

           ynuD1(1,1) = 0.d0
           ynuD1(1,2) = 0.d0
           ynuD1(1,3) = 0.d0
           ynuD1(2,1) = 0.d0
           ynuD1(2,2) = 0.d0
           ynuD1(2,3) = 0.d0
           ynuD1(3,1) = 0.d0
           ynuD1(3,2) = 0.d0
           ynuD1(3,3) = 0.d0
        endif

!----------------------------------------------------
C decoupling of Anu, Ynu (i,3) elements at MR2


         i0 = 27
         match12: do i = 1,3

         yy_d(i0 + i) = ynuD1(1,i)
         j = 3 + i
         yy_d(i0 + j) = ynuD1(2,i)
         k = 6 + i
         yy_d(i0 + k) = ynuD1(3,i)

        enddo match12


         i0 = 66

         decoup2: do i=1,6
         yy_d(i0 +i) = 0.d0


           if(MR1.eq.0.d0.and.rhn.eq.1)then

            yy_d(63 + 1)= 0.d0
            yy_d(63 + 2)= 0.d0
            yy_d(63 + 3)= 0.d0

           endif

         enddo decoup2


C-------------------------------------------------------------------
C calling RK4ROUTINE- running down from MR2 to MR1
C--------------------------------------------------------------------

        n0 = 126

        if(MR1.eq.0.d0.and.rhn.eq.1)then
        x2 = tq0
        else
         x2 = t1                !!!!upper limit
        endif

	x1 = t2                 !!!!lower limit

        if(maxval(yy_d(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
           flags = "NPERTYUK"
           return
        endif
        

	call RK4ROUTINE(yy_d,n0,x1,x2,eps,h1,hmin,nok,nbad,mssmrge,
     .	QMSRK4,check)

	if(check.eq.100)then
           flags = 'variable underflow '
           return
        endif

        if(maxval(yy_d(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
           flags = "NPERTYUK"
           return
        endif

!        print*,"***************MR2 --> MR1 done ***********"

        endif mr20
!-------------------------------------------at MR1
        mr10: if(MR1.eq.0.d0.and.rhn.eq.1)then
           continue
           else mr10

         ynuD1(1,1) = 0.d0
         ynuD1(1,2) = 0.d0
         ynuD1(1,3) = 0.d0
         ynuD1(2,1) = 0.d0
         ynuD1(2,2) = 0.d0
         ynuD1(2,3) = 0.d0
         ynuD1(3,1) = 0.d0
         ynuD1(3,2) = 0.d0
         ynuD1(3,3) = 0.d0

!----------------------------------------------------
C decoupling of Anu, Ynu (i,3) elements at MR1

         i0 = 63
         decoup1: do i=1,9
         yy_d(i0 + i) = 0.d0
         enddo decoup1

         i0 = 27
         match1MR1: do i = 1,3

         yy_d(i0 + i) = ynuD1(1,i)
         j = 3 + i
         yy_d(i0 + j) = ynuD1(2,i)
         k = 6 + i
         yy_d(i0 + k) = ynuD1(3,i)

        enddo match1MR1

C-------------------------------------------------------------------
C calling RK4ROUTINE- running down from MR1 to Msusy
C--------------------------------------------------------------------

        n0 = 126


	x1 = t1                 !!!!lower limit
        x2 = tq0

        if(maxval(yy_d(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
           flags = "NPERTYUK"
           return
        endif


	call RK4ROUTINE(yy_d,n0,x1,x2,eps,h1,hmin,nok,nbad,mssmrge,
     .	QMSRK4,check)
        
        if(maxval(yy_d(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
           flags = "NPERTYUK"
           return
        endif
      
        
	if(check.eq.100)then
           flags = 'variable underflow '
           return
        endif
!        print*,"***************MR1 --> Msusy done ***********"

        endif mr10

        endif rhn0susy
C-------------------------------------------------------------------------------------
C	Redefining: output
C------------------------------------------------------------------------------------

      newtbeta = yy_d(126)/yy_d(125)

      do i = 1,3
         do j = 1,3
            yuRG(i,j) = 0.d0
            ydRG(i,j) = 0.d0
            yeRG(i,j) = 0.d0

            auRG(i,j) = 0.d0
            adRG(i,j) = 0.d0
            aeRG(i,j) = 0.d0

            msqrg(i,j) = 0.d0
            msurg(i,j) = 0.d0
            msdrg(i,j) = 0.d0

            mslrg(i,j) = 0.d0
            mserg(i,j) = 0.d0

            msnurg(i,j) = 0.d0

            enddo
            enddo

       i0 = 0

      run56: do i = 1,3

         yuRG(1,i) = (4.d0*pi)*yy_d(i0 + i)
         j = 3 + i
         yuRG(2,i) = (4.d0*pi)*yy_d(i0 + j)
         k = 6 + i
         yuRG(3,i) = (4.d0*pi)*yy_d(i0 + k)

         enddo run56


C     Bottom Yukawa !!
C      ----------------------------------

         i0 = 9

       run57: do i = 1,3

         ydRG(1,i) = (4.d0*pi)*yy_d(i0 + i)
         j = 3 + i
         ydRG(2,i) = (4.d0*pi)*yy_d(i0 + j)
         k = 6 + i
         ydRG(3,i) = (4.d0*pi)*yy_d(i0 + k)

         enddo run57

C     Tau Yukawa !!!
C      ----------------------------------

         i0 = 18

         run58: do i = 1,3

         yeRG(1,i) = (4.d0*pi)*yy_d(i0 + i)
         j = 3 + i
         yeRG(2,i) = (4.d0*pi)*yy_d(i0 + j)
         k = 6 + i
         yeRG(3,i) = (4.d0*pi)*yy_d(i0 + k)

         enddo  run58

C     Au-matrix !!!!
C     ------------------

         if(VCKM(1,1).eq.1.d0) then

            AURG(1,1) = yy_d(37)/yy_d(1)
            AURG(1,2) = 0.d0
            AURG(1,3) = 0.d0
            AURG(2,1) = 0.d0
            AURG(2,2) = yy_d(41)/yy_d(5)
            AURG(2,3) = 0.d0
            AURG(3,1) = 0.d0
            AURG(3,2) = 0.d0
            AURG(3,3) = yy_d(45)/yy_d(9)

         else

         i0 = 36
         j0 = 0

            run59: do i = 1,3

            AURG(1,i) = yy_d(i0 + i)/yy_d(j0 + i)
            j = 3 + i
            AURG(2,i) = yy_d(i0 + j)/yy_d(j0 + j)
            k = 6 + i
            AURG(3,i) = yy_d(i0 + k)/yy_d(j0 + k)

      enddo  run59

      endif

C     Ad-matrix !!!!
C     ---------------

      if(VCKM(1,1).eq.1.d0) then

         ADRG(1,1) = yy_d(46)/yy_d(10)
         ADRG(1,2) = 0.d0
         ADRG(1,3) = 0.d0
         ADRG(2,1) = 0.d0
         ADRG(2,2) = yy_d(50)/yy_d(14)
         ADRG(2,3) = 0.d0
         ADRG(3,1) = 0.d0
         ADRG(3,2) = 0.d0
         ADRG(3,3) = yy_d(54)/yy_d(18)

      else


         i0 = 45
         j0 = 9

         run60: do i = 1,3

         ADRG(1,i) = yy_d(i0 + i)/yy_d(j0 + i)
         j = 3 + i
         ADRG(2,i) = yy_d(i0 + j)/yy_d(j0 + j)
         k = 6 + i
         ADRG(3,i) = yy_d(i0 + k)/yy_d(j0 + k)

      enddo run60

      endif


C     Ae-matrix !!
C     -----------------

      if((VCKM(1,1).eq.1.d0).OR.(yeRG(1,3).eq.0.d0)) then


         AERG(1,1) = yy_d(55)/yy_d(19)
         AERG(1,2) = 0.d0
         AERG(1,3) = 0.d0
         AERG(2,1) = 0.d0
         AERG(2,2) = yy_d(59)/yy_d(23)
         AERG(2,3) = 0.d0
         AERG(3,1) = 0.d0
         AERG(3,2) = 0.d0
         AERG(3,3) = yy_d(63)/yy_d(27)




      else

         i0 = 54
         j0 = 18

         run61: do i = 1,3

         AERG(1,i) =  yy_d(i0 + i)/yy_d(j0 + i)
         j = 3 + i
         AERG(2,i) =   yy_d(i0 + j)/yy_d(j0 + j)
         k = 6 + i
         AERG(3,i) =   yy_d(i0 + k)/yy_d(j0 + k)

      enddo run61



      endif

C     mQ-matrix
C     ---------------------

         i0 = 72

         run62: do i = 1,3
             mSQRG(1,i) = yy_d(i0 + i)
         enddo run62

         j0 = i0 + 3

         run63: do j = 1,2
            mSQRG(2,j+1) = yy_d(j0 + j)
         enddo run63

         k = i0 + 6

         mSQRG(3,3) = yy_d(k)

         mSQRG(2,1) = mSQRG(1,2)
         mSQRG(3,1) = mSQRG(1,3)
         mSQRG(3,2) = mSQRG(2,3)


C     mU-matrix
C     -----------------------

         i0 = 78

         run64: do i = 1,3
             mSURG(1,i) = yy_d(i0 + i)
         enddo run64

         j0 = i0 + 3

         run65: do j = 1,2
            mSURG(2,j+1) = yy_d(j0 + j)
         enddo run65

         k = i0 + 6

         mSURG(3,3) = yy_d(k)

         mSURG(2,1) = mSURG(1,2)
         mSURG(3,1) = mSURG(1,3)
         mSURG(3,2) = mSURG(2,3)


C     mD-matrix
C     ------------------------

         i0 = 84

         run66: do i = 1,3
             mSDRG(1,i) = yy_d(i0 + i)
         enddo run66

         j0 = i0 + 3

         run67: do j = 1,2
            mSDRG(2,j+1) = yy_d(j0 + j)
         enddo run67

         k = i0 + 6

         mSDRG(3,3) = yy_d(k)

         mSDRG(2,1) = mSDRG(1,2)
         mSDRG(3,1) = mSDRG(1,3)
         mSDRG(3,2) = mSDRG(2,3)


C     mL-matrix
C     ----------------------

         i0 = 90

         run68: do i = 1,3
             mSLRG(1,i) = yy_d(i0 + i)
         enddo run68

         j0 = i0 + 3

         run69: do j = 1,2
            mSLRG(2,j+1) = yy_d(j0 + j)
         enddo run69

         k = i0 + 6

         mSLRG(3,3) = yy_d(k)

         mSLRG(2,1) = mSLRG(1,2)
         mSLRG(3,1) = mSLRG(1,3)
         mSLRG(3,2) = mSLRG(2,3)



C     mE-matrix
C     --------------------------

         i0 = 96

         run70: do i = 1,3
             mSERG(1,i) = dabs(yy_d(i0 + i))
         enddo run70

         j0 = i0 + 3

         run71: do j = 1,2
            mSERG(2,j+1) = dabs(yy_d(j0 + j))
         enddo run71

         k = i0 + 6

         mSERG(3,3) = dabs(yy_d(k))

         mSERG(2,1) = mSERG(1,2)
         mSERG(3,1) = mSERG(1,3)
         mSERG(3,2) = mSERG(2,3)


C     mNU-matrix
C     -----------------------

         i0 = 102

         run72: do i = 1,3
             mSNURG(1,i) = yy_d(i0 + i)
         enddo run72

         j0 = i0 + 3

         run73: do j = 1,2
            mSNURG(2,j+1) = yy_d(j0 + j)
         enddo run73

         k = i0 + 6

         mSNURG(3,3) = yy_d(k)

         mSNURG(2,1) = mSNURG(1,2)
         mSNURG(3,1) = mSNURG(1,3)
         mSNURG(3,2) = mSNURG(2,3)


C     ------------------------------------------------

         mh1mz  =  yy_d(109)
         mh2mz  =  yy_d(110)
         murge  =  (yy_d(111))
         bmurge =  yy_d(112)

         alph3 = yy_d(119)
         alph2 = yy_d(120)
         alph1 = yy_d(121)

         M1tz   =  yy_d(122)
         M2tz   =  yy_d(123)
         M3tz   =  yy_d(124)

         vev1 = yy_d(125)
         vev2 = yy_d(126)


         RETURN
        end subroutine MSSMRUN

c==========================================================================================
c                              mssm running ends
c==========================================================================================
****f*SuSeFLAV/runrges.f/MSSM_MZ
*  NAME
*    SUBROUTINE MSSM_MZ
*  SYNOPSIS
*    Integrates mssm RGEs from msusy to MZ.
*  FUNCTION
*     Integrates mssm RGEs from msusy to MZ for given inputs at
*     msusy scale.
*  INPUTS
*     MX         -  Reference scale, 10^19(GeV)
*     msusy      -  susy breaking scale.
*     mu_conv    -  Converged value of \mu at msusy.
*     bmur_conv  -  Converged value of b_\mu at msusy.
*
*  RESULT
*     murgemz    -  RGE output: \mu at M_z scale
*     bmurgemz   -  RGE output: b_\mu at M_z scale
*     newtbetamz -  Ratio of vev at MZ from rge running.
*
*  EXAMPLE
*     SUBROUTINE MSSM_MZ(MX,msusy,mu_conv,bmur_conv,
*     $     murgemz,bmurgemz,newtbetamz)
*  NOTES
*     Common blocks and external routines used:
*      common/rgeyy/yy_d
*      common/rgeopt_mz/mt_mz, mb_mz, mtau_mz
*
*      common/rgeoutput_MZ/ mh1mzz,mh2mzz,M1tmz,M2tmz,M3tmz,
*     $     alph1MZ,alph2MZ,alph3MZ,vev1mz,vev2mz,yumz,ydmz,yemz
*
*      common/softout_mat_mz/mSQRGz,mSURGz,mSDRGz,AURGz,ADRGz,
*     $     AERGz,ANURGz,mSLRGz,mSERGz, mSNURGz,ONz,OCLz,OCRz,
*     $     MCharz, MNeutz
*
*      EXTERNAL RK4ROUTINE, QMSRK4,mssmrge
*
*  BUGS
*    ---
*  SEE ALSO
*    MSSMRUN
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
c==========================================================================================
c                              mssm-> mz running begins
c==========================================================================================

      SUBROUTINE MSSM_MZ(MX,msusy,mu_conv,bmur_conv,
     $     murgemz,bmurgemz,newtbetamz,flags)

      IMPLICIT NONE

      integer i0,i,n0,j0,j,nok,nbad,k
      integer check
      CHARACTER*100 flags

      DOUBLE PRECISION yy_d(126),MX,msusy, mu_conv, bmur_conv
      double precision x1,x2,h1,hmin,eps
      DOUBLE PRECISION tz, tq0


!--output
      DOUBLE PRECISION yumz(3,3), ydmz(3,3), yemz(3,3)
      DOUBLE PRECISION alph3MZ, alph2MZ, alph1MZ, vevMZ
      DOUBLE PRECISION AURGz(3,3),ADRGz(3,3),AERGz(3,3)
      DOUBLE PRECISION mSQRGz(3,3), mSURGz(3,3),mSDRGz(3,3)
      DOUBLE PRECISION mSLRGz(3,3),mSNURGz(3,3), mSERGz(3,3)
      DOUBLE PRECISION bmurgemz,M1tmz,M2tmz,M3tmz, mh1mzz,mh2mzz
      DOUBLE PRECISION murgemz, newtbetamz,vev1mz,vev2mz, ANURGz(3,3)

      DOUBLE PRECISION mt_mz, mb_mz, mtau_mz
      DOUBLE PRECISION ONz(4,4),OCLz(2,2),OCRz(2,2),
     $     MCharz(2,2), MNeutz(4,4)
      double precision VCKM(3,3),MZ,pi,mbpole, mtaupole, Mtpole,MZpole

!-------

      common/VCKMparam/ VCKM
      common/rgeyy/yy_d
      common/rgeopt_mz/mt_mz, mb_mz, mtau_mz
      common/sminputs/ mbpole, mtaupole, Mtpole,MZpole

      common/rgeoutput_MZ/ mh1mzz,mh2mzz,M1tmz,M2tmz,M3tmz,
     $     alph1MZ,alph2MZ,alph3MZ,vev1mz,vev2mz,yumz,ydmz,yemz

      common/softout_mat_mz/mSQRGz,mSURGz,mSDRGz,AURGz,ADRGz,
     $     AERGz,ANURGz,mSLRGz,mSERGz, mSNURGz,ONz,OCLz,OCRz,
     $     MCharz, MNeutz

!-----
      EXTERNAL RK4ROUTINE, QMSRK4,mssmrge

!--------
!     WARNING: no statements before these calls
      include 'stdinputs.h'
!----------------------------------------------------------------
      
      MZ = MZpole
      pi = 4.d0*datan(1.d0)

      tq0  =  dLog(MX**2.d0/msusy**2.d0)
      tz   =  dLog(MX**2.d0/MZ**2.d0)

      yy_d(111) = mu_conv
      yy_d(112) = bmur_conv



      n0 = 126

      x2 = tz
      x1 = tq0


      h1   =  1.d-4
      hmin =  1.d-10
      eps  =  1.d-6
!     
      
      if(maxval(yy_d(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
         flags = "NPERTYUK"
         return
      endif

      call RK4ROUTINE(yy_d,n0,x1,x2,eps,h1,hmin,nok,nbad,mssmrge,
     $     QMSRK4,check)

      if(maxval(yy_d(1:36)).gt.dsqrt(1.d0/(4.d0*pi)))then
         flags = "NPERTYUK"
         return
      endif

      if(check.eq.100)then
         flags = 'variable underflow '
         print*,"flags @ mz= ", flags
         return
      endif
c$$$  print*,"***************Msusy ---> Mz done ***********"


!-------------------------------

      newtbetamz = yy_d(126)/yy_d(125)

      do i = 1,3
         do j = 1,3
            yuMZ(i,j) = 0.d0
            ydMZ(i,j) = 0.d0
            yeMZ(i,j) = 0.d0

            auRGz(i,j) = 0.d0
            adRGz(i,j) = 0.d0
            aeRGz(i,j) = 0.d0

            msqrgz(i,j) = 0.d0
            msurgz(i,j) = 0.d0
            msdrgz(i,j) = 0.d0

            mslrgz(i,j) = 0.d0
            msergz(i,j) = 0.d0

            msnurgz(i,j) = 0.d0

         enddo
      enddo



      alph3MZ = yy_d(119)
      alph2MZ = yy_d(120)
      alph1MZ = yy_d(121)

      i0 = 0

      runmz1: do i = 1,3

      yuMZ(1,i) = (4.d0*pi)*yy_d(i0 + i)
      j = 3 + i
      yuMZ(2,i) = (4.d0*pi)*yy_d(i0 + j)
      k = 6 + i
      yuMZ(3,i) = (4.d0*pi)*yy_d(i0 + k)

      enddo runmz1


C     Bottom Yukawa !!
C     ----------------------------------

      i0 = 9

      runmz2: do i = 1,3

      ydMZ(1,i) = (4.d0*pi)*yy_d(i0 + i)
      j = 3 + i
      ydMZ(2,i) = (4.d0*pi)*yy_d(i0 + j)
      k = 6 + i
      ydMZ(3,i) = (4.d0*pi)*yy_d(i0 + k)

      enddo runmz2



C     Tau Yukawa !!!
C     ----------------------------------

      i0 = 18

      runmz3: do i = 1,3

      yeMZ(1,i) = (4.d0*pi)*yy_d(i0 + i)
      j = 3 + i
      yeMZ(2,i) = (4.d0*pi)*yy_d(i0 + j)
      k = 6 + i
      yeMZ(3,i) = (4.d0*pi)*yy_d(i0 + k)

      enddo  runmz3



      if(VCKM(1,1).eq.1.d0) then

         AURGz(1,1) = yy_d(37)/yy_d(1)

         AURGz(2,2) = yy_d(41)/yy_d(5)

         AURGz(3,3) = yy_d(45)/yy_d(9)

      else

         i0 = 36
         j0 = 0

         run59a: do i = 1,3

         AURGz(1,i) = yy_d(i0 + i)/yy_d(j0 + i)
         j = 3 + i
         AURGz(2,i) = yy_d(i0 + j)/yy_d(j0 + j)
         k = 6 + i
         AURGz(3,i) = yy_d(i0 + k)/yy_d(j0 + k)

      enddo  run59a

      endif



C     Ad-matrix !!!!
C     ---------------

      if(VCKM(1,1).eq.1.d0) then

         ADRGz(1,1) = yy_d(46)/yy_d(10)
         ADRGz(1,2) = 0.d0
         ADRGz(1,3) = 0.d0
         ADRGz(2,2) = yy_d(50)/yy_d(14)
         ADRGz(2,1) = 0.d0
         ADRGz(2,3) = 0.d0
         ADRGz(3,1) = 0.d0
         ADRGz(3,2) = 0.d0
         ADRGz(3,3) = yy_d(54)/yy_d(18)

      else


         i0 = 45
         j0 = 9

         run60a: do i = 1,3

         ADRGz(1,i) = yy_d(i0 + i)/yy_d(j0 + i)
         j = 3 + i
         ADRGz(2,i) = yy_d(i0 + j)/yy_d(j0 + j)
         k = 6 + i
         ADRGz(3,i) = yy_d(i0 + k)/yy_d(j0 + k)

      enddo run60a

      endif


C     Ae-matrix !!
C     -----------------

      if((VCKM(1,1).eq.1.d0).OR.(yeMZ(1,3).eq.0.d0)) then


         AERGz(1,1) = yy_d(55)/yy_d(19)
         AERGz(1,2) = 0.d0
         AERGz(1,3) = 0.d0
         AERGz(2,1) = 0.d0
         AERGz(2,2) = yy_d(59)/yy_d(23)
         AERGz(2,3) = 0.d0
         AERGz(3,1) = 0.d0
         AERGz(3,2) = 0.d0

         AERGz(3,3) = yy_d(63)/yy_d(27)




      else

         i0 = 54
         j0 = 18

         run61a: do i = 1,3

         AERGz(1,i) =  yy_d(i0 + i)/yy_d(j0 + i)
         j = 3 + i
         AERGz(2,i) =   yy_d(i0 + j)/yy_d(j0 + j)
         k = 6 + i
         AERGz(3,i) =   yy_d(i0 + k)/yy_d(j0 + k)

      enddo run61a



      endif

C     mQ-matrix
C     ---------------------

      i0 = 72

      run62a: do i = 1,3
      mSQRGz(1,i) = 0.d0
      mSQRGz(1,i) = yy_d(i0 + i)
      enddo run62a

      j0 = i0 + 3

      run63a: do j = 1,2
      mSQRGz(2,j+1) = 0.d0
      mSQRGz(2,j+1) = yy_d(j0 + j)
      enddo run63a

      k = i0 + 6

      mSQRGz(3,3) = yy_d(k)

      mSQRGz(2,1) = mSQRGz(1,2)
      mSQRGz(3,1) = mSQRGz(1,3)
      mSQRGz(3,2) = mSQRGz(2,3)


C     mU-matrix
C     -----------------------

      i0 = 78

      run64a: do i = 1,3
      mSURGz(1,i) = yy_d(i0 + i)
      enddo run64a

      j0 = i0 + 3

      run65a: do j = 1,2
      mSURGz(2,j+1) = yy_d(j0 + j)
      enddo run65a

      k = i0 + 6

      mSURGz(3,3) = yy_d(k)

      mSURGz(2,1) = mSURGz(1,2)
      mSURGz(3,1) = mSURGz(1,3)
      mSURGz(3,2) = mSURGz(2,3)


C     mD-matrix
C     ------------------------

      i0 = 84

      run66a: do i = 1,3
      mSDRGz(1,i) = yy_d(i0 + i)
      enddo run66a

      j0 = i0 + 3

      run67a: do j = 1,2
      mSDRGz(2,j+1) = yy_d(j0 + j)
      enddo run67a

      k = i0 + 6

      mSDRGz(3,3) = yy_d(k)

      mSDRGz(2,1) = mSDRGz(1,2)
      mSDRGz(3,1) = mSDRGz(1,3)
      mSDRGz(3,2) = mSDRGz(2,3)


C     mL-matrix
C     ----------------------

      i0 = 90

      run68a: do i = 1,3
      mSLRGz(1,i) = dabs(yy_d(i0 + i))
      enddo run68a

      j0 = i0 + 3

      run69a: do j = 1,2
      mSLRGz(2,j+1) = dabs(yy_d(j0 + j))
      enddo run69a

      k = i0 + 6

      mSLRGz(3,3) = dabs(yy_d(k))

      mSLRGz(2,1) = mSLRGz(1,2)
      mSLRGz(3,1) = mSLRGz(1,3)
      mSLRGz(3,2) = mSLRGz(2,3)



C     mE-matrix
C     --------------------------

      i0 = 96

      run70a: do i = 1,3
      mSERGz(1,i) = yy_d(i0 + i)
      enddo run70a

      j0 = i0 + 3

      run71a: do j = 1,2
      mSERGz(2,j+1) = yy_d(j0 + j)
      enddo run71a

      k = i0 + 6

      mSERGz(3,3) = yy_d(k)

      mSERGz(2,1) = mSERGz(1,2)
      mSERGz(3,1) = mSERGz(1,3)
      mSERGz(3,2) = mSERGz(2,3)


C     mNU-matrix
C     -----------------------

      i0 = 102

      run72a: do i = 1,3
      mSNURGz(1,i) = yy_d(i0 + i)
      enddo run72a

      j0 = i0 + 3

      run73a: do j = 1,2
      mSNURGz(2,j+1) = yy_d(j0 + j)
      enddo run73a

      k = i0 + 6

      mSNURGz(3,3) = yy_d(k)

      mSNURGz(2,1) = mSNURGz(1,2)
      mSNURGz(3,1) = mSNURGz(1,3)
      mSNURGz(3,2) = mSNURGz(2,3)


C     ------------------------------------------------

      M1tmz = 0.d0
      M2tmz = 0.d0
      M3tmz = 0.d0
      murgemz = 0.d0
      bmurgemz = 0.d0

      mh1mzz = 0.d0
      mh2mzz = 0.d0

      mh1mzz  =  yy_d(109)
      mh2mzz  =  yy_d(110)

      murgemz  =  (yy_d(111))
      bmurgemz =  yy_d(112)

      M1tmz   =  yy_d(122)
      M2tmz   =  yy_d(123)
      M3tmz   =  yy_d(124)

      vev1mz = yy_d(125)
      vev2mz = yy_d(126)


      vevMZ = dsqrt(vev1mz**2.d0 + vev2mz**2.d0)

!---------------------------------

      mt_mz = yuMZ(3,3)*vev2mz /dsqrt(2.d0)
      mb_mz = ydMZ(3,3)*vev1mz /dsqrt(2.d0)
      mtau_mz = yeMZ(3,3)*vev1mz /dsqrt(2.d0)

!-----------------------------

      RETURN
      END SUBROUTINE MSSM_MZ
!============================================================================
