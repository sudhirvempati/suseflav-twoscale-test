****f* SuSeFLAV/oneloopgauge.f 
*  NAME
*    subroutine strongcoupling,emcoupling,gluinose
*
*  SYNOPSIS
*    Computes One loop correction(analytical) alphas,alphaem and gluino 
*
*  FUNCTION
*     The routine calcultes one loop correction to alphas,alphaem and gluino.
*     hep-ph/ 9903404, hep-ph/0002213  - analytical expression for light higgs and 
*     CP- even higgs
*
*  INPUTS
*     p                              - external momenta
*     q                              - scale
*     alph3                          -  g3/(16 * pi^2)
*     mt,mb,mtau                     -  running masses of top, bottom and tau
*     M3tz                           - Gluino mass 
*     SUegg   =  6 eigenvalues of UP-Squark mass matrix.
*     SDegg   =  6 eigenvalues of Down-Squark mass matrix.
*     SLegg   =  6 eigenvalues of slepton mass matrix.
*     Ceg     =  2 singular values of the chargino Mass Matrix
*     mhpmsq  = mass of charged higgs ^2
*
*  RESULT
*  
*  EXAMPLE
*
*      subroutine emcoupling(q,alphaDR,mt,SUegg,SDegg,SLegg,Ceg,
*     $     mHpmsq,alphaem,delalphem)
*
*      subroutine strongcoupling(q,g3,mt,M3tz,SUegg,SDegg,
*     $     alphas1,delalphas)
*
*      subroutine gluinose(p,q,g3,mt,mB,SUegg,SDegg,
*     $     deltagluino)
*
*     
*  NOTES
*
*  BUGS
*    ---
*  SEE ALSO
*
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!-------------------------------------------------------------------------
!                  Gauge Coupling Radiative Correction
!-------------------------------------------------------------------------

      subroutine strongcoupling(q,g3,mt,M3tz,SUegg,SDegg,
     $     alphas1,delalphas)
      
      integer i, j
      double precision mT,alph,Gf,alphasin
      double precision delalphas,SUegg(6),SDegg(6),mgluino,M3tz,g3
      double precision mup(3,2), mdown(3,2),alphas1,q,pi
      data mup/ 6 * 0.d0/, mdown/ 6 * 0.d0/
      common/gauge/alph,Gf,alphasin

!--------------------------------------------
      include "stdinputs.h"
!--------------------------------------------

      pi = 4.d0 * datan(1.d0)
      
      mgluino = M3tz
      
c$$$      mup(1,1) = dsqrt(SUegg(6))
c$$$      mup(1,2) = dsqrt(SUegg(5))
c$$$      mup(2,1) = dsqrt(SUegg(4))
c$$$      mup(2,2) = dsqrt(SUegg(3))
c$$$      mup(3,1) = dsqrt(SUegg(2))
c$$$      mup(3,2) = dsqrt(SUegg(1))
c$$$      
c$$$      mdown(1,1) = dsqrt(SDegg(6)) 
c$$$      mdown(1,2) = dsqrt(SDegg(5))
c$$$      mdown(2,1) = dsqrt(SDegg(4))
c$$$      mdown(2,2) = dsqrt(SDegg(3))
c$$$      mdown(3,1) = dsqrt(SDegg(2))
c$$$      mdown(3,2) = dsqrt(SDegg(1))

      mup(1,1) = dsqrt(SUegg(5))
      mup(1,2) = dsqrt(SUegg(6))
      mup(2,1) = dsqrt(SUegg(3))
      mup(2,2) = dsqrt(SUegg(4))
      mup(3,1) = dsqrt(SUegg(1))
      mup(3,2) = dsqrt(SUegg(2))
      
      mdown(1,1) = dsqrt(SDegg(5)) 
      mdown(1,2) = dsqrt(SDegg(6))
      mdown(2,1) = dsqrt(SDegg(3))
      mdown(2,2) = dsqrt(SDegg(4))
      mdown(3,1) = dsqrt(SDegg(1))
      mdown(3,2) = dsqrt(SDegg(2))
                 
!------------------------------------------ 
      alphas = alphasin !g3 * g3 / ( 4.d0 * pi)
            
      delalphas = (alphas / (2.d0 * pi)) * (0.5d0 - (2.d0/ 3.d0) *
     $     dlog(mt/ q) - 2.d0 * dlog(dabs(mGluino) / q))
      
!      print*,"delalphas = ", delalphas
      
      loopi: do i = 1, 3
      loopj: do j = 1, 2
      
      delalphas = delalphas - (alphas / (12.d0 * pi)) * 
     $     (dlog(mup(i, j)/ q) + dlog(mdown(i, j)/ q))
      
      enddo loopj
      enddo loopi

!      print*,"delalphas = ", delalphas
      
!      print*,"alphas = ", alphas

      alphas1 = alphas/ (1.d0 - delalphas)

!------------------------------------------
      
      return     
      end subroutine strongcoupling

!======================================================================================

!-----------------------------------------------------------------------------------
c$$$   Does SUSY (and other) threshold corrections to alphaEm - returns alpha in
c$$$   DRbar scheme at scale Q. From hep-ph/9606211. Input empirical value of
c$$$   alpha at MZ external momentum....
!-----------------------------------------------------------------------------------
      
      subroutine emcoupling(q,alphaDR,mt,SUegg,SDegg,SLegg,Ceg,
     $     mHpmsq,alphaem,delalphem)

      
      double precision mT,  alphaDR
      double precision SUegg(6),SDegg(6),SLegg(6),Ceg(2)
      double precision mup(3,2), mdown(3,2), ml(3,2),q,mHpmsq
      data mup/ 6 * 0.d0/, mdown/ 6 * 0.d0/, ml/ 6 * 0.d0/
      double precision delalphem, mchargino(2)
      data mchargino/ 0.d0,  0.d0/
      double precision deltasm, deltasusy,alphaem

      double precision pi,MWpole,MW,alph,Gf,alphas
      
      common/mwpole/ MWpole
      common/gauge/alph,Gf,alphas

!----------------------------------------------------      
      include "stdinputs.h"
      
!----------------------------------------------------

      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      
c$$$      mup(1,2) = dsqrt(SUegg(6)) 
c$$$      mup(1,1) = dsqrt(SUegg(5))
c$$$      mup(2,2) = dsqrt(SUegg(4))
c$$$      mup(2,1) = dsqrt(SUegg(3))
c$$$      mup(3,2) = dsqrt(SUegg(2))
c$$$      mup(3,1) = dsqrt(SUegg(1))
c$$$      
c$$$      mdown(1,2) = dsqrt(SDegg(6)) 
c$$$      mdown(1,1) = dsqrt(SDegg(5))
c$$$      mdown(2,2) = dsqrt(SDegg(4))
c$$$      mdown(2,1) = dsqrt(SDegg(3))
c$$$      mdown(3,2) = dsqrt(SDegg(2))
c$$$      mdown(3,1) = dsqrt(SDegg(1))
c$$$      
c$$$      ml(1,2) = dsqrt(SLegg(6))
c$$$      ml(1,1) = dsqrt(SLegg(5))
c$$$      ml(2,2) = dsqrt(SLegg(4))
c$$$      ml(2,1) = dsqrt(SLegg(3))
c$$$      ml(3,2) = dsqrt(SLegg(2))
c$$$      ml(3,1) = dsqrt(SLegg(1))

      mup(1,1) = dsqrt(SUegg(5))
      mup(1,2) = dsqrt(SUegg(6))
      mup(2,1) = dsqrt(SUegg(3))
      mup(2,2) = dsqrt(SUegg(4))
      mup(3,1) = dsqrt(SUegg(1))
      mup(3,2) = dsqrt(SUegg(2))
      
      mdown(1,1) = dsqrt(SDegg(5)) 
      mdown(1,2) = dsqrt(SDegg(6))
      mdown(2,1) = dsqrt(SDegg(3))
      mdown(2,2) = dsqrt(SDegg(4))
      mdown(3,1) = dsqrt(SDegg(1))
      mdown(3,2) = dsqrt(SDegg(2))

      ml(1,2) = dsqrt(SLegg(6))
      ml(1,1) = dsqrt(SLegg(5))
      ml(2,2) = dsqrt(SLegg(4))
      ml(2,1) = dsqrt(SLegg(3))
      ml(3,2) = dsqrt(SLegg(2))
      ml(3,1) = dsqrt(SLegg(1))
      
      mchargino(1) = Ceg(1)
      mchargino(2) = Ceg(2)
      
      mHpm = dsqrt(mHpmsq)


!-------------------------------------
          
c$$$      deltasm =  (- 7.d0 * dlog(MW/ q) + 
c$$$     $     (16.d0/ 9.d0) * dlog(mt/ q))

      deltasm =  (- 1.d0/3.d0 + 
     $     (16.d0/ 9.d0) * dlog(mt/ q))
      
      deltasusy = (dlog(mHpm/ q)/ 3.d0) + (4.d0/ 9.d0) * 
     $     (dlog(mup(1,1)/ q) + dlog(mup(1,2)/ q) + 
     $     dlog(mup(2,1)/ q) + dlog(mup(2,2)/ q) + 
     $     dlog(mup(3,1)/ q) + dlog(mup(3,2)/ q)) + (1.d0/ 9.d0) *
     $     (dlog(mdown(1,1)/ q) + dlog(mdown(1,2)/ q) + 
     $     dlog(mdown(2,1)/ q) + dlog(mdown(2,2)/ q) + 
     $     dlog(mdown(3,1)/ q) + dlog(mdown(3,2)/ q)) + (1.d0/ 3.d0) *
     $     (dlog(ml(1,1)/ q) + dlog(ml(1,2)/ q) + 
     $     dlog(ml(2,1)/ q) + dlog(ml(2,2)/ q) + 
     $     dlog(ml(3,1)/ q) + dlog(ml(3,2)/ q)) + (4.d0/ 3.d0) *
     $     (dlog(dabs(mchargino(1))/ q) + dlog(dabs(mchargino(2))/ q))
      

      alphaem = alph
      
      
      delalphem = 0.d0 - 
     $     ((alphaem/(2.d0 * pi)) * (deltasm + deltasusy))
      
c$$$      print*,"deltasm = ", deltasm, " deltasusy = ", deltasusy
c$$$      print*,"delalphem = ", delalphem
c$$$      print*,"alphaem = ", alphaem

      alphaem = alphaem/ (1.d0 - delalphem)
      
!------------------------------------------------
      return
      end subroutine emcoupling

!-----------------------------------------------------------------------------

      
!=============================================================================

      subroutine gluinose(p,q,g3,mt,mB,SUegg,SDegg,
     $     deltagluino)

      
      implicit none
      integer i
      double precision mT, mB
      double precision p,q,g3,deltagluino
      
      DOUBLE PRECISION SUegg(6),SDegg(6)

      
      double precision thetat,thetab,thetatau
      double precision thetac,thetas,thetamu
      double precision thetau,thetad,thetae
      
      double precision msup(2),msdown(2),mscharm(2)
      double precision msstrange(2),msbot(2),mstop(2)
      double precision b1mstop(2),b1msup(2),b1mscharm(2)
      double precision b1msdown(2),b1msstrange(2),b1msbot(2)
      double precision b0mstop(2),b0msbot(2),pi
      
      double precision delta,delsusy(2)
      data delsusy/ 2 * 0.d0/, delta/ 0.d0/

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae
      
!---------------------------------------------------------------------
      
      external b0,b1
      
      include 'stdinputs.h'
      
!---------------------------------------------------------------------

      pi = 4.d0 * datan(1.d0)

      msup(2)    = dsqrt(SUegg(6)) 
      msup(1)    = dsqrt(SUegg(5))
      mscharm(2) = dsqrt(SUegg(4))
      mscharm(1) = dsqrt(SUegg(3))
      mstop(2)   = dsqrt(SUegg(2))
      mstop(1)   = dsqrt(SUegg(1))
      
      msdown(2)    = dsqrt(SDegg(6)) 
      msdown(1)    = dsqrt(SDegg(5))
      msstrange(2) = dsqrt(SDegg(4))
      msstrange(1) = dsqrt(SDegg(3))
      msbot(2)     = dsqrt(SDegg(2))
      msbot(1)     = dsqrt(SDegg(1))
      
c$$$      print*,"msup(1), msup(2) = ", msup(1), msup(2)
c$$$      print*,"mscharm(1), mscharm(2) = ", mscharm(1), mscharm(2)
c$$$      print*,"mstop(1), mstop(2) = ", mstop(1), mstop(2)

!------------------------------------------------
      loopi: do i = 1, 2
      
      call b1(p,muq,msup(i),q,b1msup(i))
      call b1(p,mc,mscharm(i),q,b1mscharm(i))
      call b1(p,mt,mstop(i),q,b1mstop(i))
      call b1(p,md,msdown(i),q,b1msdown(i))
      call b1(p,ms,msstrange(i),q,b1msstrange(i))
      call b1(p,mb,msbot(i),q,b1msbot(i))
      
      call b0(p,mt,mstop(i),q,b0mstop(i))
      call b0(p,mb,msbot(i),q,b0msbot(i))
      
      enddo loopi
      
!--------------------------------------------------

      delta = 15.d0 + 9.d0 * dlog((q*q)/(p*p))
      
!      print*,"delta1 = " , delta
     
C     Quark/squark correction
!----------------------------------
      
      delsusy(1) = (b1msup(1) + b1mscharm(1) + b1mstop(1) + 
     $     b1msdown(1) + b1msstrange(1) + b1msbot(1))

      delsusy(2) = (b1msup(2) + b1mscharm(2) + b1mstop(2) + 
     $     b1msdown(2) + b1msstrange(2) + b1msbot(2))
      
      
      delta = delta - (delsusy(1) + delsusy(2))
      

c$$$      print*,"delta, delsusy(1), delsusy(2) = ", delta, delsusy(1), 
c$$$     $     delsusy(2)

C     Third family mixing contribution
!-------------------------------------------
      
      delta = delta - mt * sin(2.d0 * thetat) * 
     $     ((b0mstop(1) - b0mstop(2))/p)
      
      delta = delta - mb * sin(2.d0 * thetab) * 
     $     ((b0msbot(1) - b0msbot(2))/p)

c$$$      delta = delta + mt * sin(2.d0 * thetat) * 
c$$$     $     ((b0mstop(1) - b0mstop(2))/p)
      
!      print*,"deltatop = ", delta

c$$$      delta = delta + mb * sin(2.d0 * thetab) * 
c$$$     $     ((b0msbot(1) - b0msbot(2))/p)

!      print*,"deltabtm = ", delta

!      deltagluino = - delta * p * g3 * g3/(16.d0 * pi * pi)

      deltagluino =  delta * g3 * g3/(16.d0 * pi * pi)
     
!---------------------------------------------- 
      return     
      end subroutine gluinose

!======================================================================================
      
      
