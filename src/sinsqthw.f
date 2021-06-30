****f* SuSeFLAV/sinsqthw.f 
*  NAME
*    subroutine 
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
*
*     alph3                          -  g3/(16 * pi^2)
*     mt,mb,mtau                     -  running masses of top, bottom and tau
*     yuRG,ydRG,yeRG                 - (3 X 3) Yukawas
*     AURG                           - (3 X 3) Trilinear couplings
*     pizzT,piwwT                    - self energy of W and Z bosons at M_z
*     modmu                          - modulus of the \mu paramter 
*     vev1,vev2                      - vacuum expectation values of the two 
*                                      higgs doublet fields
*     M3t                            - Gaugino mass at msusy
*     tanbeta                        - the ratio of the vevs of the 
*                                      two Higgs doublet fields.
*     SUegg   =  6 eigenvalues of UP-Squark mass matrix.
*
*  RESULT
*  
*  EXAMPLE
*
*     subroutine strongcoupling(q,g3,mt,mB,mTau,M3tz,SUegg,SDegg,
*     $     alphas1,delalphas)
*     
*     subroutine emcoupling(q,alphaDR,mt,mB,mTau,SUegg,SDegg,SLegg,Ceg,
*     $     mHpmsq,alphaem,delalphem)
*     
*     subroutine gluinose(p,q,g3,mt,mB,mTau,tanbeta,mSQRG,mSDRG,mSURG,
*     $     mSLRG,mSERG,AURG,ADRG,AERG,SUegg,SDegg,SLegg,M3tz,
*     $     modmu,deltagluino)
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
!====================================================================
      SUBROUTINE S2ThetaEff(pizgMZ, pizg0, sinsqthwDR, 
     $     alphaem,sinsqtheff)

      IMPLICIT NONE
      DOUBLE PRECISION sinsqthwDR, sinsqtheff,kl,sw2,cw2,alphaem
      DOUBLE PRECISION Vlmz2, alphhat,pizgMZ, pizg0
      DOUBLE PRECISION refx,gx,MWpole, MZpole, MW,MZ
      double precision mbpole, mtaupole, Mtpole,pi
!--------
      
      common/mwpole/ MWpole
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      
!----------------------------
      INCLUDE 'stdinputs.h'
!------
      
      pi = 4.d0 * datan(1.d0)
      MW = MWpole
      MZ = MZpole

      alphhat = alphaem 


      sw2 = 1.d0 - (MW*MW)/(MZ*MZ) !sinsqthwDR!
      cw2 = 1- sw2

      Vlmz2 = 0.5d0 * refx(1.d0/cw2) + 
     $     4.d0 * (1.d0 - sinsqthwDR) * gx(1.d0/cw2) -
     $     (1.d0 - 6.d0 * sinsqthwDR + 8.d0 * sinsqthwDR**2.d0) * 
     $     refx(1.d0)/(4.d0 * (1.d0 - sinsqthwDR))

      kl = 1.d0 + (alphhat * (1.d0 - sinsqthwDR) * dlog(cw2)/
     $     (pi * sinsqthwDR)) - 
     $     alphhat * Vlmz2 / (4.d0 * pi * sinsqthwDR) +
     $     dsqrt((1.d0 - sinsqthwDR)/sinsqthwDR) * 
     $     (pizgMZ - pizg0)/ (MZ*MZ)

      sinsqtheff = sinsqthwDR !* kl

      RETURN
      END SUBROUTINE
!--------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION refx(x)

      IMPLICIT NONE
      DOUBLE PRECISION x,pi
      DOUBLE PRECISION DDILOG

      pi = 4.d0 * datan(1.d0)
      
      refx = 2.d0/x + 3.5d0 - 
     $     (3.d0 + 2.d0/x)*dlog(x) +
     $     (1.d0 + 1.d0/x)**2.d0 * ( 2.d0 * DDILOG(1.d0/(1.d0 + x)) -
     $     pi*pi/3.d0 + (dlog(1.d0 + x))**2.d0)
      
      RETURN
      END

!------------------------------------------------------------------------

      DOUBLE PRECISION FUNCTION gx(x)
      IMPLICIT NONE
      DOUBLE PRECISION x, y
      
      y = dsqrt(x/(4.d0 - x))

      gx = (1.d0/x + 0.5d0) * (datan(y)/y - 1.d0) +
     $     (9.d0/8.d0) + (0.5d0/x) - 
     $     (1.d0 + 0.5d0/x) * 4.d0 * (datan(y))**2.d0/x 

      RETURN
      END
!-----------------------------------------------------------------------

      RECURSIVE SUBROUTINE S2ThetaW(q,g,gp,Neg,Ceg,SLegg,SNegg,ON,OCL,
     $     OCR,sinthwold, rhold, delalphhat, 
     $     alphaMZ, pizzTMZ, piwwTMW, piwwT0, mh0sq,trysinsq,flags2tw)


      IMPLICIT NONE
      INTEGER try2,trysinsq
      character*100 flags2tw
      DOUBLE PRECISION alphhat, pizzTMZ, piwwTMW, piWWT0
      DOUBLE PRECISION mh0sq, rhold, delr, sinthwold, sinthw_cor
      DOUBLE PRECISION alphaMZ,xt, delr1, delr2, MWpole
      DOUBLE PRECISION mbpole, mtaupole, Mtpole, rhonew, MZpole
      DOUBLE PRECISION g,q,gp,Neg(4),Ceg(2),ON(4,4),OCL(2,2),
     $     OCR(2,2) 
      DOUBLE PRECISION SLegg(6),SNegg(3),deltaztot,deltavtot,
     $     delalphhat

      double precision m0, m12, m10, m20, sgnmu, tanbeta, a0
      DOUBLE PRECISION alph,Gf,alphas,pi,MW,MZ
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
      common/gauge/alph,Gf,alphas
      common/mwpole/ MWpole


      double precision rho2, DELTAVB1
      EXTERNAL delrho, DELTAZ, deltav

!----------
      include 'stdinputs.h'
!----------

      
      pi = 4.d0 * datan(1.d0)
      
      MW = MWpole
      MZ = MZpole


!      print*, "alph = ", alph, 1.d0/alph

      alphhat = alph/(1.d0 - delalphhat)

!     alphhat = (1.d0/127.934d0)/(1.d0 - delalphhat)

      xt = 3.d0 * Gf * mtpole**2.d0 / (8.d0 * pi*pi * dsqrt(2.d0))

      try2 = 1


      call delrho(rhold,sinthwold**2.d0, alphhat, rhonew,
     $     alphaMZ, pizzTMZ, piwwTMW, piWWT0, mh0sq,try2)


      call DELTAZ(q,g,gp,Neg,Ceg,SLegg,SNegg,ON,OCL,OCR,
     $     deltaztot)
      

      call deltav(sinthwold**2.d0, alphhat, q,g,gp,Neg,Ceg,
     $     Slegg,SNegg,ON,OCL,OCR,deltavtot)


      delr1 = rhonew * (piwwT0)/MW**2.d0 - pizzTMZ/MZ**2.d0 

      delr2 = alphhat * alphas/(4.d0 * pi * pi * sinthwold**2.d0 *
     $     (dcos(dasin(sinthwold)))**2.d0 ) * (
     $     -2.145d0 * (Mtpole/MZ)**2.d0 + 0.575d0*dlog(MTpole/MZ)  - 
     $     2.24d0 - 0.144d0 * (MZ/Mtpole)**2.d0) + 
     $     xt**2.d0/3.d0 * rho2(dsqrt(mh0sq)/mtpole) * 
     $     (dcos(alphaMZ)/ dsin(datan(tanbeta)))**2.d0 *
     $     (1 - delr1)*rhonew

      delr = delr1 + delr2 + deltavtot + deltaztot + 
     $     DELTAVB1(alphhat,rhonew,sinthwold**2.d0,
     $     (1.d0 - (MW*MW)/(MZ*MZ)))

 


      sinthw_cor = dsqrt((1.d0 - (dsqrt(1.d0 - (4.d0 * pi * alphhat/
     $     (dsqrt(2.d0) * MZ*MZ *Gf*(1.d0 - (delr)))))))/2.d0)
   


      if((((1.d0 - (4.d0 * pi * alphhat/
     $     (dsqrt(2.d0) * MZ*MZ *Gf*(1.d0 - (delr))))))/2.d0)
     $     .lt.0.d0)then

         flags2tw = ' sin^2\theta_w nonconvergent'
         sinthw_cor = 0.48d0 ! dsqrt(sinsqtw) !<---nan protection
!         print*,flags2tw
         endif

         if(delr.ge.1.d0)then

            flags2tw = ' sin^2\theta_w nonconvergent'
            sinthw_cor = 0.48d0 ! dsqrt(sinsqtw) !<---nan protection
!     print*,flags2tw
         endif


 29   format(1pe11.4,2x,1pe11.4,2x,1pe11.4,2x,1pe11.4)   

      if(trysinsq.eq.20)then

!         WRITE(197,29) m0, m12,a0, tanbeta
         flags2tw = ' sin^2\theta_w nonconvergent'
         sinthwold = 0.48d0     !dsqrt(sinsqtw)
         return
      endif

      if(((1.d0 - min(sinthwold,sinthw_cor)/max(sinthwold,sinthw_cor))
     $     .le.2.d-6).and.(trysinsq.lt.20))then

      sinthwold = sinthw_cor

!      sinthwold = dsqrt(.23d0)
     
      else
      sinthwold = sinthw_cor
      trysinsq = trysinsq + 1
      rhold = rhonew
      call S2ThetaW(q,g,gp,Neg,Ceg,SLegg,SNegg,ON,OCL,
     $     OCR,sinthwold, rhold, delalphhat, 
     $     alphaMZ, pizzTMZ, piwwTMW, piwwT0, mh0sq,trysinsq,flags2tw)


      endif

      RETURN
      END SUBROUTINE


!------------------------------------------------------------------------------
       RECURSIVE SUBROUTINE delrho(rhold,s2tw, alphhat,rhonew, 
     $     alphaMZ, pizzTMZ, piwwTMW, piWWT0, mh0sq,try2)

      IMPLICIT NONE
      integer try2
      DOUBLE PRECISION alphhat, pizzTMZ, piwwTMW, piWWT0
      DOUBLE PRECISION mh0sq, rhold, rhonew, delrhohat,xt
      DOUBLE PRECISION mbpole, mtaupole, Mtpole, alphaMZ 
      DOUBLE PRECISION s2tw,alph,Gf,alphas,MZpole,MW,MZ
      double precision rho2,pi,MWpole
      double precision m0, m12, m10, m20, sgnmu, tanbeta, a0

      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mwpole/ MWpole
      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
      COMMON/GAUGE/alph,Gf,alphas

!-----------
      include 'stdinputs.h'
!------------


      pi = 4.d0 * datan(1.d0)

      MW = MWpole
      MZ = MZpole

      delrhohat = 0.d0
      
      xt = 3.d0 * Gf * mtpole**2.d0 / (8.d0 * pi*pi * dsqrt(2.d0))


      delrhohat = pizzTMZ/(rhold * MZ**2.d0) - piwwTMW/(MW**2.d0) +        !--one loop sm rho
     $     alphhat * alphas/(4.d0 * pi * pi * s2tw) *(
     $     -2.145d0 * (Mtpole/MW)**2.d0 + 1.262*dlog(MTpole/MZ)  - 
     $     2.24d0 - 0.85d0 * (MZ/Mtpole)**2.d0) + 
     $     xt**2.d0/3.d0 * rho2(dsqrt(mh0sq)/mtpole) * 
     $     (dcos(alphaMZ)/ dsin(datan(tanbeta)))**2.d0

      rhonew = 1.d0/(1.d0 - delrhohat)
      
      if(try2.eq.19) WRITE(198,*) m0, m12, m10, m20, tanbeta, a0, 1

      if(try2.lt.10)then

      if(((1.d0 - (min(rhold,rhonew)/max(rhold,rhonew))).le.2.d-6))then
      
      rhold = rhonew

      RETURN
      else
         rhold = rhonew
         try2 = try2 + 1

      call delrho(rhold, s2tw, alphhat,rhonew, 
     $     alphaMZ, pizzTMZ, piwwTMW, piWWT0, mh0sq,try2)

         ENDIF

         else
            rhold = 1.01d0
            return
         endif

      RETURN
      END SUBROUTINE


!-----------------------------------------------------------------
      DOUBLE PRECISION Function rho2(ar)

      DOUBLE PRECISION pi, ar, rho2ar

      pi = 4.d0 * datan(1.d0)

      if(ar.le.1.9d0) then
         rho2ar = 19.d0 - (33.d0/2.d0)*ar + (43.d0/12.d0)*(ar*ar) +
     $        (7.d0/120.d0)*(ar*ar*ar) - pi*dsqrt(ar)*(4.d0 - 1.5d0*ar 
     $        +(3.d0/32.d0)*(ar*ar) + (1.d0/256.d0)*(ar*ar*ar)) -
     $        pi*pi*(2.d0 -2.d0*ar + 0.5d0*(ar*ar)) - 
     $        dlog(ar)*(3.d0*ar - 0.5d0*(ar*ar))
         else
            
            rho2ar = (dlog(ar))**2.d0 * (1.5d0 - (9.d0/ar) - 
     $           (15.d0/ar**2.d0) - (48.d0/ar**3.d0) - (168.d0/ar**4.d0)
     $           - (612.d0/ar**5.d0)) - 
     $           dlog(ar) * (27.d0/2.d0 + (4.d0/ar) - 
     $           (125.d0/(4.d0*ar**2.d0)) - (558.d0/(5.d0*ar**3.d0)) -
     $           (8307.d0/(20.d0*ar**4.d0)) - (109321/(70*ar**5.d0))) +
     $           pi*pi * (1.d0 - (4.d0/ar) - (5.d0/(ar**2.d0)) - 
     $           (16.d0/(ar**3.d0)) -
     $           (56.d0/(ar**4.d0)) - (204.d0/(ar**5.d0))) +
     $           (49.d0/4.d0 + (2.d0/(3.d0*ar)) - 
     $        (1613.d0/(48.d0*ar**2.d0)) -(8757.d0/(100.d0*ar**3.d0)) -
     $           (341959.d0/(1200.d0*ar**4.d0)) - 
     $           (9737663.d0/(9800*ar**5.d0)))
               
            ENDIF
            
            rho2 = rho2ar
      RETURN 
      END

!-----------------------------------------------------------

      DOUBLE PRECISION FUNCTION DELTAVB1(alphhat,rhohat,sinsqtw,
     $     sinsqtwp)
      
      implicit none
      double precision alphhat,rhohat,sinsqtw,sinsqtwp
      double precision deltavbsm,pi
      
      pi = 4.d0 * datan(1.d0)
      

      deltavbsm = ((rhohat*alphhat)/(4.d0 * pi * sinsqtw)) * 
     $     (6.d0 + (dlog((1.d0 - sinsqtwp))/sinsqtwp*(3.5d0 - 2.5d0 * 
     $     sinsqtwp -
     $     sinsqtw * ( 5.d0 - 1.5d0 * 
     $     ((1.d0 - sinsqtwp)/(1.d0 - sinsqtw))))))

      deltavb1 = deltavbsm

      RETURN
      END FUNCTION DELTAVB1
!===========================================================================


      SUBROUTINE DELTAZ(q,g,gp,Neg,Ceg,SLegg,SNegg,ON,OCL,OCR,
     $     deltaztot)
      
      implicit none
      
      integer i
      double precision q,g,gp,sinsqthw,pi
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      double precision deltaztot,deltaznue,deltaznumu,deltaze,deltazmu
      double precision Neg(4),Ceg(2),SLegg(6),SNegg(3)
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),nmneut(4),mchar(2),msnu(3)

      double precision ON(4,4),OCL(2,2),OCR(2,2)
      double precision bcharel(2),bneutnue(4),bcharmul(2),bneutnumu(4)
      double precision acharnue(2),bneutel(4),acharnumu(2),bneutmul(4)

      data bcharel/ 2 * 0.d0/,bneutnue/ 4 * 0.d0/,bcharmul/ 2 * 0.d0/,
     $     bneutnumu/ 4 * 0.d0/
      data acharnue/ 2 * 0.d0/,bneutel/ 4 * 0.d0/,acharnumu/ 2 * 0.d0/,
     $     bneutmul/ 4 * 0.d0/

      double precision apsicnue(2),apsicnumu(2),bpsi0el(4),bpsi0mul(4)
      double precision bpsicel(2),bpsicmul(2),bpsi0nue(4),bpsi0numu(4)

      double precision b1mcharmel(2),b1mcharmmul(2),b1mneutmnue(4),
     $     b1mneutmnumu(4)
      double precision b1mcharmnue(2),b1mcharmnumu(2),b1mneutmel(4),
     $     b1mneutmmul(4)

       DOUBLE PRECISION  sinsqthw_mz
       common/sinsq_mz/sinsqthw_mz


!---------------------

      external b0,b1
      include 'stdinputs.h'

!---------------------

      pi = 4.d0 * datan(1.d0)

c$$$      meR   = dsqrt(SLegg(6))
c$$$      meL   = dsqrt(SLegg(5))
c$$$      mmuR  = dsqrt(SLegg(4))
c$$$      mmuL  = dsqrt(SLegg(3))
c$$$      mtauR = dsqrt(SLegg(2))
c$$$      mtauL = dsqrt(SLegg(1))


      meR   = dsqrt(SLegg(5))
      meL   = dsqrt(SLegg(6))
      mmuR  = dsqrt(SLegg(3))
      mmuL  = dsqrt(SLegg(4))
      mtauR = dsqrt(SLegg(1))
      mtauL = dsqrt(SLegg(2))
      
      msnu(1) = dsqrt(SNegg(1))
      msnu(2) = dsqrt(SNegg(2))
      msnu(3) = dsqrt(SNegg(3))

      mneut(1) = Neg(1) 
      mneut(2) = Neg(2)
      mneut(3) = Neg(3)
      mneut(4) = Neg(4)


      nmneut(1) = dabs(Neg(1)) 
      nmneut(2) = dabs(Neg(2))
      nmneut(3) = dabs(Neg(3))
      nmneut(4) = dabs(Neg(4))

      mchar(1) = Ceg(1)
      mchar(2) = Ceg(2)

!------------------------

      sinsqthw = sinsqthw_mz

      gnuL = 0.5d0
      guL = 0.5d0 - 2.d0*sinsqthw/ 3.d0 
      gdL = -0.5d0 + sinsqthw/3.d0
      geL = -0.5d0 + sinsqthw
      guR =  2.d0*sinsqthw/3.d0 
      gdR = -sinsqthw/3.d0
      geR = -sinsqthw 
      yuL = 1.d0/3.d0 
      yuR = -4.d0/3.d0 
      ydL = 1.d0/3.d0 
      ydR = 2.d0/3.d0 
      yeL = -1.d0 
      yeR = 2.d0 
      ynuL = -1.d0       

!--------------------

      bpsicel(1) = g
      bpsicel(2) = 0.d0
      
      bcharel(1) = OCL(1,1)*bpsicel(1) + OCL(1,2)*bpsicel(2)
      bcharel(2) = OCL(2,1)*bpsicel(1) + OCL(2,2)*bpsicel(2)
      
      bpsicmul(1) = g
      bpsicmul(2) = 0.d0

      bcharmul(1) = OCL(1,1)*bpsicmul(1) + OCL(1,2)*bpsicmul(2)
      bcharmul(2) = OCL(2,1)*bpsicmul(1) + OCL(2,2)*bpsicmul(2)

      apsicnue(1) = g
      apsicnue(2) = 0.d0

      acharnue(1) = OCR(1,1)*apsicnue(1) + OCR(1,2)*apsicnue(2)
      acharnue(2) = OCR(2,1)*apsicnue(1) + OCR(2,2)*apsicnue(2)

      apsicnumu(1) = g
      apsicnumu(2) = 0.d0

      acharnumu(1) = OCR(1,1)*apsicnumu(1) + OCR(1,2)*apsicnumu(2)
      acharnumu(2) = OCR(2,1)*apsicnumu(1) + OCR(2,2)*apsicnumu(2)

      bpsi0nue(1) = gp*ynul/dsqrt(2.d0)
      bpsi0nue(2) = dsqrt(2.d0) * g * 0.5d0
      bpsi0nue(3) = 0.d0
      bpsi0nue(4) = 0.d0

      bneutnue(1) = ON(1,1)*bpsi0nue(1) + ON(1,2)*bpsi0nue(2) + 
     $     ON(1,3)*bpsi0nue(3) + ON(1,4)*bpsi0nue(4)

      bneutnue(2) = ON(2,1)*bpsi0nue(1) + ON(2,2)*bpsi0nue(2) + 
     $     ON(2,3)*bpsi0nue(3) + ON(2,4)*bpsi0nue(4)

      bneutnue(3) = ON(3,1)*bpsi0nue(1) + ON(3,2)*bpsi0nue(2) + 
     $     ON(3,3)*bpsi0nue(3) + ON(3,4)*bpsi0nue(4)

      bneutnue(4) = ON(4,1)*bpsi0nue(1) + ON(4,2)*bpsi0nue(2) + 
     $     ON(4,3)*bpsi0nue(3) + ON(4,4)*bpsi0nue(4)

      bpsi0numu(1) = gp*ynul/dsqrt(2.d0)
      bpsi0numu(2) = dsqrt(2.d0) * g * 0.5d0
      bpsi0numu(3) = 0.d0
      bpsi0numu(4) = 0.d0

      bneutnumu(1) = ON(1,1)*bpsi0numu(1) + ON(1,2)*bpsi0numu(2) + 
     $     ON(1,3)*bpsi0numu(3) + ON(1,4)*bpsi0numu(4)

      bneutnumu(2) = ON(2,1)*bpsi0numu(1) + ON(2,2)*bpsi0numu(2) + 
     $     ON(2,3)*bpsi0numu(3) + ON(2,4)*bpsi0numu(4)

      bneutnumu(3) = ON(3,1)*bpsi0numu(1) + ON(3,2)*bpsi0numu(2) + 
     $     ON(3,3)*bpsi0numu(3) + ON(3,4)*bpsi0numu(4)

      bneutnumu(4) = ON(4,1)*bpsi0numu(1) + ON(4,2)*bpsi0numu(2) + 
     $     ON(4,3)*bpsi0numu(3) + ON(4,4)*bpsi0numu(4)
      
      bpsi0el(1) = gp*yel/dsqrt(2.d0)
      bpsi0el(2) = dsqrt(2.d0) * g * (-0.5d0)
      bpsi0el(3) = 0.d0
      bpsi0el(4) = 0.d0

      bneutel(1) = ON(1,1)*bpsi0el(1) + ON(1,2)*bpsi0el(2) + 
     $     ON(1,3)*bpsi0el(3) + ON(1,4)*bpsi0el(4) 

      bneutel(2) = ON(2,1)*bpsi0el(1) + ON(2,2)*bpsi0el(2) + 
     $     ON(2,3)*bpsi0el(3) + ON(2,4)*bpsi0el(4) 

      bneutel(3) = ON(3,1)*bpsi0el(1) + ON(3,2)*bpsi0el(2) + 
     $     ON(3,3)*bpsi0el(3) + ON(3,4)*bpsi0el(4) 

      bneutel(4) = ON(4,1)*bpsi0el(1) + ON(4,2)*bpsi0el(2) + 
     $     ON(4,3)*bpsi0el(3) + ON(4,4)*bpsi0el(4) 

      bpsi0mul(1) = gp*yel/dsqrt(2.d0)
      bpsi0mul(2) = dsqrt(2.d0) * g * (-0.5d0)
      bpsi0mul(3) = 0.d0
      bpsi0mul(4) = 0.d0

      bneutmul(1) = ON(1,1)*bpsi0mul(1) + ON(1,2)*bpsi0mul(2) + 
     $     ON(1,3)*bpsi0mul(3) + ON(1,4)*bpsi0mul(4) 

      bneutmul(2) = ON(2,1)*bpsi0mul(1) + ON(2,2)*bpsi0mul(2) + 
     $     ON(2,3)*bpsi0mul(3) + ON(2,4)*bpsi0mul(4) 

      bneutmul(3) = ON(3,1)*bpsi0mul(1) + ON(3,2)*bpsi0mul(2) + 
     $     ON(3,3)*bpsi0mul(3) + ON(3,4)*bpsi0mul(4) 

      bneutmul(4) = ON(4,1)*bpsi0mul(1) + ON(4,2)*bpsi0mul(2) + 
     $     ON(4,3)*bpsi0mul(3) + ON(4,4)*bpsi0mul(4) 

!-----------------
      do i = 1, 2

         call b1(2.d-5,mchar(i),meL,q,b1mcharmel(i))
         call b1(2.d-5,mchar(i),mmuL,q,b1mcharmmul(i))

         call b1(2.d-5,mchar(i),msnu(3),q,b1mcharmnue(i))
         call b1(2.d-5,mchar(i),msnu(2),q,b1mcharmnumu(i))

      enddo
      
      do i = 1, 4

         call b1(2.d-5,nmneut(i),msnu(3),q,b1mneutmnue(i))
         call b1(2.d-5,nmneut(i),msnu(2),q,b1mneutmnumu(i))

         call b1(2.d-5,nmneut(i),meL,q,b1mneutmel(i))
         call b1(2.d-5,nmneut(i),mmuL,q,b1mneutmmul(i))

      enddo

      deltaznue = - ( (bcharel(1)**2.d0) * b1mcharmel(1) + 
     $     (bcharel(2)**2.d0) * b1mcharmel(2) +  
     $     (bneutnue(1)**2.d0) * b1mneutmnue(1) + 
     $     (bneutnue(2)**2.d0) * b1mneutmnue(2) + 
     $     (bneutnue(3)**2.d0) * b1mneutmnue(3) + 
     $     (bneutnue(4)**2.d0) * b1mneutmnue(4)) 

      deltaznumu = - ( (bcharmul(1)**2.d0) * b1mcharmmul(1) + 
     $     (bcharmul(2)**2.d0) * b1mcharmmul(2) +  
     $     (bneutnumu(1)**2.d0) * b1mneutmnumu(1) + 
     $     (bneutnumu(2)**2.d0) * b1mneutmnumu(2) + 
     $     (bneutnumu(3)**2.d0) * b1mneutmnumu(3) + 
     $     (bneutnumu(4)**2.d0) * b1mneutmnumu(4))
 
      deltaze = - ( (acharnue(1)**2.d0) * b1mcharmnue(1) + 
     $     (acharnue(2)**2.d0) * b1mcharmnue(2) + 
     $     (bneutel(1)**2.d0) * b1mneutmel(1) +
     $     (bneutel(2)**2.d0) * b1mneutmel(2) +
     $     (bneutel(3)**2.d0) * b1mneutmel(3) +
     $     (bneutel(4)**2.d0) * b1mneutmel(4)) 

      deltazmu = - ( (acharnumu(1)**2.d0) * b1mcharmnumu(1) + 
     $     (acharnumu(2)**2.d0) * b1mcharmnumu(2) + 
     $     (bneutmul(1)**2.d0) * b1mneutmmul(1) +
     $     (bneutmul(2)**2.d0) * b1mneutmmul(2) +
     $     (bneutmul(3)**2.d0) * b1mneutmmul(3) +
     $     (bneutmul(4)**2.d0) * b1mneutmmul(4))

      deltaztot = (deltaznue + deltaznumu + deltaze + deltazmu)/
     $     (32.d0 * pi*pi)


      return
      end subroutine deltaz
!=============================================================================

      subroutine deltav(s2tw, alphhat, q,g,gp,Neg,Ceg,SLegg,
     $     SNegg,ON,OCL,OCR,deltavtot)
      implicit none
      
      integer i,j,k
      double precision q,g,gp,sinsqthw, s2tw, alphhat,a1
      DOUBLE PRECISION gnuL,guL,gdL,geL,guR,gdR,geR
      DOUBLE PRECISION yuL,ydL,ydR,yeL,yeR,ynuL,yuR
      double precision deltavtot,deltave,deltavmu
      double precision Neg(4),Ceg(2),SLegg(6),SNegg(3)
      DOUBLE PRECISION meL,meR,mmuL,mmuR,mtauL,mtauR
      DOUBLE PRECISION mneut(4),nmneut(4),mchar(2),msnu(3)

      double precision ON(4,4),OCL(2,2),OCR(2,2)
      double precision bcharel(2),bneutnue(4),bcharmul(2),bneutnumu(4)
      double precision acharnue(2),bneutel(4),acharnumu(2),bneutmul(4)

      data bcharel/ 2 * 0.d0/,bneutnue/ 4 * 0.d0/,bcharmul/ 2 * 0.d0/,
     $     bneutnumu/ 4 * 0.d0/
      data acharnue/ 2 * 0.d0/,bneutel/ 4 * 0.d0/,acharnumu/ 2 * 0.d0/,
     $     bneutmul/ 4 * 0.d0/

      double precision apsicnue(2),apsicnumu(2),bpsi0el(4),bpsi0mul(4)
      double precision bpsicel(2),bpsicmul(2),bpsi0nue(4),bpsi0numu(4)

      double precision apsi0psicw(4,2),bpsi0psicw(4,2),achi(4,2)
      double precision achinotpw(4,2),bchinotpw(4,2),bchi(4,2)
      
      data apsi0psicw/ 8 * 0.d0/,bpsi0psicw/ 8 * 0.d0/,achi/ 8 * 0.d0/
      data achinotpw/ 8 * 0.d0/,bchinotpw/ 8 * 0.d0/,bchi/ 8 * 0.d0/

      double precision ONdag(4,4),OCLdag(2,2),OCRdag(2,2)
      double precision b0mcharmneut(2,4),b0melmsnue,b0mmulmsnumu

      data b0mcharmneut/ 8 * 0.d0/,b0melmsnue/ 0.d0/,b0mmulmsnumu/ 0.d0/

      double precision term1e,term2e,term3e,term1mu,term2mu,term3mu

      DOUBLE PRECISION sinsqthw_mz,mbpole,mtaupole,Mtpole,MZpole,pi,MZ
      common/sinsq_mz/sinsqthw_mz
      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole


      DOUBLE PRECISION C0, D0, D27
!---------------------

      external b0,dag2d

!--------------------

      include 'stdinputs.h'

!---------------------

      pi = 4.d0 * datan(1.d0)

      MZ = MZpole

c$$$      meR   = dsqrt(SLegg(6))
c$$$      meL   = dsqrt(SLegg(5))
c$$$      mmuR  = dsqrt(SLegg(4))
c$$$      mmuL  = dsqrt(SLegg(3))
c$$$      mtauR = dsqrt(SLegg(2))
c$$$      mtauL = dsqrt(SLegg(1))

      meR   = dsqrt(SLegg(5))
      meL   = dsqrt(SLegg(6))
      mmuR  = dsqrt(SLegg(3))
      mmuL  = dsqrt(SLegg(4))
      mtauR = dsqrt(SLegg(1))
      mtauL = dsqrt(SLegg(2))

      
      msnu(1) = dsqrt(SNegg(1))
      msnu(2) = dsqrt(SNegg(2))
      msnu(3) = dsqrt(SNegg(3))

      mneut(1) = Neg(1) 
      mneut(2) = Neg(2)
      mneut(3) = Neg(3)
      mneut(4) = Neg(4)


      nmneut(1) = dabs(Neg(1)) 
      nmneut(2) = dabs(Neg(2))
      nmneut(3) = dabs(Neg(3))
      nmneut(4) = dabs(Neg(4))

      mchar(1) = Ceg(1)
      mchar(2) = Ceg(2)

!------------------------

      sinsqthw = sinsqthw_mz

      gnuL = 0.5d0
      guL = 0.5d0 - 2.d0*sinsqthw/ 3.d0 
      gdL = -0.5d0 + sinsqthw/3.d0
      geL = -0.5d0 + sinsqthw
      guR =  2.d0*sinsqthw/3.d0 
      gdR = -sinsqthw/3.d0
      geR = -sinsqthw 
      yuL = 1.d0/3.d0 
      yuR = -4.d0/3.d0 
      ydL = 1.d0/3.d0 
      ydR = 2.d0/3.d0 
      yeL = -1.d0 
      yeR = 2.d0 
      ynuL = -1.d0       

!--------------------

      bpsicel(1) = g
      bpsicel(2) = 0.d0
      
      bcharel(1) = OCL(1,1)*bpsicel(1) + OCL(1,2)*bpsicel(2)
      bcharel(2) = OCL(2,1)*bpsicel(1) + OCL(2,2)*bpsicel(2)
      

      bpsicmul(1) = g
      bpsicmul(2) = 0.d0

      bcharmul(1) = OCL(1,1)*bpsicmul(1) + OCL(1,2)*bpsicmul(2)
      bcharmul(2) = OCL(2,1)*bpsicmul(1) + OCL(2,2)*bpsicmul(2)


      apsicnue(1) = g
      apsicnue(2) = 0.d0

      acharnue(1) = OCR(1,1)*apsicnue(1) + OCR(1,2)*apsicnue(2)
      acharnue(2) = OCR(2,1)*apsicnue(1) + OCR(2,2)*apsicnue(2)


      apsicnumu(1) = g
      apsicnumu(2) = 0.d0

      acharnumu(1) = OCR(1,1)*apsicnumu(1) + OCR(1,2)*apsicnumu(2)
      acharnumu(2) = OCR(2,1)*apsicnumu(1) + OCR(2,2)*apsicnumu(2)



      bpsi0nue(1) = gp*ynul/dsqrt(2.d0)
      bpsi0nue(2) = dsqrt(2.d0) * g * 0.5d0
      bpsi0nue(3) = 0.d0
      bpsi0nue(4) = 0.d0

      bneutnue(1) = ON(1,1)*bpsi0nue(1) + ON(1,2)*bpsi0nue(2) + 
     $     ON(1,3)*bpsi0nue(3) + ON(1,4)*bpsi0nue(4)

      bneutnue(2) = ON(2,1)*bpsi0nue(1) + ON(2,2)*bpsi0nue(2) + 
     $     ON(2,3)*bpsi0nue(3) + ON(2,4)*bpsi0nue(4)

      bneutnue(3) = ON(3,1)*bpsi0nue(1) + ON(3,2)*bpsi0nue(2) + 
     $     ON(3,3)*bpsi0nue(3) + ON(3,4)*bpsi0nue(4)

      bneutnue(4) = ON(4,1)*bpsi0nue(1) + ON(4,2)*bpsi0nue(2) + 
     $     ON(4,3)*bpsi0nue(3) + ON(4,4)*bpsi0nue(4)


      bpsi0numu(1) = gp*ynul/dsqrt(2.d0)
      bpsi0numu(2) = dsqrt(2.d0) * g * 0.5d0
      bpsi0numu(3) = 0.d0
      bpsi0numu(4) = 0.d0

      bneutnumu(1) = ON(1,1)*bpsi0numu(1) + ON(1,2)*bpsi0numu(2) + 
     $     ON(1,3)*bpsi0numu(3) + ON(1,4)*bpsi0numu(4)

      bneutnumu(2) = ON(2,1)*bpsi0numu(1) + ON(2,2)*bpsi0numu(2) + 
     $     ON(2,3)*bpsi0numu(3) + ON(2,4)*bpsi0numu(4)

      bneutnumu(3) = ON(3,1)*bpsi0numu(1) + ON(3,2)*bpsi0numu(2) + 
     $     ON(3,3)*bpsi0numu(3) + ON(3,4)*bpsi0numu(4)

      bneutnumu(4) = ON(4,1)*bpsi0numu(1) + ON(4,2)*bpsi0numu(2) + 
     $     ON(4,3)*bpsi0numu(3) + ON(4,4)*bpsi0numu(4)

      
      bpsi0el(1) = gp*yel/dsqrt(2.d0)
      bpsi0el(2) = dsqrt(2.d0) * g * (-0.5d0)
      bpsi0el(3) = 0.d0
      bpsi0el(4) = 0.d0

      bneutel(1) = ON(1,1)*bpsi0el(1) + ON(1,2)*bpsi0el(2) + 
     $     ON(1,3)*bpsi0el(3) + ON(1,4)*bpsi0el(4) 

      bneutel(2) = ON(2,1)*bpsi0el(1) + ON(2,2)*bpsi0el(2) + 
     $     ON(2,3)*bpsi0el(3) + ON(2,4)*bpsi0el(4) 

      bneutel(3) = ON(3,1)*bpsi0el(1) + ON(3,2)*bpsi0el(2) + 
     $     ON(3,3)*bpsi0el(3) + ON(3,4)*bpsi0el(4) 

      bneutel(4) = ON(4,1)*bpsi0el(1) + ON(4,2)*bpsi0el(2) + 
     $     ON(4,3)*bpsi0el(3) + ON(4,4)*bpsi0el(4) 


      bpsi0mul(1) = gp*yel/dsqrt(2.d0)
      bpsi0mul(2) = dsqrt(2.d0) * g * (-0.5d0)
      bpsi0mul(3) = 0.d0
      bpsi0mul(4) = 0.d0

      bneutmul(1) = ON(1,1)*bpsi0mul(1) + ON(1,2)*bpsi0mul(2) + 
     $     ON(1,3)*bpsi0mul(3) + ON(1,4)*bpsi0mul(4) 

      bneutmul(2) = ON(2,1)*bpsi0mul(1) + ON(2,2)*bpsi0mul(2) + 
     $     ON(2,3)*bpsi0mul(3) + ON(2,4)*bpsi0mul(4) 

      bneutmul(3) = ON(3,1)*bpsi0mul(1) + ON(3,2)*bpsi0mul(2) + 
     $     ON(3,3)*bpsi0mul(3) + ON(3,4)*bpsi0mul(4) 

      bneutmul(4) = ON(4,1)*bpsi0mul(1) + ON(4,2)*bpsi0mul(2) + 
     $     ON(4,3)*bpsi0mul(3) + ON(4,4)*bpsi0mul(4) 


!--------------------------

      loopsii: DO i = 1, 4
      loopsij: DO j = 1, 2

      apsi0psicw(i,j) = 0.d0
      bpsi0psicw(i,j) = 0.d0
      
      achinotpw(i,j) = 0.d0
      bchinotpw(i,j) = 0.d0

      achi(i,j) = 0.d0
      bchi(i,j) = 0.d0

      ENDDO loopsij
      ENDDO loopsii


      apsi0psicw(2,1) = - 1.d0 * g
      apsi0psicw(4,2) =  g/(dsqrt(2.d0))

      bpsi0psicw(2,1) = - 1.d0 * g
      bpsi0psicw(3,2) = - 1.d0 * (g/(dsqrt(2.d0)))


      call dag2d(ON,ONdag)
      call dag2d(OCL,OCLdag)
      call dag2d(OCR,OCRdag)

!--------------------------------------------------------------

      loop4by2i: DO i = 1, 4
      loop4by2j: DO j = 1, 2

      achi(i,j) = 0.d0
      bchi(i,j) = 0.d0

      loop4by2k: DO k = 1, 4
      
      
      achi(i,j) = achi(i,j) + ON(i,k) * apsi0psicw(k,j)
      bchi(i,j) = bchi(i,j) + ON(i,k) * bpsi0psicw(k,j)
      
      
      ENDDO loop4by2k
      ENDDO loop4by2j
      ENDDO loop4by2i

!---------------------------------------------------------------

      loop2i: DO i = 1, 4
      loop2j: DO j = 1, 2

      achinotpw(i,j) = 0.d0
      bchinotpw(i,j) = 0.d0

      loop2k: DO k = 1, 2

      achinotpw(i,j) = achinotpw(i,j) + achi(i,k) * OCRdag(k,j)
      bchinotpw(i,j) = bchinotpw(i,j) + bchi(i,k) * OCLdag(k,j)
      
      ENDDO loop2k
      ENDDO loop2j
      ENDDO loop2i

!-----------------

      do i = 1, 2
         do j = 1, 4
            call b0(2.d-5,mchar(i),nmneut(j),q,b0mcharmneut(i,j))
         enddo
      enddo
      
      call b0(2.d-5,meL,msnu(3),q,b0melmsnue)
      call b0(2.d-5,mmuL,msnu(2),q,b0mmulmsnumu)

!----------------


      term1e = 0.d0
      term2e = 0.d0
      term3e = 0.d0

      term1mu = 0.d0
      term2mu = 0.d0
      term3mu = 0.d0

      loopti: do i = 1, 2
      looptj: do j = 1, 4

      term1e = term1e + bcharel(i) * bneutel(j) * (( - (dsqrt(2.d0)/g) * 
     $     achinotpw(j,i) * mneut(j)*mchar(i) * 
     $     c0(meL,mchar(i),nmneut(j))) + 
     $     (bchinotpw(j,i)/(dsqrt(2.d0) * g)) * (b0mcharmneut(i,j) + 
     $     meL*meL * c0(meL,mchar(i),nmneut(j)) - 0.5d0))
      
      term1mu = term1mu + 
     $     bcharmul(i) * bneutmul(j) * (( - (dsqrt(2.d0)/g) * 
     $     achinotpw(j,i) * mneut(j)*mchar(i) * 
     $     c0(mmuL,mchar(i),nmneut(j))) + 
     $     (bchinotpw(j,i)/(dsqrt(2.d0) * g)) * (b0mcharmneut(i,j) + 
     $     mmuL*mmuL * c0(mmuL,mchar(i),nmneut(j)) - 0.5d0))
      
      term2e = term2e + 
     $     acharnue(i) * bneutnue(j) * (( - (dsqrt(2.d0)/g) * 
     $     bchinotpw(j,i) * mneut(j)*mchar(i) * 
     $     c0(msnu(3),mchar(i),nmneut(j))) + 
     $     (achinotpw(j,i)/(dsqrt(2.d0) * g)) * (b0mcharmneut(i,j) + 
     $     msnu(3)**2.d0 * c0(msnu(3),mchar(i),nmneut(j)) - 0.5d0))
      
      term2mu = term2mu + 
     $     acharnumu(i) * bneutnumu(j) * (( - (dsqrt(2.d0)/g) * 
     $     bchinotpw(j,i) * mneut(j)*mchar(i) * 
     $     c0(msnu(2),mchar(i),nmneut(j))) + 
     $     (achinotpw(j,i)/(dsqrt(2.d0) * g)) * (b0mcharmneut(i,j) + 
     $     msnu(2)**2.d0 * c0(msnu(2),mchar(i),nmneut(j)) - 0.5d0))

      enddo looptj
      enddo loopti
      
      
      loopt3j: do j = 1, 4

      term3e = term3e + bneutel(j) * bneutnue(j) * (b0melmsnue + 
     $     mneut(j)**2.d0 * c0(nmneut(j),meL,msnu(3)) + 0.5d0)

      term3mu = term3mu + bneutmul(j) * bneutnumu(j) * (b0mmulmsnumu + 
     $     mneut(j)**2.d0 * c0(nmneut(j),mmuL,msnu(2)) + 0.5d0)
      
      enddo loopt3j

      deltave = (term1e - term2e + 0.5d0 * term3e)/(16.d0 * pi*pi)

      deltavmu = (term1mu - term2mu + 0.5d0 * term3mu)/(16.d0 * pi*pi)
!----------------

      a1 = 0.d0

      loopa1i: do i = 1, 2
      loopa1j:  do j = 1, 4

      a1 =    a1 +  0.5d0 * (acharnumu(i) * bcharel(i) * bneutnumu(j) * 
     $     bneutel(j) * mchar(i) * mneut(j) * 
     $     D0(meL,msnu(2),mchar(i), nmneut(j)) ) +
     $     0.5d0 * (acharnue(i) * bcharmul(i) * bneutnue(j) * 
     $     bneutmul(j) * mchar(i) * mneut(j) * 
     $     D0(mmuL,msnu(3), mchar(i), nmneut(j))) + 
     $     (bcharmul(i) * bcharel(i) * bneutmul(j) * bneutel(j)
     $     * D27(mmuL,meL,mchar(i), nmneut(j)) ) +
     $     (acharnumu(i) * acharnue(i) * bneutnumu(j) * bneutnue(j) * 
     $     D27(msnu(2),msnu(1),mchar(i), nmneut(j)))


      enddo loopa1j
      enddo loopa1i


      deltavtot = deltave + deltavmu -  
     $     s2tw * (1.d0 - s2tw) *MZ*MZ/alphhat * a1/(32.d0*pi*pi*pi )

      return

      end subroutine deltav

!-----------------------------------------------------------------------------------------------
