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
!-----------------------------------------------------------------------

      RECURSIVE SUBROUTINE S2ThetaWSM(q,g,gp,Neg,Ceg,SLegg,SNegg,ON,OCL,
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

!      delr = delr1 + delr2 + deltavtot + deltaztot + 
!     $     DELTAVB1(alphhat,rhonew,sinthwold**2.d0,
!     $     (1.d0 - (MW*MW)/(MZ*MZ)))

      delr = delr1 + delr2 + 
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
      call S2ThetaWSM(q,g,gp,Neg,Ceg,SLegg,SNegg,ON,OCL,
     $     OCR,sinthwold, rhold, delalphhat, 
     $     alphaMZ, pizzTMZ, piwwTMW, piwwT0, mh0sq,trysinsq,flags2tw)


      endif

      RETURN
      END SUBROUTINE

