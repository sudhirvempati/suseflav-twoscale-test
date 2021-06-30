C=========================================================================================================
* Program : SuSeFLAV: Supersymmetry seesaw and flavor violation calculator
* Version : 1.0
* Website : http://cts.iisc.ernet.in/Suseflav/main.html
* Authors : Debtosh Choudhury  debtosh@cts.iisc.ernet.in 
*           Raghuveer Garani   rgarani@cts.iisc.ernet.in
*           Sudhir Vempati     vempati@cts.iisc.ernet.in
C--------------------------------------------------------------------------
C     iteration for minimization of \mu parameter
C-------------------------------------------------------------------------
****f* SuSeFLAV/ewsbiterate.f/RECURSIVE SUBROUTINE ITERATE 
*  NAME
*    Subroutine Iterate
*  SYNOPSIS
*    Recursive subroutine for \mu parameter minimization 
*  FUNCTION
*    This suborutine minimizes \mu parameter using iterative methods 
*    and incorporates REWSB at msusy. Also, it checks for D-flat conditions and 
*    whether the determinant is negative. 
*  INPUTS
*    scale       - Scale at which minimization is done(ex. msusy or MZ)
*    mt          - One-loop corrected mass of top quark(GeV).
*    murge       - \mu at tree level
*    msusy       - scale at which susy is broken.
*    newtbeta    - the ratio of the vevs of the two–Higgs doublet fields at msusy.
*
*  RESULT
*    bmur        - b_{\mu} obtained from converged value of \mu
*    muflag      - flags unphysical \mu. muflag = 2 for \mu^2 <0 and if 
*                  iteration count exceeds the limit.
*    try1        - Counter for iteration.
*    muold       - Converged/Minimized value of |\mu|.
*    itcount     - Counter for RGE iteration.
*    flags       - Saves the error encountered in a character string.
*    mursq       - value of \mu^2
*
*  EXAMPLE
*      RECURSIVE SUBROUTINE ITERATE (scale,mt,murge,bmur,newtbeta,
*     $     msusy,muflag,try1,muold,itcount,flags, mursq)
*
*  NOTES
*     Common Blocks used:
*      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
*      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
*     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut
*      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
*     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG
*      common/sparticles_susy/SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
*     $     mh0sq,mhu0sq,mhpmsq,mA0sq
*      common/sminputs/ mbpole, mtaupole, Mtpole
*      common/sinsq_susy/sinsqthw_susy
*      common/hzV/ delta1,delta2
*
*     External routines used in this routine:
*     EXTERNAL pizz,tadpole1,piww
*     EXTERNAL tadpole2
*
*  BUGS
*    ---
*  SEE ALSO
*    -----
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      RECURSIVE SUBROUTINE ITERATE (scale,mt,murge,bmur,newtbeta,
     $     msusy,muflag,try1,muold,exitcalc,flags,mursq)

      IMPLICIT NONE 
      
      INTEGER try1,muflag 
      character(len=1):: exitcalc
      CHARACTER*100 flags
      DOUBLE PRECISION muold,munew,mt,MW,g3
      DOUBLE PRECISION tol,ratio,pi,q,pizzT,MZ,p,piwwT
      DOUBLE PRECISION alph2,alph3,alph1,CMusq
      DOUBLE PRECISION modmu,delta1,delta2,tan2beta,tanbeta,sgnmu
      double precision gp,sinbeta,cosbeta

      double precision correction

      data correction/0.d0/         

      double precision M3t
      DOUBLE PRECISION scale

      double precision bmur,newtbeta,msusy,mh1mz,mh2mz
      double precision mSQRG(3,3),mSURG(3,3),mSDRG(3,3),mSLRG(3,3),
     $     mSERG(3,3), MChar(2,2)
      double precision yuRG(3,3),ydRG(3,3),yeRG(3,3)
      double precision ADRG(3,3),AURG(3,3),AERG(3,3)
      double precision SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      double precision Neg(4),Ceg(2)
      double precision mh0sq,mhu0sq,mhpmsq,mA0sq
      double precision MNeut(4,4),ON(4,4)
      double precision OCR(2,2),OCL(2,2)
      double precision M2tz,M3tz
      double precision vev1,vev2,cos2beta,g
      double precision m0,m12,a0
      double precision m10,m20,mur

      DOUBLE PRECISION  ANURG(3,3), mSNURG(3,3)
      DOUBLE PRECISION beta,  M1tz
      DOUBLE PRECISION  murge, mursq
      DOUBLE PRECISION bmurcor, mAm3
      double precision sin2beta
     
      DOUBLE PRECISION mt_r,mb_r,mtau_r              !<---define common top etc mass @ msusy

      double precision noDflat, detisneg,Sueggit(6),Sdeggit(6)
      DOUBLE PRECISION SLeggit(6),SNeggit(6),Negit(4),Cegit(2)
      DOUBLE PRECISION mh0sqit,mhu0sqit,mhpmsqit,mA0sqit
      DOUBLE PRECISION  mbpole, mtaupole, Mtpole, sinsqthw_susy, musave

      DOUBLE PRECISION muwithsign

      double precision MWpole, MZpole

      double precision newvev,vev1n,vev2n

      double precision thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae

      double precision mtsq_r,mgsq,mst1sq,mst2sq,st,ct,vv,S1top,S2top
      double precision mbsq_r,msb1sq,msb2sq,sb,cb,S1btm,S2btm,
     $     S1stop,S2stop,qsq,cotbeta,mtausq_r,msnusq,mstau1sq,
     $     mstau2sq,stau,ctau,S1tau,S2tau
      double precision CMusq2loop,delta12loop,delta22loop

      double precision MZrun,MWrun
!----------------------------------
      common/muiterate/ SUeggit,SDeggit,SLeggit,SNeggit,Negit,Cegit,
     $     mh0sqit,mhu0sqit,mhpmsqit,mA0sqit, muwithsign

      common/mustuff/musave
      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0

      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut

      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG

      common/sparticles_susy/SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq

      common/sminputs/ mbpole, mtaupole,Mtpole,MZpole
      common/mwpole/MWpole
      common/sinsq_susy/sinsqthw_susy

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae
     
      common/hzV/ delta1,delta2

      common/gbrunning/ MZrun,MWrun


!-----------------------------------------------------------------------
      EXTERNAL pizz,tadpole1,piww
      external tadpole2

!-----------------------------------------------------------------------

      msusy = scale


      beta = datan(newtbeta)
      tan2beta = dtan(2.d0*datan(newtbeta))
      sinbeta  = dsin(datan(newtbeta))
      cosbeta  = dcos(datan(newtbeta))
      cos2beta = dcos(2.d0*datan(newtbeta))
      sin2beta = 2.d0*sinbeta*cosbeta

     
      pi = 4.d0*datan(1.d0)
      q = scale

      delta1 = 0.d0
      delta2 = 0.d0
      pizzT  = 0.d0
      piwwT  = 0.d0
 
      try1 = try1 + 1


      Sueggit(1) = dabs(Suegg(1))
      Sueggit(2) = dabs(Suegg(2))
      Sueggit(3) = dabs(Suegg(3))
      Sueggit(4) = dabs(Suegg(4))
      Sueggit(5) = dabs(Suegg(5))
      Sueggit(6) = dabs(Suegg(6))

      SDeggit(1) = dabs(SDegg(1))
      SDeggit(2) = dabs(SDegg(2))
      SDeggit(3) = dabs(SDegg(3))
      SDeggit(4) = dabs(SDegg(4))
      SDeggit(5) = dabs(SDegg(5))
      SDeggit(6) = dabs(SDegg(6))

      SLeggit(1) = dabs(SLegg(1))
      SLeggit(2) = dabs(SLegg(2))
      SLeggit(3) = dabs(SLegg(3))
      SLeggit(4) = dabs(SLegg(4))
      SLeggit(5) = dabs(SLegg(5))
      SLeggit(6) = dabs(SLegg(6))

      SNeggit(1) = dabs(SNegg(1)) 
      SNeggit(2) = dabs(SNegg(2)) 
      SNeggit(3) = dabs(SNegg(3))

      negit(1) = (Neg(1))
      negit(2) = (Neg(2))
      negit(3) = (Neg(3))
      negit(4) = (Neg(4))
      
      Cegit(1) = (Ceg(1))
      Cegit(2) = (Ceg(2))



      if(try1.eq.1)then
         muold = murge
         M3t = M3tz
      endif

      mA0sqit = dabs(mA0sq) 
      mh0sqit = dabs(mh0sq)
      mhpmsqit = dabs(mhpmsq)
      mHu0sqit = dabs(mHu0sq)

!-------------------------------------------------------
! 10          FORMAT(f8.2,1x, f8.2,1x, f8.2,1x, f8.2,1x, f8.2,1x, f8.2,
!     $     1x,f8.2,1x,I1,1x,I2)


      if(try1.ge.10)then
         if(isnan(muold).or.isnan(bmur))then
            muflag = 2
            flags = ' MUNOC'
            return
         else            
            muold = 10.d0
            bmur = 10.d0
            muflag = 2
            flags = ' MUNOC'
            return
         endif
      endif
         
 97   format(1pe11.4,2x,1pe11.4,2x,1pe11.4,2x,1pe11.4)
         
         
      if(exitcalc.eq.'T')then
          
         if(muwithsign.lt.0.d0)then
            muflag = 2
            flags = ' REWSB'    !-Invalid: mursq negative'
!            write(98,97) m0,m12,tanbeta,a0
            return
         endif
         
         if(bmur.lt.0.d0.and.flags.ne.' REWSB')then
            flags='BMUNEG'            
            muflag = 2  
!            write(199,97) m0,m12,tanbeta,a0
            return
         endif

      endif


      if(bmur.lt.0.d0) bmur = dabs(bmur)
!------------------------------------------------------

      MZ = MZpole               !91.1876d0
      MW = MWpole               !80.398d0

      g  = dsqrt(alph2*16.d0*pi*pi)
      gp = dsqrt(alph1*16.d0*pi*pi) * dsqrt(3.d0/5.d0)
      g3 = dsqrt(alph3*16.d0*pi*pi)
 
c$$$      MW = dsqrt((alph2*4.d0*pi*pi)*(vev1**2.d0 + vev2**2.d0))
c$$$      
c$$$      MZ = dsqrt(((alph2*4.d0*pi*pi) + 
c$$$     $     (alph1*4.d0*pi*pi) * (3.d0/5.d0)) *
c$$$     $     (vev1**2.d0 + vev2**2.d0))

      mt_r = yuRG(3,3)* vev2/dsqrt(2.d0)
      mb_r = ydRG(3,3)* vev1/dsqrt(2.d0)
      mtau_r = yeRG(3,3)* vev1/dsqrt(2.d0)

c$$$      print*,"mt_r = ", mt_r, "yuRG(3,3) = ", yuRG(3,3)
c$$$      print*,"vev2 = ", vev2
c$$$      print*,"vev1 = ", vev1

      sinsqthw_susy = (gp**2.d0/(g**2.d0 + gp**2.d0))

      modmu = dabs(muold)
      mAm3 = bmur/ (sinbeta*cosbeta)  
     
!-------------------------------------------------------


      p = MZ
 
      q = scale

      call pizz(p,q,g,mt_r,mb_r,mtau_r,newtbeta,SUeggit,SDeggit,
     $     SLeggit,
     $     SNeggit,Negit,Cegit,mh0sqit,mhu0sqit,mhpmsqit,mA0sqit,
     $     ON,OCL,OCR,sinsqthw_susy,pizzT)


      p = MW

      call piww(p,q,g,mt_r,mb_r,mtau_r,newtbeta,SUeggit,SDeggit,
     $     SLeggit,SNeggit,Negit,Cegit,
     $     mh0sqit,mhu0sqit,mhpmsqit,mA0sqit,ON,OCL,OCR,
     $     sinsqthw_susy,piwwT)


!------------------------------

      call vevewsb(newtbeta,MZ,pizzT,gp,g,vev1n,vev2n,newvev)

c$$$      print*,"vev2n = ", vev2n
c$$$      print*,"vev1n = ", vev1n


      vev1 = vev1n
      vev2 = vev2n

c$$$      print*,"MZ, MW = ", dsqrt(((alph2*4.d0*pi*pi) + 
c$$$     $     (alph1*4.d0*pi*pi) * (3.d0/5.d0)) *
c$$$     $     (vev1**2.d0 + vev2**2.d0)), 
c$$$     $     MZ*dcos(dasin(dsqrt(sinsqthw_susy))),
c$$$     $     dsqrt((alph2*4.d0*pi*pi)*(vev1**2.d0 + vev2**2.d0))

      MZrun = dsqrt(((alph2*4.d0*pi*pi) + 
     $     (alph1*4.d0*pi*pi) * (3.d0/5.d0)) *
     $     (vev1**2.d0 + vev2**2.d0))

      MWrun = dsqrt((alph2*4.d0*pi*pi)*(vev1**2.d0 + vev2**2.d0))

!      print*,"MWrun, MZrun = ", MWrun, MZrun

!------------------------------
      q = scale

      call tadpole1(q,g,mb_r,mtau_r,newtbeta,
     $     ADRG,AERG,yuRG,ydRG,yeRG,SUeggit,SDeggit,SLeggit,
     $     SNeggit,Negit,Cegit,ON,OCL,OCR,sgnmu,modmu,mh0sqit,mhu0sqit,
     $     mHpmsqit,mA0sqit,vev1,delta2)

      
      q = scale
      
      call tadpole2(q,g,mt_r,newtbeta,
     $     AURG,yuRG,ydRG,yeRG,SUeggit,SDeggit,SLeggit,
     $     SNeggit,Negit,Cegit,ON,OCL,OCR,sgnmu,modmu,mh0sqit,mhu0sqit,
     $     mHpmsqit,mA0sqit,vev2,delta1)
      

!------------------------------------------------------------------
!     2 loop higgs tadpole O(a_t a_s) correction from Slavich et al.
!------------------------------------------------------------------
      
      mtsq_r = mt_r*mt_r
      mbsq_r = mb_r*mb_r
      mtausq_r = mtau_r*mtau_r

      mgsq = M3tz*M3tz

      mst1sq = Sueggit(2)
      mst2sq = Sueggit(1)

      msb1sq = Sdeggit(2)
      msb2sq = Sdeggit(1)

      msnusq = SNeggit(1)

      mstau1sq = SLeggit(2)
      mstau2sq = SLeggit(1)

      st = dsin(thetat)
      ct = dcos(thetat)
      sb = dsin(thetab)
      cb = dcos(thetab)
      stau = dsin(thetatau)
      ctau = dcos(thetatau)

      qsq = scale*scale
      vv = (vev1**2.d0 + vev2**2.d0)
      cotbeta = 1.d0/newtbeta
      
      call ewsb2loop(mtsq_r,M3tz,mst1sq,mst2sq,st,ct,qsq,
     $     -1.d0*sgnmu*modmu,
     $     newtbeta,vv,g3,S1stop,S2stop)

      call ewsb2loop(mbsq_r,M3tz,msb1sq,msb2sq,sb,cb,qsq,
     $     -1.d0*sgnmu*modmu,
     $     cotbeta,vv,g3,S2btm,S1btm)

      call DDStad(mtsq_r,mbsq_r,mA0sqit,mst1sq,mst2sq,msb1sq,msb2sq,
     $     st,ct,sb,cb,qsq,-1.d0*sgnmu*modmu,newtbeta,vv,S1top,S2top)

      call tausqtad(mtausq_r,mA0sqit,msnusq,mstau1sq,mstau2sq,stau,
     $     ctau,qsq,-1.d0*sgnmu*modmu,newtbeta,vv,S1tau,S2tau)

      delta12loop = delta1 - S2stop - S2top - S2btm - S2tau
      delta22loop = delta2 - S1stop - S1top - S1btm - S1tau
 
!--------------------


      bmurcor = (0.5d0*sin2beta*((mh2mz - delta2) + (mh1mz - delta1) +
     $     2.d0*modmu*modmu))
      
      mAm3 = bmurcor/ (sinbeta*cosbeta)  


!-------
      CMusq = - 0.5d0 * MZ*MZ - 0.5d0 * pizzT + 0.5d0 *
     $     tan2beta * (((mh1mz-delta1)*newtbeta) - 
     $     ((mh2mz-delta2)/newtbeta))


      muwithsign = CMusq


      CMusq2loop = - 0.5d0 * MZ*MZ - 0.5d0 * pizzT + 0.5d0 *
     $     tan2beta * (((mh1mz-delta12loop)*newtbeta) - 
     $     ((mh2mz-delta22loop)/newtbeta))

!      print*,"mu1loop, mu2loop = ", dsqrt((CMusq)), dsqrt((CMusq2loop))
!      print*,"mu1loop, mu2loop = ", ((CMusq)), ((CMusq2loop))


!---------
      
      if(CMusq.lt.0.d0)then
         munew = dsqrt(dabs(CMusq))
      else
         munew = dsqrt((CMusq))
      ENDIF
      

      modmu = munew
      

!---------------------------------------------------------------------
!     CHEKING DFLAT AND DETERMINANT IS NEGATIVE FLAGS
!---------------------------------------------------------------------

      noDflat = 2.d0*(Cmusq - bmurcor) + mh2mz -delta2 + mh1mz - delta1 

      detisneg = (Cmusq+mh1mz-delta1)*(Cmusq+mh2mz-delta2) - 
     $     bmurcor*bmurcor

      mur = munew

C-----------------------------------------------------------------------
C     checking for convergence
C-----------------------------------------------------------------------
      
      ratio = 1.d0 - (MIN(dabs(muold),dabs(munew))/
     $     MAX(dabs(muold),dabs(munew)))
      
!      print*,ratio, muold, munew

      tol = 1.d-4     
      
      if((ratio.le.tol).and.Cmusq.gt.0.d0)then
         
         muold = munew
         mur = muold
         try1 = 0
         munew = 0.d0
         
c$$$         print*,"mh1mz = ", mh1mz, "mh2mz = ", mh2mz
c$$$         print*,"delta1 = ", delta1, "delta2 = ", delta2

         bmurcor = (0.5d0*sin2beta*((mh2mz-delta2) + (mh1mz-delta1) +
     $     2.d0*muold*muold))

         bmur = bmurcor

!         print*,"mur, bmur = ", mur, bmur
         
         RETURN
         
      else

         muold = munew
         bmur = (bmurcor)


      endif
      

!-----------------------------------------

      call iterate(scale,mt,murge,bmur,newtbeta,
     $     msusy, muflag,try1,muold,exitcalc,flags,mursq)


      RETURN
      
      END SUBROUTINE ITERATE
!==========================================================================
!--------------------------------------------------------------------------
C     Corrections to VEV
C---------------------------------------------------------------------------
      
      SUBROUTINE vevewsb(tanbeta,MZ,pizzT,gp,g,vev1n,vev2n,
     $     newvev)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION MZ,pizzT,gp,g,vev1n,vev2n
      DOUBLE PRECISION pi,tanbeta,cosbeta,sinbeta,newvev
      
      cosbeta = dcos(datan(tanbeta))
      sinbeta = dsin(datan(tanbeta))
      
      pi = 4.d0*datan(1.d0)

      newvev = dsqrt(4.d0*(MZ*MZ + pizzT)/(gp**2.d0 + g**2.d0))
      
      vev1n = newvev*cosbeta
      vev2n = newvev*sinbeta
      
      
      RETURN
      
      END SUBROUTINE vevewsb

C==============================================================================

      Subroutine parsdata(m0,m12,flagdet,flagDflat)
      
      IMPLICIT NONE 

      INTEGER num,flagdet,flagDflat
      DOUBLE PRECISION m0,m12
      
      if((flagdet.eq.1).AND.(flagDflat.eq.1))then
         
         num = 1
         WRITE(88,*) m0,m12,num
         
      else
           continue
      ENDIF

      END SUBROUTINE parsdata

!==============================================================================

****f* SuSeFLAV/ewsbiterate.f/REWSBCOR
*  NAME
*    Subroutine REWSBCOR
*  SYNOPSIS
*    One loop correction to all MSSM paramters at msusy.
*  FUNCTION
*    Calculates complete one-loop threshold corrections to all MSSM particles.
*    We closely follow BPMZ [hep-ph/9606211].
*     
*  INPUTS
*     mur       -  Converged value of \mu parameter.
*     bmur      -  Converged calue of b_{\mu} parameter.
*     sgnmu     -  Sign of \mu parameter.
*     newtbeta  -  \tan\beta at msusy scale. 
*     MT        -  One-loop corrected mass of top quark (GeV).      
*     msusy     -  Susy breaking scale (GeV).
*  RESULT
*    Corrected sfermions, gauginos and higgs masses. 
*      M3t       - One-loop corrected gluino mass(GeV).
*      mAm3      - One-loop corrected psuedo scalar higgs mass(GeV).
*      STeg      - (1 X 2) One-loop corrected (s)top masses, contains mixed state of L and R.
*      SCeg      - (1 X 2) One-loop corrected (s)charm masses, L and R components.
*      SUqeg     - (1 X 2) One-loop corrected (s)up masses, L and R components.
*      SBeg      - (1 X 2) One-loop corrected (s)bottom masses, contains mixed state of L and R.
*      SSTeg     - (1 X 2) One-loop corrected (s)strange masses, L and R components.
*      SDneg     - (1 X 2) One-loop corrected (s)down masses, L and R components.
*      STaueg    - (1 X 2) One-loop corrected (s)tau masses, contains mixed state of L and R.
*      SMUeg     - (1 X 2) One-loop corrected (s)mu masses, L and R components.
*      SEeg      - (1 X 2) One-loop corrected (s)electron masses, L and R components.
*      tsnu      - One-loop corrected tau (s)neutrino mass.
*      musnu     - One-loop corrected mu (s)neutrino mass.
*      elsnu     - One-loop corrected electron (s)neutrino mass.
*      newmA0sq  | one loop corrected higgs masses.  
*      newmh0sq  | Currently we use compact analytical expression at two-loop level for the lightest higgs boson mass. 
*      newmhpmsq | Heinenmeyer, Hollik and Weiglein [hep-ph/ 9903404]
*      newmHu0sq |
*      Cegm      - (1 X 2) One-loop corrected chargino masses.
*      negm      - (1 X 4) One-loop corrected neutralino masses.
*  EXAMPLE
*          SUBROUTINE  REWSBCOR(mur,bmur,sgnmu,newtbeta,MT,
*     $     msusy,msnew,M3t,mAm3,
*     $     STeg,SCeg,SUqeg,SBeg,SSTeg,SDneg,STaueg,SMUeg,SEeg,
*     $     tsnu,musnu,elsnu, newmA0sq, newmh0sq, newmhpmsq,
*     $     newmHu0sq,Cegm,negm,itcount) 
*
*  NOTES
*    Common blocks used in this routine:
*
*      common/sminputs/ mbpole, mtaupole, Mtpole
*      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
*     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut
*      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
*     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG
*      common/sparticles_susy/SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
*     $     mh0sq,mhu0sq,mhpmsq,mA0sq
*      common/sinsq_susy/sinsqthw_susy
*      common/opt_mixing/OCRm, OCLTm,ONew, alpha
*      common/vev_ewsb/ vev1n,vev2n
*
*    External subroutines in this routine:
*
*      EXTERNAL pizz,tadpole1,piww
*      EXTERNAL tadpole2,piaa,vevewsb,pihphm,pistop
*      EXTERNAL pisbottom,pischarm,pisupq,pistaul
*      EXTERNAL pisdown,pistausnu,higgs_analytical
*
*  BUGS
*    ---
*  SEE ALSO
*    -----
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!=============================================================================
      
      SUBROUTINE REWSBCOR(sgnmu,mur,bmur,newtbeta,MT,
     $     msusy,msnew,M3t,mAm3,tanbeta,
     $     STeg,SCeg,SUqeg,SBeg,SSTeg,SDneg,STaueg,SMUeg,SEeg,
     $     tsnu,musnu,elsnu, newmA0sq, newmh0sq, newmhpmsq,
     $     newmHu0sq,Cegm,negm) 


      IMPLICIT NONE
      DOUBLE PRECISION mur,mt,MW,pihphmT,g3
      DOUBLE PRECISION pi,q,pizzT,MZ,p,piwwT
      DOUBLE PRECISION alph2,alph3,alph1,piaaT,vev1n,vev2n
      DOUBLE PRECISION modmu,delta1,delta2,tan2beta,newtbeta
      double precision gp,sinbeta,cosbeta
      DOUBLE PRECISION newmA0sq,newmhpmsq,newmh0sq,newmHu0sq
      DOUBLE PRECISION STeg(2),SBeg(2)
      DOUBLE PRECISION SCeg(2),SUqeg(2),STaueg(2),SEeg(2),SMueg(2)
      DOUBLE PRECISION SSTeg(2),SDneg(2),tsnu,musnu,elsnu

      DOUBLE PRECISION STegtmp(2),SBegtmp(2),SSTegtmp(2),SDnegtmp(2)
      DOUBLE PRECISION SCegtmp(2),SUqegtmp(2),STauegtmp(2),SEegtmp(2),
     $     SMuegtmp(2)

      data piwwT/ 1 * 0.d0/, pizzT/ 1 * 0.d0/

      double precision M3t,deltagluino

      double precision MNeut1(4,4),Negm(4),neutmasstot(4,4),MChar(2,2)

      double precision charmasstot(2,2),Cegm(2),msnew     
      
c$$$      double precision Neutevc(4,4), msnew
c$$$      data Neutevc/ 16 * 0.d0/

!---------------------------------------------------------------------------
      integer i,j,OS,itcount,higgscorr
      
      parameter (OS = 0)

      double precision bmur,msusy,mh1mz,mh2mz,sgnmu
      double precision mSQRG(3,3),mSURG(3,3),mSDRG(3,3),mSLRG(3,3),
     $     mSERG(3,3),mSNURG(3,3)
      double precision yuRG(3,3),ydRG(3,3),yeRG(3,3)
      double precision ADRG(3,3),AURG(3,3),AERG(3,3), ANURG(3,3)
      double precision SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      double precision Neg(4),Ceg(2)
      double precision mh0sq,mhu0sq,mhpmsq,mA0sq
      double precision MNeut(4,4),ON(4,4)
      double precision OCR(2,2),OCL(2,2)
      double precision M1tz,M2tz,M3tz
      double precision vev1,vev2,cos2beta,g,alpha
      Double precision bmurcor,mAm3,mbpole, mtaupole, Mtpole
      data bmurcor/0.d0/

      DOUBLE PRECISION sinsqthw_susy
      DOUBLE PRECISION Mhiggsa(2,2),sigphi1,sigphi12,sigphi2,sig2phi2
      DOUBLE PRECISION newmha, mt_r, mb_r, mtau_r, sig2phi2yuk,newvev
      DOUBLE PRECISION delpisbtm(2,2),delpist(2,2),delpistau(2,2),
     $     ONew(4,4)
      DOUBLE PRECISION delpisc(2,2), delpisst(2,2), OCRm(2,2),OCLTm(2,2)
      DATA delpist/ 4 * 0.d0/,delpisbtm/ 4 * 0.d0/,delpistau/ 4 * 0.d0/
      DATA delpisc/ 4 * 0.d0/, delpisst/ 4 * 0.d0/

      DOUBLE PRECISION tanbeta

      double precision mhiggsev(2),mhmix(2,2)
      
      double precision pis1s1ans,pis2s2ans,pis1s2ans

      double precision thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae
      
      double precision mtsq_r,mgsq,mst1sq,mst2sq,st,ct,vv,
     $     S11stop,S22stop,S12stop,S11btm,S22btm,S12btm,p2btm
      double precision mbsq_r,msb1sq,msb2sq,sb,cb,qsq,cotbeta,
     $     mtausq_r,msnusq,mstau1sq,mstau2sq,stau,ctau,p2stop,
     $     S11top,S22top,S12top,p2top,S11tau,S22tau,S12tau,p2tau
      double precision sigh02loop(2,2), delmA, Mhiggs2loop(2,2)

      double precision Mhiggstree(2,2), heign(2),hl(2,2),alphatree
      double precision Mhiggsn(2,2),alphacorr

      double precision MWpole, MZpole
      
      double precision mh0sqcor,mhu0sqcor,mA0sqcor,mHpmsqcor

!--------

      common/higgs2loop/ mh0sqcor,mhu0sqcor,mA0sqcor,mHpmsqcor

      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole

      common/mwpole/MWpole

      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut

      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG

      common/sparticles_susy/SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq

      common/sinsq_susy/sinsqthw_susy

      common/opt_mixing/OCRm, OCLTm,ONew, alpha

      common/vev_ewsb/ vev1n,vev2n
      common/hzV/ delta1,delta2

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,
     $     thetamu,thetau,thetad,thetae

!---------------------------------------------------------------------------

      EXTERNAL pizz,tadpole1,piww
      external tadpole2,piaa,vevewsb,pihphm,pistop
      EXTERNAL topcor,pisbottom,pischarm,bottomcor,taucor,pisupq,pistaul
      EXTERNAL pisdown,pistausnu,higgs_analytical

!---------------------------------------------------------------------------

      MZ = MZpole
      MW = MWpole


      tan2beta = dtan(2.d0*datan(newtbeta))
      sinbeta  = dsin(datan(newtbeta))
      cosbeta  = dcos(datan(newtbeta))
      cos2beta = dcos(2.d0*datan(newtbeta))

     
      pi = 4.d0 * datan(1.d0)
      q = msusy

      delta1 = 0.d0
      delta2 = 0.d0
      pizzT  = 0.d0
      piwwT = 0.d0

      g  = dsqrt(alph2*16.d0*pi*pi)
      gp = dsqrt(alph1*16.d0*pi*pi*3.d0/5.d0) 
      g3 = dsqrt(alph3*16.d0*pi*pi)
      
!      MW = dsqrt((alph2*4.d0*pi*pi)*(vev1**2.d0 + vev2**2.d0))
      
c$$$      MZ = dsqrt(((alph2*4.d0*pi*pi) + 
c$$$     $     (alph1*4.d0*pi*pi) * (3.d0/5.d0)) *
c$$$     $     (vev1**2.d0 + vev2**2.d0))


c$$$      print*,"MW, MZ in EWSBCOR = ", dsqrt((alph2*4.d0*pi*pi)*
c$$$     $     (vev1**2.d0 + vev2**2.d0)), dsqrt(((alph2*4.d0*pi*pi) + 
c$$$     $     (alph1*4.d0*pi*pi) * (3.d0/5.d0)) *
c$$$     $     (vev1**2.d0 + vev2**2.d0))

      modmu = mur

      Mt_r = yuRG(3,3)* vev2/dsqrt(2.d0)
      Mb_r = ydRG(3,3)* vev1/dsqrt(2.d0)
      Mtau_r = yeRG(3,3)* vev1/dsqrt(2.d0)


!-------------------

      mA0sq = bmur/ (sinbeta*cosbeta) !mh1mz + mh2mz + 2.d0 * modmu*modmu

!      print*,"in REWSB mA0sq, bmur, MZ, MW = ", mA0sq, bmur, MZ, MW
      
      mhpmsq = mA0sq + MW*MW

      Mhiggstree(1,1) = mA0sq * sinbeta**2.d0 + 
     $     (MZ*MZ + pizzT) * cosbeta**2.d0

!      print*,"in REWSB pizzT, piwwT = ", pizzT, piwwT 

      Mhiggstree(1,2) = -1.d0*(mA0sq + (MZ*MZ) + pizzT) * 
     $     sinbeta*cosbeta
      
      Mhiggstree(2,1) = Mhiggstree(1,2)

      Mhiggstree(2,2) = mA0sq * cosbeta**2.d0 + 
     $     (MZ*MZ + pizzT) * sinbeta**2.d0

      Mhiggsn(1,1) = Mhiggstree(1,1) 

      Mhiggsn(1,2) = Mhiggstree(1,2) 

      Mhiggsn(2,1) = Mhiggstree(2,1) 

      Mhiggsn(2,2) = Mhiggstree(2,2) 


      Call CEigensystem(2,Mhiggsn,2,heign,hl,2,-1)

      mh0sq = heign(2)       
      mhu0sq = heign(1)      
      
!      alphatree = datan(hl(1,2)/hl(1,1))


!-------------------

      p = M3tz
      q = msusy
      
      
      call gluinose(p,q,g3,mt_r,mB_r,SUegg,SDegg,deltagluino)

      
      M3t = M3tz * (1.d0 + deltagluino)
      
!      print*,"m3tz, deltag, m3t = ", m3tz, deltagluino, m3t

!------------------------------------------------------------------------

      p = MZ

      q = msusy

      call pizz(p,q,g,mt_r,mb_r,mtau_r,newtbeta,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,ON,OCL,OCR,sinsqthw_susy,
     $     pizzT)


!------------------------------------------------------------------------
      q = msusy

      p = MW

      call piww(p,q,g,mt_r,mb_r,mtau_r,newtbeta,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,ON,OCL,OCR,sinsqthw_susy,
     $     piwwT)

!-----------------------------------------------------------------------------


      call vevewsb(newtbeta,MZ,pizzT,gp,g,vev1n,vev2n,
     $     newvev)
      
!-----------------------------------------------------------------------------

      q = msusy
      p  = q
      
      call chargino(p,q,g,gp,mt_r,mb_r,mtau_r,newtbeta,yuRG,ydRG,yeRG,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,MChar,ON,OCL,OCR,mh0sq,
     $     mhu0sq,mHpmsq,mA0sq,charmasstot,OCRm,OCLTm,Cegm)

!-----------------------------------------------------------------------------

      q = msusy
      p = q

      call neutralino(p,q,g,gp,mt_r,mb_r,mtau_r,tanbeta,yuRG,
     $     ydRG,yeRG,SUegg,SDegg,SLegg,SNegg,MNeut,Neg,Ceg,
     $     ON,OCL,OCR,mh0sq,mhu0sq,mHpmsq,mA0sq,ONew,Negm,
     $     neutmasstot)

!----------------------

      p = dsqrt(max(SUegg(1),SUegg(2)))
           
      call pistop(p,q,g,gp,g3,mt_r,mb_r,newtbeta,mSQRG,
     $     mSURG,yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,delpist,STegtmp)

!      print*,"pmax = ", dsqrt(steg(1)),dsqrt(steg(2))

      STeg(2) = max(STegtmp(1),STegtmp(2))

      p = dsqrt(min(SUegg(1),SUegg(2)))
           
      call pistop(p,q,g,gp,g3,mt_r,mb_r,newtbeta,mSQRG,
     $     mSURG,yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,delpist,STegtmp)

!      print*,"pmin = ", dsqrt(steg(1)),dsqrt(steg(2))

      STeg(1) = min(STegtmp(1),STegtmp(2))

      
!-----------------------------------------------------------------------------

      p = dsqrt(max(SDegg(1),SDegg(2)))
      
      call pisbottom(p,q,g,gp,g3,mt_r,mb_r,newtbeta,mSQRG,mSDRG,
     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,delpisbtm,SBegtmp)

!      print*,"pmax = ", dsqrt(sbeg(1)),dsqrt(sbeg(2))

      Sbeg(2) = max(SBegtmp(1),SBegtmp(2))

      p = dsqrt(min(SDegg(1),SDegg(2)))
      
      call pisbottom(p,q,g,gp,g3,mt_r,mb_r,newtbeta,mSQRG,mSDRG,
     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,delpisbtm,SBegtmp)

!      print*,"pmin = ", dsqrt(sbeg(1)),dsqrt(sbeg(2))

      Sbeg(1) = min(SBegtmp(1),SBegtmp(2))
      
      
!-----------------------------------------------------------------------------

      p = dsqrt(max(SUegg(3),SUegg(4)))
     
      call pischarm(p,q,g,gp,g3,newtbeta,mSQRG,
     $     mSURG,yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,delpisc,SCegtmp)

      SCeg(2) = max(SCegtmp(1),SCegtmp(2))

      p = dsqrt(min(SUegg(3),SUegg(4)))
     
      call pischarm(p,q,g,gp,g3,newtbeta,mSQRG,
     $     mSURG,yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,delpisc,SCegtmp)

      SCeg(1) = min(SCegtmp(1),SCegtmp(2))

!------------------------------------------------------------------------------

      p = dsqrt(max(SDegg(3),SDegg(4)))
      
      call pisstrange(p,q,g,gp,g3,newtbeta,mSQRG,mSDRG,
     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,delpisst,SSTegtmp)

      SSTeg(2) = max(SSTegtmp(1),SSTegtmp(2))

      p = dsqrt(min(SDegg(3),SDegg(4)))
      
      call pisstrange(p,q,g,gp,g3,newtbeta,mSQRG,mSDRG,
     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,delpisst,SSTegtmp)

      SSTeg(1) = min(SSTegtmp(1),SSTegtmp(2))
      
!----------------------------------------------------------------------------

      p = dsqrt(max(SUegg(5),SUegg(6)))

      call pisupq(p,q,g,gp,g3,tanbeta,mSQRG,
     $     mSURG,yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,SUqegtmp)


      SUqeg(2) = max(SUqegtmp(1),SUqegtmp(2))

      p = dsqrt(min(SUegg(5),SUegg(6)))

      call pisupq(p,q,g,gp,g3,tanbeta,mSQRG,
     $     mSURG,yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,SUqegtmp)

      SUqeg(1) = min(SUqegtmp(1),SUqegtmp(2))
            
!-----------------------------------------------------------------------------
      
      p = dsqrt(max(SDegg(5),SDegg(6)))

      call pisdown(p,q,g,gp,g3,newtbeta,mSQRG,mSDRG,
     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,SDnegtmp)

      SDneg(2) = max(SDnegtmp(1),SDnegtmp(2))

      p = dsqrt(min(SDegg(5),SDegg(6)))

      call pisdown(p,q,g,gp,g3,newtbeta,mSQRG,mSDRG,
     $     yuRG,ydRG,AURG,ADRG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,vev2,M3tz,SDnegtmp)

      SDneg(1) = min(SDnegtmp(1),SDnegtmp(2))

!----------------------------------------------------------------------------

      p = dsqrt(max(SLegg(1),SLegg(2)))
      
      call pistaul(p,q,g,gp,mtau_r,newtbeta,
     $     mSLRG,mSERG,yeRG,AERG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,M3tz,delpistau,STauegtmp)

c$$$      print*, "STaueg(1)M = ", dsqrt(STaueg(1))
c$$$      print*, "STaueg(2)M = ", dsqrt(STaueg(2))

      STaueg(2) = max(STauegtmp(1),STauegtmp(2))

      p = dsqrt(min(SLegg(1),SLegg(2)))
      
      call pistaul(p,q,g,gp,mtau_r,newtbeta,
     $     mSLRG,mSERG,yeRG,AERG,SUegg,SDegg,
     $     SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,M3tz,delpistau,STauegtmp)
      

      STaueg(1) = min(STauegtmp(1),STauegtmp(2))

c$$$      print*, "STaueg(1)m = ", dsqrt(STaueg(1))
c$$$      print*, "STaueg(2)m = ", dsqrt(STaueg(2))

!-----------------------------------------------------------------------------
      
      p =  dsqrt(max(SLegg(3),SLegg(4)))
      
      call pisMUl(p,q,g,gp,tanbeta,mSLRG,mSERG,
     $     yeRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,M3tz,SMuegtmp)

      SMueg(2) = max(SMuegtmp(1),SMuegtmp(2))

      p =  dsqrt(min(SLegg(3),SLegg(4)))
      
      call pisMUl(p,q,g,gp,tanbeta,mSLRG,mSERG,
     $     yeRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,
     $     sgnmu,modmu,ON,OCL,OCR,vev1,M3tz,SMuegtmp)

      SMueg(1) = min(SMuegtmp(1),SMuegtmp(2))

!-----------------------------------------------------------------------------

      p = dsqrt(max(SLegg(5),SLegg(6)))
      
      call pisEl(p,q,g,gp,tanbeta,mSLRG,mSERG,yeRG,
     $     AERG,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu,
     $     modmu,ON,OCL,OCR,vev1,M3tz,SEegtmp)

      SEeg(2) = max(SEegtmp(1),SEegtmp(2))

      p = dsqrt(min(SLegg(5),SLegg(6)))
      
      call pisEl(p,q,g,gp,tanbeta,mSLRG,mSERG,yeRG,
     $     AERG,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu,
     $     modmu,ON,OCL,OCR,vev1,M3tz,SEegtmp)
      
      SEeg(1) = min(SEegtmp(1),SEegtmp(2))
      
!------------------------------------------------------------------------------
      
      p =  dsqrt(Snegg(3))
      
      call pitausnu(p,q,g,gp,mtau_r,newtbeta,yeRG,AERG,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,
     $     mhpmsq,mA0sq,sgnmu,modmu,ON,OCL,OCR,vev1,tsnu)
      
!-------------------------------------------------------------------------------
      p = dsqrt(Snegg(2))
            
      call pimulsnu(p,q,g,gp,tanbeta,yeRG,AERG,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,
     $     mhpmsq,mA0sq,sgnmu,modmu,ON,OCL,OCR,vev1,musnu)

!------------------------------------------------------------------------------

      p =  dsqrt(Snegg(1))

      call pielsnu(p,q,g,gp,newtbeta,yeRG,AERG,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,mh0sq,mhu0sq,
     $     mhpmsq,mA0sq,sgnmu,modmu,ON,OCL,OCR,vev1,elsnu)

 
!---------------------------------------------------------------------------
C     Higgs mass matrix
C     from bpmz    

      mAm3 = bmur/ (sinbeta*cosbeta)  

!--------------------------      
      q = msusy
      
      call tadpole1(q,g,mb_r,mtau_r,newtbeta,
     $     ADRG,AERG,yuRG,ydRG,yeRG,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,ON,OCL,OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,
     $     mA0sq,vev1,delta2)

      q = msusy
      
      call tadpole2(q,g,mt_r,newtbeta,
     $     AURG,yuRG,ydRG,yeRG,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,ON,OCL,OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,
     $     mA0sq,vev2,delta1)

!----------------------------------------------------------------------------

      higgscorr = 1


      if(higgscorr>0)then
         
         itcount = 0

         mh0sqcor = 0.d0
         mhu0sqcor = 0.d0
         mhpmsqcor = 0.d0
         mA0sqcor = 0.d0

         call higgs(msusy,newtbeta,SUegg,SDegg,SLegg,
     $        SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu,modmu,
     $        bmur,piwwT,pizzT,delta1,delta2,mh0sqcor,mhu0sqcor,
     $        mhpmsqcor,mA0sqcor,alphacorr,itcount)


         newmHu0sq = mHu0sqcor

         newmh0sq = mh0sqcor

         newmA0sq = mA0sqcor

         newmHpmsq = mHpmsqcor

         alpha = alphacorr

      else

!----------------------------------------------------------------------------
C     Higgs 2 loop approximate expression

         mAm3 = bmur/ (sinbeta*cosbeta)  

         call higgs_analytical(MT,newtbeta,SUegg,AURG,sgnmu,modmu,
     $        alph3,sigphi1,sigphi2,sigphi12,sig2phi2,sig2phi2yuk)

C     C       Two loop approximation
         
         Mhiggsa(1,1) = (MZ*MZ + pizzT) * (cosbeta**2.d0)  + 
     $        (mAm3) * (sinbeta**2.d0) - sigphi1
         
         
         Mhiggsa(1,2) = -(MZ*MZ + pizzT + mAm3) * sinbeta*cosbeta -
     $        sigphi12
         
         Mhiggsa(2,1) = -(MZ*MZ + pizzT + mAm3) * sinbeta*cosbeta -
     $        sigphi12
         
         Mhiggsa(2,2) = (MZ*MZ + pizzT) * (sinbeta**2.d0) + 
     $        (mAm3) * (cosbeta**2.d0)  -  
     $        (sigphi2 + sig2phi2 + sig2phi2yuk)


         Call CEigensystem(2,Mhiggsa,2,mhiggsev,mhmix,2,-1)

         newmHu0sq = mhiggsev(1)
         newmh0sq = mhiggsev(2)
         alpha = datan(mhmix(1,2)/mhmix(1,1))


         newmA0sq = mAm3 - piaaT + (cosbeta**2.d0) * delta1 + 
     $        (sinbeta**2.d0) * delta2 
         
         newmhpmsq = newmA0sq + MW*MW + piaaT + piwwT - pihphmT

      endif

!----------      
      
      msnew = dsqrt(dsqrt(dabs(STeg(1)))*dsqrt(dabs(STeg(2))))

      RETURN
      END SUBROUTINE REWSBCOR

!================================================================================================================
C-----------------------------------------------------------
C     corrections to mt,md,mtau and gauge bosons at Mz
C------------------------------------------------------------    
****f* SuSeFLAV/ewsbiterate.f/coratmz
*  NAME
*    SUBROUTINE coratmz
*  SYNOPSIS
*    One loop correction to all standard model inputs at MZ.
*  FUNCTION
*    Calculates one-loop threshold corrections to $M_t, M_b, M_tau, \alpha_{em},
*    \alpha_s and \sin^2 \theta_w$ at MZ.
*    We closely follow BPMZ [hep-ph/9606211].
*  INPUTS
*     MW         - Mass of W boson.
*     MZ         - Mass of Z boson.
*    tanbeta     - the ratio of the vevs of the two Higgs doublet fields.
*    alphaDR     - \overbar{DR} Electromagnetic coupling constant.
*
*  RESULT
*     MWc_mz     - One loop corrected mass of W boson.
*     MZc_mz     - One loop corrected mass of Z boson.
*     MTc_mz     - One loop corrected top quark mass. (Expand on correction parameters)
*     mBc_mz     - One loop corrected bottom quark mass.
*     mTauc_mz   - One loop corrected tau lepton mass.
*     alphas1    - One loop corrected strong coupling constant.
*     alphaem    - One loop corrected electromagnetic coupling constant.
*     delalphaem - One loop correction to electromagnetic coupling constant.
*     delalphas  - One loop correction to strong coupling constant.
*     sinsqtheff - Corrected effective weak mixing angle.
*
*  EXAMPLE
*          SUBROUTINE coratMZ(MW,MZ,tanbeta,alphaDR,MWc_mz, MZc_mz,MTc_mz,
*     $     mBc_mz,mTauc_mz,alphas1,alphaem,delalphas, delalphem,
*     $     sinsqtheff,newvev,mbdrbar,itcount)
*
*  NOTES
*    Common blocks used in this routine
*
*      common/rgeopt_mz/mt_mz, mb_mz, mtau_mz
*      common/sminputs/ mbpole, mtaupole, Mtpole
*      common/higgsmixmz/ alphatree 
*      common/sinsq_mz/sinsqthw_mz
*      common/mu_mz/ murgemz
*      common/qcd_cor/mbmzdrbar
*      common/softout_mat_mz/mSQRGz,mSURGz,mSDRGz,AURGz,ADRGz,
*     $     AERGz,ANURGz,mSLRGz,mSERGz, mSNURGz,ONz,OCLz,OCRz,
*     $     MCharz, MNeutz
*      common/rgeoutput_MZ/ mh1mzz,mh2mzz,M1tmz,M2tmz,M3tmz,
*     $     alph1MZ,alph2MZ,alph3MZ,vev1mz,vev2mz,yumz,ydmz,yemz
*      common/sparticles_MZ/SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
*     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz 
* 
*    External Subroutines used in this routine
*
*      EXTERNAL topcor, bottomcor, taucor,pizz, piww
*      EXTERNAL strongcoupling, emcoupling, vevewsb, S2ThetaW,
*     $     pizgamma,S2ThetaEff
*      
*
*  BUGS
*    ---
*  SEE ALSO
*    -----
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!=================================================================================================================
      
      SUBROUTINE coratMZ(sgnmu,MW,MZ,tanbeta,alphaDR,MWc_mz, MZc_mz,
     $     MTc_mz,mBc_mz,mTauc_mz,alphas1,alphaem,delalphas, 
     $     delalphem,sinsqtheff,newvev,mbdrbar,flags)
      
 

      IMPLICIT NONE 

      INTEGER trysinsq,i
      character*100 flags,flags2tw
      DOUBLE PRECISION tanbeta,mbdrbar,sgnmu
      DOUBLE PRECISION MWc_mz,MZc_mz,MTc_mz,mBc_mz,mTauc_mz
      DOUBLE PRECISION g,gp,g3, p, q, modmu, newvev,pizzT0
      DOUBLE PRECISION pizzT, piwwT, vev1n,vev2n
      DOUBLE PRECISION alphaem, alphas1,correction,mbcor,mtaucor
      DOUBLE PRECISION delalphas, delalphem, alphatree

      double precision DeltaTZ, MT_qcd

      DOUBLE PRECISION mbpole, mtaupole, Mtpole
      DOUBLE PRECISION mt_mz, mb_mz, mtau_mz, alphaDR

      double precision pizgMZ, pizg0, sinsqthw_mz

      double precision sinthwold, rhold, piwwT0,sinsqtheff,
     $      mbmzdrbar,mbmzmsbar

      DOUBLE PRECISION yumz(3,3), ydmz(3,3), yemz(3,3)
      DOUBLE PRECISION alph3MZ, alph2MZ, alph1MZ
      DOUBLE PRECISION AURGz(3,3),ADRGz(3,3),AERGz(3,3)
      DOUBLE PRECISION mSQRGz(3,3), mSURGz(3,3),mSDRGz(3,3)
      DOUBLE PRECISION mSLRGz(3,3),mSNURGz(3,3), mSERGz(3,3)
      DOUBLE PRECISION M1tmz,M2tmz,M3tmz, mh1mzz,mh2mzz
      DOUBLE PRECISION murgemz, vev1mz,vev2mz, ANURGz(3,3)
      DOUBLE PRECISION ONz(4,4),OCLz(2,2),OCRz(2,2),
     $     MCharz(2,2), MNeutz(4,4)
      double precision mh0sqz,mhu0sqz,mhpmsqz,mA0sqz,Cegz(2),Negz(4)
      double precision SUeggz(6),SDeggz(6),SLeggz(6),SNeggz(3)
      DOUBLE PRECISION delta1z,delta2z,MWsqpole_MZ

      double precision MWpole, MZpole,pi,MW,MZ,mtaupoledrbar,
     $     mTauMZmsbar,mTauMZdrbar

!----------------
      common/rgeopt_mz/mt_mz, mb_mz, mtau_mz

      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole
      common/mtaurmass/ mtaupoledrbar,mTauMZmsbar,mTauMZdrbar
      common/mwpole/ MWpole
      common/higgsmixmz/ alphatree 
      common/sinsq_mz/sinsqthw_mz
      common/mu_mz/ murgemz
      common/qcd_cor/mbmzdrbar,mbmzmsbar

      common/softout_mat_mz/mSQRGz,mSURGz,mSDRGz,AURGz,ADRGz,
     $     AERGz,ANURGz,mSLRGz,mSERGz, mSNURGz,ONz,OCLz,OCRz,
     $     MCharz, MNeutz

      common/rgeoutput_MZ/ mh1mzz,mh2mzz,M1tmz,M2tmz,M3tmz,
     $     alph1MZ,alph2MZ,alph3MZ,vev1mz,vev2mz,yumz,ydmz,yemz

      common/sparticles_MZ/SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz


      common/finetuning/delta1z,delta2z
      common/deltarho/ piwwT0,pizzT0,pizzt,piwwt,MWsqpole_MZ
!--------------------
      EXTERNAL topcor, bottomcor, taucor,pizz, piww
      EXTERNAL strongcoupling, emcoupling, vevewsb, S2ThetaW,
     $     pizgamma,S2ThetaEff,tadpole1,tadpole2

!----------------

!-----------------------
      include 'stdinputs.h'
!------------------------

      MW = MWpole
      MZ = MZpole

      pi = 4.d0*datan(1.d0)

      g  = dsqrt(alph2MZ*16.d0*pi*pi)
      gp = dsqrt(alph1MZ*16.d0*pi*pi) * dsqrt(3.d0/5.d0) 
      g3 = dsqrt(alph3MZ*16.d0*pi*pi)

      modmu = murgemz


      sinthwold = dsqrt(gp**2.d0/(g**2.d0 + gp**2.d0))
      sinsqthw_mz = (gp**2.d0/(g**2.d0 + gp**2.d0))

!      print*,"mu at MZ = ", modmu


!------------------------
      q = MZ

      call strongcoupling(q,g3,mt_mz,M3tmz,SUeggz,SDeggz,
     $     alphas1,delalphas)


      q = MZ


      call emcoupling(q,alphaDR,mt_mz,SUeggz,SDeggz,SLeggz,Cegz,
     $     mHpmsqz,alphaem,delalphem)


!------------------------------------------------------------------------
      p = MZ

      q = MZ

      call pizz(p,q,g,mtpole,mb_mz,mtau_mz,tanbeta,SUeggz,SDeggz,
     $     SLeggz,SNeggz,Negz,Cegz,
     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz,ONz,OCLz,OCRz,
     $     sinsqthw_mz,pizzT)

!      print*,"pizzT = ", pizzT

      MZc_mz = dsqrt(MZ*MZ + pizzT)

!-----------------------------------------------------------------------

      call vevewsb(tanbeta,MZ,pizzT,gp,g,vev1n,vev2n,
     $     newvev)

c$$$      mt_mz = yuMZ(3,3) * vev2n/dsqrt(2.d0)
c$$$      mb_mz = ydMZ(3,3) * vev1n/dsqrt(2.d0)
c$$$      mtau_mz = yeMZ(3,3) * vev1n/dsqrt(2.d0)

!      print*,"mt_mz vev corr = ", mt_mz
!------------------------------------------------------------------------

      q = MZ

      p = MW

      call piww(p,q,g,mtpole,mb_mz,mtau_mz,tanbeta,SUeggz,SDeggz,
     $     SLeggz,SNeggz,Negz,Cegz,
     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz,ONz,OCLz,OCRz,
     $     sinsqthw_mz,piwwT)

      MWc_mz = dsqrt(MW*MW + piwwT)

!      print*,"MWc_MZ = ", MWc_mz, piwwt

!--------------

      q = MZ

      p = 2.d-5

      call piww(p,q,g,mtpole,mb_mz,mtau_mz,tanbeta,SUeggz,SDeggz,
     $     SLeggz,SNeggz,Negz,Cegz,
     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz,ONz,OCLz,OCRz,
     $     sinsqthw_mz,piwwT0)

!------------


      p = 2.d-5

      q = MZ

      call pizgamma(p,q,g,mt_mz,mb_mz,mtau_mz,SUeggz,SDeggz,
     $     SLeggz,Cegz,mhpmsqz,OCLz,OCRz,alphaDR/(1-delalphem),
     $     pizg0)


!---------

      p = MZ

      q = MZ

      call pizgamma(p,q,g,mt_mz,mb_mz,mtau_mz,SUeggz,SDeggz,
     $     SLeggz,Cegz,mhpmsqz,OCLz,OCRz,alphaDR/(1-delalphem),
     $     pizgMZ)

!-----------------

      p = 2.d-5

      q = MZ

      call pizz(p,q,g,mt_mz,mb_mz,mtau_mz,tanbeta,SUeggz,SDeggz,
     $     SLeggz,SNeggz,Negz,Cegz,
     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz,ONz,OCLz,OCRz,
     $     sinsqthw_mz,pizzT0)

      
!-------sinthw correction

      rhold = 1.d0
      trysinsq = 1
      flags2tw = 'AOK'
      q = MZ
      
      call S2ThetaW(q,g,gp,Negz,Cegz,SLeggz,SNeggz,ONz,OCLz,
     $     OCRz,sinthwold, rhold,delalphem,alphatree,pizzT,
     $     piwwT,piwwT0,mh0sqz,trysinsq,flags2tw)

c$$$      print*,"flags, flags2tw = ", flags, flags2tw

      if(flags2tw.ne.'AOK')then
         flags =' SW2NOC'       
      endif
      
!      print*,"sinthwold = ", sinthwold 
!----------------------

      call S2ThetaEff(pizgMZ, pizg0, 
     $     sinthwold**2.d0,alphaem, sinsqtheff)

!      sinsqtw = 0.2221d0

!      if(sinsqtheff.lt.sinsqtw) sinsqtheff = 0.2304d0
     
!      print*,"sintheff = ", dsqrt(sinsqtheff)
!-------------------------------------------------------------------------

!---------Tadpoles for fine tuning

      q = MZ

      call tadpole1(q,g,mb_mz,mtau_mz,tanbeta,
     $     ADRGz,AERGz,yuMZ,ydMZ,yeMZ,SUeggz,SDeggz,SLeggz,
     $     SNeggz,Negz,Cegz,ONz,OCLz,OCRz,sgnmu,modmu,mh0sqz,mhu0sqz,
     $     mHpmsqz,mA0sqz,vev1MZ,delta2z)

      
      q = MZ
      
      call tadpole2(q,g,mt_mz,tanbeta,AURGz,yuMZ,ydMZ,yeMZ,
     $     SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,ONz,OCLz,OCRz,
     $     sgnmu,modmu,mh0sqz,mhu0sqz,mHpmsqz,mA0sqz,vev2MZ,
     $     delta1z)
      
!----------------------------------------------------------------------------

      p = mtpole

      q = MZ 

      call topcor(p,q,g,gp,g3,M3tmz,mt_mz,mb_mz,tanbeta,
     $     yuMZ,ydMZ,SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz,ONz,OCLz,
     $     OCRz,sinsqthw_mz,correction)

!----------------------------------------------------------------------------


      p = mbmzMSbar

      q = MZ

      call bottomcor(p,q,g,gp,g3,M3tmz,mt_mz,mb_mz,tanbeta,
     $     yuMZ,ydMZ,SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz,ONz,OCLz,
     $     OCRz,sinsqthw_mz,mbcor,mbdrbar)
      

!----------------------------------------------------------------------------

      p = mtaupole

      q = MZ

      call taucor(p,q,g,gp,M3tmz,mtau_mz,tanbeta,yeMZ,
     $     SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,mh0sqz,
     $     mhu0sqz,mhpmsqz,mA0sqz,ONz,OCLz,
     $     OCRz,sinsqthw_mz,mtaucor)

!----------------------------------------------------------------------------
   

      alphas1 = (g3**2.d0/(4.d0 * pi))

      DeltaTZ = 2.d0 * dLog(Mt_mz/MZ)


      MT_qcd = -(alphas1 * (5.d0 - 3.d0 * DeltaTZ)/(3.d0 * pi)) -
     $    (alphas1**2.d0) * (0.538d0 - (43.d0 * DeltaTZ/(24.d0 * pi*pi)) 
     $     + (3.d0 * (DeltaTZ**2.d0)/(8.d0 * pi * pi)))

      MTc_mz =  correction/mtpole + MT_qcd

      mBc_mz =  mbcor

      mTauc_mz = mtaucor

!-----------------------------------------------------------------------


      call vevewsb(tanbeta,MZ,pizzT,gp,g,vev1n,vev2n,
     $     newvev)

!----------------------

      MWsqpole_MZ = ((g*g*(vev2n**2.d0+vev1n**2.d0)/4.d0) - piwwt)


      RETURN
      END SUBROUTINE coratMZ
!============================================================================

      recursive subroutine higgs(q,tanbeta,SUegg,SDegg,SLegg,
     $     SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu,mur,
     $     bmur,piwwT,pizzT,delta1,delta2,newmh0sq,newmhu0sq,
     $     newmhpmsq,newmA0sq,alphacorr,itcount)
      
      IMPLICIT NONE
      DOUBLE PRECISION mur,mt,MW,pihphmT,g3
      DOUBLE PRECISION pi,q,pizzT,MZ,p,piwwT
      DOUBLE PRECISION alph2,alph3,alph1,piaaT,vev1n,vev2n
      DOUBLE PRECISION modmu,delta1,delta2,tan2beta,newtbeta
      double precision gp,sinbeta,cosbeta
      DOUBLE PRECISION newmA0sq,newmhpmsq,newmh0sq,newmHu0sq
      DOUBLE PRECISION STeg(2),SBeg(2)
      DOUBLE PRECISION SCeg(2),SUqeg(2),STaueg(2),SEeg(2),SMueg(2)
      DOUBLE PRECISION SSTeg(2),SDneg(2),tsnu,musnu,elsnu

!      data piwwT/ 1 * 0.d0/, pizzT/ 1 * 0.d0/

      double precision M3t,deltagluino

      double precision MNeut1(4,4),Negm(4),neutmasstot(4,4),MChar(2,2)

      double precision charmasstot(2,2),Cegm(2),msnew     
      
c$$$      double precision Neutevc(4,4), msnew
c$$$      data Neutevc/ 16 * 0.d0/

!---------------------------------------------------------------------------
      integer i,j,OS,itcount
      
      parameter (OS = 0)

      double precision bmur,msusy,mh1mz,mh2mz,sgnmu
      double precision mSQRG(3,3),mSURG(3,3),mSDRG(3,3),mSLRG(3,3),
     $     mSERG(3,3),mSNURG(3,3)
      double precision yuRG(3,3),ydRG(3,3),yeRG(3,3)
      double precision ADRG(3,3),AURG(3,3),AERG(3,3), ANURG(3,3)
      double precision SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      double precision Neg(4),Ceg(2)
      double precision mh0sq,mhu0sq,mhpmsq,mA0sq
      double precision MNeut(4,4),ON(4,4)
      double precision OCR(2,2),OCL(2,2)
      double precision M1tz,M2tz,M3tz
      double precision vev1,vev2,cos2beta,g
      Double precision bmurcor,mAm3,mbpole, mtaupole, Mtpole
      data bmurcor/0.d0/

      DOUBLE PRECISION sinsqthw_susy
      DOUBLE PRECISION Mhiggsa(2,2),sigphi1,sigphi12,sigphi2,sig2phi2
      DOUBLE PRECISION newmha, mt_r, mb_r, mtau_r, sig2phi2yuk,newvev
      DOUBLE PRECISION delpisbtm(2,2),delpist(2,2),delpistau(2,2),
     $     ONew(4,4)
      DOUBLE PRECISION delpisc(2,2), delpisst(2,2), OCRm(2,2),OCLTm(2,2)
      DATA delpist/ 4 * 0.d0/,delpisbtm/ 4 * 0.d0/,delpistau/ 4 * 0.d0/
      DATA delpisc/ 4 * 0.d0/, delpisst/ 4 * 0.d0/

      DOUBLE PRECISION tanbeta

      double precision mhiggsev(2),mhmix(2,2),mhiggsevmHu(2)
      
      double precision pis1s1ans,pis2s2ans,pis1s2ans,pis1s1ansmHu,
     $     pis2s2ansmHu,pis1s2ansmHu

      double precision sigHu02loop(2,2),MhiggsmHu2loop(2,2)

      double precision thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae
      
      double precision mtsq_r,mgsq,mst1sq,mst2sq,st,ct,vv,
     $     S11stop,S22stop,S12stop,S11btm,S22btm,S12btm,p2btm
      double precision mbsq_r,msb1sq,msb2sq,sb,cb,qsq,cotbeta,
     $     mtausq_r,msnusq,mstau1sq,mstau2sq,stau,ctau,p2stop,
     $     S11top,S22top,S12top,p2top,S11tau,S22tau,S12tau,p2tau
      double precision sigh02loop(2,2), delmA, Mhiggs2loop(2,2)

      double precision Mhiggstree(2,2), heign(2),hl(2,2),alphatree
      double precision Mhiggsn(2,2)

      double precision MWpole, MZpole
      
      double precision mh0sqcor,mhu0sqcor,mA0sqcor,mHpmsqcor,rmh0,
     $     alphacorr

c$$$      double precision SUeggss,SDeggss,SLeggss,SNeggss,Negss,
c$$$     $     Cegss,mh0sqss,mhu0sqss,mhpmsqss,mA0sqss

!--------

      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole

      common/mwpole/MWpole

      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
     $     AERG,ANURG,mSLRG,mSERG,mSNURG,ON,OCL,OCR,MChar,MNeut

      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG

c$$$      common/sparticles_susy/SUeggss,SDeggss,SLeggss,SNeggss,Negss,
c$$$     $     Cegss,mh0sqss,mhu0sqss,mhpmsqss,mA0sqss

      common/sinsq_susy/sinsqthw_susy

!      common/opt_mixing/OCRm, OCLTm,ONew, alpha
!      common/vev_ewsb/ vev1n,vev2n
!      common/hzV/ delta1,delta2

      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,
     $     thetamu,thetau,thetad,thetae

!---------------------------------------------------------------------------

      EXTERNAL pizz,tadpole1,piww
      external tadpole2,piaa,vevewsb,pihphm,pistop
      EXTERNAL higgs_analytical

!---------------------------------------------------------------------------

!      print*,"itcount", itcount

      itcount = itcount + 1

!      print*,"itcount", itcount
      
      if(itcount>=10)then
         print*,"Non-Convergent Higgs Spectrum"
         return
      endif
         

      MZ = MZpole
      MW = MWpole

c$$$      mt_r = mt
c$$$      mb_r = mb
c$$$      mtau_r = mtau

      newtbeta = tanbeta
      tan2beta = dtan(2.d0*datan(newtbeta))
      sinbeta  = dsin(datan(newtbeta))
      cosbeta  = dcos(datan(newtbeta))
      cos2beta = dcos(2.d0*datan(newtbeta))

     
      pi = 4.d0 * datan(1.d0)

!      delta1 = 0.d0
!      delta2 = 0.d0
!      pizzT  = 0.d0
!      piwwT = 0.d0

      g  = dsqrt(alph2*16.d0*pi*pi)
      gp = dsqrt(alph1*16.d0*pi*pi*3.d0/5.d0) 
      g3 = dsqrt(alph3*16.d0*pi*pi)
      
      modmu = mur

      Mt_r = yuRG(3,3)* vev2/dsqrt(2.d0)
      Mb_r = ydRG(3,3)* vev1/dsqrt(2.d0)
      Mtau_r = yeRG(3,3)* vev1/dsqrt(2.d0)


!-------------------

      mA0sq = bmur/ (sinbeta*cosbeta) 

!      print*,"mA0sq, bmur, MZ, MW = ", mA0sq, bmur, MZ, MW 

      mhpmsq = mA0sq + MW*MW

!      print*,"pizzT, piwwT = ", pizzT, piwwT 


      Mhiggstree(1,1) = mA0sq * sinbeta**2.d0 + 
     $     (MZ*MZ + pizzT) * cosbeta**2.d0

      Mhiggstree(1,2) = -1.d0*(mA0sq + (MZ*MZ) + pizzT) * 
     $     sinbeta*cosbeta
      
      Mhiggstree(2,1) = Mhiggstree(1,2)

      Mhiggstree(2,2) = mA0sq * cosbeta**2.d0 + 
     $     (MZ*MZ + pizzT) * sinbeta**2.d0

      Mhiggsn(1,1) = Mhiggstree(1,1) 

      Mhiggsn(1,2) = Mhiggstree(1,2) 

      Mhiggsn(2,1) = Mhiggstree(2,1) 

      Mhiggsn(2,2) = Mhiggstree(2,2) 


      Call CEigensystem(2,Mhiggsn,2,heign,hl,2,-1)

      mh0sq = heign(2)       
      mhu0sq = heign(1)      
      
!      alphatree = datan(hl(1,2)/hl(1,1))

!-------------------


      if(itcount.eq.1)then
         
         mh0sqcor = mh0sq
         mhu0sqcor = mhu0sq
         mA0sqcor = mA0sq
         mHpmsqcor = mHpmsq
         
      else

         mh0sqcor = newmh0sq
         mhu0sqcor = newmhu0sq
         mA0sqcor = newmA0sq
         mHpmsqcor = newmHpmsq

      endif

!-----------------------------------------------------------------------

!      q = msusy

      p = dsqrt(((mh1mz - mh2mz)/(cos2beta)) - MZ*MZ) !(mA0sq)

      call piaa(p,q,g,gp,mt_r,mb_r,mtau_r,newtbeta,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,
     $     Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu,modmu,ON,OCL,OCR,
     $     piaaT)

!---------------------------------------------------------------------------    

      p = dsqrt(mHpmsqcor)
      

      call pihphm(p,q,g,gp,mt_r,mb_r,mtau_r,newtbeta,AURG,ADRG,AERG,
     $     SUegg,SDegg,SLegg,SNegg,Neg,Ceg,yuRG,
     $     yeRG,ydRG,mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu,modmu,ON,OCL,OCR,
     $     pihphmT)


!-----------------------------------------------------------------------------
C     Higgs mass matrix
C     from bpmz    
!----------------------------------------------------------------------------
      
!      q = msusy
      p = dsqrt(mh0sqcor)

      call pis1s1(p,q,g,gp,mt_r,mb_r,mtau_r,newtbeta,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,ON,OCL,
     $     OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,mA0sq,vev1,vev2,
     $     pis1s1ans)

!      print*,"pis1s1ans = ", pis1s1ans, p, q

      p = dsqrt(mhu0sqcor)

      call pis1s1(p,q,g,gp,mt_r,mb_r,mtau_r,newtbeta,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,ON,OCL,
     $     OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,mA0sq,vev1,vev2,
     $     pis1s1ansmHu)

!      print*,"pis1s1ansmHu = ", pis1s1ansmHu, p, q

!      q = msusy
      p = dsqrt(mh0sqcor)

      call pis2s2(p,q,g,gp,mt_r,mb_r,mtau_r,newtbeta,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,ON,OCL,
     $     OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,mA0sq,vev1,vev2,
     $     pis2s2ans)

!      print*,"pis2s2ans = ", pis2s2ans, p, q

      p = dsqrt(mhu0sqcor)

      call pis2s2(p,q,g,gp,mt_r,mb_r,mtau_r,newtbeta,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,ON,OCL,
     $     OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,mA0sq,vev1,vev2,
     $     pis2s2ansmHu)

!      print*,"pis2s2ansmHu = ", pis2s2ansmHu, p, q

!      q = msusy
      p = dsqrt(mh0sqcor)

      call pis1s2(p,q,g,gp,mt_r,mb_r,mtau_r,newtbeta,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,ON,OCL,
     $     OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,mA0sq,vev1,vev2,
     $     pis1s2ans)

!      print*,"pis1s2ans = ", pis1s2ans, p, q

      p = dsqrt(mhu0sqcor)

      call pis1s2(p,q,g,gp,mt_r,mb_r,mtau_r,newtbeta,yuRG,ydRG,yeRG,
     $     AURG,ADRG,AERG,SUegg,SDegg,SLegg,SNegg,Neg,Ceg,ON,OCL,
     $     OCR,sgnmu,modmu,mh0sq,mhu0sq,mHpmsq,mA0sq,vev1,vev2,
     $     pis1s2ansmHu)

!      print*,"pis1s2ansmHu = ", pis1s2ansmHu, p, q

!-------------------------------------------------------------------

!-------------------------------------------------------------------
!     2 loop higgs tadpole O(a_t a_s) correction from Slavich et al.
!-------------------------------------------------------------------
      
      mtsq_r = mt_r*mt_r
      mbsq_r = mb_r*mb_r
      mtausq_r = mtau_r*mtau_r

      mgsq = M3tz*M3tz

      mst1sq = Suegg(2)
      mst2sq = Suegg(1)

      msb1sq = Sdegg(2)
      msb2sq = Sdegg(1)

      msnusq = SNegg(1)

      mstau1sq = SLegg(2)
      mstau2sq = SLegg(1)

      st = dsin(thetat)
      ct = dcos(thetat)
      sb = dsin(thetab)
      cb = dcos(thetab)
      stau = dsin(thetatau)
      ctau = dcos(thetatau)

      qsq = q*q
      vv = (vev1**2.d0 + vev2**2.d0)
      cotbeta = 1.d0/newtbeta

      call DSZHiggs(mtsq_r,M3tz,mst1sq,mst2sq,st,ct,
     $     qsq,-1.d0*sgnmu*modmu,
     $     newtbeta,vv,g3,OS,S11stop,S22stop,
     $     S12stop)

!      print*,"S11stop,S22stop,S12stop = ", S11stop,S22stop,S12stop,os

      call DSZodd(mtsq_r,M3tz,mst1sq,mst2sq,st,ct,
     $     qsq,-1.d0*sgnmu*modmu,newtbeta,vv,g3,p2stop)

!      print*,"p2stop = ", p2stop

      call DSZHiggs(mbsq_r,M3tz,msb1sq,msb2sq,sb,cb,
     $     qsq,-1.d0*sgnmu*modmu,
     $     cotbeta,vv,g3,OS,S22btm,S11btm,S12btm)

!      print*,"S11btm,S22btm,S12btm = ", S11btm,S22btm,S12btm,os

      call DSZodd(mbsq_r,M3tz,msb1sq,msb2sq,sb,cb,
     $     qsq,-1.d0*sgnmu*modmu,cotbeta,vv,g3,p2btm)

!      print*,"p2btm = ", p2btm

      call DDSHiggs(mtsq_r,mbsq_r,mA0sq,mst1sq,mst2sq,msb1sq,
     $     msb2sq,st,ct,sb,cb,qsq,-1.d0*sgnmu*modmu,newtbeta,
     $     vv,S11top,S12top,S22top) 

!      print*,"S11top,S22top,S12top = ", S11top,S22top,S12top

      
      call DDSodd(mtsq_r,mbsq_r,mA0sq,mst1sq,mst2sq,msb1sq,msb2sq,
     $     st,ct,sb,cb,qsq,-1.d0*sgnmu*modmu,newtbeta,vv,p2top)

!      print*,"p2top = ", p2top

      call tausqHiggs(mtausq_r,mA0sq,msnusq,mstau1sq,mstau2sq,stau,
     $     ctau,qsq,-1.d0*sgnmu*modmu,newtbeta,vv,OS,S11tau,S22tau,
     $     S12tau)

!      print*,"S11tau,S22tau,S12tau = ", S11tau,S22tau,S12tau

      call tausqodd(mtausq_r,mA0sq,msnusq,mstau1sq,mstau2sq,stau,ctau,
     $     qsq,-1.d0*sgnmu*modmu,newtbeta,vv,p2tau)

!-------------------------------------------------------------------------

      sigh02loop(1,1) = pis1s1ans - S11stop - S11top - S11btm - S11tau
      sigh02loop(1,2) = pis1s2ans - S12stop - S12top - S12btm - S12tau
      sigh02loop(2,2) = pis2s2ans - S22stop - S22top - S22btm - S22tau
      sigh02loop(2,1) = sigh02loop(1,2)
      
      delmA = p2stop + p2top + p2btm + p2tau 


      sigHu02loop(1,1) = pis1s1ansmHu - S11stop - S11top - S11btm - 
     $     S11tau
      sigHu02loop(1,2) = pis1s2ansmHu - S12stop - S12top - S12btm - 
     $     S12tau
      sigHu02loop(2,2) = pis2s2ansmHu - S22stop - S22top - S22btm - 
     $     S22tau
      sigHu02loop(2,1) = sigHu02loop(1,2)

!      print*,"delmA = ", delmA

      Mhiggsa(1,1) = Mhiggstree(1,1) + delta2 
     $     + delmA * (sinbeta**2.d0)
      
      Mhiggsa(1,2) = Mhiggstree(1,2) - delmA * sinbeta*cosbeta  
      
      Mhiggsa(2,1) = Mhiggstree(2,1) - delmA * sinbeta*cosbeta  
      
      Mhiggsa(2,2) = Mhiggstree(2,2) + delta1 
     $     + delmA * (cosbeta**2.d0)


      Mhiggs2loop(1,1) = Mhiggsa(1,1) - sigh02loop(1,1)
      Mhiggs2loop(1,2) = Mhiggsa(1,2) - sigh02loop(1,2)
      Mhiggs2loop(2,1) = Mhiggsa(2,1) - sigh02loop(2,1)
      Mhiggs2loop(2,2) = Mhiggsa(2,2) - sigh02loop(2,2)

      MhiggsmHu2loop(1,1) = Mhiggsa(1,1) - sigHu02loop(1,1)
      MhiggsmHu2loop(1,2) = Mhiggsa(1,2) - sigHu02loop(1,2)
      MhiggsmHu2loop(2,1) = Mhiggsa(2,1) - sigHu02loop(2,1)
      MhiggsmHu2loop(2,2) = Mhiggsa(2,2) - sigHu02loop(2,2)


      Call CEigensystem(2,MhiggsmHu2loop,2,mhiggsevmHu,mhmix,2,0)

      Call CEigensystem(2,Mhiggs2loop,2,mhiggsev,mhmix,2,0)

      alphacorr = datan(mhmix(1,2)/mhmix(1,1))

      if(mhiggsev(2).gt.mhiggsev(1)) alphacorr = alphacorr + pi/2.d0

      newmHu0sq = max(mhiggsevmHu(1),mhiggsevmHu(2))

      newmh0sq = min(mhiggsev(1),mhiggsev(2))


!---------------------------------------

      newmA0sq = mA0sq - piaaT - pizzT + (cosbeta**2.d0) * delta1 + 
     $     (sinbeta**2.d0) * delta2 + delmA
      
      newmhpmsq = newmA0sq + MW*MW + piaaT + piwwT - pihphmT

      rmh0 = 1.d0 - (min(mh0sqcor,newmh0sq)/max(mh0sqcor,newmh0sq))


      if(rmh0.le.1.d-3)then

         return

      else
         
         call higgs(q,tanbeta,SUegg,SDegg,SLegg,
     $        SNegg,Neg,Ceg,mh0sq,mhu0sq,mhpmsq,mA0sq,sgnmu,mur,
     $        bmur,piwwT,pizzT,delta1,delta2,newmh0sq,newmhu0sq,
     $        newmhpmsq,newmA0sq,alphacorr,itcount)
      endif
      
      end subroutine higgs
      
!------------------------------------------------------------------------
!------------------------------------------------------------------------
