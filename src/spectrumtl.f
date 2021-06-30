****f* SuSeFLAV/spectrumtl.f/completerun 
*  NAME
*    Subroutine completerun
*  SYNOPSIS
*    Runs RGEs and computes tree level physical masses 
*    for MSSM parameters. 
*  FUNCTION
*    This suborutine integrates MSSM RGEs and computes low 
*    energy spectrum at msusy.  
*
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
*    Calculated tree level masses are stored in common blocks
*      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
*     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut
*      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
*     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG
*      common/sparticles_susy/SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
*     $     mh0sq,mhu0sq,mhpmsq,mA0sq
        
*  EXAMPLE
*
*   subroutine completerun(msusyold,vevin,yuin,ydin,yein,
*     $     alph1in,alph2in,alph3in,mur,bmur,murge,bmurge,prnstat,
*     $     check,newtbeta,msusynew,mursq,try,flags,runum,itcount)   
*
*  NOTES
*     Common Blocks used:
*      common/mu_rge/mufrge
*      common/loops/ lopt,rhn
*      common/sminputs/ mbpole, mtaupole, Mtpole
*      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
*      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
*      common/charinputs/case, model
*      common/rgeinput_high/ MX, M1X,M2X,M3X,mh10,mh20,
*     $     mQ0,mU0,mD0,mE0,mL0,mNU0
*      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
*     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut
*      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
*     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG
*      common/sparticles_susy/SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
*     $     mh0sq,mhu0sq,mhpmsq,mA0sq
*      common/soft_mat_susy/ mSL1, SLeg, USL, SNeg, USN, SUeg, USU,
*     $     SDeg, USD
*      common/runningmass_susy/MT_susy, MB_susy, Mtau_susy
*
*     External Routine used:
*      EXTERNAL MSSMRUN,mutreelevel,softspectrum
*
*  BUGS
*    ---
*  SEE ALSO
*    -----
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!=============================================================================================

      subroutine completerun(msusyold,vevin,yuin,ydin,yein,
     $     alph1in,alph2in,alph3in,mur,bmur,murge,bmurge,prnstat,
     $     check,newtbeta,msusynew,mursq,try,flags,runum,itcount)
   


      IMPLICIT NONE
      
      integer rhn,itcount
      DOUBLE PRECISION msusyold,msusynew,sgnmu
      double precision bmur,mur,bmurge,murge
      double precision M1X,M2X,M3X,MX,mh1mz,mh2mz,mh10,mh20,mursq
      double precision M1tz,M2tz,M3tz

      double precision MT_susy,MB_susy,mTau_susy

      double precision m0,m12,tanbeta,a0,ue3
      double precision mh0sq,mhu0sq,mhpmsq,mA0sq
      double precision mq0(3,3),mu0(3,3),MR1,MR2,MR3
      double precision md0(3,3),ml0(3,3),me0(3,3),mnu0(3,3) !,b1,b2,b3
      double precision msusy    !,alph10,alph20,alph30,tZ,tq0
      double precision ANURG(3,3)
      double precision yuRG(3,3),ydRG(3,3),yeRG(3,3),Ynui(3,3)

      double precision SUeg(6),SDeg(6),SLeg(6),SNeg(3)
      double precision mSQRG(3,3),mSURG(3,3),mSDRG(3,3),mSLRG(3,3),
     $     mSERG(3,3),mSNURG(3,3)
      double precision ADRG(3,3),AURG(3,3),AERG(3,3)
      double precision SUegg(6),SDegg(6),SLegg(6),SNegg(3)
      double precision USD(6,6),USL(6,6),USU(6,6),USN(3,3) !,mHcharsq, mHchar,Bbsg
      
      DOUBLE PRECISION alph1,alph2,alph3,newtbeta

      DOUBLE PRECISION Mtaupole,MTpole,Mbpole,MZpole
      

      double precision m10,m20
      double precision vev1,vev2
      double precision vevin
      double precision alph1in,alph2in,alph3in
      double precision yuin(3,3),ydin(3,3),yein(3,3)

      double precision MChar(2,2)
      double precision Neuevi(4),ONL(4,4),ONR(4,4),MNeut(4,4),
     $     Neg(4),ON(4,4)
      double precision OCR(2,2),OCL(2,2),Ceg(2)
      DOUBLE PRECISION pi
      DOUBLE PRECISION mSL1(6,6)
!------------------------------
      INTEGER runum, try
      integer prnstat,lopt,AOK,check
      integer flagdet, flagDflat
      
      character*3 case
      character*4 model
      CHARACTER*100 flags
      double precision mufrge, MWsq_susy, MZsq_susy
!---------------------------------

      common/mu_rge/mufrge
      common/loops/ lopt,rhn
      common/sminputs/ mbpole, mtaupole, Mtpole,MZpole
      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
      common/mssmrhn/ MR1,MR2,MR3,ue3, Ynui
      common/charinputs/case, model
      common/rgeinput_high/ MX, M1X,M2X,M3X,mh10,mh20,
     $     mQ0,mU0,mD0,mE0,mL0,mNU0

      common/softout_mat/mSQRG,mSURG,mSDRG,AURG,ADRG,
     $     AERG,ANURG,mSLRG,mSERG, mSNURG,ON,OCL,OCR,MChar, MNeut

      common/rgeoutput_susy/ mh1mz,mh2mz,M1tz,M2tz,M3tz,
     $     alph1,alph2,alph3,vev1,vev2,yuRG,ydRG,yeRG

      common/sparticles_susy/SUegg,SDegg,SLegg,SNegg,Neg,Ceg,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq

      common/soft_mat_susy/ mSL1, SLeg, USL, SNeg, USN, SUeg, USU,
     $     SDeg, USD
      common/runningmass_susy/MT_susy, MB_susy, Mtau_susy


!-----------------------------------------------------------------
      EXTERNAL MSSMRUN,mutreelevel,softspectrum
!-----------------------------------------------------------------

      PI = 4.d0 * datan(1.d0)

      msusy = msusyold

      check = 0
!--------------------------------

c$$$      print*,"rge input mu = ", mur 
c$$$      print*,"rge input bmu = ", bmur 

      call MSSMRUN(vevin,yuin,ydin,yein,alph1in,alph2in,alph3in,
     $     mur,bmur,murge,bmurge,prnstat,check,newtbeta,
     $     msusy,runum, itcount,flags)  

      if(flags.eq.'variable underflow '.or.flags.eq.'NPERTYUK')then
         return
      endif

c$$$      print*,"rge output mu = ", murge 
c$$$      print*,"rge output bmu = ", bmurge 
!      print*,"rge output of bmu = ", 2.d0 * murge**2.d0 + mh1mz + mh2mz

      mufrge = murge

      if(try.ne.1)then

         call mutreelevel(newtbeta,mh1mz,mh2mz,mursq,murge,
     $        bmurge,flagdet,flagDflat)


         if(mursq.lt.(0.d0).or.bmurge.lt.0.d0)then

            murge = dsqrt(dabs(mursq)) 
            bmurge = dabs(bmurge)

         endif

      elseif(try.eq.1)then
c$$$         if(sgnmu.lt.0.d0)then
c$$$            murge = -10.d0 
c$$$         else
c$$$            murge = 10.d0
c$$$         endif
         murge = 10.d0
         bmurge = 10.d0 
         
      endif

      MWsq_susy = (alph2*4.d0*pi*pi)*(vev1**2.d0 + vev2**2.d0)
      
      MZsq_susy = ((alph2*4.d0*pi*pi) + 
     $     (alph1*4.d0*pi*pi) * (3.d0/5.d0)) *
     $     (vev1**2.d0 + vev2**2.d0)
      
      MT_susy = yuRG(3,3) * vev2/ dsqrt(2.d0)
      MB_susy = ydRG(3,3) * vev1 / dsqrt(2.d0)
      Mtau_susy = yeRG(3,3) * vev1/ dsqrt(2.d0)
      
c$$$      print*,"MW in completerun = ", dsqrt(MWsq_susy)
c$$$      print*,"MZ in completerun = ", dsqrt(MZsq_susy)
c$$$      print*,"yt in completerun = ", yuRG(3,3)
c$$$      print*,"yb in completerun = ", ydRG(3,3) 
c$$$      print*,"ytau in completerun = ", yeRG(3,3)
c$$$      print*,"treelevel output mu = ", murge 
c$$$      print*,"treelevel output bmu = ", bmurge 


      call softspectrum(msusy,sgnmu,newtbeta,mSQRG,mSDRG,mSURG,AURG,
     $     ADRG,vev1,vev2,MWsq_susy,MZsq_susy,
     $     mSLRG,mSERG,AERG,yuRG,yeRG,ydRG,M1tz,M2tz,murge,bmurge,
     $     SUegg,USU,SUeg,
     $     SDegg,USD,SDeg,SLegg,USL,SLeg,SNegg,USN,SNeg,MNeut,ON,Neg,
     $     MChar,mSL1,OCR,OCL,Ceg,AOK,MT_susy,MB_susy,Mtau_susy,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,Neuevi,ONL,ONR)

c$$$      print*,"SUegg(1), SUegg(2) = ", SUegg(1), SUegg(2)
c$$$      print*,"SDegg(1), SDegg(2) = ", SDegg(1), SDegg(2)
c$$$      print*,"SLegg(1), SLegg(2) = ", SLegg(1), SLegg(2)
      
      
      
      msusynew = dsqrt(dsqrt(dabs(SUegg(1))*dabs(SUegg(2))))         



      RETURN

      END SUBROUTINE completerun
      
!============================================================================================
****f* SuSeFLAV/spectrumtl.f/runtomz 
*  NAME
*    Subroutine runtomz
*  SYNOPSIS
*    Runs RGEs from msusy to MZ and computes tree level 
*    physical masses for MSSM parameters. 
*  FUNCTION
*    This suborutine integrates MSSM RGEs from msusy to MZ and 
*    computes low energy spectrum at msusy.  
*
* *  INPUTS
*     MX         -  Reference scale, 10^19(GeV)
*     msusy      -  susy breaking scale. 
*     mu_conv    -  Converged value of \mu at msusy.
*     bmur_conv  -  Converged value of b_\mu at msusy.
*
*  RESULT
*     murgemz    -  RGE output: \mu at M_z scale
*     bmurgemz   -  RGE output: b_\mu at M_z scale
*     newtbetamz -  Ratio of vev at MZ from rge running.
*     flags      - flags any problem with the running 
*
*    Calculated tree level masses at MZ are stored in common blocks
*      common/rgeoutput_MZ/ mh1mzz,mh2mzz,M1tmz,M2tmz,M3tmz,
*     $     alph1MZ,alph2MZ,alph3MZ,vev1mz,vev2mz,yumz,ydmz,yemz
*      common/sparticles_MZ/SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
*     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz
*      common/softout_mat_mz/mSQRGz,mSURGz,mSDRGz,AURGz,ADRGz,
*     $     AERGz,ANURGz,mSLRGz,mSERGz, mSNURGz,ONz,OCLz,OCRz,
*     $     MCharz, MNeutz
*      common/soft_mat_mz/ mSL1z, SLegz, USLz, SNegz, USNz, SUegz, USUz,
*     $     SDegz, USDz
*      common/rgeopt_mz/mt_mz, mb_mz, mtau_mz
*        
*  EXAMPLE
*
*      SUBROUTINE runtomz(try,MX,msusy,mu_conv,bmur_conv,
*     $     murgemz,bmurgemz,newtbetamz,flags)
*
*  NOTES
*     Common Blocks used:
*      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
*      common/rgeoutput_MZ/ mh1mzz,mh2mzz,M1tmz,M2tmz,M3tmz,
*     $     alph1MZ,alph2MZ,alph3MZ,vev1mz,vev2mz,yumz,ydmz,yemz
*      common/sparticles_MZ/SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
*     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz
*      common/softout_mat_mz/mSQRGz,mSURGz,mSDRGz,AURGz,ADRGz,
*     $     AERGz,ANURGz,mSLRGz,mSERGz, mSNURGz,ONz,OCLz,OCRz,
*     $     MCharz, MNeutz
*      common/soft_mat_mz/ mSL1z, SLegz, USLz, SNegz, USNz, SUegz, USUz,
*     $     SDegz, USDz
*      common/rgeopt_mz/mt_mz, mb_mz, mtau_mz
*      common/runningew/MWsq_mz,MZsq_mz     
*
*    External routines used:
*      EXTERNAL MSSM_MZ, softspectrum
*
*  BUGS
*    ---
*  SEE ALSO
*    -----
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!----------------------------------------------------------------------------------
      SUBROUTINE runtomz(MX,msusy,mu_conv,bmur_conv,
     $     murgemz,bmurgemz,newtbetamz,flags)

      IMPLICIT NONE

      INTEGER AOK,i
      character*100 flags
      DOUBLE PRECISION MX, msusy, mu_conv, bmur_conv
      DOUBLE PRECISION yumz(3,3), ydmz(3,3), yemz(3,3)
      DOUBLE PRECISION alph3MZ, alph2MZ, alph1MZ,MZ
      DOUBLE PRECISION AURGz(3,3),ADRGz(3,3),AERGz(3,3)
      DOUBLE PRECISION mSQRGz(3,3), mSURGz(3,3),mSDRGz(3,3)
      DOUBLE PRECISION mSLRGz(3,3),mSNURGz(3,3), mSERGz(3,3)
      DOUBLE PRECISION bmurgemz,M1tmz,M2tmz,M3tmz, mh1mzz,mh2mzz
      DOUBLE PRECISION murgemz, newtbetamz,vev1mz,vev2mz, ANURGz(3,3)


      double precision mh0sqz,mhu0sqz,mhpmsqz,mA0sqz, MWsq_mz, MZsq_mz
      double precision SUegz(6),SDegz(6),SLegz(6),SNegz(3)
      double precision SUeggz(6),SDeggz(6),SLeggz(6),SNeggz(3)
      double precision USDz(6,6),USLz(6,6),USUz(6,6),USNz(3,3) 
      double precision MCharz(2,2)
      double precision Neueviz(4),ONLz(4,4),ONRz(4,4),MNeutz(4,4),
     $     Negz(4),ONz(4,4)
      double precision OCRz(2,2),OCLz(2,2),Cegz(2)
      DOUBLE PRECISION mSL1z(6,6),pi
      DOUBLE PRECISION MT_mz, MB_mz, mtau_mz
      DOUBLE PRECISION m0, m12, m10, m20, sgnmu, tanbeta, a0
!---------------

      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0

      common/rgeoutput_MZ/ mh1mzz,mh2mzz,M1tmz,M2tmz,M3tmz,
     $     alph1MZ,alph2MZ,alph3MZ,vev1mz,vev2mz,yumz,ydmz,yemz

      common/sparticles_MZ/SUeggz,SDeggz,SLeggz,SNeggz,Negz,Cegz,
     $     mh0sqz,mhu0sqz,mhpmsqz,mA0sqz


      common/softout_mat_mz/mSQRGz,mSURGz,mSDRGz,AURGz,ADRGz,
     $     AERGz,ANURGz,mSLRGz,mSERGz, mSNURGz,ONz,OCLz,OCRz,
     $     MCharz, MNeutz


      common/soft_mat_mz/ mSL1z, SLegz, USLz, SNegz, USNz, SUegz, USUz,
     $     SDegz, USDz

      common/rgeopt_mz/mt_mz, mb_mz, mtau_mz
      common/runningew/MWsq_mz,MZsq_mz     

!---------
      EXTERNAL MSSM_MZ, softspectrum
!-------

C     softspectrum at MZ

      pi = 4.d0 * datan(1.d0)



      CALL MSSM_MZ(MX,msusy,mu_conv,bmur_conv,
     $     murgemz,bmurgemz,newtbetamz,flags)

      if(flags.eq.'variable underflow ')then
         return
      endif

      mt_mz = yuMZ(3,3)*vev2mz /dsqrt(2.d0)
      mb_mz = ydMZ(3,3)*vev1mz /dsqrt(2.d0)
      mtau_mz = yeMZ(3,3)*vev1mz /dsqrt(2.d0)

c$$$      print*,"mt_mz = ", mt_mz, "yu_mz = ", yuMZ(3,3)
c$$$      print*,"mb_mz = ", mb_mz, "yb_mz = ", ydMZ(3,3)
c$$$      print*,"mtau_mz = ", mtau_mz, "ytau_mz = ", yeMZ(3,3)
c$$$      print*,"vev2mz = ", vev2mz
c$$$      print*,"vev1mz = ", vev1mz
c$$$      
c$$$      print*,"RGE output of mu at MZ = ", murgemz
c$$$      print*,"RGE output of bmu at MZ = ", bmurgemz

      MWsq_MZ = (alph2MZ*4.d0*pi*pi)*(vev1mz**2.d0 + vev2mz**2.d0)

      MZsq_MZ = ((alph2MZ*4.d0*pi*pi) + 
     $     (alph1MZ*4.d0*pi*pi) * (3.d0/5.d0)) *
     $     (vev1mz**2.d0 + vev2mz**2.d0)

      MZ = 91.2d0

      call softspectrum(MZ,sgnmu,newtbetamz,mSQRGz,mSDRGz,mSURGz,
     $      AURGz,ADRGz,vev1mz,vev2mz,MWsq_mz, MZsq_mz,
     $      mSLRGz,mSERGz,AERGz,yuMZ,yeMZ,ydMZ,M1tmz,M2tmz,
     $      murgemz,bmurgemz,
     $      SUeggz,USUz,SUegz,
     $      SDeggz,USDz,SDegz,SLeggz,USLz,SLegz,SNeggz,USNz,SNegz,
     $      MNeutz,ONz,Negz,
     $      MCharz,mSL1z,OCRz,OCLz,Cegz,AOK,MT_mz,MB_mz,Mtau_mz,
     $      mh0sqz,mhu0sqz,mhpmsqz,mA0sqz,Neueviz,ONLz,ONRz)

c$$$      print*,"After softspectrum at MZ"

c$$$ 100  format(/A/,("Eign(",I1,") = ", 1x, 1pe11.4))
c$$$      
c$$$!      print* ,"q,g,gp = ",q,g,gp
c$$$      print 100,"Negz: ", (i, Negz(i), i = 1, 4)
c$$$      print 100,"Cegz: ", (i, Cegz(i), i = 1, 2)
c$$$      print 100,"SLeggz: ", (i, SLeggz(i), i = 1, 6)
c$$$      print 100,"SNeggz: ", (i, SNeggz(i), i = 1, 3)
c$$$      print*,"==========================="

      RETURN
      END SUBROUTINE runtomz

!==========================================================
