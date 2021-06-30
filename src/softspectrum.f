****f* SuSeFLAV/softspectrum.f 
*  NAME
*    softspectrum
*  SYNOPSIS
*    Subroutine to generate Susy spectrum. 
*  FUNCTION
*     THIS PROGRAM CALCULATES THE MASSES AND MIXINGS IN THE SOFT
*     SECTOR FOR A GIVEN SPARTICLE MASSES AT LOW ENERGIES. THIS
*     PARAMETERS CAN BE EITHER OUTPUT FROM RGE OR DIRECT INPUTS.
*     AT PRESENT, IT HAS MAINLY SLEPTONIC AND SNEUTRINO PART ALONG
*     WITH CHARGINO AND NEUTRALINO PARTS IN THE BASIS WHERE CHARGED
*     LEPTONS ARE CONSIDERED DIAGONAL.
*     Higgs Spectrum added. Last Modified : 26/01/10. 
*
*  INPUTS
*     tanbeta  = value of the tanbeta being used. 
*     mSLRG    = (3 X 3) real mass matrix of left-handed sleptons (L)
*     mSERG    = (3 X 3) real mass matrix of right-handed sleptons (E^c)
*     AERG     = (3 X 3) real mass matrix of leptonic A-parameters.
*     M1tz     = low-energy (~M_Z) value of parameter M1
*     M2tz     = low-energy (~M_Z) value of parameter M2
*     mur      =  low-energy value of the mu-parameter either through REWSB
*            or as a free parameter.
*
*  RESULT
*     SUegg = 6 eigenvalues (ascending order) of UP-Squark mass matrix.
*     USU = (6 X 6) real orthogonal matrix such that,
*            USU*MSU*Transpose[USU] = Diag[MSU].
*     SUeg = 6 singular values (descending order) of UP-Squark mass
*            mass matrix. All positive.
*
*     SDegg = 6 eigenvalues (ascending order) of Down-Squark mass matrix.
*     USD = (6 X 6) real orthogonal matrix such that,
*            USD*MSD*Transpose[USD] = Diag[MSD].
*     SDeg = 6 singular values (descending order) of Down-Squark mass
*            mass matrix. All positive.
*
*     Slegg = 6 eigenvalues (ascending order) of slepton mass matrix.
*     USL = (6 X 6) real orthogonal matrix such that,
*            USL*MSL*Transpose[USL] = Diag[MSL].
*     SLeg = 6 singular values (descending order) of slepton mass
*            mass matrix. All positive.
*	
*     SNegg = 3 eigenvalues (ascending order) of sneutrino mass matrix.
*     USN = (3 X 3) real orthogonal matrix such that
*             USN*mSN*Transpose[USN] = Diag[MSN].
*     SNeg = 3 singular values (descending order) of sneutrino mass matrix.
*
*     ON = (4 X 4) orthogonal matrix such that 
*           ON.MNeut.Transpose[ON] = Diag[MNeut] 
*     Neg = 4 singular values (descending order) of the Neutralino mass matrix. 
*
*     OCR, OCL = (2 X 2) orthogonal matrices such that 
*                MChar = Transpose[OCR].Diag[MChar].OCL
*     Ceg = 2 singular values of the Neutralino Mass Matrix
*
*     AOK = Tells us whether all the diagonalising routines have run alright
*           or not. Without the Higgs spectrum, AOK should be 10 on output, 
*	    if everything goes well. 
*
*  EXAMPLE
*      subroutine softspectrum(tanbeta,mSQRG,mSDRG,mSURG,AURG,ADRG,vev1,
*     $     vev2,mSLRG,mSERG,AERG,yuRG,yeRG,ydRG,M1tz,M2tz,mur,SUegg,USU,
*     $     SUeg,SDegg,USD,SDeg,SLegg,USL,SLeg,SNegg,USN,SNeg,ON,Neg,OCR,
*     $     OCL,Ceg,AOK,MT_susy,MB_susy,Mtau_susy,mh0sq,
*     $     mhu0sq,mhpmsq,mA0sq,Neuevi,ONL,ONR)
*  
*  NOTES
*     ALL Masses and Parameters are in GeV.
*     Lapack is required.
*     Note that Lapack returns singular values as A = U Diag [A] Transpose[V]
*     or Transpose[U] A V = Diag[A]  which is transpose of the conventions of 
*     Mathematica as well as our conventions. (See notes). Our conventions
*     coincide with that of HN.  Also compare with that of Haber-Kane.  
*     
*  BUGS
*     
*     Possibly there are no bugs in this subroutine.
*     
*  SEE ALSO
*     -----
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!===========================================================================================
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
C     1. Checked on 28/06/2010
C     2. Removed the unused variables.
C     3. Removed the dependency on 'stddef.dat'.
C\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


      subroutine softspectrum(scale,sgnmu,tanbeta,mSQRG,mSDRG,mSURG,
     $     AURG,ADRG,vev1,vev2,MWsq, MZsq,
     $     mSLRG,mSERG,AERG,yuRG,yeRG,ydRG,M1tz,M2tz,mur,bmur,
     $     SUegg,USU,
     $     SUeg,SDegg,USD,SDeg,SLegg,USL,SLeg,SNegg,USN,SNeg,MNeut,ON,
     $     Neg,MChar,mSL2,OCR,OCL,Ceg,AOK,MT_susy,MB_susy,Mtau_susy,
     $     mh0sq,mhu0sq,mhpmsq,mA0sq,Neuevi,ONL,ONR)


      implicit none 

      integer m,a,i,j
      integer AOK,lwork1,lwork2,lwork3,info,lwork
      character*10 flag_bmu

      parameter(lwork=35,lwork1=15,lwork2=20,lwork3=10)

      double precision mSL(6,6),mSL2(6,6),mSLRG(3,3),AERG(3,3)
      double precision MChar(2,2),tanbeta, mur, M1tz, SLeg(6), SLegi(6)
      double precision MNeut(4,4),M2tz,mSN(3,3),USL(6,6),mSERG(3,3)
      double precision USN(3,3),SNeg(3),Neg(4),ON(4,4), OCRT(2,2)
      double precision work(lwork), USNT(3,3), OCR(2,2), OCL(2,2)
      double precision OCLT(2,2),mSL1(6,6),Ceg(2),mSN1(3,3)
      double precision work1(lwork1),USLT(6,6),work2(lwork2)
      double precision sw2,stw,ctw,beta,sbeta,cbeta,work3(lwork3)
      double precision vev1,vev2,cos2beta,SLegg(6),SNegg(3),mSU2(6,6)
      double precision mSU(6,6),mSQRG(3,3),mSURG(3,3),AURG(3,3)
      double precision mSU1(6,6),USU(6,6),USUT(6,6),SUegg(6),SUeg(6)
      double precision mSD(6,6),mSDRG(3,3),ADRG(3,3),mSD1(6,6),mSD2(6,6)
      double precision USD(6,6),USDT(6,6),SDegg(6),SDeg(6),slegs(6)

      double precision AURGn(3,3),ADRGn(3,3),AERGn(3,3)
      double precision MWsq,MZsq,scale,mSD3(6,6)
      double precision mh0sq,mhu0sq,mhpmsq,mA0sq
      DOUBLE PRECISION Mhiggstree(2,2), Mhiggsn(2,2), alphatree
      DOUBLE PRECISION heign(2), lhim(2),hl(2,2),hr(2,2)

      double precision SNeggt(3)


      double precision yuRG(3,3),ydRG(3,3),yeRG(3,3)
      double precision mT, mB, mTau !,SUeggt(6)
      double precision MNeut1(4,4),MChar1(2,2),bmur
      double precision MNeut2(4,4),Neuevi(4),ONL(4,4),ONR(4,4)

      double precision SUeggi(6)
      
      DOUBLE PRECISION sups(2,2), schs(2,2), stops(2,2),supegg(2),
     $     schegg(2), stopegg(2)

      DOUBLE PRECISION sdowns(2,2),sstrs(2,2), sbtms(2,2),sdownegg(2),
     $     sstregg(2), sbtmegg(2)

      DOUBLE PRECISION staus(2,2),smus(2,2), sels(2,2),staugg(2),
     $     smugg(2), selgg(2)


      DOUBLE PRECISION  thetaocl

      double precision  Mceg(2,2), thetaocr, OCRdum(2,2),OCLdum(2,2)
      double precision mmdag(2,2), mdagm(2,2), MCharT(2,2)
      DOUBLE PRECISION eigcgR(2,2),eigcgL(2,2)
      
      DOUBLE PRECISION supsi(2,2), supL(2,2), supR(2,2)
      DOUBLE PRECISION schsi(2,2), schL(2,2), schR(2,2)
      DOUBLE PRECISION stopsi(2,2), stopL(2,2), stopR(2,2)
      
      DOUBLE PRECISION sdownsi(2,2), sdownL(2,2), sdownR(2,2)
      DOUBLE PRECISION sstrsi(2,2), sstrL(2,2), sstrR(2,2)
      DOUBLE PRECISION sbtmsi(2,2), sbtmL(2,2), sbtmR(2,2)

      DOUBLE PRECISION selsi(2,2), selL(2,2), selR(2,2)
      DOUBLE PRECISION smusi(2,2), smuL(2,2), smuR(2,2)
      DOUBLE PRECISION stausi(2,2), stauL(2,2), stauR(2,2)

      DOUBLE PRECISION rrstaumu(2)
            
      DOUBLE PRECISION USLs(6,6),USLTs(6,6),MSL3(6,6)
!     $     USLsrt(6,6),USLTsrt(6,6)
           
      double precision sgnmu
      
      double precision slpmix(6,6),slpmixT(6,6),slpegn(6)
      double precision SDegi(6),sdmix(6,6),sdmixT(6,6),sdegn(6)
      double precision sumix(6,6),sumixT(6,6),suegn(6)
!     double precision ONT(4,4),SUeggi(6),work2(lwork2)
      
      integer ILO,IHI,IWORK
      double precision SCALE1(6),ABNRM, RCONDE(6),RCONDV(6)

!     MT at the susy scale from the running
      double precision MT_susy, MB_susy, Mtau_susy
!      common/suptree/ SUeggt
      common/higgsmixmz/ alphatree 
!      common/mutaup/rrstaumu,USLsrt,USLTsrt
      external dsyev, dgesvd , matmult2d, mat3prod2d,dgeevx
      common/mAflag/flag_bmu

!--------mixing angles 2x2

      DOUBLE PRECISION thetat,thetab,thetatau
      DOUBLE PRECISION thetac,thetas,thetamu
      DOUBLE PRECISION thetau,thetad,thetae

      DOUBLE PRECISION thetatz,thetabz,thetatauz
      DOUBLE PRECISION thetacz,thetasz,thetamuz
      DOUBLE PRECISION thetauz,thetadz,thetaez

      double precision SLegvd(6),SUegvd(6),SDegvd(6)
      double precision USLvd(6,6),USUvd(6,6),USDvd(6,6)

      double precision traceb, detb,eigb1,eigb2
      double precision traceL,DETL,eigl1,eigl2
      double precision traceT,DETT,eigt1,eigt2
 
      double precision mSN2(3,3),SNegvd(3),USNvd(3,3)
      double precision mSU3(6,6), diag(6),Vo(6,6),Wo(6,6)
      double precision stops1(2,2),usnce(3,3)
      
      double precision MW,MZ,sinsqtw,pi

      common/lfvmixing/USLvd,USUvd,USDvd,USNvd
      common/lfveigval/SLegvd,SUegvd,SDegvd,SNegvd


      common/sfmixing_susy/thetat,thetab,thetatau,thetac,thetas,thetamu,
     $     thetau,thetad,thetae

      common/sfmixing_mz/thetatz,thetabz,thetatauz,thetacz,thetasz,
     $     thetamuz,thetauz,thetadz,thetaez



C     --------------------------------------------------------------
C     Standard Inputs  and continuation of declaration statements
C     ---------------------------------------------------------------
      
      include 'stdinputs.h'
!------------------------------
      
!     redefining MT, MB, MTAU

      MT = MT_susy
      MB = MB_susy
      mTau = Mtau_susy
      
C     To be always be included ONLY after the defintion of tanbeta.
C     Here tanbeta is an input. 
      

      MW = dsqrt(MWsq)
      MZ = dsqrt(MZsq)
      
      pi = 4.d0 * datan(1.d0)
      
      sw2 = (1.d0 - (MWsq/MZsq))      
      stw = dsqrt(sw2)
      ctw = dsqrt(1.d0 - sw2)
      sinsqtw = sw2
      
      beta  = datan(tanbeta)
      sbeta = dsin(beta)
      cbeta = dcos(beta)
      cos2beta = cbeta*cbeta - sbeta*sbeta


      loopmulin: do i = 1, 3
      loopmuljn: do j = 1, 3

      AURGn(i,j) = AURG(i,j)
      
      ADRGn(i,j) = ADRG(i,j)

      AERGn(i,j) = AERG(i,j)
      
      enddo loopmuljn
      enddo loopmulin


      
      loopmuli:  do i = 1, 3
      loopmulj:  do j = 1, 3

      AURG(i,j) = AURG(i,j) * yuRG(i,j)
      
      ADRG(i,j) = ADRG(i,j) * ydRG(i,j)

      AERG(i,j) = AERG(i,j) * yeRG(i,j)
      
      enddo loopmulj
      enddo loopmuli
      
      AOK = 0




C     Mass Matrices
C     ---------------------------------------------------------------
C     Up-Squark Mass Matrix 
C     ---------------------------------------------------------------

c$$$      mSU(1,1) = mSQRG(1,1) + mUQ*mUQ + MZ*MZ*cos2beta*(0.5d0 - 
c$$$     .	         (2.d0/3.d0)*sinsqtw)
c$$$      mSU(1,2) = mSQRG(1,2) 
c$$$      mSU(1,3) = mSQRG(1,3)
c$$$      mSU(1,4) = - AURG(1,1)*vev2/dsqrt(2.d0)-
c$$$     $     mUQ*sgnmu*mur*(1.d0/tanbeta) 
c$$$      mSU(1,5) = - AURG(1,2)*vev2/dsqrt(2.d0) 
c$$$      mSU(1,6) = - AURG(1,3)*vev2/dsqrt(2.d0)
c$$$      mSU(2,1) = mSQRG(2,1)
c$$$      mSU(2,2) = mSQRG(2,2) + mC*mC + MZ*MZ*cos2beta*(0.5d0 -
c$$$     .     (2.d0/3.d0)*sinsqtw)
c$$$      mSU(2,3) = mSQRG(2,3)
c$$$      mSU(2,4) = - AURG(2,1)*vev2/dsqrt(2.d0) 
c$$$      mSU(2,5) = - AURG(2,2)*vev2/dsqrt(2.d0)-
c$$$     $     mC*sgnmu*mur*(1.d0/tanbeta)
c$$$      mSU(2,6) = - AURG(2,3)*vev2/dsqrt(2.d0) 
c$$$      mSU(3,1) = mSQRG(3,1)
c$$$      mSU(3,2) = mSQRG(3,2)
c$$$      mSU(3,3) = mSQRG(3,3) + mT*mT + MZ*MZ*cos2beta*(0.5d0 
c$$$     $     - (2.d0/3.d0)*sinsqtw)
c$$$      mSU(3,4) = - AURG(3,1)*vev2/dsqrt(2.d0) 
c$$$      mSU(3,5) = - AURG(3,2)*vev2/dsqrt(2.d0) 
c$$$      mSU(3,6) = - AURG(3,3)*vev2/dsqrt(2.d0) - 
c$$$     $     mT*sgnmu*mur*(1.d0/tanbeta)
c$$$      mSU(4,1) = - AURG(1,1)*vev2/dsqrt(2.d0) - 
c$$$     $     mUQ*sgnmu*mur*(1.d0/tanbeta)
c$$$      mSU(4,2) = - AURG(2,1)*vev2/dsqrt(2.d0)
c$$$      mSU(4,3) = - AURG(3,1)*vev2/dsqrt(2.d0)
c$$$      mSU(4,4) = mSURG(1,1)+mUQ*mUQ+(2.d0/3.d0)*MZ*MZ*cos2beta*sinsqtw
c$$$      mSU(4,5) = mSURG(1,2)
c$$$      mSU(4,6) = mSURG(1,3)
c$$$      mSU(5,1) = - AURG(1,2)*vev2/dsqrt(2.d0)
c$$$      mSU(5,2) = - AURG(2,2)*vev2/dsqrt(2.d0) - 
c$$$     $     mC*sgnmu*mur*(1.d0/tanbeta)
c$$$      mSU(5,3) = - AURG(3,2)*vev2/dsqrt(2.d0)
c$$$      mSU(5,4) = mSURG(2,1)
c$$$      mSU(5,5) = mSURG(2,2)+mC*mC+(2.d0/3.d0)*MZ*MZ*cos2beta*sinsqtw
c$$$      mSU(5,6) = mSURG(2,3)
c$$$      mSU(6,1) = - AURG(1,3)*vev2/dsqrt(2.d0)
c$$$      mSU(6,2) = - AURG(2,3)*vev2/dsqrt(2.d0)
c$$$      mSU(6,3) = - AURG(3,3)*vev2/dsqrt(2.d0) - 
c$$$     $     mT*sgnmu*mur*(1.d0/tanbeta)
c$$$      mSU(6,4) = mSURG(3,1)
c$$$      mSU(6,5) = mSURG(3,2)
c$$$      mSU(6,6) = mSURG(3,3) + mT*mT + (2.d0/3.d0)*MZ*MZ*cos2beta*sinsqtw

      mSU(1,1) = mSQRG(1,1) + mUQ*mUQ + MZ*MZ*cos2beta*(0.5d0 - 
     .	         (2.d0/3.d0)*sinsqtw)
      mSU(1,2) = mSQRG(1,2) 
      mSU(1,3) = mSQRG(1,3)
      mSU(1,4) =  AURG(1,1)*vev2/dsqrt(2.d0)-
     $     mUQ*sgnmu*mur*(1.d0/tanbeta) 
      mSU(1,5) =  AURG(1,2)*vev2/dsqrt(2.d0) 
      mSU(1,6) =  AURG(1,3)*vev2/dsqrt(2.d0)
      mSU(2,1) = mSQRG(2,1)
      mSU(2,2) = mSQRG(2,2) + mC*mC + MZ*MZ*cos2beta*(0.5d0 -
     .     (2.d0/3.d0)*sinsqtw)
      mSU(2,3) = mSQRG(2,3)
      mSU(2,4) =  AURG(2,1)*vev2/dsqrt(2.d0) 
      mSU(2,5) =  AURG(2,2)*vev2/dsqrt(2.d0)-
     $     mC*sgnmu*mur*(1.d0/tanbeta)
      mSU(2,6) =  AURG(2,3)*vev2/dsqrt(2.d0) 
      mSU(3,1) = mSQRG(3,1)
      mSU(3,2) = mSQRG(3,2)
      mSU(3,3) = mSQRG(3,3) + mT*mT + MZ*MZ*cos2beta*(0.5d0 
     $     - (2.d0/3.d0)*sinsqtw)
      mSU(3,4) =  AURG(3,1)*vev2/dsqrt(2.d0) 
      mSU(3,5) =  AURG(3,2)*vev2/dsqrt(2.d0) 
      mSU(3,6) =  AURG(3,3)*vev2/dsqrt(2.d0) - 
     $     mT*sgnmu*mur*(1.d0/tanbeta)
      mSU(4,1) =  AURG(1,1)*vev2/dsqrt(2.d0) - 
     $     mUQ*sgnmu*mur*(1.d0/tanbeta)
      mSU(4,2) =  AURG(2,1)*vev2/dsqrt(2.d0)
      mSU(4,3) =  AURG(3,1)*vev2/dsqrt(2.d0)
      mSU(4,4) = mSURG(1,1)+mUQ*mUQ+(2.d0/3.d0)*MZ*MZ*cos2beta*sinsqtw
      mSU(4,5) = mSURG(1,2)
      mSU(4,6) = mSURG(1,3)
      mSU(5,1) =  AURG(1,2)*vev2/dsqrt(2.d0)
      mSU(5,2) =  AURG(2,2)*vev2/dsqrt(2.d0) - 
     $     mC*sgnmu*mur*(1.d0/tanbeta)
      mSU(5,3) =  AURG(3,2)*vev2/dsqrt(2.d0)
      mSU(5,4) = mSURG(2,1)
      mSU(5,5) = mSURG(2,2)+mC*mC+(2.d0/3.d0)*MZ*MZ*cos2beta*sinsqtw
      mSU(5,6) = mSURG(2,3)
      mSU(6,1) =  AURG(1,3)*vev2/dsqrt(2.d0)
      mSU(6,2) =  AURG(2,3)*vev2/dsqrt(2.d0)
      mSU(6,3) =  AURG(3,3)*vev2/dsqrt(2.d0) - 
     $     mT*sgnmu*mur*(1.d0/tanbeta)
      mSU(6,4) = mSURG(3,1)
      mSU(6,5) = mSURG(3,2)
      mSU(6,6) = mSURG(3,3) + mT*mT + (2.d0/3.d0)*MZ*MZ*cos2beta*sinsqtw
      
C     Make a copy of the Mass Matrix. 
      
      do a = 1,6
         do m = 1,6 
            
           mSU1(a,m) = mSU(a,m)
           mSU2(a,m) = mSU(a,m)
           mSU3(a,m) = mSU(a,m)
           
        enddo
      enddo 


!-------------------------------------------------------------------------------      
C     find  Eigenvalues using Lapacks eigenvalue routine. 
            
      
C     find singular values and diaongalising matrices. 
      
      if(scale.gt.92.d0)then

      info = 10
      
c$$$      call dgesvd('A','A',6,6,mSU1,6,SUegvd,USUvd,6,USUT,6,
c$$$     $     work,lwork,info)
c$$$      
c$$$      if(info.eq.0) then
c$$$         AOK = AOK + 1
c$$$      endif
c$$$
c$$$C     According to definition take the transpose for the diagonalising matrix.
c$$$      
c$$$      do a = 1, 6
c$$$ 	 do m = 1, 6 
c$$$            
c$$$            USUvd(a,m) = USUT(a,m)
c$$$            
c$$$ 	 enddo  
c$$$      enddo 
c$$$
c$$$
c$$$      print*,'SUegvd (1) ', SUegvd(1)
c$$$      print*,'SUegvd (2) ', SUegvd(2)
c$$$      print*,'SUegvd (3) ', SUegvd(3)
c$$$      print*,'SUegvd (4) ', SUegvd(4)
c$$$      print*,'SUegvd (5) ', SUegvd(5)
c$$$      print*,'SUegvd (6) ', SUegvd(6)


      call SVD(6, 6, mSU1,6, SUegvd, USUvd,6, USUT,6, -1)
!      call CEigensystem(6, mSU3,6, diag, Vo,6, 0)
!      call HEigensystem(6, mSU3,6, diag, Vo,6, 0)

c$$$      print*,"daig(1) ", (diag(1))
c$$$      print*,"daig(2) ", (diag(2))
c$$$      print*,"daig(3) ", (diag(3))
c$$$      print*,"daig(4) ", (diag(4))
c$$$      print*,"daig(5) ", (diag(5))
c$$$      print*,"daig(6) ", (diag(6))
c$$$      do a = 1, 6
c$$$ 	 do m = 1, 6 
c$$$            
c$$$            USUvd(a,m) = USUT(a,m)
c$$$            print*,"USUVD (",a,m,") =",USUvd(a,m)
c$$$            print*,"Vo (",a,m,")=",Vo(a,m)
c$$$ 	 enddo  
c$$$      enddo 

      
      endif


!---------------------------------
C     squark sector, 3rd gen

      stops(1,1) = mSU2(3,3)
      stops(1,2) = mSU2(3,6)
      stops(2,1) = mSU2(6,3)
      stops(2,2) = mSU2(6,6)

c$$$  stops1(1,1) = stops(1,1)
c$$$  stops1(1,2) = stops(1,2)
c$$$  stops1(2,1) = stops(2,1)
c$$$  stops1(2,2) = stops(2,2)
c$$$  
c$$$  call DGEEVX( 'N', 'V', 'V', 'N', 2, stops, 2, stopegg, Stopsi,
c$$$  $                   stopL, 2, stopR, 2, ILO, IHI, SCALE1, ABNRM,
c$$$  $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
c$$$  
c$$$  print*,"stops ", stopegg(1),stopegg(2)
c$$$  print*,"stopl ",stopl(1,1),stopl(1,2),stopl(2,1),stopl(2,2)
      
      
c$$$      if(scale.gt.94.d0)then
c$$$         
c$$$c$$$         print*,"thetat d = ", datan((2.d0*stops(1,2))/(stops(1,1)
c$$$c$$$     $        - stops(2,2)))/2.d0
c$$$      else
c$$$         print*,"thetatz d = ", datan((2.d0*stops(1,2))/(stops(1,1)
c$$$     $        -stops(2,2)))/2.d0
c$$$      endif


      Call CEigensystem(2,stops,2,stopegg,stopl,2,0)
c$$$  
c$$$  print*,"stops-ce ", stopegg(1),stopegg(2)
c$$$  print*,"stopl ",stopl(1,1),stopl(1,2),stopl(2,1),stopl(2,2)

!---  2nd gen

         schs(1,1) = mSU2(2,2)
         schs(1,2) = mSU2(2,5)
         schs(2,1) = mSU2(5,2)
         schs(2,2) = mSU2(5,5)

c$$$
c$$$      call DGEEVX( 'N', 'V', 'V', 'N', 2, schs, 2, schegg, Schsi,
c$$$     $                   schL, 2, schR, 2, ILO, IHI, SCALE1, ABNRM,
c$$$     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )

      Call CEigensystem(2,schs,2,schegg,schl,2,0)

!---1st gen

         sups(1,1) = mSU2(1,1)
         sups(1,2) = mSU2(1,4)
         sups(2,1) = mSU2(4,1)
         sups(2,2) = mSU2(4,4)
         
c$$$
c$$$      call DGEEVX( 'N', 'V', 'V', 'N', 2, sups, 2, supegg, Supsi,
c$$$     $                   supL, 2, supR, 2, ILO, IHI, SCALE1, ABNRM,
c$$$     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )

      Call CEigensystem(2,sups,2,supegg,supl,2,0)
!--------------------------------------------------------------------------

       if(scale.gt.92.d0)then
!          print*,"stopL11, stopL12 = ",stopL(1,1),stopL(1,2) 
          thetat = datan(stopL(1,2)/stopL(1,1)) !dacos(dabs(stopL(1,1)))
          thetac = datan(schL(1,2)/schL(1,1)) !dacos(dabs(schL(1,1)))
          thetau = datan(supL(1,2)/supL(1,1)) !dacos(dabs(supL(1,1)))
       else
!          print*,"stopL11z, stopL12z = ",stopL(1,1),stopL(1,2) 
          thetatz = datan(stopL(1,2)/stopL(1,1)) !dacos(dabs(stopL(1,1)))
          thetacz = datan(schL(1,2)/schL(1,1)) !dacos(dabs(schL(1,1)))
          thetauz = datan(supL(1,2)/supL(1,1)) !dacos(dabs(supL(1,1)))
       endif



c$$$      if(dabs(stopr(1,1)).gt.dabs(stopr(2,1)))then
c$$$         SUegg(1) = stopegg(1)
c$$$         SUegg(2) = stopegg(2)
c$$$      else
c$$$         SUegg(1) = stopegg(2)
c$$$         SUegg(2) = stopegg(1)
c$$$      endif
c$$$      
c$$$
c$$$      if(dabs(schr(1,1)).gt.dabs(schr(2,1)))then
c$$$         SUegg(3) = schegg(1)
c$$$         SUegg(4) = schegg(2)
c$$$      else
c$$$         SUegg(3) = schegg(2)
c$$$         SUegg(4) = schegg(1)
c$$$      endif
c$$$
c$$$
c$$$
c$$$      if(dabs(supr(1,1)).gt.dabs(supr(2,1)))then
c$$$         SUegg(5) = supegg(1)
c$$$         SUegg(6) = supegg(2)
c$$$      else
c$$$         SUegg(5) = supegg(2)
c$$$         SUegg(6) = supegg(1)
c$$$      endif


c$$$      traceT = MSU2(6,6) + MSU2(3,3)
c$$$      DETT = MSU2(6,6)* MSU2(3,3) - MSU2(3,6)**2.D0
c$$$
c$$$      eigt1 = (traceT + dsqrt(traceT**2.d0 - 4.d0*DETT ))/2.d0
c$$$  
c$$$      eigt2 = (traceT - dsqrt(traceT**2.d0 - 4.d0*DETT ))/2.d0

c$$$      SUegg(1) = stopegg(1)
c$$$      SUegg(2) = stopegg(2)
c$$$      SUegg(3) = schegg(1)
c$$$      SUegg(4) = schegg(2)
c$$$      SUegg(5) = supegg(1)
c$$$      SUegg(6) = supegg(2)


      SUegg(1) = stopegg(2)
      SUegg(2) = stopegg(1)
      SUegg(3) = schegg(2)
      SUegg(4) = schegg(1)
      SUegg(5) = supegg(2)
      SUegg(6) = supegg(1)
      
C     ---------------------------------------------------------------
C     Down-Squark Mass Matrix 
C     ---------------------------------------------------------------

      mSD(1,1) = mSQRG(1,1) + mD*mD + MZ*MZ*cos2beta*(- 0.5d0 + 
     .	          (1.d0/3.d0)*sinsqtw) 
      mSD(1,2) = mSQRG(1,2) 
      mSD(1,3) = mSQRG(1,3)
      mSD(1,4) = ADRG(1,1)*vev1/dsqrt(2.d0) - mD*sgnmu*mur*tanbeta 
      mSD(1,5) = ADRG(1,2)*vev1/dsqrt(2.d0) 
      mSD(1,6) = ADRG(1,3)*vev1/dsqrt(2.d0)
      mSD(2,1) = mSQRG(2,1)
      mSD(2,2) = mSQRG(2,2) + mS*mS + MZ*MZ*cos2beta*(- 0.5d0 + 
     . 	          (1.d0/3.d0)*sinsqtw)
      mSD(2,3) = mSQRG(2,3)
      mSD(2,4) = ADRG(2,1)*vev1/dsqrt(2.d0) 
      mSD(2,5) = ADRG(2,2)*vev1/dsqrt(2.d0) - mS*sgnmu*mur*tanbeta
      mSD(2,6) = ADRG(2,3)*vev1/dsqrt(2.d0) 
      mSD(3,1) = mSQRG(3,1)
      mSD(3,2) = mSQRG(3,2)
      mSD(3,3) = mSQRG(3,3) + mB*mB + MZ*MZ*cos2beta*(- 0.5d0 + 
     .            (1.d0/3.d0)*sinsqtw)
      mSD(3,4) = ADRG(3,1)*vev1/dsqrt(2.d0) 
      mSD(3,5) = ADRG(3,2)*vev1/dsqrt(2.d0) 
      mSD(3,6) = ADRG(3,3)*vev1/dsqrt(2.d0) - mB*sgnmu*mur*tanbeta
      mSD(4,1) = ADRG(1,1)*vev1/dsqrt(2.d0) - mD*sgnmu*mur*tanbeta
      mSD(4,2) = ADRG(2,1)*vev1/dsqrt(2.d0)
      mSD(4,3) = ADRG(3,1)*vev1/dsqrt(2.d0)
      mSD(4,4) = mSDRG(1,1) + mD*mD - (1.d0/3.d0)*MZ*MZ*
     .                  cos2beta*sinsqtw
      mSD(4,5) = mSDRG(1,2)
      mSD(4,6) = mSDRG(1,3)
      mSD(5,1) = ADRG(1,2)*vev1/dsqrt(2.d0)
      mSD(5,2) = ADRG(2,2)*vev1/dsqrt(2.d0) - mS*sgnmu*mur*tanbeta
      mSD(5,3) = ADRG(3,2)*vev1/dsqrt(2.d0)
      mSD(5,4) = mSDRG(2,1)
      mSD(5,5) = mSDRG(2,2) + mS*mS - (1.d0/3.d0)*MZ*MZ*
     .            cos2beta*sinsqtw
      mSD(5,6) = mSDRG(2,3)
      mSD(6,1) = ADRG(1,3)*vev1/dsqrt(2.d0)
      mSD(6,2) = ADRG(2,3)*vev1/dsqrt(2.d0)
      mSD(6,3) = ADRG(3,3)*vev1/dsqrt(2.d0) - mB*sgnmu*mur*tanbeta
      mSD(6,4) = mSDRG(3,1)
      mSD(6,5) = mSDRG(3,2)
      mSD(6,6) = mSDRG(3,3) + mB*mB - (1.d0/3.d0)*MZ*MZ*
     .            cos2beta*sinsqtw


C     --------------------------------------------------------
C     Make a copy of the Mass Matrix. 

      do a = 1, 6
         do m = 1, 6 
            

            mSD1(a,m) = mSD(a,m)
            mSD2(a,m) = mSD(a,m)
            mSD3(a,m) = mSD(a,m)

         enddo
      enddo 

C     find sdown sector Eigenvalues using Lapacks eigenvalue routine. 
      

C     find singular values and diaongalising matrices. 
      
      info = 10
      
        if(scale.gt.92.d0)then

c$$$       call dgesvd('A','A',6,6,mSD3,6,SDegvd,USDvd,6,USDT,6,
c$$$     .	work,lwork,info)
c$$$
c$$$
c$$$C     According to definition take the transpose for the diagonalising matrix.
c$$$      
c$$$      do a = 1, 6
c$$$ 	 do m = 1, 6 
c$$$            
c$$$            USDvd(a,m) = USDT(a,m)
c$$$            
c$$$ 	 enddo  
c$$$      enddo 

        call SVD(6, 6, mSD3,6, SDegvd, USDvd,6, USDT,6, -1)

      endif

!--------------------
        
        

         sbtms(1,1) = mSD2(3,3)
         sbtms(1,2) = mSD2(3,6)
         sbtms(2,1) = mSD2(6,3)
         sbtms(2,2) = mSD2(6,6)


c$$$
c$$$      call DGEEVX( 'N', 'V', 'V', 'N', 2, sbtms, 2, sbtmegg, Sbtmsi,
c$$$     $                   sbtmL, 2, sbtmR, 2, ILO, IHI, SCALE1, ABNRM,
c$$$     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )

      Call CEigensystem(2,sbtms,2,sbtmegg,sbtml,2,0)

!---2nd gen

         sstrs(1,1) = mSD2(2,2)
         sstrs(1,2) = mSD2(2,5)
         sstrs(2,1) = mSD2(5,2)
         sstrs(2,2) = mSD2(5,5)


c$$$      call DGEEVX( 'N', 'V', 'V', 'N', 2, sstrs, 2, sstregg, Sstrsi,
c$$$     $                   sstrL, 2, sstrR, 2, ILO, IHI, SCALE1, ABNRM,
c$$$     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )

      Call CEigensystem(2,sstrs,2,sstregg,sstrl,2,0)

!---1st gen

         sdowns(1,1) = mSD2(1,1)
         sdowns(1,2) = mSD2(1,4)
         sdowns(2,1) = mSD2(4,1)
         sdowns(2,2) = mSD2(4,4)
         
c$$$
c$$$      call DGEEVX( 'N', 'V', 'V', 'N', 2, sdowns, 2, sdownegg, Sdownsi,
c$$$     $                   sdownL, 2, sdownR, 2, ILO, IHI, SCALE1, ABNRM,
c$$$     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
c$$$

      Call CEigensystem(2,sdowns,2,sdownegg,sdownl,2,0)
!---------------------------------------------------------------------
      
      if(scale.gt.92.d0)then
         
         thetab = datan(sbtmL(1,2)/sbtmL(1,1)) !dacos((sbtmL(1,1)))
         thetas = datan(sstrL(1,2)/sstrL(1,1)) !dacos((sstrL(1,1)))
         thetad = datan(sdownL(1,2)/sdownL(1,1)) !dacos((sdownL(1,1)))
      else
         thetabz = datan(sbtmL(1,2)/sbtmL(1,1)) !dacos((sbtmL(1,1)))
         thetasz = datan(sstrL(1,2)/sstrL(1,1)) !dacos((sstrL(1,1)))
         thetadz = datan(sdownL(1,2)/sdownL(1,1)) !dacos((sdownL(1,1)))
      endif
         
c$$$  thetabz = dacos((sbtmL(1,1)))
c$$$  thetasz = dacos((sstrL(1,1)))
c$$$  thetadz = dacos((sdownL(1,1)))
c$$$  endif
      
c$$$   
c$$$
c$$$
c$$$      if(dabs(sbtmr(1,1)).gt.dabs(sbtmr(2,1)))then
c$$$         SDegg(1) = sbtmegg(1)
c$$$         SDegg(2) = sbtmegg(2)
c$$$      else
c$$$         SDegg(1) = sbtmegg(2)
c$$$         SDegg(2) = sbtmegg(1)
c$$$      endif
c$$$      
c$$$
c$$$      if(dabs(sstrr(1,1)).gt.dabs(sstrr(2,1)))then
c$$$         SDegg(3) = sstregg(1)
c$$$         SDegg(4) = sstregg(2)
c$$$      else
c$$$         SDegg(3) = sstregg(2)
c$$$         SDegg(4) = sstregg(1)
c$$$      endif
c$$$
c$$$
c$$$      if(dabs(sdownr(1,1)).gt.dabs(sdownr(2,1)))then
c$$$         SDegg(5) = sdownegg(1)
c$$$         SDegg(6) = sdownegg(2)
c$$$      else
c$$$         SDegg(5) = sdownegg(2)
c$$$         SDegg(6) = sdownegg(1)
c$$$      endif
c$$$
c$$$
c$$$
c$$$      traceB = MSD2(6,6) + MSD2(3,3)
c$$$      DETB = MSD2(6,6)* MSD2(3,3) - MSD2(3,6)**2.D0
c$$$
c$$$      eigb1 = (traceb + dsqrt(traceb**2.d0 - 4.d0*DETB ))/2.d0
c$$$  
c$$$      eigb2 = (traceb - dsqrt(traceb**2.d0 - 4.d0*DETB ))/2.d0

c$$$      SDegg(1) = sbtmegg(1)
c$$$      SDegg(2) = sbtmegg(2)
c$$$      SDegg(3) = sstregg(1)
c$$$      SDegg(4) = sstregg(2)
c$$$      SDegg(5) = sdownegg(1)
c$$$      SDegg(6) = sdownegg(2)


      SDegg(1) = sbtmegg(2)
      SDegg(2) = sbtmegg(1)
      SDegg(3) = sstregg(2)
      SDegg(4) = sstregg(1)
      SDegg(5) = sdownegg(2)
      SDegg(6) = sdownegg(1)


C     --------------------------------------------------------
C     Slepton Mass Matrix
C     --------------------------------------------------------


      mSL(1,1) = mSLRG(1,1)+mE*mE+MZ*MZ*cos2beta*(- 0.5d0 + sinsqtw) 
      mSL(1,2) = mSLRG(1,2) 
      mSL(1,3) = mSLRG(1,3)
      mSL(1,4) = AERG(1,1)*vev1/dsqrt(2.d0) - mE*sgnmu*mur*tanbeta 
      mSL(1,5) = AERG(1,2)*vev1/dsqrt(2.d0) 
      mSL(1,6) = AERG(1,3)*vev1/dsqrt(2.d0)
      mSL(2,1) = mSLRG(2,1)
      mSL(2,2) = mSLRG(2,2)+mMU*mMU+MZ*MZ*cos2beta*(-0.5d0 + sinsqtw)
      mSL(2,3) = mSLRG(2,3)
      mSL(2,4) = AERG(2,1)*vev1/dsqrt(2.d0) 
      mSL(2,5) = AERG(2,2)*vev1/dsqrt(2.d0) - mMU*sgnmu*mur*tanbeta
      mSL(2,6) = AERG(2,3)*vev1/dsqrt(2.d0) 
      mSL(3,1) = mSLRG(3,1)
      mSL(3,2) = mSLRG(3,2)
      mSL(3,3) = mSLRG(3,3)+mTAU*mTAU+MZ*MZ*cos2beta*(-0.5d0+sinsqtw)
      mSL(3,4) = AERG(3,1)*vev1/dsqrt(2.d0) 
      mSL(3,5) = AERG(3,2)*vev1/dsqrt(2.d0) 
      mSL(3,6) = AERG(3,3)*vev1/dsqrt(2.d0) - mTAU*sgnmu*mur*tanbeta
      mSL(4,1) = AERG(1,1)*vev1/dsqrt(2.d0) - mE*sgnmu*mur*tanbeta
      mSL(4,2) = AERG(2,1)*vev1/dsqrt(2.d0)
      mSL(4,3) = AERG(3,1)*vev1/dsqrt(2.d0)
      mSL(4,4) = mSERG(1,1) + mE*mE - MZ*MZ*cos2beta*sinsqtw
      mSL(4,5) = mSERG(1,2)
      mSL(4,6) = mSERG(1,3)
      mSL(5,1) = AERG(1,2)*vev1/dsqrt(2.d0)
      mSL(5,2) = AERG(2,2)*vev1/dsqrt(2.d0) - mMU*sgnmu*mur*tanbeta
      mSL(5,3) = AERG(3,2)*vev1/dsqrt(2.d0)
      mSL(5,4) = mSERG(2,1)
      mSL(5,5) = mSERG(2,2) + mMU*mMU - MZ*MZ*cos2beta*sinsqtw
      mSL(5,6) = mSERG(2,3) 
      mSL(6,1) = AERG(1,3)*vev1/dsqrt(2.d0)
      mSL(6,2) = AERG(2,3)*vev1/dsqrt(2.d0)
      mSL(6,3) = AERG(3,3)*vev1/dsqrt(2.d0) - mTAU*sgnmu*mur*tanbeta
      mSL(6,4) = mSERG(3,1)
      mSL(6,5) = mSERG(3,2) 
      mSL(6,6) = mSERG(3,3) + mTAU*mTAU - MZ*MZ*cos2beta*sinsqtw


C	-------------------------------------------------------------

C     Make a copy of the Mass Matrix. 

      do a = 1,6
         do m = 1,6 
            
            mSL1(a,m) = mSL(a,m)
            mSL2(a,m) = mSL(a,m)
            mSL3(a,m) = mSL(a,m)

         enddo
      enddo 



C     find singular values and diaongalising matrices. 
      
      info = 10
      

C----------------------------
C     Lepton Sector 6X6
C----------------------------

      by66:  if(scale.gt.92.d0)then
         
         call SVD(6, 6, mSL1,6, SLegvd, USLvd,6, USLT,6, -1)
         
         
         call SVD(6, 6, mSL3,6, SLegs, USLs,6, USLTs,6, 0)
         
      endif by66
      
c$$$      call dgesvd('A','A',6,6,mSL1,6,SLegvd,USLvd,6,USLT,6,
c$$$     $     work,lwork,info)
c$$$
c$$$!------------------------------------
c$$$C     According to definition take the transpose for the diagonalising matrix.
c$$$
c$$$      do a = 1, 6
c$$$ 	 do m = 1,6 
c$$$            
c$$$            USLvd(a,m) = USLT(a,m)
c$$$            
c$$$ 	 enddo
c$$$      enddo 

c$$$      call DGEEVX( 'N', 'V', 'V', 'N', 6, mSL3, 6, slegs, slegi,
c$$$     $                   USLs, 6, USLTs, 6, ILO, IHI, SCALE1, ABNRM,
c$$$     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
c$$$      
c$$$
c$$$      if(info.eq.0) then
c$$$         AOK = AOK + 1
c$$$      endif



!-----2x2 slepton

         staus(1,1) = mSL2(3,3)
         staus(1,2) = mSL2(3,6)
         staus(2,1) = mSL2(6,3)
         staus(2,2) = mSL2(6,6)

c$$$         if(scale.gt.94.d0)then
c$$$            
c$$$            print*,"thetatau d = ", datan((2.d0*staus(1,2))/(staus(1,1)
c$$$     $           - staus(2,2)))/2.d0
c$$$         else
c$$$            print*,"thetatauz d = ", datan((2.d0*staus(1,2))/(staus(1,1)
c$$$     $           -staus(2,2)))/2.d0
c$$$         endif
         


C     find SLepton Eigenvalues using Lapacks eigenvalue routine. 
c$$$
c$$$      call DGEEVX( 'N', 'V', 'V', 'N', 2, staus, 2, staugg, Stausi,
c$$$     $                   stauL, 2, stauR, 2, ILO, IHI, SCALE1, ABNRM,
c$$$     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
c$$$

      Call CEigensystem(2,staus,2,staugg,staul,2,0)

         smus(1,1) = mSL2(2,2)
         smus(1,2) = mSL2(2,5)
         smus(2,1) = mSL2(5,2)
         smus(2,2) = mSL2(5,5)


C     find SLepton Eigenvalues using Lapacks eigenvalue routine. 


c$$$      call DGEEVX( 'N', 'V', 'V', 'N', 2, smus, 2, smugg, Smusi,
c$$$     $                   smuL, 2, smuR, 2, ILO, IHI, SCALE1, ABNRM,
c$$$     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
c$$$


      Call CEigensystem(2,smus,2,smugg,smul,2,0)

         sels(1,1) = mSL2(1,1)
         sels(1,2) = mSL2(1,4)
         sels(2,1) = mSL2(4,1)
         sels(2,2) = mSL2(4,4)


C     find SLepton Eigenvalues using Lapacks eigenvalue routine. 


c$$$      call DGEEVX( 'N', 'V', 'V', 'N', 2, sels, 2, selgg, Selsi,
c$$$     $                   selL, 2, selR, 2, ILO, IHI, SCALE1, ABNRM,
c$$$     $                   RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )

      Call CEigensystem(2,sels,2,selgg,sell,2,0)

      if(scale.gt.94.d0)then
         
!     print*,"thetatau d = ", datan((2.d0*sels(1,2))/(sels(1,1)
!     $         -sels(2,2)))/2.d0
         
         thetatau = datan(stauL(1,2)/stauL(1,1)) !dacos(dabs(stauL(1,1)))
         thetamu = datan(smuL(1,2)/smuL(1,1)) !dacos(dabs(smuL(1,1)))
         thetae = datan(selL(1,2)/selL(1,1)) !dacos(dabs(selL(1,1)))

!         print*,"thetatau = ", thetatau

      else
         
!     print*,"thetatauz d = ", datan((2.d0*sels(1,2))/(sels(1,1)
!     $         -sels(2,2)))/2.d0
         
         thetatauz = datan(stauL(1,2)/stauL(1,1)) !dacos(dabs(stauL(1,1)))
         thetamuz = datan(smuL(1,2)/smuL(1,1)) !dacos(dabs(smuL(1,1)))
         thetaez = datan(selL(1,2)/selL(1,1)) !dacos(dabs(selL(1,1)))
         
c$$$  thetatauz = dacos(dabs(stauL(1,1)))
c$$$  thetamuz = dacos(dabs(smuL(1,1)))
c$$$  thetaez = dacos(dabs(selL(1,1)))

!         print*,"thetatauz = ", thetatau

      endif

c$$$
c$$$      if(dabs(staur(1,1)).gt.dabs(staur(2,1)))then
c$$$         SLegg(1) = staugg(1)
c$$$         SLegg(2) = staugg(2)
c$$$      else
c$$$         SLegg(1) = staugg(2)
c$$$         SLegg(2) = staugg(1)
c$$$      endif
c$$$      
c$$$
c$$$      if(dabs(smur(1,1)).gt.dabs(smur(2,1)))then
c$$$         SLegg(3) = smugg(1)
c$$$         SLegg(4) = smugg(2)
c$$$      else
c$$$         SLegg(3) = smugg(2)
c$$$         SLegg(4) = smugg(1)
c$$$      endif
c$$$
c$$$
c$$$      if(dabs(selr(1,1)).gt.dabs(selr(2,1)))then
c$$$         SLegg(5) = selgg(1)
c$$$         SLegg(6) = selgg(2)
c$$$      else
c$$$         SLegg(5) = selgg(2)
c$$$         SLegg(6) = selgg(1)
c$$$      endif
c$$$
c$$$
c$$$      traceL = MSL2(6,6) + MSL2(3,3)
c$$$      DETL = MSL2(6,6)* MSL2(3,3) - MSL2(3,6)**2.D0
c$$$
c$$$      eigl1 = (traceL + dsqrt(traceL**2.d0 - 4.d0*DETL ))/2.d0
c$$$  
c$$$      eigl2 = (traceL - dsqrt(traceL**2.d0 - 4.d0*DETL ))/2.d0

c$$$      SLegg(1) = staugg(1)
c$$$      SLegg(2) = staugg(2)
c$$$      SLegg(3) = smugg(1)
c$$$      SLegg(4) = smugg(2)
c$$$      SLegg(5) = selgg(1)
c$$$      SLegg(6) = selgg(2)

      SLegg(1) = staugg(2)
      SLegg(2) = staugg(1)
      SLegg(3) = smugg(2)
      SLegg(4) = smugg(1)
      SLegg(5) = selgg(2)
      SLegg(6) = selgg(1)


!---------------------------------
C     ------------------------------------------------------------------	
C     Sneutrino Mass Matrix
C     -------------------------------------------------------------------

      mSN(1,1) = mSLRG(1,1) + 0.5d0*MZ*MZ*cos2beta
      mSN(1,2) = mSLRG(1,2)
      mSN(1,3) = mSLRG(1,3)
      mSN(2,1) = mSLRG(2,1)
      mSN(2,2) = mSLRG(2,2) + 0.5d0*MZ*MZ*cos2beta
      mSN(2,3) = mSLRG(2,3)
      mSN(3,1) = mSLRG(3,1)
      mSN(3,2) = mSLRG(3,2)
      mSN(3,3) = mSLRG(3,3) + 0.5d0*MZ*MZ*cos2beta

C     Make a copy of the Mass Matrix. 

      do a = 1,3
         do m = 1,3 
            
            mSN1(a,m) = mSN(a,m)
            mSN2(a,m) = mSN(a,m)

         enddo
      enddo 

c$$$ 100  format(/A/,("Eign(",I1,") = ", 1x, 1pe11.4))
c$$$      
c$$$!      print* ,"q,g,gp = ",q,g,gp
c$$$      print 100,"mSN: ", (i, mSN(1,i), i = 1, 3)
c$$$      print 100,"mSLRG: ", (i, mSLRG(1,i), i = 1, 6)
c$$$!      print 100,"SLeggz: ", (i, SLeggz(i), i = 1, 6)
c$$$!      print 100,"SNeggz: ", (i, SNeggz(i), i = 1, 3)
c$$$!      print*,"==========================="


      by33:  if(scale.gt.92.d0)then
         
c$$$         call dgesvd('A','A',3,3,mSN2,3,SNegvd,USNvd,3,USNT,3,
c$$$     $        work,lwork,info)
c$$$         
         call SVD(3, 3, mSN2,3, SNegvd, USNvd,3, USNT,3, -1)

!------------------------------------
C     According to definition take the transpose for the diagonalising matrix.
         
c$$$         do a = 1, 3
c$$$            do m = 1,3 
c$$$               
c$$$               USNvd(a,m) = USNT(a,m)
c$$$               
c$$$            enddo
c$$$         enddo 

      endif by33
      
      info = 10

C     find sneutrino eigenvalues 
c$$$
c$$$      call dsyev('V','U',3,mSN,3,SNegg,work1,lwork1,info)

      Call CEigensystem(3,mSN,3,SNeggt,usnce,3,0)

      if(info.eq.0)then
          AOK = AOK + 1
      endif 

C     find the singular values and the diagonalising matrix.

      info = 10
      
c$$$      call dgesvd('A','A',3,3,mSN1,3,SNeg,USN,3,USNT,3,
c$$$     .     work1,lwork1,info)

        call SVD(3, 3, mSN1,3, SNeg, USN,3, USNT,3, -1)
      
c$$$      if(info.eq.0) then
c$$$         AOK = AOK + 1
c$$$      endif
c$$$
c$$$C     Take the transpose as per the defintion.
c$$$
c$$$      do a = 1,3
c$$$         do m = 1,3 
c$$$            
c$$$            USN(a,m) = USNT(a,m)
c$$$
c$$$         enddo
c$$$      enddo

c$$$      print*,"Sneg(1),Sneg(2),Sneg(3)",Sneg(1),
c$$$     $     Sneg(2),Sneg(3)

        SNegg(1) = SNeggt(3) 
        SNegg(2) = SNeggt(2) 
        SNegg(3) = SNeggt(1) 


C     --------------------------------------------------------
C     Neutralino Mass Matrix 
C     ---------------------------------------------------------

      MNeut(1,1) = M1tz
      MNeut(1,2) = 0.d0
      MNeut(1,3) = - (MZ*stw*cbeta)
      MNeut(1,4) = MZ*stw*sbeta
      MNeut(2,1) = 0.d0
      MNeut(2,2) = M2tz
      MNeut(2,3) = MZ*ctw*cbeta
      MNeut(2,4) = -(MZ*ctw*sbeta)
      MNeut(3,1) = -(MZ*stw*cbeta)
      MNeut(3,2) =  MZ*ctw*cbeta
      MNeut(3,3) = 0.d0
      MNeut(3,4) = - sgnmu*mur
      MNeut(4,1) = MZ*stw*sbeta
      MNeut(4,2) = - (MZ*ctw*sbeta)
      MNeut(4,3) = - sgnmu*mur 
      MNeut(4,4) = 0.d0 


!-------------------------------------------------------------------

C     Make a copy of the Mass Matrix. 

      do a = 1, 4
         do m = 1, 4 
            
            MNeut1(a,m) = MNeut(a,m)
            MNeut2(a,m) = MNeut(a,m)

         enddo
      enddo 

      Call CEigensystem(4,MNeut1,4,Neg,ON,4,1)
      
C     find singular values and diagonalising matrix.

c$$$      call DGEEV( 'V', 'V', 4, MNeut2, 4, Neg, Neuevi, ONL, 4, 
c$$$     $     ONR, 4, WORK, LWORK, INFO )
c$$$      
c$$$      if(info.eq.0) then
c$$$         AOK = AOK + 1
c$$$      endif
c$$$      
c$$$      do a = 1,4
c$$$         print*,"Neg (",a,")", Neg(a)
c$$$         do m = 1,4 
c$$$            
c$$$            ON(a,m) = ONR(m,a)
c$$$            print*,"ON(",a,m,")", ON(a,m)
c$$$         enddo
c$$$      enddo

 
c$$$ !     call SVD(4, 4, Mneut1,4, Neg, ON,4, ONR,4, 1)
c$$$
c$$$      do a = 1,4
c$$$         print*,"Neg-CE (",a,")", Neg(a)
c$$$         do m = 1,4 
c$$$            
c$$$ !           ON(a,m) = ONR(m,a)
c$$$            print*,"ON-CE(",a,m,")", ON(a,m)
c$$$         enddo
c$$$      enddo

C     ----------------------------------------------------------------
C     Chargino Mass Matrix
C     ----------------------------------------------------------------


      MChar(1,1) = M2tz  
      MChar(1,2) = sqrt(2.d0)*MW*sbeta
      MChar(2,1) = sqrt(2.d0)*MW*cbeta
      MChar(2,2) = sgnmu*mur

      call SVD(2, 2, MChar,2, Ceg, OCL,2, OCR,2, 1)
      

c$$$ 10   format(1x,A,1x,1pe11.4,1x,1pe11.4)
c$$$
c$$$      print 10,"mchar11, mchar12 = ",MChar(1,1),Mchar(1,2)
c$$$      print 10,"mchar21, mchar22 = ",MChar(2,1),Mchar(2,2)
c$$$      
c$$$      call dgesvd('A','A',2,2,MChar,2,Ceg,OCL,2,OCRT,2,
c$$$     $     work1,lwork1,info)
c$$$
c$$$      
c$$$      print 10,"ceg1, ceg2 = ",Ceg(1),Ceg(2)
c$$$
c$$$      do a = 1, 2
c$$$         do m = 1,2 
c$$$            OCR(a,m) = OCRT(a,m)
c$$$            OCL(a,m) = OCL(m,a)
c$$$!            OCR(m,a) = OCRT(a,m)
c$$$            print*,"OCR(",a,m,")",OCR(a,m)
c$$$            print*,"OCL(",a,m,")",OCL(a,m)
c$$$
c$$$         enddo
c$$$      enddo 

c$$$      print 10,"ceg1, ceg2, svd = ",Ceg(1),Ceg(2)
c$$$      do a = 1, 2
c$$$         do m = 1,2 
c$$$        !    OCR(a,m) = OCRT(a,m)
c$$$         !   OCL(a,m) = OCL(m,a)
c$$$!            OCR(m,a) = OCRT(a,m)
c$$$            print*,"OCR(",a,m,")",OCR(a,m)
c$$$            print*,"OCL(",a,m,")",OCL(a,m)
c$$$
c$$$         enddo
c$$$      enddo 


C-----------------------------------------------------------------
C     Tree level Higgs mass spectrum 
C-------------------------------------------------------------------

      if(bmur.lt.0.d0) then

         mA0sq =   (dsqrt(dabs(SUegg(1)))*
     $        dsqrt(dabs(SUegg(2))) )/(16.d0 * pi *pi) 

         flag_bmu = 'BMUNEG'

!         print*,"flag -ve Bmu "

      else

         mA0sq = bmur/(sbeta*cbeta)

         flag_bmu = ' AOK'

!         print*,"flag AOK "

      endif

!----------
      
!      print*,"mz in higgs = ", mz

      mhpmsq = mA0sq + MW*MW


      Mhiggstree(1,1) = mA0sq*dsin(beta)**2.d0 + 
     $     (MZ*MZ)*dcos(beta)**2.d0

      Mhiggstree(1,2) = -1.d0*(mA0sq + 
     $     (MZ*MZ)) * dsin(beta)*dcos(beta)
      
      Mhiggstree(2,1) = Mhiggstree(1,2)

      Mhiggstree(2,2) = mA0sq*dcos(beta)**2.d0 + 
     $     (MZ*MZ)*dsin(beta)**2.d0

      Mhiggsn(1,1) = Mhiggstree(1,1) 

      Mhiggsn(1,2) = Mhiggstree(1,2) 

      Mhiggsn(2,1) = Mhiggstree(2,1) 

      Mhiggsn(2,2) = Mhiggstree(2,2) 


      Call CEigensystem(2,Mhiggsn,2,heign,hl,2,-1)

      mh0sq = heign(2)       
      mhu0sq = heign(1)      
      
      alphatree = datan(hl(1,2)/hl(1,1))

c$$$      print*,"mhtree = ", dsqrt(mh0sq)
c$$$      print*,"mhu0tree = ", dsqrt(mhu0sq)
      
c$$$        alphatree    = 0.d0 + 
c$$$     $       datan((2.d0*Mhiggsn(1,2))/(Mhiggsn(1,1)-Mhiggsn(2,2)))
c$$$     $       *0.5d0
c$$$      call DGEEV('V', 'V', 2, Mhiggsn, 2, heign, lhim, hl, 2, 
c$$$     $     hr, 2, WORK, LWORK, INFO)


!----------------

      loopmulin1: do i = 1, 3
      loopmuljn1: do j = 1, 3

      AURG(i,j) = AURGn(i,j)
      
      ADRG(i,j) = ADRGn(i,j)

      AERG(i,j) = AERGn(i,j)
      
      enddo loopmuljn1
      enddo loopmulin1

      return 

      end subroutine softspectrum 

!=======================================================================
!    Computes tree level \mu and b\mu. 
!-----------------------------------------------------------------------

      SUBROUTINE mutreelevel(tanbetaQ,mh1mz,mh2mz,mursq,mur,bmur,
     $     flagdet,flagDflat)

      double precision  cbeta, sbeta,beta
      double precision sin2beta, tan2beta
      double precision  cos2beta, tanbetaQ
      DOUBLE PRECISION mh1mz,mh2mz,mur,mursq,bmur,mh1sq,mh2sq
      integer flagdet, flagDflat
      double precision detisneg,nodflat
      DOUBLE PRECISION m0, m12, m10, m20, sgnmu, tanbeta, a0
      double precision mbpole, mtaupole, Mtpole,MZpole,MZ

      common/mssminputs/ m0, m12, m10, m20, sgnmu, tanbeta, a0
      common/sminputs/ mbpole, mtaupole, Mtpole,MZpole

!-----------------------------------------------------------------
!      include 'stdinputs.h'
c-----------------------------------------------------------------
C     MU and bmu Part - Tree Level
C ----------------------------------------------------------------

      MZ = MZpole
      beta  = datan(tanbetaQ)
      sbeta = dsin(beta)
      cbeta = dcos(beta)
      
      cos2beta = cbeta*cbeta - sbeta*sbeta
      sin2beta = 2.d0*cbeta*sbeta
      tan2beta = sin2beta/cos2beta
      
      mh1sq = mh1mz
      mh2sq = mh2mz

      flagdet = 0
      flagDflat = 0

C--------------------------------------------------------------------

      mursq = (((mh1sq*sbeta*sbeta - mh2sq*cbeta*cbeta)/cos2beta)
     $     - 0.5d0*MZ*MZ)

      bmur = (0.5d0*sin2beta*(mh2sq + mh1sq + 2.d0*mursq))
      
c$$$      print*,"mu in treelavelmu = ", dsqrt(mursq)
c$$$      print*,"bmu in treelavelmu = ", bmur, bmur/(sbeta*cbeta)
      
C--------------------------------------------------------------------
!     no D flat directions or detisneg (TL)?
!---------------------------------------------------------

      noDflat = 2.d0*(mursq - bmur) + mh2sq + mh1sq 

      detisneg = (mursq+mh1sq)*(mursq+mh2sq) - bmur*bmur

      sud1if: if(detisneg.lt.0.d0)then
         flagdet = 1            !<---------det is negative satisfied
!       write(*,*)'detisneg '
      endif sud1if

      sud2if: if(noDflat.gt.0.d0)then
         flagDflat = 1          !<-----------Dflat satisfied
!      write(*,*)'dflat satisfied'

      endif sud2if


C-------------------------------------------
C     choosing sign of \mu
C-------------------------------------------


      muge:  if(mursq.ge.0)then
         
         mur = dsqrt(mursq)

      else 

         mur = dsqrt(dabs(mursq))
         
      endif muge
      

      RETURN

      END SUBROUTINE mutreelevel
C==============================================================================
C==============================================================================

      
