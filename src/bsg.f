****f* SuSeFLAV/bsg.f 
*  NAME
*    subroutine bsg
*  SYNOPSIS
*    Main routine for running SuSeFLAV. 
*
*  FUNCTION
*     This subroutine calculates the b->s,gamma rate
*     The functions follow from that of Bertolini-Borzumati-Masiero
*     And A. Bartl et. al (PRD 64 076009 (2001) ).
*     The SM NLO Contibutions and MSSM LO Contributions are added.
*     The functions are combinations of various functions presented
*     in the papers above. 
*
*  INPUTS
*     tanbeta  -   ratio of vevs
*     mT       -   top quark mass 
*     mB       -   bottom quark mass
*     Ceg      -   (1x2) Chargino eigen values
*     SUeg     -   (1x6) up-type squark eigenvalue
*     USU      -   (6x6) squark mixing matrix
*     OCL,OCR  -   (2 X 2) orthogonal matrices such that 
*                 MChar = Transpose[OCR].Diag[MChar].OCL 
*     mHchar   -   mass of charged higgs
*
*  RESULT
*
*     Bbsg     - b->s gamma decay rate    
*
*  EXAMPLE
*
*  subroutine bsg(tanbeta,mT,mB,Ceg,SUeg,USU,OCL,OCR,mHchar,Bbsg)
*
*  NOTES
*	mHchar = Mass of the Charged Higgs Boson (NOT Squared !!)
*  BUGS
*    ---
*  SEE ALSO
*    -----
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!==================================================================================================     

      subroutine bsg(tanbeta,mT,mB,Ceg,SUeg,USU,OCL,OCR,mHchar,Bbsg)
      implicit none 

      double precision OCR(2,2),Ceg(2),OCL(2,2),mT,mB
      double precision SUeg(6),USU(6,6), Bbsg,z,mHchar,xt,y(2,6)
      double precision lambdat,lambdab,c7smmw,c8smmw
      double precision beta,sbeta,cbeta,tanbeta
      double precision gg1,gg2
      double precision f1,f2,f3,fF4,f2mod,gsm1
      double precision c7h,c8h,c7ch,c8ch,r7,r8
      integer a, x, alp, alp2 
      double precision VCKM(3,3)

      double precision MWpole,MW
      
      common/mwpole/ MWpole

      
      common/VCKMparam/ VCKM

      include 'stdinputs.h' 

!--------------------------------

      MW = MWpole

      beta = datan(tanbeta)
      sbeta = dsin(beta)
      cbeta = dcos(beta)

      lambdat = mT/(dsqrt(2.d0)*MW*sbeta)

C	check the lambdab definition. 

      lambdab = mB/(dSqrt(2.d0)*MW*cbeta)

C	Defining y 

	do a = 1, 2
	do x = 1, 6 
	
	y(a,x) = SUeg(x)/(Ceg(a)**2.d0)

	enddo
	enddo 

	z = (mT**2.d0/mHchar**2.d0)

	xt = (mT*mT)/(MW*MW)

       c7smmw = (xt/4.d0)*f1(xt)
       c8smmw = (xt/4.d0)*gsm1(xt)

       c7h = (1.d0/6.d0)*((1.d0/2.d0)*z*(1.d0/tanbeta**2.d0)*f1(z)
     .	      + f2(z))

       c8h = (1.d0/6.d0)*((1.d0/2.d0)*z*(1.d0/tanbeta**2.d0)*gg1(z) 
     .	+ gg2(z))

	c7ch = 0.d0


c$$$       print*,"VCKM"
c$$$       opu20a: do i = 1, 3
c$$$ 207   format(1x,3(2x,1pe25.15))
c$$$       write(*,207) (VCKM(i,j),j = 1, 3)
c$$$      enddo opu20a  
c$$$
c$$$       print*,"USU"
c$$$       opu20b: do i = 1, 6
c$$$ 208   format(1x,6(2x,1pe25.15))
c$$$       write(*,208) (USU(i,j),j = 1, 6)
c$$$      enddo opu20b  
c$$$
c$$$       print*,"OCL"
c$$$       opu20c: do i = 1, 2
c$$$ 209   format(1x,2(2x,1pe25.15))
c$$$       write(*,209) (OCL(i,j),j = 1, 2)
c$$$      enddo opu20c  

	
	do a = 1, 2
	do x = 1, 6
	do alp = 1, 3
	do alp2 = 1, 3 

	c7ch =   VCKM(alp,3)*VCKM(alp2,2)/(VCKM(3,3)*VCKM(3,2)) 
     .   *(MW*MW/(Ceg(a)**2.d0)*( -0.5d0*f1(y(a,x))*(lambdat*
     .	USU(x,alp2)*OCL(a,2) - USU(x,alp2+3)* OCL(a,1))*(lambdat*
     .	USU(x,alp)*OCL(a,2) - USU(x,alp+3)*OCL(a,1)) + (Ceg(a)/mB)*
     .	f2mod(y(a,x))*USU(x,alp)*lambdab* OCR(a,2)*(lambdat*
     .  USU(x,alp2)*OCL(a,2) - USU(x,alp2+3)*OCL(a,1) )))

	enddo
	enddo 
	enddo 
	enddo 


	c7ch = (1.d0/3.d0)*c7ch 


        c8ch = 0.d0 
	
	do a = 1, 2
	do x = 1, 6
	do alp = 1, 3
	do alp2 = 1, 3 
	
	c8ch =  (VCKM(alp,3)*VCKM(alp2,2))/(VCKM(3,3)*VCKM(3,2))
     .  *( (MW*MW)/(Ceg(a)**2.d0)*(-0.5d0*f3(y(a,x))*
     .  (lambdat*USU(x,alp2)*OCL(a,2) - USU(x,alp2+3)*OCL(a,1) )* 
     .  (lambdat*USU(x,alp)*OCL(a,2) - USU(x,alp+3)*OCL(a,1)) +
     .   (Ceg(a)/Mb)*fF4(y(a,x))*USU(x,alp)*lambdab*OCR(a,2)*
     .   (lambdat*USU(x,alp2)*OCL(a,2) - USU(x,alp2+3)*OCL(a,1) ) ) )

	enddo
	enddo
	enddo 
	enddo 

	c8ch = (1.d0/3.d0)*c8ch 

        r7 = (c7smmw + c7h + c7ch)/c7smmw

        r8 = (c8smmw + c8h + c8ch)/c8smmw


        Bbsg = (1.258d0 + 0.382d0*r7*r7 + 0.015d0*r8*r8 + 1.395d0*r7 
     . 	+ 0.161d0*r8 + 0.083d0*r7*r8 )*(10.d0**(-4.d0))

        
	return 
	end subroutine bsg

!=============================================================================
	double precision function f1(x)
	implicit none 
	double precision  x 

        f1 = (-7.d0 + 5.d0*x + 8.d0*x**2.d0)/(6.d0*(1.d0 - x)**3.d0)
     . 	- (2.d0*x - 3.d0*x**2.d0)*dLog(x)/((1.d0 - x)**4.d0)

	end function 


	double precision function f2(x)
	implicit none 
	double precision  x 

       f2 = (3.d0*x - 5.d0*x**2.d0)/(2.d0*(1.d0 - x)**2.d0) + 
     . 	(2.d0*x - 3.d0*x**2.d0)*dLog(x)/((1.d0 -x)**3.d0)
	end function 

C     f2 is modified by adding 5/2 that is the limit for x->Infinity 

	double precision function f2mod(x)
	implicit none 
	double precision  x 

       f2mod = (3.d0*x - 5.d0*x**2.d0)/(2.d0*(1.d0 - x)**2.d0) + 
     . (2.d0*x - 3.d0*x**2.d0)*dLog(x)/((1.d0 -x)**3.d0)+(5.d0/2.d0)

	end function 

C	F2(x) in the notes and paper. 

	double precision function fF2(x)
	implicit none 
	double precision  x 

	fF2 = (2.d0*x**3.d0 + 3.d0*x**2.d0 - 6.d0*x + 1.d0 - 
     .	6.d0*x**2.d0*dLog(x))/(12.d0*(x-1.d0)**4.d0)

	end function 

	double precision function F4(x) 
	implicit none
	double precision x 

        F4 = (x**2.d0 -1.d0-2.d0*x*dLog(x))/(2.d0*(x-1.d0)**3.d0)

	end function

	double precision function f3(x) 
	implicit none
	double precision x 

	f3 = -(2.d0 + x*(3.d0 + (-6.d0 + x)*x) - 6.d0*
     .	       x*dLog(1.d0/x))/(2.d0*(-1.d0 + x)**4.d0)

	end function

C	f4(x) of the math program.

	double precision function fF4(x) 
	implicit none
	double precision x 

         fF4 = (3.d0*(-1.d0 + x**2.d0 + 2.d0*x*dLog(1.d0/x)))/
     .	(2.d0*(-1.d0 + x)**3.d0)

	end function
	

	double precision function gsm1(x)
	implicit none
	double precision x

	gsm1 = (3.d0*(2.d0 + x*(3.d0 + (-6.d0 + x)*x) - 
     .	6.d0*x*dLog(1.d0/x)))/(2.d0*(-1.d0 + x)**4.d0)

	end function 

	double precision function gg1(x)
	implicit none
	double precision x

	gg1 = (3.d0*(2.d0 + x*(3.d0 + (-6.d0 + x)*x) - 
     .	6.d0*x*Log(1.d0/x)))/(2.d0*(-1.d0 + x)**4.d0)


	end function 

	double precision function gg2(x)
	implicit none
	double precision x

	gg2 =  (-3.d0*(2.d0 + x*(3.d0 + (-6.d0 + x)*x) - 
     . 6.d0*x*dLog(1.d0/x)))/(2.d0*(-1.d0 + x)**3.d0)


	end function 
