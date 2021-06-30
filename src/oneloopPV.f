****f* SuSeFLAV/oneloopPV.f 
*  NAME
*    Subroutine a0,b0,b1.b22,b22t,f,funcg,h
*    function D0,D27,C0
*
*  SYNOPSIS
*    Solves one loop scalar integrals (Passarino-Veltman Integrals)  
*
*  FUNCTION 
*    Solves one loop scalar integrals (Passarino-Veltman Integrals)
*    using analytical expressions given the scale and external momenta.    
*
*  INPUTS
*
*     p         -  external momenta
*     q         -  enegry scale
*   m1,m2,m3,m4 -  input Masses 
*
*  RESULT
*     
*     ans       - solution of the intergral
*
*  EXAMPLE
*     SUBROUTINE b0(p,m1,m2,q,ans)
*     SUBROUTINE a0(m,q,ans)
*     SUBROUTINE b1(p,m1,m2,q,ans)
*     SUBROUTINE b22(p,m1,m2,q,ans)
*     SUBROUTINE f(p,m1,m2,q,ans)
*     SUBROUTINE funcg(p,m1,m2,q,ans)
*     SUBROUTINE b22t(p,m1,m2,q,ans)
*     double complex function xlogx(x)
*     DOUBLE PRECISION FUNCTION D0(m1,m2,m3,m4)            !<-------check function and impose constraints
*     DOUBLE PRECISION FUNCTION D27(m1,m2,m3,m4)            !<-------check function and impose constraints
*     double precision function C0(m1,m2,m3)
*      
*  NOTES
*     We use the convention of BPMZ. 
*  BUGS
*    ---
*  SEE ALSO
*
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!=======================================================================================
C-----------------------------------------------------------------------------------

       SUBROUTINE b0(p,m1,m2,q,ans)

       IMPLICIT NONE

       DOUBLE PRECISION s,p,m1,m2,q,ans
       DOUBLE PRECISION mMax, mMin
       COMPLEX*8 fbxplus,fbxminus                   ! complex declaration must be *8 for the correct functioning  
                                                    ! of the intrinsic function clog
       COMPLEX*8 xplus,xminus                 
       COMPLEX*8 iEpsilon                     

       ans = 0.d0

       if(p.lt.5.d-3)then


          if(dabs(1.d0 - min(m1**2.d0,m2**2.d0)/max(m1**2.d0,m2**2.d0))
     $         .lt.2.d-6)then   !<if m1 = m2 --------protection
             ans  = - dlog((m1 / q)**2.d0)
             return
          else
             if(m1**2.d0.lt.4.d-8) ans = 1.d0 + dlog((q/m2)**2.d0)
             
             if(m2**2.d0.lt.4.d-8) ans = 1.d0 + dlog((q/m1)**2.d0)

             
             if(m1**2.d0.ge.4.d-8.and.m2**2.d0.ge.4.d-8) ans = -1.d0 *  
     $            dlog(Max(m1**2.d0,m2**2.d0)/q**2.d0) + 1.d0 + 1.d0 * 
     $            min(m1**2.d0,m2**2.d0) * 
     $            dlog(max(m1**2.d0,m2**2.d0)/min(m1**2.d0,m2**2.d0))/
     $            (min(m1**2.d0,m2**2.d0) - max(m1**2.d0,m2**2.d0))
             
          endif
          RETURN

       else

          if(max(m2**2.d0, m1**2.d0).lt.4.d-8) then           

             ans =  2.d0 - dlog((p/q)**2.d0) !(- 2.d0 * dlog(p)  - 2.d0 * q**2.d0)/q**2.d0
             return

          else

             mMax = max(m1**2.d0, m2**2.d0)
             mMin = min(m1**2.d0,m2**2.d0)

             iEpsilon = ( 0.d0 , 0.00000000001d0)

             s = p**2.d0 - mMin + mMax

             xplus = (s + sqrt(s*s - (4.d0*p*p*((mMax) - iEpsilon))))/
     $            (2.d0*p*p)
             
             xminus = (s - sqrt(s*s - (4.d0*p*p*((mMax) - iEpsilon))))/
     $            (2.d0*p*p)
             
             
             fbxplus = zlog(1.d0 - xplus) - xplus * 
     $            zlog(1.d0 - (1.d0/xplus)) - cmplx(1.d0) 
             
             fbxminus = zlog(1.d0 - xminus) - xminus * 
     $            zlog(1.d0 - (1.d0/xminus)) - cmplx(1.d0) 


             ans = REAL(-2.d0 * dlog(p / q) - fbxplus - fbxminus)

          endif
!     endif
       endif

       RETURN

       END SUBROUTINE b0

C----------------------------------------------------------------------
      
      SUBROUTINE a0(m,q,ans)
      
      IMPLICIT NONE 
      
      DOUBLE PRECISION m,q,ans
      
      if(m**2.d0.eq.0.d0)then
         ans = 0.d0
         RETURN
      else
         continue
      endif
      
      ans = m*m*(1.d0 - dlog((m/q)**2.d0))
      
      RETURN
      
      END SUBROUTINE a0
      
      
!----------------------------------------------------------------------

      SUBROUTINE b1(p,m1,m2,q,ans)
      
      IMPLICIT NONE 
      
      DOUBLE PRECISION p,m1,m2,q,ans,x
      DOUBLE PRECISION s, a0m2,a0m1,b0pm1m2

      EXTERNAL a0,b0

      x = (m2/m1)**2.d0
      if(x.eq.0.d0) x = 1.d-5

      if(x.eq.(x+1)) x = 1.d-5

      if(p.lt.5.d-3.OR.p.lt.91.2d0)then

         if(x.ge.1.d0) then
            ans = .5d0*(.5d0+1.d0/(1.d0-x)+dlog(x)/(1.d0-x)**2.d0 - 
     $           dlog(Max(m1**2.d0,m2**2.d0)/q**2.d0))
            return
         else
            ans = .5d0*(.5d0+1.d0/(1.d0-x)+dlog(x)/(1.d0-x)**2.d0 
     $           -dlog(x) - dlog(Max(m1**2.d0,m2**2.d0)/q**2.d0))
            return
         endif
         
      else
!--------------------------------


         s = p*p - m2*m2 + m1*m1

         call a0(m2,q,a0m2)
         call a0(m1,q,a0m1)
         
         call b0(p,m1,m2,q,b0pm1m2)

         ans = (a0m2 - a0m1 + (s*b0pm1m2))/(2*p*p)
      endif

      RETURN

      END SUBROUTINE b1

!----------------------------------------------------------------------


      SUBROUTINE b22(p,m1,m2,q,ans)

      IMPLICIT NONE

      DOUBLE PRECISION p,m1,m2,q,ans
      DOUBLE PRECISION  a0m2,a0m1,b0pm1m2
             
      EXTERNAL a0,b0

       ans = 0.d0

      call a0(m2,q,a0m2)
      call a0(m1,q,a0m1)
      
      call b0(p,m1,m2,q,b0pm1m2)


	if (p.lt.5.d-3.and.m1**2.d0.le.4.d-8)then     ! m1, m2 = 0

	  ans = 0.375d0 * (m2)**2.d0 - 0.25d0 * (m2)**2.d0 *  
     $          dlog((m2 / q)**2.d0)

	return
        else

        if(p.lt.5.d-3.and.m2**2.d0.le.4.d-8)then
	  ans = 0.375d0 * (m1)**2.d0 - 0.25d0 * (m1)**2.d0 * 
     $          dlog((m1 / q)**2.d0)
                
          return
          else


      if(p.lt.5.d-3.and.dabs(1.d0 - min(m1**2.d0,m2**2.d0)/
     $            max(m1**2.d0,m2**2.d0)).lt.2.d-6)then                    ! checked <m1=m2 approx--------protection
 

      ans = -(m1)**2.d0 * dlog((m1/q)**2.d0) *0.5d0 + (m1)**2.d0 * 0.5d0


      return
      else

      if(p.lt.5.d-3.and.dabs(1.d0 - min(m1**2.d0,m2**2.d0)/
     $        max(m1**2.d0,m2**2.d0)).gt.2.d-6)then


	ans = 0.375d0 * ((m1)**2.d0 + (m2)**2.d0) - 0.25d0 * 
     $     (((m2)**2.d0)**2.d0 * log((m2 / q)**2.d0) - ((m1))**4.d0 * 
     $       dlog((m1 / q)**2.d0)) / ((m2)**2.d0 - (m1)**2.d0) 



          return
          else


         if(p.gt.5.d-3) then
      ans = (0.5d0*(a0m1 + a0m2) + (m1*m1 + m2*m2 - 0.5d0*p*p)*b0pm1m2
     $     + ((m2*m2 - m1*m1)/(2.d0*p*p)) * (a0m2 - a0m1 - 
     $     (m2*m2 - m1*m1)*b0pm1m2) + m1*m1 + m2*m2 -(p*p/3.d0))/6.d0
      
      

      RETURN
      endif
      endif
      endif
      endif
      endif

      return
      END SUBROUTINE b22

!--------------------------------------------------------------------------

      
      SUBROUTINE f(p,m1,m2,q,ans)
      
      IMPLICIT NONE 
      
      DOUBLE PRECISION p,m1,m2,q,ans
      DOUBLE PRECISION  a0m2,a0m1,b0pm1m2
      
      EXTERNAL a0,b0
      
      call a0(m2,q,a0m2)
      call a0(m1,q,a0m1)
      
      call b0(p,m1,m2,q,b0pm1m2)

      ans = a0m1 - 2.d0*a0m2 -((2.d0 * p*p + 2.d0 * m1*m1 - m2*m2) * 
     $     b0pm1m2)
      
      RETURN
      
      END SUBROUTINE f
      
      
!--------------------------------------------------------------------


      SUBROUTINE funcg(p,m1,m2,q,ans)
      
      IMPLICIT NONE 
      
      DOUBLE PRECISION p,m1,m2,q,ans
      DOUBLE PRECISION  a0m2,a0m1,b0pm1m2

      EXTERNAL a0,b0
    

      call a0(m2,q,a0m2)
      call a0(m1,q,a0m1)
      
      call b0(p,m1,m2,q,b0pm1m2)

      ans = - a0m1 - a0m2  + ((p*p - m1*m1 - m2*m2)*b0pm1m2)

      RETURN

      END SUBROUTINE funcg
!------------------------------------------------------------------------


      SUBROUTINE h(p,m1,m2,q,ans)
      
      IMPLICIT NONE 
      
      DOUBLE PRECISION p,m1,m2,q,ans
      DOUBLE PRECISION gpm1m2,b22pm1m2

      EXTERNAL funcg,b22
      
      call funcg(p,m1,m2,q,gpm1m2)
      call b22(p,m1,m2,q,b22pm1m2)

      ans = 4.d0*b22pm1m2 + gpm1m2

      RETURN

      END SUBROUTINE h
!-----------------------------------------------------------------------

      SUBROUTINE b22t(p,m1,m2,q,ans)
      
      IMPLICIT NONE 
      
      DOUBLE PRECISION p,m1,m2,q,ans
      DOUBLE PRECISION a0m1,a0m2,b22pm1m2

      EXTERNAL a0,b22
      
      call a0(m1,q,a0m1)
      call a0(m2,q,a0m2)
      call b22(p,m1,m2,q,b22pm1m2)

      ans = b22pm1m2 - (a0m1 + a0m2)/4.d0

      RETURN

      END SUBROUTINE b22t

!-----------------------------------------------------------------------
!====================================================================================
!-----------------------------------------------------------------------
      double complex function xlogx(x)
      implicit none
      double complex x

      if(abs(x) .eq. 0) then
         xlogx = 0
      else
         xlogx = x*log(x)
      endif
      end
!-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION D0(m1,m2,m3,m4)            !<-------check function and impose constraints

      IMPLICIT NONE 
      DOUBLE PRECISION m1, m2, m3, m4, ans
      DOUBLE PRECISION C0
      

      if (dabs(1.d0 - min(m1**2.d0, m2**2.d0)/max(m1**2.d0,m2**2.d0))
     $     .lt.2.d-2) then  

         if (dabs(1.d0 - min(m2, m3)/max(m2,m3)).lt.2.d-2.and.
     $        dabs(1.d0 - min(m1, m2)/max(m1,m2)).lt.2.d-2) then

            ans =  1.d0 / (6.d0 * (m2)**4.d0)

             D0 = ans
             RETURN
         
         else 
         if (dabs(1.d0 - min(m2**2.d0, m3**2.d0)/max(m2**2.d0,m3**2.d0))
     $           .lt.2.d-2)then
            ans =
     $           ((m2)**4.d0 - (m4)**4.d0 + 2.d0 * m4**2.d0 * m2**2.d0 * 
     $           dlog((m4 / m2)**2.d0)) / 
     $           (2.d0 * m2**2.d0 * (m2**2.d0 - m4**2.d0)**2.d0 * 
     $           (m2**2.d0 - m4**2.d0))

              D0 = ans
              RETURN
          else 
            if (dabs(1.d0 - min(m2**2.d0, m4**2.d0)/
     $            max(m2**2.d0,m4**2.d0)).lt.2.d-2) then 
               ans =
     $        ((m2)**4.d0 - (m3)**4.d0 + 2.d0 * m3**2.d0 * m2**2.d0 * 
     $              dlog((m3 / m2)**2.d0)) / 
     $              (2.0 * m2**2.d0 * (m2**2.d0 - m3**2.d0)**2.d0 * 
     $              (m2**2.d0 - m3**2.d0))

                D0 = ans
                RETURN
            else 
               if (dabs(1.d0 - min(m3**2.d0, m4**2.d0)/
     $              max(m3**2.d0,m4**2.d0)).lt.2.d-2) then 
                  ans = 
     $                 -1.d0 / (m2**2.d0 - m3**2.d0)**2.d0 * 
     $           ((m2**2.d0 + m3**2.d0) / (m2**2.d0 - m3**2.d0) *  
     $                 dlog((m3 / m2)**2.d0) + 2.d0)
                  D0 = ans
                  RETURN
               else

                  ans = (m4**2.d0 / (m2**2.d0 - m4**2.d0)**2.d0 *  
     $                 dlog((m4 / m2)**2.d0) + 
     $                 m4**2.d0 / (m2**2.d0 * (m2**2.d0 - m4**2.d0)) -
     $                 m3**2.d0 / (m2**2.d0 - m3**2.d0)**2.d0 *
     $                 dlog((m3 / m2)**2.d0) -
     $                 m3**2.d0 / (m2**2.d0 * (m2**2.d0 - m3**2.d0))) / 
     $                 (m3**2.d0 - m4**2.d0)
                  D0 = ans
                  RETURN
                  
               endif
             endif       
            endif
         endif 
      ans = 1.d0/(m1**2.d0 - m2**2.d0)*(C0(m1,m3,m4) - C0(m2,m3,m4) )
      D0 = ans
      RETURN
         endif
      ans = 1.d0/(m1**2.d0 - m2**2.d0)*(C0(m1,m3,m4) - C0(m2,m3,m4) )
      D0 = ans
     

      END
!--------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION D27(m1,m2,m3,m4)            !<-------check function and impose constraints

      IMPLICIT NONE 
      DOUBLE PRECISION m1, m2, m3, m4, ans
      DOUBLE PRECISION C0
      
      ans = 0.25d0/(m1**2.d0 - m2**2.d0)*(m1**2.d0 * C0(m1,m3,m4) -
     $     m2**2.d0 * C0(m2,m3,m4) )
      
      D27 = ans
      RETURN
      END

!===============================================================================

      double precision function C0(m1,m2,m3)
      
      implicit none
      double precision m1,m2,m3,result


      if (dabs(1.d0 - min(m2**2.d0,m3**2.d0)/max(m2**2.d0,m3**2.d0))
     $     .lt.2.d-2) then

         if (dabs(1.d0 - min(m1**2.d0, m2**2.d0)/max(m1**2.d0,m2**2.d0))
     $        .lt.2.d-2) then 

            result = ( - 0.5d0 /(m2**2.d0)) 
            
         else 
            result = ((m1)**2.d0 / ((m1)**2.d0 - (m2)**2.d0 )**2.d0 * 
     $         dlog((m2/m1)**2.d0)  + 1.d0 / ((m1)**2.d0 - (m2)**2.d0)) 

         endif
         
         
      else
         if (dabs(1.d0 - min(m1**2.d0, m2**2.d0)/
     $        max(m1**2.d0,m2**2.d0)).lt.2.d-2) then        

           result = ( - ( 1.d0 + (m3)**2.d0 / ((m2)**2.d0 -(m3)**2.d0)* 
     $           dlog((m3/m2)**2.d0))
     $           / ((m2)**2.d0 -(m3)**2.d0) ) 
                       
         else
            if (dabs(1.d0 - min(m1**2.d0, m3**2.d0)/
     $           max(m1**2.d0,m3**2.d0)).lt.2.d-2) then   

            result = ( - (1.d0 + (m2)**2.d0 / ((m3)**2.d0 - (m2)**2.d0)* 
     $           dlog((m2/m3)**2.d0)) 
     $           / ((m3)**2.d0-(m2)**2.d0) )
           
            
         else 
            result = ( ((m2*m2/(m1*m1 - m2*m2)) * dlog((m2/m1)**2.d0)) - 
     $          ((m3*m3/(m1*m1 - m3*m3)) * dlog((m3/m1)**2.d0)))/(m2*m2- 
     $           m3*m3)

            endif
            endif
         endif

      c0 = result

      return
      end 

!==============================================================================
