****f* susyflav/smrge.f 
*  NAME
*    Subroutine smrge, smrgemt
*  SYNOPSIS
*     In this subroutine we write all the standard model renormalization 
*     group. 
*  FUNCTION
*     Computes the numerical values of all the sm rge at a given energy scale  
*  INPUTS
*     yy(31)    - Initial values for all RGEs
*     t          - energy scale 
*  RESULT
*     yy(31)    - RGE output at a scale t
*  EXAMPLE
*     subroutine smrge(t,yy,dydx)
*     subroutine smrgemt(t,yy,dydx) - SM rge running without top
*
*  NOTES
*     The notation we use closely follows that
*     of Arason, castano et al, PRD46(1192)3945.
*     Note that alpha1, alpha2, alpha3 in the following are 
*     further normalised by a "4 pi " factor. 
*     Remember all the yukawa matrices are also normalised by this "4 pi x 2"
*     factor. thus : yu(i,j) = [ 1/(4 pi x 2) ] (yu(i,j) ); where yu = yukawa
*     in the lagrangian and yu == yukawa in the program.
*     Also remember that all the A-parameters and the soft breaking masses 
*     are also scaled/normalised by the "4 pi x 2" factor. 
*	
*	dydx(1)-dydx(3) : yu(1,1) - yu(1,3) 
*	dydx(4)-dydx(6) : yu(2,1) - yu(2,3) 
*	dydx(7)-dydx(9) : yu(3,1) - yu(3,3) 
*
*	dydx(10)-dydx(12) : yd(1,1) - yd(1,3) 
*	dydx(13)-dydx(15) : yd(2,1) - yd(2,3) 
*	dydx(16)-dydx(18) : yd(3,1) - yd(3,3) 
*
*	dydx(19)-dydx(21) : ye(1,1) - ye(1,3) 
*	dydx(22)-dydx(24) : ye(2,1) - ye(2,3) 
*	dydx(25)-dydx(27) : ye(3,1) - ye(3,3) 
*	dydx(28)          : alph3 
*	dydx(39)          : alph2 
*	dydx(30)          : alph1 
*	dydx(31)          : vev
*
*  BUGS
*    ---
*  SEE ALSO
*    smrge
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C
!===========================================================================================

      subroutine smrge(t,yy_sm,dydx) 
      
      implicit none

      integer i, j, k,c

      double precision t,e, dydx(31), yy_sm(31)

      double precision tz, alph1, alph2, alph3
      double precision b1_ms, b2_ms, b3_ms,a1,a2,a3,vev

      double precision t1,t2u,t2d,t2e
      
      double precision yubu(3,3), ydbd(3,3), yebe(3,3)
      double precision yu(3,3), yd(3,3), ye(3,3)
      double precision yudag(3,3), yddag(3,3), yedag(3,3)
      double precision bu(3,3), bd(3,3), be(3,3)
      double precision yedagye(3,3),yddagyd(3,3),yudagyu(3,3)

      double precision tryedagye,tryddagyd,tryudagyu


      double precision mx 

      DOUBLE PRECISION beta1yu(9),beta2yu(9),beta1yd(9),beta2yd(9)
      DOUBLE PRECISION beta1ye(9),beta2ye(9)
!----------------------------------------------------------------------------
      double precision yuyudagyuyudagyu(3,3),yuyudagyuyddagyd(3,3)
      double precision yuyddagydyudagyu(3,3)
      double precision yuyddagydyddagyd(3,3), yuyddagyd(3,3),
     $     yuyudagyu(3,3)
      double precision lambda,ki4s,Y4s
!----------------------------------------------------------------------------
      double precision ydyddagydyddagyd(3,3),ydyddagydyudagyu(3,3)
      double precision ydyudagyuyddagyd(3,3),ydyudagyuyudagyu(3,3)
      double precision ydyddagyd(3,3),ydyudagyu(3,3)
!-------------------------------------------------------------------------------

      double precision yeyedagyeyedagye(3,3),yeyedagye(3,3)

!-------------------------------------------------------------------------------

      double precision yudagyuyudagyu(3,3),tryudagyuyudagyu
      double precision yddagydyddagyd(3,3),tryddagydyddagyd
      double precision yedagyeyedagye(3,3),tryedagyeyedagye
      double precision yudagyuyddagyd(3,3),tryudagyuyddagyd
!-----------------------------------------------------------------------------
      
      double precision mbpole, mtaupole, Mtpole, MZpole,pi,MZ

      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole

      external dag,matmult,trace,mat5pr,mat4pr,mat3prod 
      
!----------------
      include 'stdinputs.h'
!----------------

      pi = 4.d0 * datan(1.d0)

      MZ = MZpole
      
      b1_ms =  (41.d0/10.d0) 
      b2_ms = -(19.d0/6.d0)
      b3_ms = - 7.d0


      MX = 5.d0*(10**19.d0)

      tz =  dlog(MX**2.d0/MZ**2.d0)
      

      yu(1,1) = yy_sm(1)
      yu(1,2) = yy_sm(2)
      yu(1,3) = yy_sm(3)
      yu(2,1) = yy_sm(4)
      yu(2,2) = yy_sm(5)
      yu(2,3) = yy_sm(6)
      yu(3,1) = yy_sm(7)
      yu(3,2) = yy_sm(8)
      yu(3,3) = yy_sm(9) 
      
      yd(1,1) = yy_sm(10)    
      yd(1,2) = yy_sm(11)    
      yd(1,3) = yy_sm(12)    
      yd(2,1) = yy_sm(13)   
      yd(2,2) = yy_sm(14)   
      yd(2,3) = yy_sm(15)
      yd(3,1) = yy_sm(16)
      yd(3,2) = yy_sm(17)
      yd(3,3) = yy_sm(18)

      ye(1,1) = yy_sm(19)
      ye(1,2) = yy_sm(20)
      ye(1,3) = yy_sm(21)
      ye(2,1) = yy_sm(22)
      ye(2,2) = yy_sm(23)
      ye(2,3) = yy_sm(24)
      ye(3,1) = yy_sm(25)
      ye(3,2) = yy_sm(26)
      ye(3,3) = yy_sm(27)
      
      alph3 = yy_sm(28) 
      alph2 = yy_sm(29)
      alph1 = yy_sm(30)

      vev = yy_sm(31)

      call dag(yd,yddag) 
      loop3: do i=1,3
      loop4: do j=1,3
      
      enddo loop4	
      enddo loop3
      
      call matmult(yddag,yd,yddagyd) 
      call trace(yddagyd,tryddagyd)
     
      
      call dag(yu,yudag)
      call matmult(yudag,yu,yudagyu)    
      call trace(yudagyu,tryudagyu)          
      
      call dag(ye,yedag)
      call matmult(yedag,ye,yedagye)    
      call trace(yedagye,tryedagye) 
      
      t1 = (3.d0*tryudagyu + 3.d0*tryddagyd + tryedagye)


      lambda = (alph1 + alph2) / 8.d0 
!---------------------------------------------------------------------------------
      call mat4pr(yudag,yu,yudag,yu,yudagyuyudagyu)
      call trace(yudagyuyudagyu,tryudagyuyudagyu)
      call mat4pr(yddag,yd,yddag,yd,yddagydyddagyd)
      call trace(yddagydyddagyd,tryddagydyddagyd)
      call mat4pr(yedag,ye,yedag,ye,yedagyeyedagye)
      call trace(yedagyeyedagye,tryedagyeyedagye)
      call mat4pr(yudag,yu,yddag,yd,yudagyuyddagyd)
      call trace(yudagyuyddagyd,tryudagyuyddagyd)

!---------------------------------------------------------------------------------
      call mat5pr(yu,yudag,yu,yudag,yu,yuyudagyuyudagyu)
      call mat5pr(yu,yudag,yu,yddag,yd,yuyudagyuyddagyd)
      call mat5pr(yu,yddag,yd,yudag,yu,yuyddagydyudagyu)
      call mat5pr(yu,yddag,yd,yddag,yd,yuyddagydyddagyd)
      call mat3prod(yu,yddag,yd,yuyddagyd)
      call mat3prod(yu,yudag,yu,yuyudagyu)

!----------------------------------------------------------------------------------

      call mat5pr(yd,yddag,yd,yddag,yd,ydyddagydyddagyd)
      call mat5pr(yd,yddag,yd,yudag,yu,ydyddagydyudagyu)
      call mat5pr(yd,yudag,yu,yddag,yd,ydyudagyuyddagyd)
      call mat5pr(yd,yudag,yu,yudag,yu,ydyudagyuyudagyu)
      call mat3prod(yd,yddag,yd,ydyddagyd)
      call mat3prod(yd,yudag,yu,ydyudagyu)
      
!--------------------------------------------------------------------------------
      
      call mat5pr(ye,yedag,ye,yedag,ye,yeyedagyeyedagye)
      call mat3prod(ye,yedag,ye,yeyedagye)
!---------------------------------------------------------------------------------


      
      Y4s = (((17.d0/20.d0)*alph1 + (9.d0/2.d0)*alph2 +
     $     (8.d0)*alph3)*tryudagyu) + 
     $     (((1.d0/4.d0)*alph1 + (9.d0/4.d0)*alph2 +
     $     (8.d0)*alph3)*tryddagyd) +
     $     ((3.d0/4.d0)*(alph1 + alph2))*tryedagye	
      
!----------------------------------------------------------------------------------

      ki4s = ((27.d0/4.d0)*tryudagyuyudagyu +
     $     (27.d0/4.d0)*tryddagydyddagyd +
     $     (9.d0/4.d0)*tryedagyeyedagye - (3.d0/2.d0)*tryudagyuyddagyd)


!---------------------------------------------------------------------------------



      b1do: do i=1,3
      b2do: do j=1,3

      bu(i,j) = 1.5d0*(yudagyu (i,j) - yddagyd(i,j))               	
      bd(i,j) = -bu(i,j)
      be(i,j) = 1.5d0*(yedagye(i,j))

      enddo b2do

      bu(i,i) = bu(i,i) + (t1
     $     - 17.d0/20.d0*alph1 - 9.d0/4.d0*alph2 - 8.d0*alph3)
      bd(i,i) = bd(i,i) + ( t1
     $     - 1.d0/4.d0*alph1 - 9.d0/4.d0*alph2 - 8.d0*alph3)
      be(i,i) = be(i,i) + (t1
     $     - 9.d0/4.d0*(alph1 + alph2))

      enddo b1do
      
      t2u = ((17.d0/10.d0)*tryudagyu + (1.d0/2.d0)*tryddagyd
     $     +(3.d0/2.d0)*tryedagye)

      t2d = ((3.d0/2.d0)*tryudagyu + (3.d0/2.d0)*tryddagyd
     $     +(1.d0/2.d0)*tryedagye)
      
      t2e = ((2.d0)*tryudagyu + (2.d0)*tryddagyd)
      

      call matmult(yu,bu,yubu)
      call matmult(yd,bd,ydbd)
      call matmult(ye,be,yebe)



!-------------------------------------------------------------------
      c = 1
      rgesmi: do i=1,3
      rgesmj: do j=1,3
      
      beta1yu(c) =  - 0.5d0*yubu(i,j)	
      


      beta2yu(c) = -(3.d0/4.d0)*yuyudagyuyudagyu(i,j) + 
     $     (0.5d0)*yuyudagyuyddagyd(i,j) + 
     $     (1.d0/8.d0)*yuyddagydyudagyu(i,j) - 
     $     (11.d0/8.d0)*yuyddagydyddagyd(i,j) - 
     $     (0.5d0)*t1*(((5.d0/4.d0)*yuyddagyd(i,j) - 
     $     (9.d0/4.d0)*yuyudagyu(i,j))) +
     $     (0.5d0)*ki4s*yu(i,j) - (3.d0/4.d0)*lambda**2.d0*yu(i,j)+
     $     lambda*(3.d0*yuyudagyu(i,j) + yuyddagyd(i,j)) -
     $     ((223.d0/160.d0)*alph1 + (135.d0/32.d0)*alph2 + 
     $     (8.d0)*alph3)*yuyudagyu(i,j) +
     $     (((43.d0/160.d0)*alph1 -(9.d0/32.d0)*alph2+(8.d0)*alph3)*
     $     yuyddagyd(i,j)) - (5.d0/4.d0)*Y4s*yu(i,j) -
     $     (1187.d0/600.d0)*alph1**2.d0*yu(i,j) +
     $     (9.d0/40.d0)*alph1*alph2*yu(i,j) -
     $     (19.d0/30.d0)*alph1*alph3*yu(i,j) +
     $     (23.d0/8.d0)*alph2**2.d0*yu(i,j) -
     $     (4.5d0)*alph2*alph3*yu(i,j) +
     $     (54.d0)*alph3*alph3*yu(i,j)

!---------------------------------------------------------------------
      beta1yd(c) =  - 0.5d0*ydbd(i,j)

      beta2yd(c) = -(3.d0/4.d0)*ydyddagydyddagyd(i,j) + 
     $     (0.5d0)*ydyddagydyudagyu(i,j) + 
     $     (1.d0/8.d0)*ydyudagyuyddagyd(i,j) - 
     $     (11.d0/8.d0)*ydyudagyuyudagyu(i,j) - 
     $     (0.5d0)*t1*(((5.d0/4.d0)*ydyudagyu(i,j) - 
     $     (9.d0/4.d0)*ydyddagyd(i,j))) +
     $     (0.5d0)*ki4s*yd(i,j) - (3.d0/4.d0)*(lambda**2.d0)*yd(i,j)+
     $     lambda*(3.d0*ydyddagyd(i,j) + ydyudagyu(i,j)) -
     $     ((187.d0/160.d0)*alph1 + (135.d0/32.d0)*alph2 + 
     $     (8.d0)*alph3)*ydyddagyd(i,j) +
     $     (((79.d0/160.d0)*alph1 -(9.d0/32.d0)*alph2+(8.d0)*alph3)*
     $     ydyudagyu(i,j)) - (5.d0/4.d0)*Y4s*yd(i,j) -
     $     (127.d0/600.d0)*(alph1**2.d0)*yd(i,j) +
     $     (27.d0/40.d0)*alph1*alph2*yd(i,j) -
     $     (31.d0/30.d0)*alph1*alph3*yd(i,j) +
     $     (23.d0/8.d0)*(alph2**2.d0)*yd(i,j) -
     $     (4.5d0)*alph2*alph3*yd(i,j) +
     $     (54.d0)*alph3*alph3*yd(i,j)
      
!-----------------------------------------------------------------------------
      beta1ye(c) =   - 0.5d0*yebe(i,j)

      beta2ye(c) = -(3.d0/4.d0)*yeyedagyeyedagye(i,j) + 
     $     (9.d0/8.d0)*t1*yeyedagye(i,j) +
     $     (0.5d0)*ki4s*ye(i,j) - (3.d0/4.d0)*(lambda**2.d0)*ye(i,j)+
     $     lambda*(3.d0*yeyedagye(i,j)) -
     $     ((387.d0/160.d0)*alph1 +
     $     (135.d0/30.d0)*alph2)*yeyedagye(i,j) -
     $     (5.d0/4.d0)*Y4s*ye(i,j) -
     $     (1371.d0/400.d0)*(alph1*alph1)*ye(i,j) +
     $     (27.d0/40.d0)*(alph1*alph2)*ye(i,j) -
     $     (11.d0/8.d0)*(alph2*alph2)*ye(i,j)

      c= c+1

      enddo rgesmj
      enddo rgesmi
      
!-----------------------------------------------------------------


!------------------------------------------------------------------

      looprge: do k=1,9
      
      dydx(k) = beta1yu(k) + beta2yu(k)
      

      dydx(k+9) = beta1yd(k) + beta2yd(k)
      
      
      dydx(k+18) = beta1ye(k) + beta2ye(k) 
      
      enddo looprge


      dydx(28) = -(alph3**(2.d0))*(b3_ms)
     $     -((44.d0*alph1/5.d0 +
     $     12.d0*alph2/1.d0 - 26.d0*alph3)*alph3**(2.d0))+
     $     ((alph3**(2.d0))*t2e)  
      
      dydx(29) = -(alph2**(2.d0))*(b2_ms)
     $     -((27.d0*alph1/10.d0 +
     $     35.d0*alph2/6.d0 + 9.d0*alph3/2.d0)*alph2**(2.d0))+
     $     ((alph2**(2.d0))*t2d)  
      
      dydx(30) = -(alph1**(2.d0))*(b1_ms)
     $     -((199.d0*alph1/50.d0 +
     $     9.d0*alph2/10.d0 +(11.d0/10.d0)*alph3)*alph1**(2.d0))+
     $     ((alph1**(2.d0))*t2u)
      
!-------------------------------------------------------------------
C     RGE for vev
!------------------------------------------------------------------
      
      dydx(31) = ((t1/2.d0) - (9.d0/8.d0)*((alph1/5.d0) + alph2))*vev +
     $     (((3.d0/ 4.d0) * (lambda**2.d0)) + ((5.d0/4.d0) * Y4s) - !<---- 2 loop starts here
     $     (ki4s/2.d0) + (27.d0/160.d0) * alph1 * alph2 + 
     $     ((1293.d0/800.d0) * alph1 * alph1) - 
     $     ((271.d0/32.d0) * alph2 * alph2)) * vev

      

!-------------------------------------------------------------------


      e = (MX/(dexp(t/2.d0)))

      
      a1 = 1.d0/(alph1*4.d0*pi)
      a2 = 1.d0/(alph2*4.d0*pi)
      a3 = 1.d0/(alph3*4.d0*pi)
      
      return
      
      end subroutine smrge


!------------------------------------------------------------------------
!     SMRGE ENDS
!------------------------------------------------------------------------
!========================================================================
      
      subroutine smrgemt(t,yy_sm,dydx) 
      
      implicit none

!      include 'declare_smrge.f'

      integer i, j, k,c

      double precision t,e, dydx(31), yy_sm(31)

      double precision tz, alph1, alph2, alph3
      double precision b1_ms, b2_ms, b3_ms,a1,a2,a3,vev

      double precision t1,t2u,t2d,t2e
      
      double precision yubu(3,3), ydbd(3,3), yebe(3,3)
      double precision yu(3,3), yd(3,3), ye(3,3)
      double precision yudag(3,3), yddag(3,3), yedag(3,3)
      double precision bu(3,3), bd(3,3), be(3,3)
      double precision yedagye(3,3),yddagyd(3,3),yudagyu(3,3)


      double precision tryedagye,tryddagyd,tryudagyu

      double precision mx 

      DOUBLE PRECISION beta1yu(9),beta2yu(9),beta1yd(9),beta2yd(9)
      DOUBLE PRECISION beta1ye(9),beta2ye(9)
!----------------------------------------------------------------------------
      double precision yuyudagyuyudagyu(3,3),yuyudagyuyddagyd(3,3)
      double precision yuyddagydyudagyu(3,3)
      double precision yuyddagydyddagyd(3,3), yuyddagyd(3,3),
     $     yuyudagyu(3,3)
      double precision lambda,ki4s,Y4s
!----------------------------------------------------------------------------
      double precision ydyddagydyddagyd(3,3),ydyddagydyudagyu(3,3)
      double precision ydyudagyuyddagyd(3,3),ydyudagyuyudagyu(3,3)
      double precision ydyddagyd(3,3),ydyudagyu(3,3)
!-------------------------------------------------------------------------------

      double precision yeyedagyeyedagye(3,3),yeyedagye(3,3)

!-------------------------------------------------------------------------------

      double precision yudagyuyudagyu(3,3),tryudagyuyudagyu
      double precision yddagydyddagyd(3,3),tryddagydyddagyd
      double precision yedagyeyedagye(3,3),tryedagyeyedagye
      double precision yudagyuyddagyd(3,3),tryudagyuyddagyd
!-----------------------------------------------------------------------------

      double precision mbpole, mtaupole, Mtpole, MZpole,pi,MZ

      common/sminputs/ mbpole, mtaupole, Mtpole, MZpole


      external dag,matmult,trace,mat5pr,mat4pr,mat3prod 
      

!----------------------
      include 'stdinputs.h'
!-----------------------

      pi = 4.d0 * datan(1.d0)
      
      b1_ms =  (41.d0/10.d0) 
      b2_ms = -(19.d0/6.d0)
      b3_ms = - 7.d0


      MX = 5.d0*(10**19.d0)
     
      tz =  dlog(MX**2.d0/MZ**2.d0)
      

      yu(1,1) = yy_sm(1)
      yu(1,2) = yy_sm(2)
      yu(1,3) = yy_sm(3)
      yu(2,1) = yy_sm(4)
      yu(2,2) = yy_sm(5)
      yu(2,3) = yy_sm(6)
      yu(3,1) = yy_sm(7)
      yu(3,2) = yy_sm(8)
      yu(3,3) = yy_sm(9) 
      
      yd(1,1) = yy_sm(10)    
      yd(1,2) = yy_sm(11)    
      yd(1,3) = yy_sm(12)    
      yd(2,1) = yy_sm(13)   
      yd(2,2) = yy_sm(14)   
      yd(2,3) = yy_sm(15)
      yd(3,1) = yy_sm(16)
      yd(3,2) = yy_sm(17)
      yd(3,3) = yy_sm(18)

      ye(1,1) = yy_sm(19)
      ye(1,2) = yy_sm(20)
      ye(1,3) = yy_sm(21)
      ye(2,1) = yy_sm(22)
      ye(2,2) = yy_sm(23)
      ye(2,3) = yy_sm(24)
      ye(3,1) = yy_sm(25)
      ye(3,2) = yy_sm(26)
      ye(3,3) = yy_sm(27)
      
      alph3 = yy_sm(28) 
      alph2 = yy_sm(29)
      alph1 = yy_sm(30)

      vev = yy_sm(31)

      call dag(yd,yddag) 
      loop3: do i=1,3
      loop4: do j=1,3
      
      enddo loop4	
      enddo loop3
      
      call matmult(yddag,yd,yddagyd) 
      call trace(yddagyd,tryddagyd)
      
      
      call dag(yu,yudag)
      call matmult(yudag,yu,yudagyu)    
      call trace(yudagyu,tryudagyu)          
      
      call dag(ye,yedag)
      call matmult(yedag,ye,yedagye)    
      call trace(yedagye,tryedagye) 
      
      t1 = (3.d0*tryudagyu + 3.d0*tryddagyd + tryedagye)

      lambda = (alph1 + alph2) / 8.d0 
!---------------------------------------------------------------------------------
      call mat4pr(yudag,yu,yudag,yu,yudagyuyudagyu)
      call trace(yudagyuyudagyu,tryudagyuyudagyu)
      call mat4pr(yddag,yd,yddag,yd,yddagydyddagyd)
      call trace(yddagydyddagyd,tryddagydyddagyd)
      call mat4pr(yedag,ye,yedag,ye,yedagyeyedagye)
      call trace(yedagyeyedagye,tryedagyeyedagye)
      call mat4pr(yudag,yu,yddag,yd,yudagyuyddagyd)
      call trace(yudagyuyddagyd,tryudagyuyddagyd)

!---------------------------------------------------------------------------------
      call mat5pr(yu,yudag,yu,yudag,yu,yuyudagyuyudagyu)
      call mat5pr(yu,yudag,yu,yddag,yd,yuyudagyuyddagyd)
      call mat5pr(yu,yddag,yd,yudag,yu,yuyddagydyudagyu)
      call mat5pr(yu,yddag,yd,yddag,yd,yuyddagydyddagyd)
      call mat3prod(yu,yddag,yd,yuyddagyd)
      call mat3prod(yu,yudag,yu,yuyudagyu)

!----------------------------------------------------------------------------------

      call mat5pr(yd,yddag,yd,yddag,yd,ydyddagydyddagyd)
      call mat5pr(yd,yddag,yd,yudag,yu,ydyddagydyudagyu)
      call mat5pr(yd,yudag,yu,yddag,yd,ydyudagyuyddagyd)
      call mat5pr(yd,yudag,yu,yudag,yu,ydyudagyuyudagyu)
      call mat3prod(yd,yddag,yd,ydyddagyd)
      call mat3prod(yd,yudag,yu,ydyudagyu)
      
!--------------------------------------------------------------------------------
      
      call mat5pr(ye,yedag,ye,yedag,ye,yeyedagyeyedagye)
      call mat3prod(ye,yedag,ye,yeyedagye)
!---------------------------------------------------------------------------------


      
      Y4s = (((17.d0/20.d0)*alph1 + (9.d0/2.d0)*alph2 +
     $     (8.d0)*alph3)*tryudagyu) + 
     $     (((1.d0/4.d0)*alph1 + (9.d0/4.d0)*alph2 +
     $     (8.d0)*alph3)*tryddagyd) +
     $     ((3.d0/4.d0)*(alph1 + alph2))*tryedagye	
      
!----------------------------------------------------------------------------------

      ki4s = ((27.d0/4.d0)*tryudagyuyudagyu +
     $     (27.d0/4.d0)*tryddagydyddagyd +
     $     (9.d0/4.d0)*tryedagyeyedagye - (3.d0/2.d0)*tryudagyuyddagyd)

!---------------------------------------------------------------------------------



      b1do: do i=1,3
      b2do: do j=1,3

      bu(i,j) = 1.5d0*(yudagyu (i,j) - yddagyd(i,j))               	
      bd(i,j) = -bu(i,j)
      be(i,j) = 1.5d0*(yedagye(i,j))

      enddo b2do

      bu(i,i) = bu(i,i) + (t1
     $     - 17.d0/20.d0*alph1 - 9.d0/4.d0*alph2 - 8.d0*alph3)
      bd(i,i) = bd(i,i) + ( t1
     $     - 1.d0/4.d0*alph1 - 9.d0/4.d0*alph2 - 8.d0*alph3)
      be(i,i) = be(i,i) + (t1
     $     - 9.d0/4.d0*(alph1 + alph2))

      enddo b1do
      
      t2u = ((17.d0/10.d0)*tryudagyu + (1.d0/2.d0)*tryddagyd
     $     +(3.d0/2.d0)*tryedagye)

      t2d = ((3.d0/2.d0)*tryudagyu + (3.d0/2.d0)*tryddagyd
     $     +(1.d0/2.d0)*tryedagye)
      
      t2e = ((2.d0)*tryudagyu + (2.d0)*tryddagyd)
      

      call matmult(yu,bu,yubu)
      call matmult(yd,bd,ydbd)
      call matmult(ye,be,yebe)



!-------------------------------------------------------------------
      c = 1
      rgesmi: do i=1,3
      rgesmj: do j=1,3
      
      beta1yu(c) =  - 0.5d0*yubu(i,j)	
      


      beta2yu(c) = -(3.d0/4.d0)*yuyudagyuyudagyu(i,j) + 
     $     (0.5d0)*yuyudagyuyddagyd(i,j) + 
     $     (1.d0/8.d0)*yuyddagydyudagyu(i,j) - 
     $     (11.d0/8.d0)*yuyddagydyddagyd(i,j) - 
     $     (0.5d0)*t1*(((5.d0/4.d0)*yuyddagyd(i,j) - 
     $     (9.d0/4.d0)*yuyudagyu(i,j))) +
     $     (0.5d0)*ki4s*yu(i,j) - (3.d0/4.d0)*lambda**2.d0*yu(i,j)+
     $     lambda*(3.d0*yuyudagyu(i,j) + yuyddagyd(i,j)) -
     $     ((223.d0/160.d0)*alph1 + (135.d0/32.d0)*alph2 + 
     $     (8.d0)*alph3)*yuyudagyu(i,j) +
     $     (((43.d0/160.d0)*alph1 -(9.d0/32.d0)*alph2+(8.d0)*alph3)*
     $     yuyddagyd(i,j)) - (5.d0/4.d0)*Y4s*yu(i,j) -
     $     (1187.d0/600.d0)*alph1**2.d0*yu(i,j) +
     $     (9.d0/40.d0)*alph1*alph2*yu(i,j) -
     $     (19.d0/30.d0)*alph1*alph3*yu(i,j) +
     $     (23.d0/8.d0)*alph2**2.d0*yu(i,j) -
     $     (4.5d0)*alph2*alph3*yu(i,j) +
     $     (54.d0)*alph3*alph3*yu(i,j)

!---------------------------------------------------------------------
      beta1yd(c) =  - 0.5d0*ydbd(i,j)

      beta2yd(c) = -(3.d0/4.d0)*ydyddagydyddagyd(i,j) + 
     $     (0.5d0)*ydyddagydyudagyu(i,j) + 
     $     (1.d0/8.d0)*ydyudagyuyddagyd(i,j) - 
     $     (11.d0/8.d0)*ydyudagyuyudagyu(i,j) - 
     $     (0.5d0)*t1*(((5.d0/4.d0)*ydyudagyu(i,j) - 
     $     (9.d0/4.d0)*ydyddagyd(i,j))) +
     $     (0.5d0)*ki4s*yd(i,j) - (3.d0/4.d0)*(lambda**2.d0)*yd(i,j)+
     $     lambda*(3.d0*ydyddagyd(i,j) + ydyudagyu(i,j)) -
     $     ((187.d0/160.d0)*alph1 + (135.d0/32.d0)*alph2 + 
     $     (8.d0)*alph3)*ydyddagyd(i,j) +
     $     (((79.d0/160.d0)*alph1 -(9.d0/32.d0)*alph2+(8.d0)*alph3)*
     $     ydyudagyu(i,j)) - (5.d0/4.d0)*Y4s*yd(i,j) -
     $     (127.d0/600.d0)*(alph1**2.d0)*yd(i,j) +
     $     (27.d0/40.d0)*alph1*alph2*yd(i,j) -
     $     (31.d0/30.d0)*alph1*alph3*yd(i,j) +
     $     (23.d0/8.d0)*(alph2**2.d0)*yd(i,j) -
     $     (4.5d0)*alph2*alph3*yd(i,j) +
     $     (54.d0)*alph3*alph3*yd(i,j)
      
!-----------------------------------------------------------------------------
      beta1ye(c) =   - 0.5d0*yebe(i,j)

      beta2ye(c) = -(3.d0/4.d0)*yeyedagyeyedagye(i,j) + 
     $     (9.d0/8.d0)*t1*yeyedagye(i,j) +
     $     (0.5d0)*ki4s*ye(i,j) - (3.d0/4.d0)*(lambda**2.d0)*ye(i,j)+
     $     lambda*(3.d0*yeyedagye(i,j)) -
     $     ((387.d0/160.d0)*alph1 +
     $     (135.d0/30.d0)*alph2)*yeyedagye(i,j) -
     $     (5.d0/4.d0)*Y4s*ye(i,j) -
     $     (1371.d0/400.d0)*(alph1*alph1)*ye(i,j) +
     $     (27.d0/40.d0)*(alph1*alph2)*ye(i,j) -
     $     (11.d0/8.d0)*(alph2*alph2)*ye(i,j)

      c= c+1

      enddo rgesmj
      enddo rgesmi
      
!-----------------------------------------------------------------


!------------------------------------------------------------------

      looprge: do k=1,9
      
      dydx(k) = beta1yu(k) + beta2yu(k)
      

      dydx(k+9) = beta1yd(k) + beta2yd(k)
      
      
      dydx(k+18) = beta1ye(k) + beta2ye(k) 
      
      enddo looprge

      dydx(7) = 0.d0
      dydx(8) = 0.d0
      dydx(9) = 0.d0

      dydx(28) = -(alph3**(2.d0))*(b3_ms)
     $     -((44.d0*alph1/5.d0 +
     $     12.d0*alph2/1.d0 - 26.d0*alph3)*alph3**(2.d0))+
     $     ((alph3**(2.d0))*t2e)  
      
      dydx(29) = -(alph2**(2.d0))*(b2_ms)
     $     -((27.d0*alph1/10.d0 +
     $     35.d0*alph2/6.d0 + 9.d0*alph3/2.d0)*alph2**(2.d0))+
     $     ((alph2**(2.d0))*t2d)  
      
      dydx(30) = -(alph1**(2.d0))*(b1_ms)
     $     -((199.d0*alph1/50.d0 +
     $     9.d0*alph2/10.d0 +(11.d0/10.d0)*alph3)*alph1**(2.d0))+
     $     ((alph1**(2.d0))*t2u)
      
!-------------------------------------------------------------------
C     RGE for vev
!------------------------------------------------------------------
      
      dydx(31) = ((t1/2.d0) - (9.d0/8.d0)*((alph1/5.d0) + alph2))*vev +
     $     (((3.d0/ 4.d0) * (lambda**2.d0)) + ((5.d0/4.d0) * Y4s) - !<---- 2 loop starts here
     $     (ki4s/2.d0) + (27.d0/160.d0) * alph1 * alph2 + 
     $     ((1293.d0/800.d0) * alph1 * alph1) - 
     $     ((271.d0/32.d0) * alph2 * alph2)) * vev

      

!-------------------------------------------------------------------


      e = (MX/(dexp(t/2.d0)))

      
      a1 = 1.d0/(alph1*4.d0*pi)
      a2 = 1.d0/(alph2*4.d0*pi)
      a3 = 1.d0/(alph3*4.d0*pi)
      
      return
      
      end subroutine smrgemt


!------------------------------------------------------------------------
!     SMRGEmt ENDS
!------------------------------------------------------------------------
!========================================================================
      
