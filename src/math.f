
* Program : SuSeFLAV: Supersymmetry seesaw and flavor violation calculator
* Version : 1.0
* Website : http://cts.iisc.ernet.in/Suseflav/main.html
* Authors : Debtosh Choudhury  debtosh@cts.iisc.ernet.in 
*           Raghuveer Garani   rgarani@cts.iisc.ernet.in
*           Sudhir Vempati     vempati@cts.iisc.ernet.in

C     Consists of User Defined Subroutines
!================================================================================================
C	RK4ROUTINE is used to solve the differential equations using fourth order
C	Runge-Kutte Method. The various parameters are defined as follows: 
C	
C	dyinitial(neq) : Initial conditions for the differential equations.
C		  neq : Number of differential equations.
C	          dt1 : first step size 	 
C		dtmin : minimum step size
C	           t1 : start value of the integration variable 
C	           t2 : end value of the integration variable 
C	           t  : integration variable 
C	           dt : step size for integration variable 
C		 dydt : differential equation in the variable t 
C	      diffeqs : subroutine where differential equations are written 
C	       QMSRK4 : subroutine for quality managed RK 4 step size 
C	    CSHKRPRK4 : subroutine for corrected RK 4 step size 
C	        delta : overall tolerance level 
C	        nsok  : nsok are the number of good steps taken 
C	       nsbad  : nsbad are the number of bad steps taken 
C	       nsbad  : nsbad are the number of bad steps taken 
C	       check  : variable to check underflow error in numerical evaluation of RGE 
C	----------------------------------------------------------------------------------------------
C					REST OF THE PARAMETERS ARE PURELY INTERNAL PURPOSES
!================================================================================================



      SUBROUTINE RK4ROUTINE(dyinitial,neq,t1,t2,delta,dt1,dtmin,nsok,
     . nsbad,diffeqs,QMSRK4,check)
      
      INTEGER nsbad,nsok,neq,KSMXX,MXnstpsS,NSTMAX
      integer, intent(out) :: check
!      integer, save :: check = 0
      double precision  delta,dt1,dtmin,t1,t2,dyinitial(neq),MINI
      PARAMETER (MXnstpsS=50000,NSTMAX=500,KSMXX=500,MINI=2.d-31)
      INTEGER i,ksmx,qount,nstps
      double precision  dtsav,dt,dtdid,dtnext,t,tsav,dydt(NSTMAX),
     . tp(KSMXX), y(NSTMAX),yp(NSTMAX,KSMXX),yscle(NSTMAX)
      
      EXTERNAL diffeqs,QMSRK4

C     =================================================================================
C		INTEGRATION ROUTINE STARTS 
C     =================================================================================
      ksmx = ksmxx
      dtsav = 0.1d0
 
      check = 0

      t=t1
      dt=sign(dt1,t2-t1)
      nsok=0
      nsbad=0
      qount=0
      do 11 i=1,neq
        y(i)=dyinitial(i)
11    continue
      if (ksmx.gt.0) tsav=t-2.d0*dtsav
      do 16 nstps=1,MXnstpsS
        call diffeqs(t,y,dydt)
        do 12 i=1,neq
          yscle(i)=dabs(y(i))+dabs(dt*dydt(i))+MINI
12      continue
        if(ksmx.gt.0)then
          if(dabs(t-tsav).gt.dabs(dtsav)) then
            if(qount.lt.ksmx-1)then
              qount=qount+1
              tp(qount)=t
              do 13 i=1,neq
                yp(i,qount)=y(i)
13            continue
              tsav=t
            endif
          endif
        endif
        if((t+dt-t2)*(t+dt-t1).gt.0.d0) dt=t2-t
    	call QMSRK4(y,dydt,neq,t,dt,delta,yscle,dtdid,dtnext,diffeqs,
     .  check)
	if(check.eq.100)then
	goto 97
	endif
        if(dtdid.eq.dt)then
          nsok=nsok+1
        else
          nsbad=nsbad+1
        endif
        if((t-t2)*(t2-t1).ge.0.d0)then
          do 14 i=1,neq
            dyinitial(i)=y(i)
14        continue
          if(ksmx.ne.0)then
            qount=qount+1
            tp(qount)=t
            do 15 i=1,neq
              yp(i,qount)=y(i)
15          continue
          endif
          return
        endif
        if(dabs(dtnext).lt.dtmin)then!pause
!     *'stepsize smaller than minimum in RK4ROUTINE'
        check=100
        goto 97
        endif
        dt=dtnext
16      continue
!      pause 'too many steps in RK4ROUTINE'
	check = 100
	goto 97
97     return
       END SUBROUTINE RK4ROUTINE

!---------------------------------------------------------------------

        SUBROUTINE CSHKRPRK4(y,dydt,n,t,dt,dyout,dyerr,diffeqs)   

        INTEGER n,NSTMAX
        double precision dt,t,dydt(n),y(n),dyerr(n),dyout(n)
        EXTERNAL diffeqs
        PARAMETER (NSTMAX=500)

        
CU    USES diffeqs
        INTEGER i
 
        double precision zk2(NSTMAX),zk3(NSTMAX),zk4(NSTMAX),
     .   zk5(NSTMAX),WC5,WC6,T41,T42,WC1,WC3,WC4,
     $  zk6(NSTMAX),dytemp(NSTMAX),Z2,Z3,Z4,Z5,Z6,T21,T31,T32,
     .  T43,T51,T52,T53,T54,T61,T62,T63,T64,T65,X1,X3,X4,X6

        PARAMETER (Z2=.2d0,Z3=.3d0,Z4=.6d0,Z5=1.d0,Z6=.875d0,T21=.2d0,
     .  T31=3.d0/40.d0,T32=9.d0/40.d0,T41=.3d0,T42=-.9d0,T43=1.2d0,
     .  T51=-11.d0/54.d0,T52=2.5d0,T53=-70.d0/27.d0,T54=35.d0/27.d0,
     .  T61=1631.d0/55296.d0,T62=175.d0/512.d0,T63=575.d0/13824.d0,
     .  T64=44275.d0/110592.d0,T65=253.d0/4096.d0,X1=37.d0/378.d0,
     *  X3=250.d0/621.d0,X4=125.d0/594.d0,X6=512.d0/1771.d0,
     .  WC1=X1-2825.d0/27648.d0, WC3=X3-18575.d0/48384.d0,
     .  WC4=X4-13525.d0/55296.d0,WC5=-277.d0/14336.d0,WC6=X6-.25d0)

        do 11 i=1,n
           dytemp(i)=y(i)+T21*dt*dydt(i)
11    continue
      call diffeqs(t+Z2*dt,dytemp,zk2)
      do 12 i=1,n
         dytemp(i)=y(i)+dt*(T31*dydt(i)+T32*zk2(i))
12    continue
      call diffeqs(t+Z3*dt,dytemp,zk3)

      do 13 i=1,n
         dytemp(i)=y(i)+dt*(T41*dydt(i)+T42*zk2(i)+T43*zk3(i))
13    continue
      call diffeqs(t+Z4*dt,dytemp,zk4)


      do 14 i=1,n
        dytemp(i)=y(i)+dt*(T51*dydt(i)+T52*zk2(i)+T53*zk3(i)
     .   +T54*zk4(i))
14    continue
      call diffeqs(t+Z5*dt,dytemp,zk5)


      do 15 i=1,n
        dytemp(i)=y(i)+dt*(T61*dydt(i)+T62*zk2(i)+T63*zk3(i)+T64
     .   *zk4(i)+ T65*zk5(i))
15    continue
      call diffeqs(t+Z6*dt,dytemp,zk6)


      do 16 i=1,n
         dyout(i)=y(i)+dt*(X1*dydt(i)+X3*zk3(i)+X4*zk4(i)+X6*zk6(i))
16    continue
      do 17 i=1,n
         dyerr(i)=dt*(WC1*dydt(i)+WC3*zk3(i)+WC4*zk4(i)+
     .    WC5*zk5(i)+WC6*zk6(i))
17    continue
      return
      END SUBROUTINE CSHKRPRK4
      

!------------------------------------------------------------------------
      
       SUBROUTINE QMSRK4(y,dydt,n,t,dttry,delta,yscle,dtdid,dtnext,
     .	diffeqs,check)


       integer, intent(out) :: check
!       integer, save :: check = 0
       INTEGER n,NSTMAX
       DOUBLE PRECISION delta,dtdid,dtnext,dttry,t,dydt(n),y(n),
     .  yscle(n)
       EXTERNAL diffeqs
       PARAMETER (NSTMAX=500)
CU     USES diffeqs,CSHKRPRK4
       INTEGER i
       DOUBLE PRECISION maxeror,dt,xnew,dyerr(NSTMAX),dytemp(NSTMAX),
     .  SECURE,GROWY,SHRIKY, CONEROR
       PARAMETER (SECURE=0.9d0,GROWY=-.2d0,SHRIKY=-.25d0,
     . CONEROR=1.89d-4)


       check = 0

      dt=dttry
1     call CSHKRPRK4(y,dydt,n,t,dt,dytemp,dyerr,diffeqs)

      maxeror=0.d0
      do 11 i=1,n
        maxeror=max(maxeror,abs(dyerr(i)/yscle(i)))
11    continue
      maxeror=maxeror/delta
      if(maxeror.gt.1.d0)then
        dt=SECURE*dt*(maxeror**SHRIKY)
        if(dt.lt.0.1d0*dt)then
          dt=.1d0*dt
        endif
        xnew=t+dt
        if(xnew.eq.t)then       !pause 'stepsize underflow in QMSRK4'
	check = 100
	goto 33
	endif 
       goto 1
      else
        if(maxeror.gt.CONEROR)then
          dtnext=SECURE*dt*(maxeror**GROWY)
        else
          dtnext=5.d0*dt
        endif
        dtdid=dt
        t=t+dt
        do 12 i=1,n
          y(i)=dytemp(i)
12      continue
        return
      endif
33      END
      
      
c$$$C     ====================================================================
c$$$C     -------------------------------------------------------------------
c$$$C     Portable Random Number Generator. Recommended by D. Knuth. 
c$$$C     iseed is an integer seed which can be set to any negative value. 
c$$$C     Use different seeds if you calling the same random number generator
c$$$C     at different points. 
c$$$C     -------------------------------------------------------------------
c$$$      
c$$$      FUNCTION knran(iseed)
c$$$      INTEGER iseed
c$$$      INTEGER RNBIG,RNSEED,RMZ
c$$$C     REAL RNBIG,RNSEED,RMZ
c$$$      REAL knran,RNFAC 
c$$$C      REAL RNFAC ! for modification
c$$$      PARAMETER (RNBIG=1000000000,RNSEED=161803398,RMZ=0,RNFAC=1./RNBIG)
c$$$C     PARAMETER (RNBIG=4000000.,RNSEED=1618033.,RMZ=0.,RNFAC=1./RNBIG)
c$$$      INTEGER i,iff,rnii,rinxt,rinxtp,k
c$$$      INTEGER rmj,rmk,rma(55)
c$$$C     REAL rmj,rmk,rma(55)
c$$$      SAVE iff,rinxt,rinxtp,rma
c$$$      DATA iff /0/
c$$$      if(iseed.lt.0.or.iff.eq.0)then
c$$$        iff=1
c$$$        rmj=RNSEED-iabs(iseed)
c$$$        rmj=mod(rmj,RNBIG)
c$$$        rma(55)=rmj 
c$$$        rmk=1
c$$$        do 11 i=1,54
c$$$          rnii=mod(21*i,55)
c$$$          rma(rnii)=rmk
c$$$          rmk=rmj-rmk
c$$$          if(rmk.lt.RMZ)rmk=rmk+RNBIG
c$$$          rmj=rma(rnii)
c$$$11      continue
c$$$        do 13 k=1,4
c$$$          do 12 i=1,55
c$$$            rma(i)=rma(i)-rma(1+mod(i+30,55))
c$$$            if(rma(i).lt.RMZ)rma(i)=rma(i)+RNBIG
c$$$12        continue
c$$$13      continue
c$$$        rinxt=0
c$$$        rinxtp=31
c$$$        iseed=1
c$$$      endif
c$$$      rinxt=rinxt+1
c$$$      if(rinxt.eq.56)rinxt=1
c$$$      rinxtp=rinxtp+1
c$$$      if(rinxtp.eq.56)rinxtp=1
c$$$      rmj=rma(rinxt)-rma(rinxtp)
c$$$      if(rmj.lt.RMZ)rmj=rmj+RNBIG
c$$$      rma(rinxt)=rmj
c$$$      knran=rmj*RNFAC
c$$$      return
c$$$      END
c$$$
c$$$C     ========================================================================
!
!////////////////////////////////////////////////////////////////////////////////
!     USER DEFINED SUBROUTINES: MATRIX MANIPULATIONS
!-------------------------------------------------------
!     subroutine dagger
!1------------------------------------------------------
      SUBROUTINE dag(A, B)
      IMPLICIT NONE
      INTEGER  i,j              !Matrix Dimensions
      DOUBLE PRECISION A(3,3)             !Matrix A 
      DOUBLE PRECISION B(3,3)             !Matrix B
      loopdi: DO i = 1, 3
      loopdj: DO j = 1, 3
      B(i,j) = 0.d0
      B(i,j) = A(j,i)
      
      END DO loopdj
      ENDDO loopdi
      RETURN
      
      END SUBROUTINE dag

!---check orthogonality
      SUBROUTINE chkorth(A, check)
      IMPLICIT NONE
      character*1 check
      INTEGER  i,j,c              !Matrix Dimensions
      INTEGER iorth(3,3)
      DOUBLE PRECISION A(3,3),AT(3,3),orthT(3,3)             !Matrix A 
      DOUBLE PRECISION B(3,3),orth(3,3)             !Matrix B
      external dag, matmult

      c = 0
      call dag(A,AT)
      call matmult(AT,A,orth)
      call matmult(A,AT,orthT)
      
      loopdi: DO i = 1, 3
      loopdj: DO j = 1, 3
      iorth(i,j) = 0.d0
      
      iorth(i,j) = orth(i,j)-orthT(i,j)

      if(iorth(i,j).eq.0)then
         c = c+1
      endif

      END DO loopdj
      ENDDO loopdi


      if(c.eq.9)then
         check = 'T'
      else
         check = 'F'
      endif


      RETURN
      
      END SUBROUTINE chkorth

!-------------------------------------------------------
! subroutine dagger
!1------------------------------------------------------
      SUBROUTINE dag4d(A, B)
      IMPLICIT NONE
      INTEGER  i,j              !Matrix Dimensions
      DOUBLE PRECISION A(4,4)             !Matrix A 
      DOUBLE PRECISION B(4,4)             !Matrix B
      loopd4i: DO i = 1, 4
      loopd4j: DO j = 1, 4
      B(i,j) = 0.d0
      B(i,j) = A(j,i)
      
      END DO loopd4j
      ENDDO loopd4i
      RETURN
      
      END SUBROUTINE dag4d   
      
!-------------------------------------------------------
! subroutine dagger
!1--------------------------------------------------------
      SUBROUTINE dag2d(A, B)
      IMPLICIT NONE
      INTEGER  i,j              !Matrix Dimensions
      DOUBLE PRECISION A(2,2)             !Matrix A 
      DOUBLE PRECISION B(2,2)             !Matrix B
      loopd2i: DO i = 1, 2
      loopd2j: DO j = 1, 2
      B(i,j) = 0.d0
      B(i,j) = A(j,i)
      
      END DO loopd2j
      ENDDO loopd2i
      RETURN
      
      END SUBROUTINE dag2d


!2-------------------------------------------------------      
! subroutine to compute product of two matrices
!--------------------------------------------------------
      SUBROUTINE matmult(A,B,C)

      IMPLICIT NONE
      double precision A(3,3)   !Matrix A
      double precision B(3,3)   !Matrix B
      double precision C(3,3)   !Matrix C
      
      INTEGER i, j, k 
      
!     computing product 
      loopr: DO i = 1, 3
      loopc: DO j = 1, 3
      C(i,j)=0.d0
      loopk: DO k = 1, 3
      C(i,j) = C(i,j) + A(i,k)*B(k,j)
      ENDDO  loopk
      ENDDO  loopc
      ENDDO loopr
      RETURN
      END SUBROUTINE matmult
!3-----------------------------------------------------------
!PURPOSE: SUBROUTINE TO COMPUTE TRACE AF A GIVEN MATRIX
!------------------------------------------------------------
	subroutine trace(A,tr)

        IMPLICIT NONE
        DOUBLE PRECISION A(3,3)  !Matrix A
        double precision tr
        INTEGER i
        tr=0.d0
        loopt: DO i=1, 3
                   tr=tr+A(i,i)
               ENDDO loopt 
        RETURN  
	END SUBROUTINE trace


!4----------------------------------------------------------
!PURPOSE:  Subroutine for multiplying a scalar to a matrix
!-----------------------------------------------------------
      SUBROUTINE scmul(A,sc,B)
      IMPLICIT NONE
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION B(3,3)
      DOUBLE PRECISION sc
      INTEGER i,j
        loopsi: DO i=1,3
        loopsj:      DO j=1,3
                        B(i,j)=0.d0
                        B(i,j)= sc*A(i,j)
                    ENDDO loopsj
                ENDDO loopsi
                               
       RETURN
       END SUBROUTINE scmul
!-----------------------------------------------------------


      SUBROUTINE scmul4d(A,sc,B)
      IMPLICIT NONE
      DOUBLE PRECISION A(4,4)
      DOUBLE PRECISION B(4,4)
      DOUBLE PRECISION sc
      INTEGER i,j
        loopsi: DO i=1,4
        loopsj:      DO j=1,4
                        B(i,j)=0.d0
                        B(i,j)= sc*A(i,j)
                    ENDDO loopsj
                ENDDO loopsi
                               
       RETURN
       END SUBROUTINE scmul4d

!5---------------------------------------------------------
!PURPOSE: SUBROUTINE TO COMPUTE ONE LOOP RGES,YUKAWA TERMS
!----------------------------------------------------------
      subroutine rgeb1(P, Q, R, S, T, b1)
      
      IMPLICIT NONE 
      external matmult,scmul,add,addfr,trace 
      DOUBLE PRECISION P(3,3),Q(3,3),R(3,3),S(3,3),T(3,3)
      DOUBLE PRECISION m1(3,3),m2(3,3),m3(3,3),m4(3,3)
      DOUBLE PRECISION ids(3,3),Ss(3,3)   ! gt
      DOUBLE PRECISION Ps(3,3),Rst(3,3),Qs(3,3) !,gts
      DOUBLE PRECISION b1(3,3)
      DOUBLE PRECISION trm1
      integer ro,co
      double precision id(3,3),idst(3,3)
      trm1=0.d0
      loopro: do ro=1,3
      loopco: do co=1,3
      
      m1(ro,co) = 0.d0
      m2(ro,co) = 0.d0
      m3(ro,co) = 0.d0
      m4(ro,co) = 0.d0
      ids(ro,co) = 0.d0
      b1(ro,co) = 0.d0
      Ss(ro,co) = 0.d0
      Ps(ro,co) = 0.d0
      Rst(ro,co) = 0.d0
      Qs(ro,co) = 0.d0
      
      enddo  loopco
      enddo  loopro
      
      id(1,1)=1.d0
      id(1,2)=0.d0
      id(1,3)=0.d0
      id(2,1)=0.d0
      id(2,2)=1.d0
      id(2,3)=0.d0
      id(3,1)=0.d0
      id(3,2)=0.d0
      id(3,3)=1.d0
      
      call scmul(P,1.5d0,Ps)
      call scmul(R,1.5d0,Rst)
      call scmul(Q,0.5d0,Qs)
      call scmul(S,0.5d0,Ss) 
      call add(Ps,Qs,m1)
      call add(Rst,Ss,m2)
      call trace(m1,trm1)
      call scmul(id,trm1,idst)
      call add(idst,m2,m4)
      
      
      call matmult(T,m4,b1)
      RETURN
      

      end subroutine rgeb1       
      
!6-------------------------------------------------
!SUBROUTINE TO COMPUTE TWO LOOP RGES, YUKAWA TERMS
!--------------------------------------------------
      SUBROUTINE rgeb2(A, B, C, D, E, F, G, H, I, J, K, L, M, N, b2)
      
      IMPLICIT NONE
      EXTERNAL  matmult,add,addfr,trace,scmul
      DOUBLE PRECISION  A(3,3),B(3,3),C(3,3),D(3,3),E(3,3),F(3,3)
      DOUBLE PRECISION  G(3,3), H(3,3),I(3,3),J(3,3),K(3,3)
      DOUBLE PRECISION  As(3,3), Bs(3,3),Cs(3,3),Ds(3,3),Fs(3,3),Gs(3,3)
      DOUBLE PRECISION  Hs(3,3),Is(3,3),Js(3,3),Es(3,3),Ks(3,3),ids(3,3)
      DOUBLE PRECISION  L(3,3),M(3,3),N(3,3),a9(3,3),b2(3,3) 
      DOUBLE PRECISION  a1(3,3),a2(3,3),a3(3,3),a4(3,3)
      DOUBLE PRECISION  a5(3,3),a6(3,3),a7(3,3),a8(3,3),id(3,3)
      INTEGER  ro,co
      DOUBLE PRECISION tra1, tra2, tra3
      tra1=0.d0
      tra2=0.d0
      tra3=0.d0
      loopro2:   do ro=1,3
      loopco2:       do co=1,3
      As(ro,co)=0.d0
      Bs(ro,co)=0.d0
      Cs(ro,co)=0.d0
      Ds(ro,co)=0.d0
      Fs(ro,co)=0.d0
      Gs(ro,co)=0.d0
      Hs(ro,co)=0.d0
      Is(ro,co)=0.d0
      Js(ro,co)=0.d0
      Es(ro,co)=0.d0
      Ks(ro,co)=0.d0
      ids(ro,co)=0.d0
      a1(ro,co)=0.d0
      a2(ro,co)=0.d0
      a3(ro,co)=0.d0
      a4(ro,co)=0.d0
      a5(ro,co)=0.d0
      a6(ro,co)=0.d0
      a7(ro,co)=0.d0
      a8(ro,co)=0.d0
      b2(ro,co)=0.d0
      end do    loopco2
      end do       loopro2
      id(1,1)=1.d0
      id(1,2)=0.d0
      id(1,3)=0.d0
      id(2,1)=0.d0
      id(2,2)=1.d0
      id(2,3)=0.d0
      id(3,1)=0.d0
      id(3,2)=0.d0
      id(3,3)=1.d0
      
      call scmul(A,4.5d0,As)
      call scmul(B,1.5d0,Bs)
      call scmul(C,1.5d0,Cs)
      call scmul(D,0.5d0,Ds)
      call scmul(F,1.5d0,Fs)
      call scmul(G,0.5d0,Gs)
      call addfr(As,Bs,Cs,Ds,a1)
      call trace(a1,tra1)

      call add(Fs,Gs,a2)
      call trace(a2,tra2)
      call scmul(E,tra2,Es)
      
      call scmul(I,1.5d0,Is)
      call scmul(J,0.5d0,Js)
      call add(Is,Js,a3)
      call trace(a3,tra3)
      call scmul(H,3.d0,Hs)
      call scmul(Hs,tra3,a4)
      
      call scmul(K,2.d0,Ks)
      call add(Ks,L,a5)
      call add(a5,M,a6)

      call add(a6,a4,a7)
      call add(a7,Es,a8)
      call scmul(id,tra1,ids)
      call add(ids,a8,a9)
      call matmult(N,a9,b2)
      
      RETURN

      end subroutine rgeb2

!7----------------------------------------------------
!SUBROUTINE TO COMPUTE ONE LOOP RGES FOR A PARAMETERS
!-----------------------------------------------------
       
        SUBROUTINE rgeA1(A,B,C,D,E,F,G,H,I,J,b1)
       
        IMPLICIT NONE
        EXTERNAL  matmult,add,addfr,trace,scmul
 	DOUBLE PRECISION  A(3,3),B(3,3),C(3,3),D(3,3),E(3,3),F(3,3)
        DOUBLE PRECISION  G(3,3), H(3,3),I(3,3),J(3,3)
        DOUBLE PRECISION  As(3,3), Bs(3,3),Cs(3,3),Ds(3,3),Fs(3,3),
     $       Gs(3,3)
        DOUBLE PRECISION  Hs(3,3),Is(3,3),Js(3,3),Es(3,3),Ks(3,3),
     $       ids(3,3)
	DOUBLE PRECISION  idm1(3,3),idm2(3,3),b1(3,3) 
        DOUBLE PRECISION  a1(3,3),a2(3,3),a3(3,3),a4(3,3)
        DOUBLE PRECISION  a5(3,3),a6(3,3),m1(3,3),m2(3,3)
        INTEGER  ro,co
        DOUBLE PRECISION trm1, trm2 
        double precision id(3,3)
        loopro2:   do ro=1,3
        loopco2:       do co=1,3
                           As(ro,co)=0.d0
                           Bs(ro,co)=0.d0
                           Cs(ro,co)=0.d0
                           Ds(ro,co)=0.d0
                           Fs(ro,co)=0.d0
                           Gs(ro,co)=0.d0
                           Hs(ro,co)=0.d0
                           Is(ro,co)=0.d0
                           Js(ro,co)=0.d0
                           Es(ro,co)=0.d0
                           Ks(ro,co)=0.d0
                           ids(ro,co)=0.d0
                           idm1(ro,co)=0.d0
                           idm2(ro,co)=0.d0
                           b1(ro,co)=0.d0
                           a1(ro,co)=0.d0
                           a2(ro,co)=0.d0
                           a3(ro,co)=0.d0
                           a4(ro,co)=0.d0
                           a5(ro,co)=0.d0
                           a6(ro,co)=0.d0
                           m1(ro,co)=0.d0
                           m2(ro,co)=0.d0
                           
                       end do    loopco2
                    end do       loopro2
        
          id(1,1)=1.d0
          id(1,2)=0.d0
          id(1,3)=0.d0
          id(2,1)=0.d0
          id(2,2)=1.d0
          id(2,3)=0.d0
          id(3,1)=0.d0
          id(3,2)=0.d0
          id(3,3)=1.d0
        
        trm1=0.d0
        trm2=0.d0
        call scmul(A,1.5d0,As)
        call scmul(B,0.5d0,Bs)
        call scmul(C,2.5d0,Cs)
        call scmul(D,0.5d0,Ds)
        call scmul(E,3.d0,Es)
        call scmul(G,2.d0,Gs)
         
        call add(As,Bs,m1)
        call add(Es,F,m2)

        call trace(m1,trm1)
        call trace(m2,trm2)
        call scmul(id,trm1,idm1)
        call scmul(id,trm2,idm2)

        call add(idm1,Cs,a1)
        call add(a1,Ds,a2)
        call add(idm2,Gs,a3)
        call add(H,a3,a4)

        call matmult(I,a2,a5)
        call matmult(J,a4,a6)
    
        call add(a5,a6,b1)
      
        RETURN
      
        END SUBROUTINE rgeA1

!8--------------------------------------------------------------
!SUBROUTINE TO COMPUTE TWO LOOP RGE FOR A PARAMETER, SECTION1
!---------------------------------------------------------------
        
        SUBROUTINE rgeA21(A, B, C, D, E, F, G, H, I, J, K, L, M, N, b21)
        
        IMPLICIT NONE
        EXTERNAL  matmult,add,trace,scmul,addfr
 	DOUBLE PRECISION  A(3,3),B(3,3),C(3,3),D(3,3),E(3,3),F(3,3)
        DOUBLE PRECISION  G(3,3), H(3,3),I(3,3),J(3,3),K(3,3)
        DOUBLE PRECISION  As(3,3), Bs(3,3),Cs(3,3),Ds(3,3),Fs(3,3),
     $       Gs(3,3)
        DOUBLE PRECISION  Hs(3,3),Is(3,3),Js(3,3),Es(3,3),Ks(3,3),
     $       ids(3,3)
	DOUBLE PRECISION  L(3,3),M(3,3),N(3,3),a9(3,3),b21(3,3),Ms(3,3) 
        DOUBLE PRECISION  a1(3,3),a2(3,3),a3(3,3),a4(3,3)
        DOUBLE PRECISION  a5(3,3),a6(3,3),a7(3,3),a8(3,3)
        INTEGER  ro,co
        DOUBLE PRECISION tra1, tra2, tra3
        double precision id(3,3)
        tra1=0.d0
        tra2=0.d0
        tra3=0.d0
       
        loopro2:   do ro=1,3
        loopco2:       do co=1,3
                           As(ro,co)=0.d0
                           Bs(ro,co)=0.d0
                           Cs(ro,co)=0.d0
                           Ds(ro,co)=0.d0
                           Fs(ro,co)=0.d0
                           Gs(ro,co)=0.d0
                           Hs(ro,co)=0.d0
                           Is(ro,co)=0.d0
                           Js(ro,co)=0.d0
                           Es(ro,co)=0.d0
                           Ks(ro,co)=0.d0
                           ids(ro,co)=0.d0
                           b21(ro,co)=0.d0
                           a1(ro,co)=0.d0
                           a2(ro,co)=0.d0
                           a3(ro,co)=0.d0
                           a4(ro,co)=0.d0
                           a5(ro,co)=0.d0
                           a6(ro,co)=0.d0
                           a7(ro,co)=0.d0
                           a8(ro,co)=0.d0
                           a9(ro,co)=0.d0
                       end do    loopco2
                    end do       loopro2
                    
          id(1,1)=1.d0
          id(1,2)=0.d0
          id(1,3)=0.d0
          id(2,1)=0.d0
          id(2,2)=1.d0
          id(2,3)=0.d0
          id(3,1)=0.d0
          id(3,2)=0.d0
          id(3,3)=1.d0 

         call scmul(A,4.5d0,As)
         call scmul(B,1.5d0,Bs)
         call scmul(C,1.5d0,Cs)
         call scmul(D,0.5d0,Ds)
         call scmul(F,1.5d0,Fs)
         call scmul(G,0.5d0,Gs)
 	 call addfr(As,Bs,Cs,Ds,a1)
	 call trace(a1,tra1)

	 call add(Fs,Gs,a2)
	 call trace(a2,tra2)
         call scmul(E,tra2,Es)
         
         call scmul(I,1.5d0,Is)
         call scmul(J,0.5d0,Js)
	 call add(Is,Js,a3)
	 call trace(a3,tra3)
         call scmul(H,5.d0,a4)
         call scmul(a4,tra3,Hs)
    
         call scmul(K,3.d0,Ks)
         call scmul(M,2.d0,Ms)
	 call add(Ks,L,a5)
	 call add(a5,Ms,a6)

	 call add(a6,Hs,a7)
         call add(a7,Es,a8)
         call scmul(id,tra1,ids)
	 call add(ids,a8,a9)
	 call matmult(N,a9,b21)
      
         RETURN

         end subroutine rgeA21

!9--------------------------------------------------------------
!SUBROUTINE TO COMPUTE TWO LOOP RGES FOR A PARAMETERS PART 2
!---------------------------------------------------------------

        SUBROUTINE rgeA22(A, B, C, D, E, F, G, H, I, J, K, L, M, N,O,
     $  P,Q,R,S,T,U,V,W,X,Y,b22)
        
        IMPLICIT NONE
        EXTERNAL  matmult,add,addfr,trace,scmul

 	DOUBLE PRECISION  A(3,3),B(3,3),C(3,3),D(3,3),E(3,3),F(3,3)
        DOUBLE PRECISION  G(3,3), H(3,3),I(3,3),J(3,3),K(3,3),L(3,3),
     $       M(3,3)
        DOUBLE PRECISION  N(3,3),O(3,3),P(3,3),Q(3,3),R(3,3),S(3,3),
     $       T(3,3)
        DOUBLE PRECISION  U(3,3),V(3,3),W(3,3),X(3,3),Y(3,3)


        DOUBLE PRECISION  As(3,3), Bs(3,3),Cs(3,3),Ds(3,3),Fs(3,3),
     $       Gs(3,3)
        DOUBLE PRECISION  Hs(3,3),Is(3,3),Js(3,3),Es(3,3),Ks(3,3),
     $       ids(3,3)
	DOUBLE PRECISION  Ls(3,3),Ms(3,3),Ns(3,3),Os(3,3),Ps(3,3),Qs(3,3)
        DOUBLE PRECISION  Rs(3,3),Ss(3,3),Ts(3,3),Us(3,3),Vs(3,3),
     $       Ws(3,3),Xs(3,3)
         

        DOUBLE PRECISION  a1(3,3),a2(3,3),a3(3,3),a4(3,3),a9(3,3),
     $       a10(3,3)
        DOUBLE PRECISION  a5(3,3),a6(3,3),a7(3,3),a8(3,3),m1(3,3),
     $       m2(3,3),m3(3,3)
        DOUBLE PRECISION  m4(3,3),m5(3,3),m6(3,3),m7(3,3),b22(3,3),
     $       null(3,3)
        
        INTEGER  ro,co
        DOUBLE PRECISION tra3,tra4,tra5,tra6,tra7
        double precision id(3,3)
        tra3=0.d0
        tra4=0.d0
        tra5=0.d0
        tra6=0.d0
        tra7=0.d0
        loopr02:   do ro=1,3
        loopc02:       do co=1,3
                           As(ro,co)=0.d0
                           Bs(ro,co)=0.d0
                           Cs(ro,co)=0.d0
                           Ds(ro,co)=0.d0
                           Fs(ro,co)=0.d0
                           Gs(ro,co)=0.d0
                           Hs(ro,co)=0.d0
                           Is(ro,co)=0.d0
                           Js(ro,co)=0.d0
                           Es(ro,co)=0.d0
                           Ks(ro,co)=0.d0
                           Ls(ro,co)=0.d0
                           Ms(ro,co)=0.d0
                           Ns(ro,co)=0.d0
                           Os(ro,co)=0.d0
                           Ps(ro,co)=0.d0
                           Qs(ro,co)=0.d0
                           Rs(ro,co)=0.d0
                           Ss(ro,co)=0.d0
                           Ts(ro,co)=0.d0
                           Us(ro,co)=0.d0
                           Vs(ro,co)=0.d0
                           Ws(ro,co)=0.d0
                           Xs(ro,co)=0.d0
                           ids(ro,co)=0.d0
                           b22(ro,co)=0.d0
                           a1(ro,co)=0.d0
                           a2(ro,co)=0.d0
                           a3(ro,co)=0.d0
                           a4(ro,co)=0.d0
                           a5(ro,co)=0.d0
                           a6(ro,co)=0.d0
                           a7(ro,co)=0.d0
                           a8(ro,co)=0.d0
                           a9(ro,co)=0.d0 
                           a10(ro,co)=0.d0
                           m1(ro,co)=0.d0
                           m2(ro,co)=0.d0
                           m3(ro,co)=0.d0
                           m4(ro,co)=0.d0
                           m5(ro,co)=0.d0
                           m6(ro,co)=0.d0
                           m7(ro,co)=0.d0
                           null(ro,co)=0.d0
    
                       end do    loopc02
                    end do       loopr02
        
          id(1,1)=1.d0
          id(1,2)=0.d0
          id(1,3)=0.d0
          id(2,1)=0.d0
          id(2,2)=1.d0
          id(2,3)=0.d0
          id(3,1)=0.d0
          id(3,2)=0.d0
          id(3,3)=1.d0

         call scmul(A,18.d0,As)
         call scmul(B,3.d0,Bs)
         call scmul(C,3.d0,Cs)
         call scmul(D,6.d0,Ds)
         call scmul(H,3.d0,Hs)
         call scmul(G,3.d0,Gs)
         call scmul(K,3.d0,Ks)
         call scmul(M,2.d0,Ms)
         call scmul(N,3.d0,Ns)
         call scmul(Q,3.d0,Qs)
         call scmul(S,3.d0,Ss)
         call scmul(T,4.d0,Ts)
         call scmul(U,2.d0,Us)
         call scmul(V,2.d0,Vs)
         call scmul(X,2.d0,Xs)

 	 call addfr(As,Bs,Cs,Ds,a1)
         call add(E,F,a2)
         call add(a1,a2,a3)
	 call trace(a3,tra3)

	 call add(Hs,I,a4)
	 call trace(a4,tra4)
         call scmul(Gs,tra4,m1)

         call add(Ks,L,a5)
         call trace(a5,tra5)
         call scmul(J,tra5,m2)

         call add(Ns,O,a6)
         call trace(a6,tra6)
         call scmul(Ms,tra6,m3)

         call add(Qs,R,a7)
         call trace(a7,tra7)
         call scmul(P,tra7,m4)

         call addfr(Ss,Ts,Us,Vs,a8)
         call addfr(a8,W,Xs,null,a9)
         
         call addfr(m1,m2,m3,m4,a10)
         call add(a10,a9,m5)
         call scmul(id,tra3,m6)
         call add(m5,m6,m7)

         call matmult(Y,m7,b22)
        
              
         RETURN

         END SUBROUTINE rgeA22

!10---------------------------------------------
!PURPOSE:  SUBROUTINE TO ADD GIVEN TWO MATRICES
!-----------------------------------------------
      SUBROUTINE add(A,B,C)

      IMPLICIT NONE
      DOUBLE PRECISION A(3,3),B(3,3)
      DOUBLE PRECISION C(3,3)

      INTEGER i,j
      
      loopai: DO  i = 1, 3
      loopaj:     DO   j = 1, 3
      C(i,j)=0.d0
      C(i,j) = A(i,j) + B(i,j)
      
      ENDDO   loopaj
      ENDDO       loopai
      
      RETURN
      END SUBROUTINE add
!-------------------------------------------------------------------
	SUBROUTINE add4d(A,B,C)

        IMPLICIT NONE
	DOUBLE PRECISION A(4,4),B(4,4)
	DOUBLE PRECISION C(4,4)

	INTEGER i,j
        
	loopai: DO  i=1,4
	loopaj:     DO   j=1,4
	                 C(i,j)=0.d0
                         C(i,j) = A(i,j) + B(i,j)
                   
                   ENDDO   loopaj
	       ENDDO       loopai
  
        RETURN
	END SUBROUTINE add4d

!11-----------------------------------------------------------
!PURPOSE: subroutine to compute the product of three matrices
!-------------------------------------------------------------

	SUBROUTINE  mat3prod(A, B, C, P)

 	IMPLICIT  NONE
        EXTERNAL  matmult 
                                                                           
        DOUBLE PRECISION m1(3,3)          !TEMPORARY, INTERMEDIATE MULTIPLIER M1
	DOUBLE PRECISION A(3,3)           !Matrix A
	DOUBLE PRECISION B(3,3)           !Matrix B
	DOUBLE PRECISION C(3,3)           !Matrix C
        DOUBLE PRECISION P(3,3)           !Matrix PRODUCT
        integer i,j

        loop03pr1: do i=1,3
        loop03pr2: do j=1,3
                   m1(i,j)=0.d0
        enddo loop03pr2
        enddo loop03pr1 
	
        call matmult(A, B, m1)
        call matmult(m1, C, P)
        RETURN
        END SUBROUTINE mat3prod



!12----------------------------------------------------------
!PURPOSE: SUBROUTINE TO COMPUTE THE PRODUCT OF FOUR MATRICES        
!------------------------------------------------------------
       SUBROUTINE mat4pr(A,B,C,D,P)
        
        IMPLICIT  NONE
        EXTERNAL  matmult 
        
        DOUBLE PRECISION m1(3,3),m2(3,3)  !TEMPORARY, INTERMEDIATE MULTIPLIER M1
	DOUBLE PRECISION A(3,3),C(3,3)    !Matrix A
	DOUBLE PRECISION B(3,3),D(3,3)    !Matrix B
        DOUBLE PRECISION P(3,3)           !Matrix PRODUCT
        integer i,j
        loop04pr1: do i=1,3
        loop04pr2: do j=1,3
                   m1(i,j) = 0.d0
                   m2(i,j) = 0.d0
                enddo loop04pr2
             enddo loop04pr1                                    
      
        call matmult(A, B, m1)
        call matmult(m1, C, m2)
        call matmult(m2, D, P)      
        
        RETURN
        END SUBROUTINE mat4pr

!13----------------------------------------------------------
!PURPOSE: SUBROUTINE TO COMPUTE THE PRODUCT OF Five MATRICES        
!------------------------------------------------------------
       SUBROUTINE mat5pr(A,B,C,D,E,P)
        
        IMPLICIT  NONE
        EXTERNAL  matmult 
                                                                           
        DOUBLE PRECISION m1(3,3),m2(3,3),m3(3,3) 
	DOUBLE PRECISION A(3,3),C(3,3)        
	DOUBLE PRECISION B(3,3),D(3,3),E(3,3) 
        DOUBLE PRECISION P(3,3)  
        integer i,j    
                                                               
        loop05pr1: do i=1,3
        loop05pr2: do j=1,3
                   m1(i,j)=0.d0
                   m2(i,j)=0.d0
                   m3(i,j)=0.d0
        enddo loop05pr2
        enddo loop05pr1                            
	
        call matmult(A, B, m1)
        call matmult(m1, C, m2)
        call matmult(m2,D,m3)
        call matmult(m3,E,P)      
        
        RETURN
        END SUBROUTINE mat5pr


!14---------------------------------------------------
!PURPOSE: SUBROUTINE TO COMPUTE THE SUM FOUR MATRICES
!-----------------------------------------------------
	SUBROUTINE addfr(A, B, C, D, a3)
          
        IMPLICIT NONE 
        INTEGER i,j
        DOUBLE PRECISION A(3,3),C(3,3)
        DOUBLE PRECISION B(3,3),D(3,3) 
        DOUBLE PRECISION a3(3,3)
          
        loopafi: DO i=1,3
        loopafj:     DO j=1,3
                       a3(i,j) = 0.d0 
                       a3(i,j)= A(i,j)+B(i,j)+C(i,j)+D(i,j)
                   ENDDO loopafj
               ENDDO loopafi
                      
        RETURN

	END SUBROUTINE addfr


!15-------------------------------------------------
! subroutine to compute product of two 4x4 matrices
!---------------------------------------------------
      SUBROUTINE matmult4d(A,B,C)
      
      IMPLICIT NONE
      double precision A(4,4)   !Matrix A
      double precision B(4,4)   !Matrix B
      double precision C(4,4)   !Matrix C
      
      INTEGER i, j, k 
      
!     computing product 
      loopr4: DO i = 1, 4
      loopc4: DO j = 1, 4
      C(i,j)=0.d0
      loopk4: DO k = 1, 4
      C(i,j) = C(i,j) + A(i,k)*B(k,j)
      ENDDO  loopk4
      ENDDO  loopc4
      ENDDO loopr4
      RETURN
      END SUBROUTINE matmult4d
!16-------------------------------------------------
! subroutine to compute product of two 2x2 matrices
!---------------------------------------------------
      SUBROUTINE matmult2d(A,B,C)

      IMPLICIT NONE
      double precision A(2,2)   !Matrix A
      double precision B(2,2)   !Matrix B
      double precision C(2,2)   !Matrix C
      
      INTEGER i, j, k 
      
!     computing product 
      loopr2: DO i = 1, 2
      loopc2: DO j = 1, 2
      C(i,j)=0.d0
      loopk2: DO k = 1, 2
      C(i,j) = C(i,j) + A(i,k)*B(k,j)
      ENDDO  loopk2
      ENDDO  loopc2
      ENDDO loopr2
      RETURN
      END SUBROUTINE matmult2d


!17---------------------------------------------------------------
!PURPOSE: subroutine to compute the product of three 4x4 matrices
!-----------------------------------------------------------------

	SUBROUTINE  mat3prod4d(A, B, C, P)

 	IMPLICIT  NONE
        EXTERNAL  matmult4d 
                                                                           
        DOUBLE PRECISION m1(4,4)          !TEMPORARY, INTERMEDIATE MULTIPLIER M1
	DOUBLE PRECISION A(4,4)           !Matrix A
	DOUBLE PRECISION B(4,4)           !Matrix B
	DOUBLE PRECISION C(4,4)           !Matrix C
        DOUBLE PRECISION P(4,4)           !Matrix PRODUCT
        integer i,j

        loop04pr1: do i=1,4
        loop04pr2: do j=1,4
                   m1(i,j)=0.d0
        enddo loop04pr2
        enddo loop04pr1 
	
        call matmult4d(A, B, m1)
        call matmult4d(m1, C, P)
        RETURN
        END SUBROUTINE mat3prod4d

!18---------------------------------------------------------------
!PURPOSE: subroutine to compute the product of three 2x2 matrices
!-----------------------------------------------------------------

      SUBROUTINE  mat3prod2d(A, B, C, P)

      IMPLICIT  NONE
      EXTERNAL  matmult2d 
      
      DOUBLE PRECISION m1(2,2)            !TEMPORARY, INTERMEDIATE MULTIPLIER M1
      DOUBLE PRECISION A(2,2)             !Matrix A
      DOUBLE PRECISION B(2,2)             !Matrix B
      DOUBLE PRECISION C(2,2)             !Matrix C
      DOUBLE PRECISION P(2,2)             !Matrix PRODUCT
      integer i,j

      loop02pr1: do i=1,2
      loop02pr2: do j=1,2
      m1(i,j)=0.d0
      enddo loop02pr2
      enddo loop02pr1 
      
      call matmult2d(A, B, m1)
      call matmult2d(m1, C, P)
      RETURN
      END SUBROUTINE mat3prod2d


!19---------------------------------------------------------------
!PURPOSE: subroutine to compute the product of three 2x2 matrices
!-----------------------------------------------------------------

	SUBROUTINE  rmat2d(theta,rmtheta)

 	IMPLICIT  NONE
                                                                           
        DOUBLE PRECISION rmtheta(2,2),theta        

        rmtheta(1,1) = dcos(theta)
        rmtheta(1,2) = dsin(theta)
        rmtheta(2,1) = -1.d0 * dsin(theta)
        rmtheta(2,2) = dcos(theta)

        RETURN
        END SUBROUTINE rmat2d
!20----------------------------------------------..-------------
!PURPOSE: SUBROUTINE TO COMPUTE THE SUM OF THREE (4X4) MATRICES
!--------------------------------------------------.------------
      SUBROUTINE add4dmat(A, B, C, a4)
      
      IMPLICIT NONE 
      INTEGER i,j
      DOUBLE PRECISION A(4, 4), C(4, 4)
      DOUBLE PRECISION B(4, 4)
      DOUBLE PRECISION a4(4, 4)
      
      loopafi: DO i = 1, 4
      loopafj: DO j = 1, 4

      a4(i,j) = 0.d0   

      a4(i,j)= A(i,j) + B(i,j) + C(i,j)

      ENDDO loopafj
      ENDDO loopafi
      
      RETURN

      END SUBROUTINE add4dmat


!21-----------------------------------------------------------
!PURPOSE: SUBROUTINE TO COMPUTE THE SUM OF TWO (4X4) MATRICES
!-------------------------------------------------------------
	SUBROUTINE add4dmat2(A, B, a2)
          
        IMPLICIT NONE 
        INTEGER i,j
        DOUBLE PRECISION A(4, 4), B(4, 4)
        DOUBLE PRECISION a2(4, 4)
        
        loopafi: DO i = 1, 4
        loopafj: DO j = 1, 4

        a2(i,j) = 0.d0   

        a2(i,j)= A(i,j) + B(i,j) 

      ENDDO loopafj
      ENDDO loopafi
      
      RETURN

	END SUBROUTINE add4dmat2


!22------------------------------------------------------------
!PURPOSE:  Subroutine for multiplying a scalar to a 2x2 matrix
!--------------------------------------------------------------
      SUBROUTINE scmul2d(A,sc,B)
      IMPLICIT NONE
      DOUBLE PRECISION A(2,2)
      DOUBLE PRECISION B(2,2)
      DOUBLE PRECISION sc
      INTEGER i,j
        loopsi: DO i=1,2
        loopsj:      DO j=1,2
                        B(i,j)=0.d0
                        B(i,j)= sc*A(i,j)
                    ENDDO loopsj
                ENDDO loopsi
                               
       RETURN
       END SUBROUTINE scmul2d

!23-------------------------------------------------------------
!PURPOSE: SUBROUTINE TO COMPUTE THE SUM OF THREE (2X2) MATRICES
!---------------------------------------------------------------
      SUBROUTINE add2dmat(A, B, C, a2)
      
      IMPLICIT NONE 
      INTEGER i,j
      DOUBLE PRECISION A(2, 2), C(2, 2)
      DOUBLE PRECISION B(2, 2)
      DOUBLE PRECISION a2(2, 2)
      
      loopafi: DO i = 1, 2
      loopafj: DO j = 1, 2

      a2(i,j) = 0.d0   

      a2(i,j)= A(i,j) + B(i,j) + C(i,j)

      ENDDO loopafj
      ENDDO loopafi
      
      RETURN

      END SUBROUTINE add2dmat

!24-----------------------------------------------------------
!PURPOSE: SUBROUTINE TO COMPUTE THE SUM OF TWO (2X2) MATRICES
!-------------------------------------------------------------
      SUBROUTINE add2dmat2(A, B, b2)
      
      IMPLICIT NONE 
      INTEGER i,j
      DOUBLE PRECISION A(2, 2), B(2, 2)
      DOUBLE PRECISION b2(2, 2)
      
      loopafi: DO i = 1, 2
      loopafj: DO j = 1, 2

      b2(i,j) = 0.d0   

      b2(i,j)= A(i,j) + B(i,j) 

      ENDDO loopafj
      ENDDO loopafi
      
      RETURN

      END SUBROUTINE add2dmat2

!==========================================================================
      DOUBLE PRECISION FUNCTION DDILOG(X)

      DOUBLE PRECISION X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      DOUBLE PRECISION C(0:18),H,ALFA,B0,B1,B2

      DATA ZERO /0.0D0/, ONE /1.0D0/
      DATA HALF /0.5D0/, MALF /-0.5D0/, MONE /-1.0D0/, MTWO /-2.0D0/
      DATA PI3 /3.28986 81336 96453D0/, PI6 /1.64493 40668 48226D0/

      DATA C( 0) / 0.42996 69356 08137 0D0/
      DATA C( 1) / 0.40975 98753 30771 1D0/
      DATA C( 2) /-0.01858 84366 50146 0D0/
      DATA C( 3) / 0.00145 75108 40622 7D0/
      DATA C( 4) /-0.00014 30418 44423 4D0/
      DATA C( 5) / 0.00001 58841 55418 8D0/
      DATA C( 6) /-0.00000 19078 49593 9D0/
      DATA C( 7) / 0.00000 02419 51808 5D0/
      DATA C( 8) /-0.00000 00319 33412 7D0/
      DATA C( 9) / 0.00000 00043 45450 6D0/
      DATA C(10) /-0.00000 00006 05784 8D0/
      DATA C(11) / 0.00000 00000 86121 0D0/
      DATA C(12) /-0.00000 00000 12443 3D0/
      DATA C(13) / 0.00000 00000 01822 6D0/
      DATA C(14) /-0.00000 00000 00270 1D0/
      DATA C(15) / 0.00000 00000 00040 4D0/
      DATA C(16) /-0.00000 00000 00006 1D0/
      DATA C(17) / 0.00000 00000 00000 9D0/
      DATA C(18) /-0.00000 00000 00000 1D0/

      IF(X .EQ. ONE) THEN
       DDILOG=PI6
       RETURN
      ELSE IF(X .EQ. MONE) THEN
       DDILOG=MALF*PI6
       RETURN
      END IF
      T=-X
      IF(T .LE. MTWO) THEN
       Y=MONE/(ONE+T)
       S=ONE
       A=-PI3+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
      ELSE IF(T .LT. MONE) THEN
       Y=MONE-T
       S=MONE
       A=LOG(-T)
       A=-PI6+A*(A+LOG(ONE+ONE/T))
      ELSE IF(T .LE. MALF) THEN
       Y=(MONE-T)/T
       S=ONE
       A=LOG(-T)
       A=-PI6+A*(MALF*A+LOG(ONE+T))
      ELSE IF(T .LT. ZERO) THEN
       Y=-T/(ONE+T)
       S=MONE
       A=HALF*LOG(ONE+T)**2
      ELSE IF(T .LE. ONE) THEN
       Y=T
       S=ONE
       A=ZERO
      ELSE
       Y=ONE/T
       S=MONE
       A=PI6+HALF*LOG(T)**2
      END IF

      H=Y+Y-ONE
      ALFA=H+H
      B1=ZERO
      B2=ZERO
      DO 1 I = 18,0,-1
      B0=C(I)+ALFA*B1-B2
      B2=B1
    1 B1=B0
      DDILOG=-(S*(B0-H*B2)+A)
      RETURN
      END
!------------------------------
!     Newton-Raphson Method
!------------------------------------------------
      FUNCTION rt(alphas,x1,x2,xacc)

      implicit none

      INTEGER j
      REAL*8 rt,x1,x2,xacc
      REAL*8 df,dx,f,alphas,pi

!      print*, "-----"

      pi = 4.d0 * datan(1.d0)

      rt = 0.5d0 * (x1+x2)

      l1: do j = 1, 500
      
      f = (12.d0 * pi * ( 1.d0 + 
     $     (121104.d0 *
     $     (-0.32233865107676046d0 + 
     $     (-0.5d0 + dLog(rt))**2.d0))/
     $     (279841.d0 * rt**2.d0) - 
     $     (348.d0 * dLog(rt))/(529.d0 * rt)))/(23.d0 * rt) - 
     $     alphas
      
      df = (-0.5553967575594119d0 + 
     $     (-1.0782683423515227d0 - 
     $     1.6390918192642399d0 * rt) * rt + 
     $     (3.5466671374133263d0 + 
     $     2.1565366847030454d0 * rt - 
     $     2.1280002824479958d0 * dLog(rt)) * dLog(rt))/rt**4.d0

      dx = f/df
      
      rt = rt - dx
      
      if(((x1 - rt)*(rt - x2)).lt.0.d0)then
         print*,'rt jumped out of the range'
         rt = 0.d0
         return
      endif
      
      if(abs(dx).lt.xacc) return

      enddo l1

      print*, 'rt exceeded maximum iterations'

      END function rt
!=====================================================      
