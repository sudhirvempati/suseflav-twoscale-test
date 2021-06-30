!------------------------------------------------------------------------
!                              MSSMRGE BEGINS
!------------------------------------------------------------------------
****f* susyflav/mssmrge.f 
*  NAME
*    mssmrge
*  SYNOPSIS
*     In this subroutine we write all the mssm renormalization group equations, 
*     including the flavor mixing in the squark and slepton sector. 
*  FUNCTION
*     Computes the numerical values of all the mssm rge at a given energy scale
*  INPUTS
*     yy(126)    - Initial values for all RGEs
*     t          - energy scale 
*  RESULT
*     yy(126)    - RGE output at a scale t
*  EXAMPLE
*     subroutine mssmrge(t,yy,dydx)
*  NOTES
*     The notation we use closely follows that
*     of Martin and Vaughn PRD50(1994)2282 and 
*     Ibarra and Simonetto JHEP04(2008)102.
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
*	dydx(28)-dydx(30) : ynu(1,1) - ynu(1,3) 
*	dydx(31)-dydx(33) : ynu(2,1) - ynu(2,3) 
*	dydx(34)-dydx(36) : ynu(3,1) - ynu(3,3) 
*
*	dydx(37)-dydx(39) : au(1,1) - au(1,3) 
*	dydx(40)-dydx(42) : au(2,1) - au(2,3) 
*	dydx(43)-dydx(45) : au(3,1) - au(3,3) 
*
*	dydx(46)-dydx(48) : ad(1,1) - ad(1,3) 
*	dydx(49)-dydx(51) : ad(2,1) - ad(2,3) 
*	dydx(52)-dydx(54) : ad(3,1) - ad(3,3) 
*
*	dydx(55)-dydx(57) : ae(1,1) - ae(1,3) 
*	dydx(58)-dydx(60) : ae(2,1) - ae(2,3) 
*	dydx(61)-dydx(63) : ae(3,1) - ae(3,3) 
*	dydx(64)-dydx(66) : anu(1,1) - anu(1,3) 
*	dydx(67)-dydx(69) : anu(2,1) - anu(2,3) 
*	dydx(70)-dydx(72) : anu(3,1) - anu(3,3) 
*	dydx(73)-dydx(75) : mq(1,1) - mq(1,3)
*	dydx(76)-dydx(77) : mq(2,2) - mq(2,3)
*	dydx(78)          : mq(3,3)
*
*       dydx(79)-dydx(81) : mu(1,1) - mu(1,3)
*       dydx(82)-dydx(83) : mu(2,2) - mu(2,3)
*       dydx(84)          : mu(3,3)  
*
*       dydx(85)-dydx(87) : md(1,1) - md(1,3)
*       dydx(88)-dydx(89) : md(2,2) - md(2,3)
*       dydx(90)          : md(3,3)  
*
*       dydx(91)-dydx(93) : ml(1,1) - ml(1,3)
*       dydx(94)-dydx(95) : ml(2,2) - ml(2,3)
*       dydx(96)          : ml(3,3)  
*
*	dydx(97)-dydx(99)   : me(1,1) - me(1,3)
*	dydx(100)-dydx(101) : me(2,2) - me(2,3)
*	dydx(102)           : me(3,3) 
*	dydx(103)-dydx(105) : mnu(1,1) - mnu(1,3)
*	dydx(106)-dydx(107) : mnu(2,2) - mnu(2,3)
*	dydx(108)           : mnu(3,3)
*
*
*       dydx(109)          : mh1
*	dydx(110)          : mh2
*	dydx(111)          : mu
*       dydx(112)          : b_mu 
*
*	dydx(113)-dydx(118) : m_neutrino(3x3)_symmetric
*       dydx(119)-dydx(121) : aplh3-alph1 
*       dydx(122)-dydx(124) : m1t-m3t
*
*  BUGS
*    ---
*  SEE ALSO
*    smrge
******
* You can use this space for remarks that should not be included
* in the documentation.
C****C

c=================================================================================

      subroutine mssmrge(t,yy,dydx)


      implicit none



      integer i, j,k, i0,j0,c,lopt,rc,itcount,rhn,fuscale
      real*8 a1_unif, a2_unif
      double precision t,r,e
      double precision m1t,m2t,m3t,b1,b2,b3
      double precision yuyudag(3,3),ydyddag(3,3),yeyedag(3,3),ad(3,3)
      double precision tryuyudag, trydyddag,tryeyedag,yu(3,3),yd(3,3)
      double precision ye(3,3),mx,gut,gdt,get
      double precision dydx(126), yudagyu(3,3),tryudagyu,yddagyd(3,3)
      double precision tryddagyd, yedagye(3,3),tryedagye,au(3,3)
      double precision ae(3,3), auyudag(3,3),trauyudag, adyddag(3,3)
      double precision tradyddag, aeyedag(3,3),traeyedag,gaut,gadt,gaet
      double precision yddagyu(3,3),tryddagyu,ydyudag(3,3),trydyudag
      double precision yudagyd(3,3),tryudagyd,yuyddag(3,3),tryuyddag
      double precision gqsq,gusq,gdsq,mqsq(3,3),musq(3,3),mdsq(3,3)
      double precision yumuyudag(3,3),ydmdyddag(3,3),mh1sq,mh2sq
      double precision adaddag(3,3),auaudag(3,3),audagau(3,3)
      double precision aeaedag(3,3),aedagae(3,3),glsq,gesq
      double precision yddagmqyd(3,3), yedagmlye(3,3),yemlyedag(3,3)
      double precision mlsq(3,3),mesq(3,3),yedagmeye(3,3)
      double precision trauaudag,tradaddag,traeaedag,traudagau
      double precision traedagae,trydmdyddag,tryemlyedag,trydmqyddag
      double precision tryumuyudag,tryumqyudag,ydmqyddag(3,3)
      double precision tryemeyedag,mmusq,yuaudag(3,3),tryuaudag
      double precision trydaddag,yeaedag(3,3),tryeaedag,yemeyedag(3,3)
      double precision addagad(3,3),yudagmqyu(3,3),traddagad
      double precision yumqyudag(3,3), ydaddag(3,3),yy(126),bmu
      double precision trynuynudag,ynuynudag(3,3), ynu(3,3),gnut
      double precision anu(3,3),trynuanudag,ynuanudag(3,3),ganut
      double precision anuynudag(3,3),ynudagynu(3,3),yedagynu(3,3)
      double precision tranuynudag,trynudagynu,tryedagynu
      double precision mnusq(3,3),ynumnuynudag(3,3),trynumnuynudag
      double precision anuanudag(3,3),tranuanudag,anudaganu(3,3)
      double precision tranudaganu,ynudagmlynu(3,3),gh1sq,gh2sq
      double precision ynumlynudag(3,3),trynumlynudag
      double precision trmqsq,trmdsq,trmusq,trmlsq,trmesq,s
      double precision vev1,vev2    

      double precision addag(3,3),aedag(3,3),anudag(3,3),yddag(3,3),
     $     yub1(3,3)
      double precision audag(3,3),yedag(3,3),ynudag(3,3),yudag(3,3),
     $     ydb1(3,3)
      double precision tryddagmqyd,tryedagmeye,tryedagmlye,trynudagmlynu
      double precision tryudagmqyu,ynub1(3,3),yeb1(3,3),ydb2(3,3),
     $     yeb2(3,3)
           
      double precision mnu_light(3,3)
      
      double precision yuyudagyu(3,3),yddagydyu(3,3)
      double precision ydyddagydyddag(3,3),yuyddagydyudag(3,3),
     $     yeyedagyeyedag(3,3)
      double precision ynuyedagyeynudag(3,3),yudagyuyddagyd(3,3)
      double precision yddagydyddagyd(3,3), yudagyuyudagyu(3,3)
      
      double precision yddagydyudagyu(3,3),yedagyeyedagye(3,3),yub2(3,3)
      double precision ynuynudagynuynudag(3,3)
      double precision yuyudagyuyudag(3,3),ynub2(3,3)
      double precision yedagyeynudagynu(3,3),ynudagynuynudagynu(3,3)

      double precision yddagydyd(3,3),yedagyeye(3,3),yedagyeynu(3,3)
      double precision yudagyuyd(3,3),ynudagynuynu(3,3)
      double precision yddagad(3,3),yedagae(3,3),ynudaganu(3,3),
     $     yudagau(3,3)
      double precision tryddagad,tryedagae,tryudagau,trynudaganu 
      
      double precision adyudagyu(3,3),adyddagyd(3,3),ydyudagau(3,3),
     $     ydyddagad(3,3)
      double precision ydyudagyu(3,3),auyudagyu(3,3),auyddagyd(3,3),
     $     yuyudagau(3,3)
      double precision yuyddagad(3,3),yuyddagyd(3,3),aeyedagye(3,3),
     $     yeyedagae(3,3)
      double precision yeyedagye(3,3),anuynudagynu(3,3),anuyedagye(3,3)
      double precision ynuynudaganu(3,3),ynuyedagae(3,3)
      double precision ynuyedagye(3,3),ydyddagyd(3,3),ynuynudagynu(3,3)
      double precision b1au(3,3),b1ad(3,3),b1ae(3,3),b1anu(3,3) 

      double precision adyddagydyddag(3,3),aeyedagyeyedag(3,3)
      double precision adyudagyuyddag(3,3),anuyedagyeynudag(3,3)
      double precision anuynudagynuynudag(3,3),auyudagyuyudag(3,3) 
      double precision aeynudagynuyedag(3,3),auyddagydyudag(3,3)

      double precision audagauyudagyu(3,3),traudagauyudagyu
      double precision audagyuyudagau(3,3),traudagyuyudagau
      double precision  addagadyudagyu(3,3),traddagadyudagyu
      double precision  yddagydaudagau(3,3),tryddagydaudagau
      double precision  addagydyudagau(3,3),traddagydyudagau
      double precision  yddagadaudagyu(3,3),tryddagadaudagyu
      double precision  anudaganuynudagynu(3,3),tranudaganuynudagynu,
     $     id(3,3)
      
      double precision  anudagynuynudaganu(3,3),tranudagynuynudaganu
      double precision  aedagaeynudagynu(3,3),traedagaeynudagynu
      double precision  yedagaeanudaganu(3,3),tryedagaeanudaganu
      double precision  aedagyeynudaganu(3,3),traedagyeynudaganu
      double precision  yedagaeanudagynu(3,3),tryedagaeanudagynu            

      
      double precision beta1mq(9),beta2mq(9),beta1md(9),beta2md(9),
     $     beta1mu(9)
      double precision beta2mu(9),beta1ml(9),beta2ml(9),beta1me(9),
     $     beta2me(9)
      double precision beta1mnu(9),beta2mnu(9),beta1mh1sq,beta2mh1sq
      double precision beta1mh2sq,beta2mh2sq,anuyddagydyudag(3,3)
      double precision ada21(3,3),ada22(3,3),aea21(3,3),aea22(3,3),
     $     anua21(3,3)
      double precision anua22(3,3),aua21(3,3),aua22(3,3)

      double precision mlynudagynuyedagye(3,3),trmlynudagynuyedagye
      double precision mlynudagynuynudagynu(3,3),trmlynudagynuynudagynu
      double precision mqyddagyd(3,3),trmqyddagyd,aedagyeyedagae(3,3)
      double precision traedagyeyedagae,audagyu(3,3),traudagyu,
     $     audagyuyddagad(3,3)
      double precision traudagyuyddagad,mlyedagye(3,3),trmlyedagye,
     $     tranudagynu
      double precision mlyedagyeyedagye(3,3),trmlyedagyeyedagye,
     $     anudagynu(3,3)
      double precision anudaganuyedagye(3,3),tranudaganuyedagye
      double precision anudagynuyedagae(3,3),tranudagynuyedagae
      double precision audagauyddagyd(3,3),traudagauyddagyd
      double precision mqyddagydyddagyd(3,3),trmqyddagydyddagyd
      double precision addagydyddagad(3,3),traddagydyddagad,addagyd(3,3)
      double precision addagadyddagyd(3,3),traddagadyddagyd,traddagyd,
     $     trydagad
      double precision aedagaeyedagye(3,3),traedagaeyedagye,ydagad(3,3)
      double precision yddagmdyd(3,3),tryddagmdyd,aedagye(3,3),traedagye
      double precision mqyudagyuyddagyd(3,3),trmqyudagyuyddagyd,
     $     mlynudagynu(3,3)
      double precision trmlynudagynu,tryedagmeyeyedagye,tryudagyuyddagyd

      double precision mqyudagyu(3,3),trmqyudagyu
      double precision mqyudagyuyudagyu(3,3),trmqyudagyuyudagyu
      double precision yddagmdydyddagyd(3,3),tryddagmdydyddagyd
      double precision tryddagydyddagyd
      double precision yedagmeyeyedagye(3,3)
      double precision yedagyeanudaganu(3,3),tryedagyeanudaganu
      double precision tryedagyeyedagye
      double precision ynudaganuaedagye(3,3),trynudaganuaedagye        
      double precision ynudagmnuynu(3,3),trynudagmnuynu
      double precision ynudagmnuynuyedagye(3,3),trynudagmnuynuyedagye
      double precision ynudagmnuynuynudagynu(3,3),
     $     trynudagmnuynuynudagynu
      double precision ynudagynuaedagae(3,3),trynudagynuaedagae
      double precision ynudagynuyedagmeye(3,3),trynudagynuyedagmeye   
      double precision ynudagynumlyedagye(3,3),trynudagynumlyedagye
      double precision ynudagynuyedagye(3,3),trynudagynuyedagye
      double precision trynudagynuynudagynu
      double precision yudagauaddagyd(3,3),tryudagauaddagyd
      double precision yudagmuyu(3,3),tryudagmuyu
      double precision yudagmuyuyddagyd(3,3),tryudagmuyuyddagyd
      double precision yudagmuyuyudagyu(3,3),tryudagmuyuyudagyu
      double precision yudagyuaddagad(3,3),tryudagyuaddagad
      double precision yudagyumqyddagyd(3,3),tryudagyumqyddagyd
      double precision yudagyuyddagmdyd(3,3),tryudagyuyddagmdyd
      
      double precision yddagydmq(3,3),yudagyuyudagau(3,3),
     $     yudagyuyddagad(3,3)
      double precision yudagauyddagyd(3,3),ynuyedag(3,3),
     $     ynudagynuynudaganu(3,3)
      double precision ynudagynuyedagae(3,3),ynudaganuynudagynu(3,3)
      double precision yeynudag(3,3),yedagyeynudaganu(3,3),
     $     yedagyeyedagae(3,3)
      double precision yedagaeynudagynu(3,3),yedagaeyedagye(3,3)

      double precision yddagydyddagad(3,3),yddagadyudagyu(3,3),
     $     yddagadyddagyd(3,3)
      double precision yudagauyudagyu(3,3),ynudaganuyedagye(3,3)
      double precision yddagydyudagau(3,3)
      double precision yudagyumq(3,3),yuyudagmu(3,3),ydyddagmd(3,3)
      double precision yedagyeml(3,3),yeyedagme(3,3),ynuynudagmnu(3,3)
      double precision tryudagyuyudagyu
      double precision yddagadaddagyd(3,3),muyudagyu(3,3),mdydyddag(3,3)
      double precision ynudagynuml(3,3),yedagaeaedagye(3,3)
      double precision meyeyedag(3,3),mnuynuynudag(3,3)
      double precision ydaudagauyddag(3,3)
      double precision yedagyeaedagae(3,3),yeanudaganuyedag(3,3)
      double precision ynuyedagyeynudagmnu(3,3)
      double precision ynudaganuanudagynu(3,3),ynuyedagmeyeynudag(3,3)
      double precision aeynudagynuaedag(3,3),ynuyedagyemlynudag(3,3)
      double precision yddagydaddagad(3,3),adyudagyuaddag(3,3)
      double precision ynudagynuanudaganu(3,3),yeynudaganuaedag(3,3)
      
      double precision yudagauaudagyu(3,3),ydyudagauaddag(3,3)
      double precision yedagyeyedagyeml(3,3),aeanudagynuyedag(3,3) 
      double precision ynumlyedagyeynudag(3,3)
      double precision yudagyuaudagau(3,3),adaudagyuyddag(3,3)
      double precision yedagyeyedagmeye(3,3),yeaedagaeyedag(3,3) 
      double precision mnuynuyedagyeynudag(3,3)
      double precision yddagydyddagydmq(3,3),ydaddagadyddag(3,3)
      double precision yedagyemlyedagye(3,3),aeyedagyeaedag(3,3)
      double precision ynuynudagynuynudagmnu(3,3)
      double precision yddagydyddagmdyd(3,3),adyddagydaddag(3,3)
      double precision ynudagynuynudagynumq(3,3),yeyedagaeaedag(3,3)
      double precision ynuynudagynumlynudag(3,3)
      double precision yddagydmqyddagyd(3,3),adaddagydyddag(3,3)
      double precision aeaedagyeyedag(3,3),yuaudagauyudag(3,3)
      double precision ynuynudagmnuynuynudag(3,3)
      double precision yudagyuyudagyumq(3,3),auyudagyuaudag(3,3)
      double precision ydyudagyuyddagmd(3,3),auaudagyuyudag(3,3)
      double precision ynumlynudagynuynudag(3,3)
      double precision ynudagynuynudagmnuynu(3,3),ydyudagyumqyddag(3,3)
      double precision yeynudagynuyedagme(3,3),yudagyuyudagmuyu(3,3)
      double precision mnuynuynudagynuynudag(3,3),yuyudagauaudag(3,3)
      double precision yudagyumqyudagyu(3,3),ydyudagmuyuyddag(3,3)
      double precision ynudagynumlynudagynu(3,3),yeynudagynumlyedag(3,3)
      double precision ydmqyudagyuyddag(3,3),yeynudagmnuynuyedag(3,3)
      double precision ydyudagyuyddag(3,3),yemlynudagynuyedag(3,3)
      double precision yeynudagynuyedag(3,3),mdydyudagyuyddag(3,3)

      double precision ydyddagydyddagmd(3,3),meyeynudagynuyedag(3,3)
      double precision ydyddagydmqyddag(3,3),yeyedagyeyedagme(3,3)
      double precision yeyedagyemlyedag(3,3),ydyddagmdydyddag(3,3)
      double precision ydmqyddagydyddag(3,3),yeyedagmeyeyedag(3,3)
      double precision mdydyddagydyddag(3,3),yemlyedagyeyedag(3,3)
      double precision meyeyedagyeyedag(3,3),yuyddagydyudagmu(3,3)
      double precision yuyddagydmqyudag(3,3),yuyddagmdydyudag(3,3)
      double precision yumqyddagydyudag(3,3),muyuyddagydyudag(3,3)
      double precision yuyudagyuyudagmu(3,3),yuyudagyumqyudag(3,3)
      double precision yuyudagmuyuyudag(3,3),yumqyudagyuyudag(3,3)
      double precision muyuyudagyuyudag(3,3),muyuyudag(3,3),
     $     yuaddagadyudag(3,3)
      double precision auyddagydaudag(3,3),yuyddagadaudag(3,3),
     $     auaddagydyudag(3,3)
      double precision yedagad(3,3),tryedagad,ynuaedagaeynudag(3,3)
      double precision anuyedagyeanudag(3,3),ynuyedagaeanudag(3,3)
      double precision anuaedagyeynudag(3,3),ynuanudaganuynudag(3,3)
      double precision anuynudagynuanudag(3,3),ynuynudaganuanudag(3,3) 
      double precision anuanudagynuynudag(3,3),yeynudagynu(3,3)
      double precision sd,sig1,sig2,sig3
      double precision ydyudagyuyddagyd(3,3),trydyddagydyddag
      double precision ydyudagyuyudagyu(3,3),ydyddagydyddagyd(3,3)
      double precision tryuyddagydyudag
      double precision tryeyedagyeyedag,trynuyedagyeynudag
      double precision trmnusq,ynudagynuynudagynuml(3,3)
      DOUBLE PRECISION tradyddagydyddag
      DOUBLE PRECISION tradyudagyuyddag
      DOUBLE PRECISION traeyedagyeyedag
      DOUBLE PRECISION traeynudagynuyedag
      DOUBLE PRECISION tranuyedagyeynudag
      DOUBLE PRECISION tranuynudagynuynudag
      DOUBLE PRECISION trauyddagydyudag,trauyudagyuyudag
      DOUBLE PRECISION trynuynudagynuynudag,tryuyudagyuyudag
      DOUBLE PRECISION beta1yu(9),beta2yu(9),beta1yd(9),beta2yd(9)
      DOUBLE PRECISION beta1ye(9),beta2ye(9),beta1ynu(9),beta2ynu(9)
      DOUBLE PRECISION beta1au(9),beta2au(9),beta1ad(9),beta2ad(9)
      DOUBLE PRECISION beta1ae(9),beta2ae(9),beta1anu(9),beta2anu(9)
      
      double precision alph3,alph2,alph1 !,a1,a2,a3
      double precision pi
      real a1,a2,a3

      double precision yuyudagydyddag(3,3),ydyddagyuyudag(3,3)
      double precision trydyddagyuyudag,tryuyudagydyddag,yukgut(126)


      DOUBLE PRECISION tol,tol1, e1, e_next

      data tryudagyu/0.d0/, trmqyudagyu/0.d0/, tryudagmuyu/0.d0/
      data traudagau/ 0.d0/, trynudagynu/0.d0/,trmlynudagynu/0.d0/
      data trynudagmnuynu/0.d0/, tranudaganu/0.d0/
      DATA yukgut/ 126 * 0.d0/

      character*4 model
      character*3 case
      common/charinputs/case, model
      common/loops/ lopt, rhn
      common/counter/itcount
      common/unif/ e1,yukgut
      common/unifs/fuscale


!---------------------------------------------------
      
      
      pi = datan(1.d0) * 4.d0

      b1 = 33.d0/5.d0 
      b2 = 1.d0
      b3 = -3.d0

      MX  = 5.d0*(10.d0**19.d0)

!------------------------initializing matrices to zero

      yu(1,1) = 0.d0
      yu(1,2) = 0.d0
      yu(1,3) = 0.d0
      yu(2,1) = 0.d0
      yu(2,2) = 0.d0
      yu(2,3) = 0.d0
      yu(3,1) = 0.d0
      yu(3,2) = 0.d0

      yd(1,1) = 0.d0
      yd(1,2) = 0.d0
      yd(1,3) = 0.d0
      yd(2,1) = 0.d0
      yd(2,2) = 0.d0
      yd(2,3) = 0.d0
      yd(3,1) = 0.d0
      yd(3,2) = 0.d0


      ye(1,1) = 0.d0
      ye(1,2) = 0.d0
      ye(1,3) = 0.d0
      ye(2,1) = 0.d0
      ye(2,2) = 0.d0
      ye(2,3) = 0.d0
      ye(3,1) = 0.d0
      ye(3,2) = 0.d0


   
      au(1,2) = 0.d0
      au(1,3) = 0.d0
      au(2,1) = 0.d0
      au(2,3) = 0.d0
      au(3,1) = 0.d0
      au(3,2) = 0.d0

   
      ad(1,2) = 0.d0
      ad(1,3) = 0.d0
      ad(2,1) = 0.d0
      ad(2,3) = 0.d0
      ad(3,1) = 0.d0
      ad(3,2) = 0.d0

   
      ae(1,2) = 0.d0
      ae(1,3) = 0.d0
      ae(2,1) = 0.d0
      ae(2,3) = 0.d0
      ae(3,1) = 0.d0
      ae(3,2) = 0.d0


      mqsq(1,2) = 0.d0
      mqsq(1,3) = 0.d0
      mqsq(2,1) = 0.d0
      mqsq(2,3) = 0.d0
      mqsq(3,1) = 0.d0
      mqsq(3,2) = 0.d0


      musq(1,2) = 0.d0
      musq(1,3) = 0.d0
      musq(2,1) = 0.d0
      musq(2,3) = 0.d0
      musq(3,1) = 0.d0
      musq(3,2) = 0.d0


      mdsq(1,2) = 0.d0
      mdsq(1,3) = 0.d0
      mdsq(2,1) = 0.d0
      mdsq(2,3) = 0.d0
      mdsq(3,1) = 0.d0
      mdsq(3,2) = 0.d0


      mlsq(1,2) = 0.d0
      mlsq(1,3) = 0.d0
      mlsq(2,1) = 0.d0
      mlsq(2,3) = 0.d0
      mlsq(3,1) = 0.d0
      mlsq(3,2) = 0.d0


      mesq(1,2) = 0.d0
      mesq(1,3) = 0.d0
      mesq(2,1) = 0.d0
      mesq(2,3) = 0.d0
      mesq(3,1) = 0.d0
      mesq(3,2) = 0.d0

      mh1sq = 0.d0
      mh2sq = 0.d0

!-------------------------------


C     Top Yukawa !!!
C     ----------------------------------

      i0 = 0 

      mssm01: do i = 1,3
      
      yu(1,i) = yy(i0 + i)
      j = 3 + i
      yu(2,i) = yy(i0 + j)
      k = 6 + i
      yu(3,i) = yy(i0 + k)

      enddo mssm01



C     Bottom Yukawa !!!
C     ----------------------------------

      i0 = 9 
      mssm02: do i = 1,3
      
      yd(1,i) = yy(i0 + i)
      j = 3 + i
      yd(2,i) = yy(i0 + j)
      k = 6 + i
      yd(3,i) = yy(i0 + k)
      
      enddo mssm02     


      
C     Tau Yukawa !!!
C     ----------------------------------

      i0 = 18

      mssm03: do i = 1,3
      
      ye(1,i) = yy(i0 + i)
      j = 3 + i
      ye(2,i) = yy(i0 + j)
      k = 6 + i
      ye(3,i) = yy(i0 + k)

      enddo  mssm03
     
      
C     Neutrino Yukawa !!!
C     ----------------------------------

      i0 = 27 

      mssm04: do i = 1,3
      
      ynu(1,i) = yy(i0 + i)
      j = 3 + i
      ynu(2,i) = yy(i0 + j)
      k = 6 + i
      ynu(3,i) = yy(i0 + k)

      enddo  mssm04
      


C     Au-matrix !!!!
C     ------------------
      
!     print*,"before_loop yy(38)",yy(38)

      i0 = 36

      mssm05: do i = 1,3
      
!     print*, "YY(",36+i,")", yy(36+i)

      au(1,i) = yy(i0 + i)
      j = 3 + i
      au(2,i) = yy(i0 + j)
      k = 6 + i
      au(3,i) = yy(i0 + k)

      enddo  mssm05



c$$$  loop52: do i=1,3
c$$$  loop62: do j=1,3
c$$$  !         print*, "au(",i,j,")", au(i,j) 
c$$$  !         print*, "yddag(",i,j,")", yddag(i,j) 
c$$$  C         print*, "YY(",36+i,")", yy(36+i)
c$$$  enddo loop62	
c$$$  enddo loop52

      
      
C     Ad-matrix !!!!
C     --------------- 

      i0 = 45

      mssm06: do i = 1,3
      
      ad(1,i) = yy(i0 + i)
      j = 3 + i
      ad(2,i) = yy(i0 + j)
      k = 6 + i
      ad(3,i) = yy(i0 + k)

      enddo mssm06


C     Ae-matrix !!
C     -----------------

      i0 = 54

      mssm07: do i = 1,3
      
      ae(1,i) = yy(i0 + i)
      j = 3 + i
      ae(2,i) = yy(i0 + j)
      k = 6 + i
      ae(3,i) = yy(i0 + k)

      enddo mssm07
    

C     Anu-matrix !!
C     -----------------

      i0 = 63 

      mssm08: do i = 1,3
      
      anu(1,i) = yy(i0 + i)
      j = 3 + i
      anu(2,i) = yy(i0 + j)
      k = 6 + i
      anu(3,i) = yy(i0 + k)

      enddo mssm08

c$$$  loop51: do i=1,3
c$$$  loop61: do j=1,3
c$$$  print*, "anu(",i,j,")", anu(i,j)
c$$$  
c$$$  enddo loop61	
c$$$  enddo loop51
c$$$  print*,"----------------------"



C     mQ-matrix 
C     ---------------------
      
      i0 = 72

      mssm09: do i = 1,3
      mqsq(1,i) = yy(i0 + i)
      enddo mssm09
      
      j0 = i0 + 3 

      mssm10: do j = 1,2  
      mqsq(2,j+1) = yy(j0 + j)
      enddo mssm10         
      k = i0 + 6 

      mqsq(3,3) = yy(k)

      mqsq(2,1) = mqsq(1,2)
      mqsq(3,1) = mqsq(1,3)
      mqsq(3,2) = mqsq(2,3) 

     
      trmqsq = 0.d0
      
      ltrq: do i = 1,3
      trmqsq = trmqsq + mqsq(i,i)
      enddo ltrq

!     print*, "trmqsq", trmqsq

C     mU-matrix
C-----------------------------------------------

      i0 = 78

      mssm11: do i = 1,3
      musq(1,i) = yy(i0 + i)
      enddo mssm11
      
      j0 = i0 + 3 

      mssm12: do j = 1,2  
      musq(2,j+1) = yy(j0 + j)
      enddo mssm12
      
      k = i0 + 6 

      musq(3,3) = yy(k)

      musq(2,1) = musq(1,2)
      musq(3,1) = musq(1,3)
      musq(3,2) = musq(2,3) 

      trmusq = 0.d0
      
      ltru: do i = 1,3
      trmusq = trmusq + musq(i,i)
      enddo ltru

!     print*, "trmusq", trmusq


C     mD-matrix
C-------------------------------------------------

      i0 = 84

      mssm13: do i = 1,3
      mdsq(1,i) = yy(i0 + i)
      enddo mssm13 
      
      j0 = i0 + 3 

      mssm14: do j = 1,2  
      mdsq(2,j+1) = yy(j0 + j)
      enddo mssm14
      
      k = i0 + 6 

      mdsq(3,3) = yy(k)

      mdsq(2,1) = mdsq(1,2)
      mdsq(3,1) = mdsq(1,3)
      mdsq(3,2) = mdsq(2,3) 



      trmdsq = 0.d0
      
      ltrd: do i = 1,3
      trmdsq = trmdsq + mdsq(i,i)
      enddo ltrd
      
!     print*, "trmdsq", trmdsq


C     mL-matrix
C     ----------------------

      i0 = 90

      mssm15: do i = 1,3
      mlsq(1,i) = yy(i0 + i)
      enddo mssm15
      
      j0 = i0 + 3 

      mssm16: do j = 1,2  
      mlsq(2,j+1) = yy(j0 + j)
      enddo mssm16
      
      k = i0 + 6 

      mLsq(3,3) = yy(k)

      mlsq(2,1) = mlsq(1,2)
      mlsq(3,1) = mlsq(1,3)
      mlsq(3,2) = mlsq(2,3) 



      trmlsq = 0.d0
      
      ltrl: do i = 1,3
      trmlsq = trmlsq + mLsq(i,i)
      enddo ltrl
      
!     print*, "trmlsq", trmlsq
      



C     mE-matrix
C     --------------------------

      i0 = 96

      mssm17: do i = 1,3
      mesq(1,i) = yy(i0 + i)
      enddo mssm17
      
      j0 = i0 + 3 

      mssm18: do j = 1,2  
      mesq(2,j+1) = yy(j0 + j)
      enddo mssm18
      
      k = i0 + 6 

      mesq(3,3) = yy(k)

      mesq(2,1) = mesq(1,2)
      mesq(3,1) = mesq(1,3)
      mesq(3,2) = mesq(2,3)


      trmesq = 0.d0
      
      ltre: do i = 1,3
      trmesq = trmesq + mEsq(i,i)
      enddo ltre

!     print*, "trmesq", trmesq


C     ------------------------------------------------
C     mNU-matrix
C     ------------------------------------------------

      i0 = 102

      mssm19: do i = 1,3
      mnusq(1,i) = 0.d0         !yy(i0 + i)
      enddo mssm19
      
      j0 = i0 + 3 

      mssm20: do j = 1,2  
      mnusq(2,j+1) = 0.d0       !yy(j0 + j)
      enddo mssm20 
      
      k = i0 + 6 

      mnusq(3,3) = 0.d0         !yy(k)

      mnusq(2,1) = mnusq(1,2)
      mnusq(3,1) = mnusq(1,3)
      mnusq(3,2) = mnusq(2,3) 

      call trace(mNUsq,trmnusq)
!        print*, "trmnusq", trmnusq

      
C     ------------------------------------------------
      
!     print*,"inside mssmrge yy(109)",yy(109)
      
      mh1sq = yy(109)
      mh2sq = yy(110)
      mmusq = (yy(111))
      bmu   = yy(112)

c$$$  print*,"A",mh1sq,"B",mh2sq,"C",mmusq,"D",bmu
c$$$  print*,"___________________________________"

C     mnu_light NON RENORM running
C     --------------------------------

      
      i0 = 112
      
      mssm_chank_1: do i = 1,3
      mnu_light(1,i) = yy(i0 + i)
      enddo mssm_chank_1
      
      j0 = i0 + 3
      
      mssm_chank_2: do j = 1,2
      mnu_light(2,j+1) = yy(j0 + j)
      enddo mssm_chank_2
      
      k = i0 + 6
      
      mnu_light(3,3) = yy(k)
      
      mnu_light(2,1) = mnu_light(1,2)
      mnu_light(3,1) = mnu_light(1,3)
      mnu_light(3,2) = mnu_light(2,3)  


C------------------------------------------
C     Gauge Couplings
C-------------------------------------------
      
      alph3 = yy(119)           !! defination, alphi = gi
      alph2 = yy(120)
      alph1 = yy(121)

C     print*,"A",alph3,"B",alph2,"C",alph1

C-----------------------------------------------
C     Gaugino Masses
C-----------------------------------------------

      M1t = yy(122)
      M2t = yy(123)
      M3t = yy(124)

      vev1 = yy(125)
      vev2 = yy(126)
c$$$  
!     print*,"A",M1t,"B",M2t,"C",M3t


C---------------------------------------------------------------------
C     Matrix Multiplications Begin
C---------------------------------------------------------------------

C     ---------------------------------------------------------
C     Y*Dagger[Y] product matrices of the Yukawas 
C     ---------------------------------------------------------

C     Up-sector :
C     ----------
      call dag(yu,yudag)
      call matmult(yu,yudag,yuyudag)    
      call trace(yuyudag,tryuyudag)

!      print*, "tryuyudag", tryuyudag

     


C     Neutrino-sector :
C     ----------
      call dag(ynu,ynudag)
      call matmult(ynu,ynudag,ynuynudag)    
      call trace(ynuynudag,trynuynudag)
     



C     Up-sector (A-part)
C     ------------------
      call dag(au,audag)
      call matmult(au,audag,auaudag)    
      call trace(auaudag,trauaudag)

     

C     Neutrino-sector (Anu-part)
C     ------------------
      call dag(anu,anudag)
      call matmult(anu,anudag,anuanudag)    
      call trace(anuanudag,tranuanudag)
      
      call matmult(anudag,anu,anudaganu)    
      call trace(anudaganu,tranudaganu)
     



C     Combination Y. Dagger[A] 
C     ------------------------
      call dag(au,audag)
      call matmult(yu,audag,yuaudag)    
      call trace(yuaudag,tryuaudag)
      
      

C     Combination Ynu. Dagger[Anu] 
C     -----------------------------
      call dag(anu,anudag)


  
      call matmult(anudag,ynu,anudagynu)    
      call trace(anudagynu,tranudagynu)
      
      call matmult(ynu,anudag,ynuanudag)    
      call trace(ynuanudag,trynuanudag)
!      print*, "trynuanudag= " , tranudagynu    
     

C     Down-sector : 
C     ------------ 
      call dag(yd,yddag)
      call matmult(yd,yddag,ydyddag)  
      call trace(ydyddag,trydyddag)
     

C     Down-sector (A-part)
C     -------------------
      call dag(ad,addag)
      call matmult(ad,addag,adaddag)    
      call trace(adaddag,tradaddag)
     


C     Combination Y^D Dagger[A^D]
C     ---------------------------
      call dag(ad,addag)
      call matmult(yd,addag,ydaddag)    
      call trace(ydaddag,trydaddag)
     

C     Charged Lepton sector : 
C     ----------------------
      call dag(ye,yedag)
      call matmult(ye,yedag,yeyedag)    
      call trace(yeyedag,tryeyedag)
     

C     Charged Lepton-sector (A-part)
C     -----------------------------
      call dag(ae,aedag)
      call matmult(ae,aedag,aeaedag)    
      call trace(aeaedag,traeaedag)
     


C     Combination (Y^E Dagger[A^E])
C     -----------------------------
      call dag(ae,aedag)
      call matmult(ye,aedag,yeaedag)    
      call trace(yeaedag,tryeaedag)

C    ------------------------------
      call matmult(yddag,ad,yddagad)
      call trace(yddagad,tryddagad)

      call matmult(yedag,ae,yedagae)
      call trace(yedagae,tryedagae)
     
      call matmult(yudag,au,yudagau)
      call trace(yudagau,tryudagau)
     
      call matmult(ynudag,anu,ynudaganu)
      call trace(ynudaganu,trynudaganu)
     
C     ---------------------------------------------------------
C     Dagger[Y]*Y product matrices of the Yukawas 
C     ---------------------------------------------------------

C     Up-sector :
C     ----------
      call dag(yu,yudag)
      call matmult(yudag,yu,yudagyu)    
      call trace(yudagyu,tryudagyu)
      
!      print*,"tryudagyu",tryudagyu
  

C     Neutrino-sector :
C     ----------
      call dag(ynu,ynudag)
      call matmult(ynudag,ynu,ynudagynu)    
      call trace(ynudagynu,trynudagynu)
    

C     Dagger[Y^E] Y^NU :
C     -----------------
      call dag(ye,yedag)
      call matmult(yedag,ynu,yedagynu)    
      call trace(yedagynu,tryedagynu)
     


C     Up-sector (A-part):
C     -----------------
      call dag(au,audag)
      
!     print*,"----------------------"
      
      call matmult(audag,au,audagau)    
      call trace(audagau,traudagau)
      
      

C     Down-sector : 
C     ------------ 
      call dag(yd,yddag)
      call matmult(yddag,yd,yddagyd)    
      call trace(yddagyd,tryddagyd)
     



C     Down-sector (A-part)
C     -------------------
      call dag(ad,addag)
      call matmult(addag,ad,addagad)    
      call trace(addagad,traddagad)
     


C     Charged Lepton sector : 
C     ----------------------
      call dag(ye,yedag)
      call matmult(yedag,ye,yedagye)    
      call trace(yedagye,tryedagye)
     

C     Charged-lepton (A-part)
C     ----------------------
      call dag(ae,aedag)
      call matmult(aedag,ae,aedagae)    
      call trace(aedagae,traedagae)
      
      
      
C     ====================================================================
C     Other combination of products of matrices for A-equations.
C     ====================================================================
      
C     A^U Dagger[Y^U] : 
C     ---------------
      call dag(yu,yudag)
      call matmult(au,yudag,auyudag)    
      call trace(auyudag,trauyudag)
     


C     A^NU Dagger[Y^NU] : 
C     ---------------
      call dag(ynu,ynudag)
      call matmult(anu,ynudag,anuynudag)    
      call trace(anuynudag,tranuynudag)
     
    


C     A^D Dagger[Y^D] : 
C     ---------------
      call dag(yd,yddag)
      call matmult(ad,yddag,adyddag)    
      call trace(adyddag,tradyddag)
     
    

C     A^E Dagger[Y^E] : 
C     ---------------
      call dag(ye,yedag)
      call matmult(ae,yedag,aeyedag)    
      call trace(aeyedag,traeyedag)
     

C     Dagger[Y^D] Y^U :
C     --------------- 	
      call dag(yd,yddag)
      call matmult(yddag,yu,yddagyu)    
      call trace(yddagyu,tryddagyu)
     
 

C     Dagger[Y^U] Y^D :
C     --------------- 	
      call dag(yu,yudag)
      call matmult(yudag,yd,yudagyd)    
      call trace(yudagyd,tryudagyd)
     
    


C     Y^D Dagger[Y^U] :
C     --------------- 	
      call dag(yu,yudag)
      call matmult(yd,yudag,ydyudag)    
      call trace(ydyudag,trydyudag)
     
    

C     Y^U Dagger[Y^D] :
C     --------------- 	
      call dag(yd,yddag)
      call matmult(yu,yddag,yuyddag)    
      call trace(yuyddag,tryuyddag)
     
 


C     Y^U M^U Dagger[Y^U] 
C     -------------------
      call mat3prod(yu,musq,yudag,yumuyudag)
      call trace(yumuyudag,tryumuyudag)


C     Y^U M^Q Dagger[Y^U]
C     -------------------
      call mat3prod(yu,mqsq,yudag,yumqyudag)
      call trace(yumqyudag,tryumqyudag)

C     Y^D M^D Dagger[Y^D] 
C     -------------------
      call mat3prod(yd,mdsq,yddag,ydmdyddag)
      call trace(ydmdyddag,trydmdyddag)


C     Y^D M^Q Dagger[Y^D] 
C     -------------------
      call mat3prod(yd,mqsq,yddag,ydmqyddag)
      call trace(ydmqyddag,trydmqyddag)

      



C     Dagger[Y^U] M^Q Y^U
C     ------------------
      call mat3prod(yudag,mqsq,yu,yudagmqyu)
      call trace(yudagmqyu,tryudagmqyu)



C     Dagger[Y^NU] M^L Y^NU
C     ------------------
      call mat3prod(ynudag,mlsq,ynu,ynudagmlynu)
      call trace(ynudagmlynu,trynudagmlynu)

      

C     Dagger[Y^D] M^D Y^D
C     -------------------
      call mat3prod(yddag,mqsq,yd,yddagmqyd)
      call trace(yddagmqyd,tryddagmqyd)

      

C     Y^E M^L Dagger[Y^e] 
C     ---------------------
      call mat3prod(ye,mlsq,yedag,yemlyedag)
      call trace(yemlyedag,tryemlyedag)

      

C     Y^E M^e Dagger[Y^e] 
C     ---------------------
      call mat3prod(ye,mesq,yedag,yemeyedag)
      call trace(yemeyedag,tryemeyedag)

      

C     Y^NU M^NU Dagger[Y^NU] 
C     ---------------------
      call mat3prod(ynu,mnusq,ynudag,ynumnuynudag)
      call trace(ynumnuynudag,trynumnuynudag)

      


C     Y^NU M^L Dagger[Y^NU] 
C     ---------------------
      call mat3prod(ynu,mlsq,ynudag,ynumlynudag)
      call trace(ynumlynudag,trynumlynudag)

      


C     Dagger[Y^E] M^L Y^E
C     -------------------
      call mat3prod(yedag,mlsq,ye,yedagmlye)
      call trace(yedagmlye,tryedagmlye)


C     Dagger[Y^E] M^e Y^E
C     -------------------
      call mat3prod(yedag,mesq,ye,yedagmeye)
      call trace(yedagmeye,tryedagmeye)
      
C------------------------------------------------------------------------------------------
C     Remaining matrix products(product of 3,4,5 (3x3)matrices) and their respective traces
C------------------------------------------------------------------------------------------
     
      call matmult(yedag,ad,yedagad)
      call trace(yedagad,tryedagad)
      call mat3prod(yu,yudag,yu,yuyudagyu)
      call mat3prod(yddag,yd,yu,yddagydyu)
      call mat3prod(yddag,yd,yd,yddagydyd)
      call mat3prod(yedag,ye,ye,yedagyeye)
      call mat3prod(yedag,ye,ynu,yedagyeynu)
      call mat3prod(yudag,yu,yd,yudagyuyd)
      call mat3prod(ynudag,ynu,ynu,ynudagynuynu)
      call mat3prod(ad,yudag,yu,adyudagyu)
      call mat3prod(ad,yddag,yd,adyddagyd)
      call mat3prod(yd,yudag,au,ydyudagau)
      call mat3prod(yd,yddag,ad,ydyddagad)
      call mat3prod(yd,yddag,yd,ydyddagyd)
      call mat3prod(yd,yudag,yu,ydyudagyu)
      call mat3prod(au,yudag,yu,auyudagyu)
      call mat3prod(au,yddag,yd,auyddagyd)
      call mat3prod(yu,yudag,au,yuyudagau)
      call mat3prod(yu,yddag,ad,yuyddagad)
      call mat3prod(yu,yddag,yd,yuyddagyd)
      call mat3prod(ae,yedag,ye,aeyedagye)
      call mat3prod(ye,yedag,ae,yeyedagae)
      call mat3prod(ye,yedag,ye,yeyedagye)
      call mat3prod(anu,ynudag,ynu,anuynudagynu)
      call mat3prod(anu,yedag,ye,anuyedagye)
      call mat3prod(ynu,ynudag,anu,ynuynudaganu)
      call mat3prod(ynu,yedag,ae,ynuyedagae)
      call mat3prod(ynu,ynudag,ynu,ynuynudagynu)
      call mat3prod(ynu,yedag,ye,ynuyedagye)
      call mat3prod(yu,yudag,yu,yuyudagyu)
      call mat3prod(yd,yddag,yd,ydyddagyd)
      call mat3prod(ye,yedag,ye,yeyedagye)
      
     
      call mat4pr(yd,yddag,yd,yddag,ydyddagydyddag)
      call mat4pr(yu,yddag,yd,yudag,yuyddagydyudag)
      call mat4pr(ye,yedag,ye,yedag,yeyedagyeyedag)
      call mat4pr(ynu,yedag,ye,ynudag,ynuyedagyeynudag)
      call mat4pr(yddag,yd,yddag,yd,yddagydyddagyd)
      call mat4pr(yudag,yu,yudag,yu,yudagyuyudagyu)
      call mat4pr(yudag,yu,yddag,yd,yudagyuyddagyd)

      
      call mat4pr(yu,yudag,yu,yudag,yuyudagyuyudag)
      call mat4pr(ynu,ynudag,ynu,ynudag,ynuynudagynuynudag)
      call mat4pr(yddag,yd,yudag,yu,yddagydyudagyu)
      call mat4pr(yedag,ye,ynudag,ynu,yedagyeynudagynu)
      call trace(ynuynudagynuynudag,trynuynudagynuynudag)
      call trace(yuyudagyuyudag,tryuyudagyuyudag)



      call mat4pr(ad,yddag,yd,yddag,adyddagydyddag)
      call mat4pr(ad,yudag,yu,yddag,adyudagyuyddag)
      call mat4pr(ae,yedag,ye,yedag,aeyedagyeyedag)
      call mat4pr(ae,ynudag,ynu,yedag,aeynudagynuyedag)
      call mat4pr(anu,yedag,ye,ynudag,anuyedagyeynudag)
      call mat4pr(au,yddag,yd,yudag,anuyddagydyudag)
      call mat4pr(anu,ynudag,ynu,ynudag,anuynudagynuynudag)
      call mat4pr(au,yudag,yu,yudag,auyudagyuyudag)
      call mat4pr(ynu,aedag,ae,ynudag,ynuaedagaeynudag)
    
      call mat4pr(audag,au,yudag,yu,audagauyudagyu)
      call trace(audagauyudagyu,traudagauyudagyu)

      call mat4pr(audag,yu,yudag,au,audagyuyudagau)
      call trace(audagyuyudagau,traudagyuyudagau)

      call mat4pr(addag,ad,yudag,yu,addagadyudagyu)
      call trace(addagadyudagyu,traddagadyudagyu)

      call mat4pr(yddag,yd,audag,au,yddagydaudagau)
      call trace(yddagydaudagau,tryddagydaudagau)

      call mat4pr(addag,yd,yudag,au,addagydyudagau)
      call trace(addagydyudagau,traddagydyudagau)

      call mat4pr(yddag,ad,audag,yu,yddagadaudagyu)
      call trace(yddagadaudagyu,tryddagadaudagyu)

      call mat4pr(anudag,anu,ynudag,ynu,anudaganuynudagynu)
      call trace(anudaganuynudagynu,tranudaganuynudagynu)

      call mat4pr(anudag,ynu,ynudag,anu,anudagynuynudaganu)
      call trace(anudagynuynudaganu,tranudagynuynudaganu)

      call mat4pr(aedag,ae,ynudag,ynu,aedagaeynudagynu)
      call trace(aedagaeynudagynu,traedagaeynudagynu)

      call mat4pr(yedag,ae,anudag,anu,yedagaeanudaganu)
      call trace(yedagaeanudaganu,tryedagaeanudaganu)

      call mat4pr(aedag,ye,ynudag,anu,aedagyeynudaganu)
      call trace(aedagyeynudaganu,traedagyeynudaganu)

      call mat4pr(yedag,ae,anudag,ynu,yedagaeanudagynu)
      call trace(yedagaeanudagynu,tryedagaeanudagynu)  


     

      call mat4pr(addag,ad,yddag,yd,addagadyddagyd)
      call trace(addagadyddagyd,traddagadyddagyd)
      
      call mat4pr(addag,yd,yddag,ad,addagydyddagad)
      call trace(addagydyddagad,traddagydyddagad)
        
      call matmult(addag,yd,addagyd)
      call trace(addagyd,traddagyd)

      call mat4pr(aedag,ae,yedag,ye,aedagaeyedagye)
      call trace(aedagaeyedagye,traedagaeyedagye)

      call matmult(yddag,ad,yddagad)
      call trace(ydagad,trydagad)

      call mat3prod(yddag,mdsq,yd,yddagmdyd)
      call trace(yddagmdyd,tryddagmdyd)

      call matmult(aedag,ye,aedagye)
      call trace(aedagye,traedagye)
     
      call mat5pr(mqsq,yudag,yu,yddag,yd,mqyudagyuyddagyd)
      call trace(mqyudagyuyddagyd,trmqyudagyuyddagyd)
      
      call mat3prod(mlsq,ynudag,ynu,mlynudagynu)
      call trace(mlynudagynu, trmlynudagynu)

     
      call mat5pr(mlsq,ynudag,ynu,yedag,ye,mlynudagynuyedagye)
      call trace(mlynudagynuyedagye, trmlynudagynuyedagye)
      
      call mat5pr(mlsq,ynudag,ynu,ynudag,ynu,mlynudagynuynudagynu)
      call trace(mlynudagynuynudagynu,trmlynudagynuynudagynu)
     
      trmqyddagyd = 0.d0
      
      call matmult(mqsq,yddagyd,mqyddagyd)

!      call mat3prod(mqsq,yddag,yd,mqyddagyd)
      call trace(mqyddagyd,trmqyddagyd)
  
      call mat4pr(aedag,ye,yedag,ae,aedagyeyedagae)
      call trace(aedagyeyedagae,traedagyeyedagae)
 
      call matmult(audag,yu,audagyu)
      call trace(audagyu,traudagyu)
      
C      print*, "traudagyu=",traudagyu               !<------------------------

      call mat4pr(audag,yu,yddag,ad,audagyuyddagad)
      call trace(audagyuyddagad,traudagyuyddagad)
  
      call mat3prod(mlsq,yedag,ye,mlyedagye)
      call trace(mlyedagye,trmlyedagye)

      call mat5pr(mlsq,yedag,ye,yedag,ye,mlyedagyeyedagye)
      call trace(mlyedagyeyedagye,trmlyedagyeyedagye)

      call mat4pr(anudag,anu,yedag,ye,anudaganuyedagye)
      call trace(anudaganuyedagye,tranudaganuyedagye)
 
      call matmult(anudag,ynu,anudagynu)
      call trace(anudagynu,tranudagynu)

      call mat4pr(anudag,ynu,yedag,ae,anudagynuyedagae)
      call trace(anudagynuyedagae,tranudagynuyedagae)

      call mat4pr(audag,au,yddag,yd,audagauyddagyd)
      call trace(audagauyddagyd,traudagauyddagyd)

      call mat5pr(mQsq,yddag,yd,yddag,yd,mqyddagydyddagyd)
      call trace(mqyddagydyddagyd,trmqyddagydyddagyd)

      call mat3prod(mqsq,yudag,yu,mqyudagyu)
      call trace(mqyudagyu,trmqyudagyu)

      call mat5pr(mqsq,yudag,yu,yudag,yu,mqyudagyuyudagyu)
      call trace(mqyudagyuyudagyu,trmqyudagyuyudagyu)

      call mat5pr(yddag,mdsq,yd,yddag,yd,yddagmdydyddagyd)
      call trace(yddagmdydyddagyd,tryddagmdydyddagyd)

      call mat4pr(yddag,yd,yddag,yd,yddagydyddagyd) 
      call trace(yddagydyddagyd,tryddagydyddagyd)
      call mat4pr(yedag,ye,yedag,ye,yedagyeyedagye)
      call trace(yedagyeyedagye,tryedagyeyedagye)


      call mat5pr(yedag,mesq,ye,yedag,ye,yedagmeyeyedagye)
      call trace(yedagmeyeyedagye,tryedagmeyeyedagye)

      call mat4pr(yedag,ye,anudag,anu,yedagyeanudaganu)
      call trace(yedagyeanudaganu,tryedagyeanudaganu)

     
      call mat4pr(ynudag,anu,aedag,ye,ynudaganuaedagye)
      call trace(ynudaganuaedagye,trynudaganuaedagye)

      call mat3prod(ynudag,mnusq,ynu,ynudagmnuynu)
      call trace(ynudagmnuynu,trynudagmnuynu)

      call mat5pr(ynudag,mnusq,ynu,yedag,ye,ynudagmnuynuyedagye)
      call trace(ynudagmnuynuyedagye,trynudagmnuynuyedagye)

      call mat5pr(ynudag,mnusq,ynu,ynudag,ynu,ynudagmnuynuynudagynu)
      call trace(ynudagmnuynuynudagynu,trynudagmnuynuynudagynu)
 
      call mat4pr(ynudag,ynu,aedag,ae,ynudagynuaedagae)
      call trace(ynudagynuaedagae,trynudagynuaedagae)

      call mat5pr(ynudag,ynu,mlsq,yedag,ye,ynudagynumlyedagye)
      call trace(ynudagynumlyedagye,trynudagynumlyedagye)

      call mat5pr(ynudag,ynu,yedag,mesq,ye,ynudagynuyedagmeye)
      call trace(ynudagynuyedagmeye,trynudagynuyedagmeye)

      call mat4pr(ynudag,ynu,yedag,ye,ynudagynuyedagye)
      call trace(ynudagynuyedagye,trynudagynuyedagye)

      call mat4pr(ynudag,ynu,ynudag,ynu,ynudagynuynudagynu)
      call trace(ynudagynuynudagynu,trynudagynuynudagynu)

      call mat4pr(yudag,au,addag,yd,yudagauaddagyd)
      call trace(yudagauaddagyd,tryudagauaddagyd)

      call mat3prod(yudag,musq,yu,yudagmuyu) 
      call trace(yudagmuyu,tryudagmuyu)

      call mat5pr(yudag,musq,yu,yddag,yd,yudagmuyuyddagyd)
      call trace(yudagmuyuyddagyd,tryudagmuyuyddagyd)

      call mat5pr(yudag,musq,yu,yudag,yu,yudagmuyuyudagyu)
      call trace(yudagmuyuyudagyu,tryudagmuyuyudagyu)

      call mat4pr(yudag,yu,addag,ad,yudagyuaddagad)
      call trace(yudagyuaddagad,tryudagyuaddagad)

      call mat5pr(yudag,yu,mqsq,yddag,yd,yudagyumqyddagyd)
      call trace(yudagyumqyddagyd,tryudagyumqyddagyd)

      call mat5pr(yudag,yu,yddag,mdsq,yd,yudagyuyddagmdyd)
      call trace(yudagyuyddagmdyd,tryudagyuyddagmdyd)

      call mat4pr(yudag,yu,yddag,yd,yudagyuyddagyd)
      call trace(yudagyuyddagyd,tryudagyuyddagyd)
 
      call mat3prod(yddag,yd,mqsq,yddagydmq)
      call mat4pr(yudag,yu,yudag,au,yudagyuyudagau)
      call mat4pr(yudag,yu,yddag,ad,yudagyuyddagad)
      call mat4pr(yudag,au,yudag,yu,yudagauyudagyu)
      call mat4pr(yudag,au,yddag,yd,yudagauyddagyd)
      call matmult(ynu,yedag,ynuyedag)
      call mat4pr(ynudag,ynu,ynudag,anu,ynudagynuynudaganu)
      call mat4pr(ynudag,ynu,yedag,ae,ynudagynuyedagae)
      call mat4pr(ynudag,anu,ynudag,ynu,ynudaganuynudagynu) 
      call mat4pr(ynudag,anu,yedag,ye,ynudaganuyedagye)
      call matmult(ye,ynudag,yeynudag)
      call mat4pr(yedag,ye,ynudag,anu,yedagyeynudaganu)
      call mat4pr(yedag,ye,yedag,ae,yedagyeyedagae)
      call mat4pr(yedag,ae,ynudag,ynu,yedagaeynudagynu)
      call mat4pr(yedag,ae,yedag,ye,yedagaeyedagye)
      call mat4pr(yddag,yd,yudag,au,yddagydyudagau)
      call mat4pr(yddag,yd,yddag,ad,yddagydyddagad)
      call mat4pr(yddag,ad,yudag,yu,yddagadyudagyu)
      call mat4pr(yddag,ad,yddag,yd,yddagadyddagyd)
      call mat4pr(ynu,yedag,ae,anudag,ynuyedagaeanudag)
      call mat4pr(anu,aedag,ye,ynudag,anuaedagyeynudag)
      call mat4pr(ynu,ynudag,anu,anudag,ynuynudaganuanudag)
      call mat4pr(ynu,anudag,anu,ynudag,ynuanudaganuynudag)
      call mat4pr(anu,anudag,ynu,ynudag,anuanudagynuynudag)
 
      call trace(yudagyuyudagyu,tryudagyuyudagyu)

      call mat3prod(yudag,yu,mqsq,yudagyumq)
      call mat3prod(yu,yudag,musq,yuyudagmu)
      call mat3prod(yd,yddag,mdsq,ydyddagmd)
      call mat3prod(yedag,ye,mLsq,yedagyeml)

      call mat3prod(ye,yedag,mesq,yeyedagme)
      call mat3prod(ynu,ynudag,mNUsq,ynuynudagmnu)
      call mat4pr(yddag,ad,addag,yd,yddagadaddagyd)
      call mat3prod(musq,yudag,yu,muyudagyu)
      call mat3prod(mdsq,yd,yddag,mdydyddag)

      call mat3prod(ynudag,ynu,mlsq,ynudagynuml)
      call mat4pr(yedag,ae,aedag,ye,yedagaeaedagye)
      call mat3prod(mesq,ye,yedag,meyeyedag)

      call mat3prod(mnusq,ynu,ynudag,mnuynuynudag)
      call mat4pr(addag,ad,yddag,yd,addagadyddagyd)
      call mat4pr(yd,audag,au,yddag,ydaudagauyddag)
      call mat4pr(yedag,ye,aedag,ae,yedagyeaedagae)
      call mat4pr(ye,anudag,anu,yedag,yeanudaganuyedag)
      call mat5pr(ynu,yedag,ye,ynudag,mnusq,ynuyedagyeynudagmnu)
      call mat4pr(yd,audag,au,yddag,ydaudagauyddag)
      call mat4pr(ynudag,anu,anudag,ynu,ynudaganuanudagynu)
      call mat4pr(ae,ynudag,ynu,aedag,aeynudagynuaedag)
      call mat5pr(ynu,yedag,ye,mlsq,ynudag,ynuyedagyemlynudag)

      call mat4pr(yddag,yd,addag,ad,yddagydaddagad)
      call mat4pr(ad,yudag,yu,addag,adyudagyuaddag)
      call mat4pr(ynudag,ynu,anudag,anu,ynudagynuanudaganu)
      call mat4pr(ye,ynudag,anu,aedag,yeynudaganuaedag)
      call mat5pr(ynu,yedag,mesq,ye,ynudag,ynuyedagmeyeynudag)

      call mat4pr(yudag,au,audag,yu,yudagauaudagyu)
      call mat4pr(yd,yudag,au,addag,ydyudagauaddag)
      call mat5pr(yedag,ye,yedag,ye,mlsq,yedagyeyedagyeml)
      call mat4pr(ae,anudag,ynu,yedag,aeanudagynuyedag)
      call mat5pr(ynu,mlsq,yedag,ye,ynudag,ynumlyedagyeynudag)

      call mat4pr(yudag,yu,audag,au,yudagyuaudagau)
      call mat4pr(ad,audag,yu,yddag,adaudagyuyddag)
      call mat5pr(yedag,ye,yedag,mEsq,ye,yedagyeyedagmeye)
      call mat4pr(ye,aedag,ae,yedag,yeaedagaeyedag)
      call mat4pr(yu,addag,ad,yudag,yuaddagadyudag)
      call mat5pr(mnusq,ynu,yedag,ye,ynudag,mnuynuyedagyeynudag)

      call mat5pr(yddag,yd,yddag,yd,mqsq,yddagydyddagydmq)
      call mat4pr(yd,addag,ad,yddag,ydaddagadyddag)
      call mat5pr(yedag,ye,mlsq,yedag,ye,yedagyemlyedagye)
      call mat4pr(ae,yedag,ye,aedag,aeyedagyeaedag)
      call mat5pr(ynu,ynudag,ynu,ynudag,mnusq,ynuynudagynuynudagmnu)
      
      call mat5pr(yddag,yd,yddag,mdsq,yd,yddagydyddagmdyd)
      call mat4pr(ad,yddag,yd,addag,adyddagydaddag)
      call mat5pr(ynudag,ynu,ynudag,ynu,mqsq,ynudagynuynudagynumq)
      call mat4pr(ye,yedag,ae,aedag,yeyedagaeaedag)
      call mat5pr(ynu,ynudag,ynu,mlsq,ynudag,ynuynudagynumlynudag)
      call mat5pr(yddag,yd,mqsq,yddag,yd,yddagydmqyddagyd)
      call mat4pr(ad,addag,yd,yddag,adaddagydyddag)
      call mat5pr(ynudag,ynu,ynudag,ynu,mqsq,ynudagynuynudagynumq)
      call mat4pr(ae,aedag,ye,yedag,aeaedagyeyedag)
      call mat5pr(ynu,ynudag,mnusq,ynu,ynudag,ynuynudagmnuynuynudag)
      call mat5pr(yudag,yu,yudag,yu,mqsq,yudagyuyudagyumq)
      call mat4pr(ae,aedag,ye,yedag,aeaedagyeyedag)
      call mat5pr(yd,yudag,yu,yddag,mdsq,ydyudagyuyddagmd)
      
      call mat5pr(ynu,mlsq,ynudag,ynu,ynudag,ynumlynudagynuynudag)
      call mat5pr(ynudag,ynu,ynudag,mnusq,ynu,ynudagynuynudagmnuynu)
      call mat5pr(yd,yudag,yu,mqsq,yddag,ydyudagyumqyddag)
      call mat5pr(ye,ynudag,ynu,yedag,mEsq,yeynudagynuyedagme)
      
      call mat5pr(mnusq,ynu,ynudag,ynu,ynudag,mnuynuynudagynuynudag)
      call mat5pr(yudag,yu,yudag,musq,yu,yudagyuyudagmuyu)
      call mat5pr(yudag,yu,mqsq,yudag,yu,yudagyumqyudagyu)
      call mat5pr(yd,yudag,musq,yu,yddag,ydyudagmuyuyddag)
      call mat5pr(ynudag,ynu,mlsq,ynudag,ynu,ynudagynumlynudagynu)
      call mat5pr(ye,ynudag,ynu,mlsq,yedag,yeynudagynumlyedag)
      call mat5pr(yd,mqsq,yudag,yu,yddag,ydmqyudagyuyddag)
      call mat5pr(ye,ynudag,mnusq,ynu,yedag,yeynudagmnuynuyedag)
      
      call mat5pr(yu,yddag,yd,mqsq,yudag,yuyddagydmqyudag)
      call mat4pr(yd,yudag,yu,yddag,ydyudagyuyddag)
      call mat5pr(ye,mlsq,ynudag,ynu,yedag,yemlynudagynuyedag)
      call mat4pr(ye,ynudag,ynu,yedag,yeynudagynuyedag)
      call mat5pr(mdsq,yd,yudag,yu,yddag,mdydyudagyuyddag)
      call mat5pr(yd,yddag,yd,yddag,mdsq,ydyddagydyddagmd)
      call mat5pr(mesq,ye,ynudag,ynu,yedag,meyeynudagynuyedag)
      call mat5pr(yd,yddag,yd,mqsq,yddag,ydyddagydmqyddag)
      call mat5pr(ye,yedag,ye,yedag,mesq,yeyedagyeyedagme)
      call mat5pr(ye,yedag,ye,mlsq,yedag,yeyedagyemlyedag)
      call mat5pr(yd,yddag,mdsq,yd,yddag,ydyddagmdydyddag)
      call mat5pr(yd,mqsq,yddag,yd,yddag,ydmqyddagydyddag)
      call mat5pr(ye,yedag,mesq,ye,yedag,yeyedagmeyeyedag)
      call mat5pr(mdsq,yd,yddag,yd,yddag,mdydyddagydyddag)
      call mat5pr(ye,mlsq,yedag,ye,yedag,yemlyedagyeyedag)
      call mat5pr(mesq,ye,yedag,ye,yedag,meyeyedagyeyedag)
      call mat5pr(yu,yddag,yd,yudag,mUsq,yuyddagydyudagmu)
      call mat5pr(yu,yddag,mdsq,yd,yudag,yuyddagmdydyudag)
      call mat5pr(yu,mqsq,yddag,yd,yudag,yumqyddagydyudag)
      call mat5pr(musq,yu,yddag,yd,yudag,muyuyddagydyudag)
      call mat5pr(yu,yudag,yu,yudag,musq,yuyudagyuyudagmu)
      call mat5pr(yu,yudag,yu,mqsq,yudag,yuyudagyumqyudag)
      call mat5pr(yu,yudag,musq,yu,yudag,yuyudagmuyuyudag)
      call mat5pr(yu,mqsq,yudag,yu,yudag,yumqyudagyuyudag)
      call mat5pr(musq,yu,yudag,yu,yudag,muyuyudagyuyudag)
      call mat3prod(musq,yu,yudag,muyuyudag)
      call mat4pr(au,yddag,yd,audag,auyddagydaudag)
      call mat4pr(yu,yddag,ad,audag,yuyddagadaudag)
      call mat4pr(au,addag,yd,yudag,auaddagydyudag)
      call mat4pr(yu,audag,au,yudag,yuaudagauyudag)
      call mat4pr(au,yudag,yu,audag,auyudagyuaudag)
      call mat4pr(yu,yudag,au,audag,yuyudagauaudag)
      call mat4pr(au,audag,yu,yudag,auaudagyuyudag)                
      call mat4pr(anu,yedag,ye,anudag,anuyedagyeanudag)
!----------------      

      
      call mat4pr(ad,yddag,yd,yddag,adyddagydyddag)
      call trace(adyddagydyddag,tradyddagydyddag)
 
      call mat4pr(ad,yudag,yu,yddag,adyudagyuyddag)
      call trace(adyudagyuyddag,tradyudagyuyddag)

      call mat4pr(ae,yedag,ye,yedag,aeyedagyeyedag)
      call trace(aeyedagyeyedag,traeyedagyeyedag)
      
      call mat4pr(ae,ynudag,ynu,yedag,aeynudagynuyedag)
      call trace(aeynudagynuyedag,traeynudagynuyedag)

      call mat4pr(anu,yedag,ye,ynudag,anuyedagyeynudag)
      call trace(anuyedagyeynudag,tranuyedagyeynudag)

      call mat4pr(anu,ynudag,ynu,ynudag,anuynudagynuynudag)
      call trace(anuynudagynuynudag,tranuynudagynuynudag)

      call mat4pr(au,yddag,yd,yudag,auyddagydyudag)
      call trace(auyddagydyudag,trauyddagydyudag)

      call mat4pr(au,yudag,yu,yudag,auyudagyuyudag)
      call mat4pr(anu,ynudag,ynu,anudag,anuynudagynuanudag)
      call trace(auyudagyuyudag,trauyudagyuyudag)

 
      call trace(yuyddagydyudag,tryuyddagydyudag)
      call mat3prod(ye,ynudag,ynu,yeynudagynu)
      call mat5pr(yd,yudag,yu,yddag,yd,ydyudagyuyddagyd)
      call trace(ydyddagydyddag,trydyddagydyddag)
      call mat5pr(yd,yudag,yu,yudag,yu,ydyudagyuyudagyu)
      call mat5pr(yd,yudag,yu,yudag,yu,ydyudagyuyudagyu)
      call mat5pr(yd,yddag,yd,yddag,yd,ydyddagydyddagyd)
      call trace(ynuyedagyeynudag,trynuyedagyeynudag)
      call trace(yeyedagyeyedag,tryeyedagyeyedag)
      call mat4pr(ad,yudag,yu,yddag,adyudagyuyddag)
      call mat5pr(ynudag,ynu,ynudag,ynu,mlsq,ynudagynuynudagynuml)

      call mat4pr(yd,yddag,yu,yudag,ydyddagyuyudag)
      call trace(ydyddagyuyudag,trydyddagyuyudag)
      call mat4pr(yu,yudag,yd,yddag,yuyudagydyddag)
      call trace(yuyudagydyddag,tryuyudagydyddag)


C     -----------------------------------------------
C     Gauge Part of the Yukawa RGES (one loop terms)
C     ----------------------------------------------

      gut  = (8.d0/3.d0)*(alph3) + (1.5d0)*(alph2) + 
     $     (13.d0/30.d0)*(alph1)

      
      gdt  = (8.d0/3.d0)*(alph3) + (1.5d0)*(alph2) +
     $     (7.d0/30.d0)*(alph1)
      
      get  = (1.5d0)*(alph2) + (9.d0/10.d0)*(alph1) 
      
      gnut = (1.5d0)*(alph2) + (3.d0/10.d0)*(alph1) 
 
C     -------------------------------------------------
C     Gaugino Part of the A-parameters (one loop terms)
C     --------------------------------------------------

      gAut =   (16.d0/3.d0)*alph3*M3t + (3.d0)*alph2*M2t + 
     . (13.d0/15.d0)*alph1*M1t

      gAdt =   (16.d0/3.d0)*alph3*M3t + (3.d0)*alph2*M2t +
     . (7.d0/15.d0)*alph1*M1t  

      gAet =    (3.d0)*alph2*M2t + (9.d0/5.d0)*alph1*M1t

      gAnut =   3.d0*alph2*M2t + (3.d0/5.d0)*alph1*M1t 


C     ------------------------------------------------------
C     Gaugino part for the mass parameters (one loop terms)
C     ------------------------------------------------------

 
      S = trmqsq + trmdsq - (2.d0)*trmusq - trmlsq + trmesq +
     $     mh1sq - mh2sq

      
      gQsq = (16.d0/3.d0)*alph3*M3t*M3t + 
     . (3.d0)*alph2*M2t*M2t + 
     . (1.d0/15.d0)*alph1*M1t*M1t - (1.d0/10.d0)*alph1*S


      gUsq = (16.d0/3.d0)*alph3*M3t*M3t + ((16.d0/15.d0)*alph1*M1t*M1t) 
     $     + (2.d0/5.d0)*alph1*S

      gDsq = (16.d0/3.d0)*alph3*M3t*M3t + ((4.d0/15.d0)*alph1*M1t*M1t) 
     $     - (1.d0/5.d0)*alph1*S  

      gLsq = (3.d0)*alph2*M2t*M2t + (3.d0/5.d0)*alph1*M1t*M1t + 
     $     (3.d0/10.d0)*alph1*S  

      gH1sq = (3.d0)*alph2*M2t*M2t + (3.d0/5.d0)*alph1*M1t*M1t - 
     .     (3.d0/10.d0)*alph1*S 
      
      gH2sq = (3.d0)*alph2*M2t*M2t + (3.d0/5.d0)*alph1*M1t*M1t +
     .     (3.d0/10.d0)*alph1*S 

      gEsq = (12.d0/5.d0)*alph1*M1t*M1t -(3.d0/5.d0)*alph1*S


            
      Sd =(-3.d0*mh1sq*tryuyudag - trmqyudagyu + (4.d0*tryudagmuyu)+
     .     3.d0*mh2sq*trydyddag - trmqyddagyd - (2.d0*tryddagmdyd) +
     .     mh2sq*tryedagye + trmlyedagye - 2.d0*tryedagmeye - 
     .     mh1sq*trynudagynu + trmlynudagynu) + 
     .    (((1.5d0*alph2) + ((3.d0/10.d0)*alph1))*
     .    (mh1sq - mh2sq - trmlsq))+
     .    (((8.d0/3.d0)*alph3)+(1.5d0*alph2)+
     .    ((1.d0/30.d0)*alph1))*(trmqsq) -
     .    ((((16.d0/3.d0)*alph3) + 
     .    ((16.d0/15.d0)*alph1))*trmusq)+
     .    ((((8.d0/3.d0)*alph3) + 
     .     ((2.d0/15.d0)*alph1))*trmdsq) +
     .    ((6.d0/5.d0)*alph1*trmesq)   
     
      sig1 = (1.d0/5.d0)*alph1*((3.d0*(mh1sq+mh2sq)) + trmqsq + 
     .       (3.d0*trmlsq ) + (8.d0*trmusq) + (2.d0*trmdsq) + 
     .       (6.d0*trmesq))
     
      sig2 = alph2*(mh1sq + mh2sq + (3.d0*trmqsq) + trmlsq) 
 
      sig3 = alph3*((2.d0*trmqsq) + trmusq + trmdsq)

C-----------------------------------------------------------setting yb1 to zero
       loopset01: do i=1,3
       loopset02:   do j=1,3
       
       ydb1(i,j)=0.d0
       yub1(i,j)=0.d0
       yeb1(i,j)=0.d0
       ynub1(i,j)=0.d0
       ydb2(i,j)=0.d0
       yub2(i,j)=0.d0
       yeb2(i,j)=0.d0
       ynub2(i,j)=0.d0
       b1ad(i,j)=0.d0
       b1ae(i,j)=0.d0
       b1anu(i,j)=0.d0
       b1au(i,j)=0.d0
       adA21(i,j)=0.d0
       auA21(i,j)=0.d0
       aeA21(i,j)=0.d0
       anuA21(i,j)=0.d0
       adA22(i,j)=0.d0
       auA22(i,j)=0.d0
       aeA22(i,j)=0.d0
       anuA22(i,j)=0.d0
       

      enddo loopset02
      enddo loopset01
C-----------------------------------------------------------      
C identity matrix

      id(1,1)=1.d0          
      id(1,2)=0.d0
      id(1,3)=0.d0
      id(2,1)=0.d0
      id(2,2)=1.d0
      id(2,3)=0.d0
      id(3,1)=0.d0
      id(3,2)=0.d0
      id(3,3)=1.d0

!     
C--------------------------------------------------------------------
C     one loop- yukawas-matrix terms : y_b1 is the output matrix
C--------------------------------------------------------------------
      call rgeb1(yuyudag,ynuynudag,yudagyu,yddagyd,yu,yub1)
      call rgeb1(ydyddag,yeyedag,yedagye,ynudagynu,ye,yeb1)      
      call rgeb1(yuyudag,ynuynudag,ynudagynu,yedagye,ynu,ynub1)
      call rgeb1(ydyddag,yeyedag,yddagyd,yudagyu,yd,ydb1)
      

C--------------------------------------------------------------------
C     one loop- A parameters - matrix terms : b1a_ is the output matrix
C--------------------------------------------------------------------

      call rgeA1(ydyddag,yeyedag,yddagyd,yudagyu,adyddag,aeyedag,
     $     yddagad,yudagau,ad,yd,b1ad)
      call rgeA1(yuyudag,ynuynudag,yudagyu,yddagyd,auyudag,anuynudag,
     $     yudagau,yddagad,au,yu,b1au)
      call rgeA1(ydyddag,yeyedag,yedagye,ynudagynu,adyddag,aeyedag,
     $     yedagae,ynudaganu,ae,ye,b1ae)
      call rgeA1(yuyudag,ynuynudag,ynudagynu,yedagye,auyudag,anuynudag,
     $     ynudaganu,yedagae,anu,ynu,b1anu)

C--------------------------------------------------------------------------
C     two loops- yukawas - matrix terms : y_b2 is the output matrix
C---------------------------------------------------------------------------


      call rgeb2(ydyddagydyddag,yuyddagydyudag,yeyedagyeyedag,
     *     ynuyedagyeynudag,yudagyu,yuyudag,ynuynudag,yddagyd,ydyddag,
     *     yeyedag,yddagydyddagyd,yudagyuyudagyu,yudagyuyddagyd,yd,ydb2)
      
      call rgeb2(yuyudagyuyudag,yuyddagydyudag,ynuynudagynuynudag,
     *     ynuyedagyeynudag,yddagyd,ydyddag,yeyedag,yudagyu,yuyudag,
     *     ynuynudag,yudagyuyudagyu,yddagydyddagyd,yddagydyudagyu,
     *     yu,yub2)
      
      call rgeb2(ydyddagydyddag,yuyddagydyudag,yeyedagyeyedag,
     *     ynuyedagyeynudag,ynudagynu,yuyudag,ynuynudag,yedagye,yddagyd,
     *     yeyedag,yedagyeyedagye,ynudagynuynudagynu,ynudagynuyedagye,
     *     ye,yeb2)
      
      call rgeb2(yuyudagyuyudag,yuyddagydyudag,ynuynudagynuynudag,
     *     ynuyedagyeynudag,yedagye,ydyddag,yeyedag,ynudagynu,yuyudag,
     *     ynuynudag,ynudagynuynudagynu,yedagyeyedagye,yedagyeynudagynu,
     *     ynu,ynub2)

C------------------------------------------------------------------------------------------
C     two loops- A parameters - matrix terms : a_a21 is the product matrix, A_ contribution
C-------------------------------------------------------------------------------------------

      call rgeA21(ydyddagydyddag,yuyddagydyudag,yeyedagyeyedag,
     *     ynuyedagyeynudag,yudagyu,yuyudag,ynuynudag,yddagyd,ydyddag,
     *     yeyedag,yddagydyddagyd,yudagyuyudagyu,yudagyuyddagyd,ad,
     $     ada21)

      call rgeA21(yuyudagyuyudag,yuyddagydyudag,ynuynudagynuynudag,
     *     ynuyedagyeynudag,yddagyd,ydyddag,yeyedag,yudagyu,yuyudag,
     *     ynuynudag,yudagyuyudagyu,yddagydyddagyd,yddagydyudagyu,au,
     $     aua21)

      call rgeA21(ydyddagydyddag,yuyddagydyudag,yeyedagyeyedag,
     *     ynuyedagyeynudag,ynudagynu,yuyudag,ynuynudag,yedagye,ydyddag,
     *     yeyedag,yedagyeyedagye,ynudagynuynudagynu,ynudagynuyedagye,
     $     ae,aea21)

      call rgeA21(yuyudagyuyudag,yuyddagydyudag,ynuynudagynuynudag,
     *     ynuyedagyeynudag,yedagye,ydyddag,yeyedag,ynudagynu,yuyudag,
     *     ynuynudag,ynudagynuynudagynu,yedagyeyedagye,yedagyeynudagynu,
     *     anu,anua21)

C------------------------------------------------------------------------------------------
C     two loops- A parameters - matrix terms : a_a21 is the product matrix, Y_ contribution
C-------------------------------------------------------------------------------------------
        
        
        call rgeA22(adyddagydyddag,auyddagydyudag,adyudagyuyddag,
     *  aeyedagyeyedag,anuyedagyeynudag,aeynudagynuyedag,yudagyu,
     *  auyudag,anuynudag,yddagyd,adyddag,aeyedag,yudagau,yuyudag,
     *  ynuynudag,yddagad,ydyddag,yeyedag,yddagydyddagad,yddagadyddagyd,
     *  yudagauyudagyu,yudagyuyudagau,yudagauyddagyd,yudagyuyddagad,yd,
     *  ada22)

        call rgeA22(auyudagyuyudag,auyddagydyudag,adyudagyuyddag,
     *  anuynudagynuynudag,anuyedagyeynudag,aeynudagynuyedag,yudagyu,
     *  auyudag,anuynudag,yddagyd,adyddag,aeyedag,yudagau,yuyudag,
     *  ynuynudag,yddagad,ydyddag,yeyedag,yudagyuyudagau,yudagauyudagyu,
     *  yddagydyddagad,yddagadyddagyd,yddagydyudagau,yddagadyudagyu,yu,
     *  aua22)

        call rgeA22(adyddagydyddag,auyddagydyudag,adyudagyuyddag,
     *  aeyedagyeyedag,anuyedagyeynudag,aeynudagynuyedag,ynudagynu,
     *  auyudag,anuynudag,yedagye,adyddag,aeyedag,ynudaganu,yuyudag,
     *  ynuynudag,yddagad,ydyddag,yeyedag,yedagyeyedagae,yedagaeyedagye,
     *  ynudaganuynudagynu,ynudagynuynudaganu,ynudaganuyedagye,
     *  ynudagynuyedagae,ye,aea22)

        call rgeA22(auyudagyuyudag,auyddagydyudag,adyudagyuyddag,
     *  anuynudagynuynudag,anuyedagyeynudag,aeynudagynuyedag,ynudagynu,
     *  auyudag,anuynudag,yedagye,adyddag,aeyedag,ynudaganu,yuyudag,
     *  ynuynudag,yedagae,ydyddag,yeyedag,ynudagynuynudaganu,
     *  ynudaganuynudagynu,yedagyeyedagae,yedagaeyedagye,
     *  yedagyeynudaganu,yedagaeynudagynu,ynu,anua22)


C     --------------------------------------------------------------
C     BETA1 & BETA2 FOR YUKAWAS,A PARAMETERS & MASS PARAMETERS
C     --------------------------------------------------------------
        c=1
        
        rgei: do i=1,3
        rgej: do j=1,3
        
C     yu terms--------------------------------------------------------
        beta1yu(c) = + (gut*yu(i,j)) - (yub1(i,j)) !<--------------- one loop
        
        beta2yu(c) = yub2(i,j)  !<--------------- two loop begins
     $       -(8.d0*(alph3)+((2.d0/5.d0)*(alph1)))*tryuyudag*yu(i,j)-   
     .       ((3.d0*(alph2))+((1.d0/5.d0)*(alph1)))*yuyudagyu(i,j) 
     .       -((1.d0/5.d0)*(alph1)*yuyddagyd(i,j)) +
     .       ((8.d0/9.d0)*(alph3**2)*yu(i,j)) -
     .       (4.d0*(alph3)*(alph2)*yu(i,j))- 
     .       (68.d0/45.d0)*(alph3)*(alph1)*yu(i,j) - 
     .       (15.d0/4.d0)*(alph2**2)*yu(i,j)-
     .       (1.d0/2.d0)*(alph2)*(alph1)*yu(i,j)-
     .       (2743.d0/900.d0)*(alph1**2)*yu(i,j)
        
C     yd terms--------------------------------------------------------

      beta1yd(c) = (gdt*yd(i,j)) - (ydb1(i,j))                !<--------------- one loop 
      
      beta2yd(c) = ydb2(i,j)                                    !<--------------- two loops begin  
     $ -((8.d0*(alph3))-((1.d0/5.d0)*(alph1)))*trydyddag*yd(i,j) -  
     . ((3.d0*(alph2))+((2.d0/5.d0)*(alph1)))*ydyddagyd(i,j)  
     .  -((2.d0/5.d0)*(alph1)*ydyudagyu(i,j)) -                 
     $  (3.d0/5.d0)*(alph1)*tryeyedag*yd(i,j)+
     .  ((8.d0/9.d0)*(alph3**2)*yd(i,j)) -
     .  (4.d0*(alph3)*(alph2)*yd(i,j))- 
     .  (4.d0/9.d0)*(alph3)*(alph1)*yd(i,j) - 
     .  (15.d0/4.d0)*(alph2**2)*yd(i,j)-
     .  (0.5d0)*(alph2)*(alph1)*yd(i,j)-
     .  (287.d0/180.d0)*(alph1**2)*yd(i,j)

C ye terms--------------------------------------------------------------

      beta1ye(c) =  get*ye(i,j) - (yeb1(i,j))     !<--------------- one loop 

      beta2ye(c) = yeb2(i,j)                 !<--------------- two loops begin
     $  -((8.d0*(alph3))-((1.d0/5.d0)*(alph1)))*trydyddag*ye(i,j)-  
     .  (3.d0/5.d0)*(alph1)*tryeyedag*ye(i,j)-
     .  (3.d0)*(alph2)*yeyedagye(i,j) - 
     .  ((15.d0/4.d0)*(alph2**2)*ye(i,j)-
     .  (9.d0/10.d0)*((alph2*alph1))*ye(i,j)-
     .  (27.d0/4.d0)*(alph1**2)*ye(i,j))
      
C ynu terms---------------------------------------------------------------

      beta1ynu(c) = gnut*ynu(i,j)  - ynub1(i,j)   !<--------------- one loop 

      beta2ynu(c) =  ynub2(i,j)                  !<--------------- two loops begin
     $ -((8.d0*alph3)+((2.d0/5.d0)*alph1))*tryuyudag*ynu(i,j)-   
     .  ((3.d0*alph2)+((3.d0/5.d0)*(alph1)))*ynuynudagynu(i,j) -
     .  (3.d0/5.d0)*(alph1)*ynuyedagye(i,j) -
     .  (15.d0/4.d0)*(alph2**2)*ynu(i,j) -
     .  (9.d0/10.d0)*((alph2*alph1))*ynu(i,j)-
     .  (207.d0/100.d0)*(alph1**2)*ynu(i,j)


      
C au terms--------------------------------------------------------- 

      beta1au(c) = -(b1au(i,j))                    !<--------------- one loop 
     $  + gut*au(i,j) -  (gAut*yu(i,j))


      beta2au(c) =  aua21(i,j) + aua22(i,j)       !<--------------- two loops begin
     .  -((8.d0*alph3)+((2.d0/5.d0)*alph1))*tryuyudag*au(i,j)
     .  -(6.d0*alph2*auyudagyu(i,j)) -
     .  ((1.d0/5.d0)*alph1*auyddagyd(i,j))+ 
     .  ((((((8.d0/9.d0)*alph3*alph3) -
     .  (4.d0*alph3*alph2) -
     .  ((68.d0/45.d0)*alph3*alph1) -
     .  ((15.d0/4.d0)*(alph2*alph2)) -
     .  ((1.d0/2.d0)*alph2*alph1) -
     .  ((2743.d0/900.d0)*(alph1**2)))))*au(i,j)) + 
     .  ((-16.d0*alph3)-((4.d0/5.d0)*alph1))*trauyudag*yu(i,j)-
     .  (((3.d0*alph2)+((3.d0/5.d0)*alph1))*yuyudagau(i,j)) -
     .  ((2.d0/5.d0)*alph1*yuyddagad(i,j)) +
     .  ((16.d0*alph3*M3t)+((4.d0/5.d0)*alph1*M1t))*tryuyudag*yu(i,j)+
     .  ((6.d0*alph2*M2t)+((2.d0/5.d0)*alph1*M1t))*yuyudagyu(i,j) +
     .  (2.d0/5.d0)*alph1*M1t*yuyddagyd(i,j) +
     .  ((((-32.d0/9.d0)*alph3*alph3*M3t) +
     .  (8.d0*alph3*alph2*(M3t+M2t)) +
     .  ((136.d0/45.d0)*alph3*alph1*(M3t+M1t))+
     .  ((15.d0)*(alph2*alph2)*M2t) +
     .  (1.d0*alph2*alph1*(M2t+M1t))+
     .  ((2743.d0/225.d0)*(alph1*alph1)*M1t))*yu(i,j))

c$$$     .  (0.5d0*alph2*alph1*(M2t+M1t))+
c$$$     .  ((2743.d0/900.d0)*(alph1**2)*M1t))*yu(i,j))
    
C ad terms-------------------------------------------------------------- 

      beta1ad(c) =  - b1ad(i,j)                     !<--------------- one loop
     $     + gdt*ad(i,j) - (gAdt*yd(i,j))       

      beta2ad(c) =  ada21(i,j) + ada22(i,j) !<--------------- two loops begin
     . -(((8.d0*alph3)-((1.d0/5.d0)*alph1))*trydyddag*ad(i,j))
     .  -((2.d0/5.d0)*alph1*adyudagyu(i,j)) -
     .  ((3.d0/5.d0)*alph1*tryeyedag*ad(i,j)) +
     .  (((-6.d0*alph2)-(3.d0/5.d0)*(alph1))*adyddagyd(i,j)) +
     .  (((((8.d0/9.d0)*alph3**2) -
     .  (4.d0*alph3*alph2) -
     .  ((4.d0/9.d0)*alph3*alph1) -
     .  ((15.d0/4.d0)*(alph2**2)) -
     .  ((1.d0/2.d0)*alph2*alph1) -
     .  ((287.d0/180.d0)*(alph1**2)))*ad(i,j))) - 
     .  ((16.d0*alph3)-((2.d0/5.d0)*alph1))*tradyddag*yd(i,j)
     .  -(((3.d0*alph2)+(3.d0/5.d0)*(alph1))*ydyddagad(i,j)) -
     .  ((6.d0/5.d0)*alph1*traeyedag*yd(i,j)) -
     .  ((4.d0/5.d0)*alph1*ydyudagau(i,j)) +
     .  (((16.d0*alph3*M3t)+((2.d0/5.d0)*alph1*M1t))*trydyddag*yd(i,j))+
     .  ((6.d0/5.d0)*alph1*M1t*tryeyedag*yd(i,j))+
     .  (((6.d0*alph2*M2t)+((4.d0/5.d0)*alph1*M1t))*ydyddagyd(i,j)) +
     .  ((4.d0/5.d0)*alph1*M1t*ydyudagyu(i,j)) -
     .  ((((((32.d0/9.d0)*alph3**2*M3t) -
     .  (8.d0*alph3*alph2*(M3t+M2t)) -
     .  ((8.d0/9.d0)*alph3*alph1*(M3t+M1t))-
     .  ((15.d0)*(alph2**2)*M2t) -
     .  ((1.d0/1.d0)*alph2*alph1*(M2t+M1t))-
     .  ((287.d0/45.d0)*(alph1**2)*M1t))*yd(i,j))))

C ae terms-----------------------------------------------------------

      beta1ae(c) =  - b1ae(i,j)           !<--------------- one loop
     $  + get*ae(i,j) - (gAet*ye(i,j))    

      beta2ae(c) =  aea21(i,j) + aea22(i,j) !<--------------- two loops begin
     $ -(((8.d0*alph3)-((1.d0/5.d0)*alph1))*trydyddag*ae(i,j))-
     .  (((6.d0*alph2)-(3.d0/5.d0)*(alph1))*aeyedagye(i,j)) -
     .  ((3.d0/5.d0)*alph1*tryeyedag*ae(i,j)) -
     .  ((((15.d0/4.d0)*(alph2**2)) +
     .  ((9.d0/10.d0)*alph2*alph1) +
     .  ((27.d0/4.d0)*(alph1**2)))*ae(i,j)) -
     .  (((16.d0*alph3)-((2.d0/5.d0)*alph1))*tradyddag*ye(i,j))
     .  -(((3.d0*alph2)+((3.d0/5.d0)*alph1))*yeyedagae(i,j)) -
     .  ((6.d0/5.d0)*alph1*traeyedag*ye(i,j)) +
     .  (((16.d0*alph3*M3t)-((2.d0/5.d0)*alph1*M1t))*
     $  trydyddag*ye(i,j))+
     .  (((6.d0*alph2*M2t))*yeyedagye(i,j)) +
     .  ((6.d0/5.d0)*alph1*M1t*tryeyedag*ye(i,j)) +
     .  ((((((15.d0)*(alph2**2)*M2t) +
     .  ((9.d0/5.d0)*alph2*alph1*(M2t+M1t)) +
     .  (27.d0*(alph1**2)*M1t))*ye(i,j))))

C anu terms------------------------------------------------------------ 

      beta1anu(c) =  -b1anu(i,j)                  !<--------------- one loop  
     .  + (gnut*anu(i,j)) - (gAnut*ynu(i,j))  

      beta2anu(c) = + anua21(i,j) + anua22(i,j)         !<--------------- two loops begin
     .  -(((8.d0*(alph3))+
     .  ((2.d0/5.d0)*(alph1)))*tryuyudag*anu(i,j))-
     .  (6.d0*(alph2)*anuynudagynu(i,j)) -
     .  ((6.d0/5.d0)*(alph1)*anuynudagynu(i,j)) -
     .  ((3.d0/5.d0)*alph1*anuyedagye(i,j)) -
     .  ((((((15.d0/4.d0)*(alph2**2)) +
     .  ((9.d0/10.d0)*alph2*alph1) +  
     .  ((207.d0/100.d0)*(alph1**2)))*anu(i,j)))) - 
     .  (((16.d0*(alph3))+
     .  ((4.d0/5.d0)*(alph1)))*trauyudag*ynu(i,j))-
     .  (((3.d0*(alph2))+(3.d0/5.d0)*(alph1))*ynuynudaganu(i,j)) -
     .  ((6.d0/5.d0)*(alph1)*ynuyedagae(i,j))+ 
     .  (((((16.d0*(alph3)*M3t)+
     .  ((4.d0/5.d0)*(alph1)*M1t))*tryuyudag*ynu(i,j))+
     .  (((6.d0*(alph2)*M2t)+
     .  ((6.d0/5.d0)*(alph1)*M1t))*ynuynudagynu(i,j))+
     .  ((6.d0/5.d0)*(alph1)*M1t*ynuyedagye(i,j))))+ 
     .  ((((15.d0)*(alph2**2)*M2t*ynu(i,j)) +
     .  ((9.d0/5.d0)*alph2*alph1*(M2t+M1t)*ynu(i,j))+
     .  ((207.d0/25.d0)*(alph1**2)*M1t*ynu(i,j))))


     
C 73 -78  mQ

        beta1mq(c) =  -(0.5d0*mqyudagyu(i,j))
     . - (mh1sq*yudagyu(i,j))-
     .  (0.5d0*mqyddagyd(i,j))-
     .  (mh2sq*yddagyd(i,j))-
     .  (yudagmuyu(i,j))-
     .  (yddagmdyd(i,j))-
     .  (0.5d0*yudagyumq(i,j))-
     .  (0.5d0*yddagydmq(i,j))-
     .  (audagau(i,j))-
     .  (addagad(i,j)) 
     .  +(gQsq)*id(i,j)

 
       beta2mq(c) = mqyudagyuyudagyu(i,j)+
     $       (4.d0*mh1sq*yudagyuyudagyu(i,j))+
     . (2.d0*yudagmuyuyudagyu(i,j))+
     . (2.d0*yudagyumqyudagyu(i,j))+
     . (2.d0*yudagyuyudagmuyu(i,j))+
     .  yudagyuyudagyumq(i,j)+
     .  mqyddagydyddagyd(i,j)+
     . (4.d0*mh2sq*yddagydyddagyd(i,j))+
     . (2.d0*yddagmdydyddagyd(i,j))+
     . (2.d0*yddagydmqyddagyd(i,j))+
     . (2.d0*yddagydyddagmdyd(i,j))+
     . yddagydyddagydmq(i,j)+
     . (((1.5d0)*tryudagyu +(0.5d0)*trynudagynu)*(mqyudagyu(i,j)+
     . (4.d0*mh1sq*yudagyu(i,j))+(2.d0*yudagmuyu(i,j))+yudagyumq(i,j)))+
     . (((1.5d0)*tryddagyd +(0.5d0)*tryedagye)*(mqyddagyd(i,j)+
     . (4.d0*mh2sq*yddagyd(i,j))+(2.d0*yddagmdyd(i,j))+yddagydmq(i,j)))+ !
     . (yudagyu(i,j)*(((3.d0)*trmqyudagyu)+((3.d0)*tryudagmuyu)+
     . trmlynudagynu+trynudagmnuynu))+ (yddagyd(i,j)*(((3.d0)*
     . trmqyddagyd)+((3.d0)*tryddagmdyd)+ trmlyedagye + tryedagmeye))+!
     . ((2.d0)*(yudagyuaudagau(i,j)+audagauyudagyu(i,j)+
     . yudagauaudagyu(i,j)+audagyuyudagau(i,j)+yddagydaddagad(i,j)+
     . addagadyddagyd(i,j)+yddagadaddagyd(i,j)+addagydyddagad(i,j)))+ !
     . audagau(i,j)*((3.d0*tryudagyu)+trynudagynu)+
     . yudagyu(i,j)*((3.d0*traudagau)+tranudaganu)+
     . audagyu(i,j)*((3.d0*tryudagau)+trynudaganu)+
     . yudagau(i,j)*((3.d0*traudagyu)+tranudagynu)+
     . addagad(i,j)*((3.d0*tryddagyd)+tryedagye)+
     . yddagyd(i,j)*((3.d0*traddagad)+traedagae)+
     . addagyd(i,j)*((3.d0*tryddagad)+tryedagae)+
     . yddagad(i,j)*((3.d0*traddagyd)+traedagye)-
     . (((1.d0/5.d0)*alph1)*((2.d0*mqyudagyu(i,j))+
     . (4.d0*mh1sq*yudagyu(i,j))+
     . (4.d0*yudagmuyu(i,j))+
     . (2.d0*yudagyumq(i,j))+
     . (4.d0*audagau(i,j))-
     . (4.d0*M1t*audagyu(i,j))-
     . (4.d0*M1t*yudagau(i,j))+
     . (8.d0*(M1t**2)*yudagyu(i,j))+
     . (mqyddagyd(i,j)+(2.d0)*mh2sq*yddagyd(i,j))+
     . (2.d0*yddagmdyd(i,j))+
     . (yddagydmq(i,j))+
     . (2.d0*addagad(i,j))-
     . (2.d0*M1t*addagyd(i,j))-
     . (2.d0*M1t*yddagad(i,j))+
     . (4.d0*(M1t**2)*yddagyd(i,j))))
     . -((1.d0/5.d0)*alph1*Sd)*id(i,j)+
     . ((64.d0/3.d0)*(alph3*M3t)**2)*id(i,j)-
     . (16.d0*((alph3*alph2))*((M3t**2)+(M2t**2)+(M2t*M3t)))*id(i,j)-
     . ((16.d0/45.d0)*((alph3*alph1))*
     . ((M3t*M3t)+(M1t*M1t)+(M1t*M3t)))*id(i,j)-
     . ((33.d0/2.d0)*(alph2*M2t)**2)*id(i,j)-
     . ((199.d0/150.d0)*(alph1*M1t)**2)*id(i,j)-
     . ((1.d0/5.d0)*((alph2*alph1))*
     . ((M2t*M2t)+(M1t*M1t)+(M1t*M2t)))*id(i,j)-
     . ((8.d0/3.d0)*alph3*sig3)*id(i,j)-(1.5d0*alph2*sig2)*id(i,j)-
     . ((1.d0/30.d0)*alph1*sig1)*id(i,j)

C----------------------------------------------------------------------   
             
C  mU 79-84

       beta1mu(c) = -(muyuyudag(i,j))-
     .  (2.d0*mh1sq*yuyudag(i,j))-
     .  (2.d0*yumqyudag(i,j))-
     .  (yuyudagmu(i,j))-
     .  (2.d0*auaudag(i,j))
     .  +(gUsq)*id(i,j)

       beta2mu(c) = (muyuyudagyuyudag(i,j))+ 
     .  (4.d0*mh1sq*yuyudagyuyudag(i,j)) +
     . (2.d0*yumqyudagyuyudag(i,j)) +
     . (2.d0*yuyudagmuyuyudag(i,j)) +
     . (2.d0*yuyudagyumqyudag(i,j)) +
     .  yuyudagyuyudagmu(i,j) +
     .  muyuyddagydyudag(i,j) +
     . (2.d0*(mh2sq+mh1sq)*yuyddagydyudag(i,j)) +
     . (2.d0*yumqyddagydyudag(i,j)) +
     . (2.d0*yuyddagmdydyudag(i,j))+
     . (2.d0*yuyddagydmqyudag(i,j))+
     . (yuyddagydyudagmu(i,j))+
     . (((3.d0)*tryudagyu + trynudagynu)*(muyuyudag(i,j)+
     . (4.d0*mh1sq*yuyudag(i,j))+(2.d0*yumqyudag(i,j))+yuyudagmu(i,j)))+
     . (4.d0*yuyudag(i,j)*(((3.d0)*trmqyudagyu)+((3.d0)*tryudagmuyu)+
     . (0.5d0*trmlynudagynu)+(0.5d0*trynudagmnuynu)))+
     . ((2.d0)*(auaudagyuyudag(i,j)+yuyudagauaudag(i,j)+
     .  auyudagyuaudag(i,j)+yuaudagauyudag(i,j)+auaddagydyudag(i,j)+
     .  yuyddagadaudag(i,j)+auyddagydaudag(i,j)+yuaddagadyudag(i,j)))+
     . ((2.d0*auaudag(i,j))*((3.d0*tryudagyu)+trynudagynu))+
     . ((2.d0*yuyudag(i,j))*((3.d0*traudagau)+tranudaganu))+
     . ((2.d0*auyudag(i,j))*((3.d0*traudagyu)+tranudagynu))+
     . ((2.d0*yuaudag(i,j))*((3.d0*tryudagau)+trynudaganu))-
     . (((-(1.d0/5.d0)*alph1)+(3.d0*alph2))*((muyuyudag(i,j))+
     . (2.d0*mh1sq*yuyudag(i,j))+
     . (2.d0*yumqyudag(i,j))+
     . (yuyudagmu(i,j))+
     . (2.d0*auaudag(i,j))))-
     . (6.d0*alph2*(-(M2t*auyudag(i,j))-(M2t*yuaudag(i,j))+
     . (2.d0*M2t*M2t*yuyudag(i,j)))) +
     . ((2.d0/5.d0)*alph1*((2.d0*M1t*M1t*yuyudag(i,j))-
     . (M1t*auyudag(i,j))-(M1t*yuaudag(i,j)))) +
     . ((4.d0/5.d0)*alph1*Sd)*id(i,j) +
     . ((64.d0/3.d0)*(alph3*M3t)**2)*id(i,j) -
     . ((512.d0/90.d0)*((alph3*alph1))*
     . ((M3t*M3t)+(M1t*M1t)+(M1t*M3t)))*id(i,j)-
     . ((3424.d0/150.d0)*(alph1*M1t)**2)*id(i,j)-
     . ((8.d0/3.d0)*alph3*sig3)*id(i,j)-
     . ((8.d0/15.d0)*alph1*sig1)*id(i,j)

C--------------------------------------------------------------------------
C 85-90 mD

      beta1md(c) = -(mdydyddag(i,j))
     . - (2.d0*mh2sq*ydyddag(i,j))-
     .  (2.d0*ydmqyddag(i,j))-
     .  (ydyddagmd(i,j))-
     .  (2.d0*adaddag(i,j))
     . + (gDsq)*id(i,j) 
     
      beta2md(c)=(mdydyddagydyddag(i,j))+
     $     (4.d0*mh2sq*ydyddagydyddag(i,j))+
     . (2.d0*ydmqyddagydyddag(i,j))+
     . (2.d0*ydyddagmdydyddag(i,j))+
     . (2.d0*ydyddagydmqyddag(i,j))+
     .  ydyddagydyddagmd(i,j)+
     .  mdydyudagyuyddag(i,j)+
     . (2.d0*(mh2sq+mh1sq)*ydyudagyuyddag(i,j))+
     . (2.d0*ydmqyudagyuyddag(i,j))+
     . (2.d0*ydyudagmuyuyddag(i,j))+
     . (2.d0*ydyudagyumqyddag(i,j))+
     . (ydyudagyuyddagmd(i,j))+
     . (((3.d0)*tryddagyd +tryedagye)*(mdydyddag(i,j)+
     . (4.d0*mh2sq*ydyddag(i,j))+(2.d0*ydmqyddag(i,j))+ydyddagmd(i,j)))+
     . (2.d0*ydyddag(i,j)*(((3.d0)*trmqyddagyd)+((3.d0)*tryddagmdyd)+
     . (trmlyedagye)+(tryedagmeye)))+
     . ((2.d0)*(adaddagydyddag(i,j)+ydyudagauaddag(i,j)+
     .  adyddagydaddag(i,j)+ydaddagadyddag(i,j)+adaudagyuyddag(i,j)+
     .  ydyudagauaddag(i,j)+adyudagyuaddag(i,j)+ydaudagauyddag(i,j))) +
     . (2.d0*adaddag(i,j)*((3.d0*tryddagyd)+tryedagye))+
     . (2.d0*ydyddag(i,j)*((3.d0*traddagad)+traedagae))+
     . (2.d0*adyddag(i,j)*((3.d0*traddagyd)+traedagye))+
     . (2.d0*ydaddag(i,j)*((3.d0*tryddagad)+tryedagae))-
     . ((((1.d0/5.d0)*alph1)+(3.d0*alph2))*((mdydyddag(i,j))+
     . (2.d0*mh2sq*ydyddag(i,j))+
     . (2.d0*ydmqyddag(i,j))+
     . (ydyddagmd(i,j))+
     . (2.d0*adaddag(i,j)))) -
     . (6.d0*alph2*((-M2t*adyddag(i,j))-(M2t*ydaddag(i,j))+
     . (2.d0*M2t*M2t*ydyddag(i,j)))) -
     . ((2.d0/5.d0)*alph1*((2.d0*M1t*M1t*ydyddag(i,j))-
     . (M1t*adyddag(i,j))-(M1t*ydaddag(i,j)))) -
     . ((2.d0/5.d0)*alph1*Sd)*id(i,j) +
     . ((64.d0/3.d0)*((alph3*M3t)**2))*id(i,j) -
     . ((64.d0/45.d0)*((alph3*alph1))*
     . ((M3t*M3t)+(M1t*M1t)+(M1t*M3t)))*id(i,j) -
     . ((808.d0/150.d0)*(alph1*M1t)**2)*id(i,j)-
     . ((8.d0/3.d0)*alph3*sig3)*id(i,j)-
     . ((2.d0/15.d0)*alph1*sig1)*id(i,j)

C---------------------------------------------------------------------------

C  91-96 mL

      beta1ml(c) = (-0.5d0*mlynudagynu(i,j))- 
     .  (mh1sq*ynudagynu(i,j)) -
     .  (0.5d0*mlyedagye(i,j)) -
     .  (mh2sq*yedagye(i,j)) -
     .  (ynudagmnuynu(i,j)) -
     .  (yedagmeye(i,j)) -
     .  (0.5d0*ynudagynuml(i,j)) -
     .  (0.5d0*yedagyeml(i,j)) -
     .  (anudaganu(i,j)) -
     .  (aedagae(i,j)) +
     .  (gLsq)*id(i,j)

      
      beta2ml(c)= mlynudagynuynudagynu(i,j)+ 
     . (4.d0*mh1sq*ynudagynuynudagynu(i,j)) +
     . (2.d0*ynudagmnuynuynudagynu(i,j)) +
     . (2.d0*ynudagynumlynudagynu(i,j)) +
     . (2.d0*ynudagynuynudagmnuynu(i,j)) +
     .  ynudagynuynudagynuml(i,j) +
     .  mlyedagyeyedagye(i,j) +
     . (4.d0*mh2sq*yedagyeyedagye(i,j)) +
     . (2.d0*yedagmeyeyedagye(i,j)) +
     . (2.d0*yedagyemlyedagye(i,j)) + 
     . (2.d0*yedagyeyedagmeye(i,j)) +
     .  yedagyeyedagyeml(i,j) +
     . ((((1.5d0)*tryudagyu) +((0.5d0)*trynudagynu))*((mlynudagynu(i,j))
     . +(4.d0*mh1sq*ynudagynu(i,j))+(2.d0*ynudagmnuynu(i,j))+
     . (ynudagynuml(i,j)))) +
     . (((1.5d0)*tryddagyd +(0.5d0)*tryedagye)*((mlyedagye(i,j)) +
     . (4.d0*mh2sq*yedagye(i,j))+(2.d0*yedagmeye(i,j))+yedagyeml(i,j)))+ 
     . ((ynudagynu(i,j)*(((3.d0)*trmqyudagyu)+((3.d0)*tryudagmuyu) +
     . trmlynudagynu+trynudagmnuynu))) + 
     . (yedagye(i,j)*(((3.d0)*trmqyddagyd) +
     . ((3.d0)*tryddagmdyd)+ trmlyedagye + tryedagmeye)) +
     . ((2.d0)*(ynudagynuanudaganu(i,j)+anudaganuynudagynu(i,j) +
     . ynudaganuanudagynu(i,j)+anudagynuynudaganu(i,j) + 
     . yedagyeaedagae(i,j) +
     . aedagaeyedagye(i,j)+yedagaeaedagye(i,j)+aedagyeyedagae(i,j))) +
     . anudaganu(i,j)*((3.d0*tryudagyu)+trynudagynu) +
     . ynudagynu(i,j)*((3.d0*traudagau)+tranudaganu) +
     . anudagynu(i,j)*((3.d0*tryudagau)+trynudaganu) +
     . ynudaganu(i,j)*((3.d0*traudagyu)+tranudagynu) +
     . aedagae(i,j)*((3.d0*tryddagyd)+tryedagye) +
     . yedagye(i,j)*((3.d0*traddagad)+traedagae) +
     . aedagye(i,j)*((3.d0*tryddagad)+tryedagae) +
     . yedagae(i,j)*((3.d0*traddagyd)+traedagye) -
     . ((3.d0/5.d0)*alph1*((mlyedagye(i,j)) +
     . (2.d0*mh2sq*yedagye(i,j)) +
     . (2.d0*yedagmeye(i,j)) +
     . (yedagyeml(i,j)) +
     . (2.d0*aedagae(i,j)) -
     . (2.d0*M1t*yedagae(i,j)) -
     . (2.d0*M1t*aedagye(i,j)) +
     . (4.d0*(M1t**2)*yedagye(i,j)))) +
     . ((3.d0/5.d0)*alph1*Sd)*id(i,j) -
     . ((33.d0/2.d0)*(alph2*M2t)**2)*id(i,j) -
     . ((621.d0/50.d0)*(alph1*M1t)**2)*id(i,j)-
     . ((9.d0/5.d0)*((alph2*alph1))*
     $     ((M2t*M2t)+(M1t*M1t)+(M1t*M2t)))*id(i,j)-
     . (1.5d0*alph2*sig2)*id(i,j)-
     . ((0.3d0)*alph1*sig1)*id(i,j)

C-----------------------------------------------------------------------------------
      
C  97 -102 mE

      beta1me(c) = (-meyeyedag(i,j))- 
     .  (2.d0*mh2sq*yeyedag(i,j)) -
     .  (2.d0*yemlyedag(i,j)) -
     .  (yeyedagme(i,j)) -
     .  (2.d0*aeaedag(i,j)) +
     .  (gEsq)*id(i,j) 

      beta2me(c) = meyeyedagyeyedag(i,j)+
     $     (4.d0*mh2sq*yeyedagyeyedag(i,j))+
     . (2.d0*yemlyedagyeyedag(i,j))+
     . (2.d0*yeyedagmeyeyedag(i,j))+
     . (2.d0*yeyedagyemlyedag(i,j))+
     .  yeyedagyeyedagme(i,j)+
     .  meyeynudagynuyedag(i,j)+
     . (2.d0*(mh2sq+mh1sq)*yeynudagynuyedag(i,j))+
     . (2.d0*yemlynudagynuyedag(i,j))+
     . (2.d0*yeynudagmnuynuyedag(i,j))+
     . (2.d0*yeynudagynumlyedag(i,j))+
     . (yeynudagynuyedagme(i,j))+
     . (((3.d0)*tryddagyd +tryedagye)*(meyeyedag(i,j)+
     . (4.d0*mh2sq*yeyedag(i,j))+(2.d0*yemlyedag(i,j))+yeyedagme(i,j)))+
     . (4.d0*yeyedag(i,j)*(((3.d0)*trmqyddagyd)+((3.d0)*tryddagmdyd)+
     . (0.5d0*trmlyedagye)+(0.5d0*tryedagmeye))) +
     . ((2.d0)*(aeaedagyeyedag(i,j)+yeyedagaeaedag(i,j)+
     . aeyedagyeaedag(i,j)+yeaedagaeyedag(i,j)+aeanudagynuyedag(i,j)+
     . yeynudaganuaedag(i,j)+aeynudagynuaedag(i,j) +
     . yeanudaganuyedag(i,j))) +
     . (2.d0*aeaedag(i,j)*((3.d0*tryddagyd)+tryedagye)) +
     . (2.d0*yeyedag(i,j)*((3.d0*traddagad)+traedagae)) +
     . (2.d0*aeyedag(i,j)*((3.d0*traedagye)+traedagye)) +
     . (2.d0*yeaedag(i,j)*((3.d0*tryedagad)+tryedagae)) -
     . ((((3.d0/5.d0)*alph1)+(3.d0*alph2))*((meyeyedag(i,j)) +
     . (2.d0*mh2sq*ydyddag(i,j)) +
     . (2.d0*yemlyedag(i,j)) +
     . (yeyedagme(i,j)) +
     . (2.d0*aeaedag(i,j)))) +
     . (6.d0*alph2*((M2t*aeyedag(i,j))+(M2t*yeaedag(i,j)) -
     . (2.d0*M2t*M2t*yeyedag(i,j)))) +
     . (((6.d0/5.d0)*alph1)*((2.d0*M1t*M1t*yeyedag(i,j)) -
     . (M1t*aeyedag(i,j))-(M1t*yeaedag(i,j)))) +
     . ((((6.d0/5.d0)*alph1*sd)- 
     . ((2808.d0/50.d0)*(alph1*M1t)**2) -
     . ((6.d0/5.d0)*alph1*sig1)))*id(i,j) 
 
C------------------------------------------------------------------------------    
C 103-108 mnu
    

       beta1mnu(c) = ((-mnuynuynudag(i,j))- 
     .  (2.d0*mh1sq*ynuynudag(i,j)) -
     .  (ynuynudagmnu(i,j)) -
     .  (2.d0*anuanudag(i,j)) -
     .  (2.d0*ynumlynudag(i,j)))


       beta2mnu(c) = mnuynuynudagynuynudag(i,j)+
     . (4.d0*mh1sq*ynuynudagynuynudag(i,j)) +
     . (2.d0*ynumlynudagynuynudag(i,j)) +
     . (2.d0*ynuynudagmnuynuynudag(i,j)) +
     . (2.d0*ynuynudagynumlynudag(i,j)) +
     .  ynuynudagynuynudagmnu(i,j) +
     .  mnuynuyedagyeynudag(i,j) +
     . (2.d0*(mh2sq+mh1sq)*ynuyedagyeynudag(i,j)) +
     . (2.d0*ynumlyedagyeynudag(i,j)) +
     . (2.d0*ynuyedagmeyeynudag(i,j)) +
     . (2.d0*ynuyedagyemlynudag(i,j)) +
     . (ynuyedagyeynudagmnu(i,j)) +
     . (((3.d0)*tryudagyu +trynudagynu)*(mnuynuynudag(i,j) + 
     . (4.d0*mh1sq*ynuynudag(i,j))+
     $      (2.d0*ynumlynudag(i,j))+ynuynudagmnu(i,j))) 
     . + (2.d0*ynuynudag(i,j)*(((3.d0)*trmqyudagyu) +
     . ((3.d0)*tryudagmuyu) +
     . (1.d0*trmlynudagynu) + (trynudagmnuynu)))-
     . ((2.d0)*(anuanudagynuynudag(i,j)+ynuynudaganuanudag(i,j)+
     . anuynudagynuanudag(i,j)+ynuanudaganuynudag(i,j)+
     . anuaedagyeynudag(i,j)+
     . ynuyedagaeanudag(i,j)+anuyedagyeanudag(i,j) +
     . ynuaedagaeynudag(i,j))) +
     . (2.d0*anuanudag(i,j)*((3.d0*tryudagyu)+trynudagynu)) +
     . (2.d0*ynuynudag(i,j)*((3.d0*traudagau)+tranudaganu)) +
     . (2.d0*anuynudag(i,j)*((3.d0*traudagyu)+tranudagynu)) +
     . (2.d0*ynuanudag(i,j)*((3.d0*tryudagau)+trynudaganu)) -
     . ((((3.d0/5.d0)*alph1)+(3.d0*alph2))*((mnuynuynudag(i,j)) +
     . (2.d0*mh1sq*ynuynudag(i,j)) +
     . (2.d0*ynumlynudag(i,j)) +
     . (ynuynudagmnu(i,j)) +
     . (2.d0*anuanudag(i,j)))) -
     . (6.d0*alph2*((2.d0*M2t*M2t*ynuynudag(i,j))-
     . (M2t*anuynudag(i,j)) -
     . (2.d0*M2t*ynuanudag(i,j)))) -
     . (((6.d0/5.d0)*alph1)*((2.d0*M1t*M1t*ynuynudag(i,j)) -
     . (M1t*anuynudag(i,j))-(M1t*ynuanudag(i,j)))) 

       
       c= c+1
        enddo rgej
        enddo rgei

      
!------------------------------------------------------------------------------------------
!     109 mh1sq  RGE for Higgs Boson 
!------------------------------------------------------------------------------------------
      beta1mh1sq = gH1sq 
     .    - ((3.d0*((mh1sq*tryudagyu) +
     .     (trmqyudagyu)+ (tryudagmuyu)+ (traudagau)))) -
     .     (mh1sq*trynudagynu) - (trmlynudagynu) - (trynudagmnuynu) -
     .     (tranudaganu) 

       beta2mh1sq = ((18.d0)*((mh1sq*tryudagyuyudagyu) + 
     . (trmqyudagyuyudagyu) + (tryudagmuyuyudagyu) + (traudagauyudagyu)+
     .  traudagyuyudagau ))+ 
     . ((3.d0)*(((mh1sq+mh2sq)*tryudagyuyddagyd) + trmqyudagyuyddagyd +
     . tryudagmuyuyddagyd + tryudagyumqyddagyd + tryudagyuyddagmdyd +
     . traddagadyudagyu + tryddagydaudagau + traddagydyudagau +
     . tryddagadaudagyu)) + ((6.d0)*((mh1sq*trynudagynuynudagynu) + 
     . trmlynudagynuynudagynu + trynudagmnuynuynudagynu + 
     . tranudaganuynudagynu  + tranudagynuynudaganu)) +
     . ((mh1sq+mh2sq)*trynudagynuyedagye) +
     . trmlynudagynuyedagye + trynudagmnuynuyedagye + 
     . trynudagynumlyedagye + trynudagynuyedagmeye + 
     . traedagaeynudagynu + tryedagyeanudaganu +
     . traedagyeynudaganu + tryedagaeanudagynu -
     . (((16.d0*alph3) + ((4.d0/5.d0)*alph1))*((mh1sq*tryudagyu) +
     . trmqyudagyu + tryudagmuyu + traudagau)) - 
     . ((16.d0*alph3)*((2.d0*M3t*M3t*tryudagyu) - 
     . (M3t*tryudagau) - (M3t*traudagyu))) -
     . (((4.d0/5.d0)*alph1)*((2.d0*M1t*M1t*tryudagyu) -
     . (M1t*tryudagau) - (M1t*traudagyu))) -
     . ((3.d0/5.d0)*alph1*Sd) -
     . ((33.d0/2.d0)*alph2*alph2*M2t*M2t) -
     . ((9.d0/5.d0)*alph2*alph1*((M2t*M2t)+(M1t*M1t)+(M1t*M2t))) -
     . ((621.d0/50.d0)*alph1*alph1*M1t*M1t) -
     . (1.5d0*alph2*sig2) -
     . ((3.d0/10.d0)*alph1*sig1)
 
C--------------------------------------------------------------------------
C 110 mh2sq  

       beta1mh2sq = gH2sq 
     $      - ((3.d0*((mh2sq*tryddagyd) + (traddagad) +
     .      (trmqyddagyd)+ (tryddagmdyd))))  
     .      -(mh2sq*tryedagye) - (trmlyedagye) - (tryedagmeye) -
     .      (traedagae) 
    

      beta2mh2sq = ((18.d0)*((mh2sq*tryddagydyddagyd) + 
     . (trmqyddagydyddagyd) + (tryddagmdydyddagyd) + (traddagydyddagad)+ 
     .  (traddagadyddagyd)))+
     . ((3.d0)*(((mh1sq+mh2sq)*tryudagyuyddagyd) + trmqyudagyuyddagyd +
     . tryudagmuyuyddagyd + tryudagyumqyddagyd + tryudagyuyddagmdyd +
     . traudagauyddagyd + tryudagyuaddagad + traudagyuyddagad +
     . tryudagauaddagyd)) + 
     . ((6.d0)*((mh2sq*tryedagyeyedagye) + 
     . trmlyedagyeyedagye + tryedagmeyeyedagye + 
     . traedagaeyedagye  + traedagyeyedagae)) +
     . (((mh1sq+mh2sq)*trynudagynuyedagye ) +
     . trmlynudagynuyedagye + trynudagmnuynuyedagye + 
     . trynudagynumlyedagye + trynudagynuyedagmeye + 
     . tranudaganuyedagye + trynudagynuaedagae +
     . tranudagynuyedagae + trynudaganuaedagye) -
     . (((16.d0*alph3) - ((2.d0/5.d0)*alph1))*((mh2sq*tryddagyd) +
     . trmqyddagyd + tryddagmdyd + traddagad)) - 
     . ((16.d0*alph3)*((2.d0*M3t*M3t*tryddagyd) - 
     . (M3t*tryddagad) - (M3t*traddagyd))) +
     . (((2.d0/5.d0)*alph1)*((2.d0*M1t*M1t*tryddagyd) -
     . (M1t*tryddagad) - (M1t*traddagyd))) -
     . (((6.d0/5.d0)*alph1)*((mh2sq*tryedagye) + trmlyedagye +
     . tryedagmeye + traedagae + (2.d0*M1t*M1t*tryedagye) -
     . (M1t*traedagye) - (M1t*tryedagae))) +
     . ((3.d0/5.d0)*alph1*Sd) -
     . ((33.d0/2.d0)*alph2*alph2*M2t*M2t) -
     . ((9.d0/5.d0)*alph2*alph1*((M2t*M2t)+(M1t*M1t)+(M1t*M2t)))-
     . ((621.d0/50.d0)*alph1*alph1*M1t*M1t) -
     . (1.5d0*alph2*sig2) -
     . ((3.d0/10.d0)*alph1*sig1)

C----------------------------------------------------------------------
C     option : one loop or two loop
C-----------------------------------------------------------------------
 
      if(lopt.eq.1)then

         loopbset0: do c = 1, 9

         beta2yu(c) = 0.d0
         beta2yd(c) = 0.d0
         beta2ye(c) = 0.d0
         beta2ynu(c) = 0.d0
         beta2au(c) = 0.d0
         beta2ad(c) = 0.d0
         beta2ae(c) = 0.d0
         beta2anu(c) = 0.d0
         beta2mq(c) = 0.d0
         beta2mu(c) = 0.d0
         beta2md(c) = 0.d0
         beta2ml(c) = 0.d0
         beta2me(c) = 0.d0
         beta2mnu(c) = 0.d0
        
         

      enddo loopbset0

      beta2mh1sq = 0.d0
      beta2mh2sq = 0.d0
      
      else 

         continue
      endif


!=================================================================================================================
!                                       RGEs BEGIN
!==================================================================================================================

      do i = 1,126
         dydx(i) = 0.d0
      enddo


C-----------------------------------------------------
C     RGE FOR YUKAWAS
C-----------------------------------------------------
c      RGE FOR Yu
c---------------------------------------------------------
      k =1
     
      loopyuc:  do c = 1,9
     
       dydx(k) = beta1yu(c) + beta2yu(c)                 

      k = k+1
      enddo loopyuc



!      print*," dydx(9)", dydx(9)
         
      
c------------------------------------------------------
c  RGE FOR Yd
c------------------------------------------------------
      
      loopydc: do c=1,9

       dydx(k) = beta1yd(c) + beta2yd(c)
       
      k = k+1
      enddo loopydc           


c------------------------------------------------------
c RGE FOR Ye
c------------------------------------------------------      
     
      loopyec:     do c=1,9

      dydx(k) = beta1ye(c) + beta2ye(c)
       
      
      k = k+1
    
      enddo loopyec


c-----------------------------------------------------
c RGE FOR Ynu
c-----------------------------------------------------
     
      loopynuc:     do c=1,9
         
      dydx(k) = beta1ynu(c) + beta2ynu(c)
              
      k= k+1
      
      enddo loopynuc


C--------------------------------------------------
C     RGE FOR A-parameters
C--------------------------------------------------
C     RGE FOR Au
C--------------------------------------------------

      loopauc: do c=1,9
      
      dydx(k) = beta1au(c) + beta2au(c)
      
      k=k+1
                  
      enddo loopauc 



C-------------------------------------------------
C  RGE FOR Ad
C--------------------------------------------------

      loopadc: do c=1,9

      dydx(k) = beta1ad(c) + beta2ad(c)
      
      k=k+1
      
      enddo loopadc


C-------------------------------------------------------
C RGE FOR Ae
C--------------------------------------------------------

      loopaec: do c=1,9

      dydx(k) = beta1ae(c) + beta2ae(c)
      
      k=k+1

      enddo loopaec



C------------------------------------------------------------
C RGE FOR Anu
C-------------------------------------------------------------

      loopanuc: do c=1,9
     

      dydx(k) = beta1anu(c) + beta2anu(c)
      
      k=k+1
      
      enddo loopanuc 
      
      
C-----------------------------------------------------------------
C     RGE FOR Mass Parameters
C-----------------------------------------------------------------
C     RGE for mqsq
C-------------------------------------
      dydx(73) =   beta1mq(1) 
     $     +  beta2mq(1)

      dydx(74) =   beta1mq(2) 
     $     +  beta2mq(2)

      dydx(75) =   beta1mq(3) 
     $     +  beta2mq(3)
      
      dydx(76) =   beta1mq(5) 
     $     +  beta2mq(5)
      
      dydx(77) =   beta1mq(6) 
     $     +  beta2mq(6)

      dydx(78) =   beta1mq(9) 
     $     +  beta2mq(9)


C-------------------------------------------
C     RGE for musq
C------------------------------------------
      dydx(79) =   beta1mu(1) 
     $     +  beta2mu(1)

      dydx(80) =   beta1mu(2) 
     $     +  beta2mu(2)

      dydx(81) =   beta1mu(3) 
     $     +  beta2mu(3)

      dydx(82) =   beta1mu(5) 
     $     +  beta2mu(5)

      dydx(83) =   beta1mu(6) 
     $     +  beta2mu(6)

      dydx(84) =   beta1mu(9) 
     $     +  beta2mu(9)


C-------------------------------------------
C     RGE for mdsq
C------------------------------------------

      dydx(85) =   beta1md(1) 
     $     +  beta2md(1)

      dydx(86) =   beta1md(2) 
     $     +  beta2md(2)

      dydx(87) =   beta1md(3) 
     $     +  beta2md(3)

      dydx(88) =   beta1md(5) 
     $     +  beta2md(5)

      dydx(89) =   beta1md(6) 
     $     +  beta2md(6)

      dydx(90) =   beta1md(9) 
     $     +  beta2md(9)

C-------------------------------------------
C     RGE for mlsq
C------------------------------------------

      dydx(91) =   beta1ml(1) 
     $     +  beta2ml(1)

      dydx(92) =   beta1ml(2) 
     $     +  beta2ml(2)

      dydx(93) =   beta1ml(3) 
     $     +  beta2ml(3)

      dydx(94) =   beta1ml(5) 
     $     +  beta2ml(5)

      dydx(95) =   beta1ml(6) 
     $     +  beta2ml(6)

      dydx(96) =   beta1ml(9) 
     $     +  beta2ml(9)


C-------------------------------------------
C     RGE for mesq
C------------------------------------------

      dydx(97) =   beta1me(1) 
     $     +  beta2me(1)

      dydx(98) =   beta1me(2) 
     $     +  beta2me(2)
      
      dydx(99) =   beta1me(3) 
     $     +  beta2me(3)

      dydx(100) =  beta1me(5) 
     $     +  beta2me(5)

      dydx(101) =  beta1me(6) 
     $     +  beta2me(6)

      dydx(102) =  beta1me(9) 
     $     +  beta2me(9)


C-------------------------------------------
C     RGE for mnusq
C------------------------------------------

      dydx(103) = beta1mnu(1) + beta2mnu(1)

      dydx(104) = beta1mnu(2) + beta2mnu(2)

      dydx(105) = beta1mnu(3) + beta2mnu(3)

      dydx(106) = beta1mnu(5) + beta2mnu(5)

      dydx(107) = beta1mnu(6) + beta2mnu(6)

      dydx(108) = beta1mnu(9) + beta2mnu(9)

C-------------------------------------------
C     RGE for mh1sq and mh2sq
C------------------------------------------

      dydx(109) = beta1mh1sq + beta2mh1sq

      dydx(110) = beta1mh2sq + beta2mh2sq


C------------------------------------------------
C     RGE for \mu term 
C------------------------------------------------

      dydx(111) =  1.5d0*alph2*mmusq + 0.3d0*alph1*mmusq - !<-------------one loop
     $     ((1.5d0*tryuyudag) + (1.5d0*trydyddag) +
     $     (0.5d0*tryeyedag)+ (0.5d0*trynuynudag))*mmusq +
     $     ((4.5d0*tryuyudagyuyudag) + (4.5d0*trydyddagydyddag) + !<-------two loops begin
     $     (1.5d0*tryeyedagyeyedag) + (1.5d0*trynuynudagynuynudag) +
     $     (3.d0*tryuyddagydyudag) + (trynuyedagyeynudag))*mmusq -
     $     ((8.d0*alph3) + ((2.d0/5.d0)*alph1))*tryuyudag*mmusq -
     $     ((8.d0*alph3) - ((1.d0/5.d0)*alph1))*trydyddag*mmusq -
     $     ((3.d0/5.d0)*alph1*tryeyedag*mmusq) -
     $     ((15.d0/4.d0)*alph2*alph2*mmusq)-
     $     ((9.d0/10.d0)*alph1*alph2*mmusq) -
     $     ((207.d0/100.d0)*alph1*alph1*mmusq)

C-------------------------------------------------------
C     RGE for b_\mu term
C-------------------------------------------------------

      dydx(112) = ((1.5d0*alph2) + (0.3d0*alph1))*bmu - !<---------------one loop b_\mu terms
     $     ((1.5d0*tryuyudag) + (1.5d0*trydyddag) + 
     $     (0.5d0*tryeyedag) + (0.5d0*trynuynudag))*bmu -     
     $     ((3.d0*trauyudag) + (3.d0*tradyddag) + (traeyedag) + !<---------------one loop \mu terms
     $     (tranuynudag))*mmusq - 
     $     ((3.d0*alph2*M2t) + ((3.d0/5.d0)*alph1*M1t))*mmusq +     
     $     ((4.5d0*tryuyudagyuyudag) + (4.5d0*trydyddagydyddag) + !<-------two loop b_\mu terms begin
     $     (1.5d0*tryeyedagyeyedag) + (1.5d0*trynuynudagynuynudag) +
     $     (3.d0*tryuyddagydyudag) + (trynuyedagyeynudag))*bmu -
     $     ((8.d0*alph3) + (2.d0/5.d0)*alph1)*tryuyudag*bmu -
     $     ((8.d0*alph3) - (1.d0/5.d0)*alph1)*trydyddag*bmu -
     $     ((3.d0/5.d0)*alph1*tryeyedag)*bmu -
     $     ((15.d0/4.d0)*alph2*alph2*bmu)-
     $     ((9.d0/10.d0)*alph1*alph2*bmu) -
     $     ((207.d0/100.d0)*alph1*alph1*bmu) + 
     $     (2.d0*mmusq*((9.d0*trauyudagyuyudag) + !<------------two loop \mu terms begin
     $     (9.d0*tradyddagydyddag) + (3.d0*tradyudagyuyddag) +
     $     (3.d0*trauyddagydyudag) + (3.d0*tranuynudagynuynudag) +
     $     (3.d0*traeyedagyeyedag) + (tranuyedagyeynudag) +
     $     (traeynudagynuyedag))) -
     $     ((16.d0*alph3) + (4.d0/5.d0)*alph1)*trauyudag*mmusq -
     $     ((16.d0*alph3) - (2.d0/5.d0)*alph1)*tradyddag*mmusq - 
     $     ((6.d0/5.d0)*alph1*traeyedag)*mmusq +
     $     ((16.d0*alph3*M3t)+((4.d0/5.d0)*alph1*M1t))*tryuyudag*mmusq+
     $     ((16.d0*alph3*M3t)-((2.d0/5.d0)*alph1*M1t))*trydyddag*mmusq+
     $     ((6.d0/5.d0)*alph1*M1t)*tryeyedag*mmusq +
     $     (15.d0*alph2*alph2*M2t)*mmusq +
     $     ((9.d0/5.d0)*alph1*alph2*(M1t+M2t))*mmusq +
     $     ((207.d0/25.d0)*alph1*alph1*M1t)*mmusq

C------------------------------------------------------------------------------------------


      dydx(113) = 0.5d0*(4.d0*gnut - 6.d0*tryuyudag)*mnu_light(1,1)
     .     -0.5d0*((yeyedag(1,1)+ynuynudag(1,1))*mnu_light(1,1)
     .     + (yeyedag(1,2)+ynuynudag(1,2))*mnu_light(2,1)
     .     + (yeyedag(1,3)+ynuynudag(1,3))*mnu_light(3,1)
     .     + mnu_light(1,1)*(yedagye(1,1)+ynudagynu(1,1))
     .     + mnu_light(1,2)*(yedagye(2,1)+ynudagynu(2,1))
     .     + mnu_light(1,3)*(yedagye(3,1)+ynudagynu(3,1)) )

      dydx(114) = 0.5d0*(4.d0*gnut - 6.d0*tryuyudag)*mnu_light(1,2)
     .     -0.5d0*(ynuynudag(1,1)*mnu_light(1,2)
     .     + ynuynudag(1,2)*mnu_light(2,2)
     .     + ynuynudag(1,3)*mnu_light(3,2)
     .     + mnu_light(1,1)*ynudagynu(1,2)
     .     + mnu_light(1,2)*ynudagynu(2,2)
     .     + mnu_light(1,3)*ynudagynu(3,2) )
     .     -0.5d0*(yeyedag(1,1)*mnu_light(1,2)
     .     + yeyedag(1,2)*mnu_light(2,2)
     .     + yeyedag(1,3)*mnu_light(3,2)
     .     + mnu_light(1,1)*yedagye(1,2)
     .     + mnu_light(1,2)*yedagye(2,2)
     .     + mnu_light(1,3)*yedagye(3,2) )


      dydx(115) = 0.5d0*(4.d0*gnut - 6.d0*tryuyudag)*mnu_light(1,3)
     .     -0.5d0*(ynuynudag(1,1)*mnu_light(1,3)
     .     + ynuynudag(1,2)*mnu_light(2,3)
     .     + ynuynudag(1,3)*mnu_light(3,3)
     .     + mnu_light(1,1)*ynudagynu(1,3)
     .     + mnu_light(1,2)*ynudagynu(2,3)
     .     + mnu_light(1,3)*ynudagynu(3,3) )
     .     -0.5d0*(yeyedag(1,1)*mnu_light(1,3)
     .     + yeyedag(1,2)*mnu_light(2,3)
     .     + yeyedag(1,3)*mnu_light(3,3)
     .     + mnu_light(1,1)*yedagye(1,3)
     .     + mnu_light(1,2)*yedagye(2,3)
     .     + mnu_light(1,3)*yedagye(3,3) )


      dydx(116) = 0.5d0*(4.d0*gnut - 6.d0*tryuyudag)*mnu_light(2,2)
     .     -0.5d0*(ynuynudag(2,1)*mnu_light(1,2)
     .     + ynuynudag(2,2)*mnu_light(2,2)
     .     + ynuynudag(2,3)*mnu_light(3,2)
     .     + mnu_light(2,1)*ynudagynu(1,2)
     .     + mnu_light(2,2)*ynudagynu(2,2)
     .     + mnu_light(2,3)*ynudagynu(3,2) )
     .     -0.5d0*(yeyedag(2,1)*mnu_light(1,2)
     .     + yeyedag(2,2)*mnu_light(2,2)
     .     + yeyedag(2,3)*mnu_light(3,2)
     .     + mnu_light(2,1)*yedagye(1,2)
     .     + mnu_light(2,2)*yedagye(2,2)
     .     + mnu_light(2,3)*yedagye(3,2) )


      dydx(117) = 0.5d0*(4.d0*gnut - 6.d0*tryuyudag)*mnu_light(2,3)
     .     -0.5d0*(ynuynudag(2,1)*mnu_light(1,3)
     .     + ynuynudag(2,2)*mnu_light(2,3)
     .     + ynuynudag(2,3)*mnu_light(3,3)
     .     + mnu_light(2,1)*ynudagynu(1,3)
     .     + mnu_light(2,2)*ynudagynu(2,3)
     .     + mnu_light(2,3)*ynudagynu(3,3) )
     .     -0.5d0*(yeyedag(2,1)*mnu_light(1,3)
     .     + yeyedag(2,2)*mnu_light(2,3)
     .     + yeyedag(2,3)*mnu_light(3,3)
     .     + mnu_light(2,1)*yedagye(1,3)
     .     + mnu_light(2,2)*yedagye(2,3)
     .     + mnu_light(2,3)*yedagye(3,3) )

      dydx(118) = 0.5d0*(4.d0*gnut - 6.d0*tryuyudag)*mnu_light(3,3)
     .     -0.5d0*(ynuynudag(3,1)*mnu_light(1,3)
     .     + ynuynudag(3,2)*mnu_light(2,3)
     .     + ynuynudag(3,3)*mnu_light(3,3)
     .     + mnu_light(3,1)*ynudagynu(1,3)
     .     + mnu_light(3,2)*ynudagynu(2,3)
     .     + mnu_light(3,3)*ynudagynu(3,3))
     .     -0.5d0*(yeyedag(3,1)*mnu_light(1,3)
     .     + yeyedag(3,2)*mnu_light(2,3)
     .     + yeyedag(3,3)*mnu_light(3,3)
     .     + mnu_light(3,1)*yedagye(1,3)
     .     + mnu_light(3,2)*yedagye(2,3)
     .     + mnu_light(3,3)*yedagye(3,3))

!-------------------------------------------------------------------------------

      dydx(113) = 0.d0
      dydx(114) = 0.d0
      dydx(115) = 0.d0
      dydx(116) = 0.d0
      dydx(117) = 0.d0
      dydx(118) = 0.d0


!-------------------------------------------------------------------------------
C     RGEs FOR GAUGE COUPLINGS
C-------------------------------------------------------------------------------
      if(lopt.eq.2)then

      dydx(119)= -((alph3**(2.d0))*b3)
     &     - ((alph3**(2.d0)*
     &     (((11.d0/5.d0)*alph1) + 
     &     ((9.d0)*alph2)+
     &     ((14.d0)*alph3)-((4.d0)*tryudagyu)-((4.d0)*tryddagyd))))



      dydx(120)= -((alph2**(2.d0))*b2)
     &     - ((alph2**(2.d0)*
     &     (((9.d0/5.d0)*alph1) + 
     &     ((25.d0)*alph2) +
     &     ((24.d0)*alph3)
     &     -  ((6.d0)*tryudagyu)-((6.d0)*tryddagyd)-((2.d0)*tryedagye)
     &     - ((2.d0)*trynudagynu))))


      dydx(121)= -((alph1**(2.d0))*b1)
     &     -((alph1**(2.d0)*
     &     (((199.d0/25.d0)*alph1) + 
     &     ((27.d0/5.d0)*alph2) + ((88.d0/5.d0)*alph3)
     &     -((26.d0/5.d0)*tryudagyu)-((14.d0/5.d0)*tryddagyd)
     &     -((18.d0/5.d0)*tryedagye)-((6.d0/5.d0)*trynudagynu))))

      else

      dydx(119)= -((alph3**(2.d0))*b3)
c$$$     &     - ((alph3**(2.d0)*
c$$$     &     (((11.d0/5.d0)*alph1) + 
c$$$     &     ((9.d0)*alph2)+
c$$$     &     ((14.d0)*alph3)-((4.d0)*tryudagyu)-((4.d0)*tryddagyd))))



      dydx(120)= -((alph2**(2.d0))*b2)
c$$$     &     - ((alph2**(2.d0)*
c$$$     &     (((9.d0/5.d0)*alph1) + 
c$$$     &     ((25.d0)*alph2) +
c$$$     &     ((24.d0)*alph3)
c$$$     &     -  ((6.d0)*tryudagyu)-((6.d0)*tryddagyd)-((2.d0)*tryedagye)
c$$$     &     - ((2.d0)*trynudagynu))))


      dydx(121)= -((alph1**(2.d0))*b1)
c$$$     &     -((alph1**(2.d0)*
c$$$     &     (((199.d0/25.d0)*alph1) + 
c$$$     &     ((27.d0/5.d0)*alph2) + ((88.d0/5.d0)*alph3)
c$$$     &     -((26.d0/5.d0)*tryudagyu)-((14.d0/5.d0)*tryddagyd)
c$$$     &     -((18.d0/5.d0)*tryedagye)-((6.d0/5.d0)*trynudagynu))))

      endif
!-------------------------------------------
C     RGE FOR GAUGINOs
!-------------------------------------------
      
      if(lopt.eq.2)then

      dydx(122)=   -(b1*alph1*M1t)
     .     - ((alph1)*(((398.d0/25.d0)*alph1*M1t)    +
     .     ((27.d0/5.d0)*alph2*(M1t+M2t)) +
     .     ((88.d0/5.d0)*alph3*(M1t+M3t)) +
     .     ((26.d0/5.d0)*(tryudagau-(M1t*tryudagyu))) +
     .     ((14.d0/5.d0)*(tryddagad-(M1t*tryddagyd))) +
     .     ((18.d0/5.d0)*(tryedagae-(M1t*tryedagye)))
     .     +  ((6.d0/5.d0)*(trynudaganu-(M1t*trynudagynu)))))


      dydx(123)=    -(b2*alph2*M2t)
     .     - ((alph2)*(((9.d0/5.d0)*alph1*(M2t+M1t)) +
     .     (50.d0*alph2*M2t)             +
     .     (24.d0*alph3*(M2t+M3t))       +
     .     (6.d0*(tryudagau-(M2t*tryudagyu))) +
     .     (6.d0*(tryddagad-(M2t*tryddagyd))) +
     .     (2.d0*(tryedagae-(M2t*tryedagye)))
     .     + (2.d0*(trynudaganu-(M2t*trynudagynu)))))


      dydx(124)=   -(b3*alph3*M3t)
     .     - ((alph3)*(((11.d0/5.d0)*alph1*(M3t+M1t)) +
     .     (9.d0*alph2*(M3t+M2t))         +
     .     (28.d0*alph3*M3t)              +
     .     (4.d0*(tryudagau-(M3t*tryudagyu))) +
     .     (4.d0*(tryddagad-(M3t*tryddagyd)))))

      else

      dydx(122)=   -(b1*alph1*M1t)
c$$$     .     - ((alph1)*(((398.d0/25.d0)*alph1*M1t)    +
c$$$     .     ((27.d0/5.d0)*alph2*(M1t+M2t)) +
c$$$     .     ((88.d0/5.d0)*alph3*(M1t+M3t)) +
c$$$     .     ((26.d0/5.d0)*(tryudagau-(M1t*tryudagyu))) +
c$$$     .     ((14.d0/5.d0)*(tryddagad-(M1t*tryddagyd))) +
c$$$     .     ((18.d0/5.d0)*(tryedagae-(M1t*tryedagye)))
c$$$     .     +  ((6.d0/5.d0)*(trynudaganu-(M1t*trynudagynu)))))


      dydx(123)=    -(b2*alph2*M2t)
c$$$     .     - ((alph2)*(((9.d0/5.d0)*alph1*(M2t+M1t)) +
c$$$     .     (50.d0*alph2*M2t)             +
c$$$     .     (24.d0*alph3*(M2t+M3t))       +
c$$$     .     (6.d0*(tryudagau-(M2t*tryudagyu))) +
c$$$     .     (6.d0*(tryddagad-(M2t*tryddagyd))) +
c$$$     .     (2.d0*(tryedagae-(M2t*tryedagye)))
c$$$     .     + (2.d0*(trynudaganu-(M2t*trynudagynu)))))


      dydx(124)=   -(b3*alph3*M3t)
c$$$     .     - ((alph3)*(((11.d0/5.d0)*alph1*(M3t+M1t)) +
c$$$     .     (9.d0*alph2*(M3t+M2t))         +
c$$$     .     (28.d0*alph3*M3t)              +
c$$$     .     (4.d0*(tryudagau-(M3t*tryudagyu))) +
c$$$     .     (4.d0*(tryddagad-(M3t*tryddagyd)))))

      endif

!-------------------------------------------------------------------
C     RGEs for vev
!-------------------------------------------------------------------
     

      dydx(125) = ((-3.d0/8.d0) * ((alph1/5.d0) + alph2) + 
     $     1.5d0 * trydyddag + 0.5d0 * tryeyedag) * vev1   +
     $     (( - (3.d0/8.d0) *                                         !<------ 2 loop starts here
     $     (3.d0 * trydyddagydyddag + 3.d0 * trydyddagyuyudag + 
     $     tryeyedagyeyedag)) + ((0.4d0 * alph1 + 4.5d0 * alph2 + 
     $     20.d0 * alph3) * trydyddag/2.d0) + ((1.8d0 * alph1 + 
     $     1.5d0 * alph2) * tryeyedag/2.d0) + ((5967.d0/3200.d0) * 
     $     alph1 * alph1) + ((1485.d0/128.d0) * alph2 * alph2) +
     $     ((27.d0/80.d0) * alph1 * alph2)) * vev1


      dydx(126) = ((-3.d0/8.d0) * ((alph1/5.d0) + alph2) + 
     $     1.5d0*tryuyudag) * vev2 +
     $     (( - (3.d0/8.d0) * (3.d0 * tryuyudagyuyudag + 3.d0 *       !<------- 2 loop starts here
     $     tryuyudagydyddag)) + ((1.9d0 * alph1 + 4.5d0 * alph2 + 
     $     20.d0 * alph3) * tryuyudag/2.d0) + ((5967.d0/3200.d0) * 
     $     alph1 * alph1) + ((1485.d0/128.d0) * alph2 * alph2) +
     $     ((27.d0/80.d0) * alph1 * alph2)) * vev2


      
!-----------------------------------------------------------------------------
C     RGES END
C-----------------------------------------------------------------------------

      
      e = (MX/(dexp(t/2.d0)))
      
      a1 = 1.d0/(alph1*4.d0*pi)
      a2 = 1.d0/(alph2*4.d0*pi)
      a3 = 1.d0/(alph3*4.d0*pi)
      

      a1_unif = a1
      a2_unif = a2
            
      r = (1.d0 - (MIN(a1_unif,a2_unif)/
     $     MAX(a1_unif,a2_unif)))


      
!------------------------------------------------------

      tol  = 0.001d0
      tol1 = 10.d0
      rc = 0

!      fuscale = 0

      if(model.eq.'mSUG'.or.model.eq.'NUHM'
     $     .or.model.eq.'NUGM'.or.model.eq.'CNUM')then

      up: if(M1t.eq.0.d0)then
         
         unif0: if(a1.ge.a2.and.r.lt.5.d-2)then
            
!     unif0: if(a1.le.a2.and.r.lt.5.d-2.and.fuscale.eq.0)then
!     unif0: if(a1.le.a2.and.fuscale.eq.0)then
!     unif0: if(alph1.ge.alph2.and.fuscale.eq.0)then
            
!     unif0: if(e.lt.2.3d16.and.e.gt.1.8d16.and.fuscale.eq.0)then
            
!     45         format(1x,1pg25.15)
!     write(15,45) e

            e1 = 0.d0
            e1 = e 
            
!            print*," unification scale ", e1
            
c$$$        print*,"y(121), y(120), y(119) = ", yy(121)*16.d0*pi*pi, 
c$$$     $       yy(120)*16.d0*pi*pi, yy(119)*16.d0*pi*pi
c$$$
c$$$        print*,"dydx(121), dydx(120), dydx(119) = ", dydx(121)*16.d0*pi
c$$$     $       *pi, dydx(120)*16.d0*pi*pi, dydx(119)*16.d0*pi*pi
c$$$	
c$$$        print*,"factor = ", dexp((yy(120) - yy(121))/(dydx(121) - 
c$$$     $       dydx(120)))

            e_next = e 

            looppick: do i = 1, 126
            yukgut(i) =  yy(i)
            enddo looppick     
         
!         fuscale = 1
         
         return
         
      else unif0              
         
         continue 
         
      endif unif0
      
      endif up
      
      else
         continue
      endif
      
      RETURN
      
      end subroutine mssmrge


!===============================================================================================
C                                   MSSM RGE ENDS
!==============================================================================================
