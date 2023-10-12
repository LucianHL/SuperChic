*
* subtractive mitchell-moore generator
* ronald kleiss - october 2, 1987
*
* the algorithm is N(i)=[ N(i-24) - N(i-55) ]mod M,
* implemented in a cirucular array with identifcation
* of NR(i+55) and nr(i), such that effectively:
*        N(1)   <---   N(32) - N(1)
*        N(2)   <---   N(33) - N(2)  ....
*   .... N(24)  <---   N(55) - N(24)
*        N(25)  <---   N(1)  - N(25) ....
*   .... N(54)  <---   N(30) - N(54)
*        N(55)  <---   N(31) - N(55)
*
* in this version  M =2**30  and  RM=1/M=2.D0**(-30.D0)
*
* the array NR has been initialized by putting NR(i)=i
* and subsequently running the algorithm 100,000 times.
*

      subroutine R2455(Ran)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION N(55)
      DATA N/
     . 980629335, 889272121, 422278310,1042669295, 531256381,
     . 335028099,  47160432, 788808135, 660624592, 793263632,
     . 998900570, 470796980, 327436767, 287473989, 119515078,
     . 575143087, 922274831,  21914605, 923291707, 753782759,
     . 254480986, 816423843, 931542684, 993691006, 343157264,
     . 272972469, 733687879, 468941742, 444207473, 896089285,
     . 629371118, 892845902, 163581912, 861580190,  85601059,
     . 899226806, 438711780, 921057966, 794646776, 417139730,
     . 343610085, 737162282,1024718389,  65196680, 954338580,
     . 642649958, 240238978, 722544540, 281483031,1024570269,
     . 602730138, 915220349, 651571385, 405259519, 145115737/
      DATA M/1073741824/
      DATA RM/0.9313225746154785D-09/
      DATA K/55/,L/31/
 777  IF(K.EQ.55) THEN
         K=1
      ELSE
         K=K+1
      ENDIF
      IF(L.EQ.55) THEN
         L=1
      ELSE
         L=L+1
      ENDIF
      J=N(L)-N(K)
      IF(J.LT.0) J=J+M
      N(K)=J
      RAN=J*RM
      if(ran.eq.0d0)goto 777
      END
