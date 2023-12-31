      SUBROUTINE SIMPW(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C SIMPS,T,U,V
C SIMPS
C A1,B1 -THE LIMITS OF INTEGRATION
C H1 -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN AEPS1=REPS1=0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(7),P(5)
      INTEGER K
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.
    4 X0=X
      IF((X0+4.*H-B)*S) 5,5,6
    6 H=(B-X0)/4.
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*ABS(F(K))
      DI1=(F(1)+4.*F(3)+F(5))*2.*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.) 17,14,14
   17 H=2.*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
      END

      SUBROUTINE SIMPT(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C SIMPT
C A1,B1 -THE LIMITS OF INTEGRATION
C H1 -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN AEPS1=REPS1=0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(7),P(5)
      INTEGER K      
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.
    4 X0=X
      IF((X0+4.*H-B)*S) 5,5,6
    6 H=(B-X0)/4.
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*ABS(F(K))
      DI1=(F(1)+4.*F(3)+F(5))*2.*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.) 17,14,14
   17 H=2.*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
      END

      SUBROUTINE SIMPU(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C SIMPU
C A1,B1 -THE LIMITS OF INTEGRATION
C H1 -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN AEPS1=REPS1=0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER K
      DIMENSION F(7),P(5)
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.
    4 X0=X
      IF((X0+4.*H-B)*S) 5,5,6
    6 H=(B-X0)/4.
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*ABS(F(K))
      DI1=(F(1)+4.*F(3)+F(5))*2.*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.) 17,14,14
   17 H=2.*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
      END

      SUBROUTINE SIMPV(A1,B1,H1,REPS1,AEPS1,FUNCT,X,AI,AIH,AIABS)
C SIMPV
C A1,B1 -THE LIMITS OF INTEGRATION
C H1 -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN AEPS1=REPS1=0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(7),P(5)
      INTEGER K      
      H=DSIGN(H1,B1-A1)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(X)/3.
    4 X0=X
      IF((X0+4.*H-B)*S) 5,5,6
    6 H=(B-X0)/4.
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(X)/3.
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*ABS(F(K))
      DI1=(F(1)+4.*F(3)+F(5))*2.*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.) 17,14,14
   17 H=2.*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
      END

      SUBROUTINE FDSIMP (AA1,BB1,HH1,REPS1,AEPS1,FUNCT,DFUN,DFUNIN,
     + DUMMY,AI,AIH,AIABS)
*=======================================================================
C        B1
C AI=INT {FUNCT(X)*DDFUN(X)Y}
C        A1
C A1,B1 -THE LIMITS OF INTEGRATION
C H1    -AN INITIAL STEP OF INTEGRATION
C REPS1,AEPS1 - RELATIVE AND ABSOLUTE PRECISION OF INTEGRATION
C FUNCT -A NAME OF FUNCTION SUBPROGRAM FOR CALCULATION OF INTEGRAND +
C X - AN ARGUMENT OF THE INTEGRAND
C DFUNIN - INVERSE ( DFUN ). SHOULD BE DFUNINDFUN(X)Y=X.
C AI - THE VALUE OF INTEGRAL
C AIH- THE VALUE OF INTEGRAL WITH THE STEP OF INTEGRATION
C AIABS- THE VALUE OF INTEGRAL FOR MODULE OF THE INTEGRAND
C THIS SUBROGRAM CALCULATES THE DEFINITE INTEGRAL WITH THE RELATIVE OR
C ABSOLUTE PRECISION BY SIMPSON+S METHOD WITH THE AUTOMATICAL CHOICE
C OF THE STEP OF INTEGRATION
C IF AEPS1    IS VERY SMALL(LIKE 1.E-17),THEN CALCULATION OF INTEGRAL
C WITH REPS1,AND IF REPS1 IS VERY SMALL (LIKE 1.E-10),THEN CALCULATION
C OF INTEGRAL WITH AEPS1
C WHEN  AEPS1=REPS1=0. THEN CALCULATION WITH THE CONSTANT STEP H1
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION F(7),P(5)
      INTEGER K
      DUM=DUMMY
      H1=(DFUN(BB1)-DFUN(AA1))/(BB1-AA1+1.654876596E-20)*HH1
      A1=DFUN(AA1)
      B1=DFUN(BB1)
      H=DSIGN(H1,B1-A1+1.654876596E-20)
      S=DSIGN(1.D0,H)
      A=A1
      B=B1
      AI=0.D0
      AIH=0.D0
      AIABS=0.D0
      P(2)=4.D0
      P(4)=4.D0
      P(3)=2.D0
      P(5)=1.D0
      IF(B-A) 1,2,1
    1 REPS=DABS(REPS1)
      AEPS=DABS(AEPS1)
      DO 3 K=1,7
  3   F(K)=10.D16
      X=A
      C=0.D0
      F(1)=FUNCT(DFUNIN(X))/3.
    4 X0=X
      IF((X0+4.*H-B)*S) 5,5,6
    6 H=(B-X0)/4.
      IF(H) 7,2,7
    7 DO 8 K=2,7
  8   F(K)=10.D16
      C=1.D0
    5 DI2=F(1)
      DI3=DABS(F(1))
      DO 9 K=2,5
      X=X+H
      IF((X-B)*S) 23,24,24
   24 X=B
   23 IF(F(K)-10.D16) 10,11,10
   11 F(K)=FUNCT(DFUNIN(X))/3.
   10 DI2=DI2+P(K)*F(K)
    9 DI3=DI3+P(K)*ABS(F(K))
      DI1=(F(1)+4.*F(3)+F(5))*2.*H
      DI2=DI2*H
      DI3=DI3*H
      IF(REPS) 12,13,12
   13 IF(AEPS) 12,14,12
   12 EPS=DABS((AIABS+DI3)*REPS)
      IF(EPS-AEPS) 15,16,16
   15 EPS=AEPS
   16 DELTA=DABS(DI2-DI1)
      IF(DELTA-EPS) 20,21,21
   20 IF(DELTA-EPS/8.) 17,14,14
   17 H=2.*H
      F(1)=F(5)
      F(2)=F(6)
      F(3)=F(7)
      DO 19 K=4,7
  19  F(K)=10.D16
      GO TO 18
   14 F(1)=F(5)
      F(3)=F(6)
      F(5)=F(7)
      F(2)=10.D16
      F(4)=10.D16
      F(6)=10.D16
      F(7)=10.D16
   18 DI1=DI2+(DI2-DI1)/15.
      AI=AI+DI1
      AIH=AIH+DI2
      AIABS=AIABS+DI3
      GO TO 22
   21 H=H/2.
      F(7)=F(5)
      F(6)=F(4)
      F(5)=F(3)
      F(3)=F(2)
      F(2)=10.D16
      F(4)=10.D16
      X=X0
      C=0.D0
      GO TO 5
   22 IF(C) 2,4,2
    2 RETURN
      END

      FUNCTION SPENCE(X)
*
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (F1=1.64493406684822618D0)
*
      IF(X)8,1,1
1     IF(X-.5D0)2,2,3
2     SPENCE=FSPENS(X)
      RETURN
3     IF(X-1.D0)4,4,5
4     SPENCE=F1-LOG(X)*LOG(1D0-X+1D-15)-FSPENS(1D0-X)
      RETURN
5     IF(X-2.D0)6,6,7
6     SPENCE=F1-.5D0*LOG(X)*LOG((X-1D0)**2/X)+FSPENS(1D0-1D0/X)
      RETURN
7     SPENCE=2D0*F1-.5D0*LOG(X)**2-FSPENS(1D0/X)
      RETURN
8     IF(X+1.D0)10,9,9
9     SPENCE=-.5D0*LOG(1D0-X)**2-FSPENS(X/(X-1D0))
      RETURN
10    SPENCE=-.5D0*LOG(1D0-X)*LOG(X**2/(1D0-X))-F1+FSPENS(1D0/(1D0-X))
      END

      FUNCTION FSPENS(X)
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      A=1D0
      F=0D0
      AN=0D0
      TCH=1D-20
1     AN=AN+1D0
      A=A*X
      B=A/AN**2
      F=F+B
      IF(B-TCH)2,2,1
2     FSPENS=F
      END

      DOUBLE PRECISION FUNCTION DDILOG(X)
*=======================================================================

      DOUBLE PRECISION X,Y,T,S,A,PI3,PI6,ZERO,ONE,HALF,MALF,MONE,MTWO
      DOUBLE PRECISION C(0:18),H,ALFA,B0,B1,B2
      INTEGER I
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

c     FUNCTION DLI2(X)
c     IMPLICIT REAL*8(A-Z)
c     DLI2=DDILOG(X)
c     END

c     FUNCTION DLI3(X)
* This is Li3(x) for x<=1 only!
c     IMPLICIT REAL*8(A-Z)
c     DATA ZET2/1.6449340668482262D0/
c     IF(X.LT.0D0) THEN
c     X1=1D0-X
c     DLI3=+S12(1D0/X1)-TRILOG(1D0/X1)
c    &     -LOG(X1)**2*LOG(-X)/2D0+LOG(X1)**3/3D0
c    &     +(LOG(-X)-LOG(X1))*DDILOG(X)-ZET2*LOG(X1)
c    &     -LOG(X1)*LOG(-X)*LOG(X1/(-X))/2D0
c     ELSE
c     DLI3=TRILOG(X)
c     ENDIF
c     END

c     FUNCTION DS12(X)
* This is S12(x) for x<=1 only!
c     IMPLICIT REAL*8(A-Z)
c     DATA ZET2/1.6449340668482262D0/
c     DATA ZET3/1.2020569031595943D0/
c     IF(X.LT.0D0) THEN
c     X1=1D0-X
c     DS12=-TRILOG(1D0/X1)-LOG(X1)**2*LOG(-X)/2D0+LOG(X1)**3/6D0
c    &     -LOG(X1)*DDILOG(X)+ZET3-ZET2*LOG(X1)
c     ELSE
c     DS12=S12(X)
c     ENDIF
c     END

      DOUBLE PRECISION FUNCTION TRILOG(X)
C Tsuyoshi Matsuura 1987
C
C     TRILOG: Li3   between 0 and 1  !!!
C
      DOUBLE PRECISION X,S(0:10),L(0:20),U,Z,HELP,Z3
      DOUBLE PRECISION DDILOG
      INTEGER I
C
C     S: COEFFICIENTS OF S1,2
C     L: COEFFICIENTS OF Li3
C
      DATA S/0.5000000000000000D+00,2.0833333333333333D-02,
     1      -2.3148148148148148D-04,4.1335978835978837D-06,
     2      -8.2671957671957673D-08,1.7397297489890083D-09,
     3      -3.7744215276339238D-11,8.3640853316779243D-13,
     4      -1.8831557201792127D-14,4.2930310281389223D-16,
     5      -9.8857668116275541D-18/
      DATA L/1.0000000000000000D+00,-0.3750000000000000D+00,
     1       7.8703703703703705E-02,-8.6805555555555557E-03,
     2       1.2962962962962963E-04, 8.1018518518518520E-05,
     3      -3.4193571608537598E-06,-1.3286564625850340E-06,
     4       8.6608717561098512E-08, 2.5260875955320400E-08,
     5      -2.1446944683640648E-09,-5.1401106220129790E-10,
     6       5.2495821146008296E-11, 1.0887754406636318E-11,
     7      -1.2779396094493696E-12,-2.3698241773087452E-13,
     8       3.1043578879654624E-14, 5.2617586299125060E-15,
     9      -7.5384795499492655E-16,-1.1862322577752285E-16,
     1       1.8316979965491384E-17/
      DATA Z3/1.2020569031595943D+00/
C
      IF(X .LT. 0D0  .OR.  X .GT. 1D0) THEN
         WRITE(*,*)' **************************************'
         WRITE(*,*)' Li3 called with X = ',X
         WRITE(*,*)' This should lie between 0 and 1 !!!'
         WRITE(*,*)' The program stops right here !!!'
         STOP 1
      ENDIF
      IF (X.LT.0.5D0) THEN
         IF (X.GT.0.0D0) THEN
            U=-DLOG(1.0D0-X)
            HELP=U*L(20)+L(19)
            DO 10 I=18,0,-1
               HELP=U*HELP+L(I)
   10       CONTINUE
            TRILOG=U*HELP
         ELSE
            TRILOG=0.0D0
         ENDIF
      ELSE
         IF (X.LT.1.0D0) THEN
            U=-DLOG(X)
            Z=U*U
            HELP=Z*S(10)+S(9)
            DO 20 I=8,0,-1
               HELP=Z*HELP+S(I)
   20       CONTINUE
            HELP=1.0D0/2.0D0*(Z*HELP-Z*U/6.0D0)
c           LI3=-HELP+DLOG(X)*DILOG(X)+0.5D0*DLOG(1D0-X)*DLOG(X)**2+Z3
        TRILOG=-HELP+DLOG(X)*DDILOG(X)+0.5D0*DLOG(1D0-X)*DLOG(X)**2+Z3
         ELSE
          TRILOG=Z3
         ENDIF
      ENDIF
      END

      DOUBLE PRECISION FUNCTION TRIS12(X)
C Tsuyoshi Matsuura 1987
C
C     S1,2   between 0 and 1  !!!
C
      INTEGER I
      REAL*8 X,S(0:10),L(0:20),U,Z,HELP,Z3
      REAL*8 DDILOG
C
C     S: COEFFICIENTS OF S1,2
C     L: COEFFICIENTS OF Li3
C
      DATA S/0.5000000000000000D+00,2.0833333333333333D-02,
     1      -2.3148148148148148D-04,4.1335978835978837D-06,
     2      -8.2671957671957673D-08,1.7397297489890083D-09,
     3      -3.7744215276339238D-11,8.3640853316779243D-13,
     4      -1.8831557201792127D-14,4.2930310281389223D-16,
     5      -9.8857668116275541D-18/
      DATA L/1.0000000000000000D+00,-0.3750000000000000D+00,
     1       7.8703703703703705E-02,-8.6805555555555557E-03,
     2       1.2962962962962963E-04, 8.1018518518518520E-05,
     3      -3.4193571608537598E-06,-1.3286564625850340E-06,
     4       8.6608717561098512E-08, 2.5260875955320400E-08,
     5      -2.1446944683640648E-09,-5.1401106220129790E-10,
     6       5.2495821146008296E-11, 1.0887754406636318E-11,
     7      -1.2779396094493696E-12,-2.3698241773087452E-13,
     8       3.1043578879654624E-14, 5.2617586299125060E-15,
     9      -7.5384795499492655E-16,-1.1862322577752285E-16,
     1       1.8316979965491384E-17/
      DATA Z3/1.2020569031595943D+00/
C
      IF(X .LT. 0D0  .OR.  X .GT. 1D0) THEN
         WRITE(*,*)' **************************************'
         WRITE(*,*)' TRIS12 called with X = ',X
         WRITE(*,*)' This should lie between 0 and 1 !!!'
         WRITE(*,*)' The program stops right here !!!'
         STOP 1
      ENDIF
      IF (X.LT.0.5D0) THEN
         IF (X.GT.0.0D0) THEN
            U=-DLOG(1.0D0-X)
            Z=U*U
            HELP=Z*S(10)+S(9)
            DO 10 I=8,0,-1
               HELP=HELP*Z+S(I)
   10       CONTINUE
            TRIS12=1.0D0/2.0D0*(Z*HELP-Z*U/6.0D0)
         ELSE
            TRIS12=0.0D0
         ENDIF
      ELSE
         IF (X.LT.1.0D0) THEN
            U=-DLOG(X)
            HELP=U*L(20)+L(19)
            DO 20 I=18,0,-1
               HELP=HELP*U+L(I)
   20       CONTINUE
            HELP=U*HELP
            TRIS12=-HELP+DLOG(1.0D0-X)*DDILOG(1.0D0-X)+
     +           0.5D0*DLOG(X)*DLOG(1.0D0-X)**2+Z3
         ELSE
            TRIS12=Z3
         ENDIF
      ENDIF
      END

      FUNCTION DLI2(X)
      IMPLICIT REAL*8(A-Z)
      DLI2=DDILOG(X)
      END

      FUNCTION DLI3(X)
* This is Li3(x) for x<=1 only!
      IMPLICIT REAL*8(A-Z)
      DATA ZET2/1.6449340668482262D0/
      DATA ZET3/1.2020569031595943D0/
      IF(X.GT.0D0) THEN
        DLI3=TRILOG(X)
      ELSE IF(X.EQ.0D0) THEN
        DLI3=0D0
      ELSE IF(X.GT.-1D0) THEN
        X1=1D0-X
        DLOGX1=DLOG(X1)
        XX=-X/X1
        DLI3=+TRIS12(XX)-TRILOG(XX)
     &       +DLOGX1*(DDILOG(X)+1D0/3*DLOGX1*DLOGX1)
      ELSE IF(X.EQ.-1D0) THEN
        DLI3=-3D0/4*ZET3
      ELSE
        X1=1D0-X
        DLOGX1=DLOG(X1)
        DLOGXM=DLOG(-X)
        DLI3=+TRIS12(1D0/X1)-TRILOG(1D0/X1)
     &       +DLOGX1*DLOGX1*(-DLOGXM/2D0+DLOGX1/3D0)
     &       +(DLOGXM-DLOGX1)*DDILOG(X)-ZET2*DLOGX1
     &       -DLOGX1*DLOGXM*(DLOGX1-DLOGXM)/2D0
      ENDIF
      END

      FUNCTION DS12(X)
* This is S12(x) for x<=1 only!
      IMPLICIT REAL*8(A-Z)
      DATA ZET2/1.6449340668482262D0/
      DATA ZET3/1.2020569031595943D0/
      IF(X.GT.0D0) THEN
        DS12=TRIS12(X)
      ELSE IF(X.EQ.0D0) THEN
        DS12=0D0
      ELSE IF(X.GT.-1D0) THEN
        X1=1D0-X
        DLOGX1=DLOG(X1)
        DS12=TRIS12(-X/X1)+1D0/6*DLOGX1**3
      ELSE IF(X.EQ.-1D0) THEN
        DS12=ZET3/8
      ELSE
        X1=1D0-X
        DLOGX1=DLOG(X1)
        DS12=-TRILOG(1D0/X1)-DLOGX1*DLOGX1*(DLOG(-X)/2D0-DLOGX1/6D0)
     &       -DLOGX1*DDILOG(X)+ZET3-ZET2*DLOGX1
      ENDIF
      END
