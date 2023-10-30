c      SUBROUTINE VEGAS(FXN,AVGI,SD,CHI2A)
      SUBROUTINE VEGAS(AVGI,SD,CHI2A) 
C
C   SUBROUTINE PERFORMS N-DIMENSIONAL MONTE CARLO INTEG'N
C      - BY G.P. LEPAGE   SEPT 1976/(REV)APR 1978
C
      IMPLICIT double precision(A-H,O-Z)
      implicit integer(i-n)
      COMMON/BVEG1/XL(10),XU(10),ACC,NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG2/XI(50,10),SI,SI2,SWGT,SCHI,NDO,IT
      DIMENSION D(50,10),DI(50,10),XIN(50),R(50),DX(10),DT(10)
     1   ,KG(10),IA(10)
      double precision QRAN(10),x(10)
      DATA NDMX/50/,ALPH/1.5D0/,ONE/1D0/,MDS/1/
      
      include 'wmax.f'
      include 'genunw.f'
C
      NDO=1  
      DO 1 J=1,NDIM
1     XI(1,J)=ONE
C   

      ENTRY VEGAS1(AVGI,SD,CHI2A)   
c      ENTRY VEGAS1(FXN,AVGI,SD,CHI2A)
C         - INITIALIZES CUMMULATIVE VARIABLES, BUT NOT GRID
      IT=0
      SI=0.
      SI2=SI
      SWGT=SI
      SCHI=SI
C
      ENTRY VEGAS2(AVGI,SD,CHI2A)   
c      ENTRY VEGAS2(FXN,AVGI,SD,CHI2A)
C         - NO INITIALIZATION
      ND=NDMX
      NG=1
      IF(MDS.EQ.0) GO TO 2
      NG=(NCALL/2.)**(1./NDIM)
      MDS=1
      IF((2*NG-NDMX).LT.0) GO TO 2
      MDS=-1
      NPG=NG/NDMX+1
      ND=NG/NPG
      NG=NPG*ND
2     K=NG**NDIM
      NPG=NCALL/K
      IF(NPG.LT.2) NPG=2
      CALLS=NPG*K
      DXG=ONE/NG
      DV2G=(CALLS*DXG**NDIM)**2/NPG/NPG/(NPG-ONE)
      XND=ND
      NDM=ND-1
      DXG=DXG*XND
      XJAC=ONE/CALLS
      DO 3 J=1,NDIM
      DX(J)=XU(J)-XL(J)
3     XJAC=XJAC*DX(J)
C
C   REBIN PRESERVING BIN DENSITY
C
      IF(ND.EQ.NDO) GO TO 8
      RC=NDO/XND
      DO 7 J=1,NDIM
      K=0
      XN=0.
      DR=XN
      I=K
4     K=K+1
      DR=DR+ONE
      XO=XN
      XN=XI(K,J)
5     IF(RC.GT.DR) GO TO 4
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR
      IF(I.LT.NDM) GO TO 5
      DO 6 I=1,NDM
6     XI(I,J)=XIN(I)
7     XI(ND,J)=ONE
      NDO=ND
C
c8     IF(NPRN.NE.0) WRITE(6,200) NDIM,CALLS,IT,ITMX,ACC,MDS,ND
c     1                           ,(XL(J),XU(J),J=1,NDIM)
8     IF(NPRN.NE.0)THEN
c         WRITE(6,300) NDIM,CALLS,IT
      ENDIF

C
      ENTRY VEGAS3(AVGI,SD,CHI2A)
c      ENTRY VEGAS3(FXN,AVGI,SD,CHI2A)
C         - MAIN INTEGRATION LOOP
9     IT=IT+1
      TI=0.
      TSI=TI
      DO 10 J=1,NDIM
      KG(J)=1
      DO 10 I=1,ND
      D(I,J)=TI
10    DI(I,J)=TI
C
11    FB=0.
      F2B=FB
      K=0
12    K=K+1
      CALL ARAN9(QRAN,NDIM)
      WGT=XJAC
      DO 15 J=1,NDIM
      XN=(KG(J)-QRAN(J))*DXG+ONE
      IA(J)=XN
      IF(IA(J).GT.1) GO TO 13
      XO=XI(IA(J),J)
      RC=(XN-IA(J))*XO
      GO TO 14
13    XO=XI(IA(J),J)-XI(IA(J)-1,J)
      RC=XI(IA(J)-1,J)+(XN-IA(J))*XO
14    X(J)=XL(J)+RC*DX(J)
15    WGT=WGT*XO*XND
C
      F=WGT
      F=F*cs(X,WGT)
c      F=F*FXN(X,WGT)

      F2=F*F
      FB=FB+F
      F2B=F2B+F2
      DO 16 J=1,NDIM
      DI(IA(J),J)=DI(IA(J),J)+F
16    IF(MDS.GE.0) D(IA(J),J)=D(IA(J),J)+F2
      IF(K.LT.NPG) GO TO 12
C
      F2B=DSQRT(F2B*NPG)
      F2B=(F2B-FB)*(F2B+FB)
      TI=TI+FB
      TSI=TSI+F2B
      IF(MDS.GE.0) GO TO 18
      DO 17 J=1,NDIM
17    D(IA(J),J)=D(IA(J),J)+F2B
18    K=NDIM
19    KG(K)=MOD(KG(K),NG)+1
      

      IF(KG(K).NE.1) GO TO 11
      K=K-1
      IF(K.GT.0) GO TO 19
C
C   FINAL RESULTS FOR THIS ITERATION
C
      TSI=TSI*DV2G
      TI2=TI*TI
      WGT=TI2/TSI
      SI=SI+TI*WGT
      SI2=SI2+TI2
      SWGT=SWGT+WGT
      SCHI=SCHI+TI2*WGT
      AVGI=SI/SWGT
      SD=SWGT*IT/SI2
      CHI2A=SD*(SCHI/SWGT-AVGI*AVGI)/(IT-.999)
      SD=DSQRT(ONE/SD)
C
      IF(NPRN.EQ.0) GO TO 21
      TSI=DSQRT(TSI)
c      WRITE(6,201) IT,TI,TSI,AVGI,SD,CHI2A
      if(genunw)then
         WRITE(6,303) INT(CALLS),TI,TSI,AVGI,SD,SD/AVGI*100d0,wmax,
     &        CHI2A
      elseif(calcmax)then
         WRITE(6,301) IT,INT(CALLS),TI,TSI,AVGI,SD,SD/AVGI*100d0,wmax,
     &        CHI2A
      else
         WRITE(6,302) IT,INT(CALLS),TI,TSI,AVGI,SD,SD/AVGI*100d0,CHI2A
      endif   

      IF(NPRN.GE.0) GO TO 21
      DO 20 J=1,NDIM
20    WRITE(6,202) J,(XI(I,J),DI(I,J),D(I,J),I=1,ND)
C
C   REFINE GRID
C
21    DO 23 J=1,NDIM
      XO=D(1,J)
      XN=D(2,J)
      D(1,J)=(XO+XN)/2.
      DT(J)=D(1,J)
      DO 22 I=2,NDM
      D(I,J)=XO+XN
      XO=XN
      XN=D(I+1,J)
      D(I,J)=(D(I,J)+XN)/3.
22    DT(J)=DT(J)+D(I,J)
      D(ND,J)=(XN+XO)/2.
23    DT(J)=DT(J)+D(ND,J)
C
      DO 28 J=1,NDIM
      RC=0.
      DO 24 I=1,ND
      R(I)=0.
      IF(D(I,J).LE.0.) GO TO 24
      XO=DT(J)/D(I,J)
      R(I)=((XO-ONE)/XO/DLOG(XO))**ALPH
24    RC=RC+R(I)
      RC=RC/XND
      K=0
      XN=0.
      DR=XN
      I=K
25    K=K+1
      DR=DR+R(K)
      XO=XN
      XN=XI(K,J)
26    IF(RC.GT.DR) GO TO 25
      I=I+1
      DR=DR-RC
      XIN(I)=XN-(XN-XO)*DR/R(K)
      IF(I.LT.NDM) GO TO 26
      DO 27 I=1,NDM
27    XI(I,J)=XIN(I)
28    XI(ND,J)=ONE
C
      IF(IT.LT.ITMX.AND.ACC*DABS(AVGI).LT.SD) GO TO 9
C209   FORMAT(' INPUT PARAMETERS FOR VEGAS:  NDIM=',I3,'  NCALL=',F9.0
C     1    /28X,'  IT=',I5,'  ITMX=',I5/28X,'  ACC=',G9.3
C     2    /28X,'  MDS=',I3,'   ND=',I4/28X,'  (XL,XU)=',
C     3    (T40,'( ',G12.6,' , ',G12.6,' )'))
C201   FORMAT(///' INTEGRATION BY VEGAS' / ' ITERATION NO.',I3,
C     1    ':   INTEGRAL =',G14.8/21X,'STD DEV  = ',G10.4 /
C     2    ' ACCUMULATED RESULTS:   INTEGRAL =',G14.8 /
C     3    24X,'STD DEV  = ',G10.4 / 24X,'CHI**2 PER IT''N =',G10.4)
202   FORMAT(' DATA FOR AXIS',I2 / ' ',6X,'X',7X,'  DELT I  ',
     1    2X,' CONV''CE  ',11X,'X',7X,'  DELT I  ',2X,' CONV''CE  '
     2   ,11X,'X',7X,'  DELT I  ',2X,' CONV''CE  ' /
     2    (' ',3G12.4,5X,3G12.4,5X,3G12.4))

C 300  FORMAT(' INPUT PARAMETERS FOR VEGAS:  NDIM=',I3,'  NCALL=',F9.0
C     1    /28X,'  IT=',I5)
 301  FORMAT(//' ************** Integration by Vegas (iteration',I4,') '
     .      ,'*************',/,' ******************** ',1X,I7,2X
     .     ,' calls','   *************************',/,' *',63X,'*'
     .     ,/,' *',13X,'integral = ',G12.6,' +/-  ',G12.6,9X,'*',/,' *'
     .     ,6X,'accum. integral = ',G12.6,' +/-  ',G12.6,9X,'*',/,
     .     ' *',12X,'precision = ',G12.6,' %',25X,'*',/,
     .     ' *',11X,'max weight = ',G12.6,27X,'*',/,
     .     ' *',63X,'*',/,' ******************* chi**2/iteration ='
     .     ,G10.4,'*****************')
 303  FORMAT(//' ******************* Integration by Vegas *************'
     .     ,'***********',/,' ******************** ',1X,I7,2X
     .     ,' calls','   *************************',/,' *',63X,'*'
     .     ,/,' *',13X,'integral = ',G12.6,' +/-  ',G12.6,9X,'*',/,' *'
     .     ,6X,'accum. integral = ',G12.6,' +/-  ',G12.6,9X,'*',/,
     .     ' *',12X,'precision = ',G12.6,' %',25X,'*',/,
     .     ' *',11X,'max weight = ',G12.6,27X,'*',/,
     .     ' *',63X,'*',/,' ******************* chi**2/iteration ='
     .     ,G10.4,'*****************')
 302  FORMAT(//' ************** Integration by Vegas (iteration',I4,') '
     .      ,'*************',/,' ******************** ',1X,I7,2X
     .     ,' calls','   *************************',/,' *',63X,'*'
     .     ,/,' *',13X,'integral = ',G12.6,' +/-  ',G12.6,9X,'*',/,' *'
     .     ,6X,'accum. integral = ',G12.6,' +/-  ',G12.6,9X,'*',/,
     .     ' *',12X,'precision = ',G12.6,' %',25X,'*',/,
     .     ' *',63X,'*',/,' ******************* chi**2/iteration ='
     .     ,G10.4,'*****************')


      RETURN
      END
      SUBROUTINE SAVE(NDIM)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/BVEG2/XI(50,10),SI,SI2,SWGT,SCHI,NDO,IT
C
C   STORES VEGAS DATA (UNIT 7) FOR LATER RE-INITIALIZATION
C
      WRITE(7,200) NDO,IT,SI,SI2,SWGT,SCHI,
     1             ((XI(I,J),I=1,NDO),J=1,NDIM)
      RETURN
      ENTRY RESTR(NDIM)
C
C   ENTERS INITIALIZATION DATA FOR VEGAS
C
      READ(7,200) NDO,IT,SI,SI2,SWGT,SCHI,
     1            ((XI(I,J),I=1,NDO),J=1,NDIM)
200   FORMAT(2I8,4Z16/(5Z16))
      RETURN
      END
C
      SUBROUTINE ARAN9(QRAN,NDIM)
      REAL*8 QRAN(10)
      double precision ran2
      DO I=1,NDIM
         qran(i)=ran2()
      enddo
      RETURN
      END
C
