ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C     Parameters that should not be recomputed event by event.
C
      IF(READLHA) THEN

        G = 2 * DSQRT(AS*PI)  ! for the first init

        MDL_CONJG__CKM3X3 = 1.000000D+00

        MDL_CONJG__CKM1X1 = 1.000000D+00

        MDL_CKM3X3 = 1.000000D+00

        MDL_COMPLEXI = DCMPLX(0.000000D+00,1.000000D+00)

        MDL_MZ__EXP__2 = MDL_MZ**2

        MDL_MZ__EXP__4 = MDL_MZ**4

        MDL_SQRT__2 = SQRT(DCMPLX(2.000000D+00))

        MDL_MH__EXP__2 = MDL_MH**2

        MDL_AEW = 1.000000D+00/AEWM1

        MDL_MW = DBLE(SQRT(DCMPLX(MDL_MZ__EXP__2/2.000000D+00
     $   +SQRT(DCMPLX(MDL_MZ__EXP__4/4.000000D+00-(MDL_AEW*PI
     $   *MDL_MZ__EXP__2)/(MDL_GF*MDL_SQRT__2))))))

        MDL_SQRT__AEW = DBLE(SQRT(DCMPLX(MDL_AEW)))

        MDL_EE = DBLE(2.000000D+00*MDL_SQRT__AEW*SQRT(DCMPLX(PI)))

        MDL_MW__EXP__2 = MDL_MW**2

        MDL_SW2 = 1.000000D+00-MDL_MW__EXP__2/MDL_MZ__EXP__2

        MDL_CW = DBLE(SQRT(DCMPLX(1.000000D+00-MDL_SW2)))

        MDL_SQRT__SW2 = DBLE(SQRT(DCMPLX(MDL_SW2)))

        MDL_SW = MDL_SQRT__SW2

        MDL_G1 = MDL_EE/MDL_CW

        MDL_GW = MDL_EE/MDL_SW

        MDL_VEV = (2.000000D+00*MDL_MW*MDL_SW)/MDL_EE

        MDL_VEV__EXP__2 = MDL_VEV**2

        MDL_LAM = MDL_MH__EXP__2/(2.000000D+00*MDL_VEV__EXP__2)

        MDL_YB = (MDL_YMB*MDL_SQRT__2)/MDL_VEV

        MDL_YT = (MDL_YMT*MDL_SQRT__2)/MDL_VEV

        MDL_YTAU = (MDL_YMTAU*MDL_SQRT__2)/MDL_VEV

        MDL_MUH =DBLE(SQRT(DCMPLX(MDL_LAM*MDL_VEV__EXP__2)))

        MDL_I1X33 = MDL_YB*MDL_CONJG__CKM3X3

        MDL_I2X33 = MDL_YT*MDL_CONJG__CKM3X3

        MDL_I3X33 = MDL_CKM3X3*MDL_YT

        MDL_I4X33 = MDL_CKM3X3*MDL_YB

        MDL_EE__EXP__2 = MDL_EE**2

        MDL_SW__EXP__2 = MDL_SW**2

        MDL_CW__EXP__2 = MDL_CW**2

        MDL_SQRT__GF = SQRT(DCMPLX(MDL_GF))

        MDL_NB__2__EXP__0_25 = 2.000000D+00**2.500000D-01

        MDL_STH2 = 5.000000D-01-SQRT(DCMPLX(1.000000D+00-(2.000000D+00
     $   *MDL_AEW*PI*MDL_SQRT__2)/(MDL_GF*MDL_MZ__EXP__2)))/2.000000D
     $   +00

        MDL_SQRT__STH2 = SQRT(DCMPLX(MDL_STH2))

        MDL_STH = MDL_SQRT__STH2

        MDL_CTH = SQRT(DCMPLX(1.000000D+00-MDL_STH2))

        MDL_VEVHAT = 1.000000D+00/(MDL_NB__2__EXP__0_25*MDL_SQRT__GF)

        MDL_VEVHAT__EXP__2 = MDL_VEVHAT**2

        MDL_VEVHAT__EXP__3 = MDL_VEVHAT**3

        MDL_CTH__EXP__2 = MDL_CTH**2

        MDL_STH__EXP__2 = MDL_STH**2

        MDL_CTH__EXP__3 = MDL_CTH**3

        MDL_STH__EXP__3 = MDL_STH**3

        MDL_DMZ2 = (2.000000D+00*MDL_CHWB
     $   *MDL_CTH*MDL_STH*MDL_VEVHAT__EXP__2)/MDL_LAMBDASMEFT__EXP__2


        MDL_DGW = -(MDL_CTH__EXP__2*(MDL_DMZ2+(2.000000D+00
     $   *MDL_CHWB*MDL_STH__EXP__3*MDL_VEVHAT__EXP__2)/(MDL_CTH
     $   *MDL_LAMBDASMEFT__EXP__2)))/(2.000000D+00*(1.000000D+00
     $   -2.000000D+00*MDL_STH__EXP__2))

         MDL_MWSM = MDL_CTH*MDL_MZ

        MDL_DMW = MDL_DGW*MDL_MWSM

        MDL_MW_SMEFT = MDL_DMW+MDL_MWSM

      ENDIF
C
C     Parameters that should be recomputed at an event by even basis.
C
      AS = G**2/4/PI

      MDL_SQRT__AS = DBLE(SQRT(DCMPLX(AS)))

      MDL_G__EXP__2 = G**2

C
C     Parameters that should be updated for the loops.
C
C
C     Definition of the EW coupling used in the write out of aqed
C
      GAL(1) = 3.5449077018110318D0 / DSQRT(ABS(AEWM1))
      GAL(2) = 1D0

