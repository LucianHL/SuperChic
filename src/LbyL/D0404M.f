      COMPLEX*16 FUNCTION D0404M(AQ2,AP2,AM12)
**********************************************
* AQ2 = s,u    AP2 = t,u    AM12 = M^2 - i*eps
      IMPLICIT NONE
      REAL*8 AQ2,AP2,TET,ABS,zet2
      COMPLEX*16 AM12,IPI,XSPENZ,SQRT,CMPLX,LOG,C01
      COMPLEX*16 A1,A2,A3,lA3,A2A3,A1A3
      PARAMETER (zet2=1.6449340668482264364724151666460251d0)
*
      IPI=CMPLX(0d0,1d0)*4d0*ATAN(1d0)
      A1=SQRT(1d0-4d0*AM12/AQ2)        
      A2=SQRT(1d0-4d0*AM12/AP2)        
      A3=SQRT(1d0-4d0*REAL(AM12)*(AQ2+AP2)/AQ2/AP2)
      lA3  = 4d0*AM12*(AQ2+AP2)/AQ2/AP2/(1d0+A3)
      A2A3 = 4d0*AM12/AQ2/(A2+A3)
      A1A3 = 4d0*AM12/AP2/(A1+A3)
*
      D0404M =-2d0/AQ2/AP2/A3*(
     &+XSPENZ(+(1d0+A3)/(A1+A3))
     &-XSPENZ(+(lA3)/(A1A3))
*     &+XSPENZ(-(1d0+A3)/(A1A3))
     & -XSPENZ(-A1A3/(1d0+A3))-1d0/2*LOG(A1A3/(1d0+A3))**2-zet2
     &-XSPENZ(-(lA3)/(A1+A3))
     &+XSPENZ(+(1d0+A3)/(A2+A3))
     &-XSPENZ(+(lA3)/(A2A3))
*     &+XSPENZ(-(1d0+A3)/(A2A3))
     & -XSPENZ(-A2A3/(1d0+A3))-1d0/2*LOG(A2A3/(1d0+A3))**2-zet2
     &-XSPENZ(-(lA3)/(A2+A3))
     &+LOG(1d0+1d0/A3)*2d0*IPI*TET(-IMAG(A1+A3))
     &                        )
** WC = 2*t/s
*      WC = -AWC
*      if (TET(AQ2*AP2).eq.0) then
*          if ((AP2/AQ2).ge.WC) then
*      D0404M = 1d0/AQ2/REAL(AM12)*( 
** t=0 answer: 
**             (2d0-A1*LOG((A1+1d0)/(A1-1d0)))
** t small answer:
*     & - AP2/(REAL(AM12))*
*     &         (-37d0/144-7d0/48*A1**2+(5d0/48-3d0/16*A1**2)*LOG(2d0)
*     &          + 1d0/6*A1**3*LOG(-(A1+1d0)/(1d0-A1)) )
*     &  + 2d0 - A1*LOG(-(A1+1d0)/(1d0-A1))
*     &                             )
*          endif
*          if ((AP2/AQ2).le.(-1d0-WC)) then
*      D0404M = -2d0/AQ2*(
*     & C01(0d0,0d0,AQ2,AM12,AM12,AM12)-C01(0d0,0d0,-AQ2,AM12,AM12,AM12))
*          endif
*      endif
**
      if ((TET(AQ2*AP2).eq.1).or.(TET(4d0*REAL(AM12)-AQ2).eq.1)) then
      D0404M = REAL(D0404M)
      endif
*
      RETURN
      END

