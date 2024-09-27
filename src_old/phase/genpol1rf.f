ccc   generates polarization vectors for psi (in r.f.)
      subroutine genpol1rf(p,echi1)
      implicit none
      double precision m,pnorm
      double precision n1(4),n2(4)
      double precision nchi(4),p(4)
      integer i
      complex*16 echi1(3,4)

      include 'zi.f'
      include 'mom.f'

      m=dsqrt(p(4)**2-p(3)**2-p(2)**2-p(1)**2)

      pnorm=dsqrt(p(1)**2+p(2)**2+p(3)**2)
      do i=1,3
         nchi(i)=p(i)/pnorm
      enddo

      n1(1)=nchi(2)/dsqrt(nchi(1)**2+nchi(2)**2)
      n1(2)=-nchi(1)/dsqrt(nchi(1)**2+nchi(2)**2)
      n1(3)=0d0

      n2(1)=n1(2)*nchi(3)-n1(3)*nchi(2)
      n2(2)=n1(3)*nchi(1)-n1(1)*nchi(3)
      n2(3)=n1(1)*nchi(2)-n1(2)*nchi(1)

      echi1(1,4)=0d0
      echi1(2,4)=0d0

      do i=1,3
         echi1(1,i)=(n1(i)+zi*n2(i))/dsqrt(2d0)
         echi1(2,i)=-(n1(i)-zi*n2(i))/dsqrt(2d0)
      enddo

      echi1(3,4)=0d0
      do i=1,3
         echi1(3,i)=nchi(i)
      enddo

      return
      end
