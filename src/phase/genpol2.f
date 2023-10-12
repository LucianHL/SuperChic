ccc   generates polarization vectors for chi_2
      subroutine genpol2
      implicit none
      double precision pnorm
      double precision n1(4),n2(4)
      double precision nchi(4)
      integer i,k,l

      include 'polvecs.f'
      include 'zi.f'
      include 'vars.f'
      include 'mom.f'

      pnorm=dsqrt(q(1,5)**2+q(2,5)**2+q(3,5)**2)
      do i=1,3
         nchi(i)=q(i,5)/pnorm
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
         echi1(1,i)=-(n1(i)+zi*n2(i))/dsqrt(2d0)
         echi1(2,i)=(n1(i)-zi*n2(i))/dsqrt(2d0)
      enddo

      echi1(3,4)=0d0
      do i=1,3
         echi1(3,i)=nchi(i)
      enddo

      do l=1,4
         do k=1,4
            echi2(1,l,k)=echi1(1,l)*echi1(1,k)
            echi2(2,l,k)=(echi1(1,l)*echi1(3,k)+echi1(3,l)*echi1(1,k))
     &           /dsqrt(2d0)
            echi2(3,l,k)=(echi1(1,l)*echi1(2,k)+2d0*echi1(3,l)*
     &           echi1(3,k)+echi1(2,l)*echi1(1,k))/dsqrt(6d0)
            echi2(4,l,k)=(echi1(2,l)*echi1(3,k)+echi1(3,l)*echi1(2,k))
     &           /dsqrt(2d0)
            echi2(5,l,k)=echi1(2,l)*echi1(2,k)
         enddo
      enddo

            

      return
      end
