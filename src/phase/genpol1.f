ccc   generates polarization vectors for chi_1
      subroutine genpol1(in,echi1)
      implicit none
      double precision m,pnorm
      double precision n1(4),n2(4)
      double precision nchi(4)
      integer i,in,irf
      complex*16 echi1(3,4)
      double precision pb(4),pcm(4)
      double precision pboo(4),n1b(4),n2b(4)

      include 'zi.f'
      include 'mom.f'
      include 'polwrf.f'

      if(in.eq.6.or.in.eq.16)then
         irf=1
      elseif(in.eq.7.or.in.eq.17)then
         irf=2
      else
         irf=1
      endif

      m=dsqrt(q(4,in)**2-q(3,in)**2-q(2,in)**2-q(1,in)**2)

      pnorm=dsqrt(q(1,in)**2+q(2,in)**2+q(3,in)**2)
      do i=1,3
         nchi(i)=q(i,in)/pnorm
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

      echi1(3,4)=dsqrt(q(4,in)**2-m**2)/m
      do i=1,3
         echi1(3,i)=q(4,in)*nchi(i)/m
      enddo

ccccccccccccccccccc in W RF

      do i=1,4
         pcm(i)=dreal(echi1(3,i))
         pb(i)=-q(i,in)
      enddo
      pb(4)=q(4,in)

      call boost(m,pb,pcm,pboo)

      do i=1,4
         echirf(irf,3,i)=pboo(i)
      enddo

      do i=1,4
         pcm(i)=n1(i)
         pb(i)=-q(i,in)
      enddo
      pb(4)=q(4,in)

      call boost(m,pb,pcm,n1b)

      do i=1,4
         pcm(i)=n2(i)
         pb(i)=-q(i,in)
      enddo
      pb(4)=q(4,in)

      call boost(m,pb,pcm,n2b)

      echirf(irf,1,4)=0d0
      echirf(irf,2,4)=0d0

      do i=1,3
         echirf(irf,1,i)=(n1b(i)+zi*n2b(i))/dsqrt(2d0)
         echirf(irf,2,i)=-(n1b(i)-zi*n2b(i))/dsqrt(2d0)
      enddo

      return
      end
