ccc   generate three body decay for particles of mass m1,m2,m3
      subroutine threebody (j,in,i1,i2,i3,m1,m2,m3,wt)
      implicit none
      double precision m1,m2,m3,ein,wt
      integer in,i1,i2,i3,npart
      integer i,j
      double precision pcm(4),pboo(4),px(4)
      double precision am(100),pout(4,100)

      include 'mom.f'
      include 'wtinit.f'

      npart=3
      ein=dsqrt(q(4,in)**2-q(3,in)**2-q(2,in)**2-q(1,in)**2)
      am(1)=m1
      am(2)=m2
      am(3)=m3

      call rambo(npart,ein,am,pout,wt)

      do i=1,4
         px(i)=q(i,in)
      enddo

      do i=1,4
         pcm(i)=pout(i,1)
      enddo
      call boost(ein,px,pcm,pboo)
      do i=1,4
         q(i,i1)=pboo(i)
      enddo

      do i=1,4
         pcm(i)=pout(i,2)
      enddo
      call boost(ein,px,pcm,pboo)
      do i=1,4
         q(i,i2)=pboo(i)
      enddo

      do i=1,4
         pcm(i)=pout(i,3)
      enddo
      call boost(ein,px,pcm,pboo)
      do i=1,4
         q(i,i3)=pboo(i)
      enddo

      wt=wt/wt3i(j)

      return
      end
