ccc   generate six body decay for particles of mass m1
      subroutine sixbody(m1,wt)
      implicit none
      double precision m1,wt,ein
      integer npart
      integer i,k
      double precision pcm(4),pboo(4),px(4)
      double precision am(100),pout(4,100)

      include 'mom.f'
      include 'wtinit.f'

      npart=6
      ein=dsqrt(q(4,5)**2-q(3,5)**2-q(2,5)**2-q(1,5)**2)
      do i=1,6
         am(i)=m1
      enddo

      call rambo(npart,ein,am,pout,wt)

      do i=1,4
         px(i)=q(i,5)
      enddo

      do k=1,6
         do i=1,4
            pcm(i)=pout(i,k)
         enddo
         call boost(ein,px,pcm,pboo)
         do i=1,4
            q(i,5+k)=pboo(i)
         enddo
      enddo

      wt=wt/wt6i

      return
      end
