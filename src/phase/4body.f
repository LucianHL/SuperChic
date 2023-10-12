ccc   generate four body decay for particles of pair mass m1,m2
      subroutine fourbody(m1,m2,wt)
      implicit none
      double precision wt,m1,m2,ein
      integer npart
      integer i,k
      double precision pcm(4),pboo(4),px(4)
      double precision am(100),pout(4,100)
  
      include 'mom.f'
      include 'wtinit.f'

      npart=4
      ein=dsqrt(q(4,5)**2-q(3,5)**2-q(2,5)**2-q(1,5)**2)
      am(1)=m1
      am(2)=m1
      am(3)=m2
      am(4)=m2

      call rambo(npart,ein,am,pout,wt)

      do i=1,4
         px(i)=q(i,5)
      enddo 

      do k=1,4
         do i=1,4
            pcm(i)=pout(i,k)
         enddo 
         call boost(ein,px,pcm,pboo)
         do i=1,4
            q(i,5+k)=pboo(i)
         enddo
      enddo
      
      wt=wt/wt4i

      return
      end
