ccc   generates three-body phase space
      subroutine threejetps(mx,mq,ps)
      implicit none
      double precision mx,wt,mq,ps,ein
      double precision px(4),pcm(4)
      double precision pboo(4)
      integer i,j,npart
      double precision am(100),pout(4,100)

      include 'partonmom3.f'
      include 'pi.f'
      include 'mom.f'

      pa(4)=mx/2d0
      pa(3)=mx/2d0
      pb(4)=mx/2d0
      pb(3)=-mx/2d0
      do i=1,2
         pa(i)=0d0
         pb(i)=0d0
      enddo

      npart=3
      ein=mx

      do j=1,npart
         am(j)=mq
      enddo

      call rambo(npart,ein,am,pout,wt)

      do i=1,4
         px(i)=q(i,5)
         pcm(i)=pout(i,1)
      enddo 

      call boost(mx,px,pcm,pboo)

      do i=1,4
         q(i,6)=pboo(i)
         pcm(i)=pout(i,2)
      enddo 

      call boost(mx,px,pcm,pboo)
      
      do i=1,4
         q(i,7)=pboo(i)
      enddo 

      do i=1,4
         q(i,8)=q(i,5)-q(i,7)-q(i,6)
         p3(i)=pout(i,1)
         p4(i)=pout(i,2)
         p5(i)=pout(i,3)
      enddo

      ps=wt/pi**3

      return
      end
