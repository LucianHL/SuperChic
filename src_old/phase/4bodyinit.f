ccc   calculates normalization factor when using RAMBO
ccc   for four-body phase space
      subroutine fourbodyinit(ein,m1,m2)
      implicit none
      double precision ein,m1,m2,sum,wt
      integer nrun,npart,n
      double precision am(100),pout(4,100)

      include 'wtinit.f'

      nrun=100000

      npart=4
      am(1)=m1
      am(2)=m1
      am(3)=m2
      am(4)=m2

      sum=0d0

      do n=1,nrun
         call rambo(npart,ein,am,pout,wt)
         sum=sum+wt
      enddo

      wt4i=sum/dble(nrun)

      return
      end
