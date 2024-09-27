ccc   calculates normalization factor when using RAMBO
ccc   for four-body phase space
      subroutine sixbodyinit(ein,m1)
      implicit none
      double precision ein,m1,sum,wt
      integer nrun,npart,n,i
      double precision am(100),pout(4,100)

      include 'wtinit.f'

      nrun=100000

      npart=6

      do i=1,6
         am(i)=m1
      enddo

      sum=0d0

      do n=1,nrun
         call rambo(npart,ein,am,pout,wt)
         sum=sum+wt
      enddo

      wt6i=sum/dble(nrun)

      return
      end
