ccc   calculates normalization factor when using RAMBO
ccc   for two-body phase space
      subroutine twobodyinit(i,ein,m1,m2)
      implicit none
      double precision wt,sum,m1,m2,ein
      integer nrun,npart,n,i
      double precision am(100),pout(4,100)

      include 'wtinit.f'

      nrun=100000

      npart=2
      am(1)=m1
      am(2)=m2

      sum=0d0

      do n=1,nrun
         call rambo(npart,ein,am,pout,wt)
         sum=sum+wt
      enddo

      wt2i(i)=sum/dfloat(nrun)

      return
      end
