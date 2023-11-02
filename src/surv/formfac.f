ccc   Pomeron -- diffrative eignstate i form factor
      subroutine formfac(t1,t2,out)
      implicit none
      double precision t1,t2,wt,delta1,delta2,out
      integer i1,i2

      include 'nchan.f'
      include 'survpars.f'

      out=0d0

      delta1=dexp(-t1/2d0)
      delta2=dexp(-t2/2d0)

      do i1=1,nch
         do i2=1,nch

            wt=gaa(i1)*gaa(i2)*pp0(i1)*pp0(i2)/dble(nch)**2
            wt=wt*dexp(-((t1+0.08d0+bb0(i1))*bex(i1))**cc0(i1)+
     &     (bex(i1)*(bb0(i1)+0.08d0))**cc0(i1))
            wt=wt*dexp(-((t2+0.08d0+bb0(i2))*bex(i2))**cc0(i2)+
     &     (bex(i2)*(bb0(i2)+0.08d0))**cc0(i2))

            out=out+wt

         enddo
      enddo

      return
      end
