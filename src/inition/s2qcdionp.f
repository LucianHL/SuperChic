      subroutine s2qcdionp
      implicit none
      double precision hbt,btmax,wt,wt1,sum,sum1,out,bt1
      double precision tpqcdint,opacpbpint
      integer ibt1
      integer nphi,nbt

      include 'pi.f'
      include 'ion.f'
      include 'rho0.f'
      include 's2qcd.f'

      btmax=5d0*rzg
      nphi=10
      nbt=20000

      hbt=btmax/dble(nbt)

      sum=0d0
      sum1=0d0

      do 800 ibt1=1,nbt

         bt1=(dble(ibt1)-0.5d0)*hbt

         wt=opacpbpint(bt1)**2

         wt=wt*bt1*hbt*2d0*pi
         wt=wt*tpqcdint(bt1)

         wt1=bt1*hbt*2d0*pi
         wt1=wt1*tpqcdint(bt1)

         sum=sum+wt
         sum1=sum1+wt1

 800  enddo

      out=sum/sum1
      s2qcd=out

      return
      end
