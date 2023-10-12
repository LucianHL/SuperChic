      subroutine tpqcdcalc
      implicit none
      double precision sum,rzmin,rzmax,tpqcd
      integer i
      
      include 'pi.f'
      include 'ion.f'
      include 'tpqcdpars.f'
      include 'rho0.f'

      sum=0d0
      
      itpqcd=900

      rzmax=5d0*rzg
      rzmin=0d0
      
      do i=1,itpqcd+1

         rz=rzmin+(rzmax-rzmin)*dble(i-1)/dble(itpqcd)
         
         tpqcdarr(1,i)=rz
         tpqcdarr(2,i)=tpqcd(dabs(rz))
         
      enddo

      return
      end
