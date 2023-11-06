      subroutine rhoxycalc
      implicit none
      double precision r,rzmax,rzmin,sum
      double precision rhonxyc,rhozxyc,rhozxy,rhonxy
      integer i

      include 'pi.f'
      include 'ion.f'
      include 'rhoxypars.f'
      include 'rho0.f'

ccccccccccc

      sum=0d0

      irho=900

      rzmax=2d0*rzg
      rzmin=0d0

      do i=1,irho+1

         r=rzmin+(rzmax-rzmin)*dble(i-1)/dble(irho)

         rhoxyarr(1,i)=r
         rhoxyarr(2,i)=rhozxy(dabs(r))
         rhoxyarr(3,i)=rhonxy(dabs(r))
         rhoxyarr(4,i)=rhozxyc(dabs(r))
         rhoxyarr(5,i)=rhonxyc(dabs(r))

      enddo

      return
      end
