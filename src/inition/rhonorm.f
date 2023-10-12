      subroutine rhonorm
      implicit none
      double precision wtz,wtn,sumn,sumz
      double precision rmax,r,hr,rhoz,rhon
      integer i,itot

      include 'rho0.f'
      include 'ion.f'
      include 'pi.f'

      rho0z=1d0
      rho0n=1d0
      
      rmax=10d0*rzg
      itot=10000
      hr=rmax/dble(itot)

      sumz=0d0
      sumn=0d0

      do i=1,itot

         r=(dble(i)-0.5d0)*hr

         wtz=rhoz(r)
         wtz=wtz*4d0*pi
         wtz=wtz*r**2*hr

         wtn=rhon(r)
         wtn=wtn*4d0*pi
         wtn=wtn*r**2*hr
         
         sumz=sumz+wtz
         sumn=sumn+wtn
         
      enddo

      rho0z=az/sumz
      rho0n=(an-az)/sumn

      return
      end
