ccc   evolution of flavour-octet meson wave function
      subroutine wfoctet(q,n,a28,a28n)
      implicit double precision(a-y)
      integer n

      include 'anom.f'
      include 'pi.f'
      include 'ewpars.f'

      qs=q/2d0

      if(qs.gt.mb/2d0)then
         nf=5d0
      else
         nf=4d0
      endif

      b0=(33d0-2d0*nf)/(12d0*pi)
      qs=q**2
      mu0=1d0
      a28n=a28*(alphas(mu0**2)/alphas(qs*2d0))**(gamqq(n)/b0/4d0/pi) 

      return
      end
