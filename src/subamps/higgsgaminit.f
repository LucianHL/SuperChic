ccc   gamgam --> SM higgs subprocess amplitude : normalisation
      subroutine higgsgaminit(mx)
      implicit none
      double precision xt,xw,xtau,xc,xb,mx
      complex*16 zt,zb,zc,zw,ztau

      include 'pi.f'
      include 'ewpars.f'
      include 'normh.f'
      include 'quarkonia.f'

      xc=mc**2/mx**2
      xb=mb**2/mx**2
      xt=mt**2/mx**2
      xtau=mtau**2/mx**2
      xw=mw**2/mx**2

      call iq(xtau,ztau)
      call iq(xc,zc)
      call iq(xb,zb)
      call iq(xt,zt)
      call iw(xw,zw)

      normh=cdabs(3d0*(4d0*(zc+zt)+zb)/9d0+ztau+zw)
      normh=normh*alpha*mx**2/4d0/pi/v

c      print*,normh**2*2d0/(16d0*pi*mx)/2d0
      
      return
      end

      subroutine iq(x,zout)
      implicit none
      double precision x
      complex*16 zout,fx

      include 'zi.f'
      include 'pi.f'

      if(x.lt.0.25d0)then
         fx=(dlog((1d0+dsqrt(1d0-4d0*x))/(1d0-dsqrt(1d0-4d0*x)))
     &        -zi*pi)**2/2d0
      else
         fx=-2d0*dasin(1d0/2d0/dsqrt(x))**2
      endif

      zout=4d0*x*(2d0+(4d0*x-1d0)*fx)

      return
      end


      subroutine iw(x,zout)
      implicit none
      double precision x
      complex*16 zout,fx

      include 'pi.f'
      include 'zi.f'

      if(x.lt.0.25d0)then
         fx=(dlog((1d0+dsqrt(1d0-4d0*x))/(1d0-dsqrt(1d0-4d0*x)))
     &        -zi*pi)**2/2d0
      else
         fx=-2d0*dasin(1d0/2d0/dsqrt(x))**2
      endif

      zout=-2d0*(6d0*x+1d0+6d0*x*(2d0*x-1d0)*fx)

      return
      end







