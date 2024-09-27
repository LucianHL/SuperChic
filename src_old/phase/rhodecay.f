ccc   generates rho meson invariant mass distribution according to modified BW
      subroutine bwrho(mout,jrho)
      implicit none
      double precision rm1,ran2,mwidth,mminr,mmaxr,almin,almax,al1
      double precision mout,jrho

      include 'mres.f'
      include 'mpip.f'
      include 'widths.f'
      include 'pi.f'

      mminr=2d0*mpip
      mmaxr=mres+4d0*width

      almin=datan2(-(mres**2-mminr**2),width*mres)
      almax=datan2(-(mres**2-mmaxr**2),width*mres)

      rm1=ran2()
      al1=almin+(almax-almin)*rm1
      mout=dsqrt(dtan(al1)*mres*width+mres**2)

      jrho=almax-almin
      jrho=jrho*width*mres
      jrho=jrho*(1d0+dtan(al1)**2)
      jrho=jrho/2d0/mout

cccccc up to this point just change of variables
ccccc
ccccc  now BW

      mwidth=width*((1d0-4d0*mpip**2/mout**2)/
     &        (1d0-4d0*mpip**2/mres**2))**1.5d0

      jrho=jrho*2d0*mout**2*mwidth/pi
      jrho=jrho/((mres**2-mout**2)**2+mout**2*
     &     mwidth**2)

      return
      end
