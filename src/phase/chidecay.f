ccc   generates rho meson invariant mass distribution according to modified BW 
      subroutine bwchi(mout,jrho)
      implicit none
      double precision rm1,ran2,norm,mminr,mmaxr,almin,almax,al1,mout,jrho

      include 'mres.f'
      include 'mpip.f'
      include 'widths.f'
      include 'pi.f'

      mminr=mres-10d0*width
      mmaxr=mres+10d0*width
      
      almin=datan(-(mres**2-mminr**2)/width/mres)
      almax=datan(-(mres**2-mmaxr**2)/width/mres)

      norm=almax-almin
      
      rm1=ran2()
      al1=almin+(almax-almin)*rm1
      mout=dsqrt(dtan(al1)*mres*width+mres**2)

      jrho=almax-almin
      jrho=jrho*width*mres
      jrho=jrho*(1d0+dtan(al1)**2)
 
cccccc up to this point just change of variables
ccccc
ccccc  now BW

      jrho=jrho*mout*width/norm
      jrho=jrho/((mres**2-mout**2)**2+mout**2*
     &     width**2)
      
      return
      end
