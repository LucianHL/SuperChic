ccc   evolution of flavour-singlet wave function
      subroutine wfsinglet(q,a21,a2g,a21n,a2gn)
      implicit double precision(a-y)

      include 'anom.f'
      include 'pi.f'
      include 'ewpars.f'

      if(qs.gt.mb/2d0)then
         nf=5d0
      else
         nf=4d0
      endif

      b0=(33d0-2d0*nf)/(12d0*pi)

      qs=q**2
      mu0=1d0

      alphasg0=alphas(mu0**2)
      alphasg=alphas(qs)

      xcap2=(alphasg0/alphasg)**(gamp(2)/b0/4d0/pi)
      ycap2=(alphasg0/alphasg)**(gamm(2)/b0/4d0/pi)
      
      bm2=(a21-a2g/rhop(2))/(rhom(2)-1d0/rhop(2))
      bp2=(a21/rhom(2)-a2g)/(1d0/rhom(2)-rhop(2))
      
      a21n=bp2*xcap2+rhom(2)*bm2*ycap2
      a2gn=rhop(2)*bp2*xcap2+bm2*ycap2 

      return
      end
