      function rhoxyint(ia,r)
      implicit none
      double precision rmax,rhoxyint,rbin,del,r,m
      integer it,ia

      include 'rhoxypars.f'

      rmax=rhoxyarr(1,irho)

      if(r.gt.rmax)then
         rhoxyint=0d0
      else
         rbin=rhoxyarr(1,2)-rhoxyarr(1,1)
         it=nint(r/rbin)
         if(dble(it).gt.(r/rbin))then
            it=it-1
         endif
         m=(rhoxyarr(ia+1,it+2)-rhoxyarr(ia+1,it+1))
         m=m/(rhoxyarr(1,it+2)-rhoxyarr(1,it+1))
         del=r-rhoxyarr(1,it+1)
         rhoxyint=m*del+rhoxyarr(ia+1,it+1)
      endif

      return
      end
