      function tpqcdint(r)
      implicit none
      double precision rmax,m,rbin,del,tpqcdint,r
      integer it
      
      include 'tpqcdpars.f'
      
      rmax=tpqcdarr(1,itpqcd)

      if(r.gt.rmax)then
         tpqcdint=0d0
      else
         rbin=tpqcdarr(1,2)-tpqcdarr(1,1)
         it=nint(r/rbin)
         if(dble(it).gt.(r/rbin))then
            it=it-1
         endif
         m=(tpqcdarr(2,it+2)-tpqcdarr(2,it+1))
         m=m/(tpqcdarr(1,it+2)-tpqcdarr(1,it+1))
         del=r-tpqcdarr(1,it+1)
         tpqcdint=m*del+tpqcdarr(2,it+1)
      endif
         
      return
      end
