      function opacpint(r)
      implicit none
      double precision opacpint,rmax,rbin,m,del,r
      integer it

      include 'opacppars.f'

      rmax=opacparr(1,iopp)

      if(r.gt.rmax)then
         opacpint=0d0
      else
         rbin=opacparr(1,2)-opacparr(1,1)
         it=nint(r/rbin)
         if(dble(it).gt.(r/rbin))then
            it=it-1
         endif
         m=(opacparr(2,it+2)-opacparr(2,it+1))
         m=m/(opacparr(1,it+2)-opacparr(1,it+1))
         del=r-opacparr(1,it+1)
         opacpint=m*del+opacparr(2,it+1)
      endif

      return
      end
