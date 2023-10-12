      function opacpbpint(r)
      implicit double precision(a-y)
      integer it
      
      include 'opacpbppars.f'
      
      rmax=opacpbparr(1,ioppbp)
      rmin=opacpbparr(1,1)

      
      if(r.gt.rmax)then
         opacpbpint=1d0
      elseif(r.lt.rmin)then
         opacpbpint=0d0
      else
         rbin=opacpbparr(1,2)-opacpbparr(1,1)
         it=nint((r-rmin)/rbin)
         if(dble(it).gt.((r-rmin)/rbin))then
            it=it-1
         endif
         m=(opacpbparr(2,it+2)-opacpbparr(2,it+1))
         m=m/(opacpbparr(1,it+2)-opacpbparr(1,it+1))
         del=r-opacpbparr(1,it+1)
         opacpbpint=m*del+opacpbparr(2,it+1)
      endif
         
      return
      end
