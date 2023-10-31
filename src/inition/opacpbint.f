      function opacpbint(r)
      implicit none
      double precision r,del,m,rbin,rmin,rmax
      double precision opacpbint
      integer it

      include 'opacpbpars.f'

      rmax=opacpbarr(1,ioppb)
      rmin=opacpbarr(1,1)

      if(r.gt.rmax)then
         opacpbint=opacpbarr(2,ioppb)
      elseif(r.lt.rmin)then
         opacpbint=opacpbarr(2,1)
      else
         rbin=opacpbarr(1,2)-opacpbarr(1,1)
         it=nint((r-rmin)/rbin)
         if(dble(it).gt.((r-rmin)/rbin))then
            it=it-1
         endif
         m=(opacpbarr(2,it+2)-opacpbarr(2,it+1))
         m=m/(opacpbarr(1,it+2)-opacpbarr(1,it+1))
         del=r-opacpbarr(1,it+1)
         opacpbint=m*del+opacpbarr(2,it+1)
      endif

      return
      end

      function opacpbint_3(r)
      implicit none
      double precision r,del,m,rbin,rmin,rmax
      double precision opacpbint_3
      integer it

      include 'opacpbpars.f'

      rmax=opacpbarr(1,ioppb)
      rmin=opacpbarr(1,1)

      if(r.gt.rmax)then
         opacpbint_3=opacpbarr(3,ioppb)
      elseif(r.lt.rmin)then
         opacpbint_3=opacpbarr(3,1)
      else
         rbin=opacpbarr(1,2)-opacpbarr(1,1)
         it=nint((r-rmin)/rbin)
         if(dble(it).gt.((r-rmin)/rbin))then
            it=it-1
         endif
         m=(opacpbarr(3,it+2)-opacpbarr(3,it+1))
         m=m/(opacpbarr(1,it+2)-opacpbarr(1,it+1))
         del=r-opacpbarr(1,it+1)
         opacpbint_3=m*del+opacpbarr(3,it+1)
      endif


      return
      end
