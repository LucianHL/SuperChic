      function tpint(ia,ri)
      implicit none
      double precision tpint,rmin,rmax,rbin,r,ri,m,del
      integer it,ia
      
      include 'tppars.f'

      rmin=tparr(1,1)
      rmax=tparr(1,itp)

      if(ri.lt.dexp(rmin))then
         if(ia.eq.1)tpint=tparr(2,1)
         if(ia.eq.2)tpint=tparr(3,1)
         return
      endif
      
      r=dlog(ri)
      r=r-rmin
      
      if(ri.gt.dexp(rmax))then
         tpint=0d0
      else
         rbin=tparr(1,2)-tparr(1,1)
         it=nint(r/rbin)
         if(dble(it).gt.(r/rbin))then
            it=it-1
         endif
         m=(tparr(ia+1,it+2)-tparr(ia+1,it+1))
         m=m/(tparr(1,it+2)-tparr(1,it+1))
         del=r-tparr(1,it+1)+rmin
         tpint=m*del+tparr(ia+1,it+1)
      endif
         
      return
      end
