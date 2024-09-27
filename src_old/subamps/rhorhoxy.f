ccc   gg --> rhorho subprocess amplitude (x,y dependent)
      subroutine rhorhoxy(p,x,y,cost,out)
      implicit none
      double precision a,b,x,y,out,cost
      integer p

      a=(1d0-x)*(1d0-y)+x*y
      b=(1d0-x)*(1d0-y)-x*y

      if(p.eq.1)then
      out=-cost*(1d0+cost)/(a**2-b**2*cost**2)
      out=out*(4d0/3d0*b**2-3d0/2d0*a)
      elseif(p.eq.2)then
      out=cost*(1d0-cost)/(a**2-b**2*cost**2)
      out=out*(4d0/3d0*b**2-3d0/2d0*a)
      elseif(p.eq.3)then
      out=cost**2-2d0*4d0/3d0/3d0*a
      out=out*(a-b**2)/(a**2-b**2*cost**2)
      endif

      return
      end

