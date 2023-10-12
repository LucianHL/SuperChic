ccc   gg --> pipi subprocess amplitude (x,y dependent)
      subroutine pipixy(x,y,cost,out)
      implicit none
      double precision x,y,out,cost,a,b

      a=(1d0-x)*(1d0-y)+x*y
      b=(1d0-x)*(1d0-y)-x*y

      out=cost**2-2d0*4d0/3d0/3d0*a
      out=out*(a-b**2)/(a**2-b**2*cost**2)

      return
      end

