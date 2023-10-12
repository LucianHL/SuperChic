      subroutine F1F2(diss,x,q2,mx,f1,f2)
      implicit double precision(a-y)
      logical diss

      include 'mp.f'
      include 'pi.f'
      
      wcut=dsqrt(3.5d0)
      q0=1d0
      
      if(diss)then
         if(mx.lt.wcut)then
            call f2sf(1,q2,x,f2)
            if(q2.gt.q0)then
               f1=f2*(1d0+4d0*mp**2*x**2/q2)/(1d0+rfit(2,x,q2))/2d0/x
            else
               f1=f2*(1d0+4d0*mp**2*x**2/q0)/(1d0+rfit(2,x,q0))/2d0/x
            endif
         else
            if(q2.lt.q0)then
               call gd_fit_11(x,q2,mx**2,f2)
               f1=f2*(1d0+4d0*mp**2*x**2/q0)/(1d0+rfit(2,x,q0))/2d0/x
            else
               call F1F2ap(x,q2,f2,fl)
               f1=f2*(1d0+4d0*mp**2*x**2/q2)-fl
               f1=f1/2d0/x             
            endif
         endif
      else
         call F1F2el(q2,f1,f2)
      endif

      if(f1.lt.0d0)f1=0d0
      if(f2.lt.0d0)f2=0d0

      return
      end
