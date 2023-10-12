ccc   Calculates skewed unintegrated PDF according to prescpription of
ccc   LHL, Phys.Rev. D88 (2013) 3, 034029, arXiv:1306.6661
ccc   hg : skewed pdf
ccc   diffhg : dhg/dQ^2

      subroutine hpdf(x,qsq,hg,diffhg)
      implicit double precision (a-z)
      integer ytot,iy

      include 'pi.f'

      eps=1d-3
      ytot=100

      sum=0d0
      sum1=0d0

      qmin=1d0
      xmax=1d0


      if(qsq.gt.qmin.and.x.lt.xmax)then

      yinc=(1d0-x/4d0)/dble(ytot)

      do iy=1,ytot
        
         y=x/4d0+(dble(iy)-0.5d0)*yinc
         xg0=xg(x/4d0/y,qsq)
         xgp=xg(x/4d0/y,qsq+eps)
         wt=dsqrt(y**3*(1d0-y))*xg0
         diffg=(xgp-xg0)/eps
         wt1=dsqrt(y**3*(1d0-y))*diffg

         if(dabs(wt).gt.0d0)then
         else
            wt=0d0
         endif
         if(dabs(wt1).gt.0d0)then
         else
            wt1=0d0
         endif

         sum=sum+wt*yinc
         sum1=sum1+wt1*yinc
      
      enddo

      hg=sum*16d0/pi
      diffhg=sum1*16d0/pi

      else

         hg=xg(x,qsq)
         diffhg=dxg(x,qmin)

      endif

      return
      end
