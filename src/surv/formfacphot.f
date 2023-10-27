ccccc Photoproduction form factors (assume Pomeron coupling universal
ccccc between eigenstates
      subroutine formfacphot(io,t1,t2,out)
      implicit none
      double precision ww1,ww2,wt
      double precision qsq,qsqmin,out1,out2,t1,t2,out
      double precision ge,gm,fe,fm,a2e,a2m,a1m,a1e,a0e,a0m
      integer io

      include 'bpsi.f'
      include 'photo.f'
      include 'mp.f'
      include 'pi.f'
      
      if(prot.eq.1)then
         wt=dexp(-bpsi*t2)
         qsq=(xgam**2*mp**2+t1)/(1d0-xgam)
      else
         wt=dexp(-bpsi*t1) 
         qsq=(xgam**2*mp**2+t2)/(1d0-xgam)
      endif

      qsqmin=mp**2*xgam**2/(1d0-xgam)

      a0e=0.98462d0
      a1e=0.68414d0
      a2e=0.01933d0
      a0m=0.28231d0
      a1m=1.34919d0
      a2m=0.55473d0
      ge=(a0e/(1d0+qsq/a1e)**2+(1d0-a0e)/(1d0+qsq/a2e)**2)**2
      gm=(a0m/(1d0+qsq/a1m)**2+(1d0-a0m)/(1d0+qsq/a2m)**2)**2
      gm=gm*7.78d0
      
      fe=(4d0*mp**2*ge+qsq*gm)/(4d0*mp**2+qsq)
      fm=gm

      ww1=fe/qsq**2
      ww2=xgam**2/2d0*fm/qsq

      ww1=ww1/pi/(1d0-xgam)/137d0*2d0
      out1=dsqrt(wt*ww1)

      ww2=ww2/pi/(1d0-xgam)/137d0*2d0
      out2=dsqrt(wt*ww2)

      if(io.eq.1)then
         out=out1
      else
         out=out2
      endif

      return
      end
