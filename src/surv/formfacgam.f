ccccc EPA form factors (proton)
      subroutine formfacgam(io,t1,t2,out)
      implicit none
      double precision ww1pa,ww2pa,ww1p,ww2p,ww1,ww2
      double precision qsq,qsqp,f1,f2,out,t1,t2
      integer i1,i2,io
      double precision garr(-6:6)

      include 'photo.f'
      include 'mp.f'
      include 'pi.f'
      include 'x.f'
      include 'xb.f'
      include 'diss.f'
      include 'vars.f'
      include 'mom.f'
      include 'ewpars.f'
      
      qsq=(x1**2*mp**2+t1)/(1d0-x1)
      qsqp=(x2**2*mp**2+t2)/(1d0-x2) 
      
ccccccccc
      
      call F1F2(diss1,xb1,qsq,mdiss1,f1,f2)
      
      ww1p=t1/qsq**2*f2
      ww1pa=x1**2/xb1*f1/qsq

      ww1=f2/qsq**2
      ww1=ww1/pi*alpha

      ww1p=ww1p/pi*alpha
      ww1pa=ww1pa/pi*alpha

cccccccccc
      
      call F1F2(diss2,xb2,qsqp,mdiss2,f1,f2)

      ww2p=t2/qsqp**2*f2
      ww2pa=x2**2/xb2*f1/qsqp

      ww2=f2/qsqp**2
      ww2=ww2/pi*alpha

      ww2p=ww2p/pi*alpha
      ww2pa=ww2pa/pi*alpha

      if(io.eq.1)then
         out=dsqrt(ww1*ww2)
      elseif(io.eq.2)then
         out=dsqrt(ww1p*ww2pa+ww2p*ww1pa+ww1pa*ww2pa)
      endif
      
      return
      end
