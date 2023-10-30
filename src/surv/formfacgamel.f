ccccc EPA form factors (electron)
      subroutine formfacgamel(io,t1,t2,out)
      implicit none
      double precision qsq,qsqp,fe,fm,t1,t2,out
      double precision ww1,ww2,ww1p,ww2p,ww1pa,ww2pa
      integer io

      include 'photo.f'
      include 'ewpars.f'
      include 'pi.f'
      include 'x.f'
      include 'mom.f'
      
      qsq=(x1**2*me**2+t1)/(1d0-x1)
      qsqp=(x2**2*me**2+t2)/(1d0-x2)      

      fe=1d0
      fm=1d0

cccccccccc

      ww1p=t1/qsq**2*fe
      ww1pa=x1**2/2d0*fm/qsq

      ww1=fe/qsq**2
      ww1=ww1/pi*alpha

      ww1p=ww1p/pi*alpha
      ww1pa=ww1pa/pi*alpha

ccccccccccc

      ww2p=t2/qsqp**2*fe
      ww2pa=x2**2/2d0*fm/qsqp

      ww2=fe/qsqp**2
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
