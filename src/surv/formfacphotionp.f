ccccc EPA form factors (proton)
      subroutine formfacphotionp(io,t1,t2,out)
      implicit none
      double precision qsq1,qsq2,tpint,q0,out,t1,t2
      double precision wt,ww2,ww2p,x2i
      double precision f2,mion
      integer io

      include 'photo.f'
      include 'mp.f'
      include 'pi.f'
      include 'x.f'
      include 'mn.f'
      include 'ion.f'
      include 'bpsi.f'

      mion=mp*az+(an-az)*mn
      x2i=x2/an
      
      qsq1=(x1**2*mp**2+t1)/(1d0-x1)
      qsq2=(x2i**2*mion**2+t2)/(1d0-x2i)
      q0=0.71d0

ccccccccc

      wt=dexp(-bpsi*t1/2d0)

cccccccccc

      f2=1d0/(1d0+qsq2/q0)**2
      f2=f2/(t2+x2i**2*mion**2)
      f2=f2*tpint(1,dsqrt(qsq2))
      f2=f2*dsqrt((1d0-x2i)/137d0/pi*2d0)
      
      ww2=f2
      ww2p=f2*dsqrt(t2)

      if(io.eq.1)then
         out=wt*ww2
      elseif(io.eq.2)then
         out=wt*ww2p
      endif

      return
      end
