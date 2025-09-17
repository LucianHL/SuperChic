ccccc EPA form factors (proton)
      subroutine formfacgamion_inel(io,t1,t2,out)
      implicit none
      double precision t1,t2,out
      double precision x2i,qsq1,qsq2,q0,qsqmin
      double precision ww1p,ww2p,ww1,ww2,ww1pa
      double precision tpint,gm,ge,fm,fe
      double precision a2e,a2m,a1m,a1e,a0e,a0m
      integer io
      double precision f1,f2

      include 'photo.f'
      include 'mp.f'
      include 'pi.f'
      include 'x.f'
      include 'mn.f'
      include 'ion.f'
      include 'mion.f'
      include 'mom.f'
      include 'ion_inel.f'
      include 'xb.f'
      include 'diss.f'
      include 'ewpars.f'
      include 'diff.f'


c      mion=mp*az+(an-az)*mn
c      mion=an*mp
c      x2i=x2/an

      x2i=x2

c      x1=x1*mion/mp

      qsq1=(x1**2*mp**2+t1)/(1d0-x1)
      qsq2=(x2i**2*mion**2+t2)/(1d0-x2i)

      q0=0.71d0

ccccccccc

      qsqmin=mp**2*x1**2/(1d0-x1)

      a0e=0.98462d0
      a1e=0.68414d0
      a2e=0.01933d0
      a0m=0.28231d0
      a1m=1.34919d0
      a2m=0.55473d0
      ge=(a0e/(1d0+qsq1/a1e)**2+(1d0-a0e)/(1d0+qsq1/a2e)**2)**2
      gm=(a0m/(1d0+qsq1/a1m)**2+(1d0-a0m)/(1d0+qsq1/a2m)**2)**2
      gm=gm*7.78d0

      fe=(4d0*mp**2*ge+qsq1*gm)/(4d0*mp**2+qsq1)
      fm=gm

cccc

      ww1p=((1d0-x1)*(1d0-qsqmin/qsq1)*fe)/qsq1
      ww1pa=x1**2/2d0*fm/qsq1

      ww1=((1d0-x1)*fe)/qsq1**2/(1d0-x1)
      ww1=ww1/pi/(1d0-x1)/137d0

      ww1p=ww1p/pi/(1d0-x1)/137d0
      ww1pa=ww1pa/pi/(1d0-x1)/137d0

c      print*,ww1,ww1p,ww1pa

cccccccccc

ccccccccc

c      print*,xb1,qsq1,qsq2

      call F1F2(diss1,xb1,qsq1,mdiss1,f1,f2)

c      print*,f1,f2


      ww1p=t1/qsq1**2*f2
      ww1pa=x1**2/xb1*f1/qsq1

      ww1=f2/qsq1**2
      ww1=ww1/pi*alpha

      ww1p=ww1p/pi*alpha
      ww1pa=ww1pa/pi*alpha

c      print*,ww1,ww1p,ww1pa
c      print*,''

cccccccccc


      f2=1d0/(1d0+qsq2/q0)**2
      f2=f2/(t2+x2i**2*mion**2)
      f2=f2*tpint(1,dsqrt(qsq2))
      f2=f2*dsqrt((1d0-x2i)/137d0/pi)

      ww2=f2
      ww2p=f2*dsqrt(t2)

c      print*,f2
c      print*,io,ww1,ww1pa,ww2,ww2p
c      print*,''

      if(io.eq.1)then
         out=dsqrt(ww1)*ww2
      elseif(io.eq.2)then
         out=dsqrt(ww1pa)*ww2p
      endif


c      x1=x1*mp/mion
      if(neutron_inel)then
      out=out*dsqrt(an-az)
      else
      out=out*dsqrt(az)
      endif

      return
      end



ccccc EPA form factors (proton)
      subroutine formfacgamionp(io,t1,t2,out)
      implicit none
      double precision t1,t2,out
      double precision x2i,qsq1,qsq2,q0,qsqmin
      double precision ww1p,ww2p,ww1,ww2,ww1pa
      double precision tpint,gm,ge,fm,fe,f2
      double precision a2e,a2m,a1m,a1e,a0e,a0m
      integer io

      include 'photo.f'
      include 'mp.f'
      include 'pi.f'
      include 'x.f'
      include 'mn.f'
      include 'ion.f'
      include 'mion.f'
      include 'mom.f'
      include 'ion_inel.f'

c      mion=mp*az+(an-az)*mn
c      mion=an*mp
c      x2i=x2/an

      x2i=x2

      qsq1=(x1**2*mp**2+t1)/(1d0-x1)
      qsq2=(x2i**2*mion**2+t2)/(1d0-x2i)

      q0=0.71d0

ccccccccc

      qsqmin=mp**2*x1**2/(1d0-x1)

      a0e=0.98462d0
      a1e=0.68414d0
      a2e=0.01933d0
      a0m=0.28231d0
      a1m=1.34919d0
      a2m=0.55473d0
      ge=(a0e/(1d0+qsq1/a1e)**2+(1d0-a0e)/(1d0+qsq1/a2e)**2)**2
      gm=(a0m/(1d0+qsq1/a1m)**2+(1d0-a0m)/(1d0+qsq1/a2m)**2)**2
      gm=gm*7.78d0

      fe=(4d0*mp**2*ge+qsq1*gm)/(4d0*mp**2+qsq1)
      fm=gm

cccc

      ww1p=((1d0-x1)*(1d0-qsqmin/qsq1)*fe)/qsq1
      ww1pa=x1**2/2d0*fm/qsq1

      ww1=((1d0-x1)*fe)/qsq1**2/(1d0-x1)
      ww1=ww1/pi/(1d0-x1)/137d0

      ww1p=ww1p/pi/(1d0-x1)/137d0
      ww1pa=ww1pa/pi/(1d0-x1)/137d0

cccccccccc


      f2=1d0/(1d0+qsq2/q0)**2
      f2=f2/(t2+x2i**2*mion**2)
      f2=f2*tpint(1,dsqrt(qsq2))
      f2=f2*dsqrt((1d0-x2i)/137d0/pi)

      ww2=f2
      ww2p=f2*dsqrt(t2)

      if(io.eq.1)then
         out=dsqrt(ww1)*ww2
      elseif(io.eq.2)then
         out=dsqrt(ww1pa)*ww2p
      endif

      return
      end
