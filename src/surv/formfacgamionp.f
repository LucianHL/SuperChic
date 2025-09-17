      subroutine formfacgamion_ion_em(io,t1,t2,out)
      implicit none
      double precision t1,t2,out
      double precision qsq1,qsq2,q0,qsqmin
      double precision qsq1t,qsq2t
      double precision ww1p,ww2p,ww1,ww2,ww1pa
      double precision tpint,gm,ge,fm,fe
      double precision a2e,a2m,a1m,a1e,a0e,a0m
      integer io
      double precision f1,f2
      double precision x2i,x1i

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
      include 'vars.f'

c       print*,x1*x2*s,mx**2,rts
c      print*,mdiss1,mdiss2

      q0=0.71d0

      qsq1t=(q(4,3)-q(4,1))**2-(q(3,3)-q(3,1))**2-(q(2,3)-q(2,1))**2
     &     -(q(1,3)-q(1,1))**2
      qsq1t=-qsq1t

      qsq2t=(q(4,4)-q(4,2))**2-(q(3,4)-q(3,2))**2-(q(2,4)-q(2,2))**2
     &     -(q(1,4)-q(1,2))**2
      qsq2t=-qsq2t

      x1i=x1
      x2i=x2

      x1i=((q(4,1)-q(4,3))/q(4,1))
      x2i=((q(4,2)-q(4,4))/q(4,2))

ccccccccc

      if(ioninel_pbeam.eq.1)then

      qsq1=(x1**2*mp**2+x1*(mdiss1**2-mp**2)+t1)/(1d0-x1)
      qsq2=(x2**2*mion**2+t2)/(1d0-x2)

      ww1=1d0
      ww1p=1d0
      ww1pa=1d0

cccccccccc

      f2=1d0/(1d0+qsq2/q0)**2
      f2=f2/qsq2/(1d0-x2i)
      f2=f2*tpint(1,dsqrt(qsq2))
      f2=f2*dsqrt((1d0-x2i)/137d0/pi)

      ww2=f2
      ww2p=f2*dsqrt(t2)

      else

      qsq1=(x1**2*mion**2+t1)/(1d0-x1)
      qsq2=(x2**2*mp**2+x2*(mdiss2**2-mp**2)+t2)/(1d0-x2)

      ww1=1d0
      ww1p=1d0
      ww1pa=1d0

cccccccccc

      f2=1d0/(1d0+qsq1/q0)**2
      f2=f2/qsq1/(1d0-x1i)
      f2=f2*tpint(1,dsqrt(qsq1))
      f2=f2*dsqrt((1d0-x1i)/137d0/pi)

      ww2=f2
      ww2p=f2*dsqrt(t1)

      endif

ccccccccccccccc

      if(io.eq.1)then
         out=dsqrt(ww1)*ww2
         out=1d0
      elseif(io.eq.2)then
         out=dsqrt(ww1pa)*ww2p
         out=1d0
      endif

      if(neutron_inel)then
      out=out*dsqrt(an-az)
      else
      out=out*dsqrt(az)
      endif

      return
      end

ccccc EPA form factors (proton)
      subroutine formfacgamion_inel(io,t1,t2,out)
      implicit none
      double precision t1,t2,out
      double precision qsq1,qsq2,q0,qsqmin
      double precision qsq1t,qsq2t
      double precision ww1p,ww2p,ww1,ww2,ww1pa
      double precision tpint,gm,ge,fm,fe
      double precision a2e,a2m,a1m,a1e,a0e,a0m
      integer io
      double precision f1,f2
      double precision x2i,x1i,alphaem,pf,qrec,crec

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
      include 'vars.f'

      q0=0.71d0

      qsq1t=(q_rf(4,3)-q_rf(4,1))**2-(q_rf(3,3)-q_rf(3,1))**2
     &      -(q_rf(2,3)-q_rf(2,1))**2
     &     -(q_rf(1,3)-q_rf(1,1))**2
      qsq1t=-qsq1t

      qsq2t=(q_rf(4,4)-q_rf(4,2))**2-(q_rf(3,4)-q_rf(3,2))**2
     &      -(q_rf(2,4)-q_rf(2,2))**2
     &     -(q_rf(1,4)-q_rf(1,2))**2
      qsq2t=-qsq2t

      x1i=x1
      x2i=x2

      x1i=((q(4,1)-q(4,3))/q(4,1))
      x2i=((q(4,2)-q(4,4))/q(4,2))

      pf=0.25d0

ccccccccc

      if(ioninel_pbeam.eq.1)then

      qsq1=(x1**2*mp**2+x1*(mdiss1**2-mp**2)+t1)/(1d0-x1)
      qsq2=(x2**2*mion**2+t2)/(1d0-x2)

      call F1F2(diss1,xb1,qsq1,mdiss1,f1,f2)

      qrec=qsq1**2/4d0/mp**2+qsq1
      qrec=dsqrt(qrec)

      crec=0.75d0*qrec/pf*(1d0-(qrec/pf)**2/12d0)

      if(qrec.lt.2d0*pf)then
      if(diss1)then
      else
      f1=f1*crec
      f2=f2*crec
      endif
      endif

      ww1p=t1/qsq1**2*f2
      ww1pa=x1i**2/xb1*f1/qsq1

      ww1=f2/qsq1**2
      ww1=ww1/pi*alpha

      ww1p=ww1p/pi*alpha
      ww1pa=ww1pa/pi*alpha

cccccccccc

      f2=1d0/(1d0+qsq2/q0)**2
      f2=f2/qsq2/(1d0-x2i)
      f2=f2*tpint(1,dsqrt(qsq2))
      f2=f2*dsqrt((1d0-x2i)/137d0/pi)

      ww2=f2
      ww2p=f2*dsqrt(t2)

      else

      qsq1=(x1**2*mion**2+t1)/(1d0-x1)
      qsq2=(x2**2*mp**2+x2*(mdiss2**2-mp**2)+t2)/(1d0-x2)

      call F1F2(diss2,xb2,qsq2,mdiss2,f1,f2)

      ww1p=t2/qsq2**2*f2
      ww1pa=x2i**2/xb2*f1/qsq2

      ww1=f2/qsq2**2
      ww1=ww1/pi*alpha

      ww1p=ww1p/pi*alpha
      ww1pa=ww1pa/pi*alpha

cccccccccc

      f2=1d0/(1d0+qsq1/q0)**2
      f2=f2/qsq1/(1d0-x1i)
      f2=f2*tpint(1,dsqrt(qsq1))
      f2=f2*dsqrt((1d0-x1i)/137d0/pi)

      ww2=f2
      ww2p=f2*dsqrt(t1)

      endif

ccccccccccccccc

      if(io.eq.1)then
         out=dsqrt(ww1)*ww2
      elseif(io.eq.2)then
         out=dsqrt(ww1pa)*ww2p
      endif

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

