ccc   calculates 'bare' amplitude, performs integration over gluon Q_t
      subroutine bare(mu,pt1x,pt1y,pt2x,pt2y,out)
      implicit none
      double precision qt1(2),qt2(2)
      complex*16 sum(10),wt(10)
      integer n1,ntot1,n2,ntot2,ntot3,p
      complex*16 out(10),qm0,qp2,qm2
      complex*16 wts(0:1000,10)
      complex*16 zpp,zpm,zmp,zmm
      double precision yinc,y,wtt
      double precision qtsq,qtx,qty,qqphiinc,qqmin,qqmax2,
     &     qqmax1,qqmax,qqinc2,qqinc1,qqinc,qphiinc,qphi,
     &     qp0,qmin,qmax1,qmax2,qmax,qinc1,qinc2,qinc,q2sq,
     &     q2min,q1min,q1sq,lamq,jac
      double precision fg1,fg2,fg
      double precision pt1x,pt1y,pt2x,pt2y,mu

      include 'polvecs.f'
      include 'pi.f'
      include 'x.f'
      include 'zi.f'
      include 'vars.f'
      include 'mandelstam.f'
      include 'polarization.f'
      include 'proc.f'
      include 'mq.f'
      include 'forward.f'
      include 'scorr.f'
      include 'zarr.f'
      include 'quarkonia.f'
      include 'ewpars.f'
      include 'jz2.f'
      include 'gaussvars.f'
      include 'nsurv.f'

cccccccccccc

      ntot1=45
      if(jz2)ntot1=90
      ntot2=10
      ntot3=9

      ntot1=s2int*2
      ntot2=s2int
      ntot3=s2int*2

      qmin=0.4d0
      qmax=8d0
      qmax1=qmax
      qmax2=100d0

      qinc=(qmax-qmin)/dble(ntot1)
      qphiinc=2d0*pi/dble(ntot2)

      qinc1=(qmax1-qmin)/dble(ntot1)
      qinc2=(qmax2-qmax1)/dble(ntot3)

cccccc

      qqmin=dsqrt(qmin)
      qqmax=dsqrt(qmax)
      qqmax1=dsqrt(qmax1)
      qqmax2=dsqrt(qmax2)

      qqinc=(qqmax-qqmin)/dble(ntot1)
      qqphiinc=2d0*pi/dble(ntot2)

      qqinc1=(qqmax1-qqmin)/dble(ntot1)
      qqinc2=(qqmax2-qqmax1)/dble(ntot3)

ccccccc

      jac=1d0
      yinc=(1d0/qmax2-1d0/qmax1)/dble(ntot3)

      do p=1,pol
         sum(p)=0d0
      enddo

      lamq=10d0/3d0
      lamq=-1d0

      do n1=0,ntot1+ntot3+1
         do p=1,pol
            wts(n1,p)=0d0
         enddo
      enddo

      do 333 n1=1,ntot1+ntot3+1
         do 333 n2=1,ntot2

      qphi=pi*(xiphib(n2)+1d0)

      if(n1.lt.(ntot1+1))then
         qinc=qqinc1
         qtsq=(qqmax1-qqmin)*xib(n1)/2d0+(qqmax1+qqmin)/2d0
         qtsq=qtsq**2
         jac=2d0*dsqrt(qtsq)
         jac=jac*wiphib(n2)*wib(n1)*pi*(qqmax1-qqmin)/2d0
      else
         qinc=1d0
         y=(1d0/qmax2-1d0/qmax1)*xib(n1-ntot1)/2d0
     &        +(1d0/qmax2+1d0/qmax1)/2d0
         qtsq=1d0/y
         jac=qtsq**2
         jac=-jac*wiphib(n2)*wib(n1-ntot1)*pi*(1d0/qmax2-1d0/qmax1)/2d0
      endif

cccccccccccc

      qtx=dsqrt(qtsq)*dcos(qphi)
      qty=dsqrt(qtsq)*dsin(qphi)

      if(forward)then
         q1sq=qtsq
         q2sq=qtsq
         qp0=((qtx)*(qtx)+(qty)*(qty))/2d0
         qm0=((qtx)*(qty)-(qty)*(qtx))*zi/2d0
         qp2=((qty)*(qty)-(qtx)*(qtx))/2d0
         qm2=qp2+zi*((qtx)*(qty)+(qty)*(qtx))/2d0
         qp2=qp2-zi*((qtx)*(qty)+(qty)*(qtx))/2d0
      else
         q1sq=(qtx-pt1x)**2+(qty-pt1y)**2
         q2sq=(qtx+pt2x)**2+(qty+pt2y)**2
         qp0=((qtx-pt1x)*(qtx+pt2x)+(qty-pt1y)*(qty+pt2y))/2d0
         qm0=((qtx-pt1x)*(qty+pt2y)-(qty-pt1y)*(qtx+pt2x))*zi/2d0
         qp2=((qty-pt1y)*(qty+pt2y)-(qtx-pt1x)*(qtx+pt2x))/2d0
         qm2=qp2+zi*((qtx-pt1x)*(qty+pt2y)+(qty-pt1y)*(qtx+pt2x))/2d0
         qp2=qp2-zi*((qtx-pt1x)*(qty+pt2y)+(qty-pt1y)*(qtx+pt2x))/2d0
      endif

      if(qtsq.gt.q1sq)then
         q1min=q1sq
      else
         q1min=qtsq
      endif

      if(qtsq.gt.q2sq)then
         q2min=q2sq
      else
         q2min=qtsq
      endif

      if(q1min.lt.qmin) goto 334
      if(q2min.lt.qmin) goto 334

ccccccccccc

      fg1=fg(x1,q1min,mu)
      fg2=fg(x2,q2min,mu)

      if(fg1.lt.0d0)goto 334
      if(fg2.lt.0d0)goto 334

      if(hel.eq.2)then
         qt1(1)=qtx-pt1x
         qt1(2)=qty-pt1y
         qt2(1)=-qtx-pt2x
         qt2(2)=-qty-pt2y
      endif

      if(forward)then
         wtt=fg1*fg2/qtsq**3
      else
         wtt=fg1*fg2/qtsq/q1sq/q2sq
      endif

      wtt=wtt*jac
      wtt=wtt*pi**2/mx**2

ccccccc

      do p=1,pol

         if(proc.eq.21.or.proc.eq.24.or.proc.eq.29.or.proc.eq.32
     &        .or.proc.eq.35)then
            call chi0(mx,mx/2d0,qt1,qt2,zpp)
         elseif(proc.eq.22.or.proc.eq.25.or.proc.eq.27.or.proc.eq.
     &           30.or.proc.eq.33.or.proc.eq.36)then
            call chi1(p,mx,mx/2d0,mchic0,qt1,qt2,echi1,zpp)
         elseif(proc.eq.23.or.proc.eq.26.or.proc.eq.28.or.proc.eq.
     &           31.or.proc.eq.34.or.proc.eq.37)then
            call chi2(p,mx/2d0,mchic0,qt1,qt2,echi2,zpp)
         elseif(proc.eq.38)then
            call etaq(mx,qt1,qt2,zpp)
         endif

        if(proc.eq.39.or.proc.eq.42)then
            call chi0(mx,mx/2d0,qt1,qt2,zpp)
         elseif(proc.eq.40.or.proc.eq.43.or.proc.eq.45)then
            call chi1(p,mx,mx/2d0,mchib0,qt1,qt2,echi1,zpp)
         elseif(proc.eq.41.or.proc.eq.44.or.proc.eq.46)then
            call chi2(p,mx/2d0,mchib0,qt1,qt2,echi2,zpp)
         elseif(proc.eq.47)then
            call etaq(mx,qt1,qt2,zpp)
         endif


         if(hel.eq.1)then
            zpp=zarr(1,p)
            zmm=zarr(2,p)
            zpm=zarr(3,p)
            zmp=zarr(4,p)
            wt(p)=((zpp+zmm)*qp0+(zpp-zmm)*qm0+qp2*zpm+qm2*zmp)
         else
            wt(p)=zpp
         endif

       wt(p)=wt(p)*wtt

       sum(p)=sum(p)+wt(p)

      enddo

cccccc

 334  continue

 333  enddo


      do p=1,pol
         out(p)=sum(p)
      enddo

c      print*,out(1)

c      stop

      return
      end
