
ccc   generates two-body phase space for particles of equal mass mq
      subroutine twojetps(mx,mq,rtheta,rphi,ps,u,t)
      implicit none
      double precision phi,stheta,ctheta,beta,u,t,ps
      double precision rphi,mx,mq,rtheta
      double precision px(4),pcm(4)
      double precision pboo(4)
      integer i
      common/cthetatemp/ctheta

      include 'mt.f'
      include 'partonmom2.f'
      include 'pi.f'
      include 'mom.f'

      if ( 2d0*mq .GE. mx ) then
      u=-mx**2/4.0d0
      t=-mx**2/4.0d0
      ps=0.0d0
      return 
      endif

c      rphi=ran2()

c      print*,rphi,rtheta

      phi=2d0*pi*rphi
      ctheta=-1d0+2d0*rtheta
      stheta=dsqrt(1d0-ctheta**2)
      beta=dsqrt(1d0-4d0*mq**2/mx**2)

      pcm(4)=mx/2d0
      pcm(1)=stheta*dcos(phi)*pcm(4)*beta
      pcm(2)=stheta*dsin(phi)*pcm(4)*beta
      pcm(3)=ctheta*pcm(4)*beta

      u=mq**2-mx*(pcm(4)+pcm(3))
      t=mq**2-mx*(pcm(4)-pcm(3))

      do i=1,4
         px(i)=q(i,5)
         p1(i)=pcm(i)
      enddo

      p2(4)=pcm(4)
      do i=1,3
         p2(i)=-pcm(i)
      enddo

      call boost(mx,px,pcm,pboo)

      do i=1,4
         q(i,6)=pboo(i)
         q(i,7)=q(i,5)-q(i,6)
      enddo

      ps=4d0*pi*beta

      mt=dsqrt(mq**2+p1(1)**2+p1(2)**2)

      return
      end

ccc   generates two-body phase space for particles of differing mass mq1,mq2
      subroutine twojetps_rm(mx,mq1,mq2,rtheta,rphi,ps,u,t)
      implicit double precision(a-y)
      double precision px(4),pcm(4)
      double precision pboo(4)
      integer i

      include 'mt.f'
      include 'partonmom2.f'
      include 'pi.f'
      include 'mom.f'
      if ( mq1 + mq2 .GE. mx ) THEN
      u=(mq1**2+mq2**2-mx**2)/2d0
      t=(mq1**2+mq2**2-mx**2)/2d0    
      ps=0.0D0
      endif

      phi=2d0*pi*rphi
      ctheta=-1d0+2d0*rtheta
      stheta=dsqrt(1d0-ctheta**2)
      beta=dsqrt(1d0-4d0*mq**2/mx**2)

      pcm(4)=(mx+(mq1**2-mq2**2)/mx)/2d0
      beta=dsqrt(1d0-mq1**2/pcm(4)**2)

      pcm(1)=stheta*dcos(phi)*pcm(4)*beta
      pcm(2)=stheta*dsin(phi)*pcm(4)*beta
      pcm(3)=ctheta*pcm(4)*beta

      u=mq1**2-mx*(pcm(4)+pcm(3))
      t=mq1**2-mx*(pcm(4)-pcm(3))

      do i=1,4
         px(i)=q(i,5)
         p1(i)=pcm(i)
      enddo

      p2(4)=pcm(4)
      do i=1,3
         p2(i)=-pcm(i)
      enddo

      call boost(mx,px,pcm,pboo)

      do i=1,4
         q(i,6)=pboo(i)
         q(i,7)=q(i,5)-q(i,6)
      enddo

      ps=4d0*pi*dsqrt(mx**4+mq1**4+mq2**4-2d0*mx**2*mq1**2
     &     -2d0*mx**2*mq2**2-2d0*mq1**2*mq2**2)/mx**2

      mt=dsqrt(mq**2+p1(1)**2+p1(2)**2)

      return
      end
