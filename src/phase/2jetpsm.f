ccc   generates two-body phase space for particles of mass mq1,mq2
      subroutine twojetpsm(mx,mq1,mq2,ps,u,t)
      implicit none
      double precision rphi,rtheta,ran2,phi,stheta,ctheta
      double precision beta,ps,u,t,mq1,mq2,mx
      double precision px(4),pcm(4)
      double precision pboo(4)
      integer i

      include 'mt.f'
      include 'pi.f'
      include 'mom.f'
      include 'partonmom2.f'

      rphi=ran2()
      rtheta=ran2()

      phi=2d0*pi*rphi
      ctheta=-1d0+2d0*rtheta
      stheta=dsqrt(1d0-ctheta**2)

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

      mt=dsqrt((mq1+mq2)**2/4d0+p1(1)**2+p1(2)**2)

      return
      end
