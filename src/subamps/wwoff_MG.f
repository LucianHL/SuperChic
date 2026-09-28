ccc   gamma gamma --> l+l- subprocess amplitude - off-shell
      subroutine wwoff_MG(p)
      implicit none
      integer p,i,i1,i2,j
      double precision q1(4),q2(4)
      REAL*8 Pmom(0:3,6)
      integer nhel(4)
      double precision alphaem,qsq1,qsq2,beta
      complex*16 zout,AMP_aaWW_SM,ztt1

      include 'mom.f'
      include 'vars.f'
      include 'pi.f'
      include 'norm.f'
      include 'wwpars.f'
      include 'xb.f'
      include 'pol.f'
      include 'egam0.f'
      include 'mp.f'
      include 'zoutarr.f'
      include 'ewpars.f'

      beta=dsqrt(1d0-4d0*mw**2/mx**2)

      qsq1=(q(4,3)-q(4,1))**2-(q(3,3)-q(3,1))**2-(q(2,3)-q(2,1))**2
     &     -(q(1,3)-q(1,1))**2
      qsq1=-qsq1
      qsq2=(q(4,4)-q(4,2))**2-(q(3,4)-q(3,2))**2-(q(2,4)-q(2,2))**2
     &     -(q(1,4)-q(1,2))**2
      qsq2=-qsq2

      do i=1,4
         q1(i)=q(i,1)-q(i,3)
         q2(i)=q(i,2)-q(i,4)
      enddo

      do i=1,3
         pmom(i,1)=q1(i)
         pmom(i,2)=q2(i)
         pmom(i,3)=q(i,6)
         pmom(i,4)=q(i,7)
      enddo

      pmom(0,1)=q1(4)
      pmom(0,2)=q2(4)
      pmom(0,3)=q(4,6)
      pmom(0,4)=q(4,7)


      if(p.eq.3)THEN
         nhel(3)=-1  
         nhel(4)=1 
      elseif(p.eq.4)then
         nhel(3)=1
         nhel(4)=-1
      elseif(p.eq.1)then
         nhel(3)=1
         nhel(4)=1
      elseif(p.eq.2)then
         nhel(3)=-1
         nhel(4)=-1
      elseif(p.eq.5)then
         nhel(3)=0
         nhel(4)=1
      elseif(p.eq.6)then
         nhel(3)=0
         nhel(4)=-1
      elseif(p.eq.7)then
         nhel(3)=1
         nhel(4)=0
      elseif(p.eq.8)then
         nhel(3)=-1
         nhel(4)=0
      elseif(p.eq.9)then
         nhel(3)=0
         nhel(4)=0
      endif


      do i1=1,4
            do i2=1,4

              

            call egcalc(i1,i2)
            zcalc=.true.
            zout=AMP_aaWW_SM(Pmom,nhel)
            zcalc=.false.
            zout=zout*dsqrt(alphaEM(qsq1)*alphaEM(qsq2))
            zout=zout*1.325070D+02
            zout=zout*dsqrt(conv)
            zout=zout*dsqrt(beta)
c           zoutarr_mg(p,i1,i2)=zout
             zoutarr(p,i1,i2)=zout

c           print*,zout/zoutarr(p,i1,i2)

            ENDDO
      enddo

c      stop

c$$$       zout=0d0
c$$$       do i=1,4
c$$$          zout=0d0
c$$$          do j=1,4
c$$$             ztt1=zoutarr_mg(p,i,j)*q2(j)
c$$$c$$$*     q2(i)
c$$$             if(j.lt.4)ztt1=-ztt1
c$$$c$$$c            if(i.lt.4)ztt1=-ztt1
c$$$c            print*,ztt1
c$$$             zout=zout+ztt1
c$$$          enddo
c$$$          print*,i,zout
c$$$       enddo


      
      return
      end

