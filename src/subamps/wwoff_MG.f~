ccc   gamma gamma --> l+l- subprocess amplitude - off-shell
      subroutine lloff_MG(p)
      implicit none
      integer p,i,i1,i2
      double precision q1(4),q2(4)
      REAL*8 Pmom(0:3,6)
      integer nhel(4)
      double precision alphaem,qsq1,qsq2
      complex*16 zout,AMP_aall_SM

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
      include 'tau.f'


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
      endif


      do i1=1,4
            do i2=1,4

              

            call egcalc(i1,i2)
            zcalc=.true.
            zout=AMP_aall_SM(Pmom,nhel)
            zout=zout*dsqrt(alphaEM(qsq1)*alphaEM(qsq2))
            zout=zout*1.325070D+02
            zout=zout*dsqrt(conv)
            zout=zout
c            zoutarr_mg(p,i1,i2)=zout
            zoutarr(p,i1,i2)=zout



            ENDDO
      enddo

!       zout=0d0
!       do i=1,4
!          zout=0d0
!          do j=1,4
!             ztt1=zoutarr_mg(p,i,j)*q2(j)
! c$$$*     q2(i)
!             if(j.lt.4)ztt1=-ztt1
! c$$$c            if(i.lt.4)ztt1=-ztt1
! c            print*,ztt1
!             zout=zout+ztt1
!          enddo
!          print*,i,zout
!       enddo


      
      return
      end

      subroutine egcalc(i,j)
      implicit none
      integer i1,i,j
      
      include 'egam0.f'
      
      do i1=1,4
            eg1(i1)=0d0
            eg2(i1)=0d0
      ENDDO
                  
      eg1(i)=1d0
      eg2(j)=1d0
      
      if(i.lt.4)eg1(i)=-eg1(i)
      if(j.lt.4)eg2(j)=-eg2(j)


      return
      end
