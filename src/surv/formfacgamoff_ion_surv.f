ccccc EPA form factors (proton)
      subroutine formfacgamoff_ion_surv(p,q1x,q1y,q2x,q2y,zout)
      implicit none
      complex*16 zt,zout12,zout22,zout21,zout11,zout
      double precision ww2,ww2pa,ww1,ww1pa,alphaem,tpint
      double precision t1,t2,qsq,qsqp,q0,q1x,q2x,q1y,q2y
      double precision out,out22,f2
      double precision q1(2),q2(2)
      integer p,i,j
 
 
      include 'photo.f'
      include 'mp.f'
      include 'pi.f'
      include 'x.f'
      include 'xb.f'
      include 'diss.f'
      include 'vars.f'
      include 'mom.f'
      include 'ewpars.f'
      include 'proc.f'
      include 'norm.f'
      include 'ppamp.f'
      include 'zi.f'
      include 'zoutarr.f'
      include 'surv.f'
      include 'fbeam.f'
      include 'mion.f'

      q0=0.71d0
      
      t1=q1x**2+q1y**2
      t2=q2x**2+q2y**2
      
c      qsq=(q(4,3)-q(4,1))**2-(q(3,3)-q(3,1))**2-(q(2,3)-q(2,1))**2
c     &     -(q(1,3)-q(1,1))**2
c      qsq=-qsq
c      qsqp=(q(4,4)-q(4,2))**2-(q(3,4)-q(3,2))**2-(q(2,4)-q(2,2))**2
c     &     -(q(1,4)-q(1,2))**2
c      qsqp=-qsqp
      
      q1(1)=q1x
      q1(2)=q1y
      q2(1)=q2x
      q2(2)=q2y

       qsq=(x1**2*mion**2+t1)/(1d0-x1)
       qsqp=(x2**2*mion**2+t2)/(1d0-x2)

c       qsq=(q(4,3)-q(4,1))**2-(q(3,3)-q(3,1))**2-(q(2,3)-q(2,1))**2
c     &      -(q(1,3)-q(1,1))**2
c       qsq=-qsq
c       qsqp=(q(4,4)-q(4,2))**2-(q(3,4)-q(3,2))**2-(q(2,4)-q(2,2))**2
c     &      -(q(1,4)-q(1,2))**2
c       qsqp=-qsqp


c       qsq=qsq*(1d0-qsq*mp**2/4d0/s**2/xb1**2/(1d0-x1))
c       qsqp=qsqp*(1d0-qsqp*mp**2/4d0/s**2/xb2**2/(1d0-x2))
       
ccccccccc

       
       f2=1d0/(1d0+qsq/q0)**2
       f2=f2*tpint(1,dsqrt(qsq))
       f2=f2**2
       

       ww1=2d0*f2/qsq
       ww1pa=0d0
      
cccccccccc
       
       f2=1d0/(1d0+qsqp/q0)**2
       f2=f2*tpint(1,dsqrt(qsqp))
       f2=f2**2
       
       ww2=2d0*f2/qsqp
       ww2pa=0d0
       
ccccccccc

c      if(io.eq.2)goto 111
      
      zout22=0d0               

      do i=1,2
         do j=1,2
            zt=q1(i)*q2(j)*zoutarr(p,i,j)/x1/x2
            zout22=zout22+zt
         enddo
      enddo

      out22=dsqrt(ww1*ww2)
      zout22=zout22*out22      

c      if(io.eq.1)goto 222
      
c 111  zout11=0d0
c      do i=1,4
c         do j=1,4
c            zt=zoutarr(p,i,j)*conjg(zoutarr(p,i,j))
c            if(i.lt.4)zt=-zt
c            if(j.lt.4)zt=-zt
c            zout11=zout11+zt
c         enddo
c      enddo

c      out11=ww1pa*ww2pa
c      zout11=zout11*out11

      zout11=0d0

c      zout12=0d0
c      do i=1,2
c         do j=1,2
c            do k=1,4
c               zt=zoutarr(p,k,i)*conjg(zoutarr(p,k,j))*q2(i)*q2(j)/x2**2
c               if(k.lt.4)zt=-zt
c               if(i.lt.4)zt=-zt
c               if(j.lt.4)zt=-zt
c               zout12=zout12+zt
c            enddo
c         enddo
c      enddo

c      out12=ww1pa*ww2
c      zout12=zout12*out12
      
      zout12=0d0

      
c      zout21=0d0
c      do i=1,2
c         do j=1,2
c            do k=1,4
c               zt=zoutarr(p,i,k)*conjg(zoutarr(p,j,k))*q1(i)*q1(j)/x1**2
c               if(k.lt.4)zt=-zt
c               if(i.lt.4)zt=-zt
c               if(j.lt.4)zt=-zt
c               zout21=zout21+zt
c            enddo
c         enddo
c      enddo
      
c      out21=ww1*ww2pa
c      zout21=zout21*out21

      zout21=0d0

cccccccc

c 222  if(io.eq.1)then
         zout=zout22
c      elseif(io.eq.2)then
c         zout=dsqrt(cdabs(zout11)+cdabs(zout12)+cdabs(zout21))
c         zout=dsqrt(cdabs(zout11-zout12-zout21))
c      endif
      
      out=dsqrt(cdabs(zout22)**2+cdabs(zout11)+cdabs(zout12)
     &     +cdabs(zout21))

      out=out*dsqrt(alphaEM(qsq)*alphaEM(qsqp)/qsq/qsqp)
      out=out*dsqrt(4d0)        ! rho normalisation

      zout=zout*dsqrt(alphaEM(qsq)*alphaEM(qsqp)/qsq/qsqp)
      zout=zout*dsqrt(4d0)        ! rho normalisation
 
c      print*,alphaEM(qsq),1d0/137d0
     
      return
      end
