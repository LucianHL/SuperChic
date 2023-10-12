ccc   gg --> pipi subprocess amplitude
      subroutine pipi(p,mx,t,pp,mm,pm,mp)
      implicit none
      double precision phi,cphi,sphi,out1,out2,out3
      double precision norm,fpip,cost,beta,a28n
      double precision t,mx,alphas
      complex*16 pp,mm,pm,mp
      integer p

      include 'partonmom2.f'
      include 'pi.f'
      include 'zi.f'
      include 'mq.f'
      include 'mixing.f'
 
cccccc

      phi=datan(p1(2)/p1(1))
      cphi=p1(1)/dsqrt(p1(1)**2+p1(2)**2)
      sphi=p1(2)/dsqrt(p1(1)**2+p1(2)**2)

ccccccc

      fpip=fpi/(2d0*dsqrt(6d0))

      norm=64d0*pi**2*alphas(mx**2/4d0)**2
      norm=norm*fpip**2*6d0**2
      norm=norm/2d0/mx**2

      beta=dsqrt(1d0-4d0*mq**2/mx**2)
      cost=(1d0+2d0*t/mx**2-2d0*mq**2/mx**2)/beta

      call mesint(p,1,cost,out1)
      call mesint(p,2,cost,out2)
      call mesint(p,3,cost,out3)

      call wfoctet(mx,2,a28,a28n)
      
      out2=out2*(a28n/a28)
      out3=out3*(a28n/a28)**2
      
      if(dabs(a28n).lt.1d-4)then
         out2=0d0
         out3=0d0
      endif

      pp=0d0
      mm=0d0
      pm=(out1+2d0*out2+out3)*norm
      mp=pm

      pm=pm*(cphi+zi*sphi)**2
      mp=mp*(cphi-zi*sphi)**2

      return
      end

