ccc   gg --> rhorho subprocess amplitude
      subroutine rhorho(p,mx,t,pp,mm,pm,mp)
      implicit none
      double precision phi,cphi,sphi
      double precision out1,out2,out3,out4,out5,out6,out7
      double precision out1a,out2a,out3a,out4a,out5a,out6a,out7a
      double precision norm,cost,fpip,beta,a48n,a28n,mx,t,alphas
      complex*16 pp,mm,pm,mp
      integer p,pa

      include 'partonmom2.f'
      include 'pi.f'
      include 'zi.f'
      include 'mq.f'
      include 'mixing.f'

      phi=datan2(p1(2),p1(1))
      cphi=dcos(phi)
      sphi=dsin(phi)


      fpip=fpi/(2d0*dsqrt(6d0))

      norm=64d0*pi**2*alphas(mx**2/4d0)**2
      norm=norm*fpip**2*6d0**2
      norm=norm/3d0/mx**2

      beta=dsqrt(1d0-4d0*mq**2/mx**2)
      cost=(1d0+2d0*t/mx**2-2d0*mq**2/mx**2)/beta

      if(p.eq.1)pa=2
      if(p.eq.2)pa=1

      call mesint(p,1,cost,out1)
      call mesint(p,2,cost,out2)
      call mesint(p,3,cost,out3)
      call mesint(p,4,cost,out4)
      call mesint(p,5,cost,out5)
      call mesint(p,6,cost,out6)
      call mesint(p,7,cost,out7)

      call mesint(pa,1,cost,out1a)
      call mesint(pa,2,cost,out2a)
      call mesint(pa,3,cost,out3a)
      call mesint(pa,4,cost,out4a)
      call mesint(pa,5,cost,out5a)
      call mesint(pa,6,cost,out6a)
      call mesint(pa,7,cost,out7a)

      call wfoctet(mx,2,a28,a28n)
      call wfoctet(mx,4,a48,a48n)
      out2=out2*(a28n/a28)
      out3=out3*(a28n/a28)**2
      
      if(dabs(a48n).gt.1d-4)then
      out4=out4*(a48n/a48)
      out5=out5*(a48n/a48)*(a28n/a28)
      out6=out6*(a48n/a48)*(a28n/a28)
      out7=out7*(a48n/a48)**2
      endif
      
      out2a=out2a*(a28n/a28)
      out3a=out3a*(a28n/a28)**2
      if(dabs(a48n).gt.1d-4)then
      out4a=out4a*(a48n/a48)
      out5a=out5a*(a48n/a48)*(a28n/a28)
      out6a=out6a*(a48n/a48)*(a28n/a28)
      out7a=out7a*(a48n/a48)**2
      endif

      if(dabs(a28n).lt.1d-4)then
         out2=0d0
         out3=0d0
         out5=0d0
         out6=0d0
         out2a=0d0
         out3a=0d0
         out5a=0d0
         out6a=0d0
      endif

      if(dabs(a48n).lt.1d-4)then
         out4=0d0
         out5=0d0
         out6=0d0
         out7=0d0
         out4a=0d0
         out5a=0d0
         out6a=0d0
         out7a=0d0
      endif

      pp=0d0
      mm=0d0
      if(p.lt.3)then
         pm=(out1+2d0*out2+out3)*norm
         pm=pm+(out4*2d0+out5*2d0+out6)*norm
         mp=(out1a+2d0*out2a+out3a)*norm
         mp=mp+(out4a*2d0+out5a*2d0+out6a)*norm
      else
         pm=(out1+2d0*out2+out3)*norm*3d0/2d0
         pm=pm+(out4*2d0+out5*2d0+out6)*norm*3d0/2d0
         mp=pm
      endif

      pm=pm*(cphi+zi*sphi)**2
      mp=mp*(cphi-zi*sphi)**2

      return
      end

