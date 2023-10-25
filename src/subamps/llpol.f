ccc   gamma gamma --> l+l- subprocess amplitude
      subroutine llpol(p,mx,u,t,pp,mm,pm,mp)
      implicit none
      double precision sphi,cphi,sintt,costt,cost_test
      double precision normp,beta,u,t,mx
      complex*16 pp,mm,pm,mp
      integer p

      include 'ewpars.f'
      include 'pi.f'
      include 'mq.f'
      include 'partonmom2.f'
      include 'zi.f'
      include 'mom.f'
      include 'norm.f'
      include 'proc.f'

c      print*,1d0/alpha
      
      beta=dsqrt(1d0-4d0*mq**2/mx**2)
      costt=(t-u)/beta/mx**2

      cost_test=p1(3)/dsqrt(p1(1)**2+p1(2)**2+p1(3)**2)
c      print*,costt,cost_test
c      print*,'beta =',beta,mq,mq/mx
 
c      costt=cost_test
     
      if(costt.gt.1d0)costt=1d0
      if(costt.lt.-1d0)costt=-1d0
      
      sintt=dsqrt(1d0-costt**2)

      if((p1(1)**2+p1(2)**2).gt.0d0)then
         cphi=p1(1)/dsqrt(p1(1)**2+p1(2)**2)
         sphi=p1(2)/dsqrt(p1(1)**2+p1(2)**2)
      else
         cphi=1d0
         sphi=0d0
      endif

      if(p.eq.1)then !++
         pp=4d0*mq*(1d0+beta)/mx
         pp=pp*mx**4/4d0/(t-mq**2)/(u-mq**2)
         mm=4d0*mq*(-1d0+beta)/mx
         mm=mm*mx**4/4d0/(t-mq**2)/(u-mq**2)
         pm=4d0*mq/mx*beta*sintt**2/(1d0-beta**2*costt**2)
         mp=4d0*mq/mx*beta*sintt**2/(1d0-beta**2*costt**2)
         
      elseif(p.eq.2)then !--

         pp=4d0*mq*(1d0-beta)/mx
         pp=pp*mx**4/4d0/(t-mq**2)/(u-mq**2)
         mm=4d0*mq*(-1d0-beta)/mx
         mm=mm*mx**4/4d0/(t-mq**2)/(u-mq**2)
         pm=-4d0*mq/mx*beta*sintt**2/(1d0-beta**2*costt**2)
         mp=-4d0*mq/mx*beta*sintt**2/(1d0-beta**2*costt**2)

      elseif(p.eq.3)then   !+-
         
         pp=0d0
         mm=0d0
         pm=-2d0*beta*sintt*(1d0-costt)/(1d0-beta**2*costt**2)
         mp=2d0*beta*sintt*(1d0+costt)/(1d0-beta**2*costt**2)

      elseif(p.eq.4)then  !-+

         pp=0d0
         mm=0d0
         mp=-2d0*beta*sintt*(1d0-costt)/(1d0-beta**2*costt**2)
         pm=2d0*beta*sintt*(1d0+costt)/(1d0-beta**2*costt**2)

c         mp=1d0
c         pm=1d0

c         mp=1d0
c         pm=2d0*beta*sintt*(1d0+costt)/(1d0-beta**2*costt**2)

      endif

      if(proc.eq.76)then
         pp=2d0*(1d0-beta**2)/(1d0-beta**2*costt**2)
         mm=pp
         pm=2d0*(1d0-(1d0-beta**2)/(1d0-beta**2*costt**2))
         mp=pm
      endif

c      normp=16d0*pi**2*alpha**2
      normp=16d0*pi**2/137d0**2
      normp=normp/32d0/pi/mx**2*beta/4d0
      normp=normp*2d0
      normp=normp*conv
      normp=dsqrt(normp)

      if(proc.eq.61)normp=normp*dsqrt(3d0)*(2d0/3d0)**2 ! ttbar
      if(proc.eq.76)normp=normp*dsqrt(2d0) ! L + R sleptons
      
      pm=pm*(cphi+zi*sphi)**2
      mp=mp*(cphi-zi*sphi)**2

      pp=pp*normp
      mm=mm*normp
      pm=pm*normp
      mp=mp*normp
 
c      pp=0d0
c      mm=0d0
      
      if(p.eq.1)then
c         pm=0d0
c         mp=0d0
      endif

      if(p.eq.2)then
c         pm=0d0
c         mp=0d0
      endif

      if(p.eq.3)then
c         pm=0d0
c         mp=0d0
      endif

      return
      end
