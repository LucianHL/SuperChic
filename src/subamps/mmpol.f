ccc   gamma gamma --> MM subprocess amplitude
      subroutine mmpol(p,mx,u,t,pp,mm,pm,mp)
      implicit none
      double precision sphi,cphi,sintt,costt
      double precision normp,gam,beta,u,t,mx
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
      
      beta=dsqrt(1d0-4d0*mq**2/mx**2)
      gam=mx**2/mq**2
      costt=(t-u)/beta/mx**2
      sintt=dsqrt(1d0-costt**2)

      cphi=p1(1)/dsqrt(p1(1)**2+p1(2)**2)
      sphi=p1(2)/dsqrt(p1(1)**2+p1(2)**2)

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

      endif

      normp=pi**2/alpha**2
      normp=normp/32d0/pi/mx**2*beta/4d0
      normp=normp*2d0
      normp=normp*conv
      if(proc.eq.72)normp=normp*beta**4 
      normp=dsqrt(normp)
      
      pm=pm*(cphi+zi*sphi)**2
      mp=mp*(cphi-zi*sphi)**2

      pp=pp*normp
      mm=mm*normp
      pm=pm*normp
      mp=mp*normp

      return
      end
