ccc   gg --> J/psi psi(2S) subprocess amplitude
      subroutine psipsip(p,mu,mx,u,t,pp,mm,pm,mp)
      implicit none
      complex*16 pp,mm,pm,mp
      double precision sphi,sintp,s,rts,mx,phi,cphi,cost,beta
      double precision qf,psi01s,psi02s,psi0,gampsip,gampsi
      double precision ppuu,pptt,pmuu,pmtt,mpuu,mptt,mmuu,mmtt
      double precision m,denu,dent
      double precision alpha,u,t,mu,alphas
      integer p

      include 'partonmom2.f'
      include 'pi.f'
      include 'zi.f'
      include 'mq.f'
      include 'quarkonia.f'

cccccccc

      gampsi=92.9d-6*5.95d-2
      qf=2d0/3d0
      alpha=1d0/137d0
      psi01s=mpsi*dsqrt(gampsi)/2d0/alpha/qf

      gampsip=303d-6*7.82d-3
      psi02s=mpsip*dsqrt(gampsip)/2d0/alpha/qf

      psi0=dsqrt(psi01s*psi02s)

      rts=mx
      s=rts**2
      m=mq

cccccccc

      beta=dsqrt(1d0-4d0*mq**2/mx**2)
      cost=p1(3)/dsqrt(p1(1)**2+p1(2)**2+p1(3)**2)

      phi=datan2(p1(2),p1(1))
      cphi=dcos(phi)
      sphi=dsin(phi)


      sintp=dsqrt(p1(1)**2+p1(2)**2)/p1(4)
     &*dsqrt(mx**2/8d0)

      denu=(m**2-t)**2*(u+t-2d0*m**2)**3
      dent=(m**2-u)**2*(u+t-2d0*m**2)**3


cccccccccccc

      if(p.eq.1)then     ! +- f.s.

      pmuu =
     & 7./4.*s**3*m**2 + 7./2.*s**3*m**2*cost + 7./4.*s**3*m**2*cost**2
     &  - 17./16.*s**3*m**2*beta*cost - 17./8.*s**3*m**2*beta*cost**2
     &  - 17./16.*s**3*m**2*beta*cost**3 - 15./8.*s**3*m**2*beta**2*
     & cost**2 - 15./4.*s**3*m**2*beta**2*cost**3 - 15./8.*s**3*m**2*
     & beta**2*cost**4 + 15./16.*s**3*m**2*beta**3*cost**3 + 15./8.*
     & s**3*m**2*beta**3*cost**4 + 15./16.*s**3*m**2*beta**3*cost**5 +
     & 9./8.*s**4*beta + 9./8.*s**4*beta*cost - 9./8.*s**4*beta*cost**2
     &  - 9./8.*s**4*beta*cost**3 - 9./8.*s**4*beta**3*cost**2 - 9./8.*
     & s**4*beta**3*cost**3 + 9./8.*s**4*beta**3*cost**4 + 9./8.*s**4*
     & beta**3*cost**5

      pmtt =
     & 7./4.*s**3*m**2 + 7./2.*s**3*m**2*cost + 7./4.*s**3*m**2*cost**2
     &  + 17./16.*s**3*m**2*beta*cost + 17./8.*s**3*m**2*beta*cost**2
     &  + 17./16.*s**3*m**2*beta*cost**3 - 15./8.*s**3*m**2*beta**2*
     & cost**2 - 15./4.*s**3*m**2*beta**2*cost**3 - 15./8.*s**3*m**2*
     & beta**2*cost**4 - 15./16.*s**3*m**2*beta**3*cost**3 - 15./8.*
     & s**3*m**2*beta**3*cost**4 - 15./16.*s**3*m**2*beta**3*cost**5 -
     & 9./8.*s**4*beta - 9./8.*s**4*beta*cost + 9./8.*s**4*beta*cost**2
     &  + 9./8.*s**4*beta*cost**3 + 9./8.*s**4*beta**3*cost**2 + 9./8.*
     & s**4*beta**3*cost**3 - 9./8.*s**4*beta**3*cost**4 - 9./8.*s**4*
     & beta**3*cost**5

      mpuu =
     & 7./4.*s**3*m**2 - 7./2.*s**3*m**2*cost + 7./4.*s**3*m**2*cost**2
     &  - 17./16.*s**3*m**2*beta*cost + 17./8.*s**3*m**2*beta*cost**2
     &  - 17./16.*s**3*m**2*beta*cost**3 - 15./8.*s**3*m**2*beta**2*
     & cost**2 + 15./4.*s**3*m**2*beta**2*cost**3 - 15./8.*s**3*m**2*
     & beta**2*cost**4 + 15./16.*s**3*m**2*beta**3*cost**3 - 15./8.*
     & s**3*m**2*beta**3*cost**4 + 15./16.*s**3*m**2*beta**3*cost**5 -
     & 9./8.*s**4*beta + 9./8.*s**4*beta*cost + 9./8.*s**4*beta*cost**2
     &  - 9./8.*s**4*beta*cost**3 + 9./8.*s**4*beta**3*cost**2 - 9./8.*
     & s**4*beta**3*cost**3 - 9./8.*s**4*beta**3*cost**4 + 9./8.*s**4*
     & beta**3*cost**5

      mptt =
     & 7./4.*s**3*m**2 - 7./2.*s**3*m**2*cost + 7./4.*s**3*m**2*cost**2
     &  + 17./16.*s**3*m**2*beta*cost - 17./8.*s**3*m**2*beta*cost**2
     &  + 17./16.*s**3*m**2*beta*cost**3 - 15./8.*s**3*m**2*beta**2*
     & cost**2 + 15./4.*s**3*m**2*beta**2*cost**3 - 15./8.*s**3*m**2*
     & beta**2*cost**4 - 15./16.*s**3*m**2*beta**3*cost**3 + 15./8.*
     & s**3*m**2*beta**3*cost**4 - 15./16.*s**3*m**2*beta**3*cost**5 +
     & 9./8.*s**4*beta - 9./8.*s**4*beta*cost - 9./8.*s**4*beta*cost**2
     &  + 9./8.*s**4*beta*cost**3 - 9./8.*s**4*beta**3*cost**2 + 9./8.*
     & s**4*beta**3*cost**3 + 9./8.*s**4*beta**3*cost**4 - 9./8.*s**4*
     & beta**3*cost**5

      ppuu =
     & 15./16.*s**3*m**2*beta*cost - 15./16.*s**3*m**2*beta*cost**3 -
     & 15./8.*s**3*m**2*beta**2*cost**2 + 15./8.*s**3*m**2*beta**2*
     & cost**4 + 15./16.*s**3*m**2*beta**3*cost**3 - 15./16.*s**3*m**2*
     & beta**3*cost**5

      pptt =
     &  - 15./16.*s**3*m**2*beta*cost + 15./16.*s**3*m**2*beta*cost**3
     &  - 15./8.*s**3*m**2*beta**2*cost**2 + 15./8.*s**3*m**2*beta**2*
     & cost**4 - 15./16.*s**3*m**2*beta**3*cost**3 + 15./16.*s**3*m**2*
     & beta**3*cost**5

      mmuu =
     & 15./16.*s**3*m**2*beta*cost - 15./16.*s**3*m**2*beta*cost**3 -
     & 15./8.*s**3*m**2*beta**2*cost**2 + 15./8.*s**3*m**2*beta**2*
     & cost**4 + 15./16.*s**3*m**2*beta**3*cost**3 - 15./16.*s**3*m**2*
     & beta**3*cost**5

      mmtt =
     &  - 15./16.*s**3*m**2*beta*cost + 15./16.*s**3*m**2*beta*cost**3
     &  - 15./8.*s**3*m**2*beta**2*cost**2 + 15./8.*s**3*m**2*beta**2*
     & cost**4 - 15./16.*s**3*m**2*beta**3*cost**3 + 15./16.*s**3*m**2*
     & beta**3*cost**5

      mm=mmtt/dent+mmuu/denu
      pp=pptt/dent+ppuu/denu
      pm=pmtt/dent+pmuu/denu
      mp=mptt/dent+mpuu/denu

ccc   Get relative phases right

      pm=pm*(cphi+zi*sphi)**2
      mp=mp*(cphi-zi*sphi)**2
ccc


      pm=pm*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      mp=mp*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      pp=pp*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      mm=mm*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m

      pm=pm*dsqrt(2d0)
      mp=mp*dsqrt(2d0)
      pp=pp*dsqrt(2d0)
      mm=mm*dsqrt(2d0)

      elseif(p.eq.2)then    ! ++ f.s.

      mpuu =
     & 7./4.*s**3*m**2 - 7./4.*s**3*m**2*cost**2 - 17./16.*s**3*m**2*
     & beta*cost + 17./16.*s**3*m**2*beta*cost**3 + 1./4.*s**3*m**2*
     & beta**2 - 17./8.*s**3*m**2*beta**2*cost**2 + 15./8.*s**3*m**2*
     & beta**2*cost**4 + 15./16.*s**3*m**2*beta**3*cost**3 - 15./16.*
     & s**3*m**2*beta**3*cost**5 - 9./8.*s**4*beta*cost + 9./8.*s**4*
     & beta*cost**3 + 9./8.*s**4*beta**3*cost - 9./8.*s**4*beta**3*
     & cost**5 - 9./8.*s**4*beta**5*cost**3 + 9./8.*s**4*beta**5*
     & cost**5

      mptt =
     & 7./4.*s**3*m**2 - 7./4.*s**3*m**2*cost**2 + 17./16.*s**3*m**2*
     & beta*cost - 17./16.*s**3*m**2*beta*cost**3 + 1./4.*s**3*m**2*
     & beta**2 - 17./8.*s**3*m**2*beta**2*cost**2 + 15./8.*s**3*m**2*
     & beta**2*cost**4 - 15./16.*s**3*m**2*beta**3*cost**3 + 15./16.*
     & s**3*m**2*beta**3*cost**5 + 9./8.*s**4*beta*cost - 9./8.*s**4*
     & beta*cost**3 - 9./8.*s**4*beta**3*cost + 9./8.*s**4*beta**3*
     & cost**5 + 9./8.*s**4*beta**5*cost**3 - 9./8.*s**4*beta**5*
     & cost**5

      pmuu =
     & 7./4.*s**3*m**2 - 7./4.*s**3*m**2*cost**2 - 17./16.*s**3*m**2*
     & beta*cost + 17./16.*s**3*m**2*beta*cost**3 + 1./4.*s**3*m**2*
     & beta**2 - 17./8.*s**3*m**2*beta**2*cost**2 + 15./8.*s**3*m**2*
     & beta**2*cost**4 + 15./16.*s**3*m**2*beta**3*cost**3 - 15./16.*
     & s**3*m**2*beta**3*cost**5 - 9./8.*s**4*beta*cost + 9./8.*s**4*
     & beta*cost**3 + 9./8.*s**4*beta**3*cost - 9./8.*s**4*beta**3*
     & cost**5 - 9./8.*s**4*beta**5*cost**3 + 9./8.*s**4*beta**5*
     & cost**5

      pmtt =
     & 7./4.*s**3*m**2 - 7./4.*s**3*m**2*cost**2 + 17./16.*s**3*m**2*
     & beta*cost - 17./16.*s**3*m**2*beta*cost**3 + 1./4.*s**3*m**2*
     & beta**2 - 17./8.*s**3*m**2*beta**2*cost**2 + 15./8.*s**3*m**2*
     & beta**2*cost**4 - 15./16.*s**3*m**2*beta**3*cost**3 + 15./16.*
     & s**3*m**2*beta**3*cost**5 + 9./8.*s**4*beta*cost - 9./8.*s**4*
     & beta*cost**3 - 9./8.*s**4*beta**3*cost + 9./8.*s**4*beta**3*
     & cost**5 + 9./8.*s**4*beta**5*cost**3 - 9./8.*s**4*beta**5*
     & cost**5

      pmuu =
     & 7./4.*s**3*m**2 - 7./4.*s**3*m**2*cost**2 - 17./16.*s**3*m**2*
     & beta*cost + 17./16.*s**3*m**2*beta*cost**3 + 1./4.*s**3*m**2*
     & beta**2 - 17./8.*s**3*m**2*beta**2*cost**2 + 15./8.*s**3*m**2*
     & beta**2*cost**4 + 15./16.*s**3*m**2*beta**3*cost**3 - 15./16.*
     & s**3*m**2*beta**3*cost**5 - 9./8.*s**4*beta*cost + 9./8.*s**4*
     & beta*cost**3 + 9./8.*s**4*beta**3*cost - 9./8.*s**4*beta**3*
     & cost**5 - 9./8.*s**4*beta**5*cost**3 + 9./8.*s**4*beta**5*
     & cost**5

      pmtt =
     & 7./4.*s**3*m**2 - 7./4.*s**3*m**2*cost**2 + 17./16.*s**3*m**2*
     & beta*cost - 17./16.*s**3*m**2*beta*cost**3 + 1./4.*s**3*m**2*
     & beta**2 - 17./8.*s**3*m**2*beta**2*cost**2 + 15./8.*s**3*m**2*
     & beta**2*cost**4 - 15./16.*s**3*m**2*beta**3*cost**3 + 15./16.*
     & s**3*m**2*beta**3*cost**5 + 9./8.*s**4*beta*cost - 9./8.*s**4*
     & beta*cost**3 - 9./8.*s**4*beta**3*cost + 9./8.*s**4*beta**3*
     & cost**5 + 9./8.*s**4*beta**5*cost**3 - 9./8.*s**4*beta**5*
     & cost**5

      ppuu =
     &  - 1./4.*s**3*m**2 - 19./8.*s**3*m**2*cost - 1./2.*s**3*m**2*
     & beta - 53./16.*s**3*m**2*beta*cost + 1./4.*s**3*m**2*beta*
     & cost**2 + 15./16.*s**3*m**2*beta*cost**3 - 1./4.*s**3*m**2*
     & beta**2 + 17./8.*s**3*m**2*beta**2*cost**2 + 17./8.*s**3*m**2*
     & beta**2*cost**3 - 15./8.*s**3*m**2*beta**2*cost**4 + 19./16.*
     & s**3*m**2*beta**3*cost**3 + 15./16.*s**3*m**2*beta**3*cost**5 +
     & 9./8.*s**4*cost + 9./4.*s**4*beta*cost - 9./8.*s**4*beta**2*cost
     &  - 9./4.*s**4*beta**2*cost**2 - 9./8.*s**4*beta**2*cost**3 - 9./
     & 8.*s**4*beta**3*cost - 9./8.*s**4*beta**3*cost**3 + 9./8.*s**4*
     & beta**4*cost**3 + 9./4.*s**4*beta**4*cost**4 + 9./8.*s**4*
     & beta**5*cost**3 - 9./8.*s**4*beta**5*cost**5

      pptt =
     &  - 1./4.*s**3*m**2 + 19./8.*s**3*m**2*cost - 1./2.*s**3*m**2*
     & beta + 53./16.*s**3*m**2*beta*cost + 1./4.*s**3*m**2*beta*
     & cost**2 - 15./16.*s**3*m**2*beta*cost**3 - 1./4.*s**3*m**2*
     & beta**2 + 17./8.*s**3*m**2*beta**2*cost**2 - 17./8.*s**3*m**2*
     & beta**2*cost**3 - 15./8.*s**3*m**2*beta**2*cost**4 - 19./16.*
     & s**3*m**2*beta**3*cost**3 - 15./16.*s**3*m**2*beta**3*cost**5 -
     & 9./8.*s**4*cost - 9./4.*s**4*beta*cost + 9./8.*s**4*beta**2*cost
     &  - 9./4.*s**4*beta**2*cost**2 + 9./8.*s**4*beta**2*cost**3 + 9./
     & 8.*s**4*beta**3*cost + 9./8.*s**4*beta**3*cost**3 - 9./8.*s**4*
     & beta**4*cost**3 + 9./4.*s**4*beta**4*cost**4 - 9./8.*s**4*
     & beta**5*cost**3 + 9./8.*s**4*beta**5*cost**5

      mmuu =
     &  - 1./4.*s**3*m**2 + 19./8.*s**3*m**2*cost + 1./2.*s**3*m**2*
     & beta - 53./16.*s**3*m**2*beta*cost - 1./4.*s**3*m**2*beta*
     & cost**2 + 15./16.*s**3*m**2*beta*cost**3 - 1./4.*s**3*m**2*
     & beta**2 + 17./8.*s**3*m**2*beta**2*cost**2 - 17./8.*s**3*m**2*
     & beta**2*cost**3 - 15./8.*s**3*m**2*beta**2*cost**4 + 19./16.*
     & s**3*m**2*beta**3*cost**3 + 15./16.*s**3*m**2*beta**3*cost**5 -
     & 9./8.*s**4*cost + 9./4.*s**4*beta*cost + 9./8.*s**4*beta**2*cost
     &  - 9./4.*s**4*beta**2*cost**2 + 9./8.*s**4*beta**2*cost**3 - 9./
     & 8.*s**4*beta**3*cost - 9./8.*s**4*beta**3*cost**3 - 9./8.*s**4*
     & beta**4*cost**3 + 9./4.*s**4*beta**4*cost**4 + 9./8.*s**4*
     & beta**5*cost**3 - 9./8.*s**4*beta**5*cost**5

      mmtt =
     &  - 1./4.*s**3*m**2 - 19./8.*s**3*m**2*cost + 1./2.*s**3*m**2*
     & beta + 53./16.*s**3*m**2*beta*cost - 1./4.*s**3*m**2*beta*
     & cost**2 - 15./16.*s**3*m**2*beta*cost**3 - 1./4.*s**3*m**2*
     & beta**2 + 17./8.*s**3*m**2*beta**2*cost**2 + 17./8.*s**3*m**2*
     & beta**2*cost**3 - 15./8.*s**3*m**2*beta**2*cost**4 - 19./16.*
     & s**3*m**2*beta**3*cost**3 - 15./16.*s**3*m**2*beta**3*cost**5 +
     & 9./8.*s**4*cost - 9./4.*s**4*beta*cost - 9./8.*s**4*beta**2*cost
     &  - 9./4.*s**4*beta**2*cost**2 - 9./8.*s**4*beta**2*cost**3 + 9./
     & 8.*s**4*beta**3*cost + 9./8.*s**4*beta**3*cost**3 + 9./8.*s**4*
     & beta**4*cost**3 + 9./4.*s**4*beta**4*cost**4 - 9./8.*s**4*
     & beta**5*cost**3 + 9./8.*s**4*beta**5*cost**5


      mm=mmtt/dent+mmuu/denu
      pp=pptt/dent+ppuu/denu
      pm=pmtt/dent+pmuu/denu
      mp=mptt/dent+mpuu/denu


ccc   Get relative phases right

      pm=pm*(cphi+zi*sphi)**2
      mp=mp*(cphi-zi*sphi)**2

ccc

      pm=pm*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      mp=mp*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      pp=pp*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      mm=mm*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m

      pm=pm*dsqrt(2d0)
      mp=mp*dsqrt(2d0)
      pp=pp*dsqrt(2d0)
      mm=mm*dsqrt(2d0)

      elseif(p.eq.3)then  ! 00 f.s.

      ppuu =
     &  - 1./16.*s**4 + 19./32.*s**4*beta*cost - 15./32.*s**4*beta*
     & cost**3 + 1./8.*s**4*beta**2 - 1./16.*s**4*beta**2*cost**2 + 15./
     & 16.*s**4*beta**2*cost**4 - 17./16.*s**4*beta**3*cost - 17./32.*
     & s**4*beta**3*cost**3 - 15./32.*s**4*beta**3*cost**5 - 1./16.*
     & s**4*beta**4 + s**4*beta**4*cost**2 + 1./16.*s**4*beta**5*
     & cost**3 - 9./8.*s**5*m**(-2)*beta*cost + 27./16.*s**5*m**(-2)*
     & beta**2*cost**2 + 27./32.*s**5*m**(-2)*beta**3*cost + 9./32.*
     & s**5*m**(-2)*beta**3*cost**3 - 9./16.*s**5*m**(-2)*beta**4*
     & cost**2 - 27./16.*s**5*m**(-2)*beta**4*cost**4 - 9./32.*s**5*
     & m**(-2)*beta**5*cost - 9./16.*s**5*m**(-2)*beta**5*cost**3 + 27./
     & 32.*s**5*m**(-2)*beta**5*cost**5 + 9./16.*s**5*m**(-2)*beta**6*
     & cost**4 + 9./32.*s**5*m**(-2)*beta**7*cost**3 - 9./32.*s**5*
     & m**(-2)*beta**7*cost**5

      pptt =
     &  - 1./16.*s**4 - 19./32.*s**4*beta*cost + 15./32.*s**4*beta*
     & cost**3 + 1./8.*s**4*beta**2 - 1./16.*s**4*beta**2*cost**2 + 15./
     & 16.*s**4*beta**2*cost**4 + 17./16.*s**4*beta**3*cost + 17./32.*
     & s**4*beta**3*cost**3 + 15./32.*s**4*beta**3*cost**5 - 1./16.*
     & s**4*beta**4 + s**4*beta**4*cost**2 - 1./16.*s**4*beta**5*
     & cost**3 + 9./8.*s**5*m**(-2)*beta*cost + 27./16.*s**5*m**(-2)*
     & beta**2*cost**2 - 27./32.*s**5*m**(-2)*beta**3*cost - 9./32.*
     & s**5*m**(-2)*beta**3*cost**3 - 9./16.*s**5*m**(-2)*beta**4*
     & cost**2 - 27./16.*s**5*m**(-2)*beta**4*cost**4 + 9./32.*s**5*
     & m**(-2)*beta**5*cost + 9./16.*s**5*m**(-2)*beta**5*cost**3 - 27./
     & 32.*s**5*m**(-2)*beta**5*cost**5 + 9./16.*s**5*m**(-2)*beta**6*
     & cost**4 - 9./32.*s**5*m**(-2)*beta**7*cost**3 + 9./32.*s**5*
     & m**(-2)*beta**7*cost**5

      mmuu =
     &  - 1./16.*s**4 + 19./32.*s**4*beta*cost - 15./32.*s**4*beta*
     & cost**3 + 1./8.*s**4*beta**2 - 1./16.*s**4*beta**2*cost**2 + 15./
     & 16.*s**4*beta**2*cost**4 - 17./16.*s**4*beta**3*cost - 17./32.*
     & s**4*beta**3*cost**3 - 15./32.*s**4*beta**3*cost**5 - 1./16.*
     & s**4*beta**4 + s**4*beta**4*cost**2 + 1./16.*s**4*beta**5*
     & cost**3 - 9./8.*s**5*m**(-2)*beta*cost + 27./16.*s**5*m**(-2)*
     & beta**2*cost**2 + 27./32.*s**5*m**(-2)*beta**3*cost + 9./32.*
     & s**5*m**(-2)*beta**3*cost**3 - 9./16.*s**5*m**(-2)*beta**4*
     & cost**2 - 27./16.*s**5*m**(-2)*beta**4*cost**4 - 9./32.*s**5*
     & m**(-2)*beta**5*cost - 9./16.*s**5*m**(-2)*beta**5*cost**3 + 27./
     & 32.*s**5*m**(-2)*beta**5*cost**5 + 9./16.*s**5*m**(-2)*beta**6*
     & cost**4 + 9./32.*s**5*m**(-2)*beta**7*cost**3 - 9./32.*s**5*
     & m**(-2)*beta**7*cost**5

      mmtt =
     &  - 1./16.*s**4 - 19./32.*s**4*beta*cost + 15./32.*s**4*beta*
     & cost**3 + 1./8.*s**4*beta**2 - 1./16.*s**4*beta**2*cost**2 + 15./
     & 16.*s**4*beta**2*cost**4 + 17./16.*s**4*beta**3*cost + 17./32.*
     & s**4*beta**3*cost**3 + 15./32.*s**4*beta**3*cost**5 - 1./16.*
     & s**4*beta**4 + s**4*beta**4*cost**2 - 1./16.*s**4*beta**5*
     & cost**3 + 9./8.*s**5*m**(-2)*beta*cost + 27./16.*s**5*m**(-2)*
     & beta**2*cost**2 - 27./32.*s**5*m**(-2)*beta**3*cost - 9./32.*
     & s**5*m**(-2)*beta**3*cost**3 - 9./16.*s**5*m**(-2)*beta**4*
     & cost**2 - 27./16.*s**5*m**(-2)*beta**4*cost**4 + 9./32.*s**5*
     & m**(-2)*beta**5*cost + 9./16.*s**5*m**(-2)*beta**5*cost**3 - 27./
     & 32.*s**5*m**(-2)*beta**5*cost**5 + 9./16.*s**5*m**(-2)*beta**6*
     & cost**4 - 9./32.*s**5*m**(-2)*beta**7*cost**3 + 9./32.*s**5*
     & m**(-2)*beta**7*cost**5

      pmuu =
     &  - 7./8.*s**4 + 7./8.*s**4*cost**2 + 17./32.*s**4*beta*cost - 17.
     & /32.*s**4*beta*cost**3 - 3./16.*s**4*beta**2 + 9./8.*s**4*
     & beta**2*cost**2 - 15./16.*s**4*beta**2*cost**4 - 15./32.*s**4*
     & beta**3*cost**3 + 15./32.*s**4*beta**3*cost**5 + 1./16.*s**4*
     & beta**4 - 1./16.*s**4*beta**4*cost**2 + 9./16.*s**5*m**(-2)*beta
     & *cost - 9./16.*s**5*m**(-2)*beta*cost**3 - 27./32.*s**5*m**(-2)*
     & beta**3*cost + 9./32.*s**5*m**(-2)*beta**3*cost**3 + 9./16.*s**5
     & *m**(-2)*beta**3*cost**5 + 9./32.*s**5*m**(-2)*beta**5*cost + 9./
     & 16.*s**5*m**(-2)*beta**5*cost**3 - 27./32.*s**5*m**(-2)*beta**5*
     & cost**5 - 9./32.*s**5*m**(-2)*beta**7*cost**3 + 9./32.*s**5*
     & m**(-2)*beta**7*cost**5

      pmtt =
     &  - 7./8.*s**4 + 7./8.*s**4*cost**2 - 17./32.*s**4*beta*cost + 17.
     & /32.*s**4*beta*cost**3 - 3./16.*s**4*beta**2 + 9./8.*s**4*
     & beta**2*cost**2 - 15./16.*s**4*beta**2*cost**4 + 15./32.*s**4*
     & beta**3*cost**3 - 15./32.*s**4*beta**3*cost**5 + 1./16.*s**4*
     & beta**4 - 1./16.*s**4*beta**4*cost**2 - 9./16.*s**5*m**(-2)*beta
     & *cost + 9./16.*s**5*m**(-2)*beta*cost**3 + 27./32.*s**5*m**(-2)*
     & beta**3*cost - 9./32.*s**5*m**(-2)*beta**3*cost**3 - 9./16.*s**5
     & *m**(-2)*beta**3*cost**5 - 9./32.*s**5*m**(-2)*beta**5*cost - 9./
     & 16.*s**5*m**(-2)*beta**5*cost**3 + 27./32.*s**5*m**(-2)*beta**5*
     & cost**5 + 9./32.*s**5*m**(-2)*beta**7*cost**3 - 9./32.*s**5*
     & m**(-2)*beta**7*cost**5


      mpuu =
     &  - 7./8.*s**4 + 7./8.*s**4*cost**2 + 17./32.*s**4*beta*cost - 17.
     & /32.*s**4*beta*cost**3 - 3./16.*s**4*beta**2 + 9./8.*s**4*
     & beta**2*cost**2 - 15./16.*s**4*beta**2*cost**4 - 15./32.*s**4*
     & beta**3*cost**3 + 15./32.*s**4*beta**3*cost**5 + 1./16.*s**4*
     & beta**4 - 1./16.*s**4*beta**4*cost**2 + 9./16.*s**5*m**(-2)*beta
     & *cost - 9./16.*s**5*m**(-2)*beta*cost**3 - 27./32.*s**5*m**(-2)*
     & beta**3*cost + 9./32.*s**5*m**(-2)*beta**3*cost**3 + 9./16.*s**5
     & *m**(-2)*beta**3*cost**5 + 9./32.*s**5*m**(-2)*beta**5*cost + 9./
     & 16.*s**5*m**(-2)*beta**5*cost**3 - 27./32.*s**5*m**(-2)*beta**5*
     & cost**5 - 9./32.*s**5*m**(-2)*beta**7*cost**3 + 9./32.*s**5*
     & m**(-2)*beta**7*cost**5

      mptt =
     &  - 7./8.*s**4 + 7./8.*s**4*cost**2 - 17./32.*s**4*beta*cost + 17.
     & /32.*s**4*beta*cost**3 - 3./16.*s**4*beta**2 + 9./8.*s**4*
     & beta**2*cost**2 - 15./16.*s**4*beta**2*cost**4 + 15./32.*s**4*
     & beta**3*cost**3 - 15./32.*s**4*beta**3*cost**5 + 1./16.*s**4*
     & beta**4 - 1./16.*s**4*beta**4*cost**2 - 9./16.*s**5*m**(-2)*beta
     & *cost + 9./16.*s**5*m**(-2)*beta*cost**3 + 27./32.*s**5*m**(-2)*
     & beta**3*cost - 9./32.*s**5*m**(-2)*beta**3*cost**3 - 9./16.*s**5
     & *m**(-2)*beta**3*cost**5 - 9./32.*s**5*m**(-2)*beta**5*cost - 9./
     & 16.*s**5*m**(-2)*beta**5*cost**3 + 27./32.*s**5*m**(-2)*beta**5*
     & cost**5 + 9./32.*s**5*m**(-2)*beta**7*cost**3 - 9./32.*s**5*
     & m**(-2)*beta**7*cost**5

ccccccccccc

      mm=mmtt/dent+mmuu/denu
      pp=pptt/dent+ppuu/denu
      pm=pmtt/dent+pmuu/denu
      mp=mptt/dent+mpuu/denu

ccc   Get relative phases right

      pm=pm*(cphi+zi*sphi)**2
      mp=mp*(cphi-zi*sphi)**2

ccc

      pm=pm*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      mp=mp*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      pp=pp*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      mm=mm*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m

      elseif(p.eq.4)then   ! +0 f.s.

      ppuu =
     &  - 9*s**3*m**(-1)*beta**(-2)*sintp**3 + 9*s**3*m**(-1)*sintp**3
     &  + 9*s**3*m**(-1)*cost**2*sintp**3 - 9*s**3*m**(-1)*beta**2*
     & cost**2*sintp**3 + 19./8.*s**3*m*beta**(-1)*sintp + 19./8.*s**3*
     & m*sintp - 1./4.*s**3*m*cost*sintp - 15./8.*s**3*m*cost**2*sintp
     &  - 1./4.*s**3*m*beta*cost*sintp - 17./8.*s**3*m*beta*cost**2*
     & sintp + 15./4.*s**3*m*beta*cost**3*sintp - 17./8.*s**3*m*beta**2
     & *cost**2*sintp - 15./8.*s**3*m*beta**2*cost**4*sintp - 9./8.*
     & s**4*m**(-1)*beta**(-1)*sintp - 9./4.*s**4*m**(-1)*sintp - 9./8.
     & *s**4*m**(-1)*cost**2*sintp + 9./8.*s**4*m**(-1)*beta*sintp + 9./
     & 2.*s**4*m**(-1)*beta*cost*sintp + 9./8.*s**4*m**(-1)*beta*
     & cost**2*sintp + 9./8.*s**4*m**(-1)*beta**2*cost**2*sintp + 9./8.
     & *s**4*m**(-1)*beta**2*cost**4*sintp - 9./8.*s**4*m**(-1)*beta**3
     & *cost**2*sintp - 9./2.*s**4*m**(-1)*beta**3*cost**3*sintp + 9./8.
     & *s**4*m**(-1)*beta**4*cost**4*sintp

      pptt =
     & 9*s**3*m**(-1)*beta**(-2)*sintp**3 - 9*s**3*m**(-1)*sintp**3 - 9
     & *s**3*m**(-1)*cost**2*sintp**3 + 9*s**3*m**(-1)*beta**2*cost**2*
     & sintp**3 - 19./8.*s**3*m*beta**(-1)*sintp - 19./8.*s**3*m*sintp
     &  - 1./4.*s**3*m*cost*sintp + 15./8.*s**3*m*cost**2*sintp - 1./4.
     & *s**3*m*beta*cost*sintp + 17./8.*s**3*m*beta*cost**2*sintp + 15./
     & 4.*s**3*m*beta*cost**3*sintp + 17./8.*s**3*m*beta**2*cost**2*
     & sintp + 15./8.*s**3*m*beta**2*cost**4*sintp + 9./8.*s**4*m**(-1)
     & *beta**(-1)*sintp + 9./4.*s**4*m**(-1)*sintp + 9./8.*s**4*
     & m**(-1)*cost**2*sintp - 9./8.*s**4*m**(-1)*beta*sintp + 9./2.*
     & s**4*m**(-1)*beta*cost*sintp - 9./8.*s**4*m**(-1)*beta*cost**2*
     & sintp - 9./8.*s**4*m**(-1)*beta**2*cost**2*sintp - 9./8.*s**4*
     & m**(-1)*beta**2*cost**4*sintp + 9./8.*s**4*m**(-1)*beta**3*
     & cost**2*sintp - 9./2.*s**4*m**(-1)*beta**3*cost**3*sintp - 9./8.
     & *s**4*m**(-1)*beta**4*cost**4*sintp


      mmuu =
     &  - 9*s**3*m**(-1)*beta**(-2)*sintp**3 + 9*s**3*m**(-1)*sintp**3
     &  + 9*s**3*m**(-1)*cost**2*sintp**3 - 9*s**3*m**(-1)*beta**2*
     & cost**2*sintp**3 - 19./8.*s**3*m*beta**(-1)*sintp + 19./8.*s**3*
     & m*sintp + 1./4.*s**3*m*cost*sintp - 15./8.*s**3*m*cost**2*sintp
     &  - 1./4.*s**3*m*beta*cost*sintp + 17./8.*s**3*m*beta*cost**2*
     & sintp + 15./4.*s**3*m*beta*cost**3*sintp - 17./8.*s**3*m*beta**2
     & *cost**2*sintp - 15./8.*s**3*m*beta**2*cost**4*sintp + 9./8.*
     & s**4*m**(-1)*beta**(-1)*sintp - 9./4.*s**4*m**(-1)*sintp - 9./8.
     & *s**4*m**(-1)*cost**2*sintp - 9./8.*s**4*m**(-1)*beta*sintp + 9./
     & 2.*s**4*m**(-1)*beta*cost*sintp - 9./8.*s**4*m**(-1)*beta*
     & cost**2*sintp + 9./8.*s**4*m**(-1)*beta**2*cost**2*sintp + 9./8.
     & *s**4*m**(-1)*beta**2*cost**4*sintp + 9./8.*s**4*m**(-1)*beta**3
     & *cost**2*sintp - 9./2.*s**4*m**(-1)*beta**3*cost**3*sintp + 9./8.
     & *s**4*m**(-1)*beta**4*cost**4*sintp

      mmtt =
     & 9*s**3*m**(-1)*beta**(-2)*sintp**3 - 9*s**3*m**(-1)*sintp**3 - 9
     & *s**3*m**(-1)*cost**2*sintp**3 + 9*s**3*m**(-1)*beta**2*cost**2*
     & sintp**3 + 19./8.*s**3*m*beta**(-1)*sintp - 19./8.*s**3*m*sintp
     &  + 1./4.*s**3*m*cost*sintp + 15./8.*s**3*m*cost**2*sintp - 1./4.
     & *s**3*m*beta*cost*sintp - 17./8.*s**3*m*beta*cost**2*sintp + 15./
     & 4.*s**3*m*beta*cost**3*sintp + 17./8.*s**3*m*beta**2*cost**2*
     & sintp + 15./8.*s**3*m*beta**2*cost**4*sintp - 9./8.*s**4*m**(-1)
     & *beta**(-1)*sintp + 9./4.*s**4*m**(-1)*sintp + 9./8.*s**4*
     & m**(-1)*cost**2*sintp + 9./8.*s**4*m**(-1)*beta*sintp + 9./2.*
     & s**4*m**(-1)*beta*cost*sintp + 9./8.*s**4*m**(-1)*beta*cost**2*
     & sintp - 9./8.*s**4*m**(-1)*beta**2*cost**2*sintp - 9./8.*s**4*
     & m**(-1)*beta**2*cost**4*sintp - 9./8.*s**4*m**(-1)*beta**3*
     & cost**2*sintp - 9./2.*s**4*m**(-1)*beta**3*cost**3*sintp - 9./8.
     & *s**4*m**(-1)*beta**4*cost**4*sintp

      pmuu =
     & 9*s**3*m**(-1)*beta**(-2)*sintp**3 - 9*s**3*m**(-1)*sintp**3 - 9
     & *s**3*m**(-1)*cost**2*sintp**3 + 9*s**3*m**(-1)*beta**2*cost**2*
     & sintp**3 + 7./2.*s**3*m*beta**(-1)*sintp + 7./2.*s**3*m*
     & beta**(-1)*cost*sintp - 17./8.*s**3*m*cost*sintp - 17./8.*s**3*m
     & *cost**2*sintp + 1./2.*s**3*m*beta*sintp + 1./2.*s**3*m*beta*
     & cost*sintp - 15./4.*s**3*m*beta*cost**2*sintp - 15./4.*s**3*m*
     & beta*cost**3*sintp + 15./8.*s**3*m*beta**2*cost**3*sintp + 15./8.
     & *s**3*m*beta**2*cost**4*sintp - 9./8.*s**4*m**(-1)*cost*sintp -
     & 9./8.*s**4*m**(-1)*cost**2*sintp + 9./8.*s**4*m**(-1)*beta**2*
     & cost*sintp + 9./8.*s**4*m**(-1)*beta**2*cost**2*sintp + 9./8.*
     & s**4*m**(-1)*beta**2*cost**3*sintp + 9./8.*s**4*m**(-1)*beta**2*
     & cost**4*sintp - 9./8.*s**4*m**(-1)*beta**4*cost**3*sintp - 9./8.
     & *s**4*m**(-1)*beta**4*cost**4*sintp

      pmtt =
     &  - 9*s**3*m**(-1)*beta**(-2)*sintp**3 + 9*s**3*m**(-1)*sintp**3
     &  + 9*s**3*m**(-1)*cost**2*sintp**3 - 9*s**3*m**(-1)*beta**2*
     & cost**2*sintp**3 + 7./2.*s**3*m*beta**(-1)*sintp + 7./2.*s**3*m*
     & beta**(-1)*cost*sintp + 17./8.*s**3*m*cost*sintp + 17./8.*s**3*m
     & *cost**2*sintp + 1./2.*s**3*m*beta*sintp + 1./2.*s**3*m*beta*
     & cost*sintp - 15./4.*s**3*m*beta*cost**2*sintp - 15./4.*s**3*m*
     & beta*cost**3*sintp - 15./8.*s**3*m*beta**2*cost**3*sintp - 15./8.
     & *s**3*m*beta**2*cost**4*sintp + 9./8.*s**4*m**(-1)*cost*sintp +
     & 9./8.*s**4*m**(-1)*cost**2*sintp - 9./8.*s**4*m**(-1)*beta**2*
     & cost*sintp - 9./8.*s**4*m**(-1)*beta**2*cost**2*sintp - 9./8.*
     & s**4*m**(-1)*beta**2*cost**3*sintp - 9./8.*s**4*m**(-1)*beta**2*
     & cost**4*sintp + 9./8.*s**4*m**(-1)*beta**4*cost**3*sintp + 9./8.
     & *s**4*m**(-1)*beta**4*cost**4*sintp

      mpuu =
     & 9*s**3*m**(-1)*beta**(-2)*sintp**3 - 9*s**3*m**(-1)*sintp**3 - 9
     & *s**3*m**(-1)*cost**2*sintp**3 + 9*s**3*m**(-1)*beta**2*cost**2*
     & sintp**3 - 7./2.*s**3*m*beta**(-1)*sintp + 7./2.*s**3*m*
     & beta**(-1)*cost*sintp + 17./8.*s**3*m*cost*sintp - 17./8.*s**3*m
     & *cost**2*sintp - 1./2.*s**3*m*beta*sintp + 1./2.*s**3*m*beta*
     & cost*sintp + 15./4.*s**3*m*beta*cost**2*sintp - 15./4.*s**3*m*
     & beta*cost**3*sintp - 15./8.*s**3*m*beta**2*cost**3*sintp + 15./8.
     & *s**3*m*beta**2*cost**4*sintp + 9./8.*s**4*m**(-1)*cost*sintp -
     & 9./8.*s**4*m**(-1)*cost**2*sintp - 9./8.*s**4*m**(-1)*beta**2*
     & cost*sintp + 9./8.*s**4*m**(-1)*beta**2*cost**2*sintp - 9./8.*
     & s**4*m**(-1)*beta**2*cost**3*sintp + 9./8.*s**4*m**(-1)*beta**2*
     & cost**4*sintp + 9./8.*s**4*m**(-1)*beta**4*cost**3*sintp - 9./8.
     & *s**4*m**(-1)*beta**4*cost**4*sintp

      mptt =
     &  - 9*s**3*m**(-1)*beta**(-2)*sintp**3 + 9*s**3*m**(-1)*sintp**3
     &  + 9*s**3*m**(-1)*cost**2*sintp**3 - 9*s**3*m**(-1)*beta**2*
     & cost**2*sintp**3 - 7./2.*s**3*m*beta**(-1)*sintp + 7./2.*s**3*m*
     & beta**(-1)*cost*sintp - 17./8.*s**3*m*cost*sintp + 17./8.*s**3*m
     & *cost**2*sintp - 1./2.*s**3*m*beta*sintp + 1./2.*s**3*m*beta*
     & cost*sintp + 15./4.*s**3*m*beta*cost**2*sintp - 15./4.*s**3*m*
     & beta*cost**3*sintp + 15./8.*s**3*m*beta**2*cost**3*sintp - 15./8.
     & *s**3*m*beta**2*cost**4*sintp - 9./8.*s**4*m**(-1)*cost*sintp +
     & 9./8.*s**4*m**(-1)*cost**2*sintp + 9./8.*s**4*m**(-1)*beta**2*
     & cost*sintp - 9./8.*s**4*m**(-1)*beta**2*cost**2*sintp + 9./8.*
     & s**4*m**(-1)*beta**2*cost**3*sintp - 9./8.*s**4*m**(-1)*beta**2*
     & cost**4*sintp - 9./8.*s**4*m**(-1)*beta**4*cost**3*sintp + 9./8.
     & *s**4*m**(-1)*beta**4*cost**4*sintp

cccccc

      mm=mmtt/dent+mmuu/denu
      pp=pptt/dent+ppuu/denu
      pm=pmtt/dent+pmuu/denu
      mp=mptt/dent+mpuu/denu


ccc   Get relative phases right

      pm=pm*(cphi+zi*sphi)**2
      mp=mp*(cphi-zi*sphi)**2

ccc

      pm=pm*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      mp=mp*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      pp=pp*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m
      mm=mm*64d0*alphas(mu**2/4d0)**2*pi*psi0**2/9d0/m

      pm=pm*2d0
      mp=mp*2d0
      pp=pp*2d0
      mm=mm*2d0

      endif

      return
      end

