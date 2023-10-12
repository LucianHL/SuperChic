ccc   gg --> etaeta' subprocess amplitude
      subroutine etaetap(p,mx,u,t,pp,mm,pm,mp)
      implicit none
      double precision phi,cphi,sphi
      double precision u,t,mx
      double precision ppqqs,ppqg,ppgg,pmqqs,pmqqns
      complex*16 pp,mm,pm,mp
      integer p

      include 'mixing.f'
      include 'partonmom2.f'
      include 'zi.f'

cccccc

      phi=datan(p1(2)/p1(1))
      cphi=p1(1)/dsqrt(p1(1)**2+p1(2)**2)
      sphi=p1(2)/dsqrt(p1(1)**2+p1(2)**2)

cccccc

      call eta(p,mx,u,t,ppqqs,ppqg,ppgg,pmqqs,pmqqns)


      pp=-dsin(thetap1)*dcos(thetap1)*(ppqqs+ppqg+ppgg)*fpi1**2
      mm=pp

      pm=-fpi1**2*dsin(thetap1)*dcos(thetap1)*pmqqs
     &+(fpi8**2*dcos(thetap8)*dsin(thetap1)
     &-fpi1**2*dsin(thetap1)*dcos(thetap8))*pmqqns
      mp=pm

      pm=pm*(cphi+zi*sphi)**2
      mp=mp*(cphi-zi*sphi)**2

      return
      end
