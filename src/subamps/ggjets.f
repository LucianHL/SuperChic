ccc   gg --> gg subprocess amplitude
      subroutine gg(p,mx,u,t,pp,mm,pm,mp)
      implicit none
      complex*16 pp,mm,pm,mp
      double precision sphi,cphi,u,t,mx,alphas
      integer p

      include 'zi.f'
      include 'pi.f'
      include 'partonmom2.f'

      cphi=p1(1)/dsqrt(p1(1)**2+p1(2)**2)
      sphi=p1(2)/dsqrt(p1(1)**2+p1(2)**2)
      
      if(p.eq.1)then            ! ++ final state
         pp=3d0/dsqrt(8d0)*8d0*pi*alphas(mx**2)
         pp=pp*mx**4/u/t
         mm=0d0
         mp=0d0
         pm=0d0
      elseif(p.eq.2)then        ! -- final state
         pp=0d0
         mm=3d0/dsqrt(8d0)*8d0*pi*alphas(mx**2)
         mm=mm*mx**4/u/t
         mp=0d0
         pm=0d0
      elseif(p.eq.3)then        ! +- final state
         pp=0d0
         mm=0d0
         pm=3d0/dsqrt(8d0)*8d0*pi*alphas(mx**2)
         mp=pm*t/u
         pm=pm*u/t
      elseif(p.eq.4)then        ! -+ final state
         pp=0d0
         mm=0d0
         pm=3d0/dsqrt(8d0)*8d0*pi*alphas(mx**2)
         mp=pm*u/t
         pm=pm*t/u
      endif

      pm=pm*(cphi-zi*sphi)**2
      mp=mp*(cphi+zi*sphi)**2

      return
      end







