ccc   gamma gamma --> W+W- subprocess amplitude
      subroutine wwpol(p,u,t,pp,mm,pm,mp)
      implicit none
      double precision sphi,cphi,sintt,normp,m,costt
      double precision lam1,lam2,lam3,lam4,gam,u,t,beta
      complex*16 pp,mm,pm,mp
      integer p,l1,l2

      include 'ewpars.f'
      include 'pi.f'
      include 'scorr.f'
      include 'partonmom2.f'
      include 'zi.f'
      include 'vars.f'
      include 'norm.f'
      include 'mom.f'

      cphi=p1(1)/dsqrt(p1(1)**2+p1(2)**2)
      sphi=p1(2)/dsqrt(p1(1)**2+p1(2)**2)

      beta=dsqrt(1d0-4d0*mw**2/mx**2)
      gam=mx**2/mw**2
    
      costt=(t-u)/beta/mx**2
      sintt=dsqrt(1d0-costt**2)

      normp=16d0*pi**2*alpha**2/(1d0-beta**2*costt**2)**2
      normp=normp/32d0/pi/mx**2*beta/4d0
      normp=normp*2d0
      normp=normp*conv
      normp=dsqrt(normp)

      if(p.eq.1)then ! ++
         lam3=1d0
         lam4=1d0
      elseif(p.eq.2)then ! +-
         lam3=1d0
         lam4=-1d0
      elseif(p.eq.3)then ! -+
         lam3=-1d0
         lam4=1d0
      elseif(p.eq.4)then ! --
         lam3=-1d0
         lam4=-1d0
      elseif(p.eq.5)then ! 0+
         lam4=1d0
      elseif(p.eq.6)then ! 0-
         lam4=-1d0
      elseif(p.eq.7)then ! +0
         lam3=1d0
      elseif(p.eq.8)then ! -0
         lam3=-1d0
      endif

      do l1=-1,1,2
         do l2=-1,1,2

            lam1=dble(l1)
            lam2=dble(l2)

            if(p.lt.5)then
               m=beta*(lam1+lam2)*(lam3+lam4)
     &              +(-8d0*lam1*lam2*(1d0+lam3*lam4)
     &              +gam*(1d0+lam1*lam2*lam3*lam4)*(3d0+lam1*lam2)
     &              +2d0*gam*(lam1-lam2)*(lam3-lam4)*costt
     &              -4d0*(1d0-lam1*lam2)*(1d0+lam3*lam4)*costt**2
     &              +gam*(1d0-lam1*lam2)*(1d0-lam3*lam4)*costt**2)
     &              /2d0/gam
            elseif(p.lt.7)then
               m=(lam1-lam2)*(1d0-lam2*lam4*costt)*sintt*dsqrt(8d0/gam)
            elseif(p.lt.9)then
               m=(lam1-lam2)*(1d0+lam1*lam3*costt)*sintt*dsqrt(8d0/gam)
            elseif(p.eq.9)then
               m=(-8d0+(1d0-lam1*lam2)*(4d0+(4d0+gam)*sintt**2))/gam
            endif

            m=m*normp

            if(l1.eq.-1.and.l2.eq.-1)then
               mm=m
            elseif(l1.eq.1.and.l2.eq.1)then
               pp=m
            elseif(l1.eq.1.and.l2.eq.-1)then
               pm=m
            elseif(l1.eq.-1.and.l2.eq.1)then
               mp=m
            endif
            
         enddo
      enddo

      pm=pm*(cphi-zi*sphi)**2
      mp=mp*(cphi+zi*sphi)**2

      return
      end
