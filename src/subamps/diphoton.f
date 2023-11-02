      subroutine gamgam(p,mu,u,t,pp,mm,pm,mp)
      implicit none
      complex*16 pp,mm,pm,mp
      integer i,p
      complex*16 b0fqsff,b0ftsff,b0fusff
      complex*16 c0ccqsfff,c0cctsfff,c0ccusfff
      complex*16 d0ccccqstsffff,d0ccccqsusffff,d0ccccustsffff
      complex*16 wm2
      double precision stw,sh,rwm2,norm,mz,mgen2,ctw,cqfactor
      double precision u,t,mu,alphas
      double precision qfarr(12)

      include 'lbylamps.f'
      include 'lbylfac.f'
      include 'lep.f'
      include 'pvfun.f'

      include 'pi.f'
      include 'zi.f'
      include 'vars.f'
      include 'ewpars.f'
      include 'norm.f'


      mz=91.1876d0
      stw=0.48076d0
      ctw=dsqrt(1d0-stw**2)
      mw=mz*ctw
      rwm2=mw**2
      wm2 = dcmplx(rwm2,-1d-30)

cccccccccc

      norm=4d0*alpha*alphas(mu**2/4d0)

      sh=mx**2

ccccccccccccccc

      qfarr( 4) = -1d0/3d0
      qfarr( 5) =  2d0/3d0
      qfarr( 6) = -1d0/3d0
      qfarr( 7) =  2d0/3d0
      qfarr( 8) = -1d0/3d0
      qfarr( 9) =  2d0/3d0

      do i=4,9  ! No loop over leptons for gg (1-3)
         mars(i)=rmf1(i)**2
         marsi(i)=dcmplx(rmf1(i)**2,-1d-30)
         cqfa(i)=qfarr(i)**2   ! Changed from **4 as gg i.s.
c         if(i.gt.3)cqfa(i)=cqfa(i)*3d0  ! Does not apply for gg
      enddo

      if(p.eq.1)then

      qedamp(1,1,1,1) = 0d0
      qedamp(1,1,1,2) = 0d0
      qedamp(1,1,2,2) = 0d0
      qedamp(1,2,1,2) = 0d0
      qedamp(1,2,2,1) = 0d0

      do i=4,9                  ! number of fermion loop particles

      mgen2 = mars(i)
      mgen2i = marsi(i)
      cqfactor = cqfa(i)

c       mgen2=0d0
c      mgen2i=dcmplx(0d0,-1d-30)

      if (mgen2.eq.0) then

      qedamp(1,1,2,2) = qedamp(1,1,2,2) + cqfactor

      qedamp(1,1,1,2) = qedamp(1,1,1,2) + cqfactor

      qedamp(1,1,1,1) = qedamp(1,1,1,1) + cqfactor*(
     d -1d0 + (t-u)/sh*(LOG(-u/sh) - LOG(-t/sh))
     & -(1d0/2 - u*t/sh**2)*((LOG(-u/sh)-LOG(-t/sh))**2 + pi**2))

      qedamp(1,2,2,1) = qedamp(1,2,2,1) + cqfactor*(
     & -1d0 - zi*pi*(t-sh)/u
     & -((1d0 + zi*pi)*(t-sh)/u + 2d0*zi*pi*(t/u)**2)*LOG(-t/sh)
     & -(1d0/2 - sh*t/u**2)*LOG(-t/sh)**2)

      qedamp(1,2,1,2) = qedamp(1,2,1,2) + cqfactor*(
     & -1d0 - zi*pi*(u-sh)/t
     & -((1d0 + zi*pi)*(u-sh)/t + 2d0*zi*pi*(u/t)**2)*LOG(-u/sh)
     & -(1d0/2 - sh*u/t**2)*LOG(-u/sh)**2)

      else

      b0fqsff = b0f2m(sh,mgen2i)
      b0ftsff = b0f2m(t,mgen2i)
      b0fusff = b0f2m(u,mgen2i)


      c0ccqsfff = -c01(0d0,0d0,-sh,mgen2i,mgen2i,mgen2i)
      c0cctsfff = -c01(0d0,0d0,-t,mgen2i,mgen2i,mgen2i)
      c0ccusfff = -c01(0d0,0d0,-u,mgen2i,mgen2i,mgen2i)
      d0ccccqstsffff = d0404M(sh,t,mgen2i)
      d0ccccqsusffff = d0404M(sh,u,mgen2i)
      d0ccccustsffff = d0404M(u,t,mgen2i)

      qedamp(1,1,2,2) = qedamp(1,1,2,2) + cqfactor*(
     & +1d0-2d0*mgen2**2*(d0ccccqstsffff+d0ccccqsusffff+d0ccccustsffff))

      qedamp(1,1,1,2) = qedamp(1,1,1,2) + cqfactor*(
     & +1d0-mgen2*(sh**2+t**2+u**2)*(c0ccqsfff/u/t
     &                              +c0cctsfff/sh/u+c0ccusfff/sh/t)
     & -mgen2*((2d0*mgen2+sh*t/u)*d0ccccqstsffff
     &         +(2d0*mgen2+sh*u/t)*d0ccccqsusffff
     &         +(2d0*mgen2+t*u/sh)*d0ccccustsffff))

      qedamp(1,1,1,1) = qedamp(1,1,1,1) + cqfactor*(
     & -1d0+(u-t)/sh*(b0fusff-b0ftsff)
     & +(4d0*mgen2/sh+2d0*(t*u/sh**2-1d0/2d0))*(u*c0ccusfff+t*c0cctsfff)
     & -2d0*mgen2*sh*(mgen2/sh-1d0/2d0)*(d0ccccqstsffff
     &                               +d0ccccqsusffff+d0ccccustsffff)
     & -t*u*(4d0*mgen2/sh+t*u/sh**2-1d0/2d0)*d0ccccustsffff)


      qedamp(1,2,2,1) = qedamp(1,2,2,1) + cqfactor*(
     & -1d0+(sh-t)/u*(b0fqsff-b0ftsff)
     & +(4d0*mgen2/u+2d0*(t*sh/u**2-1d0/2d0))*(sh*c0ccqsfff+t*c0cctsfff)
     & -2d0*mgen2*u*(mgen2/u-1d0/2d0)*(d0ccccqstsffff
     &                        +d0ccccqsusffff+d0ccccustsffff)
     & -t*sh*(4d0*mgen2/u+t*sh/u**2-1d0/2d0)*d0ccccqstsffff)

      qedamp(1,2,1,2) = qedamp(1,2,1,2) + cqfactor*(
     & -1d0+(u-sh)/t*(b0fusff-b0fqsff)
     & +(4d0*mgen2/t+2d0*(sh*u/t**2-1d0/2d0))*(u*c0ccusfff+sh*c0ccqsfff)
     & -2d0*mgen2*t*(mgen2/t-1d0/2d0)*(d0ccccqstsffff
     &                        +d0ccccqsusffff+d0ccccustsffff)
     & -sh*u*(4d0*mgen2/t+sh*u/t**2-1d0/2d0)*d0ccccqsusffff)

      endif
      enddo

      qedamp(1,1,2,1) = qedamp(1,1,1,2)
      qedamp(1,2,1,1) = qedamp(1,1,1,2)
      qedamp(1,2,2,2) = qedamp(1,1,1,2)
      qedamp(2,1,1,1) = qedamp(1,2,2,2)
      qedamp(2,1,1,2) = qedamp(1,2,2,1)
      qedamp(2,1,2,1) = qedamp(1,2,1,2)
      qedamp(2,1,2,2) = qedamp(1,2,1,1)
      qedamp(2,2,1,1) = qedamp(1,1,2,2)
      qedamp(2,2,1,2) = qedamp(1,1,2,1)
      qedamp(2,2,2,1) = qedamp(1,1,1,2)
      qedamp(2,2,2,2) = qedamp(1,1,1,1)

      endif

      if(p.eq.1)then
         pp=qedamp(1,1,1,1)
         mm=qedamp(1,1,2,2)
         pm=qedamp(1,1,1,2)
         mp=qedamp(1,1,2,1)
      elseif(p.eq.2)then
         pp=qedamp(2,2,1,1)
         mm=qedamp(2,2,2,2)
         pm=qedamp(2,2,1,2)
         mp=qedamp(2,2,2,1)
      elseif(p.eq.3)then
         pp=qedamp(1,2,1,1)
         mm=qedamp(1,2,2,2)
         pm=qedamp(1,2,1,2)
         mp=qedamp(1,2,2,1)
      elseif(p.eq.4)then
         pp=qedamp(2,1,1,1)
         mm=qedamp(2,1,2,2)
         pm=qedamp(2,1,1,2)
         mp=qedamp(2,1,2,1)
      endif

      pp=pp*norm
      mm=mm*norm
      pm=pm*norm
      mp=mp*norm

      return
      end
