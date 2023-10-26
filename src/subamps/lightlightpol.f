      subroutine lightlightpol(p,mu,u,t,pp,mm,pm,mp)
      implicit none
      complex*16 pp,mm,pm,mp
      integer i,p
      complex*16 b0fqsww,b0ftsww,b0fusww
      complex*16 c0ccqswww,c0cctswww,c0ccuswww
      complex*16 d0ccccqstswwww,d0ccccqsuswwww,d0ccccustswwww
      complex*16 b0fqsff,b0ftsff,b0fusff
      complex*16 c0ccqsfff,c0cctsfff,c0ccusfff
      complex*16 d0ccccqstsffff,d0ccccqsusffff,d0ccccustsffff
      complex*16 wm2
      double precision stw,sh,rwm2,ql,qf,norm,mz,mgen2,efac
      double precision ctw,cqfactor,u,t,mu
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

      if(mu.gt.2d0*mb)then
         qf=35d0/81d0
      else
         qf=34d0/81d0
      endif

      qf=qf*3d0
      ql=3d0
      
      norm=(8d0*alpha**2)**2
      norm=norm/16d0/pi/mx**2
      norm=norm/4d0
      norm=norm*conv
      
      norm=dsqrt(norm)

      sh=mx**2
          
ccccccccccccccc
      
      qfarr( 1) = -1d0
      qfarr( 2) = -1d0
      qfarr( 3) = -1d0
      qfarr( 4) = -1d0/3d0
      qfarr( 5) =  2d0/3d0
      qfarr( 6) = -1d0/3d0
      qfarr( 7) =  2d0/3d0
      qfarr( 8) = -1d0/3d0
      qfarr( 9) =  2d0/3d0
 
      do i=1,9
         mars(i)=rmf1(i)**2
         marsi(i)=dcmplx(rmf1(i)**2,-1d-30)
         cqfa(i)=qfarr(i)**4
         if(i.gt.3)cqfa(i)=cqfa(i)*3d0
      enddo

      if(p.eq.1)then

      qedamp(1,1,1,1) = 0d0
      qedamp(1,1,1,2) = 0d0
      qedamp(1,1,2,2) = 0d0
      qedamp(1,2,1,2) = 0d0
      qedamp(1,2,2,1) = 0d0
         
      do i=1,9                  ! number of fermion loop particles
         
      mgen2 = mars(i)
      mgen2i = marsi(i)
      cqfactor = cqfa(i)

c      mgen2=0d0
c      mgen2i=dcmplx(0d0,-1d-30)
      
       if (ABS(mgen2).LT.1.D-12) then
         
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


cccccccc

      
      
      ewamp(1,1,1,1) = 0d0
      ewamp(1,1,1,2) = 0d0
      ewamp(1,1,2,2) = 0d0
      ewamp(1,2,1,2) = 0d0
      ewamp(1,2,2,1) = 0d0

      b0fqsww = b0f2m(sh,wm2)
      b0ftsww = b0f2m(t,wm2)
      b0fusww = b0f2m(u,wm2)
      c0ccqswww = -c01(0d0,0d0,-sh,wm2,wm2,wm2)
      c0cctswww = -c01(0d0,0d0,-t,wm2,wm2,wm2)
      c0ccuswww = -c01(0d0,0d0,-u,wm2,wm2,wm2)
      d0ccccqstswwww = d0404M(sh,t,wm2)
      d0ccccqsuswwww = d0404M(sh,u,wm2)
      d0ccccustswwww = d0404M(u,t,wm2)
 
      ewamp(1,1,2,2) = 
     & +1d0-2d0*rwm2**2*(d0ccccqstswwww+d0ccccqsuswwww+d0ccccustswwww)

      ewamp(1,1,1,2) = 
     & +1d0-rwm2*(sh**2+t**2+u**2)*( 1d0/u/t*c0ccqswww
     &                             +1d0/sh/u*c0cctswww
     &                             +1d0/sh/t*c0ccuswww)
     & -rwm2*((2d0*rwm2+sh*t/u)*d0ccccqstswwww
     &       +(2d0*rwm2+sh*u/t)*d0ccccqsuswwww
     &       +(2d0*rwm2+t*u/sh)*d0ccccustswwww)

      ewamp(1,1,1,1) = 
     & -1d0+(u-t)/sh*(b0fusww-b0ftsww)
     & +(4d0*rwm2/sh+2d0*(t*u/sh**2-4d0/3d0))*(u*c0ccuswww+t*c0cctswww)
     & -(2d0*rwm2*sh*(rwm2/sh-4d0/3d0)+2d0/3d0*sh**2)*
     &                        ( d0ccccqstswwww
     &                         +d0ccccqsuswwww
     &                         +d0ccccustswwww )
     & -t*u*(4d0*rwm2/sh+t*u/sh**2-4d0/3d0)*d0ccccustswwww

      ewamp(1,2,2,1) = 
     & -1d0+(sh-t)/u*(b0fqsww-b0ftsww)
     & +(4d0*rwm2/u+2d0*(t*sh/u**2-4d0/3d0))*(sh*c0ccqswww+t*c0cctswww)
     & -(2d0*rwm2*u*(rwm2/u-4d0/3d0)+2d0/3d0*u**2)*
     &                        ( d0ccccqstswwww
     &                         +d0ccccqsuswwww
     &                         +d0ccccustswwww )
     & -t*sh*(4d0*rwm2/u+t*sh/u**2-4d0/3d0)*d0ccccqstswwww

      ewamp(1,2,1,2) = 
     & -1d0+(u-sh)/t*(b0fusww-b0fqsww)
     & +(4d0*rwm2/t+2d0*(sh*u/t**2-4d0/3d0))*(u*c0ccuswww+sh*c0ccqswww)
     & -(2d0*rwm2*t*(rwm2/t-4d0/3d0)+2d0/3d0*t**2)*
     &                        ( d0ccccqstswwww
     &                         +d0ccccqsuswwww
     &                         +d0ccccustswwww )
     & -sh*u*(4d0*rwm2/t+sh*u/t**2-4d0/3d0)*d0ccccqsuswwww

      ewamp(1,1,2,1) = ewamp(1,1,1,2)
      ewamp(1,2,1,1) = ewamp(1,1,1,2)
      ewamp(1,2,2,2) = ewamp(1,1,1,2)
      ewamp(2,1,1,1) = ewamp(1,2,2,2)
      ewamp(2,1,1,2) = ewamp(1,2,2,1)
      ewamp(2,1,2,1) = ewamp(1,2,1,2)
      ewamp(2,1,2,2) = ewamp(1,2,1,1)
      ewamp(2,2,1,1) = ewamp(1,1,2,2)
      ewamp(2,2,1,2) = ewamp(1,1,2,1)
      ewamp(2,2,2,1) = ewamp(1,1,1,2)
      ewamp(2,2,2,2) = ewamp(1,1,1,1)
      
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

c$$$      pp=0d0
c$$$      mm=0d0
c$$$      pm=0d0
c$$$      mp=0d0
      
      efac=1.5d0
c      efac=0d0
      
      if(p.eq.1)then
         pp=pp-ewamp(1,1,1,1)*efac
         mm=mm-ewamp(1,1,2,2)*efac
         pm=pm-ewamp(1,1,1,2)*efac
         mp=mp-ewamp(1,1,2,1)*efac         
      elseif(p.eq.2)then
         pp=pp-ewamp(2,2,1,1)*efac
         mm=mm-ewamp(2,2,2,2)*efac
         pm=pm-ewamp(2,2,1,2)*efac
         mp=mp-ewamp(2,2,2,1)*efac
      elseif(p.eq.3)then         
         pp=pp-ewamp(1,2,1,1)*efac
         mm=mm-ewamp(1,2,2,2)*efac
         pm=pm-ewamp(1,2,1,2)*efac
         mp=mp-ewamp(1,2,2,1)*efac
      elseif(p.eq.4)then         
         pp=pp-ewamp(2,1,1,1)*efac
         mm=mm-ewamp(2,1,2,2)*efac
         pm=pm-ewamp(2,1,1,2)*efac
         mp=mp-ewamp(2,1,2,1)*efac
      endif

      pp=pp*norm
      mm=mm*norm
      pm=pm*norm
      mp=mp*norm
      
      return
      end
