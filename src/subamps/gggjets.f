ccc   gg --> ggg subprocess amplitude
      subroutine ggg(p,mx,pp,mm,pm,mp)
      implicit none
      complex*16 zout34,zout15,zout25,zout45,zout35
     &,zout24,zout23,zout13,zout14,zout12
      complex*16 pp,mm,pm,mp
      double precision znorm,mx,alphas
      integer p

      include 'pi.f'
      include 'zi.f'
      include 'mhvpars.f'
      include 'partonmom3.f'
      include 'ewpars.f'
      
      znorm=dsqrt(24d0)/8d0/8d0
      znorm=znorm*(4d0*dsqrt(2d0))
      znorm=znorm*(4d0*pi*alphas(mx**2))**(1.5d0)

      if(p.eq.1)then
         call denmhvc(pa,pb,p3,p4,p5,zdenp)
         call denmhvmc(pa,pb,p3,p4,p5,zdenm)
      endif

      if(p.eq.1)then            ! +++ fs
         call prodp(pa,pb,zout12)
         
         mm=0d0
         pm=0d0
         mp=0d0
         pp=zdenp*zout12**4

      elseif(p.eq.2)then        ! ++- fs

         call prodp(p3,p4,zout34)
         call prodm(pa,p5,zout15)
         call prodm(pb,p5,zout25)
    
         
         pp=zdenp*zout34**4
         mm=0d0
         pm=zdenm*zout15**4
         mp=zdenm*zout25**4

      elseif(p.eq.3)then        ! +-+ fs
         call prodp(p3,p5,zout35)
         call prodm(pa,p4,zout14)
         call prodm(pb,p4,zout24)

         pp=zdenp*zout35**4
         mm=0d0
         pm=zdenm*zout14**4
         mp=zdenm*zout24**4

      elseif(p.eq.4)then        ! -++ fs
         call prodp(p4,p5,zout45)
         call prodm(pa,p3,zout13)
         call prodm(pb,p3,zout23)
 
         pp=zdenp*zout45**4
         mm=0d0
         pm=zdenm*zout13**4
         mp=zdenm*zout23**4

      elseif(p.eq.5)then        ! --- fs
         call prodm(pa,pb,zout12)
         
         pp=0d0
         mm=zdenm*zout12**4
         pm=0d0
         mp=0d0

      elseif(p.eq.6)then        ! --+ fs
         call prodm(p3,p4,zout34)
         call prodp(pa,p5,zout15)
         call prodp(pb,p5,zout25)

         pp=0d0
         mm=zdenm*zout34**4
         pm=zdenp*zout25**4
         mp=zdenp*zout15**4

      elseif(p.eq.7)then        ! -+- fs
         call prodm(p3,p5,zout35)
         call prodp(pa,p4,zout14)
         call prodp(pb,p4,zout24)

         pp=0d0
         mm=zdenm*zout35**4
         pm=zdenp*zout24**4
         mp=zdenp*zout14**4

      elseif(p.eq.8)then       ! +-- fs
         call prodm(p4,p5,zout45)
         call prodp(pa,p3,zout13)
         call prodp(pb,p3,zout23)

         pp=0d0
         mm=zdenm*zout45**4
         pm=zdenp*zout23**4
         mp=zdenp*zout13**4

      endif

      pp=pp*znorm
      mm=mm*znorm
      pm=pm*znorm
      mp=mp*znorm

      return
      end




