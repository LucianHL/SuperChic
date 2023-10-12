ccc   gg --> qqbarg subprocess amplitude
      subroutine qqg(p,mx,pp,mm,pm,mp)
      implicit none
      complex*16 pp,mm,pm,mp
      complex*16 zpp,zpm,zout,zmp,zmm
      complex*16 zout45,zout35,zout13,zout14,zout23,zout24
      double precision mx,alphas,normp
      integer p
   
      include 'partonmom3.f'
      include 'pi.f'
      include 'zi.f'
      include 'mhvpars.f'

        if(p.eq.1)then
           call denmhvcq(pa,pb,p3,p4,p5,zdenp)
           call denmhvmcq(pa,pb,p3,p4,p5,zdenm)
        endif

        if(p.eq.1)then          ! g(+)q(-)qbar(+) final state

         call prodp(p4,p5,zout45)
         call prodp(p3,p5,zout35)
         call prodm(pa,p3,zout13)
         call prodm(pa,p4,zout14)
         call prodm(pb,p3,zout23)
         call prodm(pb,p4,zout24)

         zpp=zout45**3*zout35*zdenp
         zmm=0d0
         zpm=zout13**3*zout14*zdenm         
         zmp=zout23**3*zout24*zdenm     

        elseif(p.eq.2)then      ! g(+)q(+)qbar(-) final state

         call prodp(p4,p5,zout45)
         call prodp(p3,p5,zout35)
         call prodm(pa,p3,zout13)
         call prodm(pa,p4,zout14)
         call prodm(pb,p3,zout23)
         call prodm(pb,p4,zout24)

         zpp=zout45*zout35**3*zdenp
         zmm=0d0
         zpm=zout13*zout14**3*zdenm         
         zmp=zout23*zout24**3*zdenm    

        elseif(p.eq.3)then      ! g(-)q(-)qbar(+) final state

        call prodm(p3,p5,zout35)
        call prodm(p4,p5,zout45)
        call prodp(pb,p3,zout23)
        call prodp(pb,p4,zout24)
        call prodp(pa,p3,zout13)
        call prodp(pa,p4,zout14)
  
        zpp=0d0
        zmm=zout45*zout35**3*zdenm
        zpm=zout23*zout24**3*zdenp
        zmp=zout13*zout14**3*zdenp

        elseif(p.eq.4)then      ! g(-)q(+)qbar(-) final state

        call prodm(p3,p5,zout35)
        call prodm(p4,p5,zout45)
        call prodp(pb,p3,zout23)
        call prodp(pb,p4,zout24)
        call prodp(pa,p3,zout13)
        call prodp(pa,p4,zout14)

        zpp=0d0
        zmm=zout45**3*zout35*zdenm
        zpm=zout23**3*zout24*zdenp
        zmp=zout13**3*zout14*zdenp

        endif

        normp=(4d0*pi*alphas(mx**2))**(1.5d0)/dsqrt(2d0)
        pp=zpp*normp
        mm=zmm*normp
        pm=zpm*normp
        mp=zmp*normp

      return
      end



