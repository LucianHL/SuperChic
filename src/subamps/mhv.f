ccc   MHV amplitudes
      subroutine denmhvc(p1,p2,p3,p4,p5,zout)
      implicit none
      complex*16 zout,zouta,zoutb,zoutc,zoutd
     &     ,zoute,zoutf
      double precision p1(4),p2(4),p3(4),p4(4),p5(4)

      call zmhvg(p1,p3,p4,p5,p2,zouta)
      call zmhvg(p1,p3,p5,p4,p2,zoutb)
      call zmhvg(p1,p4,p3,p5,p2,zoutc)
      call zmhvg(p1,p4,p5,p3,p2,zoutd)
      call zmhvg(p1,p5,p3,p4,p2,zoute)
      call zmhvg(p1,p5,p4,p3,p2,zoutf)

      zout=zouta-zoutb-zoutc+zoutd+zoute-zoutf
      zout=zout*6d0

      return
      end

      subroutine denmhvmc(p1,p2,p3,p4,p5,zout)
      implicit none
      complex*16 zout,zouta,zoutb,zoutc,zoutd
     &     ,zoute,zoutf
      double precision p1(4),p2(4),p3(4),p4(4),p5(4)

      call zmhvgm(p1,p3,p4,p5,p2,zouta)
      call zmhvgm(p1,p3,p5,p4,p2,zoutb)
      call zmhvgm(p1,p4,p3,p5,p2,zoutc)
      call zmhvgm(p1,p4,p5,p3,p2,zoutd)
      call zmhvgm(p1,p5,p3,p4,p2,zoute)
      call zmhvgm(p1,p5,p4,p3,p2,zoutf)

      zout=zouta-zoutb-zoutc+zoutd+zoute-zoutf
      zout=zout*6d0

      return
      end

      subroutine denmhvcq(p1,p2,p3,p4,p5,zout)
      implicit none
      complex*16 zout,zouta,zoutb,zoutc,zoutd
     &     ,zoute,zoutf
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),c1,c2

      call zmhvg(p3,p1,p2,p5,p4,zouta)
      call zmhvg(p3,p1,p5,p2,p4,zoutb)
      call zmhvg(p3,p2,p1,p5,p4,zoutc)
      call zmhvg(p3,p2,p5,p1,p4,zoutd)
      call zmhvg(p3,p5,p1,p2,p4,zoute)
      call zmhvg(p3,p5,p2,p1,p4,zoutf)

      c1=4d0/3d0
      c2=-1d0/6d0

      zout=c1*(zouta+zoutc+zoute+zoutf)
     &+c2*(zoutb+zoutd)

      return
      end

      subroutine denmhvmcq(p1,p2,p3,p4,p5,zout)
      implicit none
      complex*16 zout,zouta,zoutb,zoutc,zoutd
     &     ,zoute,zoutf
      double precision p1(4),p2(4),p3(4),p4(4),p5(4),c1,c2

      call zmhvgm(p3,p1,p2,p5,p4,zouta)
      call zmhvgm(p3,p1,p5,p2,p4,zoutb)
      call zmhvgm(p3,p2,p1,p5,p4,zoutc)
      call zmhvgm(p3,p2,p5,p1,p4,zoutd)
      call zmhvgm(p3,p5,p1,p2,p4,zoute)
      call zmhvgm(p3,p5,p2,p1,p4,zoutf)

      c1=4d0/3d0
      c2=-1d0/6d0

      zout=c1*(zouta+zoutc+zoute+zoutf)
     &+c2*(zoutb+zoutd)

      return
      end

      subroutine zmhvg(p1,p2,p3,p4,p5,zout)
      implicit none
      complex*16 zout,zout12,zout23,zout34,zout45,zout51
      double precision p1(4),p2(4),p3(4),p4(4),p5(4)

      call prodp(p1,p2,zout12)
      call prodp(p2,p3,zout23)
      call prodp(p3,p4,zout34)
      call prodp(p4,p5,zout45)
      call prodp(p5,p1,zout51)

ccccccccccccc

      zout=1d0/(zout12*zout23*zout34*zout45*zout51)

      return
      end

      subroutine zmhvgm(p1,p2,p3,p4,p5,zout)
      implicit none
      complex*16 zout,zout12,zout23,zout34,zout45,zout51
      double precision p1(4),p2(4),p3(4),p4(4),p5(4)

      call prodm(p1,p2,zout12)
      call prodm(p2,p3,zout23)
      call prodm(p3,p4,zout34)
      call prodm(p4,p5,zout45)
      call prodm(p5,p1,zout51)

ccccccccccccc

      zout=1d0/(zout12*zout23*zout34*zout45*zout51)

      return
      end

      subroutine prodm(pa,pb,zout)
      implicit none
      double precision pa(4),pb(4),sdot
      complex*16 zout,zout1

      call prodp(pa,pb,zout1)
      zout=2d0*sdot(pa,pb)/zout1

      return
      end

      subroutine prodp(pa,pb,zout)
      implicit none
      double precision pa(4),pb(4),cphi,sphi
      double precision pap,pbp,pam,pbm,sdot
      complex*16 zout

      include 'zi.f'

      pap=pa(4)+pa(1)
      pbp=pb(4)+pb(1)
      pam=pa(4)-pa(1)
      pbm=pb(4)-pb(1)

      cphi=(pa(3)*pbp-pb(3)*pap)/dsqrt(pap*pbp*dabs(sdot(pa,pb))*2d0)
      sphi=(pa(2)*pbp-pb(2)*pap)/dsqrt(pap*pbp*dabs(sdot(pa,pb))*2d0)

      zout=dsqrt(2d0*dabs(sdot(pa,pb)))*(cphi+zi*sphi)

      return
      end
