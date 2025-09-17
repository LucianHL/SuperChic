ccc   re-initialise incoming nucleon momenta in nucleon-nucleon cms
      subroutine AAinit
      implicit none
      double precision beta

      include 'mom.f'
      include 'vars.f'
      include 'varsi.f'
      include 'mp.f'
      include 'mion.f'
      include 'ion.f'
      include 'sAA.f'

c      beta=dsqrt(1d0-4d0*mp**2/s)
      saa=an**2*s
      rtsaa=dsqrt(saa)

      beta=dsqrt(1d0-4d0*mion**2/saa)

      q(1,1)=0d0
      q(2,1)=0d0
      q(4,1)=rtsaa/2d0
      q(3,1)=rtsaa/2d0*beta

      q(1,2)=0d0
      q(2,2)=0d0
      q(4,2)=rtsaa/2d0
      q(3,2)=-rtsaa/2d0*beta

      return
      end

ccc   Generates correct photon kinematics for incoherent photon emission from ion
      subroutine gen_ion_inel(rarr,x1in,x2in,pt1sq,pt2sq,exitgen)
      implicit none
      double precision mdissmax,lmdissmax,lmdissmin
      double precision lmdiss1,rdiss1,wtdiss1
      double precision lmdiss2,rdiss2,wtdiss2
      double precision mpp1,mpp2
      double precision aa1,aa2,cc1,cc2,root1sq,root2sq
      double precision p1p,p2p,p1m,p2m
      double precision pt1sq,pt2sq
      double precision x1in,x2in,x1,x2
      logical exitgen
      double precision rarr(10)

      include 'vars.f'
      include 'mp.f'
      include 'diss.f'
      include 'proc.f'
      include 'mion.f'
      include 'ion.f'
      include 'mom.f'

      exitgen=.false.

      mdissmax=rts

      lmdissmax=dlog(mdissmax)
      lmdissmin=dlog(mp)

      x1=x1in
      x2=x2in

      if(diss1)then

         x1=x1in*an

         if(dps.eq.2)then
            rdiss1=rarr(7)
         else
            rdiss1=rarr(5)
         endif

         lmdiss1=lmdissmin+(lmdissmax-lmdissmin)*rdiss1
         mdiss1=dexp(lmdiss1)
         wtdiss1=2d0*mdiss1**2*(lmdissmax-lmdissmin)

      else
         mdiss1=mion
         wtdiss1=1d0
      endif

      if(diss2)then

         x2=x2in*an

         if(dps.eq.2)then
            rdiss2=rarr(7)
         else
            rdiss2=rarr(5)
         endif

         lmdiss2=lmdissmin+(lmdissmax-lmdissmin)*rdiss2
         mdiss2=dexp(lmdiss2)
         wtdiss2=2d0*mdiss2**2*(lmdissmax-lmdissmin)

      else
         mdiss2=mion
         wtdiss2=1d0
      endif

      mpp1=mdiss1
      mpp2=mdiss2

      aa1=(1d0-x1)*rts/dsqrt(2d0)
      aa2=(1d0-x2)*rts/dsqrt(2d0)
      cc1=0.5d0*(pt2sq+mpp2**2)
      cc2=0.5d0*(pt1sq+mpp1**2)

c     impose massive on-shell condition by solving
c                   p1+ + cc1/p2- = aa1
c                   p2- + cc2/p1+ = aa2

      root1sq=(cc1-cc2-aa1*aa2)**2-4d0*cc2*aa1*aa2
      root2sq=(cc2-cc1-aa1*aa2)**2-4d0*cc1*aa1*aa2

      if(root1sq.lt.0d0)then
      exitgen=.true.
      return
      endif
      if(root1sq.lt.0d0)then
      exitgen=.true.
      return
      endif

      p1p=(cc2-cc1+aa1*aa2+dsqrt(root1sq))/(2d0*aa2)
      p2m=(cc1-cc2+aa1*aa2+dsqrt(root2sq))/(2d0*aa1)
      p1m=(pt1sq+mpp1**2)/(2d0*p1p)
      p2p=(pt2sq+mpp2**2)/(2d0*p2m)

      if(p1m.lt.0d0)then
      exitgen=.true.
      return
      endif
      if(p2p.lt.0d0)then
      exitgen=.true.
      return
      endif


      print*,mdissmax
      stop

      return
      end
      