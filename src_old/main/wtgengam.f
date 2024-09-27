ccc   calls subprocess amplitude
      subroutine wtgengam
      implicit none
      integer p
      complex*16 pp,mm,pm,mp

      include 'polarization.f'
      include 'proc.f'
      include 'ppamp.f'
      include 'mandelstam.f'
      include 'vars.f'
c      include 'zoutarr.f'
      include 'diss.f'
      include 'eff.f'


      if(offshell)then

         do p=1,pol
            if(proc.eq.54.or.proc.eq.55)then
c               if(p.eq.1.or.p.eq.2.or.p.eq.5.or.p.eq.7.or.p.eq.9)then
                  call wwoff_axial(p)
c               endif
            endif
            if(proc.eq.56.or.proc.eq.57.or.proc.eq.58)then
               call lloff(p)
            endif
            if(proc.eq.68)then
               call alp_off(p)
            endif
         enddo
      else

      if(proc.eq.54.or.proc.eq.55)then
         do p=1,pol
            call wwpol(p,uh,th,pp,mm,pm,mp)
            ppa(p)=pp
            mma(p)=mm
            pma(p)=pm
            mpa(p)=mp
            pincarr(p)=cdabs(pp)**2+cdabs(mm)**2
     &           +cdabs(pm)**2+cdabs(mp)**2
         enddo
      elseif(proc.eq.56.or.proc.eq.57.or.proc.eq.58.or.proc.eq.61
     &   .or.proc.eq.73.or.proc.eq.74.or.proc.eq.75.or.
     &        proc.eq.76)then
         do p=1,pol
            call llpol(p,mx,uh,th,pp,mm,pm,mp)
            ppa(p)=pp
            mma(p)=mm
            pma(p)=pm
            mpa(p)=mp

c            print*,p
c            print*,'pp,mm,pm,mp=',pp,mm,pm,mp
c            print*,''

            pincarr(p)=cdabs(pp)**2+cdabs(mm)**2
     &           +cdabs(pm)**2+cdabs(mp)**2
         enddo
      elseif(proc.eq.59)then
         do p=1,pol
            call lightlightpol(p,mx,uh,th,pp,mm,pm,mp)
            ppa(p)=pp
            mma(p)=mm
            pma(p)=pm
            mpa(p)=mp
            pincarr(p)=cdabs(pp)**2+cdabs(mm)**2
     &           +cdabs(pm)**2+cdabs(mp)**2
         enddo
      elseif(proc.eq.60)then
          call higgsgam(mx,pp,mm,pm,mp)
            ppa(1)=pp
            mma(1)=mm
            pma(1)=pm
            mpa(1)=mp
            pincarr(1)=cdabs(pp)**2+cdabs(mm)**2
     &           +cdabs(pm)**2+cdabs(mp)**2
      elseif(proc.eq.68)then
         call alp(mx,pp,mm,pm,mp)
         ppa(1)=pp
         mma(1)=mm
         pma(1)=pm
         mpa(1)=mp
         pincarr(1)=cdabs(pp)**2+cdabs(mm)**2
     &        +cdabs(pm)**2+cdabs(mp)**2
      elseif(proc.eq.69.or.proc.eq.70)then
         call monop(mx,pp,mm,pm,mp)
         ppa(1)=pp
         mma(1)=mm
         pma(1)=pm
         mpa(1)=mp
         pincarr(1)=cdabs(pp)**2+cdabs(mm)**2
     &        +cdabs(pm)**2+cdabs(mp)**2
      elseif(proc.eq.71.or.proc.eq.72)then
         do p=1,pol
            call mmpol(p,mx,uh,th,pp,mm,pm,mp)
            ppa(p)=pp
            mma(p)=mm
            pma(p)=pm
            mpa(p)=mp
            pincarr(p)=cdabs(pp)**2+cdabs(mm)**2
     &           +cdabs(pm)**2+cdabs(mp)**2
         enddo
      elseif(proc.eq.82.or.proc.eq.83.or.proc.eq.84)then
         call sim_modVX(mx,pp,mm,pm,mp)
         ppa(1)=pp
         mma(1)=mm
         pma(1)=pm
         mpa(1)=mp
         pincarr(1)=cdabs(pp)**2+cdabs(mm)**2
     &        +cdabs(pm)**2+cdabs(mp)**2
      endif

      endif

      return
      end
