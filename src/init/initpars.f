ccc   sets parameters for soft survival model in=1,.,4
      subroutine initpars(in)
      implicit double precision(a-y)
      integer in,i1,i2

      include 'nchan.f'
      include 'vars.f'
      include 'survpars.f'

      if(in.eq.1)then

         ep=0.13d0
         asp=0.08d0
         ep1=0d0
         sigo=23d0/0.39d0
         gaa(3)=0.4d0
         gd2=0.3d0
         nch=2
         ntf=0
         
         cc0(1)=0.45d0
         bm(1)=3d0
         bb0(1)=0.1d0
         pp0(1)=0.92d0
         bex(1)=8.5d0
         
         cc0(2)=0.45d0
         bm(2)=1.5d0
         bb0(2)=0.5d0
         pp0(2)=0.1d0
         bex(2)=4.5d0
         
         cc0(3)=1d0
         bm(3)=0d0
         bb0(3)=0.8d0
         pp0(3)=0.5d0
         bex(3)=0.5d0

      elseif(in.eq.2)then
         
         ep=0.115d0
         asp=0.11d0
         ep1=0d0
         sigo=33d0/0.39d0
         gaa(3)=0.6d0
         gd2=0.16d0
         nch=2d0
         ntf=0d0

         cc0(1)=0.63d0
         bm(1)=3d0
         bb0(1)=0.1d0
         pp0(1)=0.5d0
         bex(1)=8d0
         
         cc0(2)=0.47d0
         bm(2)=1.5d0
         bb0(2)=0.5d0
         pp0(2)=0.1d0
         bex(2)=6d0
         
         cc0(3)=1d0
         bm(3)=0d0
         bb0(3)=0.8d0
         pp0(3)=0.5d0
         bex(3)=0.5d0
        
      elseif(in.eq.3)then

         ep=0.093d0
         asp=0.075d0
         ep1=0d0
         sigo=60d0/0.39d0
         gaa3=4.8d0
         gd2=1.03d0
         nch=2
         nga=1
         cc0(1)=0.55d0
         bm(1)=3d0
         bb0(1)=0.27d0
         pp0(1)=0.48d0
         bex(1)=5.3d0
         
         cc0(2)=0.48d0
         bm(2)=1.5d0
         bb0(2)=0.1d0
         pp0(2)=1d0
         bex(2)=3.8d0
         
         cc0(3)=0.24d0
         
      elseif(in.eq.4)then

         ep=0.11d0              !!! Capital delta
         asp=0.06d0             !!! alpha'
         ep1=0d0                !!! Zero in all models, matters (?)
         sigo=50d0/0.39d0       !!! sigma_0 (GeV^-2)
         gaa3=6d0               !!! k2/k(1.8 TeV)
         gd2=1.3d0              !!! k1/k(1.8 TeV)
         nch=2
         nga=1                  !!! A flag (?)
         cc0(1)=0.6d0           !!! d1
         bm(1)=3d0              !!! Doesn't matter
         bb0(1)=0.45d0          !!! c1-0.08 (added back later)
         pp0(1)=0.5d0           !!! 2*|a_1|^2
         bex(1)=7.2d0           !!! b1
         
         cc0(2)=0.48d0          !!! d2
         bm(2)=1.5d0            !!! Doesn't matter
         bb0(2)=0.16d0          !!! c2-0.08 (added back later)
         pp0(2)=1d0             !!! |a_2|^2 is set later
         bex(2)=4.2d0           !!! b2
         
         cc0(3)=0.12d0          !!! Beta : k^2_min ~ s^Beta
      
      endif

      if(nch.eq.3) pp0(3)=3d0-pp0(2)-pp0(1)
      if(nch.eq.2) pp0(2)=2d0-pp0(1) !!! Set |a_2|^2
      if(nch.eq.1) pp0(1)=1d0

      if(in.eq.3.or.in.eq.4)then

      gamm=(1800d0/rts)**cc0(3)
      ga1=1d0/(1d0+gamm*gd2)
      ga2=1d0/(1d0+gamm*gaa3)
      gaa(1)=2d0*ga1/(ga1+ga2)   !!! gamma_1 (+) Eq (15)
      gaa(2)=2d0*ga2/(ga1+ga2)   !!! gamma_2 (-) Eq (15)

      elseif(in.eq.1.or.in.eq.2)then

         if(nch.eq.2)then
            gaa(1)=1d0+dsqrt(gd2)
            gaa(2)=1d0-dsqrt(gd2)
            gaa(3)=0.
         elseif(nch.eq.1)then
            gaa(1)=1d0
            gaa(2)=0d0
            gaa(3)=0d0
         endif

      endif

      sum=0d0

      do i1=1,nch
         do i2=1,nch 
            sum=sum+gaa(i1)*gaa(i2)*pp0(i1)*pp0(i2)/dble(nch)**2
         enddo
      enddo
      
      norm=sum

      return
      end
