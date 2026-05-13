ccc   boosts from nucleon-nucleon c.m.s to lab frame for pA collsions
      subroutine pAboost(io)
      implicit none
      double precision spa,rtspa,beta,ein,betaa
      double precision pcm(4),pboo(4),px(4)
      integer i,j,k,jmax,io

      include 'mom.f'
      include 'leshouches.f'
      include 'ion.f'
      include 'vars.f'
      include 'varsi.f'
      include 'hepevt.f'
      include 'mp.f'
      include 'ion_inel.f'

      do j=1,20
      do i=1,4
            q_rf(i,j)=q(i,j)
      enddo
      enddo

      if(AA_frame)return

      rtspA=rtsi*dsqrt(an/az)
      spA=rtspA**2
      beta=dsqrt(1d0-4d0*mp**2/spA)
      betaa=dsqrt(1d0-4d0*mp**2/spA*an**2/az**2)
      px(1)=0d0
      px(2)=0d0

      if(ion_inel)then
      if(io.eq.1)then
      px(3)=rtsnn/2d0*(beta-an*betaa)
      px(3)=rtsnn/2d0*(1d0-an)
      else
      px(3)=-rtsnn/2d0*(beta-an*betaa)
      endif
      px(4)=rtsnn/2d0*(1d0+an)
      else
cccc  proton in +/- z direction      
      if(io.eq.1)then
      px(3)=rtspA/2d0*(beta-az*betaa)
      else
      px(3)=-rtspA/2d0*(beta-az*betaa)
      endif
      px(4)=rtspA/2d0*(1d0+az)
      endif

      ein=dsqrt(px(4)**2-px(3)**2-px(2)**2-px(1)**2)

      jmax=nup+2
      if(ion_inel.and.ion_incoh_type.eq.'inel')jmax=jmax-1

      yx=0.5d0*dlog((q(4,5)+q(3,5))/(q(4,5)-q(3,5)))

      do j=1,jmax
         do i=1,4
            pcm(i)=q(i,j)
         enddo
         call boost(ein,px,pcm,pboo)
         do i=1,4
            q(i,j)=pboo(i)
         enddo
      enddo

      yx=0.5d0*dlog((q(4,5)+q(3,5))/(q(4,5)-q(3,5)))

      do k=1,2
         do j=1,4
            pup(j,k)=q(j,k)
         enddo
         do j=1,4
            phep(j,k)=q(j,k)
         enddo
         phep(5,k)=dsqrt(q(4,k)**2-q(3,k)**2-q(2,k)**2-q(1,k)**2)
         pup(5,k)=phep(5,k)
      enddo

      return
      end

ccc   boosts from nucleon-nucleon c.m.s to lab frame for pA collsions
      subroutine pAboost_8(io)
      implicit none
      double precision spa,rtspa,beta,ein,betaa
      double precision pcm(4),pboo(4),px(4)
      integer i,j,k,jmax,io

      include 'mom.f'
      include 'leshouches.f'
      include 'ion.f'
      include 'vars.f'
      include 'varsi.f'
      include 'hepevt.f'
      include 'mp.f'
      include 'ion_inel.f'

      do j=18,18
      do i=1,4
            q_rf(i,j)=q(i,j)
      enddo
      enddo

      if(AA_frame)return

      rtspA=rtsi*dsqrt(an/az)
      spA=rtspA**2
      beta=dsqrt(1d0-4d0*mp**2/spA)
      betaa=dsqrt(1d0-4d0*mp**2/spA*an**2/az**2)
      px(1)=0d0
      px(2)=0d0

      if(ion_inel)then
      if(io.eq.1)then
      px(3)=rtsnn/2d0*(beta-an*betaa)
      else
      px(3)=-rtsnn/2d0*(beta-an*betaa)
      endif
      px(4)=rtsnn/2d0*(1d0+an)
      else
cccc  proton in +/- z direction      
      if(io.eq.1)then
      px(3)=rtspA/2d0*(beta-az*betaa)
      else
      px(3)=-rtspA/2d0*(beta-az*betaa)
      endif
      px(4)=rtspA/2d0*(1d0+az)
      endif

      ein=dsqrt(px(4)**2-px(3)**2-px(2)**2-px(1)**2)

      jmax=nup+2
      if(ion_inel.and.ion_incoh_type.eq.'inel')jmax=jmax-1

      do j=18,18
         do i=1,4
            pcm(i)=q(i,j)
         enddo
         call boost(ein,px,pcm,pboo)
         do i=1,4
            q(i,j)=pboo(i)
         enddo
      enddo

      yx=0.5d0*dlog((q(4,5)+q(3,5))/(q(4,5)-q(3,5)))

      do k=1,2
         do j=1,4
            pup(j,k)=q(j,k)
         enddo
         do j=1,4
            phep(j,k)=q(j,k)
         enddo
         phep(5,k)=dsqrt(q(4,k)**2-q(3,k)**2-q(2,k)**2-q(1,k)**2)
         pup(5,k)=phep(5,k)
      enddo

      return
      end
