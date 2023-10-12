      subroutine qinit(x,qsq,i)
      implicit double precision(a-y)
      double precision pq(4)
      integer i,j
      
      include 'mom.f'
      include 'leshouches.f'
      include 'mp.f'
      include 'diff.f'
      include 'pdg.f'
      include 'vars.f'

      mq=0.33d0

      q0=q(4,i+2)-q(4,i)
      q3=q(3,i+2)-q(3,i)

      aq=q0**2-q3**2
      bq=-2d0*q0*qsq/rts
      cq=4d0*mq**2*q3**2/s+qsq**2/s
      rootq=bq**2-4d0*aq*cq

      xit=-bq-dsqrt(rootq)
      xit=xit/2d0/aq

      betaqt=dsqrt(1d0-4d0*mq**2/s/xit**2)
      pq(4)=xit*rts/2d0
      if(i.eq.1)pq(3)=pq(4)*betaqt
      if(i.eq.2)pq(3)=-pq(4)*betaqt

      do j=1,4
         if(diff.eq.'dd')q(j,nup+i)=pq(j)
         if(diff.eq.'dd')q(j,nup+i+2)=pq(j)-q(j,i)+q(j,i+2)
         if(diff.eq.'sd')q(j,nup+3)=pq(j)
         if(diff.eq.'sd')q(j,nup+4)=pq(j)-q(j,i)+q(j,i+2)
      enddo

      if(diff.eq.'sd')then ! elastic photon assignment
         do j=1,4
            if(i.eq.1)q(j,nup+2)=q(j,2)-q(j,4)
            if(i.eq.2)q(j,nup+2)=q(j,1)-q(j,3)
         enddo
      endif

      if(diff.eq.'dd')then
         pdgid(i+2)=2
         pdgid(i+4)=2
         istup(i+2)=-1
         istup(i+4)=1
         icolup(1,i+2)=500+i
         icolup(1,i+4)=500+i
         mothup(1,i+4)=1
         mothup(2,i+4)=2
      endif
      
      if(diff.eq.'sd')then
         pdgid(4)=2
         istup(4)=-1
         icolup(1,4)=501
         pdgid(5)=2
         istup(5)=1
         icolup(1,5)=501
         pdgid(3)=22
         istup(3)=-1
         mothup(1,5)=1
         mothup(2,5)=2
      endif
      
      return
      end

      
