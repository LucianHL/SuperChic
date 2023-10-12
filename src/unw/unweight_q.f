ccc   writes event information to array for unweighted generation
      subroutine unweightq(wt,r)
      implicit double precision(a-y)
      integer i,j,k

      include 'record.f'
      include 'unweighted.f'
      include 'mom.f'
      include 'proc.f'
      include 'decay.f'
      include 'hepevt.f'
      include 'leshouches.f'
      include 'wmax.f'
      include 'diss.f'
      include 'diff.f'
      include 'pdg.f'

      if(wt/wmax.gt.r)then
         
         evnum=evnum+1

         if(diff.eq.'dd')then
            do i=1,4
               evrec(evnum,3,i)=q(i,nup+1)
               evrec(evnum,4,i)=q(i,nup+2)
               evrec(evnum,5,i)=q(i,nup+3)
               evrec(evnum,6,i)=q(i,nup+4)
            enddo
         elseif(diff.eq.'sd')then
            do i=1,4
               evrec(evnum,3,i)=q(i,nup+2)
               evrec(evnum,4,i)=q(i,nup+3)
               evrec(evnum,5,i)=q(i,nup+4)
            enddo
         elseif(diff.eq.'el')then
            pdgid(3)=22
            istup(3)=-1
            pdgid(4)=22
            istup(4)=-1 
            do i=1,4
               evrec(evnum,3,i)=q(i,1)-q(i,3)
               evrec(evnum,4,i)=q(i,2)-q(i,4)
            enddo
         endif
         
         if(diss1.and.diss2)then
            do j=7,nup+1
               do i=1,4
                  k=j-1
                  evrec(evnum,j,i)=q(i,k)
               enddo
            enddo
         elseif(diss1.or.diss2)then
            do j=6,nup+1
               do i=1,4
                  k=j
                  evrec(evnum,j,i)=q(i,k)
               enddo
            enddo
         else
            do j=5,nup+2
               do i=1,4
                  evrec(evnum,j,i)=q(i,j+1)
               enddo
            enddo
         endif
         
      endif
            
      return
      end
