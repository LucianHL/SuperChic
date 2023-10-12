ccc   writes event information to array for unweighted generation
      subroutine unweight(wt,r)
      implicit double precision(a-y)
      integer i,j

      include 'record.f'
      include 'unweighted.f'
      include 'mom.f'
      include 'proc.f'
      include 'decay.f'
      include 'hepevt.f'
      include 'leshouches.f'
      include 'wmax.f'

      if(wt/wmax.gt.r)then
         
         evnum=evnum+1
         
         do j=3,nup+2
            do i=1,4
               evrec(evnum,j,i)=q(i,j)
            enddo
         enddo
         
      endif
            
      return
      end
