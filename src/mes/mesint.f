ccc   interpolates gg-->MM amplitude in cos(theta)
      subroutine mesint(pol,j,cost,out)
      implicit double precision(a-y)
      integer i,j,p,pol

      include 'mes.f'

      i=nint((cost-cmin)/cinc)
      if(dble(i).lt.(cost-cmin)/cinc)then
         i=i+1
      endif

      if(i.eq.0)i=i+1

      m=(mesamp(pol,2,j,i+1)-mesamp(pol,2,j,i))/cinc
      out=mesamp(pol,2,j,i)+m*(cost-mesamp(pol,1,j,i))

      return
      end
