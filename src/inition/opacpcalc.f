      subroutine opacpcalc
      implicit double precision(a-y)
      integer i
      
      include 'opacppars.f'
           
      btmax=20d0
      iopp=500

      hb=btmax/dble(iopp)

      do i=1,iopp+1
         
         bt=(dble(i)-1d0)*hb
         opacparr(1,i)=bt
         opacparr(2,i)=opacp(bt)
         
      enddo

      return
      end
