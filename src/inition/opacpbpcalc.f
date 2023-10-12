      subroutine opacpbpcalc
      implicit double precision(a-y)
      integer i
      
      include 'opacpbppars.f'
      include 'ion.f'
      include 'pi.f'
      include 'ionqcd.f'
      include 'qcd.f'
      
      btmax=20d0
      ioppbp=900

      btmax=1.5d0*rzg
      btmin=0d0
      hb=(btmax-btmin)/dble(ioppbp)

      if(ionqcd.eq.'incoh')sigin=6.5d0 ! fm^2
      if(ionqcd.eq.'coh')sigin=6.74d0 ! fm^2
         
      do i=1,ioppbp+1

         bt=btmin+(dble(i)-1d0)*hb
  
         opacpbparr(1,i)=bt
         opac=opacpbp(bt)
         if(qcd)opac=opac*(1d0-1d0/nshell)
         if(opac.gt.100d0)then
            opacpbparr(2,i)=0d0
         else
            opacpbparr(2,i)=dexp(-opac/2d0)
         endif
         
      enddo
      
      return
      end
