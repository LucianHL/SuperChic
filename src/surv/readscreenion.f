ccc   read in ion form factor from file
      subroutine readscreenion
      implicit none
      integer i,outl

      include 'intag.f'
      include 'scionpars.f'

      qmin=1d-3
      qmax=2d0

      lqmin=dlog(qmin)
      lqmax=dlog(qmax)

      itot=1900
      
      call length(intag,outl)

      open(40,file='inputs/screeningion'//intag(1:outl)//'.dat')


      do i=1,itot+1

         lq=lqmin+(lqmax-lqmin)*dble(i-1)/dble(itot)
         qt=dexp(lq)

         read(40,*)scionarr(1,i),scionarr(2,i)
         
      enddo

      close(40)

      return
      end


  
