ccc   writes out sudakov factor to file
      subroutine calcsud
      implicit double precision(a-y)
      integer i,itot,j,jtot,outl

      include 'intag.f'
      include 'vars.f'

      qmax=100d0
      qmin=0.4d0

      mmax=dlog(rts)
      mmin=dlog(2d0)

      itot=300
      jtot=200

      call length(intag,outl)

      open(42,file='inputs/sudakov'//intag(1:outl)//'.dat')

      qinc=(qmax-qmin)/dble(itot)
      minc=(mmax-mmin)/dble(jtot)

      print*,'Calculating Sudakov factor...'

      do 555 i=1,itot+1
      do 555 j=1,jtot+1
         
         qsq=qmin+qinc*(dble(i)-1d0)
         mass=mmin+minc*(dble(j)-1d0)

         call Sud(qsq,dexp(mass),tg,dtg)

         write(42,*)qsq,mass,tg,dtg

 555  enddo

      close(42)

      print*,'Done!'

      return
      end
