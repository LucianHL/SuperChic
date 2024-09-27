ccc   read in skewed pdf from file
C      subroutine calchg
      subroutine readhg
      implicit double precision(a-y)
      integer i,itot,j,jtot,outl

      include 'hgpars.f'
      include 'vars.f'
      include 'intag.f'

      mmin=2d0

      call length(intag,outl)

      open(44,file='inputs/hg'//intag(1:outl)//'.dat')

      qmax=100d0
      qmin=0.4d0

      ycut=dlog(rts/mmin)

      xmin=mmin*dexp(-ycut)/rts
      xmax=1d0

      lxmax=dlog(xmax)
      lxmin=dlog(xmin)

      itot=650

      itot=900
      jtot=400


      xinc=(xmax-xmin)/dble(jtot)
      qinc=(qmax-qmin)/dble(itot)

      lxinc=(lxmax-lxmin)/dble(jtot)

      do 555 i=1,itot+1
      do 555 j=1,jtot+1

         read(44,*)hgint(1,i,j),hgint(2,i,j),hgint(3,i,j),hgint(4,i,j)

 555  enddo

      close(44)

      return
      end
