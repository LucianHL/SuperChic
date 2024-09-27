ccc   reads Sudakov factor from file
C      subroutine calcsud
      subroutine readsud
      implicit double precision(a-y)
      integer i,itot,j,jtot,outl

      include 'sudpars.f'
      include 'intag.f'
      include 'vars.f'

ccccc MUST BE SAME AS IN INITSUD.F

      qmax=100d0
      qmin=0.4d0

      mmax=dlog(rts)
      mmin=dlog(2d0)

      itot=300
      jtot=200

      qinc=(qmax-qmin)/dble(itot)
      minc=(mmax-mmin)/dble(jtot)

      call length(intag,outl)

      open(42,file='inputs/sudakov'//intag(1:outl)//'.dat')

      do 555 i=1,itot+1
      do 555 j=1,jtot+1

         read(42,*)tgint(1,i,j),tgint(2,i,j),tgint(3,i,j),tgint(4,i,j)

 555  enddo

      close(42)

      return
      end
