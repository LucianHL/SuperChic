ccc   writes out skewed PDF to file
      subroutine calchg
      implicit double precision(a-y)
      integer i,itot,j,jtot,outl

      include 'vars.f'
      include 'intag.f'

      call length(intag,outl)

      open(44,file='inputs/hg'//intag(1:outl)//'.dat')

      qmax=100d0
      qmin=0.4d0

      mmin=2d0

      ycut=dlog(rts/mmin)

      xmin=mmin*dexp(-ycut)/rts         
      xmax=1d0

      lxmax=dlog(xmax)
      lxmin=dlog(xmin)

      itot=900
      jtot=400

      xinc=(xmax-xmin)/dble(jtot)
      qinc=(qmax-qmin)/dble(itot)

      lxinc=(lxmax-lxmin)/dble(jtot)

      print*,'Calculating Skewed PDF...'

      do 555 i=1,itot+1
      do 555 j=1,jtot+1
         
         qsq=qmin+qinc*(dble(i)-1d0)
         x=xmin+xinc*(dble(j)-1d0)

         lnx=lxmin+lxinc*(dble(j)-1d0)
         x=dexp(lnx)

         if(x.gt.1d0)x=1d0

         call hpdf(x,qsq,tg,dtg)

         write(44,*)qsq,lnx,tg,dtg

 555  enddo

      close(44)

      print*,'Done!'

      return
      end
