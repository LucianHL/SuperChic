      subroutine screencalc
      implicit none
      double precision qmin,qmax,qt,lqmin,lqmax,lq
      double precision screen,screen_01
      integer i

      include 'scionpars.f'
      include 'p0Xn.f'

      qmin=1d-5
      qmax=1.5d0

      lqmin=dlog(qmin)
      lqmax=dlog(qmax)

      itot=400

c      open(10,file='testold.dat')

      do i=1,itot+1

         lq=lqmin+(lqmax-lqmin)*dble(i-1)/dble(itot)
         qt=dexp(lq)

         scionarr(1,i)=lq

         if(pAAvar)then
            scionarr(ifaa+1,i)=screen(qt)
         elseif(int_01)then
            scionarr(2,i)=screen_01(qt)
            scionarr(2,i)=scionarr(2,i)+screen(qt)
         else
            scionarr(2,i)=screen(qt)
         endif

      enddo

      return
      end
