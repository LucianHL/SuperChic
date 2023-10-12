      subroutine tpcalc
      implicit none
      double precision qtmin,qtmax,sum,qt,lqtmin,lqtmax,lqt
      double precision tpz,tpn
      integer i,j,jmax
      
      include 'tppars.f'
      include 'ion.f'
      
      itp=1900

      qtmin=1d-3
      qtmax=2d0

      lqtmin=dlog(qtmin)
      lqtmax=dlog(qtmax)

      sum=0d0
      
      do i=1,itp+1

         lqt=lqtmin+(lqtmax-lqtmin)*dble(i-1)/dble(itp)
         qt=dexp(lqt)

         tparr(1,i)=lqt
         tparr(2,i)=tpz(qt)
         tparr(3,i)=tpn(qt)

c         print*,qt,tparr(2,i),tparr(3,i)
         
      enddo

c$$$      meta=0.6d0
c$$$      rts=8.12d3
c$$$
c$$$      ymin=-5d0
c$$$      ymax=5d0
c$$$
c$$$      jmax=40
c$$$      
c$$$c      do j=1,jmax+1
c$$$
c$$$         y=ymin+(ymax-ymin)*dble(j-1)/dble(jmax)
c$$$         
c$$$c      x=meta/rts*dexp(y)
c$$$
c$$$         x=6.7d-5
c$$$         
c$$$      itp=1900
c$$$
c$$$      ypmin=dlog(x**2*0.938d0**2)
c$$$      ypmax=dlog(x**2*0.938d0**2+qtmax**2)
c$$$
c$$$      
c$$$      sum=0d0
c$$$      
c$$$      do i=1,itp
c$$$
c$$$         yp=ypmin+(ypmax-ypmin)*(dble(i)-0.5d0)/dble(itp)
c$$$
c$$$         qt=dexp(yp)-x**2*0.938d0**2
c$$$         qt=dsqrt(qt)
c$$$         
c$$$         wt=1d0/137d0/3.141d0*tpint(1,qt)**2
c$$$c         wt=1d0/137d0/3.141d0*82d0**2
c$$$         wt=wt*qt**2/(qt**2+(x*0.938d0)**2)
c$$$         wt=wt*(ypmax-ypmin)/dble(itp)
c$$$
c$$$         print*,qt,tpint(1,qt),tpint(2,qt)
c$$$         
c$$$c         if(qt.gt.1d0)wt=0d0
c$$$         
c$$$         sum=sum+wt
c$$$         
c$$$c         print*,qt,tparr(2,i)
c$$$         
c$$$       enddo
c$$$
c$$$         print*,sum
c$$$         
c$$$c      print*,y,sum
c$$$
c$$$c      enddo
c$$$      
c$$$      stop

      return
      end
