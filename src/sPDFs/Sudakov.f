ccc   calculates Sudakov factor
      subroutine Sud(qtsq,mx,tg,dtg)
      implicit double precision (a-z)
      integer i,itot

      include 'pi.f'    

      lam=0.5d0
      llam=dlog(lam**2)

      sum=0d0

      itot=100

      lnktmax=dlog(mx**2/4d0)
      lnktmin=dlog(qtsq)
      llnktmax=dlog(lnktmax-llam)
      llnktmin=dlog(lnktmin-llam)
      hllnkt=(llnktmax-llnktmin)/dble(itot)    

      do i=1,itot

         llnkt=llnktmin+hllnkt*(dble(i)-0.5d0)
         lnkt=dexp(llnkt)+llam
         ktsq=dexp(lnkt)
         delta=dsqrt(ktsq)/mx

         wt=alphas(ktsq)/2d0/pi
         pggi=(dlog(1d0/delta)-(1d0-delta)**2*(3d0*delta**2
     &        -2d0*delta+11d0)/12d0)*2d0*3d0
         pqgi=(2d0*(1d0-delta**3)/3d0-delta*(1d0-delta))*0.5d0*nf(ktsq)
         wt=wt*(pggi+pqgi)
         wt=wt*dexp(llnkt)
         wt=wt*hllnkt

         sum=sum+wt
 
      enddo

      Tg=dexp(-sum)
      delta=dsqrt(qtsq)/mx

      pggi=(dlog(1d0/delta)-(1d0-delta)**2*(3d0*delta**2
     &     -2d0*delta+11d0)/12d0)*2d0*3d0
      pqgi=(2d0*(1d0-delta**3)/3d0-delta*(1d0-delta))*0.5d0*nf(qtsq)
      ftg=alphas(qtsq)/2d0/pi*(pggi+pqgi)/qtsq

      dtg=tg*ftg

      return
      end

   
  
