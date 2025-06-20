      subroutine gdrset
      implicit none
      double precision sum1,sumx,hb,btmin,btmax,bt,acc
      double precision bti,lbtmin,lbtmax,hlb,lbt,sum
      logical recalc
      integer i

      include 'ion.f'
      include 'pgdr.f'
      include 'pi.f'
      include 'p0Xn.f'

      open(20,file='probtest1a.dat')

      btmax=rzg*500d0
      btmin=rzg*1.5d0

      ibmax=500

      lbtmax=dlog(btmax)
      lbtmin=dlog(btmin)

      acc=1d0

      hb=(btmax-btmin)/dble(ibmax)
      hlb=(lbtmax-lbtmin)/dble(ibmax)

      recalc=.true.


      do i=1,ibmax+1

         bt=btmin+hb*dble(i-1)
         lbt=lbtmin+hlb*dble(i-1)
         bt=dexp(lbt)

         call gdrcalc(bt,sum1,sumX)
         bti=bt

         gdrarr(1,i)=dlog(bti)
         gdrarr(2,i)=sum1*fracsigX
         gdrarr(3,i)=sumX*fracsigX

         if(wrho)then
         call rho_exc_calc(bt,sum)
         gdrarr(4,i)=sum
c         call twogamprob_calc(bt,sum)
c         gdrarr(5,i)=sum
         endif

      enddo


      close(20)

      stop

      return
      end

      function pgdrint(i,bt)
      implicit none
      double precision bt,pgdrint,rbin,out,m,del,btmin,btmax
      integer i,it

      include 'pgdr.f'

      btmax=gdrarr(1,ibmax+1)
      btmin=gdrarr(1,1)


      if(bt.gt.btmax)then
         it=ibmax-1
         m=(dlog(gdrarr(i,it+2))-dlog(gdrarr(i,it+1)))
         m=m/(gdrarr(1,it+2)-gdrarr(1,it+1))
         del=bt-gdrarr(1,it+1)
         out=m*del+dlog(gdrarr(i,it+1))
         out=dexp(out)
      else
         rbin=gdrarr(1,2)-gdrarr(1,1)
         it=nint((bt-btmin)/rbin)
         if(dble(it).gt.((bt-btmin)/rbin))then
            it=it-1
         endif


         if(bt.lt.btmin)it=0

         m=(gdrarr(i,it+2)-gdrarr(i,it+1))
         m=m/(gdrarr(1,it+2)-gdrarr(1,it+1))
         del=bt-gdrarr(1,it+1)
         out=m*del+gdrarr(i,it+1)
      endif

      pgdrint=out

      return
      end

      subroutine twogamprob_calc(bt,sum)
      implicit none
      double precision bt,sum,phi
      double precision btmax,wt,hbt,hphi
      double precision btx,bty,btxp,btyp,btp
      double precision bt1x,bt1y,bt2x,bt2y,bt1,bt2
      double precision x,mll,gdrint
      integer ib,nbt,iphi,nphi

      include 'ion.f'
      include 'pgdr.f'
      include 'pi.f'
      include 'p0Xn.f'
      include 'gdr.f'
      include 'sAA.f'
      include 'mion.f'
      include 'vars.f'
      include 'mp.f'


cc   0.25943604146912896    
      btmax=rzg*100d0
      btmax=dlog(btmax)

      nbt=100
      nphi=10

      hbt=btmax/dble(nbt)
      hphi=2d0*pi/dble(nphi)

      bty=bt
      btx=0d0

      mll=1d0
      x=mll/rtsaa

      sum=0d0

      do ib=1,nbt
      do iphi=1,nphi

         btp=(dble(ib)-0.5d0)*hbt
         btp=dexp(btp)
         phi=(dble(iphi)-0.5d0)*hphi

         btxp=btp*dcos(phi)
         btyp=btp*dsin(phi)

         bt1x=(btx+btxp)/2d0
         bt1y=(bty+btyp)/2d0
         bt2x=(btx-btxp)/2d0
         bt2y=(bty-btyp)/2d0

         bt1=dsqrt(bt1x**2+bt1y**2)
         bt2=dsqrt(bt2x**2+bt2y**2)

         wt=gdrint(x,bt1)*gdrint(x,bt2)
         wt=wt*hbt*hphi*btp**2/2d0

         sum=sum+wt

c         print*,btp/rzg,wt

      enddo
      enddo


      return
      end

      subroutine rho_exc_calc(bt,sum)
      implicit none
      double precision kmin,kmax,sig,htlk,lkmin,lkmax
      double precision ycut,xmin,xmax,mrho
      double precision lk,k,wt,sum,x,bt
      double precision gdrint
      double precision erho,prho,sig0,msig,wgam
      integer ik,nk

      include 'ion.f'
      include 'gdr.f'
      include 'sAA.f'
      include 'mion.f'
      include 'vars.f'
      include 'mp.f'
      include 'p0Xn.f'

      kmin=1d0
      kmax=1d4

      mrho=0.77545d0
      ycut=yrho

      xmin=mrho*dexp(-ycut)/rtsaa
      xmax=mrho*dexp(ycut)/rtsaa

      kmin=xmin/2d0/mion*(rtsaa**2-2d0*mion**2)
      kmax=xmax/2d0/mion*(rtsaa**2-2d0*mion**2)

      nk=100
      
      lkmin=dlog(kmin)
      lkmax=dlog(kmax)
      htlk=(lkmax-lkmin)/dble(nk)

      sig0=2d0/0.389379d0
      sig=sig0
      msig=0.2d0/0.389379d0/30d0

      sum=0d0

      do ik=1,nk
         
         lk=lkmin+htlk*(dble(ik)-0.5d0)
         k=dexp(lk)
         
         wgam=dsqrt(2d0*0.938d0*k)
         wt=htlk*sig*accrho

         x=2d0*mion*k/(rtsaa**2-2d0*mion**2)
         wt=wt*gdrint(x,bt)

         sum=sum+wt


      enddo

      x=mrho/rtsaa


      return
      end


      subroutine gdrcalc(bt,sum1,sumX)
      implicit none
      double precision x,wtx,wt1,sum1,sumx,bt
      double precision gdrint
      integer i

      include 'ion.f'
      include 'gdr.f'
      include 'sAA.f'
      include 'mion.f'
      include 'vars.f'
      include 'mp.f'

      sum1=0d0
      sumX=0d0

      do i=1,i0

         wt1=sneut(2,i)/sneut(1,i)

c         print*,i,sneut(1,i)

         x=2d0*mion*sneut(1,i)*1d-3/(rtsaa**2-2d0*mion**2)

         wt1=wt1*gdrint(x,bt)

         if(i.lt.i0)then
            wt1=wt1*(sneut(1,i+1)-sneut(1,i))
         else
            wt1=wt1*(sneut(1,i0)-sneut(1,i0-1))
         endif

         sum1=sum1+wt1

      enddo

      do i=1,i4

         wtX=mneut(2,i)/mneut(1,i)

         x=2d0*mion*mneut(1,i)*1d-3/(rtsaa**2-2d0*mion**2)

         wtX=wtX*gdrint(x,bt)

         if(i.lt.i4)then
            wtX=wtX*(mneut(1,i+1)-mneut(1,i))
         else
            wtX=wtX*(mneut(1,i4)-mneut(1,i4-1))
         endif

         sumX=sumX+wtX

      enddo

      do i=i4+1,i5

         x=2d0*mion*mneut(1,i)*1d-3/(rtsaa**2-2d0*mion**2)

         wtx=dlog(mneut(1,i))-dlog(mneut(1,i-1))
         wtx=wtx*gdrint(x,bt)*mneut(2,i)

         sumx=sumx+wtx

      enddo

      sum1=sum1/0.389389d0
      sumX=sumX/0.389389d0

      return
      end

      

      function gdrint(x,bt)  ! is |tilde{N}|^2 in (34) of 2303.04826
      implicit none
      double precision x,bt,gdrint,wt,sum,sum1
      double precision qt,q2min,q2max,lq2min,lq2max,qsq,q0,lqsq,hq2,f1
      double precision qtmin,hlq2
      double precision tpint
      integer i,itot
      include 'ion.f'
      include 'gdr.f'
      include 'mion.f'
      include 'pi.f'
      include 'gaussvars.f'

      q0=0.71d0

      qtmin=0d0

      q2min=(qtmin**2+x**2*mion**2)/(1d0-x)
      q2max=2d0

      if(q2min.gt.q2max)then
         gdrint=0d0
         return
      endif

      lq2min=dlog(q2min)
      lq2max=dlog(q2max)

      itot=1000

      itot=nint(bt/2d0)
      if(itot.lt.1000)itot=1000

      hlq2=(lq2max-lq2min)/dble(itot)
      hq2=(q2max-q2min)/dble(itot)

      sum=0d0
      sum1=0d0

      do i=1,itot

         lqsq=lq2min+hlq2*(dble(i)-0.5d0)
         qsq=dexp(lqsq)

         qt=(1d0-x)*qsq-x**2*mion**2
         qt=dsqrt(qt)

         f1=1d0/(1d0+qsq/q0)**2

         wt=tpint(1,dsqrt(qsq))*qt
         wt=wt*BESSEL_J1(bt*qt)*f1

         wt=wt*hlq2

         sum=sum+wt

      enddo

      sum=sum**2*(1d0-x)/4d0/pi**2/137d0

      gdrint=sum

      return
      end
