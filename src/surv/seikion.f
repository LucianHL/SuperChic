ccc   integrates bare + screened amplitude over k_t 
ccc   (QCD induced processes)
      subroutine schimcion(p1x,p1y,p2x,p2y,out)
      implicit none
      double precision wt,x00p,qtmax1,qtmax,qt,phiq
      double precision tpx,tpy,tp2,tp1,t11,t12,t22,tpint
      double precision betaionex,sc,screeningionint
      double precision p1xp,p1yp,p2xp,p2yp,mu,hphi,del
      double precision p1x,p1y,p2x,p2y
      integer jx,jy,i1,i2,p,nphi,nqt,jqt,jphi,nnphi,nk,nk1
      complex*16 out(10),x0(10),x00(10),x0p(10),outt
      complex*16 screen(2,2)
c      complex*16 outt(0:4000,10)
      integer icount
      common/icount/icount

      include 'nchan.f'
      include 'mt.f'
      include 'surv.f'
      include 'vars.f'
      include 'survpars.f'
      include 'polarization.f'
      include 'nsurv.f'
      include 'pi.f'
      include 'zarr.f'
      include 'beam.f'
      include 'rho0.f'
      include 'ionqcd.f'
      include 'gaussvars.f'
      
      call setmu(mu)

      do p=1,pol
         out(p)=0d0
      enddo
      outt=0d0
      
      nphi=s2int*4

      hphi=2d0*pi/dble(nphi)

      nk=s2int*4
      nk1=s2int*4

      del=0.5d0

      qtmax=0.5d0
      qtmax1=1.5d0

      call bare(mu,p1x,p1y,p2x,p2y,x0p)
      
      if(sfac)then
         
        do jqt=1,nk+nk1

           if(jqt.le.nk)then
              qt=qtmax/2d0*(xikt4(jqt)+1d0)
           else
              qt=((qtmax1-qtmax)*xikt4(jqt-nk)+qtmax1+qtmax)/2d0
           endif
           
           sc=screeningionint(qt)
            
           do jphi=1,nphi
              
              phiq=pi*(xiphi4(jphi)+1d0)
               
               tpx=qt*dcos(phiq)
               tpy=qt*dsin(phiq)

               if(jqt.le.nk)then
                  wt=wikt4(jqt)*qtmax/2d0*qt
               else
                  wt=wikt4(jqt-nk)*(qtmax1-qtmax)/2d0*qt
               endif
               
               wt=wt*pi*wiphi4(jphi)
               
               p1xp=p1x-tpx
               p1yp=p1y-tpy
               t12=p1xp**2+p1yp**2
               p2xp=tpx+p2x
               p2yp=tpy+p2y
               t22=p2xp**2+p2yp**2
               
               x00p=betaionex(-t12)*betaionex(-t22)              
                
            if(beam.eq.'ion')then
               tp1=tpint(1,dsqrt(t12))+tpint(2,dsqrt(t12))
               tp2=tpint(1,dsqrt(t22))+tpint(2,dsqrt(t22))
               x00p=x00p*tp1*tp2
            elseif(beam.eq.'ionp')then
               tp1=tpint(1,dsqrt(t12))+tpint(2,dsqrt(t12))
               x00p=x00p*tp1
            endif
               
           do p=1,pol
              
              x0(p)=x0p(p)
              x0(p)=x0(p)*x00p
   
              out(p)=out(p)+x0(p)*wt*sc

           enddo
           
        enddo
      enddo

      icount=icount+1
      
      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2
      tp1=tpint(1,dsqrt(t11))+tpint(2,dsqrt(t11))
      tp2=tpint(1,dsqrt(t22))+tpint(2,dsqrt(t22))
      
      x00p=betaionex(-t11)*betaionex(-t22)
      
      endif

      t11=p1x**2+p1y**2
      t22=p2x**2+p2y**2

      x00p=betaionex(-t11)*betaionex(-t22)

      if(beam.eq.'ion')then
         tp1=tpint(1,dsqrt(t11))+tpint(2,dsqrt(t11))
         tp2=tpint(1,dsqrt(t22))+tpint(2,dsqrt(t22))
         if(ionqcd.eq.'coh')x00p=x00p*tp1*tp2
      elseif(beam.eq.'ionp')then
         tp1=tpint(1,dsqrt(t11))+tpint(2,dsqrt(t11))
         if(ionqcd.eq.'coh')x00p=x00p*tp1
      endif
      
      do p=1,pol
         out(p)=out(p)+x0p(p)*x00p
         if(cdabs(x0p(p)*x00p).lt.cdabs(out(p))*del)
     &        out(p)=x0p(p)*x00p
      enddo

      return
      end
