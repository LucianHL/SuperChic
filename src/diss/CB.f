C--   Christy-Bosted fit to resonance sig(L+T) - 0712.3731
      
      function sigtlnr(x,qsq)
      implicit double precision(a-z)
      integer i
      dimension m(7),gam(7),at(7),a(7),b(7),c(7)
     &     ,al(7),d(7),e(7),l(7),betapi(7),beta2pi(7),
     &     betaeta(7),xi(7)
      dimension signrt(2),signrl(1),anr(3),bnr(3),cnr(3),
     &     dnr(3),enr(3)
      data signrt /
     &     246.1,-89.4/
      data signrl /
     &     86.7/
      data anr /
     &     0.0675,0.2098,0./
      data bnr /
     &     1.3501,1.5715,4.0294/
      data cnr /
     &     0.1205,0.0907,3.1285/
      data dnr /
     &     -0.0038,0.0104,0.334/
      data enr /
     &     0.,0.,4.9623/
      data m /
     &     1.23,1.53,1.506,1.698,1.665,1.433,1.934/
      data gam /
     &     0.136,0.22,0.083,0.096,0.109,0.379,0.38/
      data at /
     &     7.78,6.335,0.603,2.33,1.979,0.0225,3.419/
      data a /
     &     4.229,6823.2,21.24,-0.288,-0.562,462.13,0./
      data b /
     &     1.26,33521.,0.056,0.186,0.39,0.192,0.0/
      data c /
     &     2.124,2.569,2.489,0.064,0.549,1.914,1./
      data al /
     &     29.414,0.,157.92,4.216,13.764,5.5124,11./
      data d /
     &     19.91,0.,97.046,0.038,0.314,0.054,1.895/
      data e /
     &     0.226,0.,0.31,1.218,3.,1.309,0.514/
      data l/
     &     1.5,0.5,1.5,2.5,2.5,0.5,1.5/
      data betapi/
     &     1.0,0.45,0.65,0.65,0.4,0.65,0.5/
      data beta2pi/
     &     0.,0.1,0.35,0.35,0.5,0.35,0.5/
      data betaeta/
     &     0.,0.45,0.,0.,0.1,0.,0./
      data xi/
     &     0.1446,0.215,0.215,0.215,0.215,0.215,0.215/

      sum=0d0
      
      mp=0.938d0
      mpi=0.135d0
      meta=0.5479d0
      q0t=0.05d0
      m0l=4.2802d0
      q0l=0.125d0

      wsq=mp**2+qsq*(1d0-x)/x   
      
      k=(wsq-mp**2)/2d0/mp
      kcm=(wsq-mp**2)/2d0/dsqrt(wsq)
      
      do i=1,7

         ki=(m(i)**2-mp**2)/2d0/mp
         kicm=(m(i)**2-mp**2)/2d0/m(i)
         
         gamig=gam(i)*(kcm/kicm)**2*(kicm**2+xi(i)**2)
     &        /(kcm**2+xi(i)**2)

         ej1=(wsq-mp**2+mpi**2)/dsqrt(wsq)
         pj1=dsqrt(ej1**2-mpi**2)

         ej2=(wsq-mp**2+4d0*mpi**2)/dsqrt(wsq)
         pj2=dsqrt(ej2**2-mpi**2)

  
         ej1m=(m(i)**2-mp**2+mpi**2)/m(i)
         pj1m=dsqrt(ej1m**2-mpi**2)

         ej2m=(m(i)**2-mp**2+4d0*mpi**2)/m(i)
         pj2m=dsqrt(ej2m**2-mpi**2)

         gam1=gam(i)*(pj1/pj1m)**(2d0*l(i)+1d0)
         gam1=gam1*((pj1m**2+xi(i)**2)/(pj1**2+xi(i)**2))**l(i)

         gam2=gam(i)*(pj2/pj2m)**(2d0*l(i)+4d0)
         gam2=gam2*((pj2m**2+xi(i)**2)/(pj2**2+xi(i)**2))**(l(i)+2d0)
         gam2=gam2*dsqrt(wsq)/m(i)

         if((wsq-mp**2+meta**2).gt.meta*dsqrt(wsq))then

            ej3=(wsq-mp**2+meta**2)/dsqrt(wsq)
            pj3=dsqrt(ej3**2-meta**2)
            ej3m=(m(i)**2-mp**2+meta**2)/m(i)
            pj3m=dsqrt(ej3m**2-meta**2)
            gam3=gam(i)*(pj3/pj3m)**(2d0*l(i)+1d0)
            gam3=gam3*((pj3m**2+xi(i)**2)/(pj3**2+xi(i)**2))**l(i)

         else

            gam3=0d0

         endif
            
         gamitot=betapi(i)*gam1+beta2pi(i)*gam2+betaeta(i)*gam3

         
         bw=gamitot*gamig/((wsq-m(i)**2)**2+(m(i)*gamitot)**2)
         bw=bw*ki*kicm/k/kcm/gam(i)
      
         atq=at(i)/(1d0+qsq/0.91d0)**c(i)
         atq=atq*(1d0+a(i)*qsq/(1d0+b(i)*qsq))

         alq=al(i)*qsq*dexp(-e(i)*qsq)/(1d0+d(i)*qsq)     
         
         sigi=dsqrt(wsq)*bw*(atq**2+alq**2)
         
         sum=sum+sigi

      enddo

      xpt=(1d0+(wsq-(mp+mpi)**2)/(qsq+q0t))
      xpt=1d0/xpt

      snr1=signrt(1)*(dsqrt(wsq)-mpi-mp)**1.5d0
      snr1=snr1/(qsq+anr(1))**(bnr(1)+cnr(1)*qsq+dnr(1)*qsq**2)
      
      snr2=signrt(2)*(dsqrt(wsq)-mpi-mp)**2.5d0
      snr2=snr2/(qsq+anr(2))**(bnr(2)+cnr(2)*qsq+dnr(2)*qsq**2)
 
      nrt=(snr1+snr2)*xpt

      xpl=(1d0+(wsq-(mp+mpi)**2)/(qsq+q0l))
      xpl=1d0/xpl
  
      tl=dlog(dlog((qsq+m0l)/0.33d0**2)/dlog(m0l/0.33d0**2))
      
      snrl=signrl(1)*(1d0-xpl)**(anr(3)*tl+bnr(3))/(1d0-x)
      snrl=snrl*qsq**cnr(3)/(qsq+q0l)**(1d0+cnr(3))
     &     *xpl**(dnr(3)+enr(3)*tl)
      
      sum=sum+nrt+snrl
      
      sigtlnr=sum
        
      return
      end
