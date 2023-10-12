      subroutine sdcoh
      implicit double precision(a-y)
      integer ikt,nkt,nbt,ibt,ikt1,ikt2
      integer iin,i1,i2,i3,i4
      integer iphi1,iphi2,nphi,iphi
      integer ibt1,ibt2
      integer ipt,npt,ipphi,npphi
      integer n,ntotal,ntotal1
      integer ix1,nx1,outl
      integer iq,iqtot
      integer in1,in2

      include 'nchan.f'
      include 'survpars.f'
      include 'pi.f'
      include 'vars.f'
      include 'intag.f'
      include 'mp.f'

      print*,'Calculating S^2 (Elastic + Evolution)...'

      call length(intag,outl)
      open(10,file='inputs/sdcoh'//intag(1:outl)//'.dat')

cccccccc

      ntotal=1000000
      ntotal1=200
    
      ktmax=3d0
 

      crossb=0d0
      cross=0d0

      sige=sigo*dexp(dlog(rts)*2d0*ep)

   
ccccccccccc

      xmin=1d-2
      xmax=0.6d0
      
      lxmin=dlog(xmin)
      lxmax=dlog(xmax)
      
      nx1=4
      iqtot=10

      do 902 ix1=1,nx1+1

         lnx1=lxmin+(lxmax-lxmin)*dble(ix1-1)/dble(nx1)
         x=dexp(lnx1)

      do 902 iq=1,iqtot
                  
         qtmax=3d0
         qtmin=0.1d-1

         lqtmax=dlog(qtmax)
         lqtmin=dlog(qtmin)
         
         lqt=(lqtmax-lqtmin)*(dble(iq))/dble(iqtot)+lqtmin
         qt=dexp(lqt)

         pt=qt

         wt=1d0/(1d0-x)
         sum=wt*fem(x,pt**2,1)**2*pt**2

         norm0=1d0
         sum0=sum*norm0

      sum=0d0

      do 900 in1=1,ntotal1
      do 900 in2=1,ntotal1

         rphi=(dble(in1)-0.5d0)/dble(ntotal1)
         rkt=(dble(in2)-0.5d0)/dble(ntotal1)
         
         phi=2d0*pi*rphi
         kt=ktmax*rkt

         pt=qt
         ktc=(kt**2+x**2*mp**2)/(1d0-x)
  
         opac=0d0
         do i1=1,nch
            do i2=1,nch
               call screeningint(i1,i2,kt**2,sc,sc1) 
               opac=opac+sc*pp0(i1)*pp0(i2)/dble(nch)**2
     &              /(gaa(i1)*gaa(i2))
     &              *fe(ktc,i1)
            enddo
         enddo
         
         wt=2d0*pi*ktmax*kt
         wt=wt/(1d0-x)

         pt1=dsqrt(pt**2+kt**2+2d0*dcos(phi)*kt*pt)

         ptp=pt**2+pt*kt*dcos(phi)

         sum=sum+wt*opac*fem(x,pt1**2,1)*fem(x,pt**2,1)*ptp

 900  enddo

      sum=sum/dfloat(ntotal1)**2


      sum1=0d0
      sum1v=0d0


         do 810 n=1,ntotal

            rkt1=ran2()
            rkt2=ran2()
            rphi1=ran2()
            rphi2=ran2()
            
            pt=qt
            
            
            kt1=ktmax*rkt1
            kt2=ktmax*rkt2
            phi1=2d0*pi*rphi1
            phi2=2d0*pi*rphi2

 
            wt=ktmax**2*4d0*pi**2
            wt=wt*kt1*kt2

            wt=wt/(1d0-x)
       
         ktc=kt1**2+kt2**2+2d0*kt1*kt2*dcos(phi1-phi2)
         ktc=(ktc+x**2*mp**2)/(1d0-x)

         opac=0d0
         
         do i1=1,nch
            do i2=1,nch
               do i3=1,nch
                  do i4=1,nch
                     call screeningint(i1,i2,kt1**2,sca,sc1) 
                     call screeningint(i3,i4,kt2**2,scb,sc1) 
                     opac=opac+sca*scb*pp0(i1)*pp0(i2)/dble(nch)**2
     &                    *pp0(i3)*pp0(i4)/dble(nch)**2
     &                    /(gaa(i1)*gaa(i2))/(gaa(i3)*gaa(i4))
     &                    *fe(ktc,i1)
                  enddo
               enddo
            enddo
         enddo
         
         pt1=dsqrt(pt**2+kt1**2+2d0*kt1*pt*dcos(phi1))
         pt2=dsqrt(pt**2+kt2**2+2d0*kt2*pt*dcos(phi2))
         
         ptp=pt**2+pt*kt1*dcos(phi1)+pt*kt2*dcos(phi2)+kt1*
     &        kt2*dcos(phi2-phi1)

         
         sum1=sum1+wt*opac*fem(x,pt1**2,1)*fem(x,pt2**2,1)*ptp
         
 810  enddo

      sum1=sum1/dfloat(ntotal)

      surv=1d0+2d0*sum/sum0+sum1/sum0

c      print*,x,qt,surv
c      print*,sum0,sum,sum1

      if(sum0.eq.0d0.or.surv.lt.0d0)then
         write(10,*)dlog(x),dlog(qt),0d0
c         print*,x,qt,0d0
       else
          write(10,*)dlog(x),dlog(qt),1d0+2d0*sum/sum0+sum1/sum0
c          print*,x,qt,surv
      endif

 902  enddo
     
 887  close(10)

      
cccccccccccccccccc


      return
      end


      function fem(x,qtsq,i)
      implicit double precision(a-y)
      integer i

      include 'mp.f'

      qsq=(qtsq+x**2*mp**2)/(1d0-x)

      
      call F1F2el(qsq,f1,fem)
c      call F1F2el(qtsq,f1,fem)
      
      fem=dsqrt(fem)/qsq

      return
      end

      function fe(qsq,i)
      implicit double precision(a-y)
      integer i

      call F1F2el(qsq,f1,fout)
      fe=dsqrt(fout)
c      fe=fout

      return 
      end
