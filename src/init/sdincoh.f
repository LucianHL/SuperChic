      subroutine sdincoh
      implicit double precision(a-z)
      integer ikt,nkt,nbt,ibt,ikt1,ikt2
      integer iin,i1,i2,i3,i4
      integer iphi1,iphi2,nphi,iphi
      integer ibt1,ibt2
      integer ipt,npt,ipphi,npphi
      integer n,ntotal,ntotal0,ntotal1
      integer ix1,nx1,outl
      integer iq,iqtot
      integer in1,in2,in3


      include 'nchan.f'
      include 'survpars.f'
      include 'pi.f'
      include 'vars.f'
      include 'intag.f'
      include 'mp.f'

      print*,'Calculating S^2 (Inelastic + Evolution)...'
      
      call length(intag,outl)
      open(10,file='inputs/sdincoh'//intag(1:outl)//'.dat')
      
cccccccc

      ntotal=1000000    
      ntotal1=200
      ntotal0=100
      
      ktmax=2d0
      
      ptmax=3d0


      q0=4d0
      q2max=q0**2

ccccccccccc
      
      xmin=1d-6
      xmax=0.6d0
      
      lxmin=dlog(xmin)
      lxmax=dlog(xmax)

      nx1=4
      iqtot=10
     
      do 902 ix1=1,nx1+1
         
         lnx1=lxmin+(lxmax-lxmin)*dble(ix1-1)/dble(nx1)
         x=dexp(lnx1)
         
      do 902 iq=1,iqtot

         qtmin=0.01d0
         qtmax=3d0

         lqtmax=dlog(qtmax)
         lqtmin=dlog(qtmin)
         
         lqt=(lqtmax-lqtmin)*(dble(iq))/dble(iqtot)+lqtmin
         qt=dexp(lqt)

         zmin=x
         zmax=0.99d0
         umax=-dlog(zmin)
         umin=-dlog(zmax)
   
         sum=0d0
         
      do n=1,ntotal0


         ru=dble(n)-0.5d0
         ru=ru/dble(ntotal0)
         
         u=umin+(umax-umin)*ru
  
         z=dexp(-u)

         wz=(umax-umin)
 
         
ccccccccc

         q2min=x**2*mp**2/(1d0-z)

         if(q2min.gt.q2max)then
            wt=0d0
            goto 101
         endif

         q2=(qt**2+x**2*mp**2)/(1d0-z)
         wt=1d0
         wt=wt*wz/(1d0-z)
         pt=qt
         
         sum=sum+wt*fem_in(x/z,q2,1)**2*pt**2

 101     continue
         
      enddo

      sum=sum/dfloat(ntotal0)

      sum0=sum

      sum=0d0
      
      do 900 in1=1,ntotal1/2
      do 900 in2=1,ntotal1/5
      do 900 in3=1,ntotal1
            
         ru=(dble(in1)-0.5d0)/dble(ntotal1)*2d0
         
         u=umin+(umax-umin)*ru
         z=dexp(-u)

         wz=(umax-umin)

         rphi=(dble(in2)-0.5d0)/dble(ntotal1)*5d0
         rkt=(dble(in3)-0.5d0)/dble(ntotal1)     
         
         phi=2d0*pi*rphi
         kt=ktmax*rkt
         

         q2min=x**2*mp**2/(1d0-z)

         if(q2min.gt.q2max)then
            wt=0d0
            goto 102
         endif
         
         q2=(qt**2+x**2*mp**2)/(1d0-z)
         pt=qt
         
         qkt=(kt**2+x**2*mp**2)/(1d0-z)

         opac=0d0
         do i1=1,nch
            do i2=1,nch
               call screeningint(i1,i2,kt**2,sc,sc1) 
               opac=opac+sc*pp0(i1)*pp0(i2)/dble(nch)**2
     &              /(gaa(i1)*gaa(i2))
     &              *fe(qkt,i1)
            enddo
         enddo
         
         wt=2d0*pi*ktmax*kt
         wt=wt*wz/(1d0-z)

         pt1=dsqrt(pt**2+kt**2+2d0*dcos(phi)*kt*pt)
         q2_1=(pt1**2+x**2*mp**2)/(1d0-z)

         ptp=pt**2+pt*kt*dcos(phi)

         sum=sum+wt*opac*fem_in(x/z,q2,1)*
     &        fem_in(x/z,q2_1,1)*ptp

 102     continue

 900  enddo

      sum=sum/dfloat(ntotal1)**3*10d0
      
      sum1=0d0


         do 810 n=1,ntotal

            ru=ran2()
            
            u=umin+(umax-umin)*ru
            z=dexp(-u)
    
            ptmax0=ptmax

            wz=(umax-umin)

            rkt1=ran2()
            rkt2=ran2()
            rphi1=ran2()
            rphi2=ran2()
      
            pt=ptmax0*rpt
            kt1=ktmax*rkt1
            kt2=ktmax*rkt2
            phi1=2d0*pi*rphi1
            phi2=2d0*pi*rphi2

            q2min=x**2*mp**2/(1d0-z)

            if(q2min.gt.q2max)then
               wt=0d0
               goto 103
            endif
            
            q2=(qt**2+x**2*mp**2)/(1d0-z)

            pt=qt
 
            wt=ktmax**2*4d0*pi**2
            wt=wt*kt1*kt2
            wt=wt*wz/(1d0-z)
       
            ktc=kt1**2+kt2**2+2d0*kt1*kt2*dcos(phi1-phi2)
            ktc=(ktc+x**2*mp**2)/(1d0-z)

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

         q2_1=(pt1**2+x**2*mp**2)/(1d0-z)
         q2_2=(pt2**2+x**2*mp**2)/(1d0-z)

         ptp=pt**2+pt*kt1*dcos(phi1)+pt*kt2*dcos(phi2)+kt1*
     &        kt2*dcos(phi2-phi1)
         
         sum1=sum1+wt*opac*fem_in(x/z,q2_1,1)*
     &        fem_in(x/z,q2_2,1)
     &        *ptp

 103     continue
         
 810  enddo

      sum1=sum1/dfloat(ntotal)

      surv=1d0+2d0*sum/sum0+sum1/sum0

c      print*,x,qt,1d0+2d0*sum/sum0+sum1/sum0
c      print*,sum0,sum,sum1

       if(sum0.eq.0d0.or.surv.lt.0d0)then
          write(10,*)dlog(x),dlog(qt),0d0
c          print*,x,qt,0d0
       else
          write(10,*)dlog(x),dlog(qt),1d0+2d0*sum/sum0+sum1/sum0
c          print*,x,qt,surv
      endif
      
 902  enddo
     
 887  close(10)

      
cccccccccccccccccc


      return
      end
      
      function fem_in(xb,qsq,i)
      implicit double precision(a-z)
      integer i

      include 'mp.f'

      mx=qsq*(1d0-xb)/xb+mp**2
      mx=dsqrt(mx)
      call F1F2(.true.,xb,qsq,mx,f1,f2)

      fem_in=dsqrt(f2)/qsq
      
      return
      end
