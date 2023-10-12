      subroutine dd_test
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

  
      print*,'Calculating S^2 (Evolution + Evolution)...TEST'

cccccccc

      ntotal=10000000
      ntotal1=500
    
      ktmax=4d0
 


      norm0=0d0
      do i1=1,nch
         do i2=1,nch
            norm0=norm0+pp0(i1)*pp0(i2)/dble(nch)**2
c     &           *gaa(i1)
         enddo
      enddo

      sum0=norm0
      sum0=1d0


      sum=0d0

      do 900 in1=1,ntotal1
      do 900 in2=1,ntotal1

         rphi=(dble(in1)-0.5d0)/dble(ntotal1)
         rkt=(dble(in2)-0.5d0)/dble(ntotal1)
         
         phi=2d0*pi*rphi
         kt=ktmax*rkt

         ktc=kt**2
               
         opac=0d0
         do i1=1,nch
            do i2=1,nch
               call screeningint(i1,i2,kt**2,sc,sc1) 
               opac=opac+sc*pp0(i1)*pp0(i2)/dble(nch)**2
     &              /(gaa(i1)*gaa(i2))
     &              *fe(ktc,i1)**2
            enddo
         enddo
         
         wt=2d0*pi*ktmax*kt
         sum=sum+wt*opac

 900  enddo

      sum=sum/dfloat(ntotal1)**2

      sum1=0d0
      
      do 810 n=1,ntotal
         
         rkt1=ran2()
         rkt2=ran2()
         rphi1=ran2()
         rphi2=ran2()   
         
         kt1=ktmax*rkt1
         kt2=ktmax*rkt2
         phi1=2d0*pi*rphi1
         phi2=2d0*pi*rphi2

         wt=ktmax**2*4d0*pi**2
         wt=wt*kt1*kt2
         
         ktc=kt1**2+kt2**2+2d0*kt1*kt2*dcos(phi1-phi2)
 
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
     &                    *fe(ktc,i1)**2
                  enddo
               enddo
            enddo
         enddo
         
         sum1=sum1+wt*opac

 810  enddo

      sum1=sum1/dfloat(ntotal)

c      sum1=sum1/4d0/pi**2
c      sum=sum/2d0/pi


      print*,1d0+2d0*sum/sum0+sum1/sum0

     
 887  close(10)

      
cccccccccccccccccc


      return
      end

