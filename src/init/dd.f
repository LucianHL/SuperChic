      subroutine dd
      implicit double precision(a-y)
      integer nphi,nkt,nbt
      integer iphi,ikt,ibt,i1,i2,outl

      include 'nchan.f'
      include 'survpars.f'
      include 'pi.f'
      include 'vars.f'
      include 'intag.f'

      print*,'Calculating S^2 (Evolution + Evolution)...'
      
      call length(intag,outl)
      open(10,file='inputs/dd'//intag(1:outl)//'.dat')

      nphi=10

      ktmax=5d0
      nkt=20000

      btmax=100d0
      nbt=2000

      crossb=0d0
      cross=0d0

      sige=sigo*dexp(dlog(rts)*2d0*ep)

      hkt=ktmax/dble(nkt)
      hbt=btmax/dble(nbt)
      hphi=2d0*pi/dble(nphi)

      
      sum=0d0

 999  do ibt=1,nbt

         bt=(dble(ibt)-0.5d0)*hbt

         sigmab=0d0
         sigma=0d0

         opac=0d0
         do i1=1,nch
            do i2=1,nch
               call opacityint(i1,i2,bt,fr,fr1)
               tet=-0.5d0*sige*fr*gaa(i1)*gaa(i2)
               opac=opac+dexp(tet)
     &              *pp0(i1)*pp0(i2)/dble(nch)**2
            enddo
         enddo

         opac=opac**2
         
      do ikt=1,nkt

         kt=(dble(ikt)-0.5d0)*hkt
         
         wt=besj0(bt*kt)
         wt=wt/2d0/pi*hkt*kt         

         call F1F2el(kt**2,f1,f2)

         sigmab=sigmab+wt*f2

              
      enddo

      crossb=crossb+sigmab*hbt*bt
      cross=cross+sigmab*opac*hbt*bt

      enddo
 
c      print*,cross/crossb
      write(10,*)cross/crossb
     

 887  close(10)

      return
      end
