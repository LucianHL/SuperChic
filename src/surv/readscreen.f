ccc   read in screened amplitude from file
      subroutine readscreen
      implicit none
      double precision ksqma,inck,lgksq,lginck,ksqmin,
     &     ksq
      integer ns,i1,i2,ib,outl

      include 'nchan.f'
      include 'screenamp.f'
      include 'intag.f'

      ns=900
      ksqma=8.2d0
      inck=ksqma/dble(ns)
      ksqmin=0.001d0
      lginck=(dlog(ksqma/ksqmin))/dble(ns)

      call length(intag,outl)

      open(40,file='inputs/screening'//intag(1:outl)//'.dat')

      do ib=0,ns+1

         ksq=dble(ib-1)*inck
         lgksq=ksq

         if(ib.eq.0)then
            ksq=0d0
            lgksq=0d0
         else
            lgksq=dble(ib-1)*lginck+dlog(ksqmin)
            ksq=dexp(lgksq)
         endif

         do i1=1,nch
            do i2=1,nch

               read(40,*)sca(i1,i2,ib,1),sca(i1,i2,ib,2)
     &              ,sca1(i1,i2,ib,2)

               sca1(i1,i2,ib,1)=sca(i1,i2,ib,1)

            enddo
         enddo

      enddo

      close(40)

      return
      end


  
