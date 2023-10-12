ccc   reads Sudakov factor from file
      subroutine calcs2diss
      implicit none
      integer outl,ix,iqt,nx,nqt

      include 's2_diss.f'
      include 'intag.f'
      
      nx=5
      nqt=10

      call length(intag,outl)

      open(42,file='inputs/sdcoh'//intag(1:outl)//'.dat')
      open(43,file='inputs/sdincoh'//intag(1:outl)//'.dat')

      do ix=1,nx
         do iqt=1,nqt
            read(42,*)scearr(1,ix,iqt),scearr(2,ix,iqt),scearr(3,ix,iqt)
            read(43,*)siearr(1,ix,iqt),siearr(2,ix,iqt),siearr(3,ix,iqt)
         enddo
      enddo
         
      close(42)
      close(43)

      open(44,file='inputs/dd'//intag(1:outl)//'.dat')

      read(44,*)s2dd

      close(44)
      
      return
      end
