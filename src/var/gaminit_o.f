      subroutine gaminit
      implicit double precision (a-z)
      integer i,j,k
      
      include 'gmatrices.f'
      include 'zi.f'
      

      do i=1,4
         do j=1,4
            ident(i,j)=0d0
            d_(i,j)=0d0
            if(i.eq.j)then
               if(i.lt.4)d_(i,j)=-1d0
               if(i.eq.4)d_(i,j)=1d0
            endif
            do k=1,4
               gmatrix(i,j,k)=0d0
            enddo
         enddo
      enddo

      gmatrix(4,1,1)=1d0
      gmatrix(4,2,2)=1d0
      gmatrix(4,3,3)=-1d0
      gmatrix(4,4,4)=-1d0

      gmatrix(1,1,4)=1d0
      gmatrix(1,2,3)=1d0
      gmatrix(1,3,2)=-1d0
      gmatrix(1,4,1)=-1d0

      gmatrix(2,1,4)=-zi
      gmatrix(2,2,3)=zi
      gmatrix(2,3,2)=zi
      gmatrix(2,4,1)=-zi

      gmatrix(3,1,3)=1d0
      gmatrix(3,2,4)=-1d0
      gmatrix(3,3,1)=-1d0
      gmatrix(3,4,2)=1d0

      ident(1,1)=1d0
      ident(2,2)=1d0
      ident(3,3)=1d0
      ident(4,4)=1d0
      
      return
      end
