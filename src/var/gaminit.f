      subroutine gaminit
      implicit none
      integer i,j,k,l
      
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
            g5(i,j)=0d0
            do k=1,4
               gmatrix(i,j,k)=0d0
            enddo
         enddo
      enddo

      do i=1,4
         do j=1,4
            do k=1,4
               do l=1,4
                  e_(i,j,k,l)=0d0
               enddo
            enddo
         enddo
      enddo
      
      e_(4,1,2,3)=1d0
      e_(4,1,3,2)=-1d0
      e_(4,2,1,3)=-1d0
      e_(4,2,3,1)=1d0
      e_(4,3,2,1)=-1d0
      e_(4,3,1,2)=1d0

      e_(1,4,3,2)=1d0
      e_(1,4,2,3)=-1d0
      e_(1,3,4,2)=-1d0
      e_(1,3,2,4)=1d0
      e_(1,2,3,4)=-1d0
      e_(1,2,4,3)=1d0

      e_(2,1,3,4)=1d0
      e_(2,1,4,3)=-1d0
      e_(2,3,1,4)=-1d0
      e_(2,3,4,1)=1d0
      e_(2,4,3,1)=-1d0
      e_(2,4,1,3)=1d0

      e_(3,1,4,2)=1d0
      e_(3,1,2,4)=-1d0
      e_(3,2,4,1)=-1d0
      e_(3,2,1,4)=1d0
      e_(3,4,2,1)=1d0
      e_(3,4,1,2)=-1d0

      g5(4,2)=1d0
      g5(1,3)=1d0
      g5(2,4)=1d0
      g5(3,1)=1d0
      
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

      subroutine gaminit_comb
      implicit none
      integer i,j,k,i1,j1,k1,k2,l,m
      complex*16 zt
      
      include 'gmatrices.f'
      include 'gmatrices_comb.f'
      include 'zi.f'

c      do i=1,4
c         do j=1,4
c            do k=1,4
c               gmatrix5_1(i,j,k)=0d0
c            do i1=1,4
c               gmatrix5_1(i,j,k)=gmatrix5_1(i,j,k)
c     &              +g5(j,i1)*gmatrix(i,i1,k)
c            enddo
c         enddo
c      enddo
c      enddo

      do i=1,4
      do j=1,4
      do k=1,4
         do i1=1,4
            do j1=1,4
               epgam(i,j,k,i1,j1)=0d0
               do k1=1,4
                  zt=e_(i,j,k,k1)*gmatrix(k1,i1,j1)
                  if(k1.lt.4)zt=-zt
                  epgam(i,j,k,i1,j1)=epgam(i,j,k,i1,j1)+zt
               enddo
            enddo
         enddo
      enddo
      enddo
      enddo
         
      do i=1,4
      do j=1,4
      do i1=1,4
      do j1=1,4  
         gmatrix_2(i,j,i1,j1)=0d0
         do k1=1,4
            gmatrix_2(i,j,i1,j1)=gmatrix_2(i,j,i1,j1)
     &           +gmatrix(i,i1,k1)*gmatrix(j,k1,j1)
            
         enddo
      enddo
      enddo
      enddo
      enddo

      
      
      do i=1,4
      do j=1,4
      do k=1,4
      do i1=1,4
      do j1=1,4  
         gmatrix_3(i,j,k,i1,j1)=0d0
         do k1=1,4
               gmatrix_3(i,j,k,i1,j1)=gmatrix_3(i,j,k,i1,j1)
     &           +gmatrix(i,i1,k1)*gmatrix_2(j,k,k1,j1)
         enddo
      enddo
      enddo
      enddo
      enddo
      enddo

c      do i=1,4
c      do j=1,4
c      do k=1,4
c      do l=1,4    
c      do i1=1,4
c      do j1=1,4  
c         gmatrix_4(i,j,k,l,i1,j1)=0d0
c         do k1=1,4
c            gmatrix_4(i,j,k,l,i1,j1)=gmatrix_4(i,j,k,l,i1,j1)
c     &           +gmatrix(i,i1,k1)*gmatrix_3(j,k,l,k1,j1)
c         enddo
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo

c      do i=1,4
c      do j=1,4
c      do k=1,4
c      do l=1,4
c      do m=1,4     
c      do i1=1,4
c      do j1=1,4  
c         gmatrix_5(i,j,k,l,m,i1,j1)=0d0
c         do k1=1,4
c            gmatrix_5(i,j,k,l,m,i1,j1)=gmatrix_5(i,j,k,l,m,i1,j1)
c     &           +gmatrix(i,i1,k1)*gmatrix_4(j,k,l,m,k1,j1)
c         enddo
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo
c      enddo
c
               
      return
      end
