      subroutine passign(m,ip,nv,nout,barv,pdgv,momv,massv,statv)
      implicit double precision(a-y)
      integer m,n
      
      include 'vertex.f'
      include 'pvert.f'
      include 'leshouches.f'

      nout(nv)=nout(nv)+1
      barv(nv,ip,1)=m
      barv(nv,ip,2)=0
      pdgv(nv,ip)=idup(m)
      do n=1,4
         momv(nv,ip,n)=pup(n,m)
      enddo
      massv(nv,ip)=pup(5,m)
      statv(nv,ip)=istup(m)
      
      return
      end

      subroutine vertidscan(nvert,orph,nout,barv)
      implicit double precision(a-y)
      integer m,k,n

      include 'vertex.f'
      include 'pvert.f'
      include 'leshouches.f'
      
      do m=6,nup+2
         
         mv1=mothup(1,m)
         
         do n=4,nvert
            do k=1,orph(n)+nout(n)

               if(barv(n,k,1).eq.m)nv=n
               
            enddo
         enddo
         
         do n=4,nvert
            do k=1,orph(n)+nout(n)
               
               if(barv(n,k,1).eq.mv1)then
                  barv(n,k,2)=nv
               endif
               
            enddo
         enddo
         
      enddo
      

      return
      end

      subroutine vertinc(i,barv,vert,orph,nout,momv,statv,pdgv,massv)
      implicit double precision(a-y)
      integer i,m

      include 'vertex.f'
      include 'pvert.f'
      include 'leshouches.f'
      include 'mom.f'
      include 'mp.f'
      
      vert(i)=i
      orph(i)=1
      nout(i)=2
      
      barv(i,1,1)=i
      barv(i,2,1)=i+2
      barv(i,1,2)=vert(i)
      barv(i,2,2)=vert(i)
      pdgv(i,1)=idup(i)
      pdgv(i,2)=idup(i+2)
      
      do m=1,4
         momv(i,1,m)=q(m,i)
      enddo
      massv(i,1)=mp
      statv(i,1)=3
      statv(i,2)=1
      statv(i,3)=2

      do m=1,4
         momv(i,2,m)=pup(m,i+2)
         momv(i,3,m)=q(m,i)-pup(m,i+2)
      enddo
c     massv(i,2)=mp
      massv(i,2)=dsqrt(pup(4,i+2)**2-pup(3,i+2)**2
     &     -pup(2,i+2)**2-pup(1,i+2)**2)
      massv(i,3)=0d0
      
      barv(i,3,1)=-i
      barv(i,3,2)=0
      pdgv(i,3)=22
      
      return
      end

      subroutine vertcent(i,barv,vert,orph,nout,momv,statv,pdgv,massv)
      implicit double precision(a-y)
      integer i,m

      include 'vertex.f'
      include 'pvert.f'
      include 'leshouches.f'
      include 'mom.f'
      include 'mp.f'

      vert(i)=i
      orph(i)=0
      nout(i)=1
      
      barv(i,1,1)=5
      barv(i,1,2)=1
      pdgv(i,1)=idup(5)
      do m=1,4
         momv(i,1,m)=pup(m,5)
      enddo
      massv(i,1)=dsqrt(pup(4,5)**2-pup(3,5)**2
     &     -pup(2,5)**2-pup(1,5)**2)

      statv(i,1)=1
      
      return
      end
