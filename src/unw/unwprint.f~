ccc   randomizes order of VEGAS unweighted events and
ccc   prints nev events to record
      subroutine unwprint
      implicit double precision(a-y)
      integer i,j,k,l,m,n
      integer evfill(2000000)

      include 'pdg.f'
      include 'unweighted.f'
      include 'mom.f'
      include 'proc.f'
      include 'decay.f'
      include 'record.f'
      include 'hepevt.f'
      include 'leshouches.f'
      include 'ewpars.f'
      include 'mp.f'
      include 'vars.f'
      include 'beam.f'
      include 'mion.f'
      include 'x.f'
      include 'vertex.f'
      include 'pvert.f'
      include 'rech.f'
      include 'ion.f'

      
      do i=1,evnum
         evfill(i)=1
      enddo

      range=dble(evnum)

      do i=1,nev

 555     r=ran2()

         j=nint(r*range)
         if(dble(j).lt.r*range)j=j+1
         if(evfill(j).eq.0)goto 555
         evfill(j)=0
         
ccccccccccccccccccccccccccccccccccccccccccccccc
ccccc HepMC
ccccccccccccccccccccccccccccccccccccccccccccccc

         
         if(erech)then

            
           do k=3,nup+2
               do l=1,4
                  pup(l,k)=evrec(j,k,l)
               enddo
               pup(5,k)=dsqrt(dabs(pup(4,k)**2-pup(3,k)**2
     &              -pup(2,k)**2-pup(1,k)**2))
            enddo

            do k=1,nup+2
               idup(k)=pdgid(k)
            enddo

            if(beam.eq.'el')then
               pup(5,3)=me
               pup(5,4)=me
            elseif(beam.eq.'prot')then
               pup(5,3)=mp
               pup(5,4)=mp
            elseif(beam.eq.'ion')then
               pup(5,3)=dsqrt(dabs(pup(4,3)**2-pup(3,3)**2
     &              -pup(2,3)**2-pup(1,3)**2))
               pup(5,4)=dsqrt(dabs(pup(4,4)**2-pup(3,4)**2
     &              -pup(2,4)**2-pup(1,4)**2))
            elseif(beam.eq.'ionp')then
               pup(5,3)=mp
               pup(5,4)=dsqrt(dabs(pup(4,4)**2-pup(3,4)**2
     &              -pup(2,4)**2-pup(1,4)**2))
            endif

cccccc incoming proton vertices

            call vertinc(1,barv,vert,orph,nout,momv,statv,pdgv,massv)
            call vertinc(2,barv,vert,orph,nout,momv,statv,pdgv,massv)

cccccccccc Central production vertex

            call vertcent(3,barv,vert,orph,nout,momv,statv,pdgv,massv)

ccccccccc Outgoing vertices

            nvert=3
            mv1=0
            mv2=0
            
            do m=6,nup+2
               
               vert(nvert+1)=nvert+1
               orph(nvert+1)=0
               nout(nvert+1)=0
               
               if(mothup(1,m).eq.mv1.and.mothup(2,m).eq.mv2)then
                  ip=ip+1
                  call passign(m,ip,nvert,nout,barv,pdgv,momv,
     &                 massv,statv)                  
               else
                  mv1=mothup(1,m)
                  mv2=mothup(2,m)
                  ip=1
                  call passign(m,ip,nvert+1,nout,barv,pdgv,momv,
     &                 massv,statv)
                  nvert=nvert+1
               endif
                  
            enddo

ccccccccc            

            call vertidscan(nvert,orph,nout,barv)   ! Assign vertex barcodes
                        
            scalup=mx
            aqcdup=alphas(mx**2)
            aqedup=alpha
            
            write(45,51)'E',i,0,scalup,aqcdup,aqedup,proc,0
     &           ,nvert,1,2,0,0
            write(45,55)'U GEV CM'
            write(45,52)'F',nfl1,nfl2,x1,x2,scalup,0d0,0d0,0,0
            
            do n=1,nvert

               write(45,53)'V',vert(n),0,0d0,0d0,0d0,0d0,
     &           orph(n),nout(n)

               do m=1,orph(n)+nout(n)
                  write(45,54)'P',barv(n,m,1),pdgv(n,m),momv(n,m,1)
     &                 ,momv(n,m,2),momv(n,m,3),momv(n,m,4),massv(n,m),
     &                 statv(n,m),0d0,0d0,barv(n,m,2),0,0
               enddo
               
            enddo
            

            write(45,*)''
            
            goto 500


            
         endif

 51      format(1a,1x,i8,1x,i4,1x,E16.9,1x,E16.9,1x,E16.9,1x,i4,1x,i1
     &        ,1x,i4,1x,i1,1x,i1,1x,i1,1x,i1,1x)
 52      format(1a,1x,i4,1x,i4,1x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x
     &        ,E16.9,1x,i1,1x,i1)
 53      format(1a,1x,i2,1x,i1,1x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x
     &        ,i2,1x,i2)
 54      format(1a,1x,i2,1x,i10,1x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x
     &        ,E16.9,1x,i2,1x,E16.9,1x,E16.9,1x,i4,1x,i1,1x,i1)
 55      format(8a)

ccccccccccccccccccccccccccccccccccccccccccccccc
ccccc Les Houches
ccccccccccccccccccccccccccccccccccccccccccccccc

         if(erec.eq.'lhe')then
            
            do k=1,nup+2
               idup(k)=pdgid(k)
            enddo
            
           do k=3,nup+2
               do l=1,4
                  pup(l,k)=evrec(j,k,l)
               enddo
               pup(5,k)=dsqrt(dabs(pup(4,k)**2-pup(3,k)**2
     &              -pup(2,k)**2-pup(1,k)**2))
            enddo

            if(beam.eq.'el')then
               pup(5,3)=me
               pup(5,4)=me
            elseif(beam.eq.'prot')then
               pup(5,3)=mp
               pup(5,4)=mp
            elseif(beam.eq.'ion')then
               pup(5,3)=dsqrt(dabs(pup(4,3)**2-pup(3,3)**2
     &              -pup(2,3)**2-pup(1,3)**2))
               pup(5,4)=dsqrt(dabs(pup(4,4)**2-pup(3,4)**2
     &              -pup(2,4)**2-pup(1,4)**2))
            elseif(beam.eq.'ionp')then
               pup(5,3)=mp
               pup(5,4)=dsqrt(dabs(pup(4,4)**2-pup(3,4)**2
     &              -pup(2,4)**2-pup(1,4)**2))
            endif

            if(i.eq.1)then
               write(45,*)'<LesHouchesEvents version="1.0">'
               write(45,*)'<header>'
               call headerlhe
               write(45,*)'</header>'
               write(45,*)'<init>'
               if(beam.eq.'el'.or.beam.eq.'prot')then
                  write(45,302)idup(1),idup(2),pup(4,1),pup(4,2)
     &                 ,pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2)
     &                 ,idwtup,nprup
               elseif(beam.eq.'ion')then
                  write(45,202)idup(1),idup(2),pup(4,1)*an,pup(4,2)*an
     &                 ,pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2)
     &                 ,idwtup,nprup
               elseif(beam.eq.'ionp')then
                  write(45,202)idup(1),idup(2),pup(4,1),pup(4,2)*an
     &                 ,pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2)
     &                 ,idwtup,nprup
               endif
               write(45,312)xsecup(1),xerrup(1),xmaxup(1),1
               write(45,*)'</init>'
            endif

 202        format(5x,i10,2x,i10,2x,E16.9,2x,E16.9,2x,i5,2x,i5,2x,i6,2x
     &           ,i6,2x,i1,2x,i1)
            
 302        format(5x,i4,2x,i4,2x,E16.9,2x,E16.9,2x,i5,2x,i5,2x,i6,2x,i6
     &           ,2x,i1,2x,i1)

 312        format(E16.9,2x,E16.9,2x,E16.9,2x,i1)
            
            scalup=mx
            aqcdup=alphas(mx**2)
            
            write(45,*)'<event>'
            write(45,304)nup,idprup,xwgtup,scalup,aqedup,aqcdup
            if(beam.eq.'prot'.or.beam.eq.'el')then
               do m=3,nup+2
                  write(45,303)idup(m),istup(m),mothup(1,m),
     &                 mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
     &                 ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
     &                 ,spinup(m)
               enddo
            elseif(beam.eq.'ion'.or.beam.eq.'ionp')then
               do m=3,nup+2
                  write(45,203)idup(m),istup(m),mothup(1,m),
     &                 mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
     &                 ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
     &                 ,spinup(m)
               enddo
            endif
            write(45,*)'</event>'
            
c            write(45,*)''

 304        format(i2,4x,i1,3x,F2.0,3x,E16.9,3x,E16.9,3x,E16.9)
            
         endif
         
ccccccccccccccccccccccccccccccccccccccccccccccc
cccc  HEPEVT
ccccccccccccccccccccccccccccccccccccccccccccccc
         
         if(erec.eq.'hepevt')then
            
            nevhep=nev
            
            do k=1,nhep
               idhep(k)=pdgid(k)
            enddo
            
           do k=3,nhep
              do l=1,4
                 phep(l,k)=evrec(j,k,l)
              enddo
              phep(5,k)=dsqrt(dabs(phep(4,k)**2-phep(3,k)**2
     &             -phep(2,k)**2-phep(1,k)**2))
           enddo

            if(beam.eq.'el')then
               phep(5,3)=me
               phep(5,4)=me
            elseif(beam.eq.'prot')then
               phep(5,3)=mp
               phep(5,4)=mp
            elseif(beam.eq.'ion')then
               phep(5,2)=dsqrt(dabs(phep(4,3)**2-phep(3,3)**2
     &              -phep(2,3)**2-phep(1,3)**2))
               phep(5,4)=dsqrt(dabs(phep(4,4)**2-phep(3,4)**2
     &              -phep(2,4)**2-phep(1,4)**2))
            elseif(beam.eq.'ionp')then
               phep(5,3)=mp
               phep(5,4)=dsqrt(dabs(phep(4,4)**2-phep(3,4)**2
     &              -phep(2,4)**2-phep(1,4)**2))
            endif

            do k=1,2
               do m=5,nhep
                  jmohep(k,m)=mothup(k,m)
               enddo
            enddo
               
            write(45,*)i


            if(beam.eq.'el'.or.beam.eq.'prot')then
            do m=1,nhep
               write(45,300)m,idhep(m),isthep(m),jmohep(1,m),
     &              jmohep(2,m),jdahep(1,m),jdahep(2,m),
     &              phep(1,m),phep(2,m),phep(3,m),phep(4,m)
     &              ,phep(5,m),vhep(1,m),vhep(2,m),vhep(3,m),vhep(4,m)
            enddo
            elseif(beam.eq.'ion')then
            do m=1,nhep
               write(45,200)m,idhep(m),isthep(m),jmohep(1,m),
     &              jmohep(2,m),jdahep(1,m),jdahep(2,m),
     &              phep(1,m),phep(2,m),phep(3,m),phep(4,m)
     &              ,phep(5,m),vhep(1,m),vhep(2,m),vhep(3,m),vhep(4,m)
            enddo
            elseif(beam.eq.'ionp')then
            do m=1,nhep
            write(45,300)m,idhep(m),isthep(m),jmohep(1,m),
     &              jmohep(2,m),jdahep(1,m),jdahep(2,m),
     &              phep(1,m),phep(2,m),phep(3,m),phep(4,m)
     &              ,phep(5,m),vhep(1,m),vhep(2,m),vhep(3,m),vhep(4,m)
            enddo
            endif
               
            write(45,*)''
            
         endif
               
 500  enddo

      if(erec.eq.'lhe')then
         if(erech)then
         else
            write(45,*)'</LesHouchesEvents>'
         endif
      endif

 200  format(i4,1x,i10,1x,i4,1x,i4,1x,i4,1x,i4,1x,i4,1x,E18.11,1x,E18.11
     &,1x,E18.11,1x,E18.11,1x,E18.11,1x,E18.11,1x,E18.11,1x,E18.11,1x
     &,E18.11)

 300  format(i4,1x,i10,1x,i4,1x,i4,1x,i4,1x,i4,1x,i4,1x,E18.11,1x,E18.11
     &,1x,E18.11,1x,E18.11,1x,E18.11,1x,E18.11,1x,E18.11,1x,E18.11,1x
     &     ,E18.11)

 301  format(i5,1x,i4,1x,i8,1x,i4,1x,i4,1x,i4,1x,i4,1x,E16.9,1x,
     &E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x,E16.9)

 303  format(7x,i10,1x,i8,1x,i4,1x,i4,1x,i4,1x,i4,1x,E16.9,1x,
     &     E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x,F2.0,1x,F2.0)

 203  format(7x,i10,1x,i8,1x,i4,1x,i4,1x,i4,1x,i4,1x,E16.9,1x,
     &E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x,F2.0,1x,F2.0)

      return
      end
