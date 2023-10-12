ccc   randomizes order of VEGAS unweighted events and
ccc   prints nev events to record
      subroutine unwprintq
      implicit double precision(a-y)
      integer i,j,k,l,m,n,bar
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
      include 'diff.f'  
      
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

            do k=1,nup+2
               idup(k)=pdgid(k)
            enddo
            
           do k=3,nup+2
               do l=1,4
                  pup(l,k)=evrec(j,k,l)
               enddo
               pup(5,k)=dsqrt(dabs(pup(4,k)**2-pup(3,k)**2
     &              -pup(2,k)**2-pup(1,k)**2))
               if(idup(k).eq.22)pup(5,k)=0d0
            enddo

           if(beam.eq.'prot')then
              if(pup(5,3).gt.1d0)then
                 idup(3)=92
              endif
              if(pup(5,4).gt.1d0)then
                 idup(4)=92
              endif
           endif

c$$$            if(beam.eq.'prot')then
c$$$              if(pup(5,3).gt.1d0)idup(3)=92
c$$$              if(pup(5,4).gt.1d0)idup(4)=92
c$$$           endif
c$$$
c$$$cccccc incoming proton vertices
c$$$
c$$$            call vertinc(1,barv,vert,orph,nout,momv,statv,pdgv,massv)
c$$$            call vertinc(2,barv,vert,orph,nout,momv,statv,pdgv,massv)
c$$$
c$$$cccccccccc Central production vertex
c$$$
c$$$            call vertcent(3,barv,vert,orph,nout,momv,statv,pdgv,massv)
c$$$
c$$$ccccccccc Outgoing vertices
c$$$
c$$$            nvert=3
c$$$            mv1=0
c$$$            mv2=0
c$$$            
c$$$            do m=6,nup+2
c$$$               
c$$$               vert(nvert+1)=nvert+1
c$$$               orph(nvert+1)=0
c$$$               nout(nvert+1)=0
c$$$               
c$$$               if(mothup(1,m).eq.mv1.and.mothup(2,m).eq.mv2)then
c$$$                  ip=ip+1
c$$$                  call passign(m,ip,nvert,nout,barv,pdgv,momv,
c$$$     &                 massv,statv)                  
c$$$               else
c$$$                  mv1=mothup(1,m)
c$$$                  mv2=mothup(2,m)
c$$$                  ip=1
c$$$                  call passign(m,ip,nvert+1,nout,barv,pdgv,momv,
c$$$     &                 cmassv,statv)
c$$$                  nvert=nvert+1
c$$$               endif
c$$$                  
c$$$            enddo
c$$$
c$$$ccccccccc            
c$$$
c$$$            call vertidscan(nvert,orph,nout,barv)   ! Assign vertex barcodes
                        
c$$$            scalup=mx
c$$$            aqcdup=alphas(mx**2)
c$$$            aqedup=alpha
           
           if(diff.eq.'dd')then
              
              qsq1=(pup(4,5)-pup(4,3))**2-(pup(3,5)-pup(3,3))**2
     &             -(pup(2,5)-pup(2,3))**2-(pup(1,5)-pup(1,3))**2
              qsq1=-qsq1
              
              qsq2=(pup(4,6)-pup(4,4))**2-(pup(3,6)-pup(3,4))**2
     &             -(pup(2,6)-pup(2,4))**2-(pup(1,6)-pup(1,4))**2
              qsq2=-qsq2
              
              if(qsq1.gt.qsq2)then
                 scalup=dsqrt(qsq1)
                 aqcdup=alphas(scalup**2)
              else
                 scalup=dsqrt(qsq2)
              endif
              
           elseif(diff.eq.'sd')then
              
              qsq=(pup(4,5)-pup(4,4))**2-(pup(3,5)-pup(3,4))**2
     &             -(pup(2,5)-pup(2,4))**2-(pup(1,5)-pup(1,4))**2
              qsq=-qsq
               
              scalup=dsqrt(qsq)
              
           elseif(diff.eq.'el')then

              scalup=dsqrt((pup(4,5)+pup(4,6))**2-(pup(3,5)+pup(3,6))**2
     &             -(pup(2,5)+pup(2,6))**2-(pup(1,5)+pup(1,6))**2)
              
           endif
           
           aqcdup=alphas(scalup**2)
           
           write(45,51)'E',i,0,scalup,aqcdup,aqedup,proc,0
     &          ,nvert,1,2,0,0
           write(45,55)'U GEV CM'
           write(45,52)'F',nfl1,nfl2,x1,x2,scalup,0d0,0d0,0,0

           if(diff.eq.'el')then

              write(45,53)'V',1,0,0d0,0d0,0d0,0d0,
     &             2,2

              istup(3)=3
              istup(4)=3

              do m=3,nup+1
                 bar=1
                 if(m.gt.4)bar=0
                 write(45,54)'P',m-2,idup(m),pup(1,m)
     &                ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),
     &                istup(m),0d0,0d0,bar,0,0
              enddo

c$$$                            do m=3,nup+1
c$$$                     write(45,303)idup(m),istup(m),mothup(1,m),
c$$$     &                    mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
c$$$     &                    ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
c$$$     &                 ,spinup(m)
c$$$                  enddo    
              
           elseif(diff.eq.'dd')then
              
              write(45,53)'V',1,0,0d0,0d0,0d0,0d0,
     &             1,2

              istup(3)=3
              istup(4)=3

              massgam=(pup(4,5)-pup(4,3))**2-(pup(3,5)-pup(3,3))**2
     &             -(pup(2,5)-pup(2,3))**2-(pup(1,5)-pup(1,3))**2
              massgam=-dsqrt(-massgam)
              write(45,54)'P',1,idup(3),pup(1,3)
     &                ,pup(2,3),pup(3,3),pup(4,3),pup(5,3),
     &             istup(3),0d0,0d0,1,0,0
              write(45,54)'P',2,idup(5),pup(1,5)
     &                ,pup(2,5),pup(3,5),pup(4,5),pup(5,5),
     &             istup(5),0d0,0d0,0,0,0
              write(45,54)'P',3,22,pup(1,3)-pup(1,5)
     &             ,pup(2,3)-pup(2,5),pup(3,3)-pup(3,5)
     &             ,pup(4,3)-pup(4,5),massgam,
     &                2,0d0,0d0,3,0,0

              
              write(45,53)'V',2,0,0d0,0d0,0d0,0d0,
     &             1,2


              massgam=(pup(4,6)-pup(4,4))**2-(pup(3,6)-pup(3,4))**2
     &             -(pup(2,6)-pup(2,4))**2-(pup(1,6)-pup(1,4))**2
              massgam=-dsqrt(-massgam)
              write(45,54)'P',4,idup(4),pup(1,4)
     &                ,pup(2,4),pup(3,4),pup(4,4),pup(5,4),
     &             istup(4),0d0,0d0,2,0,0
              write(45,54)'P',5,idup(6),pup(1,6)
     &                ,pup(2,6),pup(3,6),pup(4,6),pup(5,6),
     &             istup(6),0d0,0d0,0,0,0
              write(45,54)'P',6,22,pup(1,4)-pup(1,6)
     &             ,pup(2,4)-pup(2,6),pup(3,4)-pup(3,6)
     &             ,pup(4,4)-pup(4,6),massgam,
     &             2,0d0,0d0,3,0,0


              write(45,53)'V',3,0,0d0,0d0,0d0,0d0,
     &             0,2
              massgam=(pup(4,5)-pup(4,3))**2-(pup(3,5)-pup(3,3))**2
     &             -(pup(2,5)-pup(2,3))**2-(pup(1,5)-pup(1,3))**2
              massgam=-dsqrt(-massgam)
              write(45,54)'P',3,22,pup(1,3)-pup(1,5)
     &             ,pup(2,3)-pup(2,5),pup(3,3)-pup(3,5)
     &             ,pup(4,3)-pup(4,5),massgam,
     &             2,0d0,0d0,3,0,0
              massgam=(pup(4,6)-pup(4,4))**2-(pup(3,6)-pup(3,4))**2
     &             -(pup(2,6)-pup(2,4))**2-(pup(1,6)-pup(1,4))**2
              massgam=-dsqrt(-massgam)
              write(45,54)'P',6,22,pup(1,4)-pup(1,6)
     &             ,pup(2,4)-pup(2,6),pup(3,4)-pup(3,6)
     &             ,pup(4,4)-pup(4,6),massgam,
     &             2,0d0,0d0,3,0,0
              write(45,54)'P',7,idup(7),pup(1,7)
     &                ,pup(2,7),pup(3,7),pup(4,7),pup(5,7),
     &             istup(7),0d0,0d0,0,0,0
              write(45,54)'P',8,idup(8),pup(1,8)
     &                ,pup(2,8),pup(3,8),pup(4,8),pup(5,8),
     &             istup(8),0d0,0d0,0,0,0

              
           elseif(diff.eq.'sd'.or.diff.eq.'sda'.or.diff.eq.'sdb')then

              istup(3)=3
              istup(4)=3
              
             write(45,53)'V',1,0,0d0,0d0,0d0,0d0,
     &             1,2

              massgam=(pup(4,5)-pup(4,4))**2-(pup(3,5)-pup(3,4))**2
     &             -(pup(2,5)-pup(2,4))**2-(pup(1,5)-pup(1,4))**2
              massgam=-dsqrt(-massgam)
              write(45,54)'P',1,idup(4),pup(1,4)
     &                ,pup(2,4),pup(3,4),pup(4,4),pup(5,4),
     &             istup(4),0d0,0d0,1,0,0
              write(45,54)'P',2,idup(5),pup(1,5)
     &                ,pup(2,5),pup(3,5),pup(4,5),pup(5,5),
     &             istup(5),0d0,0d0,0,0,0
              write(45,54)'P',3,22,pup(1,4)-pup(1,5)
     &             ,pup(2,4)-pup(2,5),pup(3,4)-pup(3,5)
     &             ,pup(4,4)-pup(4,5),massgam,
     &             2,0d0,0d0,2,0,0

              write(45,53)'V',2,0,0d0,0d0,0d0,0d0,
     &             1,2
              massgam=(pup(4,5)-pup(4,4))**2-(pup(3,5)-pup(3,4))**2
     &             -(pup(2,5)-pup(2,4))**2-(pup(1,5)-pup(1,4))**2
              massgam=-dsqrt(-massgam)
              write(45,54)'P',3,22,pup(1,4)-pup(1,5)
     &             ,pup(2,4)-pup(2,5),pup(3,4)-pup(3,5)
     &             ,pup(4,4)-pup(4,5),massgam,
     &             2,0d0,0d0,2,0,0
              write(45,54)'P',4,22,pup(1,3)
     &                ,pup(2,3),pup(3,3),pup(4,3),0d0,
     &             istup(3),0d0,0d0,2,0,0
              write(45,54)'P',5,idup(6),pup(1,6)
     &                ,pup(2,6),pup(3,6),pup(4,6),pup(5,6),
     &             istup(6),0d0,0d0,0,0,0
              write(45,54)'P',6,idup(7),pup(1,7)
     &                ,pup(2,7),pup(3,7),pup(4,7),pup(5,7),
     &             istup(7),0d0,0d0,0,0,0

              
           endif
              
c$$$              do n=1,nvert
c$$$                 
c$$$                 write(45,53)'V',vert(n),0,0d0,0d0,0d0,0d0,
c$$$     &                orph(n),nout(n)
c$$$                 
c$$$                 do m=1,orph(n)+nout(n)
c$$$                    write(45,54)'P',barv(n,m,1),pdgv(n,m),momv(n,m,1)
c$$$     &                  ,momv(n,m,2),momv(n,m,3),momv(n,m,4),massv(n,m),
c$$$     &                   statv(n,m),0d0,0d0,barv(n,m,2),0,0
c$$$                 enddo
c$$$                 
c$$$            enddo
c$$$            
c$$$         endif

            
            

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
               if(idup(k).eq.22)pup(5,k)=0d0
            enddo

           if(beam.eq.'prot')then
              if(pup(5,3).gt.1d0)then
                 idup(3)=92
              endif
              if(pup(5,4).gt.1d0)then
                 idup(4)=92
              endif
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

            if(diff.eq.'dd')then
               
               qsq1=(pup(4,5)-pup(4,3))**2-(pup(3,5)-pup(3,3))**2
     &              -(pup(2,5)-pup(2,3))**2-(pup(1,5)-pup(1,3))**2
               qsq1=-qsq1
               
               qsq2=(pup(4,6)-pup(4,4))**2-(pup(3,6)-pup(3,4))**2
     &              -(pup(2,6)-pup(2,4))**2-(pup(1,6)-pup(1,4))**2
               qsq2=-qsq2
               
               if(qsq1.gt.qsq2)then
                  scalup=dsqrt(qsq1)
                  aqcdup=alphas(scalup**2)
               else
                  scalup=dsqrt(qsq2)
               endif
               
            elseif(diff.eq.'sd')then
               
               qsq=(pup(4,5)-pup(4,4))**2-(pup(3,5)-pup(3,4))**2
     &              -(pup(2,5)-pup(2,4))**2-(pup(1,5)-pup(1,4))**2
               qsq=-qsq
               
               scalup=dsqrt(qsq)

            elseif(diff.eq.'el')then

              scalup=dsqrt((pup(4,5)+pup(4,6))**2-(pup(3,5)+pup(3,6))**2
     &              -(pup(2,5)+pup(2,6))**2-(pup(1,5)+pup(1,6))**2)
               
            endif

            aqcdup=alphas(scalup**2)
               
            write(45,*)'<event>'
            write(45,304)nup-1,idprup,xwgtup,scalup,aqedup,aqcdup
            if(beam.eq.'prot'.or.beam.eq.'el')then
               if(diff.eq.'sd')then
                  if(pup(3,3).lt.0d0)then
                  do m=4,4
                     write(45,303)idup(m),istup(m),mothup(1,m),
     &                    mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
     &                    ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
     &                 ,spinup(m)
                  enddo
                  do m=3,3
                     write(45,303)idup(m),istup(m),mothup(1,m),
     &                    mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
     &                    ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
     &                 ,spinup(m)
                  enddo
                  do m=5,5
                     write(45,303)idup(m),istup(m),mothup(1,m),
     &                    mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
     &                    ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
     &                 ,spinup(m)
                  enddo
                  do m=6,nup+1
                     write(45,303)idup(m),istup(m),mothup(1,m),
     &                    mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
     &                    ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
     &                 ,spinup(m)
                  enddo
               else
                  do m=3,nup+1
                     write(45,303)idup(m),istup(m),mothup(1,m),
     &                    mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
     &                    ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
     &                 ,spinup(m)
                  enddo
               endif
               else
                  do m=3,nup+1
                     write(45,303)idup(m),istup(m),mothup(1,m),
     &                    mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
     &                    ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
     &                 ,spinup(m)
                  enddo
               endif
            elseif(beam.eq.'ion'.or.beam.eq.'ionp')then
c$$$               do m=3,nup+2
c$$$                  write(45,203)idup(m),istup(m),mothup(1,m),
c$$$     &                 mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
c$$$     &                 ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
c$$$     &                 ,spinup(m)
c$$$               enddo
               do m=3,nup+1
                  write(45,203)idup(m),istup(m),mothup(1,m),
     &                 mothup(2,m),icolup(1,m),icolup(2,m),pup(1,m)
     &                 ,pup(2,m),pup(3,m),pup(4,m),pup(5,m),vtimup(m)
     &                 ,spinup(m)
               enddo
               
            endif
            write(45,*)'</event>'

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

           if(beam.eq.'prot')then
              if(phep(5,3).gt.1d0)idhep(3)=92
              if(phep(5,4).gt.1d0)idhep(4)=92
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
