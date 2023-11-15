      subroutine write_particle(i)
      implicit none
      integer i
      integer vert(20,4)
      integer prodv(20)
      integer endv(20)
      integer nvert
      integer ii
      COMMON /HM/nvert,vert,prodv,endv
      include 'leshouches.f'
      include 'hepevt.f'      
      ii = isthep(i)
      write(45,154)'P',i,idup(i),pup(1,i),pup(2,i),pup(3,i),pup(4,i)
     & ,0d0,ii,0d0,0d0,-endv(i),0

 154      format(A1,' ',i0,' ',i0,' ',E15.9,' ',E15.9,' ',E15.9
     &        ,' ',E15.9,' ',E15.9,' ',i0,' ',E15.9,' ',E15.9,' ',i0
     &        ,' ',i0)

      end


      subroutine write_vertex(i)
      implicit none
      integer i
      integer vert(20,4)
      integer prodv(20)
      integer endv(20)
      integer nvert
      COMMON /HM/nvert,vert,prodv,endv
      include 'leshouches.f'
      include 'hepevt.f'        
      write(45,153)'V',-i,0,0d0,0d0,0d0,0d0,vert(i,3),vert(i,4),0
 153      format(A1,' ',i0,' ',i0,' ',E15.9,' ',E15.9,' ',E15.9,' '
     &        ,E15.9,' ',i0,' ',i0,' ',i0)
      end
     
      subroutine hmout(last,enr,nfl1,nfl2,x1,x2,proc)
      implicit none
      integer enr,nfl1,nfl2,proc
      !double precision scalup,aqcdup,aqedup, x1,x2, xsecup,xerrup
      double precision x1,x2 !,xsecup,xerrup, scalup
      integer vert(20,4)
      integer mom(20,2)
      integer prodv(20)
      integer endv(20)
      integer nvert
      integer i,j
      COMMON /HM/nvert,vert,prodv,endv
      integer st,last     
      LOGICAL found
      include 'leshouches.f'      
      include 'hepevt.f'      
      nvert=0
      st=1
      do i = 1,20      
       vert(i,1)=0
       vert(i,2)=0
       vert(i,3)=0
       vert(i,4)=0
       prodv(i)=0
       endv(i)=0
       mom(i,1)=0
       mom(i,2)=0
      end do    
      do i = 1,nhep
        if (jdahep(2,i).eq.0) jdahep(2,i)=jdahep(1,i)
      enddo
      do i = 1,nhep
      do j = jdahep(1,i),jdahep(2,i)
        if (j.ne.0) then
        if (jmohep(1,j).eq.0) jmohep(1,j) = i
        if (jmohep(2,j).eq.0) jmohep(2,j) = i
         jmohep(1,j) = min(i,jmohep(1,j))
         jmohep(2,j) = max(i,jmohep(2,j))
        endif
      enddo
      enddo
      
      do i = 1,nhep
        mom(i,1)=jmohep(1,i)
        mom(i,2)=jmohep(2,i)
      end do
  
      do i = 1,nhep
      found = .FALSE.
      do j=1,nvert
      if (prodv(i) .ne. 0 ) then
        found = .TRUE.
        continue
      endif
      if (mom(i,1) .eq.vert(j,1).and.mom(i,2) .eq.vert(j,2))then
      prodv(i) = j
      found = .TRUE.
      endif
      end do
      if ( .not. found )  then
      nvert = nvert +1
      prodv(i) = nvert
      vert(nvert,1) = mom(i,1)
      vert(nvert,2) = mom(i,2)
      endif
      end do
 
      do i = 1,nhep
      do j = 1,nhep
      if  (i .ge. mom(j,1) .and. i .le. mom(j,2) ) then
      endv(i) = prodv(j)
      endif
      end do
      end do
     
      do i = 1,nhep
      if (endv(i).ne.0) vert(endv(i),3)=vert(endv(i),3)+1
      if (prodv(i).ne.0) vert(prodv(i),4)=vert(prodv(i),4)+1
      end do

      write(45,151)'E',enr,0,scalup,aqcdup,aqedup,proc,
     &0,nvert,1,2,0,1,1d0
      write(45,155)'U GEV CM'
      write(45,'(A)')'N 1 "Default"'
      write(45,156)'C',xsecup(1),xerrup(1)
      write(45,152)'F',nfl1,nfl2,x1,x2,scalup,0d0,0d0,0,0 
      do i=1,nvert
      call write_vertex(i)
      do j=1,nhep
      if ((prodv(j).eq. i) .or. (prodv(j) .eq. 0 .and. i .eq. 1))  then
      call write_particle(j)
      endif
      enddo
      end do

 151      format(A1,' ',I0,' ',I0,' ',E15.9,' ',E15.9,' ',E15.9,' ',I0
     &        ,' ',I0,' ',I0,' ',I0,' ',I0,' ',I0,' ',I0,' ',E15.9)
 152      format(A1,' ',i0,' ',i0,' ',E15.9,' ',E15.9,' ',E15.9,' '
     &        ,E15.9,' ',E15.9,' ',i0,' ',i0)
 156      format(A1,' ',E15.9,' ',E15.9)
 155      format(8a)     
      end

 




ccc   randomizes order of VEGAS unweighted events and
ccc   prints nev events to record
      subroutine unwprintq
      implicit double precision(a-y)
      integer i,j,k,l,m,nfl1,nfl2
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

            if(i.eq.1)then
               write(45,'(A)')'HepMC::Version 2.06.11'
               write(45,'(A)')'HepMC::IO_GenEvent-START_EVENT_LISTING'
            endif
            
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

       call hmout(nup+1,i,nfl1,nfl2,x1,x2,proc) 
       if(i.eq.nev)then
         write(45,'(A)')'HepMC::IO_GenEvent-END_EVENT_LISTING'
       endif

        goto 500
       endif

ccccccccccccccccccccccccccccccccccccccccccccccc
ccccc Les Houches
ccccccccccccccccccccccccccccccccccccccccccccccc

         if(erec.eq.'lhe')then
            
            do k=1,nup+1
               idup(k)=pdgid(k)
            enddo
            
           do k=3,nup+1
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
               write(45,'(A)')'<LesHouchesEvents version="1.0">'
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
               
            write(45,201)'E ',i,nhep
            isthep(1)=4
            isthep(2)=4

            if(beam.eq.'el'.or.beam.eq.'prot')then
            do m=1,nhep
               write(45,300)isthep(m),idhep(m),jmohep(1,m),
     &              jmohep(2,m),jdahep(1,m),jdahep(2,m),
     &              phep(1,m),phep(2,m),phep(3,m),phep(4,m)
     &              ,phep(5,m),vhep(1,m),vhep(2,m),vhep(3,m),vhep(4,m)
            enddo
            elseif(beam.eq.'ion')then
            do m=1,nhep
               write(45,200)isthep(m),idhep(m),jmohep(1,m),
     &              jmohep(2,m),jdahep(1,m),jdahep(2,m),
     &              phep(1,m),phep(2,m),phep(3,m),phep(4,m)
     &              ,phep(5,m),vhep(1,m),vhep(2,m),vhep(3,m),vhep(4,m)
            enddo
            elseif(beam.eq.'ionp')then
            do m=1,nhep
            write(45,300)isthep(m),idhep(m),jmohep(1,m),
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

 201  format(A1,i8,1x,i8)

 200  format(i8,i8,i8,i8,i8,i8,E19.8,E19.8
     &,E19.8,E19.8,E19.8,/,48x,E19.8,E19.8,
     & E19.8,E19.8)

 300  format(i8,i8,i8,i8,i8,i8,E19.8,E19.8
     &,E19.8,E19.8,E19.8,/,48x,E19.8,E19.8,
     & E19.8,E19.8)

 303  format(7x,i10,1x,i8,1x,i4,1x,i4,1x,i4,1x,i4,1x,E16.9,1x,
     &     E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x,F2.0,1x,F2.0)

 203  format(7x,i10,1x,i8,1x,i4,1x,i4,1x,i4,1x,i4,1x,E16.9,1x,
     &E16.9,1x,E16.9,1x,E16.9,1x,E16.9,1x,F2.0,1x,F2.0)

      return
      end
