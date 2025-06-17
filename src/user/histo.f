!cc   initialises number of histograms
      subroutine inithist(nhist)
      implicit none
      integer nhist,i

      do i=1,nhist
         call histo3(i)
      enddo

      return
      end

!cc   binning subroutine
      subroutine binit(wt)
      implicit none
      double precision wt
      integer i

      include 'vars.f'
      include 'mom.f'
      include 'range.f'
      include 'proc.f'
      include 'decay.f'
      include 'pi.f'
      include 'partonmom3.f'
      include 'x.f'
      include 'diss.f'
      include 'mp.f'
      include 'xb.f'
      include 'ewpars.f'
      include 'partonmom2.f'

!cccccccc

      mmin0=6d0

      if(dps.eq.1)then
         call histo1(1,10,ymin,ymax,yx,wt,'yx')
      else
         call histo1(1,10,mmin0,mmax,mx,wt,'mx')
         call histo1(2,10,ymin,ymax,yx,wt,'yx')
      endif

      return
      end


!cc   prints histograms
      subroutine histo1(ih,ib,x0,x1,x,w,labin)
      implicit real*8(a-h,o-z)
      integer ix,iz,iv,ic,ih,ib,il,i,ib1
      character *(*) labin
      character*256 fname
      integer outl
      character*1 regel(30),blank,star
      integer io,iu,ii
      dimension h(20,100),hx(20),io(20),iu(20),ii(20)
      dimension y0(20),y1(20),ic(20)
      data regel / 30*' ' /,blank /' ' /,star /'*'/
      save

      include 'output.f'
      include 'iteration.f'
      include 'lab.f'
      include 'nhist.f'

      if(nhist.lt.ih)nhist=ih

      lab(ih)=labin

      call length(labin,outl)

      y0(ih)=x0
      y1(ih)=x1
      ic(ih)=ib
      if(x.lt.x0) goto 11
#if defined(__NVCOMPILER)
      if(x.gt.x1 .or. ( x .ne. x) ) goto 12
#else
      if(x.gt.x1 .or. isnan(x) ) goto 12
#endif

      ix=idint((x-x0)/(x1-x0)*dble(ib))+1


      h(ih,ix)=h(ih,ix)+w
      if(h(ih,ix).gt.hx(ih)) hx(ih)=h(ih,ix)
      ii(ih)=ii(ih)+1


      return
   11 iu(ih)=iu(ih)+1
      return
   12 io(ih)=io(ih)+1
      return
      entry histo2(ih,il)

      call length(outtag,outl)
      fname='outputs/output'//outtag(1:outl)//'.dat'
      open(10,file=trim(fname),POSITION='append',Status='old')

       call length(lab(ih),outl)


      ib1=ic(ih)
      x01=y0(ih)
      x11=y1(ih)
      bsize=(x11-x01)/dble(ib1)
      hx(ih)=hx(ih)*(1.d0+1.d-06)
      write(6,*)''
      write(6,*)'***************************************'
      write(6,*)lab(ih)(1:outl)
      write(6,*)'***************************************'
      write(6,21)ii(ih),iu(ih),io(ih)
      write(6,*)'***************************************'

   21 format(' inside,under,over : ',3i6)

      if(ii(ih).eq.0) goto 28
      write(6,23)
   23 format(35(1h ),3(10h----+----i))
      do 27 iv=1,ib1
         z=(dble(iv)-0.5d0)/dble(ib1)*(x11-x01)+x01
         zl=(dble(iv)-1d0)/dble(ib1)*(x11-x01)+x01
         zh=(dble(iv))/dble(ib1)*(x11-x01)+x01
      if(il.eq.1) goto 24
      iz=idint(h(ih,iv)/hx(ih)*30.)+1
      goto 25
   24 iz=-1
      if(h(ih,iv).gt.0.d0) iz=idint(dlog(h(ih,iv))/dlog(hx(ih))*30.)+1
   25 if(iz.gt.0.and.iz.le.30) regel(iz)=star
      h(ih,iv)=h(ih,iv)
      write(6,26) z,h(ih,iv)/bsize/dble(it),(regel(i),i=1,30)
      write(10,*)z,zl,zh,h(ih,iv)/bsize/dble(it)

   26 format(1h ,2g15.6,4h   i,30a1,1hi)
      if(iz.gt.0.and.iz.le.30) regel(iz)=blank
   27 continue
      write(6,23)

      write(10,*)
      close(10)

      return
   28 write(6,29)
   29 format('  no entries inside histogram')

      return

      entry histo3(ih)
      do 31 i=1,100
   31 h(ih,i)=0.
      hx(ih)=0.
      io(ih)=0
      iu(ih)=0
      ii(ih)=0
      return
      end

