      subroutine F1F2el(q2,f1,f2)
      implicit double precision(a-y)

      include 'mp.f'
      include 'fbeam.f'

      a0e=0.98462d0
      a1e=0.68414d0
      a2e=0.01933d0
      a0m=0.28231d0
      a1m=1.34919d0
      a2m=0.55473d0
      ge=(a0e/(1d0+q2/a1e)**2+(1d0-a0e)/(1d0+q2/a2e)**2)**2
      gm=(a0m/(1d0+q2/a1m)**2+(1d0-a0m)/(1d0+q2/a2m)**2)**2
      gm=gm*7.78d0

      ge0=ge
      gm0=gm

!ccccc interpolate

      if(q2.lt.10d0)then
         gm=7.78d0*gmint(q2)**2
         ge=geint(q2)**2
      else
!c         gm=0d0
!c         ge=0d0
         gm=gm0
         ge=ge0
      endif

!cccccccc

!c      ge=1d0/(1d0+q2/0.71d0)**4
!c      gm=7.78d0*ge

      fe=(4d0*mp**2*ge+q2*gm)/(4d0*mp**2+q2)

      tau=q2/4d0*mp**2

      if(fb1)then

!ccc   F1

!c      fe=(dsqrt(ge)+tau*dsqrt(gm))/(1d0+tau)
!c      fe=fe**2

!ccc   F2
!c      fe=(dsqrt(gm)-dsqrt(ge))/(1d0+tau)
!c      fe=fe**2*tau

!cccccc

      endif

      if(fb2)then

!ccc   F1

!c      fe=(dsqrt(ge)+tau*dsqrt(gm))/(1d0+tau)
!c      fe=fe**2

!ccc   F2
!c      fe=(dsqrt(gm)-dsqrt(ge))/(1d0+tau)
!c      fe=fe**2*tau

!cccccc

      endif

      fm=gm

      f1=fm/2d0
      f2=fe

      return
      end

      function geint(t)
      implicit double precision(a-y)
      double precision gcoh(3,1000)
      common/cohar/gcoh
      integer i

      i=nint((t-0.005d0)/0.01d0)

      if(((t-0.005d0)/0.01d0).gt.dble(i))i=i+1

      if(t.lt.0.005d0)then
         m=(gcoh(2,1)-1d0)/0.005d0
         geint=1d0+m*t
      elseif(i.ge.1000)then
         geint=gcoh(3,1000)
      else
         m=(gcoh(2,i+1)-gcoh(2,i))/0.01d0
         geint=gcoh(2,i)+m*(t-gcoh(1,i))
      endif

      return
      end

      function gmint(t)
      implicit double precision(a-y)
      double precision gcoh(3,1000)
      common/cohar/gcoh
      integer i

      i=nint((t-0.005d0)/0.01d0)

      if(((t-0.005d0)/0.01d0).gt.dble(i))i=i+1


      if(t.lt.0.005d0)then
         m=(gcoh(3,1)-1d0)/0.005d0
         gmint=1d0+m*t
      elseif(i.ge.1000)then
         gmint=gcoh(3,1000)
      else
         m=(gcoh(3,i+1)-gcoh(3,i))/0.01d0
         gmint=gcoh(3,i)+m*(t-gcoh(1,i))
      endif

      return
      end


      subroutine readcoh
      implicit double precision(a-y)
      integer i
      include 'gcoh.f'
      character*500 defpath
#if defined(DATA_PATH)
      data defpath/DATA_PATH/
#else
      data defpath/'share/SuperChic'/
#endif
      integer length
      character*500 valuepath
      length = 500

      length=0
      CALL GETENV('SUPERCHIC_DATA_PATH', valuepath)
      length=len(trim(valuepath))

      ! Check if the environment variable is set
      if (length > 0) then
      write(*,*) 'Reading data from(env. var.)',valuepath(1:length)
      open(40,file=valuepath(1:length)//'/SplinesWithVariableKnots.dat')
      else
      write(*,*) 'Reading data from ',trim(defpath)
      open(40,file=trim(defpath) //'/SplinesWithVariableKnots.dat')
      endif


      do i=1,1000
         read(40,*)gcoh(1,i),gcoh(2,i),dum,dum,dum,gcoh(3,i)
      enddo

      close(40)

      return
      end
