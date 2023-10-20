      subroutine gdrin
      implicit none
      double precision gdrx_res,regge_gdr
      double precision xmax,lemax,lemin,leint,le
      double precision o0,emax,e,shad
      double precision sigt(2,62),test
      integer i,j,ir

      include 'mion.f'
      include 'ion.f'
      include 'gdr.f'
      include 'sAA.f'
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
        write(*,*) 'Reading data from(env. var.) ',valuepath(1:length)
        open(10,file=valuepath(1:length) // '/Veyssiere_singleneut.dat')
        open(11,file=valuepath(1:length) // '/Lepretre25_103.dat')
        open(12,file=valuepath(1:length) // '/Carlos_106_440.dat')
        open(13,file=valuepath(1:length) // '/gampgamn.dat')
        open(14,file=valuepath(1:length) // '/Caldwell.dat')
        open(50,file=valuepath(1:length) // '/Muccifora.dat')
      else
      write(*,*) 'Reading data from ',trim(defpath)
      open(10,file=trim(defpath)// '/Veyssiere_singleneut.dat')
      open(11,file=trim(defpath)// '/Lepretre25_103.dat')
      open(12,file=trim(defpath)// '/Carlos_106_440.dat')
      open(13,file=trim(defpath)// '/gampgamn.dat')
      open(14,file=trim(defpath)// '/Caldwell.dat')
      open(50,file=trim(defpath)// '/Muccifora.dat')
      endif
      
      i0=161 ! GDR, Veyssiere et al. Nucl. Phys. A159, 561 (1970)
      i1=189 ! 25-103 MeV, Lepretre, et al., Nucl. Phys. A367, 237 (1981)
      i2=210 ! 106-440 MeV, Carlos, et al., Nucl. Phys. A431, 573 (1984)
      i3=272 ! 440 MeV - 2 GeV, Armstrong et. al Phys. Rev. D 5, 1640 1972, Nucl. Phys. B41, 445 (1972)
      i3=i2+13
      i4=283 ! 2 - 16.4 GeV, Caldwell, Phys. Rev. D 7, 1362 (1973)
      i4=i3+11
      i5=383 ! Regge fit from 16.4 GeV to E_max
      i5=i4+100

      o0=7.4d0
      shad=0.65d0
      shad=1d0

      emax=(saa-2d0*mion**2)/2d0/mion ! in theory, but better to just have Q^2 < 10 GeV^2
      xmax=10d0/2d0/mion**2*(-1d0+dsqrt(1d0+4d0*mion**2/10d0))
      emax=emax*xmax

      lemax=dlog(emax)
      lemin=dlog(16.4d0)
      ir=100
      leint=(lemax-lemin)/dble(ir)
      
      if(nint(az).eq.82)then
         read(10,*)
         read(10,*)(sneut(2,j),j=1,i0)
      elseif(nint(az).eq.79)then
         read(10,*)(sneut(2,j),j=1,i0)
      else
         read(10,*)
         read(10,*)(sneut(2,j),j=1,i0)
         do i=1,i0
            sneut(2,i)=sneut(2,i)*(an-az)/an**(2d0/3d0)/3.397d0 ! reweight so follow TRK sum rule
         enddo
      endif

      do i=1,i0
         sneut(1,i)=o0+0.1d0*(i-1)
         sneut(2,i)=sneut(2,i)*1d3  ! convert to mb
         mneut(1,i)=o0+0.1d0*(i-1)
         mneut(2,i)=gdrx_res(mneut(1,i))
      enddo


      read(11,*)(mneut(1,j),j=i0+1,i1)
      read(11,*)(mneut(2,j),j=i0+1,i1)

      do i=i0+1,i1
         mneut(2,i)=mneut(2,i)*an/208d0
      enddo

      read(12,*)(mneut(1,j),j=i1+1,i2)
      read(12,*)(mneut(2,j),j=i1+1,i2)

      do i=i1+1,i2
         mneut(2,i)=mneut(2,i)*an/208d0
      enddo

      read(13,*)(sigt(1,j),j=1,62)
      read(13,*)(sigt(2,j),j=1,62)
      
      do i=1,13
         read(50,*)mneut(1,i+i2),mneut(2,i+i2)
         mneut(1,i+i2)=mneut(1,i+i2)*1d3
         mneut(2,i+i2)=mneut(2,i+i2)*an/1d3
      enddo



      read(14,*)(mneut(1,j),j=i3+1,i4)
      read(14,*)(mneut(2,j),j=i3+1,i4)

      do i=1,11
         mneut(2,i3+i)=mneut(2,i3+i)*an      
         e=mneut(1,i3+i)*1d-3 
      enddo

      do i=i4+1,i5
         le=lemin+dble(i-i4)*leint
         e=dexp(le)
         
         mneut(1,i)=e*1d3
         mneut(2,i)=Regge_gdr(e)*an

      enddo
         
      return
      end


      function gdrx_res(e)
      implicit none
      double precision e,e1,gam1,sig1,gdrx_res

      include 'ion.f'

      if(nint(az).eq.79)then
         sig1=540d0
         e1=13.7d0
         gam1=4.75d0
      else
         sig1=640d0
         e1=13.42d0
         gam1=4.05d0
      endif

      gdrx_res=sig1*e**2*gam1**2/((e**2-e1**2)**2+e**2*gam1**2)
      
      if(nint(az).eq.79)then
      elseif(nint(az).eq.82)then
      else
         gdrx_res=gdrx_res*(an-az)/an**(2d0/3d0)/3.397d0 ! reweight so follow TRK sum rule
      endif

      return
      end

      function Regge_gdr(e)
      implicit none
      double precision regge_gdr,shad
      double precision y,x,vec,s,pom,mn,eta,eps,e

      include 'mion.f'
      include 'ion.f'
      
      mn=0.94d0
      
      x=0.0677d0
      y=0.129d0
      eps=0.0808d0
      eta=0.4525d0

      x=0.057d0
      y=0.121d0
      eps=0.1d0
      eta=0.716d0/2d0

      shad=0.65d0

      s=2d0*mn*e+mn**2

      pom=x*s**eps
      vec=y*s**(-eta)

      Regge_gdr=shad*(pom+vec)

      if(e.gt.500d0)regge_gdr=0d0

      return
      end
