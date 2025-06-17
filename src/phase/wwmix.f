      
            
      
      subroutine wdec_pdg(j,i)
      implicit none
      double precision ranwp,ranwm,ran2,ml_wp,ml_wm
      integer i1,i,l,k,j

      include 'wdecay.f'
      include 'diff.f'
      include 'pdg.f'
      include 'record.f'
      include 'leshouches.f'

      if(diff.eq.'el')i1=5
      if(diff.eq.'sd')i1=6
      if(diff.eq.'dd')i1=7

      do k=3,nup+2
         do l=1,4
            pup(l,k)=evrec(j,k,l)
         enddo
         pup(5,k)=dsqrt(dabs(pup(4,k)**2-pup(3,k)**2
     &              -pup(2,k)**2-pup(1,k)**2))
         if(idup(k).eq.22)pup(5,k)=0d0
      enddo

      ml_wp=pup(5,i1+3)
      ml_wm=pup(5,i1+5)

      if(ml_wp.gt.0.1d0)then
         wlp='mu'
      else  
         wlp='el'
      endif

      if(ml_wm.gt.0.1d0)then
         wlm='mu'
      else  
         wlm='el'
      endif


      if(wlp.eq.'mu')then
         pdgid(i1+2)=14
         pdgid(i1+3)=-13
      elseif(wlp.eq.'el')then
         pdgid(i1+2)=12
         pdgid(i1+3)=-11
      else
         pdgid(i1+2)=16
         pdgid(i1+3)=-15
      endif

      if(wlm.eq.'mu')then
         pdgid(i1+4)=-14
         pdgid(i1+5)=13
      elseif(wlm.eq.'el')then
         pdgid(i1+4)=-12
         pdgid(i1+5)=11
      else
         pdgid(i1+4)=-16
         pdgid(i1+5)=15
      endif

      return
      end

      
      subroutine wwmix
      implicit none
      double precision ranwp,ranwm,ran2
      integer i1

      include 'wdecay.f'
      include 'diff.f'
      include 'pdg.f'


      if(diff.eq.'el')i1=5
      if(diff.eq.'sd')i1=6
      if(diff.eq.'dd')i1=7


      if(wlp_lep)then
         ranwp=ran2()
         if(ranwp.lt.0.5d0)then
            wlp='mu'
         else
            wlp='el'
         endif 
      endif
      if(wlm_lep)then
         ranwm=ran2()
         if(ranwm.lt.0.5d0)then
            wlm='mu'
         else
            wlm='el'
         endif 
      endif



      ! if(wlp.eq.'mu')then
      !    pdgid(i1+2)=14
      !    pdgid(i1+3)=-13
      ! elseif(wlp.eq.'el')then
      !    pdgid(i1+2)=12
      !    pdgid(i1+3)=-11
      ! else
      !    pdgid(i1+2)=16
      !    pdgid(i1+3)=-15
      ! endif

      ! if(wlm.eq.'mu')then
      !    pdgid(i1+4)=-14
      !    pdgid(i1+5)=13
      ! elseif(wlm.eq.'el')then
      !    pdgid(i1+4)=-12
      !    pdgid(i1+5)=11
      ! else
      !    pdgid(i1+4)=-16
      !    pdgid(i1+5)=15
      ! endif

c      print*,i1,pdgid(i1+4),pdgid(i1+5)

c      print*,wlm,wlp

      return
      end