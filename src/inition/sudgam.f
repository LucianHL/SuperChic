      subroutine sudgaminit
      implicit none
      double precision qtmin,qtmax,lqtmin,lqtmax,qt,lqt
      integer iqtot,i

      iqtot=100
      qtmin=1d-3
      qtmax=2d0

      lqtmin=dlog(qtmin)
      lqtmax=dlog(qtmax)

      do i=1,iqtot
         
         lqt=lqtmin+(lqtmax-lqtmin)*dble(i)/dble(iqtot)
         qt=dexp(lqt)



      enddo




      return
      end



      function sudgam(lt,mll)
      implicit none
      double precision btmax,hb,sum,wt,sudgam,bt,lt,mll
      double precision sudgam_bt
      integer n,ntot

      include 'pi.f'
      include 'ion.f'
      
      sum=0d0

      btmax=rzg*1000d0
      btmax=1d20
      
      ntot=1000
      hb=btmax/dble(ntot)

c      write(6,*)besj0(0d0)
c      stop
      
      do n=1,ntot

         bt=(dble(n)-0.5d0)*hb

         wt=sudgam_bt(bt,mll)
         wt=wt*bt*besj0(bt*lt)
         wt=wt/(2d0*pi)*hb

c         write(6,*)'bt',bt,sudgam_bt(bt,mll)
c         write(6,*)bt*besj0(bt*lt),bt*besj0(bt*lt)*hb

         sum=sum+wt
         
      enddo

      sudgam=sum

      

      return
      end

      function sudgam_bt(bt,mll)
      implicit none
      double precision pow,alpha,eul,mub,bt,mll
      double precision sudgam_bt,test,sud,const

      include 'mq.f'
      include 'pi.f'

      alpha=1d0/137d0
      eul=0.57721d0
      mub=2d0*dexp(-eul)/bt

      pow=alpha/2d0/pi*dlog(mll**2/mq**2)

c      write(6,*)'mub',mub,2d0*dexp(-eul)

      if(mub.gt.mq)then
         sud=alpha/2d0/pi*dlog(mll**2/mub**2)**2
      else
         sud=alpha/2d0/pi*dlog(mll**2/mq**2)
c         sud=sud*(dlog(mll**2/mub**2)+dlog(mq**2/mub**2))
         sud=sud*dlog(mll**2*mq**2/mub**4)

c         sud=mll**2*mq**2/mub**4
c         sud=sud**(pow)
c         sud=dlog(sud)
c         write(6,*)'sud',mll,mub,mq,mll**2/mub**2,dlog(mll**2/mub**2)
      endif
      
      sudgam_bt=dexp(-sud)

c      sudgam_bt=1d0/sud

c      pow=alpha/2d0/pi*dlog(mll**2/mq**2)
      const=8d0*dexp(-4d0*eul)/mll**2/mq**2
      const=const**pow

      test=const*bt**(-4d0*pow)
c      test=const*mub**(4d0*pow)
c      write(6,*)bt,sudgam_bt,test,4d0*pow

      return
      end
