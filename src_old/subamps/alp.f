ccc   gamgam --> ALP amplitude
      subroutine alp(mx,pp,mm,pm,mp)
      implicit none
      double precision norma,mx
      complex*16 pp,mm,pm,mp

      include 'pi.f'
      include 'ewpars.f'
      include 'norm.f'
      include 'gax.f'

      norma=gax/2d0*mx**2
      norma=norma*dsqrt(conv)
      norma=norma/2d0   ! to get correct M^2/4

      pp=norma
      mm=-norma

      if(alpt.eq.'ps')then
         mm=-norma
      elseif(alpt.eq.'sc')then
         mm=norma
      endif
      pm=0d0
      mp=0d0

      return
      end

      subroutine alp_off(p)
      implicit none
      complex*16 zt,zout
      double precision width,uh,th
      double precision qsq1,qsq2,q1_q2,mx
      integer i1,i2,i,j,p
      double precision q1(4),q2(4)

      include 'pi.f'
      include 'ewpars.f'
      include 'gax.f'
      include 'zoutarr.f'
      include 'gmatrices.f'
      include 'mom.f'
      include 'partonmom4.f'
      include 'zi.f'
      include 'norm.f'

      width=gax**2*malp**3/64d0/pi
      mx=dsqrt(q(4,5)**2-q(3,5)**2-q(2,5)**2-q(1,5)**2)
      th=-mx*(paa(4)-paa(3))
      uh=-mx*(paa(4)+paa(3))

      do i=1,4
         q1(i)=q(i,1)-q(i,3)
         q2(i)=q(i,2)-q(i,4)
      enddo

      qsq1=(q(4,3)-q(4,1))**2-(q(3,3)-q(3,1))**2-(q(2,3)-q(2,1))**2
     &     -(q(1,3)-q(1,1))**2
      qsq1=-qsq1
      qsq2=(q(4,4)-q(4,2))**2-(q(3,4)-q(3,2))**2-(q(2,4)-q(2,2))**2
     &     -(q(1,4)-q(1,2))**2
      qsq2=-qsq2

      q1_q2=q1(4)*q2(4)-q1(3)*q2(3)-q1(2)*q2(2)-q1(1)*q2(1)

      do i1=1,4
         do i2=1,4

            zout=0d0

            do i=1,4
               do j=1,4

                  if(alpt.eq.'ps')then
                     zt=e_(i1,i2,i,j)*q1(i)*q2(j)*gax
                  endif

                  if(i.lt.4)zt=-zt
                  if(j.lt.4)zt=-zt

                  zout=zout+zt

               enddo
            enddo

            if(alpt.eq.'sc')then
               zout=(d_(i1,i2)*q1_q2-q1(i2)*q2(i1))*gax
            endif

            zout=zout*dsqrt(conv)
            zoutarr(p,i1,i2)=zout

         enddo
      enddo

      return
      end


ccc   gamgam --> VX  amplitude
      subroutine sim_modVX(mx,pp,mm,pm,mp)
      implicit double precision (a-z)
      complex*16 pp,mm,pm,mp

      include 'mxs.f'

      tau=0.04d0

      norm=dexp(-tau*mx)  ! mx distribution

      pp=dsqrt(norm) ! square root as amplitude here

      mm=pp
      pm=0d0
      mp=0d0

      return
      end
