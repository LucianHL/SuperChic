ccc   gg --> eta(')eta(') subprocess amplitudes
      subroutine eta(p,mx,u,t,ppqqs,ppqg,ppgg,pmqqs,pmqqns)
      implicit none
      double precision qqn1,mx,pmqqns,pmqqs,ppgg,ppqg,ppqqs
      double precision ppn,mmn,pmn,mpn
      double precision out1,out2,out3,norm,cost,beta
      double precision ggn1,a2gn,a28n,a21n
      double precision alphas,u,t
      integer p

      include 'pi.f'
      include 'mq.f'
      include 'mixing.f'
      include 'partonmom2.f'

      beta=dsqrt(p1(1)**2+p1(2)**2+p1(3)**2)/p1(4)
      cost=p1(3)/p1(4)/beta

      ppn=(1d0+cost**2)/(1d0-cost**2)**2
      mmn=ppn

      pmn=(1d0+3d0*cost**2)/(1d0-cost**2)**2/2d0
      mpn=pmn

cccccccc

      norm=64d0*pi**2*alphas(mx**2/4d0)**2

      call wfoctet(mx,2,a28,a28n)
      call wfsinglet(mx,a21,a2g,a21n,a2gn)

      call mesint(p,1,cost,out1)
      call mesint(p,2,cost,out2)
      call mesint(p,3,cost,out3)

      out2=out2*(a28n/a28)
      out3=out3*(a28n/a28)**2

      qqn1=(1d0+a21n)
      ggn1=5d0/3d0*a2gn

cccccccc

      ppqqs=qqn1**2*ppn
      ppqqs=ppqqs*(6d0/(2d0*dsqrt(6d0)))**2
      ppqqs=ppqqs*norm/3d0/mx**2
      ppqqs=ppqqs*3d0

      pmqqns=(out1+2d0*out2+out3)*norm
      pmqqns=pmqqns*(6d0/(2d0*dsqrt(6d0)))**2
      pmqqns=pmqqns/2d0/mx**2

      pmqqs=qqn1**2*pmn
      pmqqs=pmqqs*(6d0/(2d0*dsqrt(6d0)))**2
      pmqqs=pmqqs*norm/3d0/mx**2
      pmqqs=pmqqs*3d0

ccccccccc

      ppqg=qqn1*ggn1*ppn
      ppqg=ppqg*norm/3d0/mx**2
      ppqg=ppqg*2d0*dsqrt(3d0**3/8d0)
      ppqg=ppqg*6d0/(2d0*dsqrt(6d0))
      ppqg=ppqg/2d0/dsqrt(6d0)*dsqrt(4d0/3d0/2d0/3d0)
      ppqg=ppqg*dsqrt(3d0)

      ppgg=ggn1**2*ppn
      ppgg=ppgg*norm/3d0/mx**2
      ppgg=ppgg*4d0*3d0**3/8d0
      ppgg=ppgg*(1d0/2d0/dsqrt(6d0)*dsqrt(4d0/3d0/2d0/3d0))**2

ccccccccc

      return
      end


