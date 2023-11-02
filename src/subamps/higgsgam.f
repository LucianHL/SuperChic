ccc   gamgam --> SM higgs subprocess amplitude
      subroutine higgsgam(mx,pp,mm,pm,mp)
      implicit none
      double precision mx
      complex*16 pp,mm,pm,mp

      include 'pi.f'
      include 'ewpars.f'
      include 'normh.f'
      include 'norm.f'

      pp=normh*dsqrt(conv/4d0)
      mm=normh*dsqrt(conv/4d0)
      mp=0d0
      pm=0d0

      return
      end







