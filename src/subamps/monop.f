ccc   gamgam --> SM higgs subprocess amplitude
      subroutine monop(mx,pp,mm,pm,mp)
      implicit none
      double precision game,normh,msq,betags,mx
      complex*16 pp,mm,pm,mp

      include 'pi.f'
      include 'norm.f'
      include 'ewpars.f'
      include 'monopar.f'
      include 'proc.f'
      
      game=(2d0-mmon/mpol)**1.5d0*mpol**3d0
      game=game*2d0/mmon**2/alpha**2
      game=game/2d0

      msq=game*32d0*pi*mmon
      msq=msq/2d0

      if(mx.gt.mmon)then
         betags=1d0-mmon**2/mx**2
      else
         betags=0d0
      endif
      
      if(proc.eq.70)msq=msq*betags**2
      
      normh=dsqrt(msq)
      
      pp=normh*dsqrt(conv/4d0)
      mm=normh*dsqrt(conv/4d0)
      mp=0d0
      pm=0d0

      return
      end







