    
      function ampWWt_fun(pp,pm,ep,em,q1,q2,ep_em,ep_q1,ep_q2,em_q1,
     &     em_q2,pp_q1,pp_q2,pm_q1,pm_q2,q1_q2,pp_pm,qsq1,qsq2,mw,d_,
     &     mu,nu,tl,ul)
      implicit none
      complex*16 ep(4),em(4),ep_em,ep_q1,ep_q2,em_q1,em_q2,
     &     ampWWt_fun,ampWWt,d_(4,4)
      double precision pp(4),pm(4),q1(4),q2(4),q1_q2,pp_pm,qsq1,
     &     qsq2,mw,pp_q1,pp_q2,pm_q1,pm_q2,tl,ul,tl2
      integer mu,nu

      ampWWt =
     & 4*pp(mu)*pp(nu)*ep_em - 4*pp(mu)*ep(nu)*em_q2 + 4*pp(mu)*em(nu)*
     & ep_q2 - 4*pp(mu)*q1(nu)*ep_em - 2*pp(mu)*q2(nu)*ep_em - 4*pp(nu)
     & *ep(mu)*em_q1 + 4*pp(nu)*em(mu)*ep_q1 - 2*pp(nu)*q1(mu)*ep_em - 
     & ep(mu)*em(nu)*mw**(-2)*qsq2*qsq1 + 2*ep(mu)*em(nu)*ul + ep(mu)*
     & em(nu)*tl - 3*ep(mu)*em(nu)*mw**2 + 4*ep(mu)*q1(nu)*em_q1 + 4*
     & ep(mu)*q1(nu)*em_q2 + 2*ep(mu)*q2(nu)*em_q1 - ep(mu)*q2(nu)*
     & em_q2*mw**(-2)*qsq1 + 2*ep(nu)*q1(mu)*em_q2 - 4*em(mu)*q1(nu)*
     & ep_q1 - 2*em(mu)*q2(nu)*ep_q1 - em(nu)*q1(mu)*ep_q1*mw**(-2)*
     & qsq2 - 2*em(nu)*q1(mu)*ep_q2 + 4*em(nu)*q2(mu)*ep_q1 + 2*q1(mu)*
     & q1(nu)*ep_em + q1(mu)*q2(nu)*ep_em - q1(mu)*q2(nu)*ep_q1*em_q2*
     & mw**(-2) - 4*d_(mu,nu)*ep_q1*em_q2

      
      
      ampWWt_fun = ampWWt
      return
      end

      function ampWWu_fun(pp,pm,ep,em,q1,q2,ep_em,ep_q1,ep_q2,em_q1,
     &     em_q2,pp_q1,pp_q2,pm_q1,pm_q2,q1_q2,pp_pm,qsq1,qsq2,mw,d_,
     &     mu,nu,tl,ul)
      implicit none
      complex*16 ep(4),em(4),ep_em,ep_q1,ep_q2,em_q1,em_q2,
     &     ampWWu_fun,ampWWu,d_(4,4)
      double precision pp(4),pm(4),q1(4),q2(4),q1_q2,pp_pm,qsq1,
     &     qsq2,mw,pp_q1,pp_q2,pm_q1,pm_q2,tl,ul
      integer mu,nu

      
      ampWWu =
     & 4*pp(mu)*pp(nu)*ep_em - 4*pp(mu)*ep(nu)*em_q2 + 4*pp(mu)*em(nu)*
     & ep_q2 - 2*pp(mu)*q2(nu)*ep_em - 4*pp(nu)*ep(mu)*em_q1 + 4*pp(nu)
     & *em(mu)*ep_q1 - 2*pp(nu)*q1(mu)*ep_em - 4*pp(nu)*q2(mu)*ep_em + 
     & 2*ep(mu)*q2(nu)*em_q1 - ep(nu)*em(mu)*mw**(-2)*qsq2*qsq1 + 
     & ep(nu)*em(mu)*ul + 2*ep(nu)*em(mu)*tl - 3*ep(nu)*em(mu)*mw**2 - 
     & ep(nu)*q1(mu)*em_q1*mw**(-2)*qsq2 + 2*ep(nu)*q1(mu)*em_q2 + 4*
     & ep(nu)*q2(mu)*em_q1 + 4*ep(nu)*q2(mu)*em_q2 + 4*em(mu)*q1(nu)*
     & ep_q2 - 2*em(mu)*q2(nu)*ep_q1 - em(mu)*q2(nu)*ep_q2*mw**(-2)*
     & qsq1 - 2*em(nu)*q1(mu)*ep_q2 - 4*em(nu)*q2(mu)*ep_q2 + q1(mu)*
     & q2(nu)*ep_em - q1(mu)*q2(nu)*ep_q2*em_q1*mw**(-2) + 2*q2(mu)*
     & q2(nu)*ep_em - 4*d_(mu,nu)*ep_q2*em_q1
      
      ampWWu_fun = ampWWu
      return
      end

      
      function ampWWs_fun(ep,em,ep_em,d_,mu,nu)
      implicit none
      complex*16 ep(4),em(4),ep_em,ampWWs_fun,ampWWs,d_(4,4)
      integer mu,nu
      ampWWs =
     &     ep(mu)*em(nu) + ep(nu)*em(mu) - 2*d_(mu,nu)*ep_em
      ampWWs_fun = ampWWs
      return
      end

      
      function ampax_fun(pp,pm,ep,em,q1,q2,ep_em,ep_q1,ep_q2,
     &              em_q1,em_q2,pp_q1,pp_q2,pm_q1,pm_q2,tl,ul,qsq1,qsq2,
     &     mw,d_,mu,nu,pp_em,pp_ep,n,n_q1,n_q2,n_pp,n_pm,n_n)
      implicit none
      complex*16 ep_em,ep_q1,ep_q2,em_q1,em_q2,pp_em,pp_ep,ep(4),em(4)
     &     ,ampax_fun,amp1ax,amp2ax,amp3ax,amp4ax,d_(4,4),
     &     amp5ax,amp3axn,amp4axn,amp3axnn,amp4axnn,pm_em,zt,ztest,
     &     zt1,zt2,zt3,zt4,zt5,zt6,ztest1,ztest2,ztest3,ztest4,ztest5,
     &     ztest6,zt7,zt8,ztest7,ztest8,amp3axnind,amp4axnind
      double precision tl,ul,qsq1,pp(4),pm(4),q1(4),q2(4),
     &     qsq2,mw,pp_q1,pp_q2,pm_q1,pm_q2,n_q1,n_q2
     &     ,den1ax,den2ax,den3ax,den4ax,den5ax,n(4),
     &     den3axn,den4axn,den3axnn,den4axnn,n_pp,n_pm,n_n
      integer mu,nu,nu2,mu2
      logical test,n_ind

      n_ind = .true.
      
      amp1ax =
     &  + n_q1*n_q2 * (  - 4*pp(mu)*pp(nu)*pp_ep*pp_em + 4*pp(mu)*
     &    pp(nu)*pp_ep*em_q1 + 4*pp(mu)*pp(nu)*pp_ep*em_q2 + 2*pp(mu)*
     &    em(nu)*pp_ep*mw**2 + 4*pp(mu)*q1(nu)*pp_ep*pp_em - 4*pp(mu)*
     &    q1(nu)*pp_ep*em_q1 - 4*pp(mu)*q1(nu)*pp_ep*em_q2 + 2*pp(mu)*
     &    q2(nu)*pp_ep*pp_em - 2*pp(mu)*q2(nu)*pp_ep*em_q1 - 2*pp(mu)*
     &    q2(nu)*pp_ep*em_q2 + 2*pp(nu)*ep(mu)*pp_em*mw**2 - 2*pp(nu)*
     &    ep(mu)*em_q1*mw**2 - 2*pp(nu)*ep(mu)*em_q2*mw**2 + 2*pp(nu)*
     &    q1(mu)*pp_ep*pp_em - 2*pp(nu)*q1(mu)*pp_ep*em_q1 - 2*pp(nu)*
     &    q1(mu)*pp_ep*em_q2 - ep(mu)*em(nu)*mw**4 - 2*ep(mu)*q1(nu)*
     &    pp_em*mw**2 + 2*ep(mu)*q1(nu)*em_q1*mw**2 + 2*ep(mu)*q1(nu)*
     &    em_q2*mw**2 - ep(mu)*q2(nu)*pp_em*mw**2 + ep(mu)*q2(nu)*em_q1
     &    *mw**2 + ep(mu)*q2(nu)*em_q2*mw**2 - em(nu)*q1(mu)*pp_ep*
     &    mw**2 - 2*q1(mu)*q1(nu)*pp_ep*pp_em + 2*q1(mu)*q1(nu)*pp_ep*
     &    em_q1 + 2*q1(mu)*q1(nu)*pp_ep*em_q2 - q1(mu)*q2(nu)*pp_ep*
     &    pp_em + q1(mu)*q2(nu)*pp_ep*em_q1 + q1(mu)*q2(nu)*pp_ep*em_q2
     &     )
      amp1ax = amp1ax + n(mu)*n(nu) * ( pp_ep*pp_em*tl**2 - 2*pp_ep*
     &    pp_em*mw**2*tl + pp_ep*pp_em*mw**4 - pp_ep*em_q1*tl**2 + 2*
     &    pp_ep*em_q1*mw**2*tl - pp_ep*em_q1*mw**4 - pp_ep*em_q2*tl**2
     &     + pp_ep*em_q2*mw**2*tl + pp_em*ep_q1*mw**2*tl - pp_em*ep_q1*
     &    mw**4 - ep_q1*em_q1*mw**2*tl + ep_q1*em_q1*mw**4 - ep_q1*
     &    em_q2*mw**2*tl )
      amp1ax = amp1ax + n(mu)*n_q2 * (  - 2*pp(nu)*pp_ep*pp_em*tl + 2*
     &    pp(nu)*pp_ep*pp_em*mw**2 + 2*pp(nu)*pp_ep*em_q1*tl - 2*pp(nu)
     &    *pp_ep*em_q1*mw**2 + 2*pp(nu)*pp_ep*em_q2*tl - 2*pp(nu)*pp_ep
     &    *em_q2*mw**2 - 2*pp(nu)*pp_em*ep_q1*mw**2 + 2*pp(nu)*ep_q1*
     &    em_q1*mw**2 + 2*pp(nu)*ep_q1*em_q2*mw**2 + em(nu)*pp_ep*mw**2
     &    *tl - em(nu)*pp_ep*mw**4 + em(nu)*ep_q1*mw**4 + 2*q1(nu)*
     &    pp_ep*pp_em*tl - 2*q1(nu)*pp_ep*pp_em*mw**2 - 2*q1(nu)*pp_ep*
     &    em_q1*tl + 2*q1(nu)*pp_ep*em_q1*mw**2 - 2*q1(nu)*pp_ep*em_q2*
     &    tl + 2*q1(nu)*pp_ep*em_q2*mw**2 + 2*q1(nu)*pp_em*ep_q1*mw**2
     &     - 2*q1(nu)*ep_q1*em_q1*mw**2 - 2*q1(nu)*ep_q1*em_q2*mw**2 + 
     &    q2(nu)*pp_ep*pp_em*tl - q2(nu)*pp_ep*pp_em*mw**2 - q2(nu)*
     &    pp_ep*em_q1*tl + q2(nu)*pp_ep*em_q1*mw**2 - q2(nu)*pp_ep*
     &    em_q2*tl + q2(nu)*pp_ep*em_q2*mw**2 + q2(nu)*pp_em*ep_q1*
     &    mw**2 - q2(nu)*ep_q1*em_q1*mw**2 - q2(nu)*ep_q1*em_q2*mw**2 )
      amp1ax = amp1ax + n(nu)*n_q1 * ( 2*pp(mu)*pp_ep*pp_em*tl - 2*
     &    pp(mu)*pp_ep*pp_em*mw**2 - 2*pp(mu)*pp_ep*em_q1*tl + 2*pp(mu)
     &    *pp_ep*em_q1*mw**2 - 2*pp(mu)*pp_ep*em_q2*tl - ep(mu)*pp_em*
     &    mw**2*tl + ep(mu)*pp_em*mw**4 + ep(mu)*em_q1*mw**2*tl - 
     &    ep(mu)*em_q1*mw**4 + ep(mu)*em_q2*mw**2*tl - q1(mu)*pp_ep*
     &    pp_em*tl + q1(mu)*pp_ep*pp_em*mw**2 + q1(mu)*pp_ep*em_q1*tl
     &     - q1(mu)*pp_ep*em_q1*mw**2 + q1(mu)*pp_ep*em_q2*tl )

      amp2ax =
     &  + n_q1*n_q2 * (  - 4*pp(mu)*pp(nu)*pp_ep*pp_em + 4*pp(mu)*
     &    pp(nu)*pp_ep*em_q1 + 4*pp(mu)*pp(nu)*pp_ep*em_q2 + 2*pp(mu)*
     &    ep(nu)*pp_em*mw**2 - 2*pp(mu)*ep(nu)*em_q1*mw**2 - 2*pp(mu)*
     &    ep(nu)*em_q2*mw**2 + 2*pp(mu)*q2(nu)*pp_ep*pp_em - 2*pp(mu)*
     &    q2(nu)*pp_ep*em_q1 - 2*pp(mu)*q2(nu)*pp_ep*em_q2 + 2*pp(nu)*
     &    em(mu)*pp_ep*mw**2 + 2*pp(nu)*q1(mu)*pp_ep*pp_em - 2*pp(nu)*
     &    q1(mu)*pp_ep*em_q1 - 2*pp(nu)*q1(mu)*pp_ep*em_q2 + 4*pp(nu)*
     &    q2(mu)*pp_ep*pp_em - 4*pp(nu)*q2(mu)*pp_ep*em_q1 - 4*pp(nu)*
     &    q2(mu)*pp_ep*em_q2 - ep(nu)*em(mu)*mw**4 - ep(nu)*q1(mu)*
     &    pp_em*mw**2 + ep(nu)*q1(mu)*em_q1*mw**2 + ep(nu)*q1(mu)*em_q2
     &    *mw**2 - 2*ep(nu)*q2(mu)*pp_em*mw**2 + 2*ep(nu)*q2(mu)*em_q1*
     &    mw**2 + 2*ep(nu)*q2(mu)*em_q2*mw**2 - em(mu)*q2(nu)*pp_ep*
     &    mw**2 - q1(mu)*q2(nu)*pp_ep*pp_em + q1(mu)*q2(nu)*pp_ep*em_q1
     &     + q1(mu)*q2(nu)*pp_ep*em_q2 - 2*q2(mu)*q2(nu)*pp_ep*pp_em + 
     &    2*q2(mu)*q2(nu)*pp_ep*em_q1 + 2*q2(mu)*q2(nu)*pp_ep*em_q2 )
      amp2ax = amp2ax + n(mu)*n(nu) * ( pp_ep*pp_em*ul**2 - 2*pp_ep*
     &    pp_em*mw**2*ul + pp_ep*pp_em*mw**4 - pp_ep*em_q1*ul**2 + 
     &    pp_ep*em_q1*mw**2*ul - pp_ep*em_q2*ul**2 + 2*pp_ep*em_q2*
     &    mw**2*ul - pp_ep*em_q2*mw**4 + pp_em*ep_q2*mw**2*ul - pp_em*
     &    ep_q2*mw**4 - ep_q2*em_q1*mw**2*ul - ep_q2*em_q2*mw**2*ul + 
     &    ep_q2*em_q2*mw**4 )
      amp2ax = amp2ax + n(mu)*n_q2 * ( 2*pp(nu)*pp_ep*pp_em*ul - 2*
     &    pp(nu)*pp_ep*pp_em*mw**2 - 2*pp(nu)*pp_ep*em_q1*ul - 2*pp(nu)
     &    *pp_ep*em_q2*ul + 2*pp(nu)*pp_ep*em_q2*mw**2 - ep(nu)*pp_em*
     &    mw**2*ul + ep(nu)*pp_em*mw**4 + ep(nu)*em_q1*mw**2*ul + 
     &    ep(nu)*em_q2*mw**2*ul - ep(nu)*em_q2*mw**4 - q2(nu)*pp_ep*
     &    pp_em*ul + q2(nu)*pp_ep*pp_em*mw**2 + q2(nu)*pp_ep*em_q1*ul
     &     + q2(nu)*pp_ep*em_q2*ul - q2(nu)*pp_ep*em_q2*mw**2 )
      amp2ax = amp2ax + n(nu)*n_q1 * (  - 2*pp(mu)*pp_ep*pp_em*ul + 2*
     &    pp(mu)*pp_ep*pp_em*mw**2 + 2*pp(mu)*pp_ep*em_q1*ul - 2*pp(mu)
     &    *pp_ep*em_q1*mw**2 + 2*pp(mu)*pp_ep*em_q2*ul - 2*pp(mu)*pp_ep
     &    *em_q2*mw**2 - 2*pp(mu)*pp_em*ep_q2*mw**2 + 2*pp(mu)*ep_q2*
     &    em_q1*mw**2 + 2*pp(mu)*ep_q2*em_q2*mw**2 + em(mu)*pp_ep*mw**2
     &    *ul - em(mu)*pp_ep*mw**4 + em(mu)*ep_q2*mw**4 + q1(mu)*pp_ep*
     &    pp_em*ul - q1(mu)*pp_ep*pp_em*mw**2 - q1(mu)*pp_ep*em_q1*ul
     &     + q1(mu)*pp_ep*em_q1*mw**2 - q1(mu)*pp_ep*em_q2*ul + q1(mu)*
     &    pp_ep*em_q2*mw**2 + q1(mu)*pp_em*ep_q2*mw**2 - q1(mu)*ep_q2*
     &    em_q1*mw**2 - q1(mu)*ep_q2*em_q2*mw**2 + 2*q2(mu)*pp_ep*pp_em
     &    *ul - 2*q2(mu)*pp_ep*pp_em*mw**2 - 2*q2(mu)*pp_ep*em_q1*ul + 
     &    2*q2(mu)*pp_ep*em_q1*mw**2 - 2*q2(mu)*pp_ep*em_q2*ul + 2*
     &    q2(mu)*pp_ep*em_q2*mw**2 + 2*q2(mu)*pp_em*ep_q2*mw**2 - 2*
     &    q2(mu)*ep_q2*em_q1*mw**2 - 2*q2(mu)*ep_q2*em_q2*mw**2 )

      amp3ax =
     &  + n_q1*n_q2*qsq1 * ( 2*pp(nu)*ep(mu)*pp_em*tl - 2*pp(nu)*ep(mu)
     &    *em_q1*tl - 2*pp(nu)*ep(mu)*em_q2*tl + ep(mu)*em(nu)*tl**2 - 
     &    ep(mu)*em(nu)*mw**2*tl - 2*ep(mu)*q1(nu)*pp_em*tl + 2*ep(mu)*
     &    q1(nu)*em_q1*tl + 2*ep(mu)*q1(nu)*em_q2*tl - ep(mu)*q2(nu)*
     &    pp_em*tl + ep(mu)*q2(nu)*em_q1*tl + ep(mu)*q2(nu)*em_q2*tl )
      amp3ax = amp3ax + n_q1*n_q2*qsq2 * ( 2*pp(mu)*em(nu)*pp_ep*tl + 
     &    ep(mu)*em(nu)*tl**2 - ep(mu)*em(nu)*mw**2*tl - em(nu)*q1(mu)*
     &    pp_ep*tl )
      amp3ax = amp3ax + n_q1*n_q2 * (  - 4*pp(mu)*pp(nu)*pp_ep*pp_em*tl
     &     + 4*pp(mu)*pp(nu)*pp_ep*em_q1*tl + 4*pp(mu)*pp(nu)*pp_ep*
     &    em_q2*tl + 4*pp(mu)*pp(nu)*ep_em*tl**2 - 4*pp(mu)*ep(nu)*
     &    em_q2*tl**2 - 2*pp(mu)*em(nu)*pp_ep*tl**2 + 2*pp(mu)*em(nu)*
     &    pp_ep*mw**2*tl + 4*pp(mu)*em(nu)*ep_q2*tl**2 + 4*pp(mu)*
     &    q1(nu)*pp_ep*pp_em*tl - 4*pp(mu)*q1(nu)*pp_ep*em_q1*tl - 4*
     &    pp(mu)*q1(nu)*pp_ep*em_q2*tl - 4*pp(mu)*q1(nu)*ep_em*tl**2 + 
     &    2*pp(mu)*q2(nu)*pp_ep*pp_em*tl - 2*pp(mu)*q2(nu)*pp_ep*em_q1*
     &    tl - 2*pp(mu)*q2(nu)*ep_em*tl**2 - 2*pp(nu)*ep(mu)*pp_em*
     &    tl**2 + 2*pp(nu)*ep(mu)*pp_em*mw**2*tl - 2*pp(nu)*ep(mu)*
     &    em_q1*tl**2 - 2*pp(nu)*ep(mu)*em_q1*mw**2*tl + 2*pp(nu)*
     &    ep(mu)*em_q2*tl**2 - 2*pp(nu)*ep(mu)*em_q2*mw**2*tl + 4*
     &    pp(nu)*em(mu)*ep_q1*tl**2 + 2*pp(nu)*q1(mu)*pp_ep*pp_em*tl - 
     &    2*pp(nu)*q1(mu)*pp_ep*em_q1*tl - 2*pp(nu)*q1(mu)*pp_ep*em_q2*
     &    tl + 2*pp(nu)*q1(mu)*pp_em*ep_q1*tl - 2*pp(nu)*q1(mu)*ep_em*
     &    tl**2 )
      amp3ax = amp3ax + n_q1*n_q2 * (  - 2*pp(nu)*q1(mu)*ep_q1*em_q1*tl
     &     - 2*pp(nu)*q1(mu)*ep_q1*em_q2*tl + 2*ep(mu)*em(nu)*tl**2*ul
     &     + ep(mu)*em(nu)*tl**3 - 2*ep(mu)*em(nu)*mw**2*tl**2 - ep(mu)
     &    *em(nu)*mw**4*tl + 2*ep(mu)*q1(nu)*pp_em*tl**2 - 2*ep(mu)*
     &    q1(nu)*pp_em*mw**2*tl + 2*ep(mu)*q1(nu)*em_q1*tl**2 + 2*
     &    ep(mu)*q1(nu)*em_q1*mw**2*tl + 2*ep(mu)*q1(nu)*em_q2*tl**2 + 
     &    2*ep(mu)*q1(nu)*em_q2*mw**2*tl + ep(mu)*q2(nu)*pp_em*tl**2 - 
     &    ep(mu)*q2(nu)*pp_em*mw**2*tl + ep(mu)*q2(nu)*em_q1*tl**2 + 
     &    ep(mu)*q2(nu)*em_q1*mw**2*tl + 2*ep(nu)*q1(mu)*em_q2*tl**2 - 
     &    4*em(mu)*q1(nu)*ep_q1*tl**2 - 2*em(mu)*q2(nu)*ep_q1*tl**2 + 
     &    em(nu)*q1(mu)*pp_ep*tl**2 - em(nu)*q1(mu)*pp_ep*mw**2*tl + 
     &    em(nu)*q1(mu)*ep_q1*tl**2 - em(nu)*q1(mu)*ep_q1*mw**2*tl - 2*
     &    em(nu)*q1(mu)*ep_q2*tl**2 + 4*em(nu)*q2(mu)*ep_q1*tl**2 - 2*
     &    q1(mu)*q1(nu)*pp_ep*pp_em*tl + 2*q1(mu)*q1(nu)*pp_ep*em_q1*tl
     &     + 2*q1(mu)*q1(nu)*pp_ep*em_q2*tl - 2*q1(mu)*q1(nu)*pp_em*
     &    ep_q1*tl )
      amp3ax = amp3ax + n_q1*n_q2 * ( 2*q1(mu)*q1(nu)*ep_em*tl**2 + 2*
     &    q1(mu)*q1(nu)*ep_q1*em_q1*tl + 2*q1(mu)*q1(nu)*ep_q1*em_q2*tl
     &     - q1(mu)*q2(nu)*pp_ep*pp_em*tl + q1(mu)*q2(nu)*pp_ep*em_q1*
     &    tl - q1(mu)*q2(nu)*pp_em*ep_q1*tl + q1(mu)*q2(nu)*ep_em*tl**2
     &     + q1(mu)*q2(nu)*ep_q1*em_q1*tl + q1(mu)*q2(nu)*ep_q1*em_q2*
     &    tl - 4*d_(mu,nu)*ep_q1*em_q2*tl**2 )
      amp3ax = amp3ax + n(mu)*n(nu) * ( pp_ep*pp_em*tl**3 - 2*pp_ep*
     &    pp_em*mw**2*tl**2 + pp_ep*pp_em*mw**4*tl - pp_ep*em_q1*tl**3
     &     + 2*pp_ep*em_q1*mw**2*tl**2 - pp_ep*em_q1*mw**4*tl - pp_em*
     &    ep_q1*tl**3 + 2*pp_em*ep_q1*mw**2*tl**2 - pp_em*ep_q1*mw**4*
     &    tl - ep_em*tl**4 + 2*ep_em*mw**2*tl**3 - ep_em*mw**4*tl**2 + 
     &    ep_q1*em_q1*tl**3 - 2*ep_q1*em_q1*mw**2*tl**2 + ep_q1*em_q1*
     &    mw**4*tl )
      amp3ax = amp3ax + n(mu)*n_q2*qsq2 * ( em(nu)*pp_ep*tl**2 - em(nu)
     &    *pp_ep*mw**2*tl - em(nu)*ep_q1*tl**2 + em(nu)*ep_q1*mw**2*tl
     &     )
      amp3ax = amp3ax + n(mu)*n_q2 * (  - 2*pp(nu)*pp_ep*pp_em*tl**2 + 
     &    2*pp(nu)*pp_ep*pp_em*mw**2*tl + 2*pp(nu)*pp_ep*em_q1*tl**2 - 
     &    2*pp(nu)*pp_ep*em_q1*mw**2*tl + 2*pp(nu)*pp_ep*em_q2*tl**2 - 
     &    2*pp(nu)*pp_ep*em_q2*mw**2*tl + 2*pp(nu)*pp_em*ep_q1*tl**2 - 
     &    2*pp(nu)*pp_em*ep_q1*mw**2*tl + 2*pp(nu)*ep_em*tl**3 - 2*
     &    pp(nu)*ep_em*mw**2*tl**2 - 2*pp(nu)*ep_q1*em_q1*tl**2 + 2*
     &    pp(nu)*ep_q1*em_q1*mw**2*tl - 2*pp(nu)*ep_q1*em_q2*tl**2 + 2*
     &    pp(nu)*ep_q1*em_q2*mw**2*tl - 2*ep(nu)*em_q2*tl**3 + 2*ep(nu)
     &    *em_q2*mw**2*tl**2 - em(nu)*pp_ep*tl**3 + 2*em(nu)*pp_ep*
     &    mw**2*tl**2 - em(nu)*pp_ep*mw**4*tl + em(nu)*ep_q1*tl**3 - 2*
     &    em(nu)*ep_q1*mw**2*tl**2 + em(nu)*ep_q1*mw**4*tl + 2*em(nu)*
     &    ep_q2*tl**3 - 2*em(nu)*ep_q2*mw**2*tl**2 + 2*q1(nu)*pp_ep*
     &    pp_em*tl**2 - 2*q1(nu)*pp_ep*pp_em*mw**2*tl - 2*q1(nu)*pp_ep*
     &    em_q1*tl**2 + 2*q1(nu)*pp_ep*em_q1*mw**2*tl - 2*q1(nu)*pp_ep*
     &    em_q2*tl**2 + 2*q1(nu)*pp_ep*em_q2*mw**2*tl - 2*q1(nu)*pp_em*
     &    ep_q1*tl**2 )
      amp3ax = amp3ax + n(mu)*n_q2 * ( 2*q1(nu)*pp_em*ep_q1*mw**2*tl - 
     &    2*q1(nu)*ep_em*tl**3 + 2*q1(nu)*ep_em*mw**2*tl**2 + 2*q1(nu)*
     &    ep_q1*em_q1*tl**2 - 2*q1(nu)*ep_q1*em_q1*mw**2*tl + 2*q1(nu)*
     &    ep_q1*em_q2*tl**2 - 2*q1(nu)*ep_q1*em_q2*mw**2*tl + q2(nu)*
     &    pp_ep*pp_em*tl**2 - q2(nu)*pp_ep*pp_em*mw**2*tl - q2(nu)*
     &    pp_ep*em_q1*tl**2 + q2(nu)*pp_ep*em_q1*mw**2*tl - q2(nu)*
     &    pp_em*ep_q1*tl**2 + q2(nu)*pp_em*ep_q1*mw**2*tl - q2(nu)*
     &    ep_em*tl**3 + q2(nu)*ep_em*mw**2*tl**2 + q2(nu)*ep_q1*em_q1*
     &    tl**2 - q2(nu)*ep_q1*em_q1*mw**2*tl )
      amp3ax = amp3ax + n(nu)*n_q1*qsq1 * (  - ep(mu)*pp_em*tl**2 + 
     &    ep(mu)*pp_em*mw**2*tl + ep(mu)*em_q1*tl**2 - ep(mu)*em_q1*
     &    mw**2*tl )
      amp3ax = amp3ax + n(nu)*n_q1 * ( 2*pp(mu)*pp_ep*pp_em*tl**2 - 2*
     &    pp(mu)*pp_ep*pp_em*mw**2*tl - 2*pp(mu)*pp_ep*em_q1*tl**2 + 2*
     &    pp(mu)*pp_ep*em_q1*mw**2*tl - 2*pp(mu)*ep_em*tl**3 + 2*pp(mu)
     &    *ep_em*mw**2*tl**2 + ep(mu)*pp_em*tl**3 - 2*ep(mu)*pp_em*
     &    mw**2*tl**2 + ep(mu)*pp_em*mw**4*tl + ep(mu)*em_q1*tl**3 - 
     &    ep(mu)*em_q1*mw**4*tl - 2*em(mu)*ep_q1*tl**3 + 2*em(mu)*ep_q1
     &    *mw**2*tl**2 - q1(mu)*pp_ep*pp_em*tl**2 + q1(mu)*pp_ep*pp_em*
     &    mw**2*tl + q1(mu)*pp_ep*em_q1*tl**2 - q1(mu)*pp_ep*em_q1*
     &    mw**2*tl - q1(mu)*pp_em*ep_q1*tl**2 + q1(mu)*pp_em*ep_q1*
     &    mw**2*tl + q1(mu)*ep_em*tl**3 - q1(mu)*ep_em*mw**2*tl**2 + 
     &    q1(mu)*ep_q1*em_q1*tl**2 - q1(mu)*ep_q1*em_q1*mw**2*tl )

      amp4ax =
     &  + n_q1*n_q2*qsq1 * ( 2*pp(nu)*em(mu)*pp_ep*ul + ep(nu)*em(mu)*
     &    ul**2 - ep(nu)*em(mu)*mw**2*ul - em(mu)*q2(nu)*pp_ep*ul )
      amp4ax = amp4ax + n_q1*n_q2*qsq2 * ( 2*pp(mu)*ep(nu)*pp_em*ul - 2
     &    *pp(mu)*ep(nu)*em_q1*ul - 2*pp(mu)*ep(nu)*em_q2*ul + ep(nu)*
     &    em(mu)*ul**2 - ep(nu)*em(mu)*mw**2*ul - ep(nu)*q1(mu)*pp_em*
     &    ul + ep(nu)*q1(mu)*em_q1*ul + ep(nu)*q1(mu)*em_q2*ul - 2*
     &    ep(nu)*q2(mu)*pp_em*ul + 2*ep(nu)*q2(mu)*em_q1*ul + 2*ep(nu)*
     &    q2(mu)*em_q2*ul )
      amp4ax = amp4ax + n_q1*n_q2 * (  - 4*pp(mu)*pp(nu)*pp_ep*pp_em*ul
     &     + 4*pp(mu)*pp(nu)*pp_ep*em_q1*ul + 4*pp(mu)*pp(nu)*pp_ep*
     &    em_q2*ul + 4*pp(mu)*pp(nu)*ep_em*ul**2 - 2*pp(mu)*ep(nu)*
     &    pp_em*ul**2 + 2*pp(mu)*ep(nu)*pp_em*mw**2*ul + 2*pp(mu)*
     &    ep(nu)*em_q1*ul**2 - 2*pp(mu)*ep(nu)*em_q1*mw**2*ul - 2*
     &    pp(mu)*ep(nu)*em_q2*ul**2 - 2*pp(mu)*ep(nu)*em_q2*mw**2*ul + 
     &    4*pp(mu)*em(nu)*ep_q2*ul**2 + 2*pp(mu)*q2(nu)*pp_ep*pp_em*ul
     &     - 2*pp(mu)*q2(nu)*pp_ep*em_q1*ul - 2*pp(mu)*q2(nu)*pp_ep*
     &    em_q2*ul + 2*pp(mu)*q2(nu)*pp_em*ep_q2*ul - 2*pp(mu)*q2(nu)*
     &    ep_em*ul**2 - 2*pp(mu)*q2(nu)*ep_q2*em_q1*ul - 2*pp(mu)*
     &    q2(nu)*ep_q2*em_q2*ul - 4*pp(nu)*ep(mu)*em_q1*ul**2 - 2*
     &    pp(nu)*em(mu)*pp_ep*ul**2 + 2*pp(nu)*em(mu)*pp_ep*mw**2*ul + 
     &    4*pp(nu)*em(mu)*ep_q1*ul**2 + 2*pp(nu)*q1(mu)*pp_ep*pp_em*ul
     &     - 2*pp(nu)*q1(mu)*pp_ep*em_q2*ul - 2*pp(nu)*q1(mu)*ep_em*
     &    ul**2 + 4*pp(nu)*q2(mu)*pp_ep*pp_em*ul - 4*pp(nu)*q2(mu)*
     &    pp_ep*em_q1*ul )
      amp4ax = amp4ax + n_q1*n_q2 * (  - 4*pp(nu)*q2(mu)*pp_ep*em_q2*ul
     &     - 4*pp(nu)*q2(mu)*ep_em*ul**2 + 2*ep(mu)*q2(nu)*em_q1*ul**2
     &     + ep(nu)*em(mu)*ul**3 + 2*ep(nu)*em(mu)*tl*ul**2 - 2*ep(nu)*
     &    em(mu)*mw**2*ul**2 - ep(nu)*em(mu)*mw**4*ul + ep(nu)*q1(mu)*
     &    pp_em*ul**2 - ep(nu)*q1(mu)*pp_em*mw**2*ul + ep(nu)*q1(mu)*
     &    em_q2*ul**2 + ep(nu)*q1(mu)*em_q2*mw**2*ul + 2*ep(nu)*q2(mu)*
     &    pp_em*ul**2 - 2*ep(nu)*q2(mu)*pp_em*mw**2*ul + 2*ep(nu)*
     &    q2(mu)*em_q1*ul**2 + 2*ep(nu)*q2(mu)*em_q1*mw**2*ul + 2*
     &    ep(nu)*q2(mu)*em_q2*ul**2 + 2*ep(nu)*q2(mu)*em_q2*mw**2*ul + 
     &    4*em(mu)*q1(nu)*ep_q2*ul**2 + em(mu)*q2(nu)*pp_ep*ul**2 - 
     &    em(mu)*q2(nu)*pp_ep*mw**2*ul - 2*em(mu)*q2(nu)*ep_q1*ul**2 + 
     &    em(mu)*q2(nu)*ep_q2*ul**2 - em(mu)*q2(nu)*ep_q2*mw**2*ul - 2*
     &    em(nu)*q1(mu)*ep_q2*ul**2 - 4*em(nu)*q2(mu)*ep_q2*ul**2 - 
     &    q1(mu)*q2(nu)*pp_ep*pp_em*ul + q1(mu)*q2(nu)*pp_ep*em_q2*ul
     &     - q1(mu)*q2(nu)*pp_em*ep_q2*ul + q1(mu)*q2(nu)*ep_em*ul**2
     &     + q1(mu)*q2(nu)*ep_q2*em_q1*ul )
      amp4ax = amp4ax + n_q1*n_q2 * ( q1(mu)*q2(nu)*ep_q2*em_q2*ul - 2*
     &    q2(mu)*q2(nu)*pp_ep*pp_em*ul + 2*q2(mu)*q2(nu)*pp_ep*em_q1*ul
     &     + 2*q2(mu)*q2(nu)*pp_ep*em_q2*ul - 2*q2(mu)*q2(nu)*pp_em*
     &    ep_q2*ul + 2*q2(mu)*q2(nu)*ep_em*ul**2 + 2*q2(mu)*q2(nu)*
     &    ep_q2*em_q1*ul + 2*q2(mu)*q2(nu)*ep_q2*em_q2*ul - 4*d_(mu,nu)
     &    *ep_q2*em_q1*ul**2 )
      amp4ax = amp4ax + n(mu)*n(nu) * ( pp_ep*pp_em*ul**3 - 2*pp_ep*
     &    pp_em*mw**2*ul**2 + pp_ep*pp_em*mw**4*ul - pp_ep*em_q2*ul**3
     &     + 2*pp_ep*em_q2*mw**2*ul**2 - pp_ep*em_q2*mw**4*ul - pp_em*
     &    ep_q2*ul**3 + 2*pp_em*ep_q2*mw**2*ul**2 - pp_em*ep_q2*mw**4*
     &    ul - ep_em*ul**4 + 2*ep_em*mw**2*ul**3 - ep_em*mw**4*ul**2 + 
     &    ep_q2*em_q2*ul**3 - 2*ep_q2*em_q2*mw**2*ul**2 + ep_q2*em_q2*
     &    mw**4*ul )
      amp4ax = amp4ax + n(mu)*n_q2*qsq2 * (  - ep(nu)*pp_em*ul**2 + 
     &    ep(nu)*pp_em*mw**2*ul + ep(nu)*em_q2*ul**2 - ep(nu)*em_q2*
     &    mw**2*ul )
      amp4ax = amp4ax + n(mu)*n_q2 * ( 2*pp(nu)*pp_ep*pp_em*ul**2 - 2*
     &    pp(nu)*pp_ep*pp_em*mw**2*ul - 2*pp(nu)*pp_ep*em_q2*ul**2 + 2*
     &    pp(nu)*pp_ep*em_q2*mw**2*ul - 2*pp(nu)*ep_em*ul**3 + 2*pp(nu)
     &    *ep_em*mw**2*ul**2 + ep(nu)*pp_em*ul**3 - 2*ep(nu)*pp_em*
     &    mw**2*ul**2 + ep(nu)*pp_em*mw**4*ul + ep(nu)*em_q2*ul**3 - 
     &    ep(nu)*em_q2*mw**4*ul - 2*em(nu)*ep_q2*ul**3 + 2*em(nu)*ep_q2
     &    *mw**2*ul**2 - q2(nu)*pp_ep*pp_em*ul**2 + q2(nu)*pp_ep*pp_em*
     &    mw**2*ul + q2(nu)*pp_ep*em_q2*ul**2 - q2(nu)*pp_ep*em_q2*
     &    mw**2*ul - q2(nu)*pp_em*ep_q2*ul**2 + q2(nu)*pp_em*ep_q2*
     &    mw**2*ul + q2(nu)*ep_em*ul**3 - q2(nu)*ep_em*mw**2*ul**2 + 
     &    q2(nu)*ep_q2*em_q2*ul**2 - q2(nu)*ep_q2*em_q2*mw**2*ul )
      amp4ax = amp4ax + n(nu)*n_q1*qsq1 * ( em(mu)*pp_ep*ul**2 - em(mu)
     &    *pp_ep*mw**2*ul - em(mu)*ep_q2*ul**2 + em(mu)*ep_q2*mw**2*ul
     &     )
      amp4ax = amp4ax + n(nu)*n_q1 * (  - 2*pp(mu)*pp_ep*pp_em*ul**2 + 
     &    2*pp(mu)*pp_ep*pp_em*mw**2*ul + 2*pp(mu)*pp_ep*em_q1*ul**2 - 
     &    2*pp(mu)*pp_ep*em_q1*mw**2*ul + 2*pp(mu)*pp_ep*em_q2*ul**2 - 
     &    2*pp(mu)*pp_ep*em_q2*mw**2*ul + 2*pp(mu)*pp_em*ep_q2*ul**2 - 
     &    2*pp(mu)*pp_em*ep_q2*mw**2*ul + 2*pp(mu)*ep_em*ul**3 - 2*
     &    pp(mu)*ep_em*mw**2*ul**2 - 2*pp(mu)*ep_q2*em_q1*ul**2 + 2*
     &    pp(mu)*ep_q2*em_q1*mw**2*ul - 2*pp(mu)*ep_q2*em_q2*ul**2 + 2*
     &    pp(mu)*ep_q2*em_q2*mw**2*ul - 2*ep(mu)*em_q1*ul**3 + 2*ep(mu)
     &    *em_q1*mw**2*ul**2 - em(mu)*pp_ep*ul**3 + 2*em(mu)*pp_ep*
     &    mw**2*ul**2 - em(mu)*pp_ep*mw**4*ul + 2*em(mu)*ep_q1*ul**3 - 
     &    2*em(mu)*ep_q1*mw**2*ul**2 + em(mu)*ep_q2*ul**3 - 2*em(mu)*
     &    ep_q2*mw**2*ul**2 + em(mu)*ep_q2*mw**4*ul + q1(mu)*pp_ep*
     &    pp_em*ul**2 - q1(mu)*pp_ep*pp_em*mw**2*ul - q1(mu)*pp_ep*
     &    em_q2*ul**2 + q1(mu)*pp_ep*em_q2*mw**2*ul - q1(mu)*pp_em*
     &    ep_q2*ul**2 + q1(mu)*pp_em*ep_q2*mw**2*ul - q1(mu)*ep_em*
     &    ul**3 )
      amp4ax = amp4ax + n(nu)*n_q1 * ( q1(mu)*ep_em*mw**2*ul**2 + 
     &    q1(mu)*ep_q2*em_q2*ul**2 - q1(mu)*ep_q2*em_q2*mw**2*ul + 2*
     &    q2(mu)*pp_ep*pp_em*ul**2 - 2*q2(mu)*pp_ep*pp_em*mw**2*ul - 2*
     &    q2(mu)*pp_ep*em_q1*ul**2 + 2*q2(mu)*pp_ep*em_q1*mw**2*ul - 2*
     &    q2(mu)*pp_ep*em_q2*ul**2 + 2*q2(mu)*pp_ep*em_q2*mw**2*ul - 2*
     &    q2(mu)*pp_em*ep_q2*ul**2 + 2*q2(mu)*pp_em*ep_q2*mw**2*ul - 2*
     &    q2(mu)*ep_em*ul**3 + 2*q2(mu)*ep_em*mw**2*ul**2 + 2*q2(mu)*
     &    ep_q2*em_q1*ul**2 - 2*q2(mu)*ep_q2*em_q1*mw**2*ul + 2*q2(mu)*
     &    ep_q2*em_q2*ul**2 - 2*q2(mu)*ep_q2*em_q2*mw**2*ul )

      amp5ax =
     &  + n_q1*n_q2 * ( ep(mu)*em(nu)*mw**2 + ep(nu)*em(mu)*mw**2 + 2*
     &    d_(mu,nu)*pp_ep*pp_em - 2*d_(mu,nu)*pp_ep*em_q1 - 2*d_(mu,nu)
     &    *pp_ep*em_q2 - 2*d_(mu,nu)*ep_em*mw**2 )
      amp5ax = amp5ax + n(mu)*n(nu) * (  - pp_ep*pp_em*ul - pp_ep*pp_em
     &    *tl + 2*pp_ep*pp_em*mw**2 + pp_ep*em_q1*ul + pp_ep*em_q1*tl
     &     - 2*pp_ep*em_q1*mw**2 + pp_ep*em_q2*ul + pp_ep*em_q2*tl - 2*
     &    pp_ep*em_q2*mw**2 + ep_em*mw**2*ul + ep_em*mw**2*tl - 2*ep_em
     &    *mw**4 + ep_q1*em_q2*mw**2 + ep_q2*em_q1*mw**2 )
      amp5ax = amp5ax + n(mu)*n_q2 * (  - ep(nu)*em_q1*mw**2 - em(nu)*
     &    ep_q1*mw**2 - 2*q1(nu)*pp_ep*pp_em + 2*q1(nu)*pp_ep*em_q1 + 2
     &    *q1(nu)*pp_ep*em_q2 + 2*q1(nu)*ep_em*mw**2 )
      amp5ax = amp5ax + n(nu)*n_q1 * (  - ep(mu)*em_q2*mw**2 - em(mu)*
     &    ep_q2*mw**2 - 2*q2(mu)*pp_ep*pp_em + 2*q2(mu)*pp_ep*em_q1 + 2
     &    *q2(mu)*pp_ep*em_q2 + 2*q2(mu)*ep_em*mw**2 )

      if (.not.n_ind) then
      
      amp3axn =
     &  + n_pp*n_q1*n_q2*qsq1 * ( 2*pp(nu)*ep(mu)*pp_em*tl - 2*pp(nu)*
     &    ep(mu)*em_q1*tl - 2*pp(nu)*ep(mu)*em_q2*tl + ep(mu)*em(nu)*
     &    tl**2 - ep(mu)*em(nu)*mw**2*tl - 2*ep(mu)*q1(nu)*pp_em*tl + 2
     &    *ep(mu)*q1(nu)*em_q1*tl + 2*ep(mu)*q1(nu)*em_q2*tl - ep(mu)*
     &    q2(nu)*pp_em*tl + ep(mu)*q2(nu)*em_q1*tl + ep(mu)*q2(nu)*
     &    em_q2*tl )
      amp3axn = amp3axn + n_pp*n_q1*n_q2*qsq2 * ( 2*pp(mu)*em(nu)*pp_ep
     &    *tl + ep(mu)*em(nu)*tl**2 - ep(mu)*em(nu)*mw**2*tl - em(nu)*
     &    q1(mu)*pp_ep*tl )
      amp3axn = amp3axn + n_pp*n_q1*n_q2 * ( 2*pp(mu)*q2(nu)*pp_ep*
     &    em_q2*tl + 2*pp(nu)*q1(mu)*pp_em*ep_q1*tl - 2*pp(nu)*q1(mu)*
     &    ep_q1*em_q1*tl - 2*pp(nu)*q1(mu)*ep_q1*em_q2*tl + ep(mu)*
     &    q2(nu)*em_q2*tl**2 - ep(mu)*q2(nu)*em_q2*mw**2*tl + em(nu)*
     &    q1(mu)*ep_q1*tl**2 - em(nu)*q1(mu)*ep_q1*mw**2*tl - 2*q1(mu)*
     &    q1(nu)*pp_em*ep_q1*tl + 2*q1(mu)*q1(nu)*ep_q1*em_q1*tl + 2*
     &    q1(mu)*q1(nu)*ep_q1*em_q2*tl - q1(mu)*q2(nu)*pp_ep*em_q2*tl
     &     - q1(mu)*q2(nu)*pp_em*ep_q1*tl + q1(mu)*q2(nu)*ep_q1*em_q1*
     &    tl + q1(mu)*q2(nu)*ep_q1*em_q2*tl )
      amp3axn = amp3axn + n_q1*n_q2**2*qsq1 * (  - 2*ep(mu)*em(nu)*
     &    tl**2 )
      amp3axn = amp3axn + n_q1*n_q2**2 * (  - 2*em(nu)*q1(mu)*ep_q1*
     &    tl**2 )
      amp3axn = amp3axn + n_q1**2*n_q2*qsq1 * (  - 2*pp(nu)*ep(mu)*
     &    pp_em*tl + 2*pp(nu)*ep(mu)*em_q1*tl + 2*pp(nu)*ep(mu)*em_q2*
     &    tl - ep(mu)*em(nu)*tl**2 + ep(mu)*em(nu)*mw**2*tl + 2*ep(mu)*
     &    q1(nu)*pp_em*tl - 2*ep(mu)*q1(nu)*em_q1*tl - 2*ep(mu)*q1(nu)*
     &    em_q2*tl + ep(mu)*q2(nu)*pp_em*tl - ep(mu)*q2(nu)*em_q1*tl - 
     &    ep(mu)*q2(nu)*em_q2*tl )
      amp3axn = amp3axn + n_q1**2*n_q2*qsq2 * (  - 2*pp(mu)*em(nu)*
     &    pp_ep*tl + ep(mu)*em(nu)*tl**2 + ep(mu)*em(nu)*mw**2*tl + 
     &    em(nu)*q1(mu)*pp_ep*tl )
      amp3axn = amp3axn + n_q1**2*n_q2 * (  - 2*pp(mu)*q2(nu)*pp_ep*
     &    em_q2*tl - 2*pp(nu)*q1(mu)*pp_em*ep_q1*tl + 2*pp(nu)*q1(mu)*
     &    ep_q1*em_q1*tl + 2*pp(nu)*q1(mu)*ep_q1*em_q2*tl + ep(mu)*
     &    q2(nu)*em_q2*tl**2 + ep(mu)*q2(nu)*em_q2*mw**2*tl - em(nu)*
     &    q1(mu)*ep_q1*tl**2 + em(nu)*q1(mu)*ep_q1*mw**2*tl + 2*q1(mu)*
     &    q1(nu)*pp_em*ep_q1*tl - 2*q1(mu)*q1(nu)*ep_q1*em_q1*tl - 2*
     &    q1(mu)*q1(nu)*ep_q1*em_q2*tl + q1(mu)*q2(nu)*pp_ep*em_q2*tl
     &     + q1(mu)*q2(nu)*pp_em*ep_q1*tl - q1(mu)*q2(nu)*ep_q1*em_q1*
     &    tl - q1(mu)*q2(nu)*ep_q1*em_q2*tl )
      amp3axn = amp3axn + n(mu)*n_pp*n_q2*qsq2 * ( em(nu)*pp_ep*tl**2
     &     - em(nu)*pp_ep*mw**2*tl - em(nu)*ep_q1*tl**2 + em(nu)*ep_q1*
     &    mw**2*tl )
      amp3axn = amp3axn + n(mu)*n_pp*n_q2 * ( q2(nu)*pp_ep*em_q2*tl**2
     &     - q2(nu)*pp_ep*em_q2*mw**2*tl - q2(nu)*ep_q1*em_q2*tl**2 + 
     &    q2(nu)*ep_q1*em_q2*mw**2*tl )
      amp3axn = amp3axn + n(mu)*n_q1*n_q2*qsq2 * (  - em(nu)*pp_ep*
     &    tl**2 + em(nu)*pp_ep*mw**2*tl - em(nu)*ep_q1*tl**2 - em(nu)*
     &    ep_q1*mw**2*tl )
      amp3axn = amp3axn + n(mu)*n_q1*n_q2 * (  - q2(nu)*pp_ep*em_q2*
     &    tl**2 + q2(nu)*pp_ep*em_q2*mw**2*tl - q2(nu)*ep_q1*em_q2*
     &    tl**2 - q2(nu)*ep_q1*em_q2*mw**2*tl )
      amp3axn = amp3axn + n(nu)*n_pp*n_q1*qsq1 * (  - ep(mu)*pp_em*
     &    tl**2 + ep(mu)*pp_em*mw**2*tl + ep(mu)*em_q1*tl**2 - ep(mu)*
     &    em_q1*mw**2*tl )
      amp3axn = amp3axn + n(nu)*n_pp*n_q1 * (  - q1(mu)*pp_em*ep_q1*
     &    tl**2 + q1(mu)*pp_em*ep_q1*mw**2*tl + q1(mu)*ep_q1*em_q1*
     &    tl**2 - q1(mu)*ep_q1*em_q1*mw**2*tl )
      amp3axn = amp3axn + n(nu)*n_q1*n_q2*qsq1 * ( 2*ep(mu)*em_q2*tl**2
     &     )
      amp3axn = amp3axn + n(nu)*n_q1*n_q2 * ( 2*q1(mu)*ep_q1*em_q2*
     &    tl**2 )
      amp3axn = amp3axn + n(nu)*n_q1**2*qsq1 * ( ep(mu)*pp_em*tl**2 - 
     &    ep(mu)*pp_em*mw**2*tl - ep(mu)*em_q1*tl**2 + ep(mu)*em_q1*
     &    mw**2*tl )
      amp3axn = amp3axn + n(nu)*n_q1**2 * ( q1(mu)*pp_em*ep_q1*tl**2 - 
     &    q1(mu)*pp_em*ep_q1*mw**2*tl - q1(mu)*ep_q1*em_q1*tl**2 + 
     &    q1(mu)*ep_q1*em_q1*mw**2*tl )

      amp4axn =
     &  + n_pp*n_q1*n_q2*qsq1 * (  - 2*pp(nu)*em(mu)*pp_ep*ul - ep(nu)*
     &    em(mu)*ul**2 + ep(nu)*em(mu)*mw**2*ul + em(mu)*q2(nu)*pp_ep*
     &    ul )
      amp4axn = amp4axn + n_pp*n_q1*n_q2*qsq2 * (  - 2*pp(mu)*ep(nu)*
     &    pp_em*ul + 2*pp(mu)*ep(nu)*em_q1*ul + 2*pp(mu)*ep(nu)*em_q2*
     &    ul - ep(nu)*em(mu)*ul**2 + ep(nu)*em(mu)*mw**2*ul + ep(nu)*
     &    q1(mu)*pp_em*ul - ep(nu)*q1(mu)*em_q1*ul - ep(nu)*q1(mu)*
     &    em_q2*ul + 2*ep(nu)*q2(mu)*pp_em*ul - 2*ep(nu)*q2(mu)*em_q1*
     &    ul - 2*ep(nu)*q2(mu)*em_q2*ul )
      amp4axn = amp4axn + n_pp*n_q1*n_q2 * (  - 2*pp(mu)*q2(nu)*pp_em*
     &    ep_q2*ul + 2*pp(mu)*q2(nu)*ep_q2*em_q1*ul + 2*pp(mu)*q2(nu)*
     &    ep_q2*em_q2*ul - 2*pp(nu)*q1(mu)*pp_ep*em_q1*ul - ep(nu)*
     &    q1(mu)*em_q1*ul**2 + ep(nu)*q1(mu)*em_q1*mw**2*ul - em(mu)*
     &    q2(nu)*ep_q2*ul**2 + em(mu)*q2(nu)*ep_q2*mw**2*ul + q1(mu)*
     &    q2(nu)*pp_ep*em_q1*ul + q1(mu)*q2(nu)*pp_em*ep_q2*ul - q1(mu)
     &    *q2(nu)*ep_q2*em_q1*ul - q1(mu)*q2(nu)*ep_q2*em_q2*ul + 2*
     &    q2(mu)*q2(nu)*pp_em*ep_q2*ul - 2*q2(mu)*q2(nu)*ep_q2*em_q1*ul
     &     - 2*q2(mu)*q2(nu)*ep_q2*em_q2*ul )
      amp4axn = amp4axn + n_q1*n_q2**2*qsq1 * ( 2*pp(nu)*em(mu)*pp_ep*
     &    ul - ep(nu)*em(mu)*ul**2 - ep(nu)*em(mu)*mw**2*ul - em(mu)*
     &    q2(nu)*pp_ep*ul )
      amp4axn = amp4axn + n_q1*n_q2**2*qsq2 * ( 2*pp(mu)*ep(nu)*pp_em*
     &    ul - 2*pp(mu)*ep(nu)*em_q1*ul - 2*pp(mu)*ep(nu)*em_q2*ul + 
     &    ep(nu)*em(mu)*ul**2 - ep(nu)*em(mu)*mw**2*ul - ep(nu)*q1(mu)*
     &    pp_em*ul + ep(nu)*q1(mu)*em_q1*ul + ep(nu)*q1(mu)*em_q2*ul - 
     &    2*ep(nu)*q2(mu)*pp_em*ul + 2*ep(nu)*q2(mu)*em_q1*ul + 2*
     &    ep(nu)*q2(mu)*em_q2*ul )
      amp4axn = amp4axn + n_q1*n_q2**2 * ( 2*pp(mu)*q2(nu)*pp_em*ep_q2*
     &    ul - 2*pp(mu)*q2(nu)*ep_q2*em_q1*ul - 2*pp(mu)*q2(nu)*ep_q2*
     &    em_q2*ul + 2*pp(nu)*q1(mu)*pp_ep*em_q1*ul - ep(nu)*q1(mu)*
     &    em_q1*ul**2 - ep(nu)*q1(mu)*em_q1*mw**2*ul + em(mu)*q2(nu)*
     &    ep_q2*ul**2 - em(mu)*q2(nu)*ep_q2*mw**2*ul - q1(mu)*q2(nu)*
     &    pp_ep*em_q1*ul - q1(mu)*q2(nu)*pp_em*ep_q2*ul + q1(mu)*q2(nu)
     &    *ep_q2*em_q1*ul + q1(mu)*q2(nu)*ep_q2*em_q2*ul - 2*q2(mu)*
     &    q2(nu)*pp_em*ep_q2*ul + 2*q2(mu)*q2(nu)*ep_q2*em_q1*ul + 2*
     &    q2(mu)*q2(nu)*ep_q2*em_q2*ul )
      amp4axn = amp4axn + n_q1**2*n_q2*qsq2 * ( 2*ep(nu)*em(mu)*ul**2 )
      amp4axn = amp4axn + n_q1**2*n_q2 * ( 2*em(mu)*q2(nu)*ep_q2*ul**2
     &     )
      amp4axn = amp4axn + n(mu)*n_pp*n_q2*qsq2 * ( ep(nu)*pp_em*ul**2
     &     - ep(nu)*pp_em*mw**2*ul - ep(nu)*em_q2*ul**2 + ep(nu)*em_q2*
     &    mw**2*ul )
      amp4axn = amp4axn + n(mu)*n_pp*n_q2 * ( q2(nu)*pp_em*ep_q2*ul**2
     &     - q2(nu)*pp_em*ep_q2*mw**2*ul - q2(nu)*ep_q2*em_q2*ul**2 + 
     &    q2(nu)*ep_q2*em_q2*mw**2*ul )
      amp4axn = amp4axn + n(mu)*n_q1*n_q2*qsq2 * (  - 2*ep(nu)*em_q1*
     &    ul**2 )
      amp4axn = amp4axn + n(mu)*n_q1*n_q2 * (  - 2*q2(nu)*ep_q2*em_q1*
     &    ul**2 )
      amp4axn = amp4axn + n(mu)*n_q2**2*qsq2 * (  - ep(nu)*pp_em*ul**2
     &     + ep(nu)*pp_em*mw**2*ul + ep(nu)*em_q2*ul**2 - ep(nu)*em_q2*
     &    mw**2*ul )
      amp4axn = amp4axn + n(mu)*n_q2**2 * (  - q2(nu)*pp_em*ep_q2*ul**2
     &     + q2(nu)*pp_em*ep_q2*mw**2*ul + q2(nu)*ep_q2*em_q2*ul**2 - 
     &    q2(nu)*ep_q2*em_q2*mw**2*ul )
      amp4axn = amp4axn + n(nu)*n_pp*n_q1*qsq1 * (  - em(mu)*pp_ep*
     &    ul**2 + em(mu)*pp_ep*mw**2*ul + em(mu)*ep_q2*ul**2 - em(mu)*
     &    ep_q2*mw**2*ul )
      amp4axn = amp4axn + n(nu)*n_pp*n_q1 * (  - q1(mu)*pp_ep*em_q1*
     &    ul**2 + q1(mu)*pp_ep*em_q1*mw**2*ul + q1(mu)*ep_q2*em_q1*
     &    ul**2 - q1(mu)*ep_q2*em_q1*mw**2*ul )
      amp4axn = amp4axn + n(nu)*n_q1*n_q2*qsq1 * ( em(mu)*pp_ep*ul**2
     &     - em(mu)*pp_ep*mw**2*ul + em(mu)*ep_q2*ul**2 + em(mu)*ep_q2*
     &    mw**2*ul )
      amp4axn = amp4axn + n(nu)*n_q1*n_q2 * ( q1(mu)*pp_ep*em_q1*ul**2
     &     - q1(mu)*pp_ep*em_q1*mw**2*ul + q1(mu)*ep_q2*em_q1*ul**2 + 
     &    q1(mu)*ep_q2*em_q1*mw**2*ul )

      amp3axnn =
     &  + n_n*n_q1*n_q2*qsq1 * ( ep(mu)*q2(nu)*em_q2*tl**2 )
      amp3axnn = amp3axnn + n_n*n_q1*n_q2*qsq2 * ( em(nu)*q1(mu)*ep_q1*
     &    tl**2 )
      amp3axnn = amp3axnn + n_n*n_q1*n_q2*qsq2*qsq1 * ( ep(mu)*em(nu)*
     &    tl**2 )
      amp3axnn = amp3axnn + n_n*n_q1*n_q2 * ( q1(mu)*q2(nu)*ep_q1*em_q2
     &    *tl**2 )

      amp4axnn =
     &  + n_n*n_q1*n_q2*qsq1 * ( em(mu)*q2(nu)*ep_q2*ul**2 )
      amp4axnn = amp4axnn + n_n*n_q1*n_q2*qsq2 * ( ep(nu)*q1(mu)*em_q1*
     &    ul**2 )
      amp4axnn = amp4axnn + n_n*n_q1*n_q2*qsq2*qsq1 * ( ep(nu)*em(mu)*
     &    ul**2 )
      amp4axnn = amp4axnn + n_n*n_q1*n_q2 * ( q1(mu)*q2(nu)*ep_q2*em_q1
     &    *ul**2 )
      endif

      if (n_ind) then

      amp3axnind =
     &  + qsq1 * ( n(nu)*ep(mu)*n_q1*pp_em*tl**2 - n(nu)*ep(mu)*n_q1*
     &    pp_em*mw**2*tl - n(nu)*ep(mu)*n_q1*em_q1*tl**2 + n(nu)*ep(mu)
     &    *n_q1*em_q1*mw**2*tl - n(nu)*ep(mu)*n_q1*em_q2*tl**2 - 2*
     &    pp(nu)*ep(mu)*n_q1*n_q2*pp_em*tl + 2*pp(nu)*ep(mu)*n_q1*n_q2*
     &    em_q1*tl + 2*pp(nu)*ep(mu)*n_q1*n_q2*em_q2*tl + ep(mu)*em(nu)
     &    *n_q1*n_q2*mw**2*tl + 2*ep(mu)*q1(nu)*n_q1*n_q2*pp_em*tl - 2*
     &    ep(mu)*q1(nu)*n_q1*n_q2*em_q1*tl - 2*ep(mu)*q1(nu)*n_q1*n_q2*
     &    em_q2*tl + ep(mu)*q2(nu)*n_q1*n_q2*pp_em*tl - ep(mu)*q2(nu)*
     &    n_q1*n_q2*em_q1*tl - ep(mu)*q2(nu)*n_q1*n_q2*em_q2*tl )
      amp3axnind = amp3axnind + qsq2 * (  - n(mu)*em(nu)*n_q2*pp_ep*
     &    tl**2 + n(mu)*em(nu)*n_q2*pp_ep*mw**2*tl - n(mu)*em(nu)*n_q2*
     &    ep_q1*mw**2*tl - 2*pp(mu)*em(nu)*n_q1*n_q2*pp_ep*tl + ep(mu)*
     &    em(nu)*n_q1*n_q2*mw**2*tl + em(nu)*q1(mu)*n_q1*n_q2*pp_ep*tl
     &     )
      amp3axnind = amp3axnind - n(mu)*q2(nu)*n_q2*pp_ep*em_q2*tl**2 + 
     &    n(mu)*q2(nu)*n_q2*pp_ep*em_q2*mw**2*tl - n(mu)*q2(nu)*n_q2*
     &    ep_q1*em_q2*mw**2*tl + n(nu)*q1(mu)*n_q1*pp_em*ep_q1*tl**2 - 
     &    n(nu)*q1(mu)*n_q1*pp_em*ep_q1*mw**2*tl - n(nu)*q1(mu)*n_q1*
     &    ep_q1*em_q1*tl**2 + n(nu)*q1(mu)*n_q1*ep_q1*em_q1*mw**2*tl - 
     &    n(nu)*q1(mu)*n_q1*ep_q1*em_q2*tl**2 - 2*pp(mu)*q2(nu)*n_q1*
     &    n_q2*pp_ep*em_q2*tl - 2*pp(nu)*q1(mu)*n_q1*n_q2*pp_em*ep_q1*
     &    tl + 2*pp(nu)*q1(mu)*n_q1*n_q2*ep_q1*em_q1*tl + 2*pp(nu)*
     &    q1(mu)*n_q1*n_q2*ep_q1*em_q2*tl + ep(mu)*q2(nu)*n_q1*n_q2*
     &    em_q2*mw**2*tl + em(nu)*q1(mu)*n_q1*n_q2*ep_q1*mw**2*tl + 2*
     &    q1(mu)*q1(nu)*n_q1*n_q2*pp_em*ep_q1*tl - 2*q1(mu)*q1(nu)*n_q1
     &    *n_q2*ep_q1*em_q1*tl - 2*q1(mu)*q1(nu)*n_q1*n_q2*ep_q1*em_q2*
     &    tl + q1(mu)*q2(nu)*n_q1*n_q2*pp_ep*em_q2*tl + q1(mu)*q2(nu)*
     &    n_q1*n_q2*pp_em*ep_q1*tl - q1(mu)*q2(nu)*n_q1*n_q2*ep_q1*
     &    em_q1*tl - q1(mu)*q2(nu)*n_q1*n_q2*ep_q1*em_q2*tl

      amp4axnind =
     &  + qsq1 * (  - n(nu)*em(mu)*n_q1*pp_ep*ul**2 + n(nu)*em(mu)*n_q1
     &    *pp_ep*mw**2*ul - n(nu)*em(mu)*n_q1*ep_q2*mw**2*ul - 2*pp(nu)
     &    *em(mu)*n_q1*n_q2*pp_ep*ul + ep(nu)*em(mu)*n_q1*n_q2*mw**2*ul
     &     + em(mu)*q2(nu)*n_q1*n_q2*pp_ep*ul )
      amp4axnind = amp4axnind + qsq2 * ( n(mu)*ep(nu)*n_q2*pp_em*ul**2
     &     - n(mu)*ep(nu)*n_q2*pp_em*mw**2*ul - n(mu)*ep(nu)*n_q2*em_q1
     &    *ul**2 - n(mu)*ep(nu)*n_q2*em_q2*ul**2 + n(mu)*ep(nu)*n_q2*
     &    em_q2*mw**2*ul - 2*pp(mu)*ep(nu)*n_q1*n_q2*pp_em*ul + 2*
     &    pp(mu)*ep(nu)*n_q1*n_q2*em_q1*ul + 2*pp(mu)*ep(nu)*n_q1*n_q2*
     &    em_q2*ul + ep(nu)*em(mu)*n_q1*n_q2*mw**2*ul + ep(nu)*q1(mu)*
     &    n_q1*n_q2*pp_em*ul - ep(nu)*q1(mu)*n_q1*n_q2*em_q1*ul - 
     &    ep(nu)*q1(mu)*n_q1*n_q2*em_q2*ul + 2*ep(nu)*q2(mu)*n_q1*n_q2*
     &    pp_em*ul - 2*ep(nu)*q2(mu)*n_q1*n_q2*em_q1*ul - 2*ep(nu)*
     &    q2(mu)*n_q1*n_q2*em_q2*ul )
      amp4axnind = amp4axnind + n(mu)*q2(nu)*n_q2*pp_em*ep_q2*ul**2 - 
     &    n(mu)*q2(nu)*n_q2*pp_em*ep_q2*mw**2*ul - n(mu)*q2(nu)*n_q2*
     &    ep_q2*em_q1*ul**2 - n(mu)*q2(nu)*n_q2*ep_q2*em_q2*ul**2 + 
     &    n(mu)*q2(nu)*n_q2*ep_q2*em_q2*mw**2*ul - n(nu)*q1(mu)*n_q1*
     &    pp_ep*em_q1*ul**2 + n(nu)*q1(mu)*n_q1*pp_ep*em_q1*mw**2*ul - 
     &    n(nu)*q1(mu)*n_q1*ep_q2*em_q1*mw**2*ul - 2*pp(mu)*q2(nu)*n_q1
     &    *n_q2*pp_em*ep_q2*ul + 2*pp(mu)*q2(nu)*n_q1*n_q2*ep_q2*em_q1*
     &    ul + 2*pp(mu)*q2(nu)*n_q1*n_q2*ep_q2*em_q2*ul - 2*pp(nu)*
     &    q1(mu)*n_q1*n_q2*pp_ep*em_q1*ul + ep(nu)*q1(mu)*n_q1*n_q2*
     &    em_q1*mw**2*ul + em(mu)*q2(nu)*n_q1*n_q2*ep_q2*mw**2*ul + 
     &    q1(mu)*q2(nu)*n_q1*n_q2*pp_ep*em_q1*ul + q1(mu)*q2(nu)*n_q1*
     &    n_q2*pp_em*ep_q2*ul - q1(mu)*q2(nu)*n_q1*n_q2*ep_q2*em_q1*ul
     &     - q1(mu)*q2(nu)*n_q1*n_q2*ep_q2*em_q2*ul + 2*q2(mu)*q2(nu)*
     &    n_q1*n_q2*pp_em*ep_q2*ul - 2*q2(mu)*q2(nu)*n_q1*n_q2*ep_q2*
     &    em_q1*ul - 2*q2(mu)*q2(nu)*n_q1*n_q2*ep_q2*em_q2*ul


      endif
      
      den1ax =
     &  + mw**2*tl

      den2ax =
     &  + mw**2*ul

      den3ax =
     &  + tl**3 - mw**2*tl**2

      den4ax =
     &  + ul**3 - mw**2*ul**2

      den5ax =
     &  + mw**2


      den3axnn =
     &  + n_pp*n_q1 * (  - 2*tl**3 + 2*mw**2*tl**2 )
      den3axnn = den3axnn + n_pp**2 * ( tl**3 - mw**2*tl**2 )
      den3axnn = den3axnn + n_q1**2 * ( tl**3 - mw**2*tl**2 )

      
      den3axn = (n_pp-n_q1)*tl**2*(mw**2-tl);

      
      den4axn = (n_pp-n_q2)*ul**2*(mw**2-ul);


      den4axn =
     &  + n_pp * ( ul**3 - mw**2*ul**2 )
      den4axn = den4axn + n_q2 * (  - ul**3 + mw**2*ul**2 )

      den4axnn =
     &  + n_pp*n_q2 * (  - 2*ul**3 + 2*mw**2*ul**2 )
      den4axnn = den4axnn + n_pp**2 * ( ul**3 - mw**2*ul**2 )
      den4axnn = den4axnn + n_q2**2 * ( ul**3 - mw**2*ul**2 )

      if (.not.n_ind) then
      ampax_fun = amp1ax/den1ax+amp2ax/den2ax+amp3ax/den3ax+
     &     amp4ax/den4ax+amp5ax/den5ax+amp3axn/den3axn+
     &     amp4axn/den4axn+amp3axnn/den3axnn+amp4axnn/den4axnn
      else
      ampax_fun = amp1ax/den1ax+amp2ax/den2ax+amp3ax/den3ax+
     &     amp4ax/den4ax+amp5ax/den5ax+amp3axnind/den3ax+
     &     amp4axnind/den4ax
      endif
      ampax_fun = ampax_fun/(n_q1*n_q2)

      return
      end



      function ampaxial_fun(pp,pm,ep,em,q1,q2,ep_em,ep_q1,ep_q2,
     &              em_q1,em_q2,pp_q1,pp_q2,pm_q1,pm_q2,tl,ul,qsq1,qsq2,
     &     mw,d_,mu,nu,pp_em,pp_ep,n,n_q1,n_q2,n_pp,n_pm,n_n)
      implicit none
      complex*16 ep_em,ep_q1,ep_q2,em_q1,em_q2,pp_em,pp_ep,ep(4),em(4)
     &     ,d_(4,4),amp1,amp2,amp3,amp4,amp5,ampaxial_fun
      double precision tl,ul,qsq1,pp(4),pm(4),q1(4),q2(4),
     &     qsq2,mw,pp_q1,pp_q2,pm_q1,pm_q2,n_q1,n_q2
     &     ,n(4),n_pp,n_pm,n_n,den1,den2,den3,den4,den5
      integer mu,nu,nu2,mu2
      logical test,n_ind

      amp1 =
     &  - 4*pp(mu)*pp(nu)*pp_ep*pp_em + 4*pp(mu)*pp(nu)*pp_ep*em_q1 + 4
     &    *pp(mu)*pp(nu)*pp_ep*em_q2 + 2*pp(mu)*em(nu)*pp_ep*mw**2 + 4*
     &    pp(mu)*q1(nu)*pp_ep*pp_em - 4*pp(mu)*q1(nu)*pp_ep*em_q1 - 4*
     &    pp(mu)*q1(nu)*pp_ep*em_q2 + 2*pp(mu)*q2(nu)*pp_ep*pp_em - 2*
     &    pp(mu)*q2(nu)*pp_ep*em_q1 - 2*pp(mu)*q2(nu)*pp_ep*em_q2 + 2*
     &    pp(nu)*ep(mu)*pp_em*mw**2 - 2*pp(nu)*ep(mu)*em_q1*mw**2 - 2*
     &    pp(nu)*ep(mu)*em_q2*mw**2 + 2*pp(nu)*q1(mu)*pp_ep*pp_em - 2*
     &    pp(nu)*q1(mu)*pp_ep*em_q1 - 2*pp(nu)*q1(mu)*pp_ep*em_q2 - 
     &    ep(mu)*em(nu)*mw**4 - 2*ep(mu)*q1(nu)*pp_em*mw**2 + 2*ep(mu)*
     &    q1(nu)*em_q1*mw**2 + 2*ep(mu)*q1(nu)*em_q2*mw**2 - ep(mu)*
     &    q2(nu)*pp_em*mw**2 + ep(mu)*q2(nu)*em_q1*mw**2 + ep(mu)*
     &    q2(nu)*em_q2*mw**2 - em(nu)*q1(mu)*pp_ep*mw**2 - 2*q1(mu)*
     &    q1(nu)*pp_ep*pp_em + 2*q1(mu)*q1(nu)*pp_ep*em_q1 + 2*q1(mu)*
     &    q1(nu)*pp_ep*em_q2 - q1(mu)*q2(nu)*pp_ep*pp_em + q1(mu)*
     &    q2(nu)*pp_ep*em_q1
      amp1 = amp1 + q1(mu)*q2(nu)*pp_ep*em_q2

      den1 =
     &  + mw**2*tl

      amp2 =
     &  - 4*pp(mu)*pp(nu)*pp_ep*pp_em + 4*pp(mu)*pp(nu)*pp_ep*em_q1 + 4
     &    *pp(mu)*pp(nu)*pp_ep*em_q2 + 2*pp(mu)*ep(nu)*pp_em*mw**2 - 2*
     &    pp(mu)*ep(nu)*em_q1*mw**2 - 2*pp(mu)*ep(nu)*em_q2*mw**2 + 2*
     &    pp(mu)*q2(nu)*pp_ep*pp_em - 2*pp(mu)*q2(nu)*pp_ep*em_q1 - 2*
     &    pp(mu)*q2(nu)*pp_ep*em_q2 + 2*pp(nu)*em(mu)*pp_ep*mw**2 + 2*
     &    pp(nu)*q1(mu)*pp_ep*pp_em - 2*pp(nu)*q1(mu)*pp_ep*em_q1 - 2*
     &    pp(nu)*q1(mu)*pp_ep*em_q2 + 4*pp(nu)*q2(mu)*pp_ep*pp_em - 4*
     &    pp(nu)*q2(mu)*pp_ep*em_q1 - 4*pp(nu)*q2(mu)*pp_ep*em_q2 - 
     &    ep(nu)*em(mu)*mw**4 - ep(nu)*q1(mu)*pp_em*mw**2 + ep(nu)*
     &    q1(mu)*em_q1*mw**2 + ep(nu)*q1(mu)*em_q2*mw**2 - 2*ep(nu)*
     &    q2(mu)*pp_em*mw**2 + 2*ep(nu)*q2(mu)*em_q1*mw**2 + 2*ep(nu)*
     &    q2(mu)*em_q2*mw**2 - em(mu)*q2(nu)*pp_ep*mw**2 - q1(mu)*
     &    q2(nu)*pp_ep*pp_em + q1(mu)*q2(nu)*pp_ep*em_q1 + q1(mu)*
     &    q2(nu)*pp_ep*em_q2 - 2*q2(mu)*q2(nu)*pp_ep*pp_em + 2*q2(mu)*
     &    q2(nu)*pp_ep*em_q1
      amp2 = amp2 + 2*q2(mu)*q2(nu)*pp_ep*em_q2

      
      den2 =
     &  + mw**2*ul

      amp3 =
     &  + n_n*qsq1 * ( ep(mu)*q2(nu)*em_q2*tl**2 )
      amp3 = amp3 + n_n*qsq2 * ( em(nu)*q1(mu)*ep_q1*tl**2 )
      amp3 = amp3 + n_n*qsq2*qsq1 * ( ep(mu)*em(nu)*tl**2 )
      amp3 = amp3 + n_n * ( q1(mu)*q2(nu)*ep_q1*em_q2*tl**2 )
      amp3 = amp3 + n_pp*n_q1*qsq2 * (  - 2*ep(mu)*em(nu)*tl**2 )
      amp3 = amp3 + n_pp*n_q1 * ( 8*pp(mu)*pp(nu)*pp_ep*pp_em*tl - 8*
     &    pp(mu)*pp(nu)*pp_ep*em_q1*tl - 8*pp(mu)*pp(nu)*pp_ep*em_q2*tl
     &     - 8*pp(mu)*pp(nu)*ep_em*tl**2 + 8*pp(mu)*ep(nu)*em_q2*tl**2
     &     + 4*pp(mu)*em(nu)*pp_ep*tl**2 - 4*pp(mu)*em(nu)*pp_ep*mw**2*
     &    tl - 8*pp(mu)*em(nu)*ep_q2*tl**2 - 8*pp(mu)*q1(nu)*pp_ep*
     &    pp_em*tl + 8*pp(mu)*q1(nu)*pp_ep*em_q1*tl + 8*pp(mu)*q1(nu)*
     &    pp_ep*em_q2*tl + 8*pp(mu)*q1(nu)*ep_em*tl**2 - 4*pp(mu)*
     &    q2(nu)*pp_ep*pp_em*tl + 4*pp(mu)*q2(nu)*pp_ep*em_q1*tl + 4*
     &    pp(mu)*q2(nu)*pp_ep*em_q2*tl + 4*pp(mu)*q2(nu)*ep_em*tl**2 + 
     &    4*pp(nu)*ep(mu)*pp_em*tl**2 - 4*pp(nu)*ep(mu)*pp_em*mw**2*tl
     &     + 4*pp(nu)*ep(mu)*em_q1*tl**2 + 4*pp(nu)*ep(mu)*em_q1*mw**2*
     &    tl - 4*pp(nu)*ep(mu)*em_q2*tl**2 + 4*pp(nu)*ep(mu)*em_q2*
     &    mw**2*tl - 8*pp(nu)*em(mu)*ep_q1*tl**2 - 4*pp(nu)*q1(mu)*
     &    pp_ep*pp_em*tl + 4*pp(nu)*q1(mu)*pp_ep*em_q1*tl + 4*pp(nu)*
     &    q1(mu)*pp_ep*em_q2*tl + 4*pp(nu)*q1(mu)*ep_em*tl**2 - 4*
     &    ep(mu)*em(nu)*tl**2*ul )
      amp3 = amp3 + n_pp*n_q1 * (  - 2*ep(mu)*em(nu)*tl**3 + 4*ep(mu)*
     &    em(nu)*mw**2*tl**2 + 2*ep(mu)*em(nu)*mw**4*tl - 4*ep(mu)*
     &    q1(nu)*pp_em*tl**2 + 4*ep(mu)*q1(nu)*pp_em*mw**2*tl - 4*
     &    ep(mu)*q1(nu)*em_q1*tl**2 - 4*ep(mu)*q1(nu)*em_q1*mw**2*tl - 
     &    4*ep(mu)*q1(nu)*em_q2*tl**2 - 4*ep(mu)*q1(nu)*em_q2*mw**2*tl
     &     - 2*ep(mu)*q2(nu)*pp_em*tl**2 + 2*ep(mu)*q2(nu)*pp_em*mw**2*
     &    tl - 2*ep(mu)*q2(nu)*em_q1*tl**2 - 2*ep(mu)*q2(nu)*em_q1*
     &    mw**2*tl - 2*ep(mu)*q2(nu)*em_q2*mw**2*tl - 4*ep(nu)*q1(mu)*
     &    em_q2*tl**2 + 8*em(mu)*q1(nu)*ep_q1*tl**2 + 4*em(mu)*q2(nu)*
     &    ep_q1*tl**2 - 2*em(nu)*q1(mu)*pp_ep*tl**2 + 2*em(nu)*q1(mu)*
     &    pp_ep*mw**2*tl + 4*em(nu)*q1(mu)*ep_q2*tl**2 - 8*em(nu)*
     &    q2(mu)*ep_q1*tl**2 + 4*q1(mu)*q1(nu)*pp_ep*pp_em*tl - 4*
     &    q1(mu)*q1(nu)*pp_ep*em_q1*tl - 4*q1(mu)*q1(nu)*pp_ep*em_q2*tl
     &     - 4*q1(mu)*q1(nu)*ep_em*tl**2 + 2*q1(mu)*q2(nu)*pp_ep*pp_em*
     &    tl - 2*q1(mu)*q2(nu)*pp_ep*em_q1*tl - 2*q1(mu)*q2(nu)*pp_ep*
     &    em_q2*tl )
      amp3 = amp3 + n_pp*n_q1 * (  - 2*q1(mu)*q2(nu)*ep_em*tl**2 + 8*
     &    d_(mu,nu)*ep_q1*em_q2*tl**2 )
      amp3 = amp3 + n_pp*n_q2*qsq1 * ( 2*ep(mu)*em(nu)*tl**2 )
      amp3 = amp3 + n_pp*n_q2 * ( 2*em(nu)*q1(mu)*ep_q1*tl**2 )
      amp3 = amp3 + n_pp**2 * (  - 4*pp(mu)*pp(nu)*pp_ep*pp_em*tl + 4*
     &    pp(mu)*pp(nu)*pp_ep*em_q1*tl + 4*pp(mu)*pp(nu)*pp_ep*em_q2*tl
     &     + 4*pp(mu)*pp(nu)*ep_em*tl**2 - 4*pp(mu)*ep(nu)*em_q2*tl**2
     &     - 2*pp(mu)*em(nu)*pp_ep*tl**2 + 2*pp(mu)*em(nu)*pp_ep*mw**2*
     &    tl + 4*pp(mu)*em(nu)*ep_q2*tl**2 + 4*pp(mu)*q1(nu)*pp_ep*
     &    pp_em*tl - 4*pp(mu)*q1(nu)*pp_ep*em_q1*tl - 4*pp(mu)*q1(nu)*
     &    pp_ep*em_q2*tl - 4*pp(mu)*q1(nu)*ep_em*tl**2 + 2*pp(mu)*
     &    q2(nu)*pp_ep*pp_em*tl - 2*pp(mu)*q2(nu)*pp_ep*em_q1*tl - 2*
     &    pp(mu)*q2(nu)*pp_ep*em_q2*tl - 2*pp(mu)*q2(nu)*ep_em*tl**2 - 
     &    2*pp(nu)*ep(mu)*pp_em*tl**2 + 2*pp(nu)*ep(mu)*pp_em*mw**2*tl
     &     - 2*pp(nu)*ep(mu)*em_q1*tl**2 - 2*pp(nu)*ep(mu)*em_q1*mw**2*
     &    tl + 2*pp(nu)*ep(mu)*em_q2*tl**2 - 2*pp(nu)*ep(mu)*em_q2*
     &    mw**2*tl + 4*pp(nu)*em(mu)*ep_q1*tl**2 + 2*pp(nu)*q1(mu)*
     &    pp_ep*pp_em*tl - 2*pp(nu)*q1(mu)*pp_ep*em_q1*tl - 2*pp(nu)*
     &    q1(mu)*pp_ep*em_q2*tl - 2*pp(nu)*q1(mu)*ep_em*tl**2 + 2*
     &    ep(mu)*em(nu)*tl**2*ul )
      amp3 = amp3 + n_pp**2 * ( ep(mu)*em(nu)*tl**3 - 2*ep(mu)*em(nu)*
     &    mw**2*tl**2 - ep(mu)*em(nu)*mw**4*tl + 2*ep(mu)*q1(nu)*pp_em*
     &    tl**2 - 2*ep(mu)*q1(nu)*pp_em*mw**2*tl + 2*ep(mu)*q1(nu)*
     &    em_q1*tl**2 + 2*ep(mu)*q1(nu)*em_q1*mw**2*tl + 2*ep(mu)*
     &    q1(nu)*em_q2*tl**2 + 2*ep(mu)*q1(nu)*em_q2*mw**2*tl + ep(mu)*
     &    q2(nu)*pp_em*tl**2 - ep(mu)*q2(nu)*pp_em*mw**2*tl + ep(mu)*
     &    q2(nu)*em_q1*tl**2 + ep(mu)*q2(nu)*em_q1*mw**2*tl - ep(mu)*
     &    q2(nu)*em_q2*tl**2 + ep(mu)*q2(nu)*em_q2*mw**2*tl + 2*ep(nu)*
     &    q1(mu)*em_q2*tl**2 - 4*em(mu)*q1(nu)*ep_q1*tl**2 - 2*em(mu)*
     &    q2(nu)*ep_q1*tl**2 + em(nu)*q1(mu)*pp_ep*tl**2 - em(nu)*
     &    q1(mu)*pp_ep*mw**2*tl - 2*em(nu)*q1(mu)*ep_q2*tl**2 + 4*
     &    em(nu)*q2(mu)*ep_q1*tl**2 - 2*q1(mu)*q1(nu)*pp_ep*pp_em*tl + 
     &    2*q1(mu)*q1(nu)*pp_ep*em_q1*tl + 2*q1(mu)*q1(nu)*pp_ep*em_q2*
     &    tl + 2*q1(mu)*q1(nu)*ep_em*tl**2 - q1(mu)*q2(nu)*pp_ep*pp_em*
     &    tl + q1(mu)*q2(nu)*pp_ep*em_q1*tl + q1(mu)*q2(nu)*pp_ep*em_q2
     &    *tl )
      amp3 = amp3 + n_pp**2 * ( q1(mu)*q2(nu)*ep_em*tl**2 - 4*d_(mu,nu)
     &    *ep_q1*em_q2*tl**2 )
      amp3 = amp3 + n_q1*n_q2*qsq1 * (  - 2*ep(mu)*em(nu)*tl**2 )
      amp3 = amp3 + n_q1*n_q2 * (  - 2*em(nu)*q1(mu)*ep_q1*tl**2 )
      amp3 = amp3 + n_q1**2*qsq2 * ( 2*ep(mu)*em(nu)*tl**2 )
      amp3 = amp3 + n_q1**2 * (  - 4*pp(mu)*pp(nu)*pp_ep*pp_em*tl + 4*
     &    pp(mu)*pp(nu)*pp_ep*em_q1*tl + 4*pp(mu)*pp(nu)*pp_ep*em_q2*tl
     &     + 4*pp(mu)*pp(nu)*ep_em*tl**2 - 4*pp(mu)*ep(nu)*em_q2*tl**2
     &     - 2*pp(mu)*em(nu)*pp_ep*tl**2 + 2*pp(mu)*em(nu)*pp_ep*mw**2*
     &    tl + 4*pp(mu)*em(nu)*ep_q2*tl**2 + 4*pp(mu)*q1(nu)*pp_ep*
     &    pp_em*tl - 4*pp(mu)*q1(nu)*pp_ep*em_q1*tl - 4*pp(mu)*q1(nu)*
     &    pp_ep*em_q2*tl - 4*pp(mu)*q1(nu)*ep_em*tl**2 + 2*pp(mu)*
     &    q2(nu)*pp_ep*pp_em*tl - 2*pp(mu)*q2(nu)*pp_ep*em_q1*tl - 2*
     &    pp(mu)*q2(nu)*pp_ep*em_q2*tl - 2*pp(mu)*q2(nu)*ep_em*tl**2 - 
     &    2*pp(nu)*ep(mu)*pp_em*tl**2 + 2*pp(nu)*ep(mu)*pp_em*mw**2*tl
     &     - 2*pp(nu)*ep(mu)*em_q1*tl**2 - 2*pp(nu)*ep(mu)*em_q1*mw**2*
     &    tl + 2*pp(nu)*ep(mu)*em_q2*tl**2 - 2*pp(nu)*ep(mu)*em_q2*
     &    mw**2*tl + 4*pp(nu)*em(mu)*ep_q1*tl**2 + 2*pp(nu)*q1(mu)*
     &    pp_ep*pp_em*tl - 2*pp(nu)*q1(mu)*pp_ep*em_q1*tl - 2*pp(nu)*
     &    q1(mu)*pp_ep*em_q2*tl - 2*pp(nu)*q1(mu)*ep_em*tl**2 + 2*
     &    ep(mu)*em(nu)*tl**2*ul )
      amp3 = amp3 + n_q1**2 * ( ep(mu)*em(nu)*tl**3 - 2*ep(mu)*em(nu)*
     &    mw**2*tl**2 - ep(mu)*em(nu)*mw**4*tl + 2*ep(mu)*q1(nu)*pp_em*
     &    tl**2 - 2*ep(mu)*q1(nu)*pp_em*mw**2*tl + 2*ep(mu)*q1(nu)*
     &    em_q1*tl**2 + 2*ep(mu)*q1(nu)*em_q1*mw**2*tl + 2*ep(mu)*
     &    q1(nu)*em_q2*tl**2 + 2*ep(mu)*q1(nu)*em_q2*mw**2*tl + ep(mu)*
     &    q2(nu)*pp_em*tl**2 - ep(mu)*q2(nu)*pp_em*mw**2*tl + ep(mu)*
     &    q2(nu)*em_q1*tl**2 + ep(mu)*q2(nu)*em_q1*mw**2*tl + ep(mu)*
     &    q2(nu)*em_q2*tl**2 + ep(mu)*q2(nu)*em_q2*mw**2*tl + 2*ep(nu)*
     &    q1(mu)*em_q2*tl**2 - 4*em(mu)*q1(nu)*ep_q1*tl**2 - 2*em(mu)*
     &    q2(nu)*ep_q1*tl**2 + em(nu)*q1(mu)*pp_ep*tl**2 - em(nu)*
     &    q1(mu)*pp_ep*mw**2*tl - 2*em(nu)*q1(mu)*ep_q2*tl**2 + 4*
     &    em(nu)*q2(mu)*ep_q1*tl**2 - 2*q1(mu)*q1(nu)*pp_ep*pp_em*tl + 
     &    2*q1(mu)*q1(nu)*pp_ep*em_q1*tl + 2*q1(mu)*q1(nu)*pp_ep*em_q2*
     &    tl + 2*q1(mu)*q1(nu)*ep_em*tl**2 - q1(mu)*q2(nu)*pp_ep*pp_em*
     &    tl + q1(mu)*q2(nu)*pp_ep*em_q1*tl + q1(mu)*q2(nu)*pp_ep*em_q2
     &    *tl )
      amp3 = amp3 + n_q1**2 * ( q1(mu)*q2(nu)*ep_em*tl**2 - 4*d_(mu,nu)
     &    *ep_q1*em_q2*tl**2 )
      amp3 = amp3 + n(mu)*n_pp*qsq2 * ( 2*em(nu)*ep_q1*tl**2 )
      amp3 = amp3 + n(mu)*n_pp * ( 2*q2(nu)*ep_q1*em_q2*tl**2 )
      amp3 = amp3 + n(mu)*n_q1*qsq2 * (  - 2*em(nu)*ep_q1*tl**2 )
      amp3 = amp3 + n(mu)*n_q1 * (  - 2*q2(nu)*ep_q1*em_q2*tl**2 )
      amp3 = amp3 + n(nu)*n_pp*qsq1 * (  - 2*ep(mu)*em_q2*tl**2 )
      amp3 = amp3 + n(nu)*n_pp * (  - 2*q1(mu)*ep_q1*em_q2*tl**2 )
      amp3 = amp3 + n(nu)*n_q1*qsq1 * ( 2*ep(mu)*em_q2*tl**2 )
      amp3 = amp3 + n(nu)*n_q1 * ( 2*q1(mu)*ep_q1*em_q2*tl**2 )

      
      den3 =
     &  + n_pp*n_q1 * (  - 2*tl**3 + 2*mw**2*tl**2 )
      den3 = den3 + n_pp**2 * ( tl**3 - mw**2*tl**2 )
      den3 = den3 + n_q1**2 * ( tl**3 - mw**2*tl**2 )

      
      amp4 =
     &  + n_n*qsq1 * ( em(mu)*q2(nu)*ep_q2*ul**2 )
      amp4 = amp4 + n_n*qsq2 * ( ep(nu)*q1(mu)*em_q1*ul**2 )
      amp4 = amp4 + n_n*qsq2*qsq1 * ( ep(nu)*em(mu)*ul**2 )
      amp4 = amp4 + n_n * ( q1(mu)*q2(nu)*ep_q2*em_q1*ul**2 )
      amp4 = amp4 + n_pp*n_q1*qsq2 * ( 2*ep(nu)*em(mu)*ul**2 )
      amp4 = amp4 + n_pp*n_q1 * ( 2*em(mu)*q2(nu)*ep_q2*ul**2 )
      amp4 = amp4 + n_pp*n_q2*qsq1 * (  - 2*ep(nu)*em(mu)*ul**2 )
      amp4 = amp4 + n_pp*n_q2 * ( 8*pp(mu)*pp(nu)*pp_ep*pp_em*ul - 8*
     &    pp(mu)*pp(nu)*pp_ep*em_q1*ul - 8*pp(mu)*pp(nu)*pp_ep*em_q2*ul
     &     - 8*pp(mu)*pp(nu)*ep_em*ul**2 + 4*pp(mu)*ep(nu)*pp_em*ul**2
     &     - 4*pp(mu)*ep(nu)*pp_em*mw**2*ul - 4*pp(mu)*ep(nu)*em_q1*
     &    ul**2 + 4*pp(mu)*ep(nu)*em_q1*mw**2*ul + 4*pp(mu)*ep(nu)*
     &    em_q2*ul**2 + 4*pp(mu)*ep(nu)*em_q2*mw**2*ul - 8*pp(mu)*
     &    em(nu)*ep_q2*ul**2 - 4*pp(mu)*q2(nu)*pp_ep*pp_em*ul + 4*
     &    pp(mu)*q2(nu)*pp_ep*em_q1*ul + 4*pp(mu)*q2(nu)*pp_ep*em_q2*ul
     &     + 4*pp(mu)*q2(nu)*ep_em*ul**2 + 8*pp(nu)*ep(mu)*em_q1*ul**2
     &     + 4*pp(nu)*em(mu)*pp_ep*ul**2 - 4*pp(nu)*em(mu)*pp_ep*mw**2*
     &    ul - 8*pp(nu)*em(mu)*ep_q1*ul**2 - 4*pp(nu)*q1(mu)*pp_ep*
     &    pp_em*ul + 4*pp(nu)*q1(mu)*pp_ep*em_q1*ul + 4*pp(nu)*q1(mu)*
     &    pp_ep*em_q2*ul + 4*pp(nu)*q1(mu)*ep_em*ul**2 - 8*pp(nu)*
     &    q2(mu)*pp_ep*pp_em*ul + 8*pp(nu)*q2(mu)*pp_ep*em_q1*ul + 8*
     &    pp(nu)*q2(mu)*pp_ep*em_q2*ul + 8*pp(nu)*q2(mu)*ep_em*ul**2 - 
     &    4*ep(mu)*q2(nu)*em_q1*ul**2 )
      amp4 = amp4 + n_pp*n_q2 * (  - 2*ep(nu)*em(mu)*ul**3 - 4*ep(nu)*
     &    em(mu)*tl*ul**2 + 4*ep(nu)*em(mu)*mw**2*ul**2 + 2*ep(nu)*
     &    em(mu)*mw**4*ul - 2*ep(nu)*q1(mu)*pp_em*ul**2 + 2*ep(nu)*
     &    q1(mu)*pp_em*mw**2*ul - 2*ep(nu)*q1(mu)*em_q1*mw**2*ul - 2*
     &    ep(nu)*q1(mu)*em_q2*ul**2 - 2*ep(nu)*q1(mu)*em_q2*mw**2*ul - 
     &    4*ep(nu)*q2(mu)*pp_em*ul**2 + 4*ep(nu)*q2(mu)*pp_em*mw**2*ul
     &     - 4*ep(nu)*q2(mu)*em_q1*ul**2 - 4*ep(nu)*q2(mu)*em_q1*mw**2*
     &    ul - 4*ep(nu)*q2(mu)*em_q2*ul**2 - 4*ep(nu)*q2(mu)*em_q2*
     &    mw**2*ul - 8*em(mu)*q1(nu)*ep_q2*ul**2 - 2*em(mu)*q2(nu)*
     &    pp_ep*ul**2 + 2*em(mu)*q2(nu)*pp_ep*mw**2*ul + 4*em(mu)*
     &    q2(nu)*ep_q1*ul**2 + 4*em(nu)*q1(mu)*ep_q2*ul**2 + 8*em(nu)*
     &    q2(mu)*ep_q2*ul**2 + 2*q1(mu)*q2(nu)*pp_ep*pp_em*ul - 2*
     &    q1(mu)*q2(nu)*pp_ep*em_q1*ul - 2*q1(mu)*q2(nu)*pp_ep*em_q2*ul
     &     - 2*q1(mu)*q2(nu)*ep_em*ul**2 + 4*q2(mu)*q2(nu)*pp_ep*pp_em*
     &    ul - 4*q2(mu)*q2(nu)*pp_ep*em_q1*ul - 4*q2(mu)*q2(nu)*pp_ep*
     &    em_q2*ul )
      amp4 = amp4 + n_pp*n_q2 * (  - 4*q2(mu)*q2(nu)*ep_em*ul**2 + 8*
     &    d_(mu,nu)*ep_q2*em_q1*ul**2 )
      amp4 = amp4 + n_pp**2 * (  - 4*pp(mu)*pp(nu)*pp_ep*pp_em*ul + 4*
     &    pp(mu)*pp(nu)*pp_ep*em_q1*ul + 4*pp(mu)*pp(nu)*pp_ep*em_q2*ul
     &     + 4*pp(mu)*pp(nu)*ep_em*ul**2 - 2*pp(mu)*ep(nu)*pp_em*ul**2
     &     + 2*pp(mu)*ep(nu)*pp_em*mw**2*ul + 2*pp(mu)*ep(nu)*em_q1*
     &    ul**2 - 2*pp(mu)*ep(nu)*em_q1*mw**2*ul - 2*pp(mu)*ep(nu)*
     &    em_q2*ul**2 - 2*pp(mu)*ep(nu)*em_q2*mw**2*ul + 4*pp(mu)*
     &    em(nu)*ep_q2*ul**2 + 2*pp(mu)*q2(nu)*pp_ep*pp_em*ul - 2*
     &    pp(mu)*q2(nu)*pp_ep*em_q1*ul - 2*pp(mu)*q2(nu)*pp_ep*em_q2*ul
     &     - 2*pp(mu)*q2(nu)*ep_em*ul**2 - 4*pp(nu)*ep(mu)*em_q1*ul**2
     &     - 2*pp(nu)*em(mu)*pp_ep*ul**2 + 2*pp(nu)*em(mu)*pp_ep*mw**2*
     &    ul + 4*pp(nu)*em(mu)*ep_q1*ul**2 + 2*pp(nu)*q1(mu)*pp_ep*
     &    pp_em*ul - 2*pp(nu)*q1(mu)*pp_ep*em_q1*ul - 2*pp(nu)*q1(mu)*
     &    pp_ep*em_q2*ul - 2*pp(nu)*q1(mu)*ep_em*ul**2 + 4*pp(nu)*
     &    q2(mu)*pp_ep*pp_em*ul - 4*pp(nu)*q2(mu)*pp_ep*em_q1*ul - 4*
     &    pp(nu)*q2(mu)*pp_ep*em_q2*ul - 4*pp(nu)*q2(mu)*ep_em*ul**2 + 
     &    2*ep(mu)*q2(nu)*em_q1*ul**2 )
      amp4 = amp4 + n_pp**2 * ( ep(nu)*em(mu)*ul**3 + 2*ep(nu)*em(mu)*
     &    tl*ul**2 - 2*ep(nu)*em(mu)*mw**2*ul**2 - ep(nu)*em(mu)*mw**4*
     &    ul + ep(nu)*q1(mu)*pp_em*ul**2 - ep(nu)*q1(mu)*pp_em*mw**2*ul
     &     - ep(nu)*q1(mu)*em_q1*ul**2 + ep(nu)*q1(mu)*em_q1*mw**2*ul
     &     + ep(nu)*q1(mu)*em_q2*ul**2 + ep(nu)*q1(mu)*em_q2*mw**2*ul
     &     + 2*ep(nu)*q2(mu)*pp_em*ul**2 - 2*ep(nu)*q2(mu)*pp_em*mw**2*
     &    ul + 2*ep(nu)*q2(mu)*em_q1*ul**2 + 2*ep(nu)*q2(mu)*em_q1*
     &    mw**2*ul + 2*ep(nu)*q2(mu)*em_q2*ul**2 + 2*ep(nu)*q2(mu)*
     &    em_q2*mw**2*ul + 4*em(mu)*q1(nu)*ep_q2*ul**2 + em(mu)*q2(nu)*
     &    pp_ep*ul**2 - em(mu)*q2(nu)*pp_ep*mw**2*ul - 2*em(mu)*q2(nu)*
     &    ep_q1*ul**2 - 2*em(nu)*q1(mu)*ep_q2*ul**2 - 4*em(nu)*q2(mu)*
     &    ep_q2*ul**2 - q1(mu)*q2(nu)*pp_ep*pp_em*ul + q1(mu)*q2(nu)*
     &    pp_ep*em_q1*ul + q1(mu)*q2(nu)*pp_ep*em_q2*ul + q1(mu)*q2(nu)
     &    *ep_em*ul**2 - 2*q2(mu)*q2(nu)*pp_ep*pp_em*ul + 2*q2(mu)*
     &    q2(nu)*pp_ep*em_q1*ul + 2*q2(mu)*q2(nu)*pp_ep*em_q2*ul + 2*
     &    q2(mu)*q2(nu)*ep_em*ul**2 )
      amp4 = amp4 + n_pp**2 * (  - 4*d_(mu,nu)*ep_q2*em_q1*ul**2 )
      amp4 = amp4 + n_q1*n_q2*qsq2 * (  - 2*ep(nu)*em(mu)*ul**2 )
      amp4 = amp4 + n_q1*n_q2 * (  - 2*em(mu)*q2(nu)*ep_q2*ul**2 )
      amp4 = amp4 + n_q2**2*qsq1 * ( 2*ep(nu)*em(mu)*ul**2 )
      amp4 = amp4 + n_q2**2 * (  - 4*pp(mu)*pp(nu)*pp_ep*pp_em*ul + 4*
     &    pp(mu)*pp(nu)*pp_ep*em_q1*ul + 4*pp(mu)*pp(nu)*pp_ep*em_q2*ul
     &     + 4*pp(mu)*pp(nu)*ep_em*ul**2 - 2*pp(mu)*ep(nu)*pp_em*ul**2
     &     + 2*pp(mu)*ep(nu)*pp_em*mw**2*ul + 2*pp(mu)*ep(nu)*em_q1*
     &    ul**2 - 2*pp(mu)*ep(nu)*em_q1*mw**2*ul - 2*pp(mu)*ep(nu)*
     &    em_q2*ul**2 - 2*pp(mu)*ep(nu)*em_q2*mw**2*ul + 4*pp(mu)*
     &    em(nu)*ep_q2*ul**2 + 2*pp(mu)*q2(nu)*pp_ep*pp_em*ul - 2*
     &    pp(mu)*q2(nu)*pp_ep*em_q1*ul - 2*pp(mu)*q2(nu)*pp_ep*em_q2*ul
     &     - 2*pp(mu)*q2(nu)*ep_em*ul**2 - 4*pp(nu)*ep(mu)*em_q1*ul**2
     &     - 2*pp(nu)*em(mu)*pp_ep*ul**2 + 2*pp(nu)*em(mu)*pp_ep*mw**2*
     &    ul + 4*pp(nu)*em(mu)*ep_q1*ul**2 + 2*pp(nu)*q1(mu)*pp_ep*
     &    pp_em*ul - 2*pp(nu)*q1(mu)*pp_ep*em_q1*ul - 2*pp(nu)*q1(mu)*
     &    pp_ep*em_q2*ul - 2*pp(nu)*q1(mu)*ep_em*ul**2 + 4*pp(nu)*
     &    q2(mu)*pp_ep*pp_em*ul - 4*pp(nu)*q2(mu)*pp_ep*em_q1*ul - 4*
     &    pp(nu)*q2(mu)*pp_ep*em_q2*ul - 4*pp(nu)*q2(mu)*ep_em*ul**2 + 
     &    2*ep(mu)*q2(nu)*em_q1*ul**2 )
      amp4 = amp4 + n_q2**2 * ( ep(nu)*em(mu)*ul**3 + 2*ep(nu)*em(mu)*
     &    tl*ul**2 - 2*ep(nu)*em(mu)*mw**2*ul**2 - ep(nu)*em(mu)*mw**4*
     &    ul + ep(nu)*q1(mu)*pp_em*ul**2 - ep(nu)*q1(mu)*pp_em*mw**2*ul
     &     + ep(nu)*q1(mu)*em_q1*ul**2 + ep(nu)*q1(mu)*em_q1*mw**2*ul
     &     + ep(nu)*q1(mu)*em_q2*ul**2 + ep(nu)*q1(mu)*em_q2*mw**2*ul
     &     + 2*ep(nu)*q2(mu)*pp_em*ul**2 - 2*ep(nu)*q2(mu)*pp_em*mw**2*
     &    ul + 2*ep(nu)*q2(mu)*em_q1*ul**2 + 2*ep(nu)*q2(mu)*em_q1*
     &    mw**2*ul + 2*ep(nu)*q2(mu)*em_q2*ul**2 + 2*ep(nu)*q2(mu)*
     &    em_q2*mw**2*ul + 4*em(mu)*q1(nu)*ep_q2*ul**2 + em(mu)*q2(nu)*
     &    pp_ep*ul**2 - em(mu)*q2(nu)*pp_ep*mw**2*ul - 2*em(mu)*q2(nu)*
     &    ep_q1*ul**2 - 2*em(nu)*q1(mu)*ep_q2*ul**2 - 4*em(nu)*q2(mu)*
     &    ep_q2*ul**2 - q1(mu)*q2(nu)*pp_ep*pp_em*ul + q1(mu)*q2(nu)*
     &    pp_ep*em_q1*ul + q1(mu)*q2(nu)*pp_ep*em_q2*ul + q1(mu)*q2(nu)
     &    *ep_em*ul**2 - 2*q2(mu)*q2(nu)*pp_ep*pp_em*ul + 2*q2(mu)*
     &    q2(nu)*pp_ep*em_q1*ul + 2*q2(mu)*q2(nu)*pp_ep*em_q2*ul + 2*
     &    q2(mu)*q2(nu)*ep_em*ul**2 )
      amp4 = amp4 + n_q2**2 * (  - 4*d_(mu,nu)*ep_q2*em_q1*ul**2 )
      amp4 = amp4 + n(mu)*n_pp*qsq2 * (  - 2*ep(nu)*em_q1*ul**2 )
      amp4 = amp4 + n(mu)*n_pp * (  - 2*q2(nu)*ep_q2*em_q1*ul**2 )
      amp4 = amp4 + n(mu)*n_q2*qsq2 * ( 2*ep(nu)*em_q1*ul**2 )
      amp4 = amp4 + n(mu)*n_q2 * ( 2*q2(nu)*ep_q2*em_q1*ul**2 )
      amp4 = amp4 + n(nu)*n_pp*qsq1 * ( 2*em(mu)*ep_q2*ul**2 )
      amp4 = amp4 + n(nu)*n_pp * ( 2*q1(mu)*ep_q2*em_q1*ul**2 )
      amp4 = amp4 + n(nu)*n_q2*qsq1 * (  - 2*em(mu)*ep_q2*ul**2 )
      amp4 = amp4 + n(nu)*n_q2 * (  - 2*q1(mu)*ep_q2*em_q1*ul**2 )

      
      den4 =
     &  + n_pp*n_q2 * (  - 2*ul**3 + 2*mw**2*ul**2 )
      den4 = den4 + n_pp**2 * ( ul**3 - mw**2*ul**2 )
      den4 = den4 + n_q2**2 * ( ul**3 - mw**2*ul**2 )

      
      amp5 =
     &  + ep(mu)*em(nu)*mw**2 + ep(nu)*em(mu)*mw**2 + 2*d_(mu,nu)*pp_ep
     &    *pp_em - 2*d_(mu,nu)*pp_ep*em_q1 - 2*d_(mu,nu)*pp_ep*em_q2 - 
     &    2*d_(mu,nu)*ep_em*mw**2

      den5 =
     &  + mw**2
      

      
      ampaxial_fun = amp1/den1+amp2/den2+amp3/den3+amp4/den4
     &     +amp5/den5
      return
      end


      
