      subroutine MGcross(pf,out)
      implicit double precision(a-y)
      double precision pq1(4),pq2(4),q1(4),q2(4)
      double precision garr(-6:6),phot
      REAL*8 P(0:3,6)
      integer i,pf
        
      include 'mom.f'
      include 'vars.f'
      include 'pi.f'
      include 'norm.f'
      include 'wwpars.f'
      include 'xb.f'
      include 'pol.f'

c      print*,higgs,sfonly,addnsf
      
      pol=pf

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

      
      xi1=-qsq1/rts/(q1(4)-q1(3))
      xi2=-qsq2/rts/(q2(4)+q2(3))
      
      do i=1,4
         pq1(i)=0d0
         pq2(i)=0d0
      enddo
      
      pq1(4)=rts/2d0*xi1
      pq1(3)=pq1(4)
      pq2(4)=rts/2d0*xi2
      pq2(3)=-pq2(4)
     
      do i=1,3
         p(i,2)=pq1(i)
         p(i,1)=pq2(i)
         p(i,4)=pq1(i)-q1(i)
         p(i,3)=pq2(i)-q2(i)
         p(i,5)=q(i,6)
         p(i,6)=q(i,7)
      enddo
      
      p(0,2)=pq1(4)
      p(0,1)=pq2(4)
      p(0,4)=pq1(4)-q1(4)
      p(0,3)=pq2(4)-q2(4)
      p(0,5)=q(4,6)
      p(0,6)=q(4,7)

cccccccccccc
      
      CALL SMATRIX_uc(P,MATELEM)                                                                                                                                 
      out_uu=matelem                                                                            
      out_uu=out_uu/(4d0*pi)**2   
         
      CALL SMATRIX_ds(P,MATELEM)
      
      out_dd=matelem      
      out_dd=out_dd/(4d0*pi)**2 ! normalize so consistent with me
            
      CALL SMATRIX_du(P,MATELEM)

      out_du=matelem
      out_du=out_du/(4d0*pi)**2 ! normalize so consistent with me

      CALL SMATRIX_ubcb(P,MATELEM)
      
      out_ubub=matelem
      out_ubub=out_ubub/(4d0*pi)**2 ! normalize so consistent with me

      CALL SMATRIX_dbsb(P,MATELEM)
     
      out_dbdb=matelem
      out_dbdb=out_dbdb/(4d0*pi)**2 ! normalize so consistent with me

      CALL SMATRIX_dbub(P,MATELEM)

      out_dbub=matelem
      out_dbub=out_dbub/(4d0*pi)**2 ! normalize so consistent with me
      
      CALL SMATRIX_udb(P,MATELEM)

      out_dbu=matelem
      out_dbu=out_dbu/(4d0*pi)**2 ! normalize so consistent with me

c$$$cccccc swap quarks

      do i=1,3
         p(i,1)=pq1(i)
         p(i,2)=pq2(i)
         p(i,3)=pq1(i)-q1(i)
         p(i,4)=pq2(i)-q2(i)
      enddo
      
      p(0,1)=pq1(4)
      p(0,2)=pq2(4)
      p(0,3)=pq1(4)-q1(4)
      p(0,4)=pq2(4)-q2(4)
      
      CALL SMATRIX_du(P,MATELEM)

      out_ud=matelem
      out_ud=out_ud/(4d0*pi)**2 ! normalize so consistent with me
      
      CALL SMATRIX_dbub(P,MATELEM)

      out_ubdb=matelem
      out_ubdb=out_ubdb/(4d0*pi)**2 ! normalize so consistent with me

      CALL SMATRIX_udb(P,MATELEM)

      out_udb=matelem
      out_udb=out_udb/(4d0*pi)**2 ! normalize so consistent with me

ccccccdifferent ordering in theses cases
      
      do i=1,3
         p(i,1)=pq1(i)
         p(i,2)=pq2(i)
         p(i,4)=pq1(i)-q1(i)
         p(i,3)=pq2(i)-q2(i)
      enddo
      
      p(0,1)=pq1(4)
      p(0,2)=pq2(4)
      p(0,4)=pq1(4)-q1(4)
      p(0,3)=pq2(4)-q2(4)


      call SMATRIX_cub(P,MATELEM)

      out_uub=matelem
      out_uub=out_uub/(4d0*pi)**2 ! normalize so consistent with me

      CALL SMATRIX_sdb(P,MATELEM)   

      out_ddb=matelem
      out_ddb=out_ddb/(4d0*pi)**2 ! normalize so consistent with me 
      
      CALL SMATRIX_dub(P,MATELEM)
      
      out_dub=matelem
      out_dub=out_dub/(4d0*pi)**2 ! normalize so consistent with me
      
      
ccccccccccc  swap momenta
      
      do i=1,3
         p(i,2)=pq1(i)
         p(i,1)=pq2(i)
         p(i,3)=pq1(i)-q1(i)
         p(i,4)=pq2(i)-q2(i)
      enddo
      
      p(0,2)=pq1(4)
      p(0,1)=pq2(4)
      p(0,3)=pq1(4)-q1(4)
      p(0,4)=pq2(4)-q2(4)

      call SMATRIX_cub(P,MATELEM)

      out_ubu=matelem
      out_ubu=out_ubu/(4d0*pi)**2 ! normalize so consistent with me
            
      CALL SMATRIX_sdb(P,MATELEM)   

      out_dbd=matelem
      out_dbd=out_dbd/(4d0*pi)**2 ! normalize so consistent with me 
      
      CALL SMATRIX_dub(P,MATELEM)

      out_ubd=matelem
      out_ubd=out_ubd/(4d0*pi)**2 ! normalize so consistent with me
      
ccccccccccc

 444  call evolvePDFphoton(xi1,dsqrt(qsq1),garr,phot)
      
      f2_1u=(garr(2)+garr(4))  ! note charge weighting included in SMATRIX
      f2_1d=(garr(1)+garr(3)+garr(5))
      f2_1ub=(garr(-2)+garr(-4))
      f2_1db=(garr(-1)+garr(-3)+garr(-5))
      
      call evolvePDFphoton(xi2,dsqrt(qsq2),garr,phot)
      
      f2_2u=(garr(2)+garr(4))
      f2_2d=(garr(1)+garr(3)+garr(5))
      f2_2ub=(garr(-2)+garr(-4))
      f2_2db=(garr(-1)+garr(-3)+garr(-5))
      
cccccccccccc      

      out_uu=out_uu*f2_1u*f2_2u
      out_dd=out_dd*f2_1d*f2_2d
      out_du=out_du*f2_1d*f2_2u
      out_ud=out_ud*f2_2d*f2_1u
      out_ubub=out_ubub*f2_1ub*f2_2ub
      out_dbdb=out_dbdb*f2_1db*f2_2db
      out_dbub=out_dbub*f2_1db*f2_2ub
      out_ubdb=out_ubdb*f2_2db*f2_1ub
      out_uub=out_uub*f2_1u*f2_2ub
      out_ubu=out_ubu*f2_1ub*f2_2u
      out_ddb=out_ddb*f2_1d*f2_2db
      out_dbd=out_dbd*f2_1db*f2_2d
      out_dub=out_dub*f2_1d*f2_2ub
      out_ubd=out_ubd*f2_1ub*f2_2d
      out_udb=out_udb*f2_1u*f2_2db
      out_dbu=out_dbu*f2_1db*f2_2u


ccccc
      
      out=out_uu+out_dd+out_du+out_ud+out_ubub+out_dbdb+out_ubdb+out_uub
     &     +out_ubu+out_ddb+out_dbd+out_dub+out_ubd+out_udb+out_dbu
     &     +out_dbub

      
      out=out*conv/xi1/xi2

cccccccc

      call F1F2ev(xi1,qsq1,f1lo,f2lo)      
      call F1F2ap(xi1,qsq1,f2,fl)
      f1=f2-fl
      f1=f1/2d0/xi1

c      fk1=f1/f1lo
      fk1=f2/f2lo

      call F1F2ev(xi2,qsq2,f1lo,f2lo)      
      call F1F2ap(xi2,qsq2,f2,fl)
      f1=f2-fl
      f1=f1/2d0/xi2

c      fk2=f1/f1lo
      fk2=f2/f2lo

      out=out*fk1*fk2  ! reweight so that pure SF part is ~ NNLO

      return
      end
