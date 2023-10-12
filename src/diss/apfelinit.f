      subroutine apfelinit
      implicit double precision(a-y)
      integer order
      
      include 'pdfinf.f'
      include 'ewsf.f'

c      call initpdfsetm(2,'SF_MMHT2015qed_nnlo')
      call initpdfsetm(2,'SF_MSHT20qed_nnlo')
      call initpdfm(2,0)

      Sin2ThetaW=0.23126d0
      ve = - 0.5d0 + 2d0 * Sin2ThetaW
      ae = - 0.5d0

      eq(1) = - 1d0 / 3d0
      eq(2) = 2d0 / 3d0
      eq(3) = eq(1) ! - 1d0 / 3d0
      eq(4) = eq(2) ! 2d0 / 3d0
      eq(5) = eq(1) ! - 1d0 / 3d0
      eq(6) = eq(2)             ! 2d0 / 3d0

      eq2(1) = eq(1) * eq(1)
      eq2(2) = eq(2) * eq(2)
      eq2(3) = eq2(1) ! eq(3) * eq(3)
      eq2(4) = eq2(2) ! eq(4) * eq(4)
      eq2(5) = eq2(1) ! eq(5) * eq(5)
      eq2(6) = eq2(2) ! eq(6) * eq(6)
      
      vq(1) = - 0.5d0 + 2d0 / 3d0 * Sin2ThetaW
      vq(2) = + 0.5d0 - 4d0 / 3d0 * Sin2ThetaW
      vq(3) = vq(1) ! - 0.5d0 + 2d0 / 3d0 * Sin2ThetaW
      vq(4) = vq(2) ! + 0.5d0 - 4d0 / 3d0 * Sin2ThetaW
      vq(5) = vq(1) ! - 0.5d0 + 2d0 / 3d0 * Sin2ThetaW
      vq(6) = vq(2) ! + 0.5d0 - 4d0 / 3d0 * Sin2ThetaW

      aq(1) = - 0.5d0
      aq(2) = + 0.5d0
      aq(3) = aq(1) ! - 0.5d0
      aq(4) = aq(2) ! + 0.5d0
      aq(5) = aq(1) ! - 0.5d0
      aq(6) = aq(2)             ! + 0.5d0

      do i=1,6
         agam(i)=eq2(i) 
         azgam(i)=-2d0*eq(i)*vq(i)
         az(i)=vq(i)**2+aq(i)**2
         bzgam(i)=-2d0*eq(i)*aq(i)
         bz(i)=2d0*vq(i)*aq(i)
      enddo
      
c$$$      return
c$$$      
c$$$      call SetMassScheme("ZM-VFNS")
c$$$      call SetMaxFlavourPDFs(5)    
c$$$      call SetQLimits(1d0,1d5)
c$$$      call SetNumberOfGrids(3)
c$$$      call SetGridParameters(1,100,3,1d-9)
c$$$      call SetGridParameters(2,100,3,1d-4)
c$$$      call SetGridParameters(3,100,3,1d-2)
c$$$
c$$$      
c$$$      
c$$$      call GetOrderAs(order)
c$$$      call SetPerturbativeOrder(order)         
c$$$      call SetPDFSet(PDFname)
c$$$      alphas=alphasPDF(91.1876d0)
c$$$      call SetAlphaQCDRef(alphas,91.1876d0) 
c$$$
c$$$      call InitializeAPFEL_DIS
c$$$      call CacheStructureFunctionsAPFEL(1d0)
      
      return
      end

