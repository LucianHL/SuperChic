      program apfelsf
      implicit none
      double precision alphas,alphasPDF
      integer order, argc, N
      character pdfname*100
      argc=iargc()
      if (argc.LT.1) THEN
        write(*,*)'Not enought arguments',argc
        stop 1
      end if
      CALL GETARG(1, pdfname)
      call SetMassScheme("ZM-VFNS")
      call SetMaxFlavourPDFs(5)    
      call SetQLimits(1d0,1d5)
      call SetNumberOfGrids(3)
      call SetGridParameters(1,100,3,1d-9)
      call SetGridParameters(2,100,3,1d-4)
      call SetGridParameters(3,100,3,1d-2)

      call InitPDFsetByName(trim(pdfname)) !AV
      call InitPDF(0) !AV
      call numberPDF(N)
      call GetOrderAs(order)
      call SetPerturbativeOrder(order) 

      
      call SetPDFSet(trim(pdfname))
      alphas=alphasPDF(91.1876d0)
      call SetAlphaQCDRef(alphas,91.1876d0) 

      call SetReplica(0)
      
      call InitializeAPFEL_DIS
      call CacheStructureFunctionsAPFEL(1d0)
       
      call EXECUTE_COMMAND_LINE('mkdir -p SF_'//trim(pdfname)) ! AV
      call LHAPDFgridStructureFunctions(N,1d0,'SF_'//trim(pdfname))
      end