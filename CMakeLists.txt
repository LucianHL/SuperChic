cmake_minimum_required(VERSION 3.15.0 FATAL_ERROR)
#----------------------------------------------------------------------------
project(SuperChic LANGUAGES  Fortran)
set(PROJECT_VERSION 5.0)
option(SUPERCHIC_ENABLE_TESTS        "Enables building of tests." OFF)
option(SUPERCHIC_ENABLE_FPES         "Enables FPEs." OFF)
option(SUPERCHIC_ENABLE_DOCS         "Enables compilation of manual." OFF)
option(SUPERCHIC_ENABLE_SHARED       "Enables compilation of shared library." OFF)
option(SUPERCHIC_INSTALL_LIBRARIES   "Enables installation of libraries." OFF)

include("GNUInstallDirs")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})
find_package(LHAPDF REQUIRED)
message(STATUS "SuperChic: LHAPDF_VERSION=${LHAPDF_VERSION} LHAPDF_LIBRARIES=${LHAPDF_LIBRARIES} LHAPDF_INCLUDE_DIRS=${LHAPDF_INCLUDE_DIRS}")
find_package(APFEL REQUIRED)
message(STATUS "SuperChic: APFEL_VERSION=${APFEL_VERSION} APFEL_LIBRARIES=${APFEL_LIBRARIES} APFEL_INCLUDE_DIRS=${APFEL_INCLUDE_DIRS}")
#GNU
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-std=legacy;-fbacktrace;-fcheck=all;-fno-automatic;-fimplicit-none>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wall;-Wextra;-Wno-unused-dummy-argument;-Wno-compare-reals;-Wno-do-subscript>")
#Classical Intel
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-noauto;-traceback;-check>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-warn;-diag-disable=10448;-diag-disable 8889;-diag-disable 6717;-diag-disable 8291;-diag-disable 7712;-diag-disable 5194>")
#Intel LLVM
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-noauto;-traceback>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-warn;-diag-disable=10448;-diag-disable 8889;-diag-disable 6717;-diag-disable 8291;-diag-disable 7712;-diag-disable 5194>")
#NVHPC
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Msave;-Mbounds;-traceback>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Minform=inform>")

#Tier 2 support compilers
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,PGI>:-Msave>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,XL>:-qsave>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,LLVMFlang>:-Msave>")
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,Flang>:-Msave>")
#Tier 3 support compilers (will not compile,probably)
add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-dcfuns;-dusty>")
if (SUPERCHIC_ENABLE_FPES)
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-ffpe-trap=zero,overflow,underflow,invalid,denormal>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-fpe0>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-fpe0>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-Ktrap=divz,denorm,fp,inexact,inv,ovf,unf>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,Flang>:-ffp-exception-behavior=strict>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,XL>:-qflttrap=overflow:underflow:invalid:zerodivide:enable>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,NAG>:-ieee=stop>")
  add_compile_options("$<$<COMPILE_LANG_AND_ID:Fortran,ARM>:-ffp-contract=on>")
endif()

add_definitions("-DDATA_PATH=\'${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_DATADIR}/SuperChic\'")

message(STATUS "SuperChic: CMAKE_Fortran_COMPILER_ID ${CMAKE_Fortran_COMPILER_ID} " )
message(STATUS "SuperChic: CMAKE_Fortran_FLAGS ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE}}  " )
message(STATUS "SuperChic: SUPERCHIC_ENABLE_FPES=${SUPERCHIC_ENABLE_FPES}" )
message(STATUS "SuperChic: SUPERCHIC_ENABLE_TESTS=${SUPERCHIC_ENABLE_TESTS}" )
message(STATUS "SuperChic: SUPERCHIC_ENABLE_DOCS=${SUPERCHIC_ENABLE_DOCS}" )

set(Dissf CLAS.f gd-fit-11.f CB.f Elastic.f F1F2.f F2ev.f rfit.f 
          alphaEM.f qinit.f apfelinit.f )
list(TRANSFORM Dissf PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/diss/")

set(LbyLf B0F2M.f D0404M.f xspenz.f eett_aux.f C01_gen.f )
list(TRANSFORM LbyLf PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/LbyL/")

set(Mesf calcmes.f mesint.f wfinit.f wfoctet.f wfsinglet.f )
list(TRANSFORM Mesf PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/mes/")

set(PDFsfLHA alphas.f inpdf.f )
list(TRANSFORM PDFsfLHA PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/PDFs/")

set(PDFsfUSER alphasuser.f inpdfuser.f )
list(TRANSFORM PDFsfUSER PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/PDFsfUSER/")

set(Intf vegas.f rann.f )
list(TRANSFORM Intf PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/int/")

set(Phasef 2body.f 2bodyw.f 2jetps.f 2jetpsm.f 3jetps.f boost.f 
           chic0decay3.f chic1decay3.f chic1decay2s.f chic1decay2f.f 
           chic2decay3.f chic2decay2s.f chic2decay2f.f jpsidecayphot.f 
           genpol1.f genpol1rf.f genpol2.f rambo.f 6body.f 6bodyinit.f 
           4body.f 4bodyinit.f 3body.f 3bodyinit.f 2bodyinit.f wwcorr.f 
           jpsidecay.f rhodecay.f chidecay.f monow.f alpdecay.f 
           pAboost.f pAinit.f AAinit.f genpol_axial.f )
list(TRANSFORM Phasef PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/phase/")

set(Subampsf chi0.f chi1.f chi2.f etaq.f higgs.f higgsinit.f pipi.f 
             qqjets.f diphoton.f etaeta.f etapetap.f etaetap.f eta.f 
             gggjets.f pipixy.f rhorho.f djpsi.f djpsip.f djpsipp.f 
             ggjets.f qqgjets.f rhorhoxy.f wwpol.f llpol.f mhv.f 
             lightlightpol.f higgsgam.f higgsgaminit.f alp.f mmpol.f 
             monop.f lloff.f spinors.f wwoff.f wwoff_amp.f )
list(TRANSFORM Subampsf PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/subamps/")

set(Survf initparsr.f formfac.f formfacphot.f formfacgam.f seik.f 
          seikphot.f seikgam.f screeningint.f readscreen.f 
          formfacgamel.f formfacgamion.f formfacgamionp.f tpint.f 
          seikgamion.f screeningionint.f seikion.f betaionex.f 
          #betaion.f
          seikphotionp.f formfacphotionp.f formfacgamoff.f 
          formfacgamoff_ion.f formfacgamoff_ion_surv.f 
          formfacgamoff_surv.f s2_diss_int.f s2_diss_init.f )
list(TRANSFORM Survf PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/surv/")

set(Userf cuts.f histo.f )
list(TRANSFORM Userf PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/user/")

set(Mainf bare.f header.f main.f process.f wtgen.f wtgengam.f 
          header_out.f )
list(TRANSFORM Mainf PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/main/")

set(sPDFsf calchg.f hpdfint.f calcsud.f sudint.f sPDF.f )
list(TRANSFORM sPDFsf PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/sPDFs/")

set(Unwf unweight.f headerlhe.f  unweight_q.f unwprint_q.f 
         unwprint.f unweight.f )
list(TRANSFORM Unwf PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/unw/")
         
set(Varf mu.f nf.f string.f varfuncs.f gaminit.f gauss.f )
list(TRANSFORM Varf PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/var/")

set(InitfLHA PDFs/alphas.f init/initsud.f var/nf.f 
             var/string.f 
             sPDFs/hg.f init/inithg.f init/initpars.f 
             init/calcop.f init/calcscreen.f init/opacityint.f surv/screeningint.f 
             init/screening.f init/opacity.f PDFs/PDF.f PDFs/PDFlha.f sPDFs/Sudakov.f PDFs/inpdf.f 
             init/dd.f diss/Elastic.f init/sdcoh.f int/rann.f surv/readscreen.f init/sdincoh.f 
             diss/F1F2.f diss/F2ev.f diss/rfit.f diss/gd-fit-11.f diss/CLAS.f diss/apfelinit.f 
             init/dd_test.f )
list(TRANSFORM InitfLHA PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/")

set(InitfION rho.f ionpars.f rhonorm.f rhoxy.f rhoxycalc.f tpcalc.f tp.f
             rhoxyint.f betaion.f opacpcalc.f opacp.f opacpint.f opacpb.f
             opacpbcalc.f opacpbp.f opacpbpcalc.f opacpbint.f 
             opacpbpint.f screencalc.f screen.f 
             #string.f 
             s2qcdion.f 
             s2qcdionp.f ioninit.f tpqcd.f tpqcdcalc.f tpqcdint.f 
             gdrin.f gdr_int.f sudgam.f screen_01.f )
list(TRANSFORM InitfION PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/inition/")

set(InitfUSER alphasuser.f init.f initsud.f nf.f 
              #string.f 
              hg.f inithg.f 
              initpars.f calcop.f calcscreen.f opacityint.f 
              screeningint.f screening.f opacity.f PDF.f PDFuser.f 
              Sudakov.f inpdfuser.f dd.f Elastic.f sdcoh.f rann.f 
              readscreen.f sdincoh.f F1F2.f F2ev.f rfit.f gd-fit-11.f 
              CLAS.f apfelinit.f )
list(TRANSFORM InitfUSER PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/init/")

set(Helas FFV2_0.f FFV2_1.f FFV2_2.f FFV2_3.f FFV1_0.f FFV1_1.f FFV1_2.f
          FFV1P0_3.f FFV3_0.f FFV3_1.f FFV3_2.f FFV3_3.f VVS1_0.f 
          VVS1_3.f VVS1_3.f VVVV5_0.f FFV5_0.f FFV5_1.f FFV5_2.f FFV5_3.f
          VVVV2_0.f VVV1_0.f VVV1_1.f VVV1P0_1.f VVV1_2.f VVV1_3.f 
          aloha_functions.f VVVV5_0.f VVVV2P0_1.f VVVV5_4.f )
list(TRANSFORM Helas PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/DHELAS/")

set(Model rw_para.f lha_read.f couplings.f couplings1.f couplings2.f )
list(TRANSFORM Model PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/MODEL/")

set(MG wwpol_match.f is_born.f matrix_uc.f matrix_ds.f matrix_ubcb.f 
       matrix_dbsb.f matrix_dbub.f matrix_udb.f matrix_cub.f matrix_sdb.f 
       matrix_dub.f matrix_du.f mgcross.f cs_SF_SD.f matrix_au.f 
       matrix_aub.f matrix_ad.f matrix_adb.f )
list(TRANSFORM MG PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/src/MG/")

set( sCODELHA ${Mainf} ${Mesf} ${EW} ${PDFsfLHA} ${Intf} ${Phasef} 
              ${Subampsf} ${Survf} ${Userf} ${sPDFsf} ${Unwf} ${Varf} 
              ${LbyLf} ${InitfION} ${Dissf} ${Helas} ${Model} ${MG})

set(scincludes ${CMAKE_CURRENT_SOURCE_DIR}/src/main ${CMAKE_CURRENT_SOURCE_DIR}/src/inc)
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(sccompile_options )
else()
  set(sccompile_options "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-mcmodel=large>" 
                        "$<$<COMPILE_LANG_AND_ID:Fortran,Intel>:-mcmodel=large>"
                        "$<$<COMPILE_LANG_AND_ID:Fortran,IntelLLVM>:-mcmodel=large>" 
                        "$<$<COMPILE_LANG_AND_ID:Fortran,NVHPC>:-mcmodel=medium>"
                        "$<$<COMPILE_LANG_AND_ID:Fortran,PGI>:-mcmodel=medium>")
endif()
set_source_files_properties( ${CMAKE_CURRENT_SOURCE_DIR}/src/diss/Elastic.f ${CMAKE_CURRENT_SOURCE_DIR}/src/inition/gdrin.f 
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/user/histo.f 
                             PROPERTIES Fortran_PREPROCESS ON Fortran_FORMAT FREE)


if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin" OR  ${CMAKE_VERSION} VERSION_LESS 3.18)
set_source_files_properties( ${CMAKE_CURRENT_SOURCE_DIR}/src/user/histo.f ${CMAKE_CURRENT_SOURCE_DIR}/src/diss/Elastic.f ${CMAKE_CURRENT_SOURCE_DIR}/src/inition/gdrin.f 
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/user/histo.f 
                             PROPERTIES COMPILE_FLAGS "-traditional -cpp" )  
endif()

if (SUPERCHIC_ENABLE_SHARED)
  add_library(superchicLib SHARED ${sCODELHA})
  target_include_directories(superchicLib PRIVATE  ${scincludes})
  target_link_libraries(superchicLib PRIVATE LHAPDF::LHAPDF APFEL::APFEL APFEL::APFELevol)
  target_compile_options(superchicLib PRIVATE ${sccompile_options})
  set_target_properties(superchicLib PROPERTIES 
    ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
  )
  if (SUPERCHIC_ENABLE_INSTALL_LIBRARIES)
    install(TARGETS superchicLib   DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libs)
  endif()
endif()

add_library(superchicLibstatic STATIC ${sCODELHA})
target_include_directories(superchicLibstatic PRIVATE  ${scincludes})
target_link_libraries(superchicLibstatic PRIVATE LHAPDF::LHAPDF APFEL::APFEL APFEL::APFELevol)
target_compile_options(superchicLibstatic PRIVATE ${sccompile_options})
set_target_properties(superchicLibstatic PROPERTIES  OUTPUT_NAME superchicLib-static
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}/$<0:>
  )
  
add_executable(superchic ${CMAKE_CURRENT_SOURCE_DIR}/src/main/superchic.f  ${InitfLHA} )
target_include_directories(superchic PRIVATE  ${scincludes})
target_link_libraries(superchic PRIVATE LHAPDF::LHAPDF APFEL::APFEL APFEL::APFELevol superchicLibstatic)
target_compile_options(superchic PRIVATE ${sccompile_options})
set_target_properties(superchic PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/$<0:>)

add_executable(init ${CMAKE_CURRENT_SOURCE_DIR}/src/init/init.f ${InitfLHA})
target_include_directories(init PRIVATE  ${scincludes})
target_link_libraries(init PRIVATE LHAPDF::LHAPDF APFEL::APFEL  APFEL::APFELevol )
target_compile_options(init PRIVATE ${sccompile_options})
set_target_properties(init PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}/$<0:>)


## MANUAL
if (SUPERCHIC_ENABLE_DOCS)
find_package(LATEX COMPONENTS PDFLATEX BIBTEX)
message(STATUS "SUPERCHIC: LATEX_FOUND=${LATEX_FOUND} PDFLATEX_COMPILER=${PDFLATEX_COMPILER} BIBTEX_COMPILER=${BIBTEX_COMPILER}")
if ( LATEX_FOUND AND LATEX_PDFLATEX_FOUND AND LATEX_BIBTEX_FOUND)
    file(COPY ${PROJECT_SOURCE_DIR}/doc/h-physrev.bst DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc/)
    file(COPY ${PROJECT_SOURCE_DIR}/doc/references.bib DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/doc/)
    configure_file(${PROJECT_SOURCE_DIR}/doc/superchic.tex ${PROJECT_BINARY_DIR}/doc/superchic${PROJECT_VERSION}.tex IMMEDIATE @ONLY)
    ADD_CUSTOM_TARGET(   latexcompileall
                       COMMAND ${PDFLATEX_COMPILER} superchic${PROJECT_VERSION}
                       COMMAND ${BIBTEX_COMPILER} superchic${PROJECT_VERSION}
                       COMMAND ${PDFLATEX_COMPILER} superchic${PROJECT_VERSION}
                       COMMAND ${PDFLATEX_COMPILER} superchic${PROJECT_VERSION}
                       COMMENT "Creating manual "
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc
                       DEPENDS ${PROJECT_BINARY_DIR}/doc/superchic${PROJECT_VERSION}.tex
                      )
    ADD_CUSTOM_TARGET(all-formats ALL) # Entry point of execution.
    ADD_DEPENDENCIES( all-formats latexcompileall)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/doc/superchic${PROJECT_VERSION}.pdf DESTINATION  ${CMAKE_INSTALL_DOCDIR} COMPONENT doc)
  endif()
endif()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/bin/input.DAT
              DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/bin/)

install(TARGETS superchic init DESTINATION ${CMAKE_INSTALL_BINDIR} COMPONENT bin)
if (SUPERCHIC_ENABLE_INSTALL_LIBRARIES)
 install(TARGETS superchicLibstatic   DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT libs)
endif()

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Cards/ DESTINATION  ${CMAKE_INSTALL_DOCDIR}/Cards COMPONENT doc)
if ( NOT "${CMAKE_CURRENT_SOURCE_DIR}" STREQUAL "${CMAKE_CURRENT_BINARY_DIR}")
   file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/SplinesWithVariableKnots.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Carlos_106_440.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/gampgamn.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Veyssiere_singleneut.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Muccifora.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Lepretre25_103.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Caldwell.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Lepretre_25_103.dat 
              DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/share/SuperChic/)
endif()

install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/SplinesWithVariableKnots.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Carlos_106_440.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/gampgamn.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Veyssiere_singleneut.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Muccifora.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Lepretre25_103.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Caldwell.dat
              ${CMAKE_CURRENT_SOURCE_DIR}/share/SuperChic/Lepretre_25_103.dat 
              DESTINATION ${CMAKE_INSTALL_DATADIR}/SuperChic COMPONENT data)

if (SUPERCHIC_ENABLE_TESTS)
  ENABLE_TESTING()
  add_subdirectory(test)
endif()
