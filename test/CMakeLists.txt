
enable_language(CXX)
find_package(HepMC3 3.1.0 REQUIRED)
message(STATUS "SuperChic test: HEPMC3_VERSION=${HEPMC3_VERSION} HEPMC3_LIBRARIES=${HEPMC3_LIBRARIES} HEPMC3_INCLUDE_DIR=${HEPMC3_INCLUDE_DIR}")
add_executable(validator ${CMAKE_CURRENT_SOURCE_DIR}/validator.cxx)
target_link_libraries(validator PRIVATE ${HEPMC3_LIBRARIES})
target_include_directories(validator PRIVATE ${HEPMC3_INCLUDE_DIR})
target_compile_features(validator PRIVATE cxx_std_11)

find_package(Pythia8 8.3 REQUIRED)
add_executable(shower ${CMAKE_CURRENT_SOURCE_DIR}/shower.cxx)
target_link_libraries(shower PRIVATE ${HEPMC3_LIBRARIES} Pythia8::Pythia8 )
target_include_directories(shower PRIVATE ${HEPMC3_INCLUDE_DIRS})
target_compile_features(shower PRIVATE cxx_std_11)
message(STATUS "SuperChic test: PYTHIA8_VERSION=${PYTHIA8_VERSION} PYTHIA8_LIBRARIES=${PYTHIA8_LIBRARIES} PYTHIA8_INCLUDE_DIRS=${PYTHIA8_INCLUDE_DIRS}")


set(STANDARDENVIRONMENT "SUPERCHIC_DATA_PATH=${PROJECT_SOURCE_DIR}/share/SuperChic;PYTHIA8_XMLDOC_DIR=${PYTHIA8_XMLDOC_DIR}")

if (SUPERCHIC_DOWNLOAD_PDFS)
  include(FetchContent)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/PDFS)
  FetchContent_Declare(SF_MSHT20qed_nnlo URL https://superchic.hepforge.org/SF_MSHT20qed_nnlo.tar.gz URL_HASH   MD5=956fe924fc72dd947769078b1b68d7d5)
  FetchContent_GetProperties(SF_MSHT20qed_nnlo)
  if (NOT SF_MSHT20qed_nnlo_POPULATED)
    FetchContent_Populate(SF_MSHT20qed_nnlo)
    file(CREATE_LINK ${sf_msht20qed_nnlo_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/PDFS/SF_MSHT20qed_nnlo SYMBOLIC)
  endif()

  FetchContent_Declare(MSHT20qed_nnlo URL https://lhapdfsets.web.cern.ch/current/MSHT20qed_nnlo.tar.gz URL_HASH   MD5=cce676b0bc21bb985c40fceac67e97ce)
  FetchContent_GetProperties(MSHT20qed_nnlo)
  if (NOT MSHT20qed_nnlo_POPULATED)
    FetchContent_Populate(MSHT20qed_nnlo)
    file(CREATE_LINK ${msht20qed_nnlo_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/PDFS/MSHT20qed_nnlo SYMBOLIC)
  endif()

  FetchContent_Declare(cteq6l1 URL https://lhapdfsets.web.cern.ch/current/cteq6l1.tar.gz URL_HASH   MD5=5611f1e9235151d9f67254aeb13bb65f)
  FetchContent_GetProperties(cteq6l1)
  if (NOT cteq6l1_POPULATED)
    FetchContent_Populate(cteq6l1)
    file(CREATE_LINK ${cteq6l1_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/PDFS/cteq6l1 SYMBOLIC)
  endif()
  LIST(APPEND STANDARDENVIRONMENT "LHAPDF_DATA_PATH=${CMAKE_CURRENT_BINARY_DIR}/PDFS")
endif()
message(STATUS "SuperChic test: STANDARDENVIRONMENT=${STANDARDENVIRONMENT}")


macro(inittest NN)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  add_test(NAME INIT_${NN} COMMAND sh -c "${PROJECT_BINARY_DIR}/bin/init < ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_${NN}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  set_tests_properties( INIT_${NN} PROPERTIES TIMEOUT 1000 ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS init)
endmacro()

macro(sctesthepmc NN MM)
  add_test(NAME SUPERCHIC_${NN}_${MM} COMMAND sh -c "${PROJECT_BINARY_DIR}/bin/superchic < ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_${MM}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  set_tests_properties( SUPERCHIC_${NN}_${MM} PROPERTIES TIMEOUT 500 ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS INIT_${NN})  
  add_test(NAME VALIDATE_${NN}_${MM} COMMAND validator ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat)
  set_tests_properties( VALIDATE_${NN}_${MM} PROPERTIES TIMEOUT 500 DEPENDS "validator;SUPERCHIC_${NN}_${MM}")
endmacro()

macro(sctestlhe NN MM df pr)
  add_test(NAME SUPERCHIC_${NN}_${MM} COMMAND sh -c "${PROJECT_BINARY_DIR}/bin/superchic < ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_${MM}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  set_tests_properties( SUPERCHIC_${NN}_${MM} PROPERTIES TIMEOUT 500 ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS INIT_${NN})
  add_test(NAME SHOWER_${NN}_${MM} COMMAND shower ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.hepmc ${df} ${pr} real)
  #add_test(NAME SHOWER_${NN}_${MM} COMMAND shower ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.hepmc ${df} ${pr} dummy)
  set_tests_properties( SHOWER_${NN}_${MM} PROPERTIES TIMEOUT 500 ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS "shower;SUPERCHIC_${NN}_${MM}")
  add_test(NAME VALIDATE_SHOWERED_HEPMC_${NN}_${MM} COMMAND validator ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.hepmc)
  set_tests_properties( VALIDATE_SHOWERED_HEPMC_${NN}_${MM} PROPERTIES TIMEOUT 500 DEPENDS "validator;SHOWER_${NN}_${MM}")
  add_test(NAME VALIDATE_LHEF_${NN}_${MM} COMMAND validator ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat)
  set_tests_properties( VALIDATE_LHEF_${NN}_${MM} PROPERTIES TIMEOUT 500 DEPENDS "validator;SUPERCHIC_${NN}_${MM}")

endmacro()


set(diffs  el sd dd)
set(diffs_dd   sd dd)
set(diffs_el   el)
#set(diffs  el  dd) 
#sd does not work so far
set(fmts  hepevt lhe hepmc)
if (SUPERCHIC_ENABLE_ALL_TESTS)
set(procs  1 2 3 4 5 6 7 8 9  
           10 11 12 13 14 15 16 17 18 19 
           20 21 22 23 24 25 26 27 28 29
           30 31 32 33 34 35 36 37 38 39
           40 41 42 43 44 45 46 47 48 49
           50 51 52 53 54 55 56 57 58 59
           60 61                   68 69
           70 71 72 73 74 75 76 77 78 79
                 82 83 84 
           )
set(procs_dd 54 55 56 57 58
             68 )
else()
set(procs 56 )
set(procs_dd 56 )
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_LHC13ppunw ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_LHC13ppunw) # for init
foreach ( diff ${diffs} )
foreach ( process ${procs} )
foreach ( format ${fmts} )
set (gencuts ".false.")
if (${process} EQUAL 3  OR  ${process} EQUAL 7 OR  ${process} EQUAL 8 OR  ${process} EQUAL 14 OR  ${process} EQUAL 15 OR  ${process} EQUAL 16 OR  ${process} EQUAL 70 OR  ${process} EQUAL 71 OR  ${process} EQUAL 72)
set (gencuts ".true.")
endif()
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_LHC13ppunw ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_LHC13ppunw_${process}_${diff}_${format})
endforeach()
endforeach()
endforeach()


inittest(LHC13ppunw)
foreach ( process ${procs} )
foreach ( diff ${diffs_el} )
sctestlhe(LHC13ppunw LHC13ppunw_${process}_${diff}_lhe ${diff} ${process})
sctesthepmc(LHC13ppunw LHC13ppunw_${process}_${diff}_hepmc)
# FIXME sctesthepmc(LHC13ppunw LHC13ppunw_${process}_${diff}_hepevt)
endforeach()
endforeach()
foreach ( process ${procs_dd} )
foreach ( diff ${diffs_dd} )
sctestlhe(LHC13ppunw LHC13ppunw_${process}_${diff}_lhe ${diff} ${process})
endforeach()
endforeach()

