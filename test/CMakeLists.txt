
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


set(STANDARDENVIRONMENT "SUPERCHIC_DATA_PATH=${PROJECT_SOURCE_DIR}/data;PYTHIA8_XMLDOC_DIR=${PYTHIA8_XMLDOC_DIR}")

if (SUPERCHIC_DOWNLOAD_PDFS)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/PDFS)
  file(DOWNLOAD https://superchic.hepforge.org/SF_MSHT20qed_nnlo.tar.gz ${CMAKE_CURRENT_BINARY_DIR}/PDFS/SF_MSHT20qed_nnlo.tar.gz)
  file(DOWNLOAD https://lhapdfsets.web.cern.ch/current/MSHT20qed_nnlo.tar.gz ${CMAKE_CURRENT_BINARY_DIR}/PDFS/MSHT20qed_nnlo.tar.gz)
  add_custom_target( UNZIPPDFS ALL)
  add_custom_command(TARGET UNZIPPDFS PRE_BUILD 
         COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/PDFS/SF_MSHT20qed_nnlo
         COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/PDFS/MSHT20qed_nnlo
         COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_BINARY_DIR}/PDFS/SF_MSHT20qed_nnlo.tar.gz
         COMMAND ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_BINARY_DIR}/PDFS/MSHT20qed_nnlo.tar.gz
         WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/PDFS
         DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/PDFS/SF_MSHT20qed_nnlo.tar.gz  ${CMAKE_CURRENT_BINARY_DIR}/PDFS/MSHT20qed_nnlo.tar.gz
         COMMENT "Unpacking PDFs"
         VERBATIM)
  LIST(APPEND STANDARDENVIRONMENT "LHAPDF_DATA_PATH=${CMAKE_CURRENT_BINARY_DIR}/PDFS")
endif()
message(STATUS "SuperChic test: STANDARDENVIRONMENT=${STANDARDENVIRONMENT}")


macro(inittest NN)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  add_test(NAME INIT_${NN} COMMAND sh -c "${PROJECT_BINARY_DIR}/init < ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_${NN}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  if(SUPERCHIC_DOWNLOAD_PDFS)
    set_tests_properties( INIT_${NN} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS UNZIPPDFS)
  else()
    set_tests_properties( INIT_${NN} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}")
  endif()
endmacro()

macro(sctesthepmc NN MM)
  add_test(NAME SUPERCHIC_${NN}_${MM} COMMAND sh -c "${PROJECT_BINARY_DIR}/superchic < ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_${MM}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  if(SUPERCHIC_DOWNLOAD_PDFS)
    set_tests_properties( SUPERCHIC_${NN}_${MM} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS "UNZIPPDFS;INIT_${NN}")
  else()
    set_tests_properties( SUPERCHIC_${NN}_${MM} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS "INIT_${NN}")  
  endif()
  add_test(NAME VALIDATE_${NN}_${MM} COMMAND validator ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat)
  set_tests_properties( VALIDATE_${NN}_${MM} PROPERTIES DEPENDS "validator;SUPERCHIC_${NN}_${MM}")
endmacro()




macro(sctestlhe NN MM)
  add_test(NAME SUPERCHIC_${NN}_${MM} COMMAND sh -c "${PROJECT_BINARY_DIR}/superchic < ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_${MM}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  if(SUPERCHIC_DOWNLOAD_PDFS)
    set_tests_properties( SUPERCHIC_${NN}_${MM} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS "UNZIPPDFS;INIT_${NN}")
  else()
    set_tests_properties( SUPERCHIC_${NN}_${MM} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS "INIT_${NN}")  
  endif()
  add_test(NAME SHOWER_${NN}_${MM} COMMAND shower ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.hepmc)
  set_tests_properties( SHOWER_${NN}_${MM} PROPERTIES DEPENDS "shower;SUPERCHIC_${NN}_${MM}")
  add_test(NAME VALIDATE_${NN}_${MM} COMMAND validator ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.hepmc)
  set_tests_properties( VALIDATE_${NN}_${MM} PROPERTIES DEPENDS "validator;SHOWER_${NN}_${MM}")
endmacro()


set(diffs  el sd dd)
set(diffs  el  dd) 
#sd does not work so far
set(fmts  hepevt lhe hepmc)
set(procs  1 2 3 4 5 6 7 8 9 10 11)
set(procs 56)

foreach ( diff ${diffs} )
foreach ( process ${procs} )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_PROCESS ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_PROCESS_${process}_${diff}) # for init
foreach ( format ${fmts} )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_PROCESS ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_PROCESS_${process}_${diff}_${format})
endforeach()
endforeach()
endforeach()


foreach ( diff ${diffs} )
foreach ( process ${procs} )
inittest(PROCESS_${process}_${diff})
sctestlhe(PROCESS_${process}_${diff} PROCESS_${process}_${diff}_lhe)
sctesthepmc(PROCESS_${process}_${diff} PROCESS_${process}_${diff}_hepmc)
# FIXME sctesthepmc(PROCESS_${process}_${diff} PROCESS_${process}_${diff}_hepevt)
endforeach()
endforeach()
