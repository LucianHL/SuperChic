
enable_language(CXX)
find_package(HepMC3 3.1.0 REQUIRED)
message(STATUS "SuperChic test: HEPMC3_VERSION=${HEPMC3_VERSION} HEPMC3_LIBRARIES=${HEPMC3_LIBRARIES} HEPMC3_INCLUDE_DIRS=${HEPMC3_INCLUDE_DIRS}")
add_executable(validator ${CMAKE_CURRENT_SOURCE_DIR}/validator.cxx)
target_link_libraries(validator PRIVATE ${HEPMC3_LIBRARIES})
target_include_directories(validator PRIVATE ${HEPMC3_INCLUDE_DIRS})
target_compile_features(validator PRIVATE cxx_std_11)

find_package(Pythia8 8.3 REQUIRED)
add_executable(shower ${CMAKE_CURRENT_SOURCE_DIR}/shower.cxx)
target_link_libraries(shower PRIVATE ${HEPMC3_LIBRARIES} Pythia8::Pythia8 )
target_include_directories(shower PRIVATE ${HEPMC3_INCLUDE_DIRS})
target_compile_features(shower PRIVATE cxx_std_11)
message(STATUS "SuperChic test: PYTHIA8_VERSION=${PYTHIA8_VERSION} PYTHIA8_LIBRARIES=${PYTHIA8_LIBRARIES} PYTHIA8_INCLUDE_DIRS=${PYTHIA8_INCLUDE_DIRS}")


set(STANDARDENVIRONMENT "SUPERCHIC_SOURCE_PATH=${PROJECT_SOURCE_DIR};PYTHIA8_XMLDOC_DIR=${PYTHIA8_XMLDOC_DIR}")

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
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/inputs)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/outputs)
  add_test(NAME INIT_${NN} COMMAND sh -c "${PROJECT_BINARY_DIR}/init < ${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_${NN}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  if(SUPERCHIC_DOWNLOAD_PDFS)
    set_tests_properties( INIT_${NN} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS UNZIPPDFS)
  else()
    set_tests_properties( INIT_${NN} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}")
  endif()
endmacro()

macro(sctesthepmc NN MM)
  add_test(NAME SUPERCHIC_${MM} COMMAND sh -c "${PROJECT_BINARY_DIR}/superchic < ${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_${MM}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  set_property(TEST SUPERCHIC_${MM} PROPERTY  DEPENDS INIT_${NN}) 
  if(SUPERCHIC_DOWNLOAD_PDFS)
    set_tests_properties( SUPERCHIC_${MM} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS UNZIPPDFS)
  else()
    set_tests_properties( SUPERCHIC_${MM} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}")  
  endif()
  add_test(NAME VALIDATE_${MM} COMMAND validator ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat)
  set_tests_properties( VALIDATE_${MM} PROPERTIES DEPENDS validator)
  set_property(TEST VALIDATE_${MM} PROPERTY  DEPENDS SUPERCHIC_${MM}) 
endmacro()




macro(sctestlhe NN MM)
  add_test(NAME SUPERCHIC_${MM} COMMAND sh -c "${PROJECT_BINARY_DIR}/superchic < ${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_${MM}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  set_property(TEST SUPERCHIC_${MM} PROPERTY  DEPENDS INIT_${NN}) 
  if(SUPERCHIC_DOWNLOAD_PDFS)
    set_tests_properties( SUPERCHIC_${MM} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS UNZIPPDFS)
  else()
    set_tests_properties( SUPERCHIC_${MM} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}")  
  endif()
  add_test(NAME SHOWER_${MM} COMMAND shower ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.hepmc)
  set_tests_properties( SHOWER_${MM} PROPERTIES DEPENDS shower)
  set_property(TEST SHOWER_${MM} PROPERTY  DEPENDS SUPERCHIC_${MM}) 
  add_test(NAME VALIDATE_${MM} COMMAND validator ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.hepmc)
  set_tests_properties( VALIDATE_${MM} PROPERTIES DEPENDS validator)
  set_property(TEST VALIDATE_${MM} PROPERTY  DEPENDS SHOWER_${MM}) 
endmacro()



inittest(1)
sctesthepmc(1 1)
sctestlhe(1 2)

inittest(3)
sctesthepmc(3 3)

inittest(9)
sctesthepmc(9 9)
