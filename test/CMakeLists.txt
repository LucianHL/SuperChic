option(SUPERCHIC_ENABLE_PROFILE      "Enables callgrind profiling." OFF)
option(SUPERCHIC_ENABLE_ALL_TESTS    "Enables building of all tests." OFF)
option(SUPERCHIC_DOWNLOAD_PDFS       "Enables downloading of PDFs for tests." OFF)
option(SUPERCHIC_ENABLE_RIVET        "Enables testing with Rivet." OFF)
option(SUPERCHIC_ENABLE_PP           "Enables testing with pp beams." ON)
option(SUPERCHIC_ENABLE_EE           "Enables testing with ee beams." OFF)
option(SUPERCHIC_ENABLE_AA           "Enables testing with AA beams." OFF)
option(SUPERCHIC_ENABLE_PA           "Enables testing with pA beams." OFF)

message(STATUS "SuperChic test: SUPERCHIC_ENABLE_ALL_TESTS=${SUPERCHIC_ENABLE_ALL_TESTS}" )
message(STATUS "SuperChic test: SUPERCHIC_DOWNLOAD_PDFS=${SUPERCHIC_DOWNLOAD_PDFS}" )
message(STATUS "SuperChic test: SUPERCHIC_ENABLE_RIVET=${SUPERCHIC_ENABLE_RIVET}" )
message(STATUS "SuperChic test: SUPERCHIC_ENABLE_PROFILE=${SUPERCHIC_ENABLE_PROFILE}" )
message(STATUS "SuperChic test: SUPERCHIC_ENABLE_PP=${SUPERCHIC_ENABLE_PP}" )
message(STATUS "SuperChic test: SUPERCHIC_ENABLE_EE=${SUPERCHIC_ENABLE_EE}" )
message(STATUS "SuperChic test: SUPERCHIC_ENABLE_AA=${SUPERCHIC_ENABLE_AA}" )
message(STATUS "SuperChic test: SUPERCHIC_ENABLE_PA=${SUPERCHIC_ENABLE_PA}" )

enable_language(CXX)
find_package(HepMC3 3.1.0 REQUIRED)
message(STATUS "SuperChic test: HEPMC3_VERSION=${HEPMC3_VERSION} HEPMC3_LIBRARIES=${HEPMC3_LIBRARIES} HEPMC3_INCLUDE_DIR=${HEPMC3_INCLUDE_DIR}")
add_executable(validator ${CMAKE_CURRENT_SOURCE_DIR}/validator.cxx)
target_link_libraries(validator PRIVATE ${HEPMC3_LIBRARIES})
target_include_directories(validator PRIVATE ${HEPMC3_INCLUDE_DIR})
target_compile_features(validator PRIVATE cxx_std_11)


find_package(HepMC2)
find_package(Pythia8 8.3 REQUIRED)
add_executable(shower ${CMAKE_CURRENT_SOURCE_DIR}/shower.cxx)
target_link_libraries(shower PRIVATE HepMC2::HepMC2 Pythia8::Pythia8 )
target_compile_features(shower PRIVATE cxx_std_11)
target_compile_definitions(shower PRIVATE -DHEPMC2=1)
message(STATUS "SuperChic test: PYTHIA8_VERSION=${PYTHIA8_VERSION} PYTHIA8_LIBRARIES=${PYTHIA8_LIBRARIES} PYTHIA8_INCLUDE_DIRS=${PYTHIA8_INCLUDE_DIRS}")


set(STANDARDENVIRONMENT "SUPERCHIC_DATA_PATH=${PROJECT_SOURCE_DIR}/share/SuperChic;PYTHIA8_XMLDOC_DIR=${PYTHIA8_XMLDOC_DIR}")

if (SUPERCHIC_ENABLE_RIVET)
find_package (Python 3 COMPONENTS  Development Interpreter)
find_package(Rivet REQUIRED)
find_package(YODA REQUIRED)
message(STATUS "SuperChic: RIVET_VERSION=${RIVET_VERSION} RIVET_FOUND=${RIVET_FOUND}  RIVET_LIBRARIES=${RIVET_LIBRARIES} \
  RIVET_INCLUDE_DIRS=${RIVET_INCLUDE_DIRS} RIVET_DATA_PATH=${RIVET_DATA_PATH} RIVET_ANALYSIS_PATH=${RIVET_ANALYSIS_PATH} \
  Rivet_HEPMC2_FOUND=${Rivet_HEPMC2_FOUND}  Rivet_HEPMC3_FOUND=${Rivet_HEPMC3_FOUND} \
  RIVET_CONFIG_LIBS_STRING=${RIVET_CONFIG_LIBS_STRING}  RIVET_CONFIG_LIBS=${RIVET_CONFIG_LIBS}  RIVET_CONFIG_LIB_DIRS=${RIVET_CONFIG_LIB_DIRS} \
  RIVET_CONFIG_CPPFLAGS_STRING=${RIVET_CONFIG_CPPFLAGS_STRING}  RIVET_CONFIG_CPPFLAGS_DIRS=${RIVET_CONFIG_CPPFLAGS_DIRS} ")
message(STATUS "SuperChic: YODA_VERSION=${YODA_VERSION} YODA_LIBRARIES=${YODA_LIBRARIES} YODA_INCLUDE_DIRS=${YODA_INCLUDE_DIRS} YODA_CONFIG_PYTHONPATH_STRING=${YODA_CONFIG_PYTHONPATH_STRING}")
  list(APPEND STANDARDENVIRONMENT "RIVET_ANALYSIS_PATH=${RIVET_ANALYSIS_PATH};RIVET_DATA_PATH=${RIVET_DATA_PATH}")
  list(APPEND STANDARDENVIRONMENT "PYTHONPATH=${YODA_CONFIG_PYTHONPATH_STRING}:${RIVET_CONFIG_PYTHONPATH_STRING}:$ENV{PYTHONPATH}")
endif()

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
if (SUPERCHIC_ENABLE_PROFILE)
  add_test(NAME INIT_${NN} COMMAND sh -c "valgrind --tool=callgrind  --callgrind-out-file=CALLGRIND_INIT_${NN}.log ${PROJECT_BINARY_DIR}/bin/init < ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_${NN}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  set_tests_properties( INIT_${NN} PROPERTIES TIMEOUT 5000 ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS init)
else()
  add_test(NAME INIT_${NN} COMMAND sh -c "${PROJECT_BINARY_DIR}/bin/init < ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_${NN}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  set_tests_properties( INIT_${NN} PROPERTIES TIMEOUT 1500 ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS init)
endif()

endmacro()

macro(sctesthepmc NN MM)
  add_test(NAME SUPERCHIC_${NN}_${MM} COMMAND sh -c "${PROJECT_BINARY_DIR}/bin/superchic < ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_${MM}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  set_tests_properties( SUPERCHIC_${NN}_${MM} PROPERTIES TIMEOUT 500 ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS INIT_${NN})  
  add_test(NAME VALIDATE_${NN}_${MM} COMMAND validator ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat)
  set_tests_properties( VALIDATE_${NN}_${MM} PROPERTIES TIMEOUT 500 DEPENDS "validator;SUPERCHIC_${NN}_${MM}")
endmacro()

macro(sctestlhe NN MM df pr)
if (SUPERCHIC_ENABLE_PROFILE)
  add_test(NAME SUPERCHIC_${NN}_${MM} COMMAND sh -c "valgrind --tool=callgrind  --callgrind-out-file=CALLGRIND_SUPERCHIC_${NN}_${MM}.log ${PROJECT_BINARY_DIR}/bin/superchic < ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_${MM}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  set_tests_properties( SUPERCHIC_${NN}_${MM} PROPERTIES TIMEOUT 1500 ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS INIT_${NN})
else()
  add_test(NAME SUPERCHIC_${NN}_${MM} COMMAND sh -c "${PROJECT_BINARY_DIR}/bin/superchic < ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_${MM}" WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN})
  set_tests_properties( SUPERCHIC_${NN}_${MM} PROPERTIES TIMEOUT 500 ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS INIT_${NN})
endif()  

  add_test(NAME SHOWER_${NN}_${MM} COMMAND shower ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.hepmc ${df} ${pr} real)
  #add_test(NAME SHOWER_${NN}_${MM} COMMAND shower ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.hepmc ${df} ${pr} dummy)
  set_tests_properties( SHOWER_${NN}_${MM} PROPERTIES TIMEOUT 500 ENVIRONMENT "${STANDARDENVIRONMENT}" DEPENDS "shower;SUPERCHIC_${NN}_${MM}")
  add_test(NAME VALIDATE_SHOWERED_HEPMC_${NN}_${MM} COMMAND validator ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.hepmc)
  set_tests_properties( VALIDATE_SHOWERED_HEPMC_${NN}_${MM} PROPERTIES TIMEOUT 500 DEPENDS "validator;SHOWER_${NN}_${MM}")
  add_test(NAME VALIDATE_LHEF_${NN}_${MM} COMMAND validator ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/evrectest${MM}.dat)
  set_tests_properties( VALIDATE_LHEF_${NN}_${MM} PROPERTIES TIMEOUT 500 DEPENDS "validator;SUPERCHIC_${NN}_${MM}")

endmacro()

macro(sctestrivet NN  process DD)
if (SUPERCHIC_ENABLE_RIVET)
string(REPLACE " " ";" D ${DD})
#FCC350eeunw_${process}_${diff}_lhe ${diff}_ee ${process}
  set(fo "")
  foreach( diff ${D})
    set(MM ${NN}_${process}_${diff}_lhe)
    set(out ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.yoda)
    add_test(NAME RIVET_LHEF_${NN}_${MM} COMMAND ${Python_EXECUTABLE} ${RIVET_EXE} --ignore-beams -a MC_XS ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/evrecs/showered${MM}.hepmc -o  ${out})
    set_tests_properties( RIVET_LHEF_${NN}_${MM} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}" TIMEOUT 500 DEPENDS "SUPERCHIC_${NN}_${MM}")
    list(APPEND fo ${out})
  endforeach()
  add_test(NAME PLOT_LHEF_${NN}_${process} COMMAND rivet-mkhtml ${fo} -o ${CMAKE_CURRENT_BINARY_DIR}/RUN_${NN}/rivet/${process}  )
  foreach( diff ${D})
    set(MM ${NN}_${process}_${diff}_lhe)
    set_tests_properties( PLOT_LHEF_${NN}_${process} PROPERTIES ENVIRONMENT "${STANDARDENVIRONMENT}" TIMEOUT 500 DEPENDS "RIVET_LHEF_${NN}_${MM}")
  endforeach()
  endif()
endmacro()


##NUMBER OF EVENTS
if (NOT NEV)
  set(NEV 10)
endif()
##


set(fmts  hepevt lhe hepmc)
if (SUPERCHIC_ENABLE_PP)

foreach( i RANGE 1 100)
set( proc${i} "")
endforeach()

if (SUPERCHIC_ENABLE_ALL_TESTS)
foreach( i IN ITEMS 1 2 3 4 5 6 7 8 9  
           10 11 12 13 14 15 16 17 18 19 
           20 21 22 23 24 25 26 27 28 29
           30 31 32 33 34 35 36 37 38 39
           40 41 42 43 44 45 46 47 48 49
           50 51 52 53 54 55 56 57 58 59
           60 61                   68 69
           70 71 72 73 74 75 76 77 78 79
                 82 83 84 
           )
    list(APPEND proc${i} "el")
endforeach()           
foreach( i IN ITEMS 54 55 56 57 58
             68 )
    list(APPEND proc${i} "sda")
    list(APPEND proc${i} "dd")
endforeach()  
else()
  set(proc68 "el;sda;dd")
endif()



set(process 100)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_LHC13ppunw ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_LHC13ppunw) # for init
inittest(LHC13ppunw)
foreach( process RANGE 1 100)
  LIST(LENGTH proc${process} LISTCOUNT)
  if (${LISTCOUNT}  GREATER 0)
  foreach(diff ${proc${process}}) 
    foreach ( format ${fmts} )
      set (gencuts ".false.")
      if (${process} EQUAL 3  OR  ${process} EQUAL 7 OR  ${process} EQUAL 8 OR  ${process} EQUAL 14 OR  ${process} EQUAL 15 OR  ${process} EQUAL 16 OR  ${process} EQUAL 70 OR  ${process} EQUAL 71 OR  ${process} EQUAL 72)
        set (gencuts ".true.")
      endif()
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_LHC13ppunw ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_LHC13ppunw_${process}_${diff}_${format})
    endforeach()
    sctestlhe(LHC13ppunw LHC13ppunw_${process}_${diff}_lhe ${diff}_pp ${process})
    sctesthepmc(LHC13ppunw LHC13ppunw_${process}_${diff}_hepmc)
  endforeach()
  list(JOIN proc${process} " " X)
  sctestrivet(LHC13ppunw ${process} ${X})
  endif()
endforeach()
endif()

if (SUPERCHIC_ENABLE_EE)
foreach( i RANGE 1 100)
set( proc${i} "")
endforeach()
if (SUPERCHIC_ENABLE_ALL_TESTS)
foreach( i IN  ITEMS
                          55 56 57 58 59
           60 61                   68 69
           70 71 72 73 74 75 76 
                 82 83 84 
           )
    list(APPEND proc${i} "el")
endforeach()                
else()
set(proc68 "el")
endif()


set(process 100)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_FCC350eeunw ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_FCC350eeunw) # for init
inittest(FCC350eeunw)
foreach( process RANGE 1 100)
  LIST(LENGTH proc${process} LISTCOUNT)
  if (${LISTCOUNT}  GREATER 0)
  foreach(diff ${proc${process}}) 
    foreach ( format ${fmts} )
      set (gencuts ".false.")
      if (${process} EQUAL 3  OR  ${process} EQUAL 7 OR  ${process} EQUAL 8 OR  ${process} EQUAL 14 OR  ${process} EQUAL 15 OR  ${process} EQUAL 16 OR  ${process} EQUAL 70 OR  ${process} EQUAL 71 OR  ${process} EQUAL 72)
        set (gencuts ".true.")
      endif()
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_FCC350eeunw ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_FCC350eeunw_${process}_${diff}_${format})
    endforeach()
    sctestlhe(FCC350eeunw FCC350eeunw_${process}_${diff}_lhe ${diff}_ee ${process})
    sctesthepmc(FCC350eeunw FCC350eeunw_${process}_${diff}_hepmc)
  endforeach()
 # sctestrivet(FCC350eeunw ${process} ${proc${process}})
   list(JOIN proc${process} " " X)
  sctestrivet(FCC350eeunw ${process} ${X})
  endif()
endforeach()
endif()



if (SUPERCHIC_ENABLE_AA)

foreach( i RANGE 1 100)
set( proc${i} "")
endforeach()

if (SUPERCHIC_ENABLE_ALL_TESTS)
foreach( i IN ITEMS 1 2 3 4 5 6 7 8 9  
           10 11 12 13 14 15 16 17 18 19 
           20 21 22 23 24 25 26 27 28 29
           30 31 32 33 34 35 36 37 38 39
           40 41 42 43 44 45 46 47 48 49
           50 51 52 53 54 55 56 57 58 59
           60 61                   68 69
           70 71 72 73 74 75 76 77 78 79
                 82 83 84 
           )
    list(APPEND proc${i} "el")
endforeach()           
foreach( i IN ITEMS 54 55 56 57 58
             68 )
    list(APPEND proc${i} "sda")
    list(APPEND proc${i} "dd")
endforeach()  
else()
  set(proc68 "el;sda;dd")
endif()



set(process 100)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_LHC13AAunw ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_LHC13AAunw) # for init
inittest(LHC13AAunw)
foreach( process RANGE 1 100)
  LIST(LENGTH proc${process} LISTCOUNT)
  if (${LISTCOUNT}  GREATER 0)
  foreach(diff ${proc${process}}) 
    foreach ( format ${fmts} )
      set (gencuts ".false.")
      if (${process} EQUAL 3  OR  ${process} EQUAL 7 OR  ${process} EQUAL 8 OR  ${process} EQUAL 14 OR  ${process} EQUAL 15 OR  ${process} EQUAL 16 OR  ${process} EQUAL 70 OR  ${process} EQUAL 71 OR  ${process} EQUAL 72)
        set (gencuts ".true.")
      endif()
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_LHC13AAunw ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_LHC13AAunw_${process}_${diff}_${format})
    endforeach()
    sctestlhe(LHC13AAunw LHC13AAunw_${process}_${diff}_lhe ${diff}_AA ${process})
    sctesthepmc(LHC13AAunw LHC13AAunw_${process}_${diff}_hepmc)
  endforeach()
  list(JOIN proc${process} " " X)
  sctestrivet(LHC13AAunw ${process} ${X})
  endif()
endforeach()
endif()



if (SUPERCHIC_ENABLE_PA)

foreach( i RANGE 1 100)
set( proc${i} "")
endforeach()

if (SUPERCHIC_ENABLE_ALL_TESTS)
foreach( i IN ITEMS 1 2 3 4 5 6 7 8 9  
           10 11 12 13 14 15 16 17 18 19 
           20 21 22 23 24 25 26 27 28 29
           30 31 32 33 34 35 36 37 38 39
           40 41 42 43 44 45 46 47 48 49
           50 51 52 53 54 55 56 57 58 59
           60 61                   68 69
           70 71 72 73 74 75 76 77 78 79
                 82 83 84 
           )
    list(APPEND proc${i} "el")
endforeach()           
foreach( i IN ITEMS 54 55 56 57 58
             68 )
    list(APPEND proc${i} "sda")
    list(APPEND proc${i} "dd")
endforeach()  
else()
  set(proc68 "el;sda;dd")
endif()



set(process 100)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_LHC13pAunw ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_LHC13pAunw) # for init
inittest(LHC13pAunw)
foreach( process RANGE 1 100)
  LIST(LENGTH proc${process} LISTCOUNT)
  if (${LISTCOUNT}  GREATER 0)
  foreach(diff ${proc${process}}) 
    foreach ( format ${fmts} )
      set (gencuts ".false.")
      if (${process} EQUAL 3  OR  ${process} EQUAL 7 OR  ${process} EQUAL 8 OR  ${process} EQUAL 14 OR  ${process} EQUAL 15 OR  ${process} EQUAL 16 OR  ${process} EQUAL 70 OR  ${process} EQUAL 71 OR  ${process} EQUAL 72)
        set (gencuts ".true.")
      endif()
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_LHC13pAunw ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_LHC13pAunw_${process}_${diff}_${format})
    endforeach()
    sctestlhe(LHC13pAunw LHC13pAunw_${process}_${diff}_lhe ${diff}_pA ${process})
    sctesthepmc(LHC13pAunw LHC13pAunw_${process}_${diff}_hepmc)
  endforeach()
  list(JOIN proc${process} " " X)
  sctestrivet(LHC13pAunw ${process} ${X})
  endif()
endforeach()
endif()
