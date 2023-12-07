if ( ${CMAKE_VERSION} VERSION_GREATER 3.24)
  cmake_policy(SET CMP0135 OLD)
endif()
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
target_link_libraries(validator PRIVATE HepMC3::HepMC3)
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

set(TABLE_1 " $H(5) \\to b(6)\\,+\\,\\overline{b}(7)$; ")
set(TABLE_2 " $\\gamma(6)\\,+\\, \\gamma(7)$; ")
set(TABLE_3 "$g(6)\\,+\\, g(7)$; ")
set(TABLE_4 "$q(6)\\,+\\, \\overline{q}(7)$; ")
set(TABLE_5 "$c(6)\\,+\\, \\overline{c}(7)$; ")
set(TABLE_6 "$b(6)\\,+\\, \\overline{b}(7)$; ")
set(TABLE_7 "$g(6)\\,+\\, g(7)\\,+\\, g(8)$; ")
set(TABLE_8 "$q(6)\\,+\\, \\overline{q}(7)\\,+\\, g(8)$; ")
set(TABLE_9 "$\\pi^+(6)\\,+\\,\\pi^-(7)$; ")
set(TABLE_10 "$ \\pi^0(6)\\,+\\,\\pi^0(7)$; ")
set(TABLE_11 "$K^+(6)\\,+\\,K^-(7)$; ")
set(TABLE_12 "$K_0(6)\\,+\\,K_0(7)$; ")
set(TABLE_13 "$\\rho_0(6)\\,+\\,\\rho_0(7)$; ")
set(TABLE_14 "$\\eta(6)\\,+\\,\\eta(7)$; ")
set(TABLE_15 "$\\eta(6)\\,+\\,\\eta'(7)$; ")
set(TABLE_16 "$\\eta'(6)\\,+\\,\\eta'(7)$; ")
set(TABLE_17 "$\\phi(6)\\,+\\,\\phi(7)$; ")
set(TABLE_18 "$J/\\psi(6)(\\to \\mu^+(8)+\\mu^-(9))\\,+\\,J/\\psi(7)(\\to \\mu^+(10)+\\mu^-(11))$; ")
set(TABLE_19 "$J/\\psi(6)(\\to \\mu^+(8)+\\mu^-(9))\\,+\\,\\psi_{2S}(7)(\\to \\mu^+(10)+\\mu^-(11))$; ")
set(TABLE_20 "$\\psi_{2S}(6)(\\to \\mu^+(8)+\\mu^-(9))\\,+\\,\\psi_{2S}(7)(\\to \\mu^+(10)+\\mu^-(11))$; ")
set(TABLE_21 "$\\chi_{c0}(5)\\to \\gamma(6)\\,+\\,J/\\psi(7)(\\to \\mu^+(8)+\\mu^-(9))$; ")
set(TABLE_22 "$\\chi_{c1}(5)\\to \\gamma(6)\\,+\\,J/\\psi(7)(\\to \\mu^+(8)+\\mu^-(9))$; ")
set(TABLE_23 "$\\chi_{c2}(5)\\to \\gamma(6)\\,+\\,J/\\psi(7)(\\to \\mu^+(8)+\\mu^-(9))$; ")
set(TABLE_24 "$\\chi_{c0}(5)\\to S(6)\\,+\\,S(7)$; ")
set(TABLE_25 "$\\chi_{c1}(5)\\to S(6)\\,+\\,S(7)$; ")
set(TABLE_26 "$\\chi_{c2}(5)\\to S(6)\\,+\\,S(7)$; ")
set(TABLE_27 "$\\chi_{c1}(5)\\to f(6)\\,+\\,\\overline{f}(7)$; ")
set(TABLE_28 "$\\chi_{c2}(5)\\to f(6)\\,+\\,\\overline{f}(7)$; ")
set(TABLE_29 "$\\chi_{c0}(5)\\to \\pi^+(6)+\\pi^-(7)+\\pi^+(8)+\\pi^-(9)$; ")
set(TABLE_30 "$\\chi_{c1}(5)\\to \\pi^+(6)+\\pi^-(7)+\\pi^+(8)+\\pi^-(9)$; ")
set(TABLE_31 "$\\chi_{c2}(5)\\to \\pi^+(6)+\\pi^-(7)+\\pi^+(8)+\\pi^-(9)$; ")
set(TABLE_32 "$\\chi_{c0}(5)\\to \\pi^+(6)+\\pi^-(7)+K^+(8)+K^-(9)$; ")
set(TABLE_33 "$\\chi_{c1}(5)\\to \\pi^+(6)+\\pi^-(7)+K^+(8)+K^-(9)$; ")
set(TABLE_34 "$\\chi_{c2}(5)\\to \\pi^+(6)+\\pi^-(7)+K^+(8)+K^-(9)$; ")
set(TABLE_35 "$\\chi_{c0}(5)\\to 3(\\pi^+(6,8,10)\\,+\\,\\pi^-(7,9,11))$; ")
set(TABLE_36 "$\\chi_{c1}(5)\\to 3(\\pi^+(6,8,10)\\,+\\,\\pi^-(7,9,11))$; ")
set(TABLE_37 "$\\chi_{c2}(5)\\to 3(\\pi^+(6,8,10)\\,+\\,\\pi^-(7,9,11))$; ")
set(TABLE_38 "$\\eta_c(5)$; ")
set(TABLE_39 "$\\chi_{b0}(5)\\to \\gamma(6)\\,+\\,\\Upsilon_{1S}(7)(\\to \\mu^+(8)+\\mu^-(9))$; ")
set(TABLE_40 "$\\chi_{b1}(5)\\to \\gamma(6)\\,+\\,\\Upsilon_{1S}(7)(\\to \\mu^+(8)+\\mu^-(9))$; ")
set(TABLE_41 "$\\chi_{b2}(5)\\to \\gamma(6)\\,+\\,\\Upsilon_{1S}(7)(\\to \\mu^+(8)+\\mu^-(9))$; ")
set(TABLE_42 "$\\chi_{b0}(5)\\to S(6)\\,+\\,S(7)$; ")
set(TABLE_43 "$\\chi_{b1}(5)\\to S(6)\\,+\\,S(7)$; ")
set(TABLE_44 "$\\chi_{b2}(5)\\to S(6)\\,+\\,S(7)$; ")
set(TABLE_45 "$\\chi_{b1}(5)\\to f(6)\\,+\\,\\overline{f}(7)$; ")
set(TABLE_46 "$\\chi_{b2}(5)\\to f(6)\\,+\\,\\overline{f}(7)$; ")
set(TABLE_47 " $\\eta_b(5)$; ")
set(TABLE_48 "$\\rho_0(5)\\to \\pi^+(6)\\pi^-(7)$; ")
set(TABLE_49 "$\\phi(5)\\to K^+(6)K^-(7)$; ")
set(TABLE_50 "$J/\\psi(5)\\to \\mu^+(6)\\mu^-(7)$; ")
set(TABLE_51 "$\\Upsilon_{1S}(5)\\to \\mu^+(6)\\mu^-(7)$; ")
set(TABLE_52 "$\\psi_{2S}(5)\\to \\mu^+(6)\\mu^-(7)$; ")
set(TABLE_53 "$\\psi_{2S}(5)\\to J/\\psi(6)(\\to\\mu^+(9)+\\mu^-(10))\\,+\\,\\pi^+(7)\\,+\\,\\pi^-(8)$; ")
set(TABLE_54 "Not documented; ")
set(TABLE_55 "$W^+(\\to \\nu_l(8)+l^+(9))\\,+\\, W^-(\\to \\overline{\\nu}_l(10)+l^-(11))$; ")
set(TABLE_56 "$e^+(6)\\,+\\,e^-(7)$; ")
set(TABLE_57 "$\\mu^+(6)\\,+\\,\\mu^-(7)$; ")
set(TABLE_58 "$\\tau^+(6)\\,+\\,\\tau^-(7)$; ")
set(TABLE_59 "$\\gamma(6)\\,+\\,\\gamma(7)$; ")
set(TABLE_60 "$H(5) \\to b(6)\\,+\\,\\overline{b}(6)$; ")
set(TABLE_61 "Not documented; ")
set(TABLE_68 "$a(5)\\to \\gamma(6)\\,+\\,\\gamma(7)$; ")
set(TABLE_69 "$M(5) \\to \\gamma(6)\\,+\\,\\gamma(7)$ (Dirac Coupling); ")
set(TABLE_70 "$M(5) \\to \\gamma(6)\\,+\\,\\gamma(7)$ ($\\beta g$ Coupling); ")
set(TABLE_71 "$m(6)\\,+\\,\\overline{m}(7)$ (Dirac Coupling); ")
set(TABLE_72 "$m(6)\\,+\\,\\overline{m}(7)$ ($\\beta g$ Coupling); ")
set(TABLE_73 "$\\tilde{\\chi}^-(6)(\\to \\tilde{\\chi}_0^1 (8)+\\mu^-(9)+\\overline{\\nu}_\\mu (10)) \\,+\\,\\tilde{\\chi}^+(7)(\\to \\tilde{\\chi}_0^1 (11)+\\mu^+(12)+\\nu_\\mu (13))$; ")
set(TABLE_74 "$\\tilde{\\chi}^-(6)(\\to \\tilde{\\chi}_0^1 (8)+\\overline{u}(9)+d(10)) \\,+\\,\\tilde{\\chi}^+(7)(\\to \\tilde{\\chi}_0^1 (11)+u(12)+\\overline{d} (13))$; ")
set(TABLE_75 "$\\tilde{\\chi}^-(6)(\\to \\tilde{\\chi}_0^1 (8)+\\mu^-(9)+\\overline{\\nu}_\\mu (10)) \\,+\\,\\tilde{\\chi}^+(7)(\\to \\tilde{\\chi}_0^1 (11)+u(12)+\\overline{d} (13))$; ")
set(TABLE_76 "$\\tilde{l}^-(5))(\\to \\tilde{\\chi}_0^1 (8)+\\mu^-(9)) \\,+\\,\\tilde{l}^+(6) (\\to \\tilde{\\chi}_0^1 (10)+\\mu^+(11))$; ")
set(TABLE_77 "$\\phi(5)\\to \\mu^+(6)\\mu^-(7)$; ")
set(TABLE_78 "$J/\\psi(5)\\to e^+(6)e^-(7)$; ")
set(TABLE_79 "$\\psi_{2S}(5)\\to e^+(6)e^-(7)$; ")
set(TABLE_82 "$X(6) + \\gamma(7)$; ")
set(TABLE_83 "$X(6) + Z(7)(\\to \\mu^+ (8)+\\mu^- (9))$; ")
set(TABLE_84 "Not documented; ")
set(TABLE_43 "$X(6) + Z(7)(\\to e^+ (8)+ e^- (9))$; ")



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
           50 51 52 53    55 56 57 58 59
           60                      68 69
           70 71 72 73 74 75 76 77 78 79
                 82 83 84 
           )
    list(APPEND proc${i} "el")
endforeach()           
foreach( i IN ITEMS    55 56 57 58
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
    list(TRANSFORM TABLE_${process} PREPEND " $pp$-${diff} " AT 1)
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
           60                      68 69
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
    list(TRANSFORM TABLE_${process} PREPEND " $ee$-${diff} " AT 1)
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
foreach( i IN ITEMS  
           55 56 57 58 59
           60 68 69
           70 71 72 73 74 75 76
           )
    list(APPEND proc${i} "el")
endforeach()           
else()
  set(proc68 "el")
endif()



set(process 100)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/input.DAT_LHC13AAunw ${CMAKE_CURRENT_BINARY_DIR}/input.DAT_LHC13AAunw) # for init
inittest(LHC13AAunw)
foreach( process RANGE 1 100)
  LIST(LENGTH proc${process} LISTCOUNT)
  if (${LISTCOUNT}  GREATER 0)
  foreach(diff ${proc${process}}) 
    list(TRANSFORM TABLE_${process} PREPEND " $AA$-${diff} " AT 1)
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
           50 51 52 53    55 56 57 58 59
           60                      68 69
           70 71 72 73 74 75 76 77 78 79
                 82 83 84 
           )
    list(APPEND proc${i} "el")
endforeach()           
foreach( i IN ITEMS    55 56 57 58
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
    list(TRANSFORM TABLE_${process} PREPEND " $pA$-${diff} " AT 1)
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


write_file(${CMAKE_CURRENT_BINARY_DIR}/tested.tex 
"\\documentclass{article}\n\
\\usepackage[margin=0.7in]{geometry}\n\
\\usepackage[parfill]{parskip}\n\
\\usepackage[utf8]{inputenc}\n\  
%\\usepackage{table}\n\  
    % Related to math\n\
\\usepackage{amsmath,amssymb,amsfonts,amsthm}\n\
\\begin{document}\n\
\\begin{center}\\tiny\n\
\\begin{tabular}{|c|l|l|}\\hline\n\
Number&Final--State& Beams and types\\\\\\hline\\hline\n\
")
foreach( process RANGE 1 100)
  if (TABLE_${process})
  list(JOIN TABLE_${process} " & " X)
  write_file(${CMAKE_CURRENT_BINARY_DIR}/tested.tex "${process} & ${X}\\\\\n" APPEND)
  endif()
endforeach()
write_file(${CMAKE_CURRENT_BINARY_DIR}/tested.tex 
"\\hline\n\
\\end{tabular}\n\
\\end{center}\n\ 
\\end{document}\n" 
APPEND)
find_package(LATEX COMPONENTS PDFLATEX)
if (LATEX_FOUND)
    ADD_CUSTOM_TARGET(   latexcompiletest
                       COMMAND ${PDFLATEX_COMPILER} tested.tex 
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                      )
    ADD_CUSTOM_TARGET(all-test ALL) # Entry point of execution.
    ADD_DEPENDENCIES( all-test latexcompiletest)

endif()


#cat doc/superchic.tex | grep ^[1-9].*  | cut -f 1,2 -d'&' --output-delimiter='";"' | tr -s ' ' | sed 's/$/\";\" \")/g' | sed 's/^/set(TABLE_/g' | sed 's/\";/ /1'
