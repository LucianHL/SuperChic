# - Locate pythia8 library
# Defines:
#
#  PYTHIA8_FOUND
#  PYTHIA8_VERSION
#  PYTHIA8_INCLUDE_DIR
#  PYTHIA8_XMLDOC_DIR
#  PYTHIA8_INCLUDE_DIRS (not cached)
#  PYTHIA8_LIBRARY
#  PYTHIA8_hepmcinterface_LIBRARY
#  PYTHIA8_lhapdfdummy_LIBRARY
#  PYTHIA8_LIBRARIES (not cached) : includes 3 libraries above; not to be used if lhapdf is used

if (PYTHIA8_ROOT_DIR OR PYTHIA8_DIR OR (DEFINED ENV{PYTHIA8_ROOT_DIR}) OR (DEFINED ENV{PYTHIA8_DIR}) )
  set(PYTHIA8_SEARCH_DIRS "" CACHE STRING "" FORCE)
  if (PYTHIA8_ROOT_DIR)
    list (APPEND PYTHIA8_SEARCH_DIRS "${PYTHIA8_ROOT_DIR}" )
  endif()
  if (PYTHIA8_DIR)
    list (APPEND PYTHIA8_SEARCH_DIRS "${PYTHIA8_DIR}" )
  endif()
  if (DEFINED ENV{PYTHIA8_ROOT_DIR})
    list (APPEND PYTHIA8_SEARCH_DIRS "$ENV{PYTHIA8_ROOT_DIR}" )
  endif()
  if (DEFINED ENV{PYTHIA8_DIR})
    list (APPEND PYTHIA8_SEARCH_DIRS "$ENV{PYTHIA8_DIR}" )
endif()
endif()

if (PYTHIA8_SEARCH_DIRS)
  find_path(PYTHIA8_INCLUDE_DIR Pythia.h Pythia8/Pythia.h PATHS  ${PYTHIA8_SEARCH_DIRS} PATH_SUFFIXES include  NO_DEFAULT_PATH)
  find_path(PYTHIA8_XMLDOC_DIR Version.xml PATHS ${PYTHIA8_SEARCH_DIRS} PATH_SUFFIXES xmldoc  share/pythia8/xmldoc share/Pythia8/xmldoc share/pythia8-data/xmldoc  share/doc/packages/pythia/xmldoc NO_DEFAULT_PATH)
  find_library(PYTHIA8_LIBRARY NAMES pythia8 Pythia8 PATHS ${PYTHIA8_SEARCH_DIRS} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
  find_library(PYTHIA8_lhapdfdummy_LIBRARY NAMES lhapdfdummy PATHS ${PYTHIA8_SEARCH_DIRS} PATH_SUFFIXES lib lib64 NO_DEFAULT_PATH)
else()
  find_path(PYTHIA8_INCLUDE_DIR Pythia.h Pythia8/Pythia.h PATH_SUFFIXES include)
  find_path(PYTHIA8_XMLDOC_DIR Version.xml PATH_SUFFIXES xmldoc  share/pythia8/xmldoc share/Pythia8/xmldoc share/pythia8-data/xmldoc  share/doc/packages/pythia/xmldoc)
  find_library(PYTHIA8_LIBRARY NAMES pythia8 Pythia8 PATH_SUFFIXES lib lib64 ../lib ../lib64)
  find_library(PYTHIA8_lhapdfdummy_LIBRARY NAMES lhapdfdummy PATH_SUFFIXES lib lib64 ../lib ../lib64)
endif()

if(PYTHIA8_INCLUDE_DIR AND PYTHIA8_XMLDOC_DIR)
  file(READ ${PYTHIA8_XMLDOC_DIR}/Version.xml versionstr)
  string(REGEX REPLACE ".*Pythia:versionNumber.*default.*[0-9][.]([0-9]+).*" "\\1" PYTHIA8_VERSION "${versionstr}")
  set(PYTHIA8_VERSION "8.${PYTHIA8_VERSION}")
  set(PYTHIA8_INCLUDE_DIRS ${PYTHIA8_INCLUDE_DIR} ${PYTHIA8_INCLUDE_DIR}/Pythia8 ${PYTHIA8_INCLUDE_DIR}/Pythia8Plugins)
  set(PYTHIA8_LIBRARIES ${PYTHIA8_LIBRARY})
  if(PYTHIA8_VERSION VERSION_LESS 8.200)
    #Is this library needed?
    #set(PYTHIA8_LIBRARIES ${PYTHIA8_LIBRARY} ${PYTHIA8_lhapdfdummy_LIBRARY})
  endif()
  find_file(resHEPMC3 HepMC3.h PATHS  ${PYTHIA8_INCLUDE_DIR}/Pythia8Plugins NO_DEFAULT_PATH)
  if (resHEPMC3)
    set(Pythia8_HepMC3_FOUND TRUE)
  endif()
  
endif()

# handle the QUIETLY and REQUIRED arguments and set PYTHIA8_FOUND to TRUE if
# all listed variables are TRUE

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Pythia8 REQUIRED_VARS PYTHIA8_INCLUDE_DIR PYTHIA8_LIBRARIES PYTHIA8_XMLDOC_DIR VERSION_VAR PYTHIA8_VERSION HANDLE_COMPONENTS)


if(Pythia8_FOUND AND NOT TARGET Pythia8::Pythia8)
    add_library(Pythia8::Pythia8 UNKNOWN IMPORTED)
    set_target_properties(Pythia8::Pythia8 PROPERTIES
        IMPORTED_LOCATION "${PYTHIA8_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${PYTHIA8_INCLUDE_DIR}"
    )
endif()
mark_as_advanced(PYTHIA8_FOUND PYTHIA8_INCLUDE_DIR PYTHIA8_LIBRARY PYTHIA8_LIBRARIES PYTHIA8_XMLDOC_DIR)

