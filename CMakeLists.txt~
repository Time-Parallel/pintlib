#
# Contains a CMakeLists.txt file that should work for standalone libTS.
# Jay Sitaraman 11/22/2014
#
PROJECT(LIBTS)
ENABLE_LANGUAGE(Fortran C)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
#
# Include a file that replaces functionality of Make.defs.
#
INCLUDE("standalone.cmake")
#
# Each directory includes a CMake include file that lists
# any source file which should be compiled in the directory.
# Dependencies between Fortran modules are computed automatically,
# so there is no need to build one directory before another.
#
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/lhs/lhs.cmake")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/math/math.cmake")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/src/partition/partition.cmake")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/interface/interface.cmake")
#
# Set the top level source directory for LIBTS.
#
SET(LIBTS_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}) 
#
# Set the build directory - where object files are place when compiled.
#
SET(LIBTS_BUILD ${CMAKE_CURRENT_BINARY_DIR})
#
# create a default installation dir if not provided
#
IF (NOT DEFINED INSTALL_DIR)
 SET(INSTALL_DIR ${LIBTS_BUILD}/libTS)
ENDIF()
#
# Fortran module files will be built here, so make sure preprocessor can see
# them during compile.
#
INCLUDE_DIRECTORIES(BEFORE ${LIBTS_BUILD})
#
# LIBTS source files without driver but 
# includes interface
#
SET(LIBTS_SRCS ${LHSSRCS} ${MATHSRCS} ${PARTITIONSRCS} ${INTERFACESRCS})
#
# Build the basic libTS.a one expects.
#
ADD_LIBRARY(LIBTS STATIC ${LIBTS_SRCS})
SET_TARGET_PROPERTIES(LIBTS PROPERTIES
                      ARCHIVE_OUTPUT_DIRECTORY ${LIBTS_BUILD}/lib)
#
# include interface directory to 
# build the python shared object now
#
ADD_SUBDIRECTORY(interface)
#
# install all the required files in 
# the installation directory
#
INSTALL(FILES ${LIBTS_SOURCE}/README.md
  DESTINATION ${INSTALL_DIR}/README)
INSTALL(FILES ${LIBTS_SOURCE}/interface/python/libTSInterface.py 
  DESTINATION ${INSTALL_DIR}/bin)
INSTALL(FILES ${LIBTS_BUILD}/python/libTS.so 
  DESTINATION ${INSTALL_DIR}/bin)
INSTALL(DIRECTORY ${LIBTS_BUILD}/lib
  DESTINATION ${INSTALL_DIR})
