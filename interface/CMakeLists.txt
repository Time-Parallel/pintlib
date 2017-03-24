#
# build the python shared object now
#
SET(LIBTS_SRC ${CMAKE_CURRENT_SOURCE_DIR})
SET(LIBTS_INTERFACEFILE ${LIBTS_SRC}/tsLibInterface.f90)
SET(LIBTS_F2PYFILE ${LIBTS_BUILD}/interface/libTS.pyf)
SET(LIBTS_FWRAPFILE ${LIBTS_BUILD}/interface/libTS-f2pywrappers2.f90)
SET(LIBTS_CWRAPFILE ${LIBTS_BUILD}/interface/libTSmodule.c)
SET(NUMPY_SRCFILE ${NUMPY_SRC}/fortranobject.c)
#
# All of these set in standalone.cmake automatically.  PYTHON_INCLUDE_DIRS is
# found through FIND_PACKAGE(PythonLibs).
#
include_directories(AFTER ${NUMPY_SRC} ${NUMPY_INCLUDE} ${PYTHON_INCLUDE})
#
# Include the build directory for this source directory.
#
include_directories(BEFORE ${LIBTS_BUILD}/interface)
#
# add custom commands to generate the .pyf file from cartInterface.f90
#
ADD_CUSTOM_COMMAND(OUTPUT ${LIBTS_F2PYFILE} 
                   COMMAND ${F2PYBIN} -m libTS --overwrite-signature -h ${LIBTS_F2PYFILE} ${LIBTS_INTERFACEFILE}
                   MAIN_DEPENDENCY ${LIBTS_INTERFACEFILE})
#
# add custom command to build the C and f90 source files that
# come out of the f2py wrapper
#
ADD_CUSTOM_COMMAND(OUTPUT ${LIBTS_FWRAPFILE} ${LIBTS_CWRAPFILE}
                   COMMAND ${F2PYBIN} ARGS ${LIBTS_F2PYFILE} 
                   MAIN_DEPENDENCY ${LIBTS_F2PYFILE})
#
# build the libTS.so shared object
#
ADD_LIBRARY(libTS SHARED ${LIBTS_FWRAPFILE} ${LIBTS_CWRAPFILE} 
		  ${NUMPY_SRCFILE})
SET_TARGET_PROPERTIES(libTS PROPERTIES
                      LIBRARY_OUTPUT_DIRECTORY ${LIBTS_BUILD}/python)
#
# Make shared object link against Python and libTS.a
#
TARGET_LINK_LIBRARIES(libTS LIBTS ${PYTHON_LIB_FLAGS})
#
# Sometimes linker has issues finding Fortran runtimes unless Fortran compiler
# invokes the linker directly, so make it call the linker instead of C.
#
SET_TARGET_PROPERTIES(libTS PROPERTIES LINKER_LANGUAGE Fortran)
#
# Make the format libTS.so for the output library.  On OS X we want to force
# the so extension instead of dylib.
#
SET_TARGET_PROPERTIES(libTS PROPERTIES PREFIX "")
SET_TARGET_PROPERTIES(libTS PROPERTIES SUFFIX ".so")
