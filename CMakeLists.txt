cmake_minimum_required(VERSION 3.12)
PROJECT(varikn)

SET(ENABLE_EXECUTABLES ON CACHE BOOL "Set false to skip building the executables")
SET(ENABLE_TESTING OFF CACHE BOOL "Set true to enable testing")
SET(ENABLE_PYTHON_WRAPPER ON CACHE BOOL "Set false to skip building the Python wrapper")
SET(REQUIRE_PYTHON_WRAPPER OFF CACHE BOOL "Set true to fail if Python wrapper cannot be installed")

set( CMAKE_CXX_STANDARD 11)

cmake_policy(SET CMP0054 NEW)
cmake_policy(SET CMP0078 OLD) # Should fix to work with NEW policy
cmake_policy(SET CMP0086 NEW)

include(cmake/clang-cxx-dev-tools.cmake)

macro( check_environment )
    if(WIN32)
        if (${CMAKE_GENERATOR} MATCHES "NMake")
            if ("$ENV{VSINSTALLDIR}" STREQUAL "")
                message(FATAL_ERROR "MSVC environment not set.  You need to run \"%VS90COMNTOOLS%\\vsvars32.bat\" manually before running cmake.")
            endif ("$ENV{VSINSTALLDIR}" STREQUAL "")
        endif(${CMAKE_GENERATOR} MATCHES "NMake")
    endif(WIN32)
endmacro( check_environment )

check_environment()

# detect whether we're in-source or out-of-source
#
IF (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    set(IN_SOURCE_BUILD 1)
    message("In-source builds aren't recommended!")
ENDIF(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
message("Configuring for platform ${CMAKE_SYSTEM_NAME}")
message("Build type is ${CMAKE_BUILD_TYPE}")
message("Generator is ${CMAKE_GENERATOR}")

IF(APPLE)
   set( CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}" )
ELSEIF(UNIX)
   SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
   SET( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
   # suppress -rdynamic (see http://www.itk.org/Bug/view.php?id=9985)
   SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
   SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
ELSEIF(WIN32)
   #add_debug_flags(/Zi)
   #add_release_flags(/O2 /Ox)

   # we always want to statically link to the CRT, so
   # override the default dynamic CRT linkage (/MD) with
   # static linkage (/MT)
   #replace_flags(/MD /MT)
ENDIF(APPLE)

# setup custom output dirs (unless defined)
#
IF(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
ENDIF()
IF(NOT DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
  set(PYTHON_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/python/varikn)
ELSE()
  set(PYTHON_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
ENDIF()
IF(NOT DEFINED CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
  set(PYTHON_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib/python/varikn)
ELSE()
  set(PYTHON_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
ENDIF()
file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
file(MAKE_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src )
file(GLOB c_srcs   RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/[^.]*.c")
file(GLOB cpp_srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/[^.]*.cc")
file(GLOB chdr_srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/[^.]*.h")
file(GLOB cpphdr_srcs RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "src/[^.]*.hh")
IF(NOT WIN32)
  SET_SOURCE_FILES_PROPERTIES(${c_srcs} PROPERTIES LANGUAGE C)
ELSE(NOT WIN32)
  SET_SOURCE_FILES_PROPERTIES(${c_srcs} PROPERTIES LANGUAGE CXX)
ENDIF(NOT WIN32)
set(srcs ${c_srcs} ${cpp_srcs})
set(hdrs ${chdr_srcs} ${cpphdr_srcs})
set(LIBVARIKN_SRCS ${srcs} ${hdrs})
add_library(varikn STATIC ${LIBVARIKN_SRCS} )
set_target_properties(varikn PROPERTIES PREFIX "lib")
install(TARGETS varikn  ARCHIVE DESTINATION lib/${OUTPUT_DIR_NAME})

# other targets
#
IF(ENABLE_EXECUTABLES)
  add_subdirectory(app)
ENDIF(ENABLE_EXECUTABLES)

IF(ENABLE_PYTHON_WRAPPER)
  add_subdirectory(python-wrapper)
ENDIF(ENABLE_PYTHON_WRAPPER)

IF(ENABLE_TESTING)
  enable_testing()
  add_subdirectory(tests)
  IF(ENABLE_PYTHON_WRAPPER)
    add_subdirectory(python-wrapper/tests)
  ENDIF(ENABLE_PYTHON_WRAPPER)
ENDIF(ENABLE_TESTING)

# These are so that package name and install prefix can be overriden
SET(DEB_INSTALL_PREFIX CACHE STRING "/usr")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "variKN statistical language modeling toolkit")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_GENERATOR "DEB")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Vesa Siivola")
SET(CPACK_PACKAGING_INSTALL_PREFIX ${DEB_INSTALL_PREFIX})

INCLUDE(CPack)
