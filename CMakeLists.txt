cmake_minimum_required(VERSION 2.6)

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

# setup custom output dirs
#
set( EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_BINARY_DIR}/bin/${OUTPUT_DIR_NAME} )
set( LIBRARY_OUTPUT_PATH    ${CMAKE_CURRENT_BINARY_DIR}/lib/${OUTPUT_DIR_NAME} )
file(MAKE_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
file(MAKE_DIRECTORY ${LIBRARY_OUTPUT_PATH})

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
add_subdirectory(app)
add_subdirectory(python-wrapper)   
IF(CMAKE_ENABLE_TESTING)
  enable_testing()
  add_subdirectory(tests)  
  add_subdirectory(python-wrapper/tests)
ENDIF(CMAKE_ENABLE_TESTING)

