FIND_PACKAGE(PythonInterp) 
FIND_PACKAGE(PythonLibs)
FIND_PACKAGE(SWIG 2.0)

#################################################
# An ugly hack to force cmake to use python3.
# Not needed with a recent version of cmake
# You need to set the paths manually corresponding
# to where your py3 is installed.
#if ("${CMAKE_USE_PYTHON_VERSION}" VERSION_GREATER "3.2")
#set(PYTHON_INCLUDE_PATH /opt/python3/include/python3.3m)
#set(PYTHON_LIBRARIES /opt/python3/lib/libpython3.3m.so)
#set(PYTHON_EXECUTABLE /opt/python3/bin/python3)
#set(PYTHON_VERSION ${CMAKE_USE_PYTHON_VERSION})
#endif()

# Cmake still detects the old python2 executable first
set(PYTHON_EXECUTABLE "${PYTHON_EXECUTABLE}3")

if (PYTHONLIBS_FOUND AND SWIG_FOUND)
    INCLUDE(${SWIG_USE_FILE})
    set(LIBRARY_OUTPUT_PATH ${LIBRARY_OUTPUT_PATH}/python)
    
    IF (APPLE AND MACPORTS_PYTHON_VERSION)
        #set( PLATFORM_LINKFLAGS "-framework Foundation")
        set( MACPORTS_ROOT /opt )
    
        set( MACPORTS_PYTHON ${MACPORTS_ROOT}/local/Library/Frameworks/Python.framework )
        IF( EXISTS ${MACPORTS_PYTHON} )
          set(PYTHON_VERSION ${MACPORTS_PYTHON_VERSION})
    
          # override the includes
          if ( MACPORTS_PYTHON_VERSION VERSION_GREATER 3.1 )
            set( PYTHON_INCLUDE_PATH ${MACPORTS_PYTHON}/Versions/${PYTHON_VERSION}/include/python${PYTHON_VERSION}m)
          else()
            set( PYTHON_INCLUDE_PATH ${MACPORTS_PYTHON}/Versions/${PYTHON_VERSION}/include/python${PYTHON_VERSION})
          endif()
          # dodge the issue of what framework to bind to
          # note: this may cause binary portability problems,
          # but it should be ok if someone builds and uses this
          # on the same machine
          #
          set( PYTHON_LIBRARIES "-undefined dynamic_lookup" )
          
          message("overriding system python for ${MACPORTS_PYTHON}")
        ENDIF( EXISTS ${MACPORTS_PYTHON} )
    ENDIF (APPLE AND MACPORTS_PYTHON_VERSION)

    ###########################
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
          
   # swig settings
   message("Python version ${PYTHON_VERSION}")
   if(PYTHON_INCLUDE_PATH MATCHES "(python3)")
      message("Using python ver >2, disabling FILE * typemaps")
      SET(CMAKE_SWIG_FLAGS "-shadow" "-DPYTHON3" "-py3")
    else( )
      SET(CMAKE_SWIG_FLAGS "-shadow")    
    endif()

    SET(CMAKE_SWIG_OUTDIR ${LIBRARY_OUTPUT_PATH})
    SET_SOURCE_FILES_PROPERTIES(VariKN.i PROPERTIES CPLUSPLUS ON)
    # create lib    
    swig_add_module(varikn_python python VariKN.i )
    swig_link_libraries(varikn_python varikn ${PYTHON_LIBRARIES} ${PLATFORM_LINKFLAGS})
    set_target_properties(${SWIG_MODULE_varikn_python_REAL_NAME} PROPERTIES OUTPUT_NAME "_varikn")
    # make nice top-level target name
    add_custom_target(varikn_python ALL 
                      DEPENDS ${SWIG_MODULE_varikn_REAL_NAME})

    message(STATUS "detecting Python")
    message(STATUS "\t PYTHON_INCLUDE_PATH= ${PYTHON_INCLUDE_PATH}")
    message(STATUS "\t PYTHON_LIBRARIES= ${PYTHON_LIBRARIES}")
    message(STATUS "\t PYTHON_EXECUTABLE= ${PYTHON_EXECUTABLE}")
    message(STATUS "\t Swig ${SWIG_VERSION}")    
    if (NOT PYTHONLIBS_FOUND)
      message(FATAL_ERROR "Unable to find Python libs")
    elseif (NOT SWIG_FOUND)
      message(FATAL_ERROR "Unable to find Swig")
    endif (NOT PYTHONLIBS_FOUND)
endif(PYTHONLIBS_FOUND AND SWIG_FOUND)
