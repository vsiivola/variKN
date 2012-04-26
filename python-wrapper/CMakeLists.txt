FIND_PACKAGE(PythonLibs)
FIND_PACKAGE(SWIG 2.0)

if (PYTHONLIBS_FOUND AND SWIG_FOUND)
    INCLUDE(${SWIG_USE_FILE})
    set(LIBRARY_OUTPUT_PATH ${LIBRARY_OUTPUT_PATH}/python)
    
    IF (APPLE AND MACPORTS_PYTHON_VER)
        #set( PLATFORM_LINKFLAGS "-framework Foundation")
        set( MACPORTS_ROOT /opt )
    
        set( MACPORTS_PYTHON ${MACPORTS_ROOT}/local/Library/Frameworks/Python.framework )
        IF( EXISTS ${MACPORTS_PYTHON} )
          set(PYTHON_VER ${MACPORTS_PYTHON_VER})
    
          # override the includes
          set( PYTHON_INCLUDE_PATH ${MACPORTS_PYTHON}/Versions/${PYTHON_VER}/include/python${PYTHON_VER})
             
          # dodge the issue of what framework to bind to
          # note: this may cause binary portability problems,
          # but it should be ok if someone builds and uses this
          # on the same machine
          #
          set( PYTHON_LIBRARIES "-undefined dynamic_lookup" )
          
          message("overriding system python for ${MACPORTS_PYTHON}")
        ENDIF( EXISTS ${MACPORTS_PYTHON} )
    ELSEIF(UNIX)
        #set( PLATFORM_LINKFLAGS ${GLIB2_LDFLAGS})
    ENDIF (APPLE AND MACPORTS_PYTHON_VER)
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)
          
   # swig settings
    SET(CMAKE_SWIG_FLAGS "-shadow")
    SET(CMAKE_SWIG_OUTDIR ${LIBRARY_OUTPUT_PATH})
    SET_SOURCE_FILES_PROPERTIES(Perplexity.i PROPERTIES CPLUSPLUS ON)

    # create lib    
    swig_add_module(varikn_python python Perplexity.i )
    swig_link_libraries(varikn_python varikn ${PYTHON_LIBRARIES} ${PLATFORM_LINKFLAGS})
    set_target_properties(${SWIG_MODULE_varikn_python_REAL_NAME} PROPERTIES OUTPUT_NAME "_Perplexity")
    
    # the swig target name is horrid, so make nice top-level name for it
    add_custom_target(varikn_python ALL 
                      DEPENDS ${SWIG_MODULE_varikn_REAL_NAME})
endif(PYTHONLIBS_FOUND AND SWIG_FOUND)
