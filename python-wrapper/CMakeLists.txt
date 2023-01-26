# manylinux dockers contain only Development.Module, but building on
# macOS-11 fails without full Development
IF (APPLE)
  FIND_PACKAGE(Python 3.2 COMPONENTS Interpreter Development)
ELSE()
  FIND_PACKAGE(Python 3.2 COMPONENTS Interpreter Development.Module)
ENDIF()

FIND_PACKAGE(SWIG 2.0)

if (Python_FOUND AND SWIG_FOUND)
    INCLUDE(${SWIG_USE_FILE})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PYTHON_LIBRARY_OUTPUT_DIRECTORY})
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PYTHON_ARCHIVE_OUTPUT_DIRECTORY})
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/varikn/__init__.py  DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

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
    INCLUDE_DIRECTORIES(${Python_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/..)

    # swig settings
    message("Python version ${Python_VERSION}")
    SET(CMAKE_SWIG_FLAGS "-shadow" "-DPYTHON3" "-py3")

    SET(CMAKE_SWIG_OUTDIR ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
    SET_SOURCE_FILES_PROPERTIES(VariKN.i PROPERTIES CPLUSPLUS ON)
    # create lib
    swig_add_library(varikn_python LANGUAGE python SOURCES VariKN.i )
    swig_link_libraries(varikn_python varikn ${Python_LIBRARIES} ${PLATFORM_LINKFLAGS})
    set_target_properties(${SWIG_MODULE_varikn_python_REAL_NAME} PROPERTIES OUTPUT_NAME "_varikn")
    install(TARGETS ${SWIG_MODULE_varikn_python_REAL_NAME} DESTINATION lib/site-packages)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/../lib/python/varikn/wrapper.py DESTINATION lib/site-packages)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/../lib/python/varikn/__init__.py DESTINATION lib/site-packages)

    message(STATUS "detecting Python")
    message(STATUS "\t Python_INCLUDE_DIRS= ${Python_INCLUDE_DIRS}")
    message(STATUS "\t Python_LIBRARIES= ${Python_LIBRARIES}")
    message(STATUS "\t Python_EXECUTABLE= ${Python_EXECUTABLE}")
    message(STATUS "\t Swig ${SWIG_VERSION}")
elseif(REQUIRE_PYTHON_WRAPPER)
  if (NOT Python_FOUND)
    message(FATAL_ERROR "Unable to find Python libs")
  endif()
  if (NOT SWIG_FOUND)
    message(FATAL_ERROR "Unable to find Swig")
  endif()
endif()
