FIND_PACKAGE(Python 3.2 COMPONENTS Interpreter)

ADD_CUSTOM_COMMAND(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/UnitTests.py" 
  COMMAND cmake -E copy "${CMAKE_CURRENT_SOURCE_DIR}/UnitTests.py" 
  "${CMAKE_CURRENT_BINARY_DIR}" 
  DEPENDS   "${CMAKE_CURRENT_SOURCE_DIR}/UnitTests.py" 
  ) 

FILE ( COPY "${CMAKE_CURRENT_SOURCE_DIR}/UnitTests.py" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}" )
add_test( pytest ${Python_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/UnitTests.py )

# We should make sure that the c unit tests run before trying to run python unit tests.
# However, the following does not work. FIXME !
# add_dependencies( pytest UnitTests )

