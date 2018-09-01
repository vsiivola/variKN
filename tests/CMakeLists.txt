FIND_PACKAGE(Boost)
set( TESTPATH ${CMAKE_CURRENT_BINARY_DIR} )
message("Building tests at ${TESTPATH}")

include_directories(${Boost_INCLUDE_DIRS})

add_executable(UnitTests UnitTests.cpp )
target_link_libraries(UnitTests varikn)
set_target_properties( UnitTests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TESTPATH} )

set( TESTDATAPATH ${CMAKE_CURRENT_BINARY_DIR}/data )
file(MAKE_DIRECTORY ${TESTDATAPATH})
add_custom_command( OUTPUT ${TESTDATAPATH}/ax10.txt
                           ${TESTDATAPATH}/bx10.txt
                    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/generate_data.sh ${TESTDATAPATH} )
add_custom_target(do_data DEPENDS ${TESTDATAPATH}/ax10.txt)
add_test( UnitTests  ${TESTPATH}/UnitTests ${TESTDATAPATH})
add_dependencies( UnitTests do_data )
