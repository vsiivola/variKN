FIND_PACKAGE(Boost)
message("Running tests")

add_executable(UnitTests UnitTests.cpp )
target_link_libraries(UnitTests varikn)
SET_TARGET_PROPERTIES( UnitTests PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}  
)
message("Building tests to  ${CMAKE_CURRENT_BINARY_DIR} ")
add_test( UnitTests  ${CMAKE_CURRENT_BINARY_DIR}/UnitTests )