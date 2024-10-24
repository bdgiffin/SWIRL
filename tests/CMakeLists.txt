set(test_dependencies SWIRL ${SWIRL_dependencies} gtest)
set(test_include_dirs ${SWIRL_include_dirs} ../src)

blt_add_executable(NAME       SpatialHash_tests
                   SOURCES    ${CMAKE_CURRENT_LIST_DIR}/SpatialHash_tests.cpp
                   INCLUDES   ${test_include_dirs}
                   DEPENDS_ON ${test_dependencies}
		   )
		 
blt_add_test(NAME    SpatialHash_tests
	     COMMAND SpatialHash_tests)
