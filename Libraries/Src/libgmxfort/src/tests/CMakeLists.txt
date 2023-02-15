# Include the module files from the built libgmxfort.so
include_directories( ${CMAKE_BINARY_DIR}/include )

# Create a library for the tests and link it to the built libgmxfort.so
add_library ( gmxfort_tests SHARED "gmxfort_tests.f90" )
target_link_libraries ( gmxfort_tests ${CMAKE_PROJECT_NAME} )

# Create the test programs and link them to the test library
foreach ( TESTNAME xtcfile ndxfile ndxfile_grp read_next utils )
    add_executable ( ${TESTNAME} ${TESTNAME}.f90 )
    target_link_libraries ( ${TESTNAME} gmxfort_tests )
endforeach ( TESTNAME )

add_test ( "Read_xtcfile" xtcfile )
add_test ( "Read_xtcfile_and_ndxfile" ndxfile )
add_test ( "Read_ndx_grp" ndxfile_grp )
add_test ( "Read_next" read_next )
add_test ( "Utils" utils )

# Copy over files needed for running the tests
configure_file ( test.xtc ${CMAKE_CURRENT_BINARY_DIR}/test.xtc COPYONLY )
configure_file ( test.ndx ${CMAKE_CURRENT_BINARY_DIR}/test.ndx COPYONLY )
