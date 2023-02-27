cmake_minimum_required ( VERSION 2.8.8 )
enable_language (C)

project ( xdrfile NONE )

# Set the type/configuration of build to perform
set ( CMAKE_CONFIGURATION_TYPES "Debug" "Release" "MinSizeRel" "RelWithDebInfo" )
set ( CMAKE_BUILD_TYPE "Release"
  CACHE STRING "Select which configuration to build." )

set ( PROJECT_DESCRIPTION "Reads/writes compressed trajectory files" )
set ( PROJECT_URL "https://github.com/wesbarnett/libxdrfile" )
set ( LIBDIR "lib")
set ( INCLUDEDIR "include")
set ( INSTALL_LIBDIR ${CMAKE_INSTALL_PREFIX}/${LIBDIR} )
set ( INSTALL_INCLUDEDIR ${CMAKE_INSTALL_PREFIX}/${INCLUDEDIR} )

file ( STRINGS "${CMAKE_SOURCE_DIR}/.VERSION" VERSION )
string( REPLACE "." ";" VERSION_LIST ${VERSION} )
list(GET VERSION_LIST 0 VERSION_MAJOR)
list(GET VERSION_LIST 1 VERSION_MINOR)
list(GET VERSION_LIST 2 VERSION_PATCH)
set(PROJECT_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
message ( STATUS "CMake build configuration for ${CMAKE_PROJECT_NAME} ${PROJECT_VERSION}" )

string ( TOLOWER ${CMAKE_PROJECT_NAME} PACKAGE_NAME )
set ( PACKAGE_VERSION "${PACKAGE_NAME}-${VERSION}" )

file ( GLOB SRCS "src/*.c" )
add_library ( ${CMAKE_PROJECT_NAME} SHARED ${SRCS} )

set_target_properties ( ${CMAKE_PROJECT_NAME}
  PROPERTIES
  OUTPUT_NAME ${CMAKE_PROJECT_NAME}
  PREFIX lib
  SOVERSION ${VERSION_MAJOR}.${VERSION_MINOR}
  VERSION ${VERSION}
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${LIBDIR}"
)

install ( TARGETS ${CMAKE_PROJECT_NAME} EXPORT ${PACKAGE_NAME}-targets LIBRARY DESTINATION ${LIBDIR} )
install ( DIRECTORY include/ DESTINATION ${INCLUDEDIR} FILES_MATCHING PATTERN "*.h")

# pkg-config
configure_file(
   "${CMAKE_CURRENT_SOURCE_DIR}/lib${CMAKE_PROJECT_NAME}.pc.cmake.in"
   "${CMAKE_CURRENT_BINARY_DIR}/lib${CMAKE_PROJECT_NAME}.pc"
   @ONLY
)
install ( FILES
   "${CMAKE_CURRENT_BINARY_DIR}/lib${CMAKE_PROJECT_NAME}.pc"
   DESTINATION "${LIBDIR}/pkgconfig"
)

#-----------------------------------------------------
# Publicize installed location to other CMake projects
#-----------------------------------------------------

set ( EXPORT_INSTALL_DIR "${INSTALL_LIBDIR}/cmake/${PACKAGE_VERSION}" )

include ( CMakePackageConfigHelpers ) # Standard CMake module
write_basic_package_version_file(
"${CMAKE_BINARY_DIR}/${PACKAGE_NAME}-config-version.cmake"
  VERSION ${VERSION}
  COMPATIBILITY SameMajorVersion )

# install package config file
configure_package_config_file (
  "${CMAKE_SOURCE_DIR}/cmake/pkg/${CMAKE_PROJECT_NAME}-config.cmake.in"
  "${CMAKE_BINARY_DIR}/pkg/${PACKAGE_NAME}-config.cmake"
  INSTALL_DESTINATION "${EXPORT_INSTALL_DIR}"
  PATH_VARS EXPORT_INSTALL_DIR)

# Install the config and version files so that we can find this project with
# others
install ( FILES
  "${CMAKE_BINARY_DIR}/pkg/${PACKAGE_NAME}-config.cmake"
  "${CMAKE_BINARY_DIR}/${PACKAGE_NAME}-config-version.cmake"
  DESTINATION "${EXPORT_INSTALL_DIR}" )

#-----------------------------------------------------
# Tests
#-----------------------------------------------------
enable_testing ()
include_directories ( include )
set ( TEST_NAME "${CMAKE_PROJECT_NAME}_test")
add_executable ( ${TEST_NAME} src/tests/test.c )
target_link_libraries ( ${TEST_NAME} m )
target_link_libraries ( ${TEST_NAME} ${CMAKE_PROJECT_NAME} )
add_test ( "basic_xdrfile_library" ${TEST_NAME} 0)
add_test ( "xtc_functionality"  ${TEST_NAME} 1)
add_test ( "trr_functionality" ${TEST_NAME} 2)
