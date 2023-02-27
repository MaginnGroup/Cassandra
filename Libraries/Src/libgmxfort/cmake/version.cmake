#
# This file is part of libgmxfort
# https://github.com/wesbarnett/libgmxfort
#
# Copyright (c) 2016,2017 by James W. Barnett
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# Version number is stored in a file
file ( STRINGS "${CMAKE_SOURCE_DIR}/.VERSION" VERSION )
string ( REPLACE "." ";" VERSION_LIST ${VERSION} )
list ( GET VERSION_LIST 0 VERSION_MAJOR )
list ( GET VERSION_LIST 1 VERSION_MINOR )
list ( GET VERSION_LIST 2 VERSION_PATCH )
set ( PROJECT_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}" )
message ( STATUS "CMake build configuration for ${CMAKE_PROJECT_NAME} ${PROJECT_VERSION}" )
string ( TOLOWER ${CMAKE_PROJECT_NAME} PACKAGE_NAME )
set ( PACKAGE_VERSION "${PACKAGE_NAME}-${VERSION}" )
