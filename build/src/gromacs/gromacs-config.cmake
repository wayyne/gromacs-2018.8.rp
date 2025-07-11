#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2014,2016, by the GROMACS development team, led by
# Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
# and including many others, as listed in the AUTHORS file in the
# top-level source directory and at http://www.gromacs.org.
#
# GROMACS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# GROMACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with GROMACS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to GROMACS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official GROMACS. Details are found
# in the README & COPYING files - if they are missing, get the
# official version at http://www.gromacs.org.
#
# To help us fund GROMACS development, we humbly ask that you cite
# the research papers on the package. Check out http://www.gromacs.org.

set(_gmx_cmake_dir ${CMAKE_CURRENT_LIST_DIR})
get_filename_component(_gmx_root_dir "${_gmx_cmake_dir}" PATH)
get_filename_component(_gmx_root_dir "${_gmx_root_dir}" PATH)
get_filename_component(_gmx_root_dir "${_gmx_root_dir}" PATH)

# Find the exported targets (file name depends on whether shared or static
# libraries were built to allow both to coexist in the same prefix), and
# import them.
set(_gmx_import_file ${_gmx_cmake_dir}/libgromacs.cmake)
if (GROMACS_PREFER_STATIC OR NOT EXISTS ${_gmx_import_file})
    set(_gmx_import_file_static ${_gmx_cmake_dir}/libgromacs_static.cmake)
    if (EXISTS ${_gmx_import_file_static})
        set(_gmx_import_file ${_gmx_import_file_static})
    endif()
    unset(_gmx_import_file_static)
endif()
if (NOT EXISTS ${_gmx_import_file})
    message(FATAL_ERROR
        "The GROMACS installation at ${_gmx_root_dir} does not contain "
        "libgromacs.cmake or libgromacs_static.cmake to define the imported "
        "targets.")
endif()
include(${_gmx_import_file})
unset(_gmx_import_file)

get_target_property(_libs libgromacs INTERFACE_LINK_LIBRARIES)
if (_libs MATCHES "tng_io::tng_io")
    include(CMakeFindDependencyMacro)
    find_dependency(TNG_IO)
endif()
unset(_libs)

set(GROMACS_INCLUDE_DIRS)
set(_include_dirs "include")
foreach (_dir ${_include_dirs})
    if (IS_ABSOLUTE ${_dir})
        list(APPEND GROMACS_INCLUDE_DIRS ${_dir})
    else()
        list(APPEND GROMACS_INCLUDE_DIRS ${_gmx_root_dir}/${_dir})
    endif()
endforeach()
set(GROMACS_LIBRARIES libgromacs )
set(GROMACS_DEFINITIONS -DGMX_DOUBLE=0)
set(GROMACS_IS_DOUBLE OFF)
if (DEFINED GROMACS_SUFFIX AND NOT "${GROMACS_SUFFIX}" STREQUAL "_rp")
    message(FATAL_ERROR "GROMACS_SUFFIX is set inconsistently, expected '_rp'")
endif()
set(GROMACS_SUFFIX "_rp")
set(GROMACS_CXX_COMPILER "/usr/bin/c++")
set(GROMACS_CXX_COMPILER_ID "GNU")
set(GROMACS_CXX_COMPILER_VERSION "4.8.5")
set(GROMACS_CXX_FLAGS "-std=c++11 ")

# Produce a message, since find_package() prints nothing on success.
include(FindPackageMessage)
# The version info is set by CMake when it determines whether this file
# is suitable (by calling the version file, which sets PACKAGE_VERSION).
# TODO: Make also the full version string available from somewhere.
set(_gmx_info "${GROMACS_VERSION}")
if (GROMACS_SUFFIX)
    set(_gmx_info "${_gmx_info} (suffix: ${GROMACS_SUFFIX})")
endif()
find_package_message(GROMACS "Found GROMACS: ${_gmx_info}" "${CMAKE_CURRENT_LIST_FILE}")

unset(_gmx_cmake_dir)
unset(_gmx_root_dir)
unset(_gmx_info)

#####################################################################
# Macros for use in calling code

# This does not work as a function if called as gromacs_check_double(GMX_DOUBLE)
# (i.e., with the parameter value equal to the formal parameter name) because
# of scoping rules.
macro (gromacs_check_double GMX_DOUBLE)
    if (${GMX_DOUBLE} AND NOT GROMACS_IS_DOUBLE)
        message(FATAL_ERROR
            "The found GROMACS installation is compiled in mixed precision, "
            "but double-precision compilation was requested with ${GMX_DOUBLE}=${${GMX_DOUBLE}}")
    elseif (NOT ${GMX_DOUBLE} AND GROMACS_IS_DOUBLE)
        message(FATAL_ERROR
            "The found GROMACS installation is compiled in double precision, "
            "but mixed-precision compilation was requested with ${GMX_DOUBLE}=${${GMX_DOUBLE}}")
    endif()
endmacro()

function (gromacs_check_compiler LANG)
    if (NOT LANG STREQUAL CXX)
        message(FATAL_ERROR
            "gromacs_check_compiler(CXX) is currently the only supported call")
    endif()
    # Deal with possible symlinks (it is fine if one of the used compilers was
    # a symlink to another one).
    get_filename_component(_cmake_compiler_realpath ${CMAKE_${LANG}_COMPILER} REALPATH)
    if (NOT "${_cmake_compiler_realpath}" STREQUAL "${GROMACS_${LANG}_COMPILER}" OR
        NOT "${CMAKE_${LANG}_COMPILER_ID}" STREQUAL "${GROMACS_${LANG}_COMPILER_ID}" OR
        NOT "${CMAKE_${LANG}_COMPILER_VERSION}" STREQUAL "${GROMACS_${LANG}_COMPILER_VERSION}")
        message(WARNING
            "You are compiling with a different C++ compiler from the one that was "
            "used to compile GROMACS. This may lead to linking or runtime problems. "
            "GROMACS was compiled with "
            "${GROMACS_${LANG}_COMPILER_ID} ${GROMACS_${LANG}_COMPILER_VERSION} "
            "(${GROMACS_${LANG}_COMPILER}).")
    endif()
endfunction()
