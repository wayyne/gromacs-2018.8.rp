#
# This file is part of the GROMACS molecular simulation package.
#
# Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

set(DOXYGEN_EXECUTABLE   "/bin/doxygen")
set(DOXYGEN_VERSION      "1.8.5")
set(DOXYGEN_MSCGEN_FOUND )

set(EXPECTED_VERSION     "1.8.5")

if (NOT DOXYGEN_VERSION VERSION_EQUAL EXPECTED_VERSION)
    message("NOTE: You are using Doxygen version ${DOXYGEN_VERSION}. "
            "The documentation is designed for ${EXPECTED_VERSION}. "
            "Other versions may or may not work, but very likely produce extra warnings.")
endif()

message("Running Doxygen...")
# The standard output shows a lot of progress information, but obscures errors
# that may appear.  CMake-introduced buffering also makes it appear sluggish,
# so disable it.
execute_process(COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile-${DOCTYPE}
                OUTPUT_QUIET RESULT_VARIABLE result)
if (result)
    message(FATAL_ERROR "Doxygen failed. "
            "Please run '${DOXYGEN_EXECUTABLE} Doxyfile-${DOCTYPE}' "
            "in the docs/doxygen/ subdirectory in the build tree to see the details.")
endif()

file(READ doxygen-${DOCTYPE}.log DOXYGEN_WARNINGS)
if (DOXYGEN_WARNINGS)
    string(STRIP "${DOXYGEN_WARNINGS}" STRIPPED_WARNINGS)
    message("The following warnings were produced by Doxygen:")
    message("${STRIPPED_WARNINGS}")
    # Remove some useless/hard-to-suppress warnings from the file to avoid
    # Jenkins complaining.
    string(REGEX REPLACE
           "warning:[^\n]*Consider increasing DOT_GRAPH_MAX_NODES.\n" ""
           DOXYGEN_WARNINGS "${DOXYGEN_WARNINGS}")
    # Add a note to the warnings that identify the documentation type that
    # produces them.  This makes it easier to see in Jenkins if a warning type
    # is produced only from some documentation types.
    string(REGEX REPLACE "\n" " (in doxygen-${DOCTYPE})\n"
           DOXYGEN_WARNINGS "${DOXYGEN_WARNINGS}")
    file(WRITE doxygen-${DOCTYPE}.log "${DOXYGEN_WARNINGS}")
endif ()
if (NOT DOXYGEN_MSCGEN_FOUND)
    message("NOTE: mscgen was not available. "
            "Please install it to produce graphs in the documentation.")
endif ()
