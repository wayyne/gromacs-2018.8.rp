/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \file
 * \brief
 * Version information for software that links to \Gromacs.
 *
 * \if libapi
 * This include file will be configured by CMake and contains version
 * information.  It is not used by \Gromacs, but intended for software that
 * links to \Gromacs.
 * The values come from the main CMakeLists.txt.
 * \endif
 *
 * This file exists from 4.6 onward, and can be included as
 * `<gromacs/version.h>`.  In 4.6, it is also included by
 * `<gromacs/typedefs.h>, but that header has already moved in 5.0.
 *
 * This header defines two values, the \Gromacs version, and the API version.
 * The versions are in numerical form, where, for example, version
 * 4.6.1 would be 40601.
 *
 * The API version is defined in ::GMX_API_VERSION, and denotes the
 * version of the programmer interface, i.e. the installed header files
 * and compatible library.
 *
 * Programs written against the \Gromacs library can use this file
 * to provide some backward compatibility even though parts of the API
 * change.  For example:
 * \code
   #include <gromacs/version.h>
   #if (GMX_API_VERSION < 50000)
       .... <do pre-5.0 stuff>
   #else
       .... <do post-5.0 stuff>
   #endif
   \endcode
 * where version.h is included directly. For code that must be compatible
 * between 4.5 and 4.6, an interim solution is to include typedefs.h, which
 * includes this file:
 * \code
   #include <gromacs/typedefs.h>
   #if !defined(GMX_API_VERSION) || (GMX_API_VERSION < 40600)
       ....  <do 4.5 specific stuff>
   #elif (GMX_API_VERSION < 40700)
       ....  <do 4.6 specific stuff>
   #endif
   \endcode
 *
 * \inpublicapi
 */
#ifndef GMX_VERSION_H
#define GMX_VERSION_H

/*! \brief
 * API version of this set of \Gromacs headers.
 *
 * If there are multiple versions of \Gromacs that work with the same set of
 * headers, then this version is not updated between the versions, even though
 * ::GMX_VERSION is.
 * For 4.6 and 5.0 (and likely for some time in the future as well), this
 * tracks the exact \Gromacs version.
 */
#define GMX_API_VERSION 20180008

/*! \brief
 * Exact \Gromacs version of this set of headers.
 *
 * This specifies the version number of the actual \Gromacs library that
 * installed these headers.
 */
#define GMX_VERSION 20180008

#endif
