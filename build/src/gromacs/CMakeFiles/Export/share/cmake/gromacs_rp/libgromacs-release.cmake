#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "libgromacs" for configuration "Release"
set_property(TARGET libgromacs APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(libgromacs PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libgromacs_rp.so.3.5.0"
  IMPORTED_SONAME_RELEASE "libgromacs_rp.so.3"
  )

list(APPEND _IMPORT_CHECK_TARGETS libgromacs )
list(APPEND _IMPORT_CHECK_FILES_FOR_libgromacs "${_IMPORT_PREFIX}/lib64/libgromacs_rp.so.3.5.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
