#----------------------------------------------------------------
# Generated CMake target import file for configuration "None".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "pigsl.e" for configuration "None"
set_property(TARGET pigsl.e APPEND PROPERTY IMPORTED_CONFIGURATIONS NONE)
set_target_properties(pigsl.e PROPERTIES
  IMPORTED_LOCATION_NONE "${_IMPORT_PREFIX}/bin/pigsl.e"
  )

list(APPEND _IMPORT_CHECK_TARGETS pigsl.e )
list(APPEND _IMPORT_CHECK_FILES_FOR_pigsl.e "${_IMPORT_PREFIX}/bin/pigsl.e" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
