include(FindPackageHandleStandardArgs)

find_path(
  CNPY_HEADER
  NAMES cnpy.h
  PATH_SUFFIXES include)

find_library(CNPY_LIB cnpy PATH_SUFFIXES lib)

find_package_handle_standard_args(CNPY DEFAULT_MSG CNPY_HEADER CNPY_LIB)

if(CNPY_HEADER AND CNPY_LIB)
  add_library(CNPY STATIC IMPORTED)
  set_target_properties(
    CNPY PROPERTIES IMPORTED_LOCATION ${CNPY_LIB} INTERFACE_INCLUDE_DIRECTORIES
                                                  ${CNPY_HEADER})
endif(CNPY_HEADER AND CNPY_LIB)
