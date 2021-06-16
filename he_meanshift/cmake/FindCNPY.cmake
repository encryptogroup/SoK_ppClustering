include(FindPackageHandleStandardArgs)

find_package(ZLIB REQUIRED)

find_path(
  CNPY_HEADER
  NAMES cnpy.h
  PATH_SUFFIXES include)

find_library(CNPY_LIB libcnpy.a PATH_SUFFIXES lib)

find_package_handle_standard_args(CNPY DEFAULT_MSG CNPY_HEADER CNPY_LIB)

if(CNPY_HEADER AND CNPY_LIB)
  add_library(CNPY STATIC IMPORTED)
  target_link_libraries(CNPY INTERFACE ${ZLIB_LIBRARIES})
  set_target_properties(
    CNPY PROPERTIES IMPORTED_LOCATION ${CNPY_LIB} INTERFACE_INCLUDE_DIRECTORIES
                                                  ${CNPY_HEADER})
endif(CNPY_HEADER AND CNPY_LIB)
