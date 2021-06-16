include(FindPackageHandleStandardArgs)

find_path(
  GMP_HEADERS
  NAMES gmp.h
  PATHS ${GMP_DIR} ${INCLUDE_INSTALL_DIR})

find_library(GMP_LIB gmp PATHS ${GMP_DIR} ${LIB_INSTALL_DIR})

find_package_handle_standard_args(GMP DEFAULT_MSG GMP_HEADERS GMP_LIB)

if(GMP_FOUND)
  add_library(GMP STATIC IMPORTED)
  set_target_properties(
    GMP PROPERTIES IMPORTED_LOCATION ${GMP_LIB} INTERFACE_INCLUDE_DIRECTORIES
                                                ${GMP_HEADERS})
endif(GMP_FOUND)
