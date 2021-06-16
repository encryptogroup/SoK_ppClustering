find_library(
  PAILLIER_LIB
  NAMES paillier libpaillier
  PATH_SUFFIXES lib)

find_path(
  PAILLIER_HEADER
  NAMES paillier.h
  PATH_SUFFIXES include)

if(PAILLIER_LIB AND PAILLIER_HEADER)
  message(STATUS "Found Paillier")
  add_library(paillier STATIC IMPORTED)
  set_target_properties(
    paillier PROPERTIES IMPORTED_LOCATION ${PAILLIER_LIB}
                        INTERFACE_INCLUDE_DIRECTORIES ${PAILLIER_HEADER})
else()
  message(SEND_ERROR "Paillier not found")
endif()
