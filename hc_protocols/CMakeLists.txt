cmake_minimum_required(VERSION 3.14)
project(mpo19 C CXX)

set(TEST TRUE)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Path to custom cmake files for external libs
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(ABY REQUIRED)
find_package(GMP REQUIRED)
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(Boost REQUIRED COMPONENTS python system thread program_options)
find_package(Threads REQUIRED)
find_package(Paillier REQUIRED)
find_package(RELIC REQUIRED)
find_package(CNPY REQUIRED)
find_package(nlohmann_json 3 REQUIRED)

if(RELIC_FOUND)
  message(STATUS "Found Relic")
  add_library(RELIC::relic STATIC IMPORTED)
  set_target_properties(
    RELIC::relic PROPERTIES IMPORTED_LOCATION ${RLC_LIBRARY}
                            INTERFACE_INCLUDE_DIRECTORIES ${RLC_INCLUDE_DIR})

else()
  message(SEND_ERROR "Relic not found")
endif()

if(NOT ${Python3_VERSION_MINOR} GREATER_EQUAL 6)
  message(ERROR "Python >= 3.6 is required for scikit")
endif()

add_library(mpo19_source_interface INTERFACE)
target_include_directories(
  mpo19_source_interface INTERFACE "${PROJECT_SOURCE_DIR}/include" Boost::boost)
target_compile_options(mpo19_source_interface INTERFACE "-Wall" "-Wextra" "-g"
                                                        "-fno-strict-aliasing")

macro(add_single_source_object_library source_name)
  add_library(${source_name} INTERFACE)
  target_sources(${source_name} INTERFACE src/${source_name}.cpp)
  target_link_libraries(${source_name} INTERFACE mpo19_source_interface)
endmacro(add_single_source_object_library)

add_single_source_object_library(paillier_wrapper)
target_link_libraries(paillier_wrapper INTERFACE paillier GMP::GMP)

add_single_source_object_library(channel)
target_link_libraries(channel INTERFACE Boost::system)

add_single_source_object_library(cluster)
target_link_libraries(cluster INTERFACE Boost::boost)

add_single_source_object_library(linkage)
target_link_libraries(linkage INTERFACE ABY::aby)

add_single_source_object_library(random)

add_single_source_object_library(python)
target_include_directories(python INTERFACE "${Python3_INCLUDE_DIRS}")
target_link_libraries(python INTERFACE Boost::python cluster "${Python3_LIBRARIES}")

add_single_source_object_library(pca)
target_link_libraries(pca INTERFACE paillier_wrapper channel cluster linkage
                                    random)

add_executable(mpo19 src/mpo19.cpp)
target_link_libraries(mpo19 pca python mpo19_source_interface)

add_executable(hc_benchmark benchmarks/hierarchical.cpp)
target_link_libraries(hc_benchmark pca python CNPY nlohmann_json::nlohmann_json Boost::program_options mpo19_source_interface)

if(TEST)
  enable_testing()
  add_subdirectory(tests)
endif()
