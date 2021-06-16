cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(ckp19 VERSION 0.1)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Path to custom cmake files for external libs
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set(THREADS_PREFER_PTHREAD_FLAG ON)

# Find external dependencies
find_package(GMP REQUIRED)
find_package(NTL 11.0.0 REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost 1.72.0 COMPONENTS program_options unit_test_framework system
                                     REQUIRED)
find_package(CNPY REQUIRED)
find_package(nlohmann_json 3 REQUIRED)

# Add HEAAN
add_subdirectory(extern/HEAAN)

# Add project local directories
add_subdirectory(src)
add_subdirectory(benchmark)
add_subdirectory(test)
