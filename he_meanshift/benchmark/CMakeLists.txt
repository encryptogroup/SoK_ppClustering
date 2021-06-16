set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/benchmarks)

add_executable(cluster cluster.cpp utils.cpp)
target_link_libraries(
  cluster
  CKP19
  HEAAN
  NTL
  nlohmann_json::nlohmann_json
  Boost::program_options
  CNPY)

add_executable(plaintext_cluster plaintext_cluster.cpp utils.cpp)
target_link_libraries(
  plaintext_cluster
  CKP19
  HEAAN
  NTL
  nlohmann_json::nlohmann_json
  Boost::program_options
  CNPY)

add_custom_target(benchmarks)
add_dependencies(benchmarks cluster plaintext_cluster)
