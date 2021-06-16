add_library(
  HEAAN
  src/BootContext.cpp
  src/Ciphertext.cpp
  src/EvaluatorUtils.cpp
  src/Key.cpp
  src/Plaintext.cpp
  src/Ring.cpp
  src/RingMultiplier.cpp
  src/Scheme.cpp
  src/SchemeAlgo.cpp
  src/SecretKey.cpp
  src/SerializationUtils.cpp
  src/StringUtils.cpp
  src/TestScheme.cpp
  src/TimeUtils.cpp)

target_compile_features(HEAAN PRIVATE cxx_std_11)
target_compile_options(HEAAN PRIVATE "-O3")

target_include_directories(HEAAN PUBLIC src)
target_link_libraries(HEAAN PUBLIC NTL)
target_link_libraries(HEAAN PRIVATE Threads::Threads GMP)

add_executable(heaan_test test.cpp)
target_link_libraries(heaan_test NTL HEAAN)
set_target_properties(heaan_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                            ${PROJECT_BINARY_DIR}/tests)
