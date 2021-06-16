set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tests)

add_executable(ckp19_test testrunner.cpp mean_shift.cpp)
target_link_libraries(ckp19_test CKP19 HEAAN NTL)

add_custom_target(tests)
add_dependencies(tests ckp19_test)
