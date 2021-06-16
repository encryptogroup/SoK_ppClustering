
find_package(ABY REQUIRED)
find_package(GMP REQUIRED)
find_package(Threads REQUIRED)

add_library(mpo19_test_interface INTERFACE)
target_include_directories(mpo19_test_interface INTERFACE
    "${PROJECT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_LIST_DIR}")
target_compile_options(mpo19_test_interface INTERFACE "-Wall" "-Wextra" "-g" "-fno-strict-aliasing")
target_link_libraries(mpo19_test_interface INTERFACE Boost::boost)

add_executable(utility_test utility_test.cpp)
target_link_libraries(utility_test mpo19_test_interface)

add_executable(random_test random_test.cpp)
target_link_libraries(random_test paillier_wrapper random mpo19_test_interface)
target_compile_options(random_test PRIVATE "-O2")

add_executable(point_test point_test.cpp)
target_link_libraries(point_test mpo19_test_interface)

add_executable(paillier_wrapper_test paillier_wrapper_test.cpp)
target_link_libraries(paillier_wrapper_test paillier_wrapper random mpo19_test_interface)

add_executable(channel_test channel_test.cpp)
target_link_libraries(channel_test paillier_wrapper channel Threads::Threads mpo19_test_interface)

add_executable(linkage_test linkage_test.cpp)
target_link_libraries(linkage_test linkage paillier_wrapper random python mpo19_test_interface)

add_executable(cluster_test cluster_test.cpp)
target_link_libraries(cluster_test cluster mpo19_test_interface)

add_executable(clustering_test clustering_test.cpp)
target_link_libraries(clustering_test cluster python mpo19_test_interface)

add_executable(pca_test pca_test.cpp)
target_link_libraries(pca_test pca mpo19_test_interface python)
target_compile_definitions(pca_test PUBLIC "MPO19_TESTING_PCA_")
#target_compile_options(random_test PRIVATE "-O2")

add_test(NAME random_test COMMAND random_test)
add_test(NAME point_test COMMAND point_test)
add_test(NAME paillier_wrapper_test COMMAND paillier_wrapper_test)
add_test(NAME channel_test COMMAND channel_test)
add_test(NAME linkage_test COMMAND linkage_test)
add_test(NAME utility_test COMMAND utility_test)
add_test(NAME clustering_test COMMAND clustering_test)
add_test(NAME pca_test COMMAND pca_test)
