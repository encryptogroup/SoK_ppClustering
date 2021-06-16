add_library(CKP19 ckp19/mean_shift.cpp ckp19/plaintext_clustering.cpp
                  ckp19/dataset.cpp)

target_include_directories(CKP19 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(CKP19 PUBLIC HEAAN CNPY)
