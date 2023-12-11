cmake_minimum_required(VERSION 3.14...3.19)

project(patchyAnalysisTools-user VERSION 2.0.0)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(patchyAnalysisTools REQUIRED)

include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.11.0
)
FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)

set(TEST_LIST
    test_rotation
    test_frame
    test_patches
    test_clusters
    test_grid
    )

foreach (CUR_TEST ${TEST_LIST})
    add_executable(${CUR_TEST} "${PROJECT_SOURCE_DIR}/${CUR_TEST}.cpp")
    
    target_link_libraries(${CUR_TEST}
                          PRIVATE
                            GTest::GTest
                            patchyAnalysisTools
                         )

endforeach (CUR_TEST)

foreach (CUR_TEST ${TEST_LIST})
    add_test(NAME ${CUR_TEST} COMMAND ${CUR_TEST})
endforeach (CUR_TEST)
