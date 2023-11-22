cmake_minimum_required(VERSION 3.14...3.19)

project(patchyAnalysisTools-user VERSION 2.0.0)

find_package(patchyAnalysisTools REQUIRED)


set(TEST_LIST
    test1
    test2
    )

foreach (CUR_TEST ${TEST_LIST})
    add_executable(${CUR_TEST} "${PROJECT_SOURCE_DIR}/${CUR_TEST}.cpp")
    target_link_libraries(${CUR_TEST} PRIVATE patchyAnalysisTools)
endforeach (CUR_TEST)

include(CTest)

foreach (CUR_TEST ${TEST_LIST})
    add_test(NAME ${CUR_TEST} COMMAND ${CUR_TEST})
endforeach (CUR_TEST)