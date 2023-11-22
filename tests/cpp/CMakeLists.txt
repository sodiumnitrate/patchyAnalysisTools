cmake_minimum_required(VERSION 3.14...3.19)

project(patchyAnalysisTools-user VERSION 2.0.0)

find_package(patchyAnalysisTools REQUIRED)

add_executable(patchyAnalysisTools-user "${PROJECT_SOURCE_DIR}/test.cpp")
target_link_libraries(patchyAnalysisTools-user PRIVATE patchyAnalysisTools)

include(CTest)

add_test(NAME test1 COMMAND patchyAnalysisTools-user)