cmake_minimum_required(VERSION 3.15...3.26)

set(CMAKE_CXX_STANDARD 20)

project(
    ${SKBUILD_PROJECT_NAME}
    VERSION ${SKBUILD_PROJECT_VERSION}
    LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

file (GLOB SOURCE_FILES "src/*.cpp")
file (GLOB HEADER_FILES "src/include/*.hpp")

python_add_library(patchyAnalysisTools_cpp MODULE ${SOURCE_FILES} ${HEADER_FILES}  WITH_SOABI)

target_link_libraries(patchyAnalysisTools_cpp PRIVATE pybind11::headers )
target_compile_definitions(patchyAnalysisTools_cpp PRIVATE VERSION_INFO=${PROJECT_VERSION})

install(TARGETS patchyAnalysisTools_cpp DESTINATION patchyAnalysisTools)