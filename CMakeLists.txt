cmake_minimum_required(VERSION 3.15..3.26)

set(CMAKE_CXX_STANDARD 20)

project(
    patchyAnalysisTools
    VERSION "2.0.0"
    LANGUAGES CXX)

include(GNUInstallDirs)
add_subdirectory(eigen)

file (GLOB SOURCE_FILES "${PROJECT_SOURCE_DIR}/src/*.cpp")


# define the C++ library "patchyAnalysisTools"
add_library(patchyAnalysisTools SHARED ${SOURCE_FILES})

target_include_directories(
    patchyAnalysisTools PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                 $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(patchyAnalysisTools Eigen3::Eigen)

install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/"
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})


install(
  TARGETS patchyAnalysisTools
  EXPORT patchyAnalysisToolsTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

set(PATCHYANALYSISTOOLS_CMAKE_PACKAGE_INSTALL_SUBDIR "share/patchyAnalysisTools/cmake")

install(
  EXPORT patchyAnalysisToolsTargets
  DESTINATION ${PATCHYANALYSISTOOLS_CMAKE_PACKAGE_INSTALL_SUBDIR})

include(CMakePackageConfigHelpers)

write_basic_package_version_file(
  patchyAnalysisToolsConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMinorVersion)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/patchyAnalysisToolsConfig.cmake.in" patchyAnalysisToolsConfig.cmake
  INSTALL_DESTINATION ${PATCHYANALYSISTOOLS_CMAKE_PACKAGE_INSTALL_SUBDIR})

install(FILES "${PROJECT_BINARY_DIR}/patchyAnalysisToolsConfig.cmake"
              "${PROJECT_BINARY_DIR}/patchyAnalysisToolsConfigVersion.cmake"
        DESTINATION ${PATCHYANALYSISTOOLS_CMAKE_PACKAGE_INSTALL_SUBDIR})

# if called by scikit-build
if(DEFINED SKBUILD)

  # call pybind11-config to obtain the root of the cmake package
  execute_process(COMMAND ${PYTHON_EXECUTABLE} -m pybind11 --cmakedir
                  OUTPUT_VARIABLE pybind11_ROOT_RAW)
  string(STRIP ${pybind11_ROOT_RAW} pybind11_ROOT)
  find_package(pybind11)

  pybind11_add_module(_patchyAnalysisTools MODULE
                      "${PROJECT_SOURCE_DIR}/src/patchyAnalysisTools/patchyAnalysisTools_py.cpp")

  target_link_libraries(_patchyAnalysisTools PRIVATE patchyAnalysisTools)

  # Installing the extension module to the root of the package
  install(TARGETS _patchyAnalysisTools DESTINATION .)

  configure_file("${PROJECT_SOURCE_DIR}/src/patchyAnalysisTools/__main__.py.in"
                 "${PROJECT_BINARY_DIR}/src/patchyAnalysisTools/__main__.py")

  install(FILES "${PROJECT_BINARY_DIR}/src/patchyAnalysisTools/__main__.py" DESTINATION .)


  if(APPLE)
    set_target_properties(
        _patchyAnalysisTools PROPERTIES INSTALL_RPATH "@loader_path/${CMAKE_INSTALL_LIBDIR}")
  else()
    set_target_properties(_patchyAnalysisTools PROPERTIES INSTALL_RPATH
                                            "$ORIGIN/${CMAKE_INSTALL_LIBDIR}")
  endif()

endif()
