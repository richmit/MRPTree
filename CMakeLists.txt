######################################################################################################################################################
######################################################################################################################################################
cmake_minimum_required(VERSION 3.30.0 FATAL_ERROR)

######################################################################################################################################################
# Proejct meta data

project(MRPTree 
  VERSION 0.5.0.0
  DESCRIPTION "Generalized bitree/quadtree/octree library"
  HOMEPAGE_URL "https://github.com/richmit/MRPTree/"
  LANGUAGES CXX)

set(MRPTREE_VERSION_DATE  2024-09-16)
mark_as_advanced(MRPTREE_VERSION_DATE)

######################################################################################################################################################
# Include ctest
include(CTest) 

######################################################################################################################################################
# Here are the options we use to modify the build -- these options are passed to cmake like so: cmake -DO_GTEST=YES or -DO_GTEST=NO

OPTION(O_DOXYGEN     "Include support for doxygen"             ON  )
OPTION(O_GTEST       "Include support for Google Test"         ON  )
OPTION(OWNLOAD_GTEST "Download latest Google Test from Github" OFF )

######################################################################################################################################################
# Find tools & libraries used to build targets

if(O_DOXYGEN)
  set(MRPTREE_OPT_DOXYGEN "1")
  find_package(Doxygen)
  if(Doxygen_FOUND)
    set(MRPTREE_FOUND_DOXYGEN "1")
  else()
    set(MRPTREE_FOUND_DOXYGEN "0")
  endif()
else()
  set(MRPTREE_OPT_DOXYGEN "0")
endif()

# CODE GEN: ~/MJR/core/codeBits/bin/src2orgListing ~/core/codeBits/cmake-template/include_gtest.cmake
if(O_GTEST)
  if(OWNLOAD_GTEST)
    message(STATUS "GTest: Downloading")
    include(FetchContent)
    FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/refs/heads/main.zip
      DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
    set(GTest_FOUND "YES") # If FetchContent_MakeAvailable didn't work, then cmake will error out.
  else()
    message(STATUS "GTest: Not downloading GTest")
    find_package(GTest)
  endif()
else()
  if(OWNLOAD_GTEST)
    message(FATAL_ERROR "-DO_GTEST must be YES if -DOWNLOAD_GTEST is YES!")
  endif()
endif()

# CODE GEN: ~/MJR/core/codeBits/bin/src2orgListing ../MRMathCPP/IncludeMRMathCPPLib.cmake 
find_file(MRMathCPP_PATH "MRMathCPPLib.cmake" PATHS "MRMathCPP/build"  "../MRMathCPP/build")
if(NOT MRMathCPP_PATH STREQUAL "MRMathCPP_PATH-NOTFOUND")
  message(STATUS "Found MRMathCPP: ${MRMathCPP_PATH}")
  include("${MRMathCPP_PATH}")
  get_target_property(MRMathCPP_INCLUDE MRMathCPP INTERFACE_INCLUDE_DIRECTORIES)
else()
  message(FATAL_ERROR " MRMathCPP Search Failed!\n"
                      "     The MRMathCPP repository must be located in the same directory as this repository,\n"
                      "     and must be configured with a build directory named 'build' at the root of the \n"
                      "     MRMathCPP repository.  That is to say, do the following in the same directory \n"
                      "     where you cloned this repository: \n"
                      "        git clone 'https://github.com/richmit/MRMathCPP.git'\n"
                      "        cd MRMathCPP/build                                  \n"
                      "        cmake ..                                            \n"
                      "        cd ../..                                            \n"
                      "     Then return to this repository, and try to configure it again.\n")
endif()

# CODE GEN: ~/MJR/core/codeBits/bin/src2orgListing ../mraster/IncludeMRasterLib.cmake | sed 's/FATAL_ERROR/WARNING/'
find_file(MRaster_PATH "MRasterLib.cmake" PATHS "mraster/build"  "../mraster/build")
if(NOT MRaster_PATH STREQUAL "MRaster_PATH-NOTFOUND")
  message(STATUS "Found MRaster: ${MRaster_PATH}")
  include("${MRaster_PATH}")
  get_target_property(MRaster_INCLUDE MRaster INTERFACE_INCLUDE_DIRECTORIES)
else()
  message(WARNING " MRaster Search Failed!\n"
                      "     The MRaster repository must be located in the same directory as this repository,\n"
                      "     and must be configured with a build directory named 'build' at the root of the \n"
                      "     MRaster repository.  That is to say, do the following in the same directory \n"
                      "     where you cloned this repository: \n"
                      "        git clone 'https://github.com/richmit/mraster.git' \n"
                      "        cd mraster/build                                   \n"
                      "        cmake ..                                           \n"
                      "        cd ../..                                           \n"
                      "     Then return to this repository, and try to configure it again.\n")
endif()

######################################################################################################################################################
# Create interface target for the entire project

set(MRPTREE_INCLUDES "lib/MR_rect_tree.hpp")
add_library(MRPTree INTERFACE ${MRPTREE_INCLUDES})
target_include_directories(MRPTree INTERFACE ${MRMathCPP_INCLUDE})
target_include_directories(MRPTree INTERFACE "${PROJECT_SOURCE_DIR}/lib")

######################################################################################################################################################
# Export the library interface cmake file
export(TARGETS MRPTree FILE "${CMAKE_BINARY_DIR}/MRPTreeLib.cmake")

######################################################################################################################################################
# Create variable for project includes
get_target_property(MRPTree_INCLUDE MRPTree INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "Combined Project Includes: ${MRPTree_INCLUDE}")

######################################################################################################################################################
# Create Targets for examples, functional tests, & unit tests

# CODE GEN: echo 'set(TARGETS_REQ_CELL '$(basename -s.cpp $(grep -El '#include "(MR_cell_cplx.hpp|MR_rt_to_cc.hpp)"' */*.cpp))')'
set(TARGETS_REQ_CELL hello_world_cell hello_world_mraster recipe-surf-plot-adapt recipe-surf-plot-norm recipe-surf-plot-rs-quad recipe-surf-plot-rs-tri complex_magnitude_surface curve_plot ear_surface ear_surface_glue holy_wave_surf implicit_curve_2d implicit_surface parametric_curve_3d parametric_surface_with_defects performance_with_large_surface surface_branch_glue surface_plot_annular_edge surface_plot_corner surface_plot_edge surface_plot_step surface_with_normals trefoil vector_field_3d flat_test_tree_01 nan_solver rect_fix_dup rect_fix_nan segment_folder triangle_folder check_cell_hexahedron check_cell_pyramid check_cell_quad check_cell_segment check_cell_triangle geomi_pnt_line_distance geomi_seg_isect_type geomr_pnt_line_distance geomr_pnt_pln_distance geomr_pnt_tri_distance)

# CODE GEN: echo 'set(TARGETS_REQ_BRIDGE '$(basename -s.cpp $(grep -El '#include "(ramCanvas.hpp|MRcolor.hpp)"' */*.cpp))')'
set(TARGETS_REQ_BRIDGE hello_world_mraster complex_color_image complex_magnitude_surface test_interp_scale)

# CODE GEN: echo 'set(TARGETS_REQ_TREE '$(basename -s.cpp $(grep -El '#include "(MR_rect_tree.hpp)"' */*.cpp))')'
set(TARGETS_REQ_TREE hello_world_cell hello_world_mraster hello_world_tree_adaptive hello_world_tree_regular recipe-surf-plot-adapt recipe-surf-plot-norm recipe-surf-plot-rs-quad recipe-surf-plot-rs-tri complex_magnitude_surface curve_plot ear_surface ear_surface_glue implicit_curve_2d implicit_surface parametric_curve_3d parametric_surface_with_defects performance_with_large_surface surface_branch_glue surface_plot_annular_edge surface_plot_corner surface_plot_edge surface_plot_step surface_with_normals trefoil vector_field_3d flat_test_tree_01 nan_solver rect_fix_dup rect_fix_nan segment_folder triangle_folder tree_basics_15b1 tree_basics_15b3 tree_basics_7b1 tree_basics_7b2 tree_basics_7b3 tree_basics_7b4 tree_basics_7b5 tree_children tree_corners tree_neighbors)

# CODE GEN: echo 'set(TARGETS_REQ_MRASTER '$(basename -s.cpp $(grep -El '#include "(ramCanvas.hpp|MRcolor.hpp)"' */*.cpp))')'
set(TARGETS_REQ_MRASTER hello_world_mraster complex_color_image complex_magnitude_surface test_interp_scale)

# CODE GEN: echo 'set(TARGETS_REQ_MRASTER '$(basename -s.cpp $(grep -El '#include <gtest/gtest.h>' */*.cpp))')'
set(TARGETS_REQ_GTEST check_cell_hexahedron check_cell_pyramid check_cell_quad check_cell_segment check_cell_triangle geomi_pnt_line_distance geomi_seg_isect_type geomr_pnt_line_distance geomr_pnt_pln_distance geomr_pnt_tri_distance tree_basics_15b1 tree_basics_15b3 tree_basics_7b1 tree_basics_7b2 tree_basics_7b3 tree_basics_7b4 tree_basics_7b5 tree_children tree_corners tree_neighbors)

# Construct list of targets we can build
set(COMBINED_TARGETS ${TARGETS_REQ_CELL} ${TARGETS_REQ_BRIDGE} ${TARGETS_REQ_TREE} ${TARGETS_REQ_MRASTER} ${TARGETS_REQ_GTEST})

if((MRaster_PATH STREQUAL "MRaster_PATH-NOTFOUND") AND TARGETS_REQ_MRASTER)
  list(REMOVE_ITEM COMBINED_TARGETS ${TARGETS_REQ_MRASTER})
  message(WARNING "Example programs and/or unit tests requiring MRaster will NOT be built: ${TARGETS_REQ_MRASTER}!")
endif()

if((NOT GTest_FOUND) AND TARGETS_REQ_GTEST)
  list(REMOVE_ITEM COMBINED_TARGETS ${TARGETS_REQ_GTEST})
  message(WARNING "Unit tests requiring GTest will NOT be built: ${TARGETS_REQ_GTEST}!")
endif()

list(REMOVE_DUPLICATES COMBINED_TARGETS)

# Add a target for each one we can build
foreach(CURTGT IN LISTS COMBINED_TARGETS)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/${CURTGT}.cpp")
    add_executable(${CURTGT} EXCLUDE_FROM_ALL "examples/${CURTGT}.cpp")
    list(APPEND TARGETS_EXAMPLES_MIN ${CURTGT})
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples-func-viz/${CURTGT}.cpp")
    add_executable(${CURTGT} EXCLUDE_FROM_ALL "examples-func-viz/${CURTGT}.cpp")
    list(APPEND TARGETS_EXAMPLES_VIZ ${CURTGT})
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples-cookbook/${CURTGT}.cpp")
    add_executable(${CURTGT} EXCLUDE_FROM_ALL "examples-cookbook/${CURTGT}.cpp")
    list(APPEND TARGETS_EXAMPLES_CB ${CURTGT})
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ftests/${CURTGT}.cpp")
    add_executable(${CURTGT} EXCLUDE_FROM_ALL "ftests/${CURTGT}.cpp")
    list(APPEND TARGETS_FTEST ${CURTGT})
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/utests/${CURTGT}.cpp")
    add_executable(${CURTGT} EXCLUDE_FROM_ALL "utests/${CURTGT}.cpp")
    list(APPEND TARGETS_UTEST ${CURTGT})
  else()
    message("Warning: Unable to find source for target ${CURTGT}!")
    continue()
  endif()

  target_compile_features(${CURTGT} PUBLIC cxx_std_23)

  # Add build requirements for CELL & BRIDGE
  list (FIND TARGETS_REQ_CELL   ${CURTGT} TMPR1)
  list (FIND TARGETS_REQ_BRIDGE ${CURTGT} TMPR2)
  if((${TMPR1} GREATER -1) OR (${TMPR2} GREATER -1))
    target_include_directories(${CURTGT} PRIVATE "${PROJECT_SOURCE_DIR}/examples-lib3d")
    target_link_libraries(${CURTGT} PRIVATE MRMathCPP)
  endif()

  # Add build requirements for MRPTree
  list (FIND TARGETS_REQ_TREE ${CURTGT} TMPR1)
  if(${TMPR1} GREATER -1)
    target_link_libraries(${CURTGT} PRIVATE MRPTree)
  endif()

  # Add build requirements for MRaster
  list (FIND TARGETS_REQ_MRASTER ${CURTGT} TMPR1)
  if(${TMPR1} GREATER -1)
    target_link_libraries(${CURTGT} PRIVATE MRaster)
  endif()

  # Add build requirements for GTest
  list (FIND TARGETS_REQ_GTEST ${CURTGT} TMPR1)
  if(${TMPR1} GREATER -1)
    target_link_libraries(${CURTGT} PRIVATE GTest::gtest_main)
    add_test(NAME "${CURTGT}" COMMAND "${CURTGT}")
  endif()

  # Compiler specific stuff
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${CURTGT} PUBLIC -Wall -Wconversion -Wno-unknown-pragmas -Wextra -Wno-deprecated-copy)
    target_compile_options(${CURTGT} PUBLIC -O4)
  elseif((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") OR (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang"))
    target_compile_options(${CURTGT} PUBLIC -Wall -Wconversion -Wno-unknown-pragmas -Wextra -Wno-sign-conversion)
    target_compile_options(${CURTGT} PUBLIC -O3)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    target_compile_options(${CURTGT} PUBLIC -Wall)
    target_compile_options(${CURTGT} PUBLIC -O3)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    #message("Warning: MSVC support is currently experimental")
  endif()

endforeach(CURTGT)

######################################################################################################################################################

add_custom_target(clean-data
  COMMAND  rm -f *.vtu *.vtp *.vtk *.ply
  COMMENT "Cleaning VTK & PLY files"
)

if(TARGETS_FTEST)
  add_custom_target(ftests
    DEPENDS ${TARGETS_FTEST}
    COMMENT "Building Functional Tests"
  )
  add_custom_target(run-ftests
    DEPENDS ${TARGETS_FTEST}
    COMMENT "Running Functional Tests"
    COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/ftests/run_tests.sh"
  )
endif()

if(TARGETS_EXAMPLES_MIN OR TARGETS_EXAMPLES_CB OR TARGETS_EXAMPLES_VIZ)
  add_custom_target(examples-all
    DEPENDS ${TARGETS_EXAMPLES_MIN} ${TARGETS_EXAMPLES_CB} ${TARGETS_EXAMPLES_VIZ}
    COMMENT "Building Minimal MRPtree Examples"
  )
endif()

if(TARGETS_EXAMPLES_MIN)
  add_custom_target(examples-min
    DEPENDS ${TARGETS_EXAMPLES_MIN}
    COMMENT "Building Minimal MRPtree Examples"
  )
endif()

if(TARGETS_EXAMPLES_CB)
  add_custom_target(examples-cb
    DEPENDS ${TARGETS_EXAMPLES_CB}
    COMMENT "Building Cookbook Examples"
  )
endif()

if(TARGETS_EXAMPLES_VIZ)
  add_custom_target(examples-viz
    DEPENDS ${TARGETS_EXAMPLES_VIZ}
    COMMENT "Building Func Viz Examples"
  )
endif()

if(TARGETS_UTEST)
  add_custom_target(utests
    DEPENDS ${TARGETS_UTEST}
    COMMENT "Building Unit Tests"
  )
  add_custom_target(run-utests
    DEPENDS ${TARGETS_UTEST}
    COMMENT "Running Unit Tests"
    COMMAND cmake --build . -t test
  )
endif()

# CODE GEN: ~/MJR/core/codeBits/bin/src2orgListing ~/core/codeBits/cmake-template/include_doxygen-targets.cmake | sed 's/PROJECT_NAME/MRPTree/g; s/SRCP/"lib" "examples-lib3d" "examples" "examples-func-viz" "examples-cookbook"/'
if(Doxygen_FOUND)
  foreach(DOXINPUT IN ITEMS "lib" "examples-lib3d" "examples" "examples-func-viz" "examples-cookbook")
    message(STATUS "Generateing doxygen generation target for ${DOXINPUT}")
    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/Doxyfile.cmake ${CMAKE_BINARY_DIR}/doc-${DOXINPUT}/Doxyfile)
    set(DOXSRCF "doc-${DOXINPUT}")
    file(GLOB TMP LIST_DIRECTORIES false "${CMAKE_SOURCE_DIR}/${DOXINPUT}/[a-zA-Z0-9]*.?pp")
    list(APPEND DOXSRCF ${TMP})
    file(GLOB TMP LIST_DIRECTORIES false "${CMAKE_SOURCE_DIR}/${DOXINPUT}/[a-zA-Z0-9]*.f90")
    list(APPEND DOXSRCF ${TMP})
    list(APPEND TARGETS_DOX_GENERATE doc-${DOXINPUT})
    add_custom_target(doc-${DOXINPUT}
      DEPENDS ${DOXSRCF}
      COMMAND ${DOXYGEN_EXECUTABLE} > dox.out
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc-${DOXINPUT}
      COMMENT "Generating ${DOXINPUT} documentation with Doxygen"
      VERBATIM)
    if(EXISTS "c:/Users/richmit/MJR/WWW/site/SS/MRPTree/doc-${DOXINPUT}/")
      message(STATUS "Generateing doxygen publish target for ${DOXINPUT}")
      list(APPEND TARGETS_DOX_INSTALL "install-doc-${DOXINPUT}")
      add_custom_target("install-doc-${DOXINPUT}"
        DEPENDS ${DOXSRCF}
        COMMAND  rsync -rlt --log-format=%f --stats --delete --delete-excluded --modify-window=2 "doc-${DOXINPUT}/" "/c/Users/richmit/MJR/WWW/site/SS/MRPTree/doc-${DOXINPUT}"
        COMMENT "Put ${DOXINPUT} data in web site directory")
    endif()
  endforeach(DOXINPUT)
  if(TARGETS_DOX_GENERATE)
    add_custom_target(doc-all
      DEPENDS ${TARGETS_DOX_GENERATE}
      COMMENT "Generateing all doxygen documentation"
    )
  endif()
  if(TARGETS_DOX_INSTALL)
    add_custom_target(install-doc-all
      DEPENDS ${TARGETS_DOX_INSTALL}
      COMMENT "Running all doxygen web site targets"
    )
  endif()
else()
  message(STATUS "Warning: Doxygen not found.  No documentation targets!")
endif()
