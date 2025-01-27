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

OPTION(O_DOXYGEN "Include support for doxygen"            ON )
OPTION(O_GTEST   "Include support for Google Test"        ON )

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

if(O_GTEST)
  find_package(GTest)
endif()

find_file(MRMathCPP_PATH "MRMathCPPLib.cmake" PATHS "MRMathCPP/build"  "../MRMathCPP/build")
if(NOT MRMathCPP_PATH STREQUAL "MRMathCPP_PATH-NOTFOUND")
  message(STATUS "Found MRMathCPP: ${MRMathCPP_PATH}")
  include("${MRMathCPP_PATH}")
else()
  message(FATAL_ERROR " MRMathCPP Search Failed!\n"
                      "     The MRMathCPP repository must be located in the same directory as this one, and\n"
                      "     must be configured -- i.e. a generated MRMathCPPLib.cmake file must exist in the.\n"
                      "     build directory.\n")
endif()
get_target_property(MRMathCPP_INCLUDE MRMathCPP INTERFACE_INCLUDE_DIRECTORIES)

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
# Add targets for documentation

if(Doxygen_FOUND)
  foreach(DOXINPUT IN ITEMS "lib" "examples")
    message(STATUS "Generateing doxygen target for ${DOXINPUT}")
    CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/Doxyfile.cmake ${CMAKE_BINARY_DIR}/doc-${DOXINPUT}/Doxyfile)
    add_custom_target(doc-${DOXINPUT}
      COMMAND ${DOXYGEN_EXECUTABLE} > dox.out
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/doc-${DOXINPUT}
      COMMENT "Generating ${DOXINPUT} documentation with Doxygen"
      VERBATIM)
  endforeach(DOXINPUT)
else()
  message(STATUS "Warning: Doxygen not found.  No documentation targets!")
endif()

######################################################################################################################################################
# Create Targets for examples, functional tests, & unit tests

# Master list of targets
set(TARGETS "hello_world_adaptive"
            "hello_world_regular"
            ### Unit Tests
            "tree_basics_7b1" "tree_basics_7b2" "tree_basics_7b3" "tree_basics_7b4" "tree_basics_7b5"
            "tree_basics_15b1" "tree_basics_15b3"
            "tree_corners"
            "tree_children" "tree_neighbors"
            ### Functional Tests
            # None yet
           )

# Add a target for each one we can build
foreach(CURTGT IN LISTS TARGETS)
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/${CURTGT}.cpp")
    add_executable(${CURTGT} EXCLUDE_FROM_ALL "examples/${CURTGT}.cpp")
    list(APPEND TARGETS_EXAMPLES ${CURTGT})
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/ftests/${CURTGT}.cpp")
    add_executable(${CURTGT} EXCLUDE_FROM_ALL "ftests/${CURTGT}.cpp")
    list(APPEND TARGETS_FTEST ${CURTGT})
  elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/utests/${CURTGT}.cpp")
    if(NOT GTest_FOUND)
      message(NOTICE "Target ${CURTGT} will not be built (no GTest support)")
      continue()
    endif()
    add_executable(${CURTGT} EXCLUDE_FROM_ALL "utests/${CURTGT}.cpp")
    list(APPEND TARGETS_UTEST ${CURTGT})
  else()
    message("Warning: Unable to find source for target ${CURTGT}!")
    continue()
  endif()

  target_compile_features(${CURTGT} PUBLIC cxx_std_23)
  target_link_libraries(${CURTGT} PRIVATE MRPTree)
  # target_include_directories(${CURTGT} PRIVATE "${MRMathCPP_INCLUDE}")
  # target_include_directories(${CURTGT} PRIVATE "${MRPtree_INCLUDE}")

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

  # Unit Test Specific Stuff
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/utests/${CURTGT}.cpp")
    target_link_libraries(${CURTGT} PRIVATE GTest::gtest_main)
    add_test(NAME "${CURTGT}" COMMAND "${CURTGT}")
  endif()

endforeach(CURTGT)

######################################################################################################################################################

add_custom_target(ftests
  DEPENDS ${TARGETS_FTEST}
  COMMENT "Building Functional Tests"
)

add_custom_target(run-ftests
  DEPENDS ${TARGETS_FTEST}
  COMMENT "Running Functional Tests"
  COMMAND "${CMAKE_CURRENT_SOURCE_DIR}/ftests/run_tests.sh"
)

if(GTest_FOUND)
  add_custom_target(utests
    DEPENDS ${TARGETS_UTEST}
    COMMENT "Building Unit Tests"
  )
  add_custom_target(run-utests
    DEPENDS ${TARGETS_UTEST}
    COMMENT "Running Unit Tests"
    COMMAND "ctest"
  )
endif()
 
add_custom_target(examples
  DEPENDS ${TARGETS_EXAMPLES}
  COMMENT "Building Examples"
)

add_custom_target(clean-data
  COMMAND  rm -f *.vtu *.vtp *.vtk *.ply
  COMMENT "Cleaning VTK & PLY files"
)

if(EXISTS "c:/Users/richmit/Documents/WWW/site/SS/MRPTree/doc-lib/")
  add_custom_target(install-doc-lib
    DEPENDS doc-lib
    COMMAND  rsync -rlt --log-format=%f --stats --delete --delete-excluded --modify-window=2 doc-lib/ /c/Users/richmit/Documents/WWW/site/SS/MRPTree/doc-lib
    COMMENT "Put doc-lib data in web site directory"
  )
endif()

if(EXISTS "c:/Users/richmit/Documents/WWW/site/SS/MRPTree/doc-examples/")
  add_custom_target(install-doc-examples
    DEPENDS doc-examples/
    COMMAND  rsync -rlt --log-format=%f --stats --delete --delete-excluded --modify-window=2 doc-examples/ /c/Users/richmit/Documents/WWW/site/SS/MRPTree/doc-examples
    COMMENT "Put doc-examples data in web site directory"
  )
endif()


