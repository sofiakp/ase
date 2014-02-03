cmake_minimum_required (VERSION 2.6)
cmake_policy(VERSION 2.6)

project (my_project)

# Set default build type
set (CMAKE_BUILD_TYPE, Release)

set (CMAKE_CURRENT_BINARY_DIR, builds)

# ==================== UTILITY FUNCTIONS ==================== 

#add_custom_target(core)
#add_custom_target(test)
#
#function(add_core_executable PROG)
        #add_executable(${PROG} ${ARGN})
        #add_dependencies(core ${PROG})
#endfunction(add_core_executable PROG)
#
#function(add_test_executable PROG)
        #add_executable(${PROG} ${ARGN})
        #add_dependencies(test ${PROG})
#endfunction(add_test_executable PROG)

# ==================== FLAGS ==================== 

set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D_FILE_OFFSET_BITS=64 -Wall -Wextra -Wno-sign-compare")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_C_FLAGS} -Wno-deprecated")

set (CMAKE_C_FLAGS_RELEASE "-O2")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}")

set (CMAKE_C_FLAGS_DEBUG "-g -DDEBUG=1")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")

# ==================== CHECKS ==================== 

# verify that we have openmp
#INCLUDE(CheckCCompilerFlag)
# Make sure the C compile can take -fopenmp as a flag
#CHECK_C_COMPILER_FLAG(-fopenmp HAVE_OPENMP)

# ==================== PROJECT DIRECTORIES ==================== 

add_subdirectory(deps/yaml-cpp)
add_subdirectory(src)
add_subdirectory(test)
