cmake_minimum_required(VERSION 3.2)
project(KronMat)


# Additional 'find.cmakes'
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Modules)


# Google test framework
include("test/gtest.cmake")

# Thing to get BLAS found
find_package(BLAS REQUIRED)
find_package(CBLAS REQUIRED)
find_path(BLAS_INCLUDE_DIRS cblas.h
         HINTS
        "/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers/"
        "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers/")

# Find Boost
set(Boost_MULTITHREADED ON)
find_package(Boost COMPONENTS timer system REQUIRED)


# Compiler flags - if GNU add in openmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -msse4.2")
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wa,-q")
endif()


# Includes
include_directories(${BLAS_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GTEST_INCLUDE_DIR})

# Source files for executable
set(SOURCE_FILES blas_wrap.h kronvecstack_dot_vec.h basic_kronmat.h basic_kronmat.cpp play.cpp
                    add_full_kron.h sq_dist.h kernel.h dataset.h distances.h parameter.h)
add_executable(KronMat ${SOURCE_FILES})
target_link_libraries(KronMat ${CBLAS_LIBRARIES} ${Boost_LIBRARIES} ${GTEST_BOTH_LIBRARIES})

# Tests
enable_testing()
add_executable(matrix_test test/testmatrix.cpp)
target_link_libraries(matrix_test ${CBLAS_LIBRARIES} ${Boost_LIBRARIES} ${GTEST_LIBRARY})
add_test(matrix_test matrix_test)