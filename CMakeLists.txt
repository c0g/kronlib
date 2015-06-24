cmake_minimum_required(VERSION 3.2)
project(KronMat)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Modules)

find_package(BLAS REQUIRED)
find_package(CBLAS REQUIRED)
find_path(BLAS_INCLUDE_DIRS cblas.h
         HINTS
        "/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers/"
        "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers/")




set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O0 -msse4.2")

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  # using Clang
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wa,-q")
endif()

set(Boost_MULTITHREADED ON)
find_package(Boost COMPONENTS timer system REQUIRED)

include_directories(${BLAS_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

set(SOURCE_FILES blas_wrap.h kronvecstack_dot_vec.h basic_kronmat.h basic_kronmat.cpp play.cpp)
add_executable(KronMat ${SOURCE_FILES} add_full_kron.h sq_dist.h kernel.h dataset.h distances.h parameter.h)
target_link_libraries(KronMat ${CBLAS_LIBRARIES} ${Boost_LIBRARIES})