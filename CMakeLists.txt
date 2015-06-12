cmake_minimum_required(VERSION 3.2)
project(KronMat)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMake/Modules)

include(/Users/tom/Code/dlib/dlib/cmake)

find_package(BOOST REQUIRED timer)
find_package(BLAS REQUIRED)
find_package(CBLAS REQUIRED)
find_path(BLAS_INCLUDE_DIRS cblas.h
         HINTS
        "/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Versions/Current/Headers/"
        "/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk/System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/Headers/")


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O0")

include_directories(${BLAS_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

set(SOURCE_FILES play.cpp)
add_executable(KronMat ${SOURCE_FILES} kronecker_matrix.cpp kronecker_matrix.h static_functions.h matrix.h blas_wrap.h)
target_link_libraries(KronMat ${CBLAS_LIBRARIES} dlib)