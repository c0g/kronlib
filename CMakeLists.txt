cmake_minimum_required(VERSION 3.2)
project(KronMat)
add_definitions(-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA)

# Additional 'FindX.cmake'
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake/Modules)

include("dlib/cmake")
# Google test framework
include("test/gtest.cmake")

find_package(OpenBLAS)
find_package(CUDA)
find_package(TBB)
message("CUDA version ${CUDA_VERSION}")
if(${CUDA_VERSION} LESS 7) 
    message("CUDA version < 7 - using homebrew Cholesky")
    add_definitions("-DKRONMAT_CHOL_VERSION=CUSTOM")
else()
    message("CUDA version >= 7 - using CUSolver Cholesky")
    add_definitions("-DKRONMAT_CHOL_VERSION=CUSOLVER")
endif()

# Compiler flags - if GNU add in openmp
# set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -G -g")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3") 

# Includes
include_directories(${CUDA_INCLUDE_DIRS} ${OpenBLAS_INCLUDE_DIR} ${GTEST_INCLUDE_DIR} thrust src inc)

## Source files for executable
set(SOURCE_FILES src/play.cu)
if(CUDA_FOUND)
	CUDA_ADD_EXECUTABLE(KronMat ${SOURCE_FILES})
	CUDA_ADD_CUBLAS_TO_TARGET(KronMat)
else()
	add_executable(KronMat ${SOURCE_FILES})
endif()

target_link_libraries(KronMat ${OpenBLAS_LIB} ${TBB_LIBRARY} ${GTEST_BOTH_LIBRARIES}  ${CUDA_LIBRARIES} dlib)

# Tests
enable_testing()

if(CUDA_FOUND)
	CUDA_ADD_EXECUTABLE(matrix_test test/testmatrix.cu)
	CUDA_ADD_CUBLAS_TO_TARGET(matrix_test)
    target_link_libraries(matrix_test ${OpenBLAS_LIB}  ${TBB_LIBRARY} ${GTEST_LIBRARY} ${CUDA_LIBRARIES} dlib)
	add_test(matrix_test matrix_test)
endif()

if(CUDA_FOUND)
	CUDA_ADD_EXECUTABLE(chol_test test/testchol.cu)
	CUDA_ADD_CUBLAS_TO_TARGET(chol_test)
    target_link_libraries(chol_test ${OpenBLAS_LIB}  ${TBB_LIBRARY} ${GTEST_LIBRARY} ${CUDA_LIBRARIES} dlib)
endif()

if(CUDA_FOUND)
	CUDA_ADD_EXECUTABLE(kron_test test/testkron.cu)
	CUDA_ADD_CUBLAS_TO_TARGET(kron_test)
    target_link_libraries(kron_test ${OpenBLAS_LIB}  ${TBB_LIBRARY} ${GTEST_LIBRARY} ${CUDA_LIBRARIES} dlib)
    add_test(kron_test kron_test)
endif()

if(FALSE) # fake multiline comment
if(CUDA_FOUND)
	CUDA_ADD_EXECUTABLE(kvs_test test/testkvs.cpp)
	CUDA_ADD_CUBLAS_TO_TARGET(kvs_test)
else()
	add_executable(kvs_test test/testkvs.cpp)
endif()
target_link_libraries(kvs_test ${CBLAS_LIBRARIES}  ${GTEST_LIBRARY}  ${LAPACKE_LIBRARIES})
add_test(kvs_test kvs_test)

if(CUDA_FOUND)
	CUDA_ADD_EXECUTABLE(sqexp_test test/testsqexp.cpp)
	CUDA_ADD_CUBLAS_TO_TARGET(sqexp_test)
else()
	add_executable(sqexp_test test/testsqexp.cpp)
endif()
target_link_libraries(sqexp_test ${CBLAS_LIBRARIES}  ${GTEST_LIBRARY}  ${LAPACKE_LIBRARIES})
add_test(sqexp_test sqexp_test)
endif()
