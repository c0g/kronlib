cmake_minimum_required(VERSION 3.2)
project(KronMat)

add_definitions(-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA)

# Additional 'FindX.cmake'
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake/Modules)
message("${CMAKE_MODULE_PATH}")

# Google test framework
include("test/gtest.cmake")

find_package(OpenBLAS)
find_package(CUDA)
find_package(TBB)

# Compiler flags - if GNU add in openmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

# Includes
include_directories(${CUDA_INCLUDE_DIRS} ${OpenBLAS_INCLUDE_DIR} ${GTEST_INCLUDE_DIR} src inc)

## Source files for executable
set(SOURCE_FILES src/play.cpp)
if(CUDA_FOUND)
	CUDA_ADD_EXECUTABLE(KronMat ${SOURCE_FILES})
	CUDA_ADD_CUBLAS_TO_TARGET(KronMat)
else()
	add_executable(KronMat ${SOURCE_FILES})
endif()

target_link_libraries(KronMat ${OpenBLAS_LIB} ${TBB_LIBRARY} ${GTEST_BOTH_LIBRARIES}  ${CUDA_LIBRARIES})

# Tests
enable_testing()

if(CUDA_FOUND)
	CUDA_ADD_EXECUTABLE(matrix_test_cuda test/testmatrix.cu)
	CUDA_ADD_CUBLAS_TO_TARGET(matrix_test_cuda)
    target_link_libraries(matrix_test_cuda ${OpenBLAS_LIB}  ${TBB_LIBRARY} ${GTEST_LIBRARY} ${CUDA_LIBRARIES})
	add_test(matrix_test_cuda matrix_test_cuda)
endif()
#add_executable(matrix_test test/testmatrix.cpp)
#target_link_libraries(matrix_test ${CBLAS_LIBRARIES}  ${GTEST_LIBRARY}  ${LAPACKE_LIBRARIES})
#add_test(matrix_test matrix_test)

if(FALSE) # fake multiline comment
if(CUDA_FOUND)
	CUDA_ADD_EXECUTABLE(kron_test test/testkron.cpp)
	CUDA_ADD_CUBLAS_TO_TARGET(kron_test)
else()
	add_executable(kron_test test/testkron.cpp)
endif()
target_link_libraries(kron_test ${CBLAS_LIBRARIES}  ${GTEST_LIBRARY}  ${LAPACKE_LIBRARIES})
add_test(kron_test kron_test)

if(CUDA_FOUND)
	CUDA_ADD_EXECUTABLE(chol_test test/testchol.cpp)
	CUDA_ADD_CUBLAS_TO_TARGET(chol_test)
else()
	add_executable(chol_test test/testchol.cpp)
endif()
target_link_libraries(chol_test ${CBLAS_LIBRARIES}  ${GTEST_LIBRARY} ${LAPACKE_LIBRARIES} )
add_test(chol_test chol_test)

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
