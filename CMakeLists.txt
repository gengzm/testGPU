cmake_minimum_required(VERSION 3.13)
cmake_policy(VERSION 3.13)

set(DEFAULT_CXX_STANDARD ${CMAKE_CXX_STANDARD})
if (NOT DEFAULT_CXX_STANDARD)
    set(DEFAULT_CXX_STANDARD 14)
endif()
set(CMAKE_CXX_STANDARD ${DEFAULT_CXX_STANDARD} CACHE STRING "C++ standard(min 14)" FORCE)

project(CudaTestMe LANGUAGES CUDA CXX C)
# …Ë÷√CUAD±‡“Î≈‰÷√
set(CUDA_NVCC_FLAGS "-g -G")

if (WIN32)
    set(OpenCV_DIR "D:/opencv/build/" CACHE FILEPATH "Set to find OpenCVConfig.cmake")
else()
endif()
find_package(OpenCV REQUIRED NO_DEFAULT_PATH)
INCLUDE("C:/Program Files/CMake/share/cmake-3.14/Modules/FindCUDA.cmake")

include(CheckLanguage)
enable_language(CUDA)

#set(CUDA_PATH "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v10.2")
#set(CUDA_INCLUDE ${CUDA_PATH}/include)
#set(CUDA_LIB_DIR ${CUDA_PATH}/lib/x64/)
#file(GLOB CUDA_LIBS "${CUDA_LIB_DIR}/*.lib" ABSOLUTE)
message("------------------------ CUDA_LIBRARIES ${CUDA_LIBRARIES}---------------------")
message("------------------------ CUDA_INCLUDE_DIRS ${CUDA_INCLUDE_DIRS}---------------------")

set(CUDA_LIBRARIES PUBLIC ${CUDA_LIBRARIES})

cuda_add_executable(${PROJECT_NAME} Src/cuda_test.cu)
target_include_directories(${PROJECT_NAME} PRIVATE ${OpenCV_INCLUDE_DIRS} ${CUDA_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${OpenCV_LIBS} ${CUDA_LIBRARIES})

cuda_add_executable(${PROJECT_NAME}_2 Src/cuda_main.cpp Src/cuda_fun.cuh Src/cuda_fun.cu)
target_include_directories(${PROJECT_NAME}_2 PRIVATE ${OpenCV_INCLUDE_DIRS} ${CUDA_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}_2 PRIVATE ${OpenCV_LIBS} ${CUDA_LIBRARIES})


cuda_add_executable(${PROJECT_NAME}_cpu GPU/Test_Example.cpp)
target_include_directories(${PROJECT_NAME}_cpu PRIVATE )
target_link_libraries(${PROJECT_NAME}_cpu PRIVATE )