cmake_minimum_required(VERSION 3.19)

add_library(stp_algorithms INTERFACE)
find_package(Eigen3 REQUIRED)
target_include_directories(stp_algorithms
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/lib
)
target_link_libraries(stp_algorithms
    INTERFACE
        kitty
        mockturtle
        Eigen3::Eigen
)
target_compile_features(stp_algorithms INTERFACE cxx_std_17)

if(STP_BUILD_TOOL)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

    set(SOURCES stp.cpp)
    add_executable(stp ${SOURCES} ${FILENAMES})
    find_package(OpenMP REQUIRED)
    if (OpenMP_CXX_FOUND)
        target_link_libraries(stp PRIVATE OpenMP::OpenMP_CXX)
endif()

target_include_directories(stp PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/commands)
target_link_libraries(stp PRIVATE
alice
fmt
kitty
mockturtle
stp_algorithms
pthread
Eigen3::Eigen
)

if(STP_ENABLE_CUDA AND CMAKE_CUDA_COMPILER)
add_subdirectory(cuda)
target_compile_definitions(stp PUBLIC ENABLE_CUDA)
endif()
endif()
