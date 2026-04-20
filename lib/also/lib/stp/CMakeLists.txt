cmake_minimum_required(VERSION 3.19)

# Set C++ standard
set(STP_CXX_STANDARD "17" CACHE STRING "C++ standard")
set(CMAKE_CXX_STANDARD ${STP_CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(stp LANGUAGES CXX)

# Set Eigen include directory
include_directories(/usr/local/include/eigen3)

# Options
option(TEST "Build tests" ON)

add_compile_options(-g -O3)


# Build tests if enabled and available
if(TEST)
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/test/CMakeLists.txt")
        add_subdirectory(test)
    else()
        message(WARNING "Tests requested for stp but no test directory with CMakeLists.txt found; skipping.")
    endif()
endif()

add_subdirectory(lib)
add_subdirectory(src)