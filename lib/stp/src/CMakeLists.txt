
# add_library(stp_lib INTERFACE)

# # Include stp headers
# target_include_directories(stp_lib SYSTEM INTERFACE ${PROJECT_SOURCE_DIR}/include)
# target_include_directories(stp_lib SYSTEM INTERFACE ${PROJECT_SOURCE_DIR}/commands)

# #set output directory
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# set(SOURCES stp.cpp)
# add_executable(stp ${SOURCES} ${FILENAMES})
# target_link_libraries(stp alice)


# if(STP_ENABLE_CUDA AND CMAKE_CUDA_COMPILER)
#     add_subdirectory(cuda)
#     target_compile_definitions(stp PUBLIC ENABLE_CUDA)
# endif()


# # Try to find Eigen3 on the system
# find_package(Eigen3 QUIET)
# message(STATUS "Eigen version: ${EIGEN3_VERSION}")

# # If Eigen3 is not found, use ExternalProject to clone and build Eigen3
# if(NOT Eigen3_FOUND)
#   message("-- Eigen3 is not found. Please clone and build from GitLab...")
#   message("git clone https://gitlab.com/libeigen/eigen.git")
#   message("cd eigen")
#   message("mkdir build")
#   message("cd build")
#   message("cmake ..")
#   message("make install")
#   message(FATAL_ERROR "Please visit https://eigen.tuxfamily.org/ for more details.")
# else()
# message("found")
# message(STATUS "Eigen3 found at: ${EIGEN3_INCLUDE_DIRS}")
#  # Add Eigen3 include directory to stp target
#   target_include_directories(stp INTERFACE ${EIGEN3_INCLUDE_DIRS})
# endif()

# find_package(OpenMP REQUIRED)
# target_link_libraries(stp OpenMP::OpenMP_CXX)
# target_link_libraries(stp  fmt Eigen3::Eigen pthread stp_lib)

# ===============================================
#  STP main CMakeLists.txt (src/)
# ===============================================

add_library(stp_lib INTERFACE)

# Include STP headers
target_include_directories(stp_lib SYSTEM INTERFACE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/commands
)

# Set output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Main executable
set(SOURCES stp.cpp)
add_executable(stp ${SOURCES} ${FILENAMES})


# -----------------------------------------------
# OpenMP  🔥 重点就这几行
# -----------------------------------------------
find_package(OpenMP REQUIRED)
if (OpenMP_CXX_FOUND)
    target_link_libraries(stp PRIVATE OpenMP::OpenMP_CXX)
endif()

# -----------------------------------------------
# Link core libraries
# -----------------------------------------------
target_link_libraries(stp PRIVATE
    alice
    fmt
    kitty
    mockturtle       
    stp_lib
)

# -----------------------------------------------
# CUDA (optional)
# -----------------------------------------------
if(STP_ENABLE_CUDA AND CMAKE_CUDA_COMPILER)
    add_subdirectory(cuda)
    target_compile_definitions(stp PUBLIC ENABLE_CUDA)
endif()

# -----------------------------------------------
# Eigen
# -----------------------------------------------
find_package(Eigen3 QUIET)
message(STATUS "Eigen version: ${EIGEN3_VERSION}")

if(NOT Eigen3_FOUND)
    message("-- Eigen3 not found. Please install manually:")
    message("git clone https://gitlab.com/libeigen/eigen.git")
    message("cd eigen && mkdir build && cd build")
    message("cmake .. && make install")
    message(FATAL_ERROR "Eigen3 is required but not found.")
else()
    message(STATUS "Eigen3 found at: ${EIGEN3_INCLUDE_DIRS}")
    target_include_directories(stp INTERFACE ${EIGEN3_INCLUDE_DIRS})
endif()

# -----------------------------------------------
# OpenMP
# -----------------------------------------------
find_package(OpenMP REQUIRED)
target_link_libraries(stp PRIVATE
    OpenMP::OpenMP_CXX
    pthread
    Eigen3::Eigen
)
