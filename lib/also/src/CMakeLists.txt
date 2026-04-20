include_directories(commands)

file(GLOB_RECURSE FILENAMES *.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(also also.cpp ${FILENAMES})
target_link_libraries(also alice mockturtle)

# STP bi-decomposition resynthesis headers live under lib/stp/src/include
# Add the include directory so the also executable can find them when the
# lut_resyn command is built with STP support.
target_include_directories(also PRIVATE 
  ${CMAKE_SOURCE_DIR}/lib/stp/src/include
)

# Link against STP algorithms library to get other dependencies
if(ALSO_HAS_STP)
  message(STATUS "STP support is enabled, adding Eigen3 includes")
  # Eigen3 is required for STP. Try to find it; if not found, use standard location.
  # Users can override by setting CMAKE_PREFIX_PATH or Eigen3_DIR.
  find_package(Eigen3 QUIET NO_MODULE)
  if(TARGET Eigen3::Eigen)
    # Try to get include dirs from target
    get_target_property(_eigen_includes Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    if(_eigen_includes)
      target_include_directories(also PRIVATE ${_eigen_includes})
      message(STATUS "Using Eigen3 from target: ${_eigen_includes}")
    endif()
  endif()
  
  # Also add standard location as fallback
  if(EXISTS "/usr/include/eigen3")
    target_include_directories(also PRIVATE /usr/include/eigen3)
    message(STATUS "Adding Eigen3 fallback location: /usr/include/eigen3")
  endif()
  
  target_link_libraries(also alice mockturtle stp_algorithms)
else()
  message(STATUS "STP support is NOT enabled")
  target_link_libraries(also alice mockturtle)
endif()