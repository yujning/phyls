# lib/CMakeLists.txt
add_subdirectory(alice)
add_subdirectory(nbulsi-mockturtle)

if(ALSO_ENABLE_STP)
  add_subdirectory(stp)
else()
  message(FATAL_ERROR "STP is required but ALSO_ENABLE_STP is OFF")
endif()
