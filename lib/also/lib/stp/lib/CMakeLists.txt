if (NOT TARGET fmt)
  add_library(fmt INTERFACE) # v6.3.0
  target_include_directories(fmt SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/fmt)
  target_compile_definitions(fmt INTERFACE FMT_HEADER_ONLY)
endif()
if (NOT TARGET kitty)
  add_library(kitty INTERFACE)
  target_include_directories(kitty SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/kitty/include)
endif()

