cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CAMKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBARRY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

project("MD")

#message(STATUS, "${CMAKE_CURRENT_BINARY_DIR}")
file(GLOB sources ./*.cpp
  ./*.h)


file(GLOB data "../Resources/*")
include_directories("../../eigen")
include_directories("../../libigl")

#"../Utils"

add_library(mdlib STATIC ${data} ${sources})
target_link_libraries(mdlib utillib)

add_subdirectory(Test)
add_subdirectory(Benchmark)
