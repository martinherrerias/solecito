cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CAMKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBARRY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

project("EC")

#message(STATUS, "${CMAKE_CURRENT_BINARY_DIR}")
file(GLOB sources ./*.cpp
  ./*.h)


file(GLOB data "../Resources/*")
include_directories("../../EIGENDIR")
include_directories("../../LIBIGLDIR")

#"../Utils"

add_library(eclib STATIC ${data} ${sources})
target_link_libraries(eclib utillib)

add_subdirectory(Test)
add_subdirectory(Benchmark)
