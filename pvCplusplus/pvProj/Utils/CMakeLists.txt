cmake_minimum_required(VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CAMKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBARRY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
include_directories("../../EIGENDIR")
include_directories("../../LIBIGLDIR")

project("Utils")

file(GLOB utilsources ../Utils/*.cpp ../Utils/*.h)


add_library(utillib STATIC  ${utilsources})
