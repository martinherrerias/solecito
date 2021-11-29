cmake_minimum_required(VERSION 2.8)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#message(STATUS, "${CMAKE_CURRENT_SOURCE_DIR}")
set(EXECUTABLE_OUTPUT_PATH "./exec_test")

file(GLOB sources_test "./*.cpp")
file(GLOB test_header "../../test.h")

#add_executable(testClip ./testClip.cpp ../test.h)
#target_link_libraries(testClip mdlib)
#add_test(NAME testClip COMMAND testClip)
foreach(testfile ${sources_test})
  get_filename_component(_dir ${testfile} NAME_WE)
  add_executable(${_dir} ${testfile} ${test_header})
  target_link_libraries(${_dir} mdlib)
endforeach(testfile ${sources_test})

