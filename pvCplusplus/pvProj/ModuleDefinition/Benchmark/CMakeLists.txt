cmake_minimum_required(VERSION 2.8)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#message(STATUS, "${CMAKE_CURRENT_SOURCE_DIR}")
set(EXECUTABLE_OUTPUT_PATH "./exec_bench")

file(GLOB sources_benchmark "./*.cpp")
file(GLOB benchmark_header "../../testUtil.h")

#add_executable(testClip ./testClip.cpp ../test.h)
#target_link_libraries(testClip mdlib)
#add_test(NAME testClip COMMAND testClip)
foreach(benchmarkfile ${sources_benchmark})
  get_filename_component(_dir ${benchmarkfile} NAME_WE)
  add_executable(${_dir} ${benchmarkfile} ${benchmark_header})
  target_link_libraries(${_dir} mdlib)
endforeach(benchmarkfile ${sources_benchmark})

