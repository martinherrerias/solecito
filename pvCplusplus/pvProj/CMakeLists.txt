cmake_minimum_required(VERSION 2.8.12)
SET(BOOST_MIN_VERSION "1.70.0")

# Prevent FindBoost.cmake from looking for system Boost.cmake files
set(Boost_NO_BOOST_CMAKE true)
find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)

project(PV CXX)
set(PROJECT_DESCRIPTION "A solar simulation.")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#set(CMAKE_CXX_COMPILER "/opt/cray/pe/craype/2.5.18/bin/CC")
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-DNDEBUG -DNMEX -fopenmp -funroll-loops -ffast-math -O3")

if(NOT CMAKE_CXX_PLATFORM)
	set(CMAKE_CXX_PLATFORM "Local")
endif()

if(NOT CMAKE_CXX_COMPILER_ID)
	set(CMAKE_CXX_COMPILER_ID "GNU")
endif()

set(BOOST_INCLUDE_DIRS ${CMAKE_CXX_BOOST_LIB})
include_directories(${BOOST_INCLUDE_DIRS})

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
	if(CMAKE_CXX_PLATFORM MATCHES AMD)
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=znver2")
	else()
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
	endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "7.5.0")
    	message(FATAL_ERROR "gcc version >= 7.5.0 required.\n")
	endif()
endif()

#set(CMAKE_BUILD_FILES_DIRECTORY ${CAMKE_CURRENT_SOURCE_DIR}/build)
#set(CMAKE_BUILD_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build)
#set(CMAKE_BINARY_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/build)

add_subdirectory(Utils)
add_subdirectory(ElectricalCalculation)
add_subdirectory(ModuleDefinition)


