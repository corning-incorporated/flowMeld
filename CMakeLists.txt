#************************************************************************************#
#  Copyright 2025. Corning Incorporated. All rights reserved.                        #                                                                                     #
#  This software may only be used in accordance with the identified license(s).      #  
#                                                                                    #                                                                                      
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR        # 
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,          #
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL           #
#  CORNING BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN        #
#  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN                 #
#  CONNECTION WITH THE SOFTWARE OR THE USE OF THE SOFTWARE.                          #
#  Authors:                                                                          #
# Hamed Haddadi Staff Scientist                                                      #
#              haddadigh@corning.com                                                 #
# David Heine   Principal Scientist and Manager                                      #
#               heinedr@corning.com                                                  #
#*********************************************************************************** #



if(CI)
###############################################################################
# NOTE: run all compilation of flowMeld from a build directory to avoid expansions

set(EXECUTABLE_NAME "mpflow")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

###############################################################################

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
# add_executable(${PROJECT_NAME} "./lbmImplementations/*.cpp")

# Link with the following libraries
target_link_libraries(${PROJECT_NAME} palabos)
if(ENABLE_MPI)
    target_link_libraries(${PROJECT_NAME} ${MPI_CXX_LIBRARIES})
endif()

###############################################################################
return()

###############################################################################
endif()

cmake_minimum_required(VERSION 2.8.12)

project(mpflow)
enable_language(CXX)

set(EXECUTABLE_NAME "mpflow")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
    #set(CMAKE_CONFIGURATION_TYPES "Debug;Release")
endif()
message("Generated with config types: ${CMAKE_CONFIGURATION_TYPES}")
message(${CMAKE_BUILD_TYPE})

# Compiler flags
# Append flags: set(CMAKE_XXX_FLAGS "${CMAKE_XXX_FLAGS} ...")
if(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
    message("GCC.")
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL Clang)
    message("Clang.")
    set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wnon-virtual-dtor")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -DPLB_DEBUG -O0")
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
    message("MSVC.")
    set(CMAKE_CXX_FLAGS_RELEASE "/Ox /Ot /GS- /GL /DNDEBUG")
    set(CMAKE_CXX_FLAGS_DEBUG "/DPLB_DEBUG")
    set(CMAKE_EXE_LINKER_FLAGS_RELEASE "/LTCG /INCREMENTAL:NO /OPT:REF")
else()
    message( FATAL_ERROR "CXX compiler not recognized. CMake will quit." )
endif()

option(ENABLE_MPI "Enable MPI" ON)
if(ENABLE_MPI)
    message("Enabling MPI")
    find_package(MPI REQUIRED)
    if(MPI_CXX_FOUND)
        #set(CMAKE_CXX_COMPILER "${MPI_CXX_COMPILER}")
        include_directories(${MPI_CXX_INCLUDE_PATH})
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS}")
        add_definitions(-DPLB_MPI_PARALLEL)
    endif()
endif()

if(WIN32)
    option(ENABLE_POSIX "Enable POSIX" OFF)
else()
    option(ENABLE_POSIX "Enable POSIX" ON)
endif()

if(ENABLE_POSIX)
    message("Enabling POSIX")
    add_definitions(-DPLB_USE_POSIX)
endif()

if(APPLE)
    add_definitions(-DPLB_MAC_OS_X)
endif()

if(WIN32 OR CYGWIN)
    add_definitions(-DPLB_WINDOWS)
endif()

###############################################################################
# include Palabos Library to the build files of corningLBM
###############################################################################

include_directories("./plbLibrary/src")
include_directories("./plbLibrary/externalLibraries")
include_directories("./plbLibrary/externalLibraries/Eigen3")

file(GLOB_RECURSE PALABOS_SRC "./plbLibrary/src/*.cpp")
file(GLOB_RECURSE EXT_SRC "./plbLibrary/externalLibraries/tinyxml/*.cpp")

add_library(palabos STATIC ${PALABOS_SRC} ${EXT_SRC})

###############################################################################
# NOTE: define the path to TOMA implementation files here 
file(GLOB_RECURSE LBM_SOURCE "./lbmImplementations/*.cpp")
#file(GLOB_RECURSE HELPER_SOURCE "./helpersImplementations/*.cpp")
# iclude this in the command below ${HELPER_SOURCE}
add_executable(${EXECUTABLE_NAME} ${LBM_SOURCE})

# Link with the following libraries
target_link_libraries(${EXECUTABLE_NAME} palabos)
if(ENABLE_MPI)
    target_link_libraries(${EXECUTABLE_NAME} ${MPI_CXX_LIBRARIES})
endif()
