cmake_minimum_required(VERSION 3.10)

# set the project name
project(TBD VERSION 1.0)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# use mpicc and mpicxx to compile
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# compiler flags (turn on all the warnings)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wno-reorder")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# folders for executable and libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

# find packages
if(NOT DEFINED ENV{PETSC_DIR})
  message(FATAL_ERROR "ERROR: Could not find Petsc...")
endif()

#find_package(Eigen3 3.3 REQUIRED)

# -----------------------------
# for version control
# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force version.cmake to always be re-run before the build
#add_custom_target(VersionHeader
#  ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/version.cmake
#)

#ADD_CUSTOM_COMMAND(
#    OUTPUT ${CMAKE_SOURCE_DIR}/version.h
#           ${CMAKE_SOURCE_DIR}/_version.h
#    COMMAND cmake -P
#            ${CMAKE_SOURCE_DIR}/version.h)
# -----------------------------

# add include and link directories 
include_directories("${CMAKE_SOURCE_DIR}")
#include_directories(${CMAKE_SOURCE_DIR}/GeoTools)
#include_directories(${CMAKE_SOURCE_DIR}/MathTools)
link_directories(${CMAKE_SOURCE_DIR}/parser)

# compile the parser
add_custom_target(
   extern_lib
   COMMAND $(MAKE)
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/parser
)

# add the executable
add_executable(erick
main.cpp
defs.cpp
Utils.cpp
IoData.cpp
MeshGenerator.cpp
SpaceOperator.cpp
SpaceVariable.cpp)


# link to libraries
target_link_libraries(erick petsc mpi parser)
#target_link_libraries(erick ${CMAKE_DL_LIBS}) #linking to the dl library (-ldl)
add_dependencies(erick extern_lib)
#add_dependencies(erick VersionHeader)
