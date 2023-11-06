cmake_minimum_required(VERSION 3.28)
project(Riemann)

include(CTest)
enable_testing()

include_directories(${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
set(CMAKE_CXX_STANDARD 17)
set( CMAKE_BUILD_TYPE "Debug" )
# FILEPATH: /c:/Users/YUJJ/Desktop/riemann/CMakeLists.txt

# Add all source files in this directory and its subdirectories to the SOURCES variable
file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# Add executable target

add_executable(${PROJECT_NAME} ${SOURCES})
