# cmake_minimum_required(VERSION 3.22)
# project(Riemann)

# include(CTest)
# enable_testing()

# include_directories(${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
# # Eigen库https://blog.csdn.net/hongge_smile/article/details/107296658
# # find_package(Eigen3 REQUIRED)
# # include_directories(${EIGEN3_INCLUDE_DIR})



# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
# set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
# include(CPack)
# set(CMAKE_CPP_STANDARD 17)
# set( CMAKE_BUILD_TYPE "Debug" )

# # FILEPATH: /c:/Users/YUJJ/Desktop/riemann/CMakeLists.txt

# # Add all source files in this directory and its subdirectories to the SOURCES variable
 file(GLOB_RECURSE SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/*.h")

# # Add executable target

# add_executable(${PROJECT_NAME} ${SOURCES})


cmake_minimum_required(VERSION 3.22)
project(Riemann)
include(CTest)
enable_testing()
include(CPack)
set(CMAKE_CPP_STANDARD 17)
set( CMAKE_BUILD_TYPE "Debug" )
# set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
#并行
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")


# 子目录
add_subdirectory(Fvm)
# add_subdirectory(GaussSolver)

#外部库
include_directories(${EIGEN3_INCLUDE_DIR})
find_package(Eigen3 REQUIRED)
# 查找 MPI
find_package(MPI REQUIRED)

# 添加 MPI 头文件路径
include_directories(${MPI_INCLUDE_PATH})

# 查找 SuiteSparse 库
# include_directories(${LD_INCLUDE_DIR})
# find_package(SuiteSparse REQUIRED)


#查找目录下所有的源文件
aux_source_directory(. SOURCES)
#生成可执行文件
add_executable(${PROJECT_NAME} ${SOURCES})

#生成链接库
target_link_libraries(${PROJECT_NAME} Fvm)
# target_link_libraries(${PROJECT_NAME} GaussSolver)
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen)
# 链接 MPI 库
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})
# target_link_libraries(${PROJECT_NAME} PRIVATE SuiteSparse::umfpack)

# export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH
