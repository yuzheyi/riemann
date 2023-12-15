# # 查找当前目录下的所有源文件
# # 并将名称保存到 DIR_LIB_SRCS 变量
# aux_source_directory(. DIR_LIB_Gauss)




# # 生成链接库
# add_library (GaussSolver ${DIR_LIB_Gauss})

# # target_link_libraries(${PROJECT_NAME} Fvm)
# target_include_directories(GaussSolver PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Fvm)