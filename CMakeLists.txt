# CMake最低版本要求
cmake_minimum_required(VERSION 3.12)

# 项目名称
project(MyFirstVTKApp)

# 关键步骤：告诉CMake去哪里找VTK的配置信息
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
# 这个路径需要指向你VTK安装目录下的 /lib/cmake/vtk-X.Y 目录
set(VTK_DIR "D:/vtk/VTK-8.2.0/lib/cmake/vtk-8.2")

# 查找VTK包，并设置为必需
find_package(VTK REQUIRED)

# 如果找到了VTK，则包含它的CMake配置
if (VTK_FOUND)
  include(${VTK_USE_FILE})
endif()

# 添加你的可执行文件
add_executable(my_app main.cpp)

# 将你的可执行文件与VTK库链接起来
target_link_libraries(my_app PRIVATE ${VTK_LIBRARIES})