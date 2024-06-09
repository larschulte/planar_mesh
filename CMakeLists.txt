cmake_minimum_required(VERSION 3.8)
project(eye_patch)

list(APPEND CMAKE_PREFIX_PATH "/usr/local/")

# settings for complier version
find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif()

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)

# pcl library
find_package(PCL 1.14 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

# delaunay triangulations
#include <delaunator.hpp>
find_package(delaunator REQUIRED)
include_directories(${delaunator_INCLUDE_DIRS})

# matplotlib-cpp
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

# ceres solver
find_package(Ceres REQUIRED)

# add executable
add_executable(main src/main.cpp)
target_link_libraries (main ${PCL_LIBRARIES} ${PYTHON_LIBRARIES} delaunator)
target_include_directories(main
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# install(TARGETS main DESTINATION lib/${PROJECT_NAME})
install(TARGETS 
  main
  DESTINATION lib/${PROJECT_NAME}
)

# add executable
add_executable(temp src/temp.cpp src/TriangleBVH.cpp src/flann3d.cpp)
target_link_libraries (temp ${PCL_LIBRARIES} ${PYTHON_LIBRARIES} delaunator Ceres::ceres)
target_include_directories(temp
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# install(TARGETS main DESTINATION lib/${PROJECT_NAME})
install(TARGETS 
  temp
  DESTINATION lib/${PROJECT_NAME}
)

# add executable
add_executable(RRSTree_test src/RRSTree_test.cpp)
target_link_libraries (RRSTree_test ${PCL_LIBRARIES})
target_include_directories(RRSTree_test
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)


# add executable
add_executable(MeshObject_test src/MeshObject_test.cpp src/MeshObject/Vertex.cpp src/MeshObject/Edge.cpp src/MeshObject/Face.cpp src/MeshObject/Surface.cpp src/MeshObject/Storage.cpp src/MeshObject/GenericPoint.cpp src/MeshObject/InteriorPoint.cpp src/TriangleBVH.cpp)
target_link_libraries (MeshObject_test ${PCL_LIBRARIES})
target_include_directories(MeshObject_test
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)


# add executable
add_executable(optimization src/optimization.cpp)
target_link_libraries (optimization ${PCL_LIBRARIES} Ceres::ceres)
target_include_directories(optimization
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)









# add add_executable
add_executable(add_mesh_to_mesh src/add_mesh_to_mesh.cpp)
target_link_libraries (add_mesh_to_mesh ${PCL_LIBRARIES} ${PYTHON_LIBRARIES} delaunator)
target_include_directories(add_mesh_to_mesh
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# install(TARGETS main DESTINATION lib/${PROJECT_NAME})
install(TARGETS 
  add_mesh_to_mesh
  DESTINATION lib/${PROJECT_NAME}
)


# call ament build system
ament_package()