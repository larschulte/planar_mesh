cmake_minimum_required(VERSION 3.8)
project(mesh_object)

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

# ceres solver
find_package(Ceres REQUIRED)

# MeshObject
set(MeshObject_SOURCE 
  src/MeshObject/Vertex.cpp
  src/MeshObject/Edge.cpp
  src/MeshObject/Face.cpp
  src/MeshObject/Surface.cpp
  src/MeshObject/Storage.cpp
  src/MeshObject/GenericPoint.cpp
  src/MeshObject/InteriorPoint.cpp
  src/MeshObject/TriangleBVH.cpp
  src/MeshObject/EdgeBVH.cpp
)

# add executable
add_executable(main src/main.cpp ${MeshObject_SOURCE})
target_link_libraries (main ${PCL_LIBRARIES})
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
add_executable(Surface_test src/test/Surface_test.cpp ${MeshObject_SOURCE})
target_link_libraries (Surface_test ${PCL_LIBRARIES})
target_include_directories(Surface_test
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# add executable
add_executable(RRSTree_test src/test/RRSTree_test.cpp ${MeshObject_SOURCE})
target_link_libraries (RRSTree_test ${PCL_LIBRARIES})
target_include_directories(RRSTree_test
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)


# add executable
add_executable(MeshObject_test src/test/MeshObject_test.cpp ${MeshObject_SOURCE})
target_link_libraries (MeshObject_test ${PCL_LIBRARIES})
target_include_directories(MeshObject_test
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# add executable
add_executable(TriangleBVH_test src/test/TriangleBVH_test.cpp ${MeshObject_SOURCE})
target_link_libraries (TriangleBVH_test ${PCL_LIBRARIES})
target_include_directories(TriangleBVH_test
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# add executable
add_executable(optimization src/utilities/optimization.cpp)
target_link_libraries (optimization ${PCL_LIBRARIES} Ceres::ceres)
target_include_directories(optimization
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# call ament build system
ament_package()