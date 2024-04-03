cmake_minimum_required(VERSION 3.0.2)
project(eye_patch)

## Find catkin macros and libraries
find_package(catkin REQUIRED COMPONENTS
  eigen_conversions
  pcl_conversions
  roscpp
  rospy
  std_msgs
  tf2_ros
  tf_conversions
)

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

catkin_package(
   INCLUDE_DIRS include
   CATKIN_DEPENDS roscpp rospy std_msgs tf2_ros tf_conversions eigen_conversions pcl_conversions
   DEPENDS
)

# build
include_directories(
    include
    ${catkin_INCLUDE_DIRS}
)

add_executable(         ${PROJECT_NAME}_node src/eye_patch_ros.cpp)
target_link_libraries(  ${PROJECT_NAME}_node ${catkin_LIBRARIES} ${PYTHON_LIBRARIES})

# # install scripts
# catkin_install_python(PROGRAMS
#   python/read_vilens_output.py
#   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )