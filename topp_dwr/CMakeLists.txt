cmake_minimum_required(VERSION 3.5)
project(topp_dwr)
# Default to C++17
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic -g -Wno-unused -Wno-pedantic)
endif()

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(mosek64 REQUIRED)
find_package(ecos REQUIRED)
find_package(geometry_msgs REQUIRED)
find_package(visualization_msgs REQUIRED)
find_package(yaml_cpp_vendor REQUIRED)
find_package(ament_index_cpp REQUIRED)
find_package(Eigen3)

include_directories(
  include
   ${mosek64_INCLUDE_DIRS}
   ${ecos_INCLUDE_DIRS}
)

add_library(${PROJECT_NAME} SHARED
  src/topp_dwr.cpp
  src/tools.cpp
  src/non_uniform_bspline.cpp
)
ament_target_dependencies(${PROJECT_NAME}
  rclcpp
  visualization_msgs
  yaml_cpp_vendor
  ament_index_cpp
)
file(GLOB SO_LIBERARIES "${mosek64_INCLUDE_DIRS}/../lib/*.so")
target_link_libraries(${PROJECT_NAME} ${SO_LIBERARIES})
target_link_libraries(${PROJECT_NAME}
  -lglog
  ${cvte_log_LIBRARIES}
  ecos
)

add_executable(${PROJECT_NAME}_node
  src/main.cpp
)
ament_target_dependencies(${PROJECT_NAME}_node
  rclcpp
)
target_link_libraries(${PROJECT_NAME}_node
  ${PROJECT_NAME}
  -lglog
  ${cvte_log_LIBRARIES}
)

install(TARGETS
  ${PROJECT_NAME}
  LIBRARY DESTINATION lib
)
install(TARGETS
  ${PROJECT_NAME}_node
  DESTINATION lib/${PROJECT_NAME}
)

install(
  DIRECTORY include/
  DESTINATION include/
)
install(
  DIRECTORY files/  
  DESTINATION files/
)
install(
  DIRECTORY params/  
  DESTINATION params/
)
ament_export_include_directories(include ${ecos_INCLUDE_DIRS})
ament_export_dependencies(visualization_msgs)
ament_export_dependencies(yaml_cpp_vendor)
ament_export_dependencies(ecos)
ament_export_libraries(${PROJECT_NAME})
ament_export_libraries(mosek64)
ament_export_libraries(ecos)
ament_package()
