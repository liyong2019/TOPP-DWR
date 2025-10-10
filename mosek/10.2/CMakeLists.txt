cmake_minimum_required(VERSION 3.5)
project(mosek64)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)  
endif()
message( "CMAKE_SYSTEM_PROCESSOR: ${CMAKE_SYSTEM_PROCESSOR}")
SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

 if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
   add_compile_options(-Wall -Wextra -Wpedantic -g -Wno-unused -Wno-pedantic)
 endif()

set(MOSEK_LICENSE_FILE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/mosek.lic")

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)

if(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
include_directories(
  tools/platform/linux64x86/h
)
install(
    DIRECTORY tools/platform/linux64x86/h/
    DESTINATION include
  )
install(
    DIRECTORY tools/platform/linux64x86/bin/
    DESTINATION lib/
)
elseif(${CMAKE_SYSTEM_PROCESSOR} STREQUAL "aarch64")
include_directories(
  tools/platform/linuxaarch64/h
)
install(
    DIRECTORY tools/platform/linuxaarch64/h/
    DESTINATION include
  )
install(
    DIRECTORY tools/platform/linuxaarch64/bin/
    DESTINATION lib/
)
else()
message(WARNING "Unrecognized processor architecture.")
endif()
ament_export_include_directories(include)
ament_package()
