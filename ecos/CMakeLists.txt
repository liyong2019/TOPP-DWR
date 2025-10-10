cmake_minimum_required(VERSION 3.5)
project(ecos)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)  
endif()

 if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
   add_compile_options(-Wall -Wextra -Wpedantic -g -Wno-unused -Wno-pedantic )
 endif()

 SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall ")

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)

set(INCLUDE_DIRS
    include/
    ${GLOG_INCLUDE_DIRS}/
    external/amd/include/
    external/ldl/include/
    external/SuiteSparse_config/
)
include_directories( ${INCLUDE_DIRS})
add_library(ecos SHARED
   # AMD
   external/amd/src/amd_1.c
   external/amd/src/amd_2.c
   external/amd/src/amd_aat.c
   external/amd/src/amd_control.c
   external/amd/src/amd_defaults.c
   external/amd/src/amd_dump.c
   external/amd/src/amd_global.c
   external/amd/src/amd_info.c
   external/amd/src/amd_order.c
   external/amd/src/amd_post_tree.c
   external/amd/src/amd_postorder.c
   external/amd/src/amd_preprocess.c
   external/amd/src/amd_valid.c
   # LDL
   external/ldl/src/ldl.c
   # ECOS
   src/cone.c
   src/ctrlc.c
   src/ecos.c
   src/equil.c
   src/expcone.c
   src/kkt.c
   src/preproc.c
   src/spla.c
   src/splamm.c
   src/timer.c
   src/wright_omega.c
   # ECOS BB
   ecos_bb/ecos_bb.c
   ecos_bb/ecos_bb_preproc.c
)
target_link_libraries(ecos -lglog )
install(TARGETS
  ecos
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
)

foreach(DIR ${INCLUDE_DIRS})
    install(
        DIRECTORY ${DIR}
        DESTINATION ${DIR}
    )
endforeach()
ament_export_include_directories(include)
ament_export_include_directories(external)
ament_export_include_directories(external/SuiteSparse_config)
ament_export_libraries(ecos)
ament_package()
