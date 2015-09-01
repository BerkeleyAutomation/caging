# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

project(caging)

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)
set (CAGING_LIB_TYPE "SHARED" CACHE STRING "Library type defaults to shared, options are: SHARED STATIC")
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-g -O3 -fPIC") 

if ( COMMAND cmake_policy )
  cmake_policy( SET CMP0003 NEW )  
endif()
 
# find libs
find_package(CGAL REQUIRED COMPONENTS Core )
find_package(Boost REQUIRED)
find_package(OpenCV REQUIRED)
find_package(OpenGL REQUIRED)
find_package(Eigen REQUIRED)

# hardcode important libraries / include dirs
set(PRM_INCLUDE_DIR include/)
set(3RDPARTY_INCLUDE_DIR 3rdparty/)
set(PHAT_INCLUDE_DIR 3rdparty/phat/include)

set(3RDPARTY_LIBRARIES freeglut_static glui)
set(BAYESOPT_LIBRARIES bayesopt)
set(BOX2D_LIBRARIES Box2D)
set(CCD_LIBRARIES ccd)
set(CGAL_LIBRARIES CGAL CGAL_Core CGAL_ImageIO CGAL_Qt3 CGAL_Qt4)
set(CGAL_3RD_PARTY_LIBRARIES ${CGAL_3RD_PARTY_LIBRARIES} gmp mpfr)
set(FCL_LIBRARIES fcl)
set(FCL_3RD_PARTY_LIBRARIES ccd)
set(GLOG_LIBRARIES glog)
set(OMPL_LIBRARIES ompl)
set(SOLID_LIBRARIES solid)
set(YAML_LIBRARIES yaml-cpp)

set(POLYDEPTH_LIBRARIES polydepth_Core polydepth_C2A polydepth_PQP)

set(DEPENDENCY_LIBS
  ${3RDPARTY_LIBRARIES}
  ${BAYESOPT_LIBRARIES}
  ${BOX2D_LIBRARIES}
  ${CCD_LIBRARIES}
  ${CGAL_LIBRARIES}
  ${CGAL_3RD_PARTY_LIBRARIES}   
  ${FCL_LIBRARIES}
  ${FCL_3RD_PARTY_LIBRARIES}
  ${GLOG_LIBRARIES}
  ${OMPL_LIBRARIES}
  ${OpenCV_LIBS}
  ${OPENGL_LIBRARIES}
  ${POLYDEPTH_LIBRARIES}
  ${SOLID_LIBRARIES}
  ${YAML_LIBRARIES}
)

# include stuff
include_directories (
  ${OPENGL_INCLUDE_DIR}
  ${CGAL_INCLUDE_DIR}
  ${Boost_INCLUDE_DIR}
  ${3RDPARTY_INCLUDE_DIR}
  ${FREEGLUT_INCLUDE_DIR}
  ${PHAT_INCLUDE_DIR}
  ${PRM_INCLUDE_DIR}
  ${Eigen_INCLUDE_DIR}
)

add_subdirectory(3rdparty)
add_subdirectory(src)

