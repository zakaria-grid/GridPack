#
#     Copyright (c) 2013 Battelle Memorial Institute
#     Licensed under modified BSD License. A copy of this license can be
#     found
#     in the LICENSE file in the top level directory of this distribution.
#
# -*- mode: cmake -*-
# -------------------------------------------------------------
# file: CMakeLists.txt
# -------------------------------------------------------------
# -------------------------------------------------------------
# Created May  6, 2013 by William A. Perkins
# Last Change: 2017-12-07 12:58:51 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_components
    gridpack_stream
    gridpack_partition
    gridpack_math
    gridpack_configuration
    gridpack_timer
    gridpack_parallel
    ${PARMETIS_LIBRARY} ${METIS_LIBRARY} 
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
    ${PETSC_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

add_executable(hello_world.x
   hw_app.cpp
   hw_components.cpp
   hw_main.cpp
)

target_link_libraries(hello_world.x ${target_libraries})

add_custom_target(hello_world.x.input

   COMMAND ${CMAKE_COMMAND} -E copy
   ${CMAKE_CURRENT_SOURCE_DIR}/10x10.raw
   ${CMAKE_CURRENT_BINARY_DIR} 

   DEPENDS
   ${CMAKE_CURRENT_SOURCE_DIR}/10x10.raw

)
add_dependencies(hello_world.x hello_world.x.input)

# -------------------------------------------------------------
# run application as test
# -------------------------------------------------------------
gridpack_add_run_test("hello_world" hello_world.x 10x10.raw)

# -------------------------------------------------------------
# install as a sample application
# -------------------------------------------------------------

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.install.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt @ONLY)

install(FILES 
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeLists.txt
  ${CMAKE_CURRENT_SOURCE_DIR}/hw_app.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hw_app.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hw_components.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hw_components.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hw_factory.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hw_main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/10x10.raw
  DESTINATION share/gridpack/example/hello
)

