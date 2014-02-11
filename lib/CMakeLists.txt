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
# Created October 17, 2013 by Bruce Palmer
# Last Change: 2014-02-11 07:37:16 d3g096
# -------------------------------------------------------------

add_library(gridpack 
  gridpack.cpp
  $<TARGET_OBJECTS:gridpack_parallel>
  $<TARGET_OBJECTS:gridpack_configuration>
  $<TARGET_OBJECTS:gridpack_math>
  $<TARGET_OBJECTS:gridpack_factory>
  $<TARGET_OBJECTS:gridpack_partition>
  $<TARGET_OBJECTS:gridpack_components>
  $<TARGET_OBJECTS:gridpack_timer>
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/GridPACK.cmake.in 
  ${CMAKE_CURRENT_BINARY_DIR}/GridPACK.cmake @ONLY)


set(GRIDPACK_TP_INCLUDE_FLAGS)
  
set(GRIDPACK_TP_LIBS)


# -------------------------------------------------------------
# installation
# -------------------------------------------------------------

install(TARGETS 
  gridpack
  DESTINATION lib
)

install(FILES 
  ${CMAKE_CURRENT_BINARY_DIR}/GridPACK.cmake
  DESTINATION lib
)
