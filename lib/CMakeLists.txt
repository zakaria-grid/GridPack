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
# Last Change: 2013-11-06 08:36:18 d3g096
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
