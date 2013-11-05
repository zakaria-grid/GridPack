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
# Last Change: 2013-11-05 09:23:51 d3g096
# -------------------------------------------------------------

add_library(gridpack_timer OBJECT
  coarse_timer.cpp
)

set(target_libraries
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

set(target_objects
  $<TARGET_OBJECTS:gridpack_timer>
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(AFTER ${GA_INCLUDE_DIRS})

# -------------------------------------------------------------
# TEST: test_timer
# -------------------------------------------------------------

add_executable(test_timer 
  test/test_timer.cpp
  ${target_objects}
)
target_link_libraries(test_timer ${target_libraries})
add_test(test_timer test_timer)
