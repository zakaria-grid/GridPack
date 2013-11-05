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
# Last Change: 2013-11-05 09:16:45 d3g096
# -------------------------------------------------------------

set(target_libraries
  ${GA_LIBRARIES} 
  ${PETSC_LIBRARIES}
  ${Boost_LIBRARIES} 
  ${MPI_CXX_LIBRARIES}
)

set(target_objects 
  $<TARGET_OBJECTS:gridpack_parallel>
  $<TARGET_OBJECTS:gridpack_configuration>
  $<TARGET_OBJECTS:gridpack_math>
  $<TARGET_OBJECTS:gridpack_components>
)

include_directories(AFTER ${GA_INCLUDE_DIRS})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()

# -------------------------------------------------------------
# TEST: test_serial_io
# A simple Serial IO test program
# -------------------------------------------------------------
add_executable(test_serial_io test/test_serial_io.cpp ${target_objects})
target_link_libraries(test_serial_io ${target_libraries})

