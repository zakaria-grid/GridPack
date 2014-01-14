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
# Last Change: 2014-01-14 07:49:12 d3g096
# -------------------------------------------------------------

set(target_libraries
  ${GA_LIBRARIES} 
  ${Boost_LIBRARIES} 
  ${MPI_CXX_LIBRARIES}
  ${PETSC_LIBRARIES}
)

set(target_objects 
  $<TARGET_OBJECTS:gridpack_parallel>
  $<TARGET_OBJECTS:gridpack_configuration>
  $<TARGET_OBJECTS:gridpack_math>
  $<TARGET_OBJECTS:gridpack_components>
  $<TARGET_OBJECTS:gridpack_timer>
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()
# -------------------------------------------------------------
# TEST: greetings
# A simple MPI program
# -------------------------------------------------------------
add_executable(test_mapper test/test_mapper.cpp ${target_objects})
target_link_libraries(test_mapper ${target_libraries})

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  bus_vector_map.hpp
  full_map.hpp
  DESTINATION include/gridpack/mapper
)

