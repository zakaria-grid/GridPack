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
# Last Change: 2013-11-05 09:13:43 d3g096
# -------------------------------------------------------------

set(target_libraries ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

add_library(gridpack_components OBJECT
  base_component.cpp
  data_collection.cpp
)

set(target_objects 
  $<TARGET_OBJECTS:gridpack_parallel>
  $<TARGET_OBJECTS:gridpack_components>
)
# -------------------------------------------------------------
# component serialization tests
# -------------------------------------------------------------
add_executable(component_serialization_test 
  test/component_serialization.cpp 
  ${target_objects}
)
target_link_libraries(component_serialization_test ${target_libraries})
gridpack_add_unit_test(component_serialization component_serialization_test)


