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
# Last Change: 2014-03-17 08:13:41 d3g096
# -------------------------------------------------------------

set(target_libraries
    gridpack_parallel
    ${Boost_LIBRARIES}
    ${GA_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(AFTER ${GA_INCLUDE_DIRS})

add_library(gridpack_dsmatrix_components
  ds_components.cpp
)

# -------------------------------------------------------------
# component serialization tests
# -------------------------------------------------------------
target_link_libraries(gridpack_dsmatrix_components ${target_libraries})

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  ds_components.hpp  
  DESTINATION include/gridpack/applications/components/ds_matrix
)

install(TARGETS 
  gridpack_dsmatrix_components
  DESTINATION lib
)

