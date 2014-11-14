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
# Last Change: 2014-03-24 11:58:24 d3g096
# -------------------------------------------------------------

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(AFTER ${GA_INCLUDE_DIRS})

add_library(gridpack_dsmatrix_components OBJECT
  ds_components.cpp
)

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  ds_components.hpp  
  DESTINATION include/gridpack/applications/components/ds_matrix
)


