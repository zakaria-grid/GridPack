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
# Last Change: 2013-11-05 09:06:01 d3g096
# -------------------------------------------------------------

add_library(gridpack_factory OBJECT
  base_factory.cpp
)

include_directories(AFTER ${GA_INCLUDE_DIRS})


