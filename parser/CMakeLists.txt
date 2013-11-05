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
# Last Change: 2013-11-05 13:34:30 d3g096
# -------------------------------------------------------------

set(target_libraries
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES})

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()


# -------------------------------------------------------------
# TEST: test_parser
# -------------------------------------------------------------
add_executable(test_parser test/PTI23_test.cpp)
target_link_libraries(test_parser ${target_libraries})

add_test(test_parser_serial test_parser)

