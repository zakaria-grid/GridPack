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
# Last Change: 2014-03-24 11:44:46 d3g096
# -------------------------------------------------------------

set(target_libraries
    ${GA_LIBRARIES}
    ${Boost_LIBRARIES}
    ${MPI_CXX_LIBRARIES}
)

set(target_objects 
  $<TARGET_OBJECTS:gridpack_parallel>
  $<TARGET_OBJECTS:gridpack_timer>
  $<TARGET_OBJECTS:gridpack_configuration>
  $<TARGET_OBJECTS:gridpack_components>
)

include_directories(BEFORE ${CMAKE_CURRENT_SOURCE_DIR})

if (GA_FOUND)
  include_directories(AFTER ${GA_INCLUDE_DIRS})
endif()


# -------------------------------------------------------------
# TEST: test_parser
# -------------------------------------------------------------
add_executable(test_parser 
  test/PTI23_test.cpp
  ${target_objects}
)
target_link_libraries(test_parser ${target_libraries})

add_custom_target(test_parser_input

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${CMAKE_CURRENT_SOURCE_DIR}/test/test.xml
  ${CMAKE_CURRENT_BINARY_DIR}

  COMMAND ${CMAKE_COMMAND} -E copy 
  ${GridPACK_SOURCE_DIR}/applications/powerflow/IEEE14.raw
  ${CMAKE_CURRENT_BINARY_DIR}

  DEPENDS
  ${CMAKE_CURRENT_SOURCE_DIR}/test/test.xml
  ${GridPACK_SOURCE_DIR}/applications/powerflow/IEEE14.raw
)
add_dependencies(test_parser test_parser_input)

add_test(test_parser_serial test_parser)

# -------------------------------------------------------------
# installation
# -------------------------------------------------------------
install(FILES 
  dictionary.hpp
  PTI23_parser.hpp
  DESTINATION include/gridpack/parser
)
