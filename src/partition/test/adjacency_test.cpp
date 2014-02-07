/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   partition_test.cpp
 * @author William A. Perkins
 * @date   2014-02-07 10:07:27 d3g096
 * 
 * @brief  Unit test suite for various partition classes.
 * 
 * @test
 */
// -------------------------------------------------------------

#include <ctime>
#include <iostream>
#include <iterator>
#include <ga++.h>
#include <boost/mpi/collectives.hpp>
#include <boost/serialization/vector.hpp>
#include "gridpack/utilities/exception.hpp"
#include "gridpack/parallel/parallel.hpp"

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#include <boost/test/included/unit_test.hpp>

#include "simple_adjacency.hpp"

BOOST_AUTO_TEST_SUITE( PartitionTest )


// -------------------------------------------------------------
// adjacency unit test
// -------------------------------------------------------------
/// 
/** 
 * @test
 * 
 * 
 */
BOOST_AUTO_TEST_CASE( adjacency )
{
  gridpack::parallel::Communicator world;
  const int global_nodes(4*world.size());
  const int global_edges(global_nodes - 1);
  
  using gridpack::network::AdjacencyList;

  std::auto_ptr<AdjacencyList> 
    adlist(simple_adjacency_list(world, global_nodes));

  for (int p = 0; p < world.size(); ++p) {
    if (world.rank() == p) {
      AdjacencyList::IndexVector nbr;
      std::cout << p << ": Results" << std::endl;
      for (size_t i = 0; i < adlist->nodes(); ++i) {
        adlist->node_neighbors(i, nbr);
        std::cout << p << ": node " << adlist->node_index(i) << ": ";
        std::copy(nbr.begin(), nbr.end(),
                  std::ostream_iterator<AdjacencyList::Index>(std::cout, ","));
        std::cout << std::endl;
      }
    }
    world.barrier();
  }

  for (size_t i = 0; i < adlist->nodes(); ++i) {
    AdjacencyList::Index node(adlist->node_index(i));
    AdjacencyList::IndexVector nbr;
    adlist->node_neighbors(i, nbr);

    // all nodes should connect to at least one other, but no more than two
    BOOST_CHECK(nbr.size() >= 1 && nbr.size() <= 2);

    // node n should connect to node n-1, except for 0
    AdjacencyList::IndexVector::iterator p;
    if (node > 0) {
      p = std::find(nbr.begin(), nbr.end(), node - 1);
      BOOST_CHECK(p != nbr.end());
    }

    // node n should connect to node n+1 except for the last
    if (node < global_nodes - 1) {
      p = std::find(nbr.begin(), nbr.end(), node + 1);
      BOOST_CHECK(p != nbr.end());
    }
  }
    
}

BOOST_AUTO_TEST_SUITE_END()


// -------------------------------------------------------------
// init_function
// -------------------------------------------------------------
/**
 * @test
 * 
 */
bool init_function()
{
  return true;
}

// -------------------------------------------------------------
//  Main Program
// -------------------------------------------------------------
/**
 * @test
 * 
 */
int
main(int argc, char **argv)
{
  gridpack::parallel::Environment env(argc, argv);
  GA_Initialize();
  MA_init(MT_C_CHAR, 1024*1024, 1024*1024);
  int result = ::boost::unit_test::unit_test_main( &init_function, argc, argv );
  GA_Terminate();
}




