// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_graph_partitioner_impl.hpp
 * @author William A. Perkins
 * @date   2014-02-06 09:29:41 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _petsc_graph_partitioner_impl_hpp_
#define _petsc_graph_partitioner_impl_hpp_

#include "graph_partitioner_implementation.hpp"

namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class PETScGraphPartitionerImpl
// -------------------------------------------------------------
class PETScGraphPartitionerImpl 
  : public GraphPartitionerImplementation
{
public:

  /// Default constructor.
  PETScGraphPartitionerImpl(const parallel::Communicator& comm);

  /// Construct w/ known local sizes (guesses to size containers, maybe)
  PETScGraphPartitionerImpl(const parallel::Communicator& comm,
                            const int& local_nodes, const int& local_edges);
  

  /// Destructor
  ~PETScGraphPartitionerImpl(void);

protected:

  /// Partition the graph (specialized)
  void p_partition(void);

};


} // namespace network
} // namespace gridpack


#endif
