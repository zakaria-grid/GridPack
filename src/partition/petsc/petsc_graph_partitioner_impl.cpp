// -------------------------------------------------------------
/*
 *     Copyright (c) 2013 Battelle Memorial Institute
 *     Licensed under modified BSD License. A copy of this license can be found
 *     in the LICENSE file in the top level directory of this distribution.
 */
// -------------------------------------------------------------
/**
 * @file   petsc_graph_partitioner_impl.cpp
 * @author William A. Perkins
 * @date   2014-02-14 11:51:41 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------


#include <petscmat.h>
#include <boost/mpi/collectives.hpp>

#include "gridpack/math/petsc/petsc_exception.hpp"
#include "petsc/petsc_graph_partitioner_impl.hpp"


namespace gridpack {
namespace network {

// -------------------------------------------------------------
//  class PETScGraphPartitionerImpl
// -------------------------------------------------------------

// -------------------------------------------------------------
// PETScGraphPartitionerImpl:: constructors / destructor
// -------------------------------------------------------------
PETScGraphPartitionerImpl::PETScGraphPartitionerImpl(const parallel::Communicator& comm)
  : GraphPartitionerImplementation(comm)
{
  
}

PETScGraphPartitionerImpl::PETScGraphPartitionerImpl(const parallel::Communicator& comm,
                                                     const int& local_nodes, 
                                                     const int& local_edges)
  : GraphPartitionerImplementation(comm, local_nodes, local_nodes)
{
  
}

PETScGraphPartitionerImpl::~PETScGraphPartitionerImpl(void)
{
}

// -------------------------------------------------------------
// PETScGraphPartitionerImpl::p_partition
// -------------------------------------------------------------
void
PETScGraphPartitionerImpl::p_partition(void)
{
  int me(communicator().rank());
  int locnodes(p_adjacency_list.nodes());
  int locedges(p_adjacency_list.edges());

  std::vector<int> pnodeoffset(communicator().size(), 0);

  all_gather(communicator().getCommunicator(), locnodes, pnodeoffset);
  int allnodes(0);
  for (size_t p = 0; p < pnodeoffset.size(); ++p) {
    int tmp(pnodeoffset[p]);
    pnodeoffset[p] = allnodes;
    allnodes += tmp;
  }

  Mat adjacency;
  PetscErrorCode ierr(0);
  try {
    ierr = MatCreate(this->communicator(), &adjacency); CHKERRXX(ierr);
    ierr = MatSetSizes(adjacency, PETSC_DECIDE, PETSC_DECIDE, allnodes, allnodes); CHKERRXX(ierr);
    if (this->communicator().size() == 1) {
      ierr = MatSetType(adjacency, MATSEQAIJ); CHKERRXX(ierr);
    } else {
      ierr = MatSetType(adjacency, MATMPIAIJ); CHKERRXX(ierr);
    }
    ierr = MatSetFromOptions(adjacency); CHKERRXX(ierr);
    ierr = MatSetUp(adjacency); CHKERRXX(ierr);

    for (int l = 0; l < locedges; ++l) {
      AdjacencyList::Index n1, n2;
      PetscScalar one(1.0);
      p_adjacency_list.edge(l, n1, n2);
      ierr = MatSetValue(adjacency, n1, n2, one, INSERT_VALUES); CHKERRXX(ierr);
      ierr = MatSetValue(adjacency, n2, n1, one, INSERT_VALUES); CHKERRXX(ierr);
    }

    ierr = MatAssemblyBegin(adjacency, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);
    ierr = MatAssemblyEnd(adjacency, MAT_FINAL_ASSEMBLY); CHKERRXX(ierr);

    MatPartitioning part;

    ierr = MatPartitioningCreate(this->communicator(), &part); CHKERRXX(ierr);
    ierr = MatPartitioningSetAdjacency(part, adjacency); CHKERRXX(ierr);
    ierr = MatPartitioningSetFromOptions(part); CHKERRXX(ierr);

    IS gblis;

    ierr = MatPartitioningApply(part, &gblis); CHKERRXX(ierr);
    ierr = ISView(gblis, PETSC_VIEWER_STDOUT_(this->communicator())); CHKERRXX(ierr);

    IS locis;
    
    ierr = ISAllGather(gblis, &locis);  CHKERRXX(ierr);

    const PetscInt *iptr;
    ierr = ISGetIndices(locis, &iptr); CHKERRXX(ierr);

    p_node_destinations.clear();
    p_node_destinations.resize(locnodes);
    for (int n = 0; n < locnodes; ++n) {
      AdjacencyList::Index gidx(p_adjacency_list.node_index(n));
      p_node_destinations[n] = iptr[gidx];
    }    
    ierr = ISRestoreIndices(locis, &iptr); CHKERRXX(ierr);

    ierr = ISDestroy(&locis); CHKERRXX(ierr);
    ierr = ISDestroy(&gblis); CHKERRXX(ierr);
    ierr = MatPartitioningDestroy(&part); CHKERRXX(ierr);
    ierr = MatDestroy(&adjacency); CHKERRXX(ierr);
    
  } catch (const PETSc::Exception& e) {
    throw math::PETScException(ierr, e);
  }
}


} // namespace network
} // namespace gridpack
