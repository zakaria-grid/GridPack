// -------------------------------------------------------------
/**
 * @file   pf_factory.cpp
 * @author Bruce Palmer
 * @date   2013-08-08 10:23:25 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <vector>
#include "boost/smart_ptr/shared_ptr.hpp"
#include "gridpack/network/base_network.hpp"
#include "gridpack/factory/base_factory.hpp"
#include "gridpack/applications/powerflow/pf_components.hpp"
#include "gridpack/applications/powerflow/pf_factory.hpp"

namespace gridpack {
namespace powerflow {

// Powerflow factory class implementations

/**
 * Basic constructor
 * @param network: network associated with factory
 */
PFFactory::PFFactory(PFFactory::NetworkPtr network)
  : gridpack::factory::BaseFactory<PFNetwork>(network)
{
  p_network = network;
}

/**
 * Basic destructor
 */
gridpack::powerflow::PFFactory::~PFFactory()
{
}

/**
 * Create the admittance (Y-Bus) matrix
 */
void gridpack::powerflow::PFFactory::setYBus(void)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;

  // Invoke setYBus method on all bus objects
  for (i=0; i<numBus; i++) {
    dynamic_cast<PFBus*>(p_network->getBus(i).get())->setYBus();
  }

  // Invoke setYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
    dynamic_cast<PFBranch*>(p_network->getBranch(i).get())->setYBus();
  }
}

/**
 * Get values from the admittance (Y-Bus) matrix
 */
#if 0
gridpack::ComplexType gridpack::powerflow::PFFactory::calMis(
    gridpack::math::Vector V,
    gridpack::math::Vector SBUS)
{
  int numBus = p_network->numBuses();
  int numBranch = p_network->numBranches();
  int i;
  gridpack::ComplexType ibus;
  gridpack::ComplexType mis;

  // MIS = V * conj (YBus * V) - SBUS

  // Invoke getYBus method on all bus objects
  /* for (i=0; i<numBus; i++) {
     ibus =
     (dynamic_cast<gridpack::powerflow::PFBus*>(p_network->getBus(i).get()))->getYBus()
   * V(i) ;
   mis(i) = V * conj(ibus
   }

  // Invoke getYBus method on all branch objects
  for (i=0; i<numBranch; i++) {
  (dynamic_cast<gridpack::powerflow::PFBranch*>(p_network->getBranch(i).get()))->getYBus();
  }*/
}
#endif
} // namespace powerflow
} // namespace gridpack
