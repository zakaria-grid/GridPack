/**
 * @file   linear_solver_implementation.hpp
 * @author Kevin A. Glass
 * @date   Mon Apr  19 13:51 2013
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  1, 2013 by Kevin A. Glass
// Last Change:
// -------------------------------------------------------------

#ifndef _pf_bus_component_hpp_
#define _pf_bus_component_hpp_

#include "gridpack/network/matrix_tnterface.hpp"
#include "gridpack/math/matrix.hpp"
#include "gridpack/math/vector.hpp"
#include "gridpack/utility/noncopyable.hpp"
#include "gridpack/utility/configurable.hpp"
#include "gridpack/parallel/distributed.hpp"

namespace gridpack {
namespace math {
class MapSize;
class MapData;

// -------------------------------------------------------------
//  class BranchComponent
// -------------------------------------------------------------

class PFBranchComponent : public PFComponent
{
public:
    PFBranchComponent(const reader::ParameterReader, int size) : PFComponent(size);
    virtual ~PFBranchComponent(void);

    virtual void connectInputBus(PFBusComponent * bus) {input = bus;};
    virtual void connectOutputBus(PFBusComponent * bus){output = bus;};
    void increment(MapSize * map);
    void mapData(MapData * map){
       ;
    }
protected:
    PFBusComponent    * input;
    PFBusComponent    * output;

    /*
     * BRANCH SPECIFIC DATA
     */
} // namespace math
} // namespace gridpack

#endif
