// -------------------------------------------------------------
/**
 * @file   network.hpp
 * @author Kevin A. Glass
 * @date   Fri Apr  19 13:36:28 2013
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  19, 2013 by Kevin A. Glass
// Last Change:
// -------------------------------------------------------------

#ifndef _network_hpp_
#define _network_hpp_

namespace gridpack {
namespace network {

class MatrixInterface;
#include <iostream>
#include <vector>

// -------------------------------------------------------------
//  class Network
// -------------------------------------------------------------
/**
 * 
 * To be safe, these should be used simultaneously on all processes.  
 */

class Network {
public:
  Network(void){};
  virtual ~Network(void){};

  // this constructs the objects required by the analysis
  void networkMap(MapConstructor * map) = 0;

protected:
private:
};

} // namespace math
} // namespace gridpack

#endif
