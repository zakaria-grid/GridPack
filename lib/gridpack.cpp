// -------------------------------------------------------------
/**
 * @file   gridpack.cpp
 * @author William A. Perkins
 * @date   2013-11-06 08:41:41 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#include <ga++.h>
#include "gridpack/gridpack.hpp"


namespace gridpack {

// -------------------------------------------------------------
// Initialize
// -------------------------------------------------------------
void
Initialize(int argc, char **argv) 
{
  int ierr = MPI_Init(&argc, &argv);
  gridpack::math::Initialize();
  GA_Initialize();
  MA_init(MT_C_CHAR, 1024*1024, 1024*1024);
}

// -------------------------------------------------------------
// Finalize
// -------------------------------------------------------------
void
Finalize(void)
{
  GA_Terminate();
  gridpack::math::Finalize();
  int ierr = MPI_Finalize();
}


} // namespace gridpack
