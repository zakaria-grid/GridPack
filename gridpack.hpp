// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   gridpack.hpp
 * @author William A. Perkins
 * @date   2013-11-06 08:35:18 d3g096
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------

#ifndef _gridpack_hpp_
#define _gridpack_hpp_

#include <gridpack/utilities/exception.hpp>
#include <gridpack/timer/coarse_timer.hpp>
#include <gridpack/serial_io/serial_io.hpp>
#include <gridpack/math/math.hpp>
#include <gridpack/mapper/bus_vector_map.hpp>
#include <gridpack/mapper/full_map.hpp>
#include <gridpack/parser/PTI23_parser.hpp>

namespace gridpack {

/// Do whatever is necessary to initialize the GridPACK library.
extern void Initialize(int argc, char **argv);

/// Do whatever is necessary to shut down the GridPACK library.
extern void Finalize(void);

} // namespace gridpack

#endif
