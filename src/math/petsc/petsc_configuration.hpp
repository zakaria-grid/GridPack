// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   petsc_configuration.hpp
 * @author William A. Perkins
 * @date   2013-10-02 11:35:12 d3g096
 * 
 * @brief Declaration of routines for handling PETSc options through
 * Configuration
 * 
 * 
 */
// -------------------------------------------------------------


#ifndef _petsc_configuration_hpp_
#define _petsc_configuration_hpp_

#include <string>
#include "gridpack/configuration/configuration.hpp"
#include "gridpack/parallel/parallel.hpp"

namespace gridpack {
namespace math {

/// Process any PETSc options in the configuration and get the option prefix to use with PETSc
extern std::string petscProcessOptions(const parallel::Communicator& comm,
                                       utility::Configuration::Cursor *props);

} // namespace math
} // namespace gridpack

#endif
