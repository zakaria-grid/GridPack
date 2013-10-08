// Emacs Mode Line: -*- Mode:c++;-*-
// -------------------------------------------------------------
/**
 * @file   linear_solver_implementation.hpp
 * @author William A. Perkins
 * @date   2013-10-03 13:41:55 d3g096
 * 
 * @brief  
 * 
 * 
 */

// -------------------------------------------------------------
// -------------------------------------------------------------
// Created April  1, 2013 by William A. Perkins
// Last Change: Thu Jun  3 06:45:08 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#ifndef _linear_solver_implementation_hpp_
#define _linear_solver_implementation_hpp_

#include <boost/scoped_ptr.hpp>
#include <gridpack/math/matrix.hpp>
#include <gridpack/parallel/distributed.hpp>
#include <gridpack/utilities/uncopyable.hpp>
#include <gridpack/configuration/configurable.hpp>

namespace gridpack {
namespace math {

// -------------------------------------------------------------
//  class LinearSolverImplementation
// -------------------------------------------------------------
class LinearSolverImplementation 
  : public parallel::Distributed,
    public utility::Configurable,
    private utility::Uncopyable
{
public:
  
  /// Default constructor.
  LinearSolverImplementation(const Matrix& A);

  /// Destructor
  ~LinearSolverImplementation(void);
  
  /// Solve w/ the specified RHS, put result in specified vector
  void solve(const Vector& b, Vector& x) const
  {
    this->p_solve(b, x);
  }

  /// Use different coefficient matrix (or A w/ new values)
  void set_matrix(const Matrix& A)
  {
    p_A.reset(A.clone());
    this->p_set_matrix(A);
  }
  /// Allow visits by implemetation visitor
  void accept(ImplementationVisitor& visitor)
  {
    this->p_accept(visitor);
  }

  /// Allow visits by const implemetation visitor
  void accept(ConstImplementationVisitor& visitor) const
  {
    this->p_accept(visitor);
  }

protected:

  /// The coefficient matrix (may not need to remember)
  boost::scoped_ptr<Matrix> p_A;
  
  /// Solve w/ the specified RHS, put result in specified vector
  /** 
   * Can be called repeatedly with different @c b and @c x vectors
   * 
   * @param b Vector containing right hand side of linear system
   * @param x Vector containing initial solution estimate; final
   * solution is put into this.
   */
  virtual void p_solve(const Vector& b, Vector& x) const = 0;

  /// Use different coefficient matrix (or A w/ new values) (specialized)
  /** 
   * The matrix must have the same parallel environment and same
   * nonzero layout as that used to construct the solver instance
   * 
   * @param A new coefficient matrix 
   */
  virtual void p_set_matrix(const Matrix& A) = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ImplementationVisitor& visitor) = 0;

  /// Allow visits by implementation visitors
  virtual void p_accept(ConstImplementationVisitor& visitor) const = 0;

};

} // namespace math
} // namespace gridpack



#endif
