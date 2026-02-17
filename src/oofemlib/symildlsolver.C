/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "symildlsolver.h"
#include "symildlmtrx.h"

#include "symildlwrapper.h" //symildl stuff through sanitized wrapper

#include "engngm.h"
#include "floatarray.h"
#include "verbose.h"
#include "timer.h"
#include "error.h"
#include "classfactory.h"
//

namespace oofem {
REGISTER_SparseLinSolver( SymildlSolver, ST_Symildl );

SymildlSolver ::SymildlSolver( Domain *d, EngngModel *m ) :
    SparseLinearSystemNM( d, m ) 
{
}

SymildlSolver ::~SymildlSolver() {}

void SymildlSolver ::initializeFrom( InputRecord &ir )
{
    this->dControl = ir.hasField(_IFT_SymildlSolver_updateD);
}
  

ConvergedReason SymildlSolver ::solve( SparseMtrx &A, FloatArray &b, FloatArray &x )
{
    int neqs = b.giveSize(); // Number of equations

    SymildlMtrx *As = dynamic_cast<SymildlMtrx *>( &A );
    if (As == nullptr)
      OOFEM_ERROR("Wrong matrix type for sym-ildl solver");

    CSCMatrix A_data = As->giveMatrix();

    // construct solver
    symildl::solver<double> solver;
    solver.set_message_level("none");
    solver.set_pivot("bunch");
    solver.set_equil("none");
    solver.set_reorder_scheme("none");
    solver.save_sol = false;
    solver.set_solver("full");
    solver.load( A_data.Ap, A_data.Ai, A_data.Ax );

    //debug
    //lilc_matrix Asymildl = solver.A;
    //Asymildl.save("debugmatrix.mtx");
    //As->printYourself();

    std::vector<double> rhs( b.begin(), b.end() );
    if ( !dControl ) {
      // no control of D is required, we can use the full sym-ildl solving with rhs
      solver.set_rhs( std::move( rhs ) );
    }

    // time the solution
    Timer timer;
    timer.startTimer();

    solver.solve( 10, 1.e-6, 1 );
    // /\ this line LDL factorizes the matrix and, if rhs has been set, also solves the system

    std::vector<double> sol;
    int nchanges = 0;
    if ( !dControl ) {
      // extract solution from solver
      sol = std::vector<double>(solver.sol_vec.begin(), solver.sol_vec.end());
    }else{
      // extract L, D, P from solver, modify D, then solve for solution
      lilc_matrix<double> A = solver.A;
      lilc_matrix<double> L = solver.L;
      block_diag_matrix<double> D = solver.D;
      std::vector<int> perm = solver.perm;

      int nchanges = this->modifyD(D);
      sol = this->solveModifiedSystem(A,L,D,perm,rhs);

    }


    // Copy/move values to FloatArray x
    x = FloatArray( sol.begin(), sol.end() );

    timer.stopTimer();
    if ( !dControl ) {
      OOFEM_LOG_INFO( "SymildlSolver:  User time consumed by solution: %.2fs\n", timer.getUtime() );
    } else {
      OOFEM_LOG_INFO( "SymildlSolver:  User time consumed by solution: %.2fs, changed %i diagonal blocks\n", timer.getUtime(), nchanges ); 
    }

    return CR_CONVERGED;

}

std::vector<double> SymildlSolver::solveModifiedSystem( lilc_matrix<double> &A, lilc_matrix<double> &L, block_diag_matrix<double> &D, std::vector<int> &perm, std::vector<double> &rhs ) const
{
  std::vector<double> sol_vec;
  //using the implementation from symildl's full solver
  //this could be improved if other faster solvers from symildl were implemented here 

  // we've permuted and equilibrated the matrix, so we gotta apply
  // the same permutation and equilibration to the right hand side.
  // i.e. rhs = P'S*rhs
  // 0. apply S
  for ( int i = 0; i < A.n_cols(); i++ ) {
    rhs[i] = A.S[i] * rhs[i];
  }

  // 1. apply P' (takes rhs[perm[i]] to rhs[i], i.e. inverse of perm,
  //    where perm takes i to perm[i])
  std::vector<double> tmp( A.n_cols() );
  for ( int i = 0; i < A.n_cols(); i++ ) {
    tmp[i] = rhs[perm[i]];
  }
  rhs = tmp;

  sol_vec.resize( A.n_cols(), 0 );
  // MINRES uses the preconditioned solver that
  // splits the block D into |D|^(1/2).
  // For the full solver we'll just solve D directly.
  L.backsolve( rhs, sol_vec );
  D.solve( sol_vec, tmp );
  L.forwardsolve( tmp, sol_vec );

  // 1. apply P
  for ( int i = 0; i < A.n_cols(); i++ ) {
    tmp[perm[i]] = sol_vec[i];
  }
  sol_vec = tmp;

  // 2. apply S
  for ( int i = 0; i < A.n_cols(); i++ ) {
    sol_vec[i] = A.S[i] * sol_vec[i];
  }

  return sol_vec;

}

int SymildlSolver::modifyD( block_diag_matrix<double> &D ) const
{
  // iterate through D
  double tol = 1e-12;
  int changes = 0;
  for ( int ii = 0; ii < D.n_cols(); ii++ ) {
    int blocksize = D.block_size( ii );
    if ( blocksize == 1 ) {
      // this is a standard diagonal entry
      double a = D[ii];
      if ( a < tol) {
        a *= -1.;
        changes++;
      }
    } else if ( blocksize == 2 ) {
      // this is a top left corner of a block diagonal entry
      double a = D[ii];
      double b = D[ii+1];
      double c = D.off_diagonal(ii);
      bool block_changed = makeBlockPositive(a,b,c,tol);
      if ( block_changed ) changes++;
    } else if ( blocksize == -2 ) {
      // this is a bottom right corder of a block diagonal entry
      // this entry has already been processed
      continue;
    }
  }
  return changes;
}

bool SymildlSolver::makeBlockPositive( double &a, double &b, double &c, double eps = 1e-12 ) const
{
  // ---- 1. Compute eigenvalues ----
  double trace_half = 0.5 * ( a + c );
  double diff_half = 0.5 * ( a - c );

  double rad = std::sqrt( diff_half * diff_half + b * b );

  double lambda1 = trace_half + rad;
  double lambda2 = trace_half - rad;

  // ---- 2. If already positive definite, do nothing ----
  if ( lambda1 > eps && lambda2 > eps )
    return false;

  // ---- 3. Clamp eigenvalues ----
  lambda1 = std::max( lambda1, fabs( lambda1 ) );
  lambda2 = std::max( lambda2, fabs( lambda2 ) );

  // ---- 4. Compute eigenvector for lambda1 ----
  double v1x, v1y;

  if ( std::abs( b ) > std::abs( a - lambda1 ) ) {
    v1x = 1.0;
    v1y = ( lambda1 - a ) / b;
  } else {
    v1x = b / ( lambda1 - c );
    v1y = 1.0;
  }

  // Normalize
  double norm = std::sqrt( v1x * v1x + v1y * v1y );
  v1x /= norm;
  v1y /= norm;

  // Orthogonal second eigenvector
  double v2x = -v1y;
  double v2y = v1x;

  // ---- 5. Reconstruct D = Q Λ Qᵀ ----
  a = lambda1 * v1x * v1x + lambda2 * v2x * v2x;
  b = lambda1 * v1x * v1y + lambda2 * v2x * v2y;
  c = lambda1 * v1y * v1y + lambda2 * v2y * v2y;

  return true;
}

} // end namespace oofem

