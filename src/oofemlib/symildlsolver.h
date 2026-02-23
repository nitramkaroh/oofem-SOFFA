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

#pragma once

#include "sparselinsystemnm.h"
#include "symildlwrapper.h"

#define _IFT_SymildlSolver_Name "symildl"
#define _IFT_SymildlSolver_updateD "dcontrol"
#define _IFT_SymildlSolver_controlledIncrementReductionFactor "controlledincrementreductionfactor"
#define _IFT_SymildlSolver_postcontrolledIncrementReductionFactor "postcontrolledincrementreductionfactor"
#define _IFT_SymildlSolver_flipMultiplicationFactor "flipfactor"

namespace oofem {
/**
 * Implements the solution of linear system of equation in the form @f$ A\cdot x=b @f$ using solvers
 * from the sym-ildl library.
 */
class OOFEM_EXPORT SymildlSolver : public SparseLinearSystemNM
{

  private:
  bool dControl = false;
  double controlledIncrementReductionFactor = 1.0; //increments are multiplied by this after dcontrol happens
  double postcontrolledIncrementReductionFactor = 1.0; //increments are multiplied by this after dcontrol does not happen
  double flipMultiplicationFactor = 1.0; //controlled values in D are multiplied by this number

  bool DrecentlyModified = false; //internal variable remembers whether D was modified or not at the last call of "solve"

  public:
  /**
   * Constructor.
   * @param d Domain which solver belongs to.
   * @param m Engineering model which solver belongs to.
   */
  SymildlSolver( Domain *d, EngngModel *m );

  virtual ~SymildlSolver();

  /// Initializes receiver from given record.
  void initializeFrom( InputRecord &ir ) override;

  ConvergedReason solve( SparseMtrx &A, FloatArray &b, FloatArray &x ) override;

  const char *giveClassName() const override { return "SymildlSolver"; }
  LinSystSolverType giveLinSystSolverType() const override { return ST_Symildl; }
  SparseMtrxType giveRecommendedMatrix( bool symmetric ) const override { return SMT_SymildlMtrx; }

  private:

  /// <summary>
  /// Modifies the block diagonal matrix D of the LDL^T decomposition by flipping its negative diagonal members and
  ///   also the negative eigenvalues of its 2x2 diagonal blocks
  /// This forces the solution towards stable solution branches 
  /// </summary>
  std::pair<int, int> modifyD( block_diag_matrix<double> &D ) const;

  bool makeBlockPositive( double &a, double &b, double &c, double tol, double eps = 1.e-12 ) const;

  std::vector<double> SymildlSolver::solveModifiedSystem( lilc_matrix<double> &A, lilc_matrix<double> &L, block_diag_matrix<double> &D, std::vector<int> &perm, std::vector<double> &rhs ) const;

};
} // end namespace oofem
