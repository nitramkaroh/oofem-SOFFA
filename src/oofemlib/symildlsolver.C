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

    if ( !dControl ) {
      // no control of D is required, we can use the full sym-ildl solving with rhs
      std::vector<double> rhs( b.begin(), b.end() );
      solver.set_rhs( std::move( rhs ) );
    }

    // time the solution
    Timer timer;
    timer.startTimer();

    solver.solve( 10, 1.e-6, 1 );
    // /\ this line LDL factorizes the matrix and, if rhs has been set, also solves the system

    std::vector<double> sol;

    if ( !dControl ) {
      // extract solution from solver
      sol = std::vector<double>(solver.sol_vec.begin(), solver.sol_vec.end());
    }else{
      // extract L, D, P from solver, modify D, then solve for solution
      OOFEM_ERROR("Dcontrol not supported yet");
    }


    // Copy/move values to FloatArray x
    x = FloatArray( sol.begin(), sol.end() );

    timer.stopTimer();
    OOFEM_LOG_INFO( "SymildlSolver:  User time consumed by solution: %.2fs\n", timer.getUtime() );

    return CR_CONVERGED;

}

} // end namespace oofem

