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

#include "eigensolverstability.h"
#include "eigensolver.h"
#include "eigenmtrx.h"
// #include <iostream>

#include "compcol.h"
#include "symcompcol.h"
#include "engngm.h"
#include "floatarray.h"
#include "verbose.h"
#include "timer.h"
#include "error.h"
#include "classfactory.h"

//
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

namespace oofem {
REGISTER_SparseLinSolver( EigenSolverStability, ST_EigenStability );

EigenSolverStability ::EigenSolverStability( Domain *d, EngngModel *m ) :
    EigenSolver( d, m )
{
}

EigenSolverStability ::~EigenSolverStability() {}

void EigenSolverStability ::initializeFrom( InputRecord &ir )
{
    EigenSolver ::initializeFrom( ir );
}


void EigenSolverStability::solveLDLT( Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b, Eigen::VectorXd &x )
{
    SimplicialLDLTderived<Eigen::SparseMatrix<double> > A_factorization( A );

    double minLam = 0.;
    bool negEig   = false;

     //negEig = A_factorization.updateD2( minLam, !isSingular );
    negEig = A_factorization.updateD( minLam, true );

    OOFEM_LOG_INFO( "lam_min = %.4f\n", minLam );

    if ( negEig ) {
        isSingular = true;
        OOFEM_LOG_INFO( "Negative eigenvalue\n" );
    }

    x = A_factorization.solve( b ); // Solve the system
}

} // end namespace oofem
