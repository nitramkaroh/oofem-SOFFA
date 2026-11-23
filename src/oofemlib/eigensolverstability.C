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
//#include "unsupported/Eigen/src/SparseExtra/MarketIO.h"

// for spectra
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <include/Spectra/GenEigsSolver.h>
#include <include/Spectra/SymEigsSolver.h>
#include <include/Spectra/MatOp/SparseGenMatProd.h>
#include <include/Spectra/MatOp/SparseSymMatProd.h>

#include <include/Spectra/SymEigsShiftSolver.h>
#include <include/Spectra/MatOp/SparseSymShiftSolve.h>

#include <exception>

namespace oofem {
REGISTER_SparseLinSolver( EigenSolverStability, ST_EigenStability );

EigenSolverStability ::EigenSolverStability( Domain *d, EngngModel *m ) :
    EigenSolver( d, m ),
    BifurcationInterface()
{}

EigenSolverStability ::~EigenSolverStability() {}

void EigenSolverStability ::initializeFrom( InputRecord &ir )
{
    EigenSolver ::initializeFrom( ir );
    //BifurcationInterface ::initializeFrom( ir );
}

ConvergedReason EigenSolverStability ::solve( SparseMtrx &A, FloatArray &b, FloatArray &x )
{
    //Timer timer;
    //timer.startTimer();
    this->setMethod( "lu" );
    EigenSolver::solve( A,b,x);

    this->performBifurcation( A, b, x );

    //timer.stopTimer();
    //OOFEM_LOG_INFO( "EigenSolver:  User time consumed by solution: %.2fs\n", timer.getUtime() );
    return CR_CONVERGED;
}



bool EigenSolverStability::checkPD( SparseMtrx &A )
{
    EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &A );
    if ( !Ae ) {
        OOFEM_ERROR( "definiteness check requires an EigenMtrx, got %s", A.giveClassName() );
    }

    const int n = Ae->giveNumberOfRows();
    if ( n == 0 ) {
        // no equations: there is no direction in which the matrix could fail to be positive
        this->storePDStatus( true, 0, 0., 0 );
        return true;
    }

    SimplicialLDLTderived<Eigen::SparseMatrix<double> >& A_factorization = Ae->giveLDLTFactorization();

    // The success of the factorization is checked before D is used: after a breakdown the
    // entries of D past the zero pivot are undefined (stale or uninitialized memory).
    double minLam  = 0.;
    int numNegEigs = 0;
    bool isPD;
    if ( A_factorization.giveInertia( minLam, numNegEigs ) ) {
        isPD = numNegEigs == 0; // all pivots finite and nonzero: positive definite iff none is negative
    } else {
        // LDL^T broke down on a zero pivot (or produced non-finite pivots). A positive definite matrix
        // has no zero pivot, so the matrix is reported as not positive definite; the number of negative
        // pivots counts only those computed before the breakdown and is a lower bound.
        OOFEM_WARNING( "LDL^T factorization of the checked matrix (size %d) did not succeed: it is reported as not positive definite, the number of negative pivots is only a lower bound", n );
        isPD = false;
    }

    // Keep the inertia of the factorization, so that the definiteness of the converged
    // solution can be reported in the .out file after the step has been completed.
    // The values recorded are those the verdict returned to the solver is based on.
    this->storePDStatus( isPD, numNegEigs, minLam, n );

    return isPD;
}

int EigenSolverStability::computeEigenValues( SparseMtrx &Ag, FloatArray &evaluesFA )
{
    EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &Ag );
    if ( !Ae ) {
        OOFEM_ERROR( "eigenvalue computation requires an EigenMtrx, got %s", Ag.giveClassName() );
    }

    // A local copy is kept, the Spectra operator below only holds a reference to it
    Eigen::SparseMatrix<double> A = Ae->giveMatrix();
    int neigs = 1, maxTrials = 5, ncvMult = 2, ncv = 30;
    int nconv = 0, count = 0;

    int Asize = (int)A.rows();
    if ( Asize < 2 ) {
        // Spectra needs more equations than requested eigenvalues; a 1x1 matrix is its own eigenvalue
        if ( Asize == 1 ) {
            evaluesFA.resize( 1 );
            evaluesFA.at( 1 ) = A.coeff( 0, 0 );
            return 1;
        }
        return 0;
    }

    //Spectra::SparseSymShiftSolve<double> op( A );
    Spectra::SparseSymMatProd<double> op( A );

    try {
        while ( nconv == 0 && count < maxTrials ) {
            int ncv2 = std::min( ncv, Asize );
            //Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double> > eigs( op, neigs, ncv2, 1e-8 );
            Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double> > eigs( op, neigs, ncv2 );
            eigs.init();

            nconv = eigs.compute( Spectra::SortRule::SmallestAlge );
            if ( eigs.info() == Spectra::CompInfo::Successful ) {
                Eigen::VectorXd evalues = eigs.eigenvalues();
                evaluesFA               = FloatArray( evalues.begin(), evalues.end() );
            } else {
                count++;
                OOFEM_LOG_INFO( "eigenvalues not found\n" );
                if ( ncv2 == Asize ) {
                    // ncv cannot grow any more and init() always starts from the same vector: a retry would repeat this trial
                    break;
                }
                ncv = (int)ncvMult * ncv;
                OOFEM_LOG_INFO( "ncv =  %i\n", ncv );
            }
        }
    } catch ( const std::exception &e ) {
        // e.g. sizes Spectra does not accept: reported as not computed instead of stopping the analysis
        OOFEM_WARNING( "eigenvalue computation failed: %s", e.what() );
        return 0;
    }

    return nconv;
}

//

int EigenSolverStability ::computeEigenValuesVectors( SparseMtrx &Ag, FloatArray &evaluesFA, FloatMatrix &evectorsFM )
{
    EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &Ag );
    if ( !Ae ) {
        OOFEM_ERROR( "eigenpair computation requires an EigenMtrx, got %s", Ag.giveClassName() );
    }

    // A local copy is kept, the Spectra operator below only holds a reference to it
    Eigen::SparseMatrix<double> A = Ae->giveMatrix();
    int neigs = 1, maxTrials = 5, ncvMult = 2, ncv = 200; // Predefined values, should be in input file
    int nconv = 0, count = 0;

    int Asize = (int)A.rows();
    if ( Asize < 2 ) {
        // Spectra needs more equations than requested eigenvalues; a 1x1 matrix is its own eigenpair
        if ( Asize == 1 ) {
            evaluesFA.resize( 1 );
            evaluesFA.at( 1 ) = A.coeff( 0, 0 );
            evectorsFM.resize( 1, 1 );
            evectorsFM.at( 1, 1 ) = 1.;
            return 1;
        }
        return 0;
    }

    //Spectra::SparseSymMatProd<double> op( A );
    Spectra::SparseSymShiftSolve<double> op( A );

    evectorsFM.resize( Asize, neigs );
    evaluesFA.resize( neigs );

    try {
        while ( nconv == 0 && count < maxTrials ) { // if ncv is too small, increase it ti find solution
            int ncv2 = std::min( ncv, Asize );
            //Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double> > eigs( op, neigs, ncv2 );
            Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double> > eigs( op, neigs, ncv2, 0.0 );

            eigs.init();

            nconv = eigs.compute( Spectra::SortRule::SmallestAlge ); // Compute eigenvalues
            if ( eigs.info() == Spectra::CompInfo::Successful ) { // If succesfull, retrieve results
                // allocate some objects
                this->setNcv0( ncv );
                Eigen::VectorXd evector;
                Eigen::MatrixXd evectors;
                FloatArray evectorFA;
                evectorFA.resize( Asize * neigs );

                Eigen::VectorXd evalues = eigs.eigenvalues();
                evaluesFA               = FloatArray( evalues.begin(), evalues.end() ); // eigenvalues in floatarray

                evectors  = eigs.eigenvectors(); // vectors in matrix
                evector   = Eigen::Map<Eigen::VectorXd>( evectors.data(), evectors.cols() * evectors.rows() ); // map to vector
                evectorFA = FloatArray( evector.begin(), evector.end() ); // eigenvectors in floatarray

                // eigenvectors in float matrix
                for ( size_t j = 0; j < neigs; j++ ) {
                    for ( size_t k = 0; k < evectors.rows(); k++ ) {
                        evectorsFM.at( k + 1, j + 1 ) = evectorFA.at( k + 1 + evectors.rows() * j ); // Eigenvectors to float matrix
                    }
                }
            } else { // if not succesfull
                count++;
                OOFEM_LOG_INFO( "eigenvalues not found\n" );
                if ( ncv2 == Asize ) {
                    // ncv cannot grow any more and init() always starts from the same vector: a retry would repeat this trial
                    break;
                }
                ncv = (int)ncvMult * ncv; // Increase ncv
                OOFEM_LOG_INFO( "ncv =  %i\n", ncv );
            }
        }
    } catch ( const std::exception &e ) {
        // e.g. the shift-invert factorization fails because the matrix is singular (sigma = 0):
        // reported as not computed instead of stopping the analysis
        OOFEM_WARNING( "eigenpair computation failed: %s", e.what() );
        return 0;
    }
    return nconv;
}

void EigenSolverStability::CholeskyUpdate( SparseMtrx &Ag, FloatArray &b, FloatArray &x )
{
    EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &Ag );
    if ( !Ae ) {
        OOFEM_ERROR( "Cholesky update requires an EigenMtrx, got %s", Ag.giveClassName() );
    }

    Eigen::VectorXd b_eig                                                = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>( b.givePointer(), b.giveSize() );
    double minLam                                                        = 0.;
    int numNegEigs                                                       = 0;
    SimplicialLDLTderived<Eigen::SparseMatrix<double> > &A_factorization = Ae->giveLDLTFactorization(); // compute factorization and get it
    bool complete                                                        = A_factorization.giveInertia( minLam, numNegEigs );
    OOFEM_LOG_INFO( "numNegEigs = %i\n", numNegEigs );
    if ( !complete ) {
        OOFEM_WARNING( "LDL^T factorization did not succeed, the solution is left unchanged by the Cholesky update" );
        return;
    }
    // Solve with |D| (the matrix made positive definite). The cached factorization is not modified,
    // so that it keeps describing the matrix itself, e.g. for a later checkPD of the same matrix.
    Eigen::VectorXd xeig = A_factorization.solveWithAbsolutePivots( b_eig ); // Solve the system
    x                    = FloatArray( xeig.begin(), xeig.end() );
    // this->setBifurcation( false );
}



//ConvergedReason EigenSolverStability ::solveBifurcation( SparseMtrx &A, FloatArray &b, FloatArray &x )
//{
//    int neqs = b.giveSize(); // Number of equations
//    EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &A );
//    Eigen::VectorXd x_eig; // Allocate vector of RHS
//
//    // Construct right hand side vetor
//    Eigen::VectorXd b_eig = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>( b.givePointer(), neqs );
//
//    Timer timer;
//    timer.startTimer();
//
//    // Create factorization
//    if ( method.compare( "ldlt" ) == 0 ) {
//        this->setBifurcation( true );
//        this->solveLDLT( *Ae, b_eig, x_eig );
//
//    } else {
//        OOFEM_LOG_INFO( "Wrong factorization used, bifucration cannot be performed\n" );
//    }
//
//    // Copy/move values to FloatArray x
//    x = FloatArray( x_eig.begin(), x_eig.end() );
//
//    timer.stopTimer();
//    OOFEM_LOG_INFO( "EigenSolver:  User time consumed by solution: %.2fs\n", timer.getUtime() );
//
//    return CR_CONVERGED;
//}



} // end namespace oofem
