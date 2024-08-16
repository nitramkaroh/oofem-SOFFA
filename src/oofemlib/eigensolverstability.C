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

// for spectra
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <include/Spectra/GenEigsSolver.h>
#include <include/Spectra/SymEigsSolver.h>
#include <include/Spectra/MatOp/SparseGenMatProd.h>
#include <include/Spectra/MatOp/SparseSymMatProd.h>

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
    Eigen::SparseMatrix<double> A_eig = Ae->giveMatrix();

    SimplicialLDLTderived<Eigen::SparseMatrix<double> >& A_factorization = Ae->giveLDLTFactorization();

    double minLam  = 0.;
    bool negEig;
    int numNegEigs = 0;
    negEig = A_factorization.updateD( minLam, false, numNegEigs );

    return !negEig;
}

//

int EigenSolverStability ::computeEigenValuesVectors( SparseMtrx &Ag, FloatArray &evaluesFA, FloatMatrix &evectorsFM )
{
    Eigen::SparseMatrix<double> A = dynamic_cast<EigenMtrx *>( &Ag )->giveMatrix();
    int neigs = 1, maxTrials = 5, ncvMult = 2, ncv = 200; // Predefined values, should be in input file
    int nconv = 0, count = 0;

    int Asize = (int)A.rows();
    Spectra::SparseSymMatProd<double> op( A );

    evectorsFM.resize( Asize, neigs );
    evaluesFA.resize( neigs );

    while ( nconv == 0 && count < maxTrials ) { // if ncv is too small, increase it ti find solution
        int ncv2 = std::min( ncv, Asize );
        Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double> > eigs( op, neigs, ncv2 );
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
            ncv = (int)ncvMult * ncv; // Increase ncv
            OOFEM_LOG_INFO( "eigenvalues not found\n" );
            OOFEM_LOG_INFO( "ncv =  %i\n", ncv );
        }
    }
    return nconv;
}

void EigenSolverStability::CholeskyUpdate( SparseMtrx &Ag, FloatArray &b, FloatArray &x )
{
    EigenMtrx *Ae                                                        = dynamic_cast<EigenMtrx *>( &Ag );
    Eigen::VectorXd b_eig                                                = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>( b.givePointer(), b.giveSize() );
    double minLam                                                        = 0.;
    int numNegEigs                                                       = 0;
    SimplicialLDLTderived<Eigen::SparseMatrix<double> > &A_factorization = Ae->giveLDLTFactorization(); // compute factorization and get it
    bool negEig                                                          = A_factorization.updateD( minLam, true, numNegEigs );
    OOFEM_LOG_INFO( "numNegEigs = %i\n", numNegEigs );
    Eigen::VectorXd xeig = A_factorization.solve( b_eig ); // Solve the system
    x                    = FloatArray( xeig.begin(), xeig.end() );
    // this->setBifurcation( false );
}

//int EigenSolverStability ::computeEigenValuesVectors( Eigen::SparseMatrix<double> &A, FloatArray &evaluesFA, FloatMatrix &evectorsFM )
//{
//    int neigs = 1, maxTrials = 5, ncvMult = 2, ncv = 200; // Predefined values, should be in input file
//    int nconv = 0, count = 0;
//
//    int Asize = (int)A.rows();
//    Spectra::SparseSymMatProd<double> op( A );
//
//    evectorsFM.resize( Asize, neigs );
//    evaluesFA.resize( neigs );
//
//    while ( nconv == 0 && count < maxTrials ) { // if ncv is too small, increase it ti find solution
//        int ncv2 = std::min( ncv, Asize );
//        Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double> > eigs( op, neigs, ncv2 );
//        eigs.init();
//
//        nconv = eigs.compute( Spectra::SortRule::SmallestAlge ); // Compute eigenvalues
//        if ( eigs.info() == Spectra::CompInfo::Successful ) { // If succesfull, retrieve results
//            // allocate some objects
//            this->setNcv0( ncv );
//            Eigen::VectorXd evector;
//            Eigen::MatrixXd evectors;
//            FloatArray evectorFA;
//            evectorFA.resize( Asize * neigs );
//
//            Eigen::VectorXd evalues = eigs.eigenvalues();
//            evaluesFA               = FloatArray( evalues.begin(), evalues.end() ); // eigenvalues in floatarray
//
//            evectors  = eigs.eigenvectors(); // vectors in matrix
//            evector   = Eigen::Map<Eigen::VectorXd>( evectors.data(), evectors.cols() * evectors.rows() ); // map to vector
//            evectorFA = FloatArray( evector.begin(), evector.end() ); // eigenvectors in floatarray
//
//            // eigenvectors in float matrix
//            for ( size_t j = 0; j < neigs; j++ ) {
//                for ( size_t k = 0; k < evectors.rows(); k++ ) {
//                    evectorsFM.at( k + 1, j + 1 ) = evectorFA.at( k + 1 + evectors.rows() * j ); // Eigenvectors to float matrix
//                }
//            }
//        } else { // if not succesfull
//            count++;
//            ncv = (int)ncvMult * ncv; // Increase ncv
//            OOFEM_LOG_INFO( "eigenvalues not found\n" );
//            OOFEM_LOG_INFO( "ncv =  %i\n", ncv );
//        }
//    }
//    return nconv;
//}

// void EigenSolverStability::solveLDLT( EigenMtrx &Ae, const Eigen::VectorXd &b, Eigen::VectorXd &x)
//{
//     // Get the matrix and LDLT factorization
//     Eigen::SparseMatrix<double> A = Ae.giveMatrix();
//
//     if (this->bifurcation ) { // if bifurcatuion should be performed
//         if ( this->choleskyBif ) { // using cholesky
//             double minLam = 0.;
//             bool negEig = false;
//             int numNegEigs = 0;
//             SimplicialLDLTderived<Eigen::SparseMatrix<double> > &A_factorization = Ae.giveLDLTFactorization(); // compute factorization and get it
//             negEig = A_factorization.updateD( minLam, true, numNegEigs );
//             OOFEM_LOG_INFO( "numNegEigs = %i\n", numNegEigs );
//             x      = A_factorization.solve( b ); // Solve the system
//             //this->setBifurcation( false );
//         }
//         else if ( this->deflationBifurcation ) {
//             Eigen::SparseLU<Eigen::SparseMatrix<double> > A_LU_factorization( A );
//             x = A_LU_factorization.solve( b ); // Compute x using LU factorization
//
//             FloatArray xoldFA  = FloatArray( x.begin(), x.end() );
//             int p = 2;
//             double alph_defl = 1.;
//             double dx_norm  = this->dx_Defl.computeNorm();
//             double gamma = p / ( dx_norm * dx_norm  + alph_defl * pow( dx_norm, p + 2 ) );
//
//             // apply Sherman-Morrison
//             double vt_xold = this->dx_Defl.dotProduct( xoldFA );
//             x = x*( 1 - gamma * vt_xold / ( 1 + gamma * vt_xold ) );
//         }
//         else { // Spectra eigenvector bifurcation
//             FloatArray evectorFAmax, xFA;
//             this->getEigenVectors().copyColumn( evectorFAmax, this->getEigenVectors().giveNumberOfColumns() ); // Chosen eigenvector in array
//
//             // perturbate the solution by the eigenvector
//             xFA = this->alpha * evectorFAmax;
//             x   = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>( xFA.givePointer(), xFA.giveSize() );
//             this->setBifurcation( false );
//         }
//     } else {
//         //SimplicialLDLTderived<Eigen::SparseMatrix<double> > &A_factorization = Ae.giveLDLTFactorization(); // compute factorization and get it
//         Eigen::SparseLU<Eigen::SparseMatrix<double> > A_factorization( A );
//         x = A_factorization.solve( b );
//     }
// }

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



//void EigenSolverStability::setEigenValuesVectors( Eigen::SparseMatrix<double> &A, FloatArray &Xeigs )
//{
//    FloatArray evaluesFA;
//    FloatMatrix evectorsFM;
//    int tmp = this->computeEigenValuesVectors( A, evaluesFA, evectorsFM );
//    this->storeEigenValuesVectors( evaluesFA, evectorsFM, Xeigs );
//}

} // end namespace oofem
