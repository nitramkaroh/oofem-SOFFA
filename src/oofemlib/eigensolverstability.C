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
    EigenSolver( d, m )
{
}

EigenSolverStability ::~EigenSolverStability() {}

void EigenSolverStability ::initializeFrom( InputRecord &ir )
{
    EigenSolver ::initializeFrom( ir );
}

ConvergedReason EigenSolverStability ::solve( SparseMtrx &A, FloatArray &b, FloatArray &x )
{
    int neqs = b.giveSize(); // Number of equations

    EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &A );

    Eigen::SparseMatrix<double> A_eig = Ae->giveMatrix();

    // Construct right hand side vetor
    Eigen::VectorXd b_eig = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>( b.givePointer(), neqs );

    Timer timer;
    timer.startTimer();

    Eigen::VectorXd x_eig; // Allocate vector of RHS

    // Create factorization
    SimplicialLDLTderived<Eigen::SparseMatrix<double> >& A_factorization_cast = Ae->giveLDLTFactorization(); // compute factorization and get it
   
    //// Why doesnt work?
    //Eigen::SparseSolverBase<SimplicialLDLTderived<Eigen::SparseMatrix<double> > >* A_factorization = Ae->giveFactorization<SimplicialLDLTderived<Eigen::SparseMatrix<double> > >( FT_LDLT );
    //SimplicialLDLTderived<Eigen::SparseMatrix<double> > *A_factorization_cast   = dynamic_cast<SimplicialLDLTderived<Eigen::SparseMatrix<double> > *>( A_factorization );

    this->solveLDLT( A_eig, b_eig, x_eig, A_factorization_cast );


    // Copy/move values to FloatArray x
    x = FloatArray( x_eig.begin(), x_eig.end() );

    timer.stopTimer();
    OOFEM_LOG_INFO( "EigenSolver:  User time consumed by solution: %.2fs\n", timer.getUtime() );

    return CR_CONVERGED;
}

void EigenSolverStability::solveLDLT( Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b, Eigen::VectorXd &x, SimplicialLDLTderived<Eigen::SparseMatrix<double> > &A_factorization )
{
    double minLam = 0.;
    bool negEig   = false;
    int numNegEigs = 0;
    negEig = A_factorization.updateD( minLam, false, numNegEigs );

    x = A_factorization.solve( b ); // Solve the system

    //OOFEM_LOG_INFO( "lam_min = %.4f\n", minLam );

    if ( negEig ) {
        isSingular = true;
        OOFEM_LOG_INFO( "Negative eigenvalue\n" );
    }

    if ( minLam < 0. && this->bifurcation ) { // if negative eigen value, compute eigenvectors and perturbe the solution
        if ( this->choleskyBif ) { // using cholesky
            negEig = A_factorization.updateD( minLam, true, numNegEigs );
            x      = A_factorization.solve( b ); // Solve the system
            //this->setBifurcation( false );
        }
        else { // Spectra eigenvector bifurcation
            Spectra::SparseSymMatProd<double> op( A );
            // int neigs = numNegEigs;
            int neigs      = 1;
            int ncv        = this->ncv0;
            int nconv      = 0; // number of converged eigenvalues
            int maxTrials  = 5;
            int count      = 0;
            double ncvMult = 2;
            while ( nconv == 0 && count < maxTrials ) { // if ncv is too small, increase it ti find solution
                Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double> > eigs( op, neigs, ncv );
                eigs.init();
                nconv = eigs.compute( Spectra::SortRule::SmallestAlge ); // Compute eigenvalues
                if ( eigs.info() == Spectra::CompInfo::Successful ) { // If succesfull, retrieve results
                    this->setNcv0( ncv );
                    // allocate some objects
                    Eigen::VectorXd evalues, evector, xmod;
                    Eigen::MatrixXd evectors;
                    FloatArray evaluesFA, evectorFA, xFA, evectorFAmax;
                    FloatMatrix evectorsFM;
                    double xdotx2, vdotx2;

                    evalues   = eigs.eigenvalues();
                    evaluesFA = FloatArray( evalues.begin(), evalues.end() ); // eigenvalues in floatarray

                    int indMin = evaluesFA.giveSize();

                    OOFEM_LOG_INFO( "lam = " ); // print eigenvalues
                    evaluesFA.printYourself();
                    OOFEM_LOG_INFO( "\n" );

                    evectors  = eigs.eigenvectors(); // vectors in matrix
                    evector   = Eigen::Map<Eigen::VectorXd>( evectors.data(), evectors.cols() * evectors.rows() ); // map to vector
                    evectorFA = FloatArray( evector.begin(), evector.end() ); // eigenvectors in floatarray

                    // eigenvectors in float matrix
                    evectorsFM.resize( evectors.rows(), neigs );
                    for ( size_t j = 0; j < neigs; j++ ) {
                        for ( size_t k = 0; k < evectors.rows(); k++ ) {
                            evectorsFM.at( k + 1, j + 1 ) = evectorFA.at( k + 1 + evectors.rows() * j ); // Eigenvectors to float matrix
                        }
                    }
                    evectorsFM.copyColumn( evectorFAmax, indMin ); // Chosen 1 eigenvector in array


                    // perturbate the solution by the eigenvector
                    xFA = this->alpha * evectorFAmax;
                    x   = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>( xFA.givePointer(), xFA.giveSize() );

                    //// incorporate the instable solution
                    // xFA    = FloatArray( x.begin(), x.end() ); // the instable solution
                    // xdotx2 = xFA.dotProduct( xFA );
                    // vdotx2 = evectorFAmax.dotProduct( xFA );
                    //  if ( abs( evaluesFA.at( indMin ) ) > 1e-6 ) {
                    //  double alpha = 1 / sqrt( xdotx2 - vdotx2 * vdotx2 );
                    //  xFA          = alpha * ( vdotx2 * xFA - xdotx2 * evectorFAmax );
                    //  double alpha = 3e1;
                    // xFA = this->alpha * evectorFAmax;
                    // x   = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>( xFA.givePointer(), xFA.giveSize() );
                    // }

                    this->setBifurcation( false );

                } else { // if not succesfull
                    count++;
                    ncv = (int)ncvMult * ncv; // Increase ncv
                    OOFEM_LOG_INFO( "eigenvalues not found\n" );
                    OOFEM_LOG_INFO( "ncv =  %i\n", ncv );
                    this->setBifurcation( false );
                }
            }
        }        
    } 
}

bool EigenSolverStability::checkPD( SparseMtrx &A )
{
    EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &A );
    Eigen::SparseMatrix<double> A_eig = Ae->giveMatrix();

    SimplicialLDLTderived<Eigen::SparseMatrix<double> >& A_factorization = Ae->giveLDLTFactorization();

    double minLam  = 0.;
    bool negEig    = false;
    int numNegEigs = 0;
    negEig         = A_factorization.updateD( minLam, false, numNegEigs );

    return !negEig;
}

ConvergedReason EigenSolverStability ::solveBifurcation( SparseMtrx &A, FloatArray &b, FloatArray &x )
{
    int neqs = b.giveSize(); // Number of equations

    EigenMtrx *Ae = dynamic_cast<EigenMtrx *>( &A );

    Eigen::SparseMatrix<double> A_eig = Ae->giveMatrix();

    // Construct right hand side vetor
    Eigen::VectorXd b_eig = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>( b.givePointer(), neqs );

    Timer timer;
    timer.startTimer();

    Eigen::VectorXd x_eig; // Allocate vector of RHS

    SimplicialLDLTderived<Eigen::SparseMatrix<double> >& A_factorization_cast = Ae->giveLDLTFactorization(); // compute factorization and get it

    // Create factorization
    if ( method.compare( "ldlt" ) == 0 ) {
        this->setBifurcation( true );
        this->solveLDLT( A_eig, b_eig, x_eig, A_factorization_cast );

    } else {
        OOFEM_LOG_INFO( "Wrong factorization used, bifucration cannot be performed\n" );
    }

    // Copy/move values to FloatArray x
    x = FloatArray( x_eig.begin(), x_eig.end() );


    timer.stopTimer();
    OOFEM_LOG_INFO( "EigenSolver:  User time consumed by solution: %.2fs\n", timer.getUtime() );

    return CR_CONVERGED;
}

} // end namespace oofem
