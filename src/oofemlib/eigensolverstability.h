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

#ifndef eigensolverstability_h
#define eigensolverstability_h

#include "eigensolver.h"
#include "eigenmtrx.h"
#include "bifurcationinterface.h"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include "Eigen/src/Core/util/Constants.h"
#include "Eigen/src/OrderingMethods/Ordering.h"


#define _IFT_EigenSolverStability_Name "eigenstability"

namespace oofem {
/**
 * Implements the solution of linear system of equation in the form @f$ A\cdot x=b @f$ using solvers
 * from eigen.tuxfamily.org.
 */
class OOFEM_EXPORT EigenSolverStability : public EigenSolver, public BifurcationInterface
{
protected:
    int ncv0 = 5; // multiplier for eigenvector analysis
public:
    /**
     * Constructor.
     * @param d Domain which solver belongs to.
     * @param m Engineering model which solver belongs to.
     */
    EigenSolverStability( Domain *d, EngngModel *m );

    virtual ~EigenSolverStability();

    /// Initializes receiver from given record.
    void initializeFrom( InputRecord &ir ) override;
    ConvergedReason solve( SparseMtrx &A, FloatArray &b, FloatArray &x ) override;

    //void solveLDLT( EigenMtrx &A, const Eigen::VectorXd &b, Eigen::VectorXd &x);
    bool checkPD( SparseMtrx &A ) override;
    void CholeskyUpdate( SparseMtrx &Ae, FloatArray &b, FloatArray &x ) override;
    void setNcv0( double ncv0new ) { this->ncv0 = ncv0new; }
    const char *giveClassName() const override { return "EigenSolverStability"; }
    LinSystSolverType giveLinSystSolverType() const override { return ST_EigenStability; }
    SparseMtrxType giveRecommendedMatrix( bool symmetric ) const override { return SMT_EigenMtrx; }

    //int computeEigenValuesVectors( Eigen::SparseMatrix<double> &A, FloatArray &evaluesFA, FloatMatrix &evectorsFM );
    int computeEigenValuesVectors( SparseMtrx &Ag, FloatArray &evaluesFA, FloatMatrix &evectorsFM ) override;
    // ConvergedReason solveBifurcation( SparseMtrx &A, FloatArray &b, FloatArray &x );
    //void setEigenValuesVectors( Eigen::SparseMatrix<double> &A, FloatArray &Xeigs );
    bool canCholeskyBifurcation() override { return true; };
    bool canComputeEigenVectors() override { return true; };
};
} // end namespace oofem
#endif // eigensolver_hj