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

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
#include "Eigen/src/Core/util/Constants.h"
#include "Eigen/src/OrderingMethods/Ordering.h"


#define _IFT_EigenSolverStability_Name "eigenstability"

namespace oofem {

template <typename _MatrixType, int _UpLo, typename _Ordering>
//template <typename _MatrixType>
class SimplicialLDLTderived : public Eigen::SimplicialLDLT<_MatrixType, _UpLo, _Ordering>
{
public:
    typedef _MatrixType MatrixType;
    typedef Eigen::SimplicialLDLT<_MatrixType, _UpLo, _Ordering> Base;


public:
    /** Default constructor */
    SimplicialLDLTderived() :
        Base() {}

    /** Constructs and performs the LLT factorization of \a matrix */
    explicit SimplicialLDLTderived( const MatrixType &matrix ) :
        Base( matrix ) {}

    bool updateD( double &minEig, const bool update, int &numNegEigs )
    {
        numNegEigs = 0;
        bool answ = false;
        eigen_assert( this->m_factorizationIsOk && "Simplicial LDLT not factorized" );
        minEig = m_diag.coeffRef( 0 );
        for ( auto &Di : m_diag ) {
            if ( Di < minEig ) minEig = Di;
            if ( Di < 0 ) {
                if ( update ) Di *= ( -1. );
                answ = true;
                numNegEigs++;
            };
        };
        return answ;
    }
};

template <typename _MatrixType, int _UpLo = Eigen ::Lower, typename _Ordering = Eigen ::AMDOrdering<typename _MatrixType::StorageIndex> >
class SimplicialLDLTderived;


/**
 * Implements the solution of linear system of equation in the form @f$ A\cdot x=b @f$ using solvers
 * from eigen.tuxfamily.org.
 */
class OOFEM_EXPORT EigenSolverStability : public EigenSolver
{
protected:


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

    virtual void solveLDLT( Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b, Eigen::VectorXd &x, bool doBifurcation ) override;
    bool checkPD( SparseMtrx &A );

    const char *giveClassName() const override { return "EigenSolverStability"; }
    LinSystSolverType giveLinSystSolverType() const override { return ST_EigenStability; }
    SparseMtrxType giveRecommendedMatrix( bool symmetric ) const override { return SMT_EigenMtrx; }
};
} // end namespace oofem
#endif // eigensolver_h