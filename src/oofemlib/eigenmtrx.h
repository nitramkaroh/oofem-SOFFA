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



#ifndef eigenmtrx_h
#define eigenmtrx_h

#include "sparsemtrx.h"
#include "intarray.h"
#include <Eigen/Dense>

#include <Eigen/Sparse>

#define _IFT_EigenMtrx_Name "EigenMtrx"

namespace oofem {
/**
 * Implementation of sparse matrix stored using Eigen library.
 */
template <typename _MatrixType, int _UpLo, typename _Ordering>
// template <typename _MatrixType>
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
        bool answ  = false;
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


enum FactorizationType {
    FT_LLT,
    FT_LU,
    FT_QR,
    FT_LDLT
};

class OOFEM_EXPORT EigenMtrx : public SparseMtrx
{
    

protected:
    Eigen::SparseMatrix<double> EigMat;

    std::vector<bool> areFactorized = { false, false, false, false };
    int versionUpdate;

    // factorizations
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > LLT_factorization; 
    Eigen::SparseLU<Eigen::SparseMatrix<double> > LU_factorization; 
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > QR_factorization; 
    SimplicialLDLTderived<Eigen::SparseMatrix<double> > LDLT_factorization; 

    //bool factorized = false;

    //template <typename Derived>
    //Eigen::SparseSolverBase<Derived> factorization;


public:
    /** Constructor. Before any operation an internal profile must be built.
     * @see buildInternalStructure
     */
    EigenMtrx( int n=0 );

    /// Destructor
    virtual ~EigenMtrx() {}

    // Overloaded methods:
    /*std::unique_ptr<SparseMtrx> clone() const override;
    void times( const FloatArray &x, FloatArray &answer ) const override;
    void timesT( const FloatArray &x, FloatArray &answer ) const override;
    void times( double x ) override;*/
    int buildInternalStructure( EngngModel *, int, const UnknownNumberingScheme &s ) override;
    int assemble( const IntArray &loc, const FloatMatrix &mat ) override;
    int assemble( const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat ) override;
    bool canBeFactorized() const override { return false; }
    void zero() override;
    double &at( int i, int j ) override;
    double at( int i, int j ) const override;
    /*void toFloatMatrix( FloatMatrix &answer ) const override;*/
    void printYourself() const override;
    const char *giveClassName() const override { return "EigenMtrx"; }
    SparseMtrxType giveType() const override { return SMT_EigenMtrx; }
    bool isAsymmetric() const override { return true; }

    Eigen::SparseMatrix<double> giveMatrix();

    //template <typename Derived>
    //std::shared_ptr<Eigen::SparseSolverBase<Derived> > giveFactorization( FactorizationType factorizationType );

    template <typename Derived>
    Eigen::SparseSolverBase<Derived>& giveFactorization( FactorizationType Factorization );


    //template <typename Derived>
    //void setFactorization( FactorizationType factorizationType, Eigen::SparseSolverBase<Derived> factorization );

    //void setLDLTFactorization(Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& factorization );

    //std::shared_ptr<SimplicialLDLTderived<Eigen::SparseMatrix<double> > > giveLDLTFactorization();
    SimplicialLDLTderived<Eigen::SparseMatrix<double> >& EigenMtrx::giveLDLTFactorization();

    void computeFactorization( FactorizationType factorizationType );
    
    bool isFactorized( FactorizationType factorizationType );

    void times( const FloatArray &x, FloatArray &answer ) const override;
};
} // end namespace oofem
#endif // eigenmtrx_h
