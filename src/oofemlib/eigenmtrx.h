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

#include <fstream>
#include <string>

namespace oofem {
/**
 * Eigen's SimplicialLDLT (P A P^T = L D L^T with a fill-reducing permutation P) with access to the
 * diagonal D, whose signs give the inertia of the symmetric matrix A.
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

    /** Constructs and performs the LDL^T factorization of \a matrix */
    explicit SimplicialLDLTderived( const MatrixType &matrix ) :
        Base( matrix ) {}

    /**
     * Inspects the diagonal D of the factorization. By Sylvester's law of inertia the number of its
     * negative entries is the number of negative eigenvalues of the matrix.
     * @param minPivot Smallest entry of D (a pivot: only its sign means the same as for an eigenvalue).
     * @param update If true, the negative entries of D are replaced by their absolute values (used by
     * the Cholesky bifurcation update of EigenSolverStability).
     * @param numNegEigs Number of negative entries of D.
     * @return True if D has a negative entry, i.e. the matrix is not positive semidefinite.
     */
    bool updateD( double &minPivot, const bool update, int &numNegEigs )
    {
        numNegEigs = 0;
        bool hasNegative = false;
        eigen_assert( this->m_factorizationIsOk && "Simplicial LDLT not factorized" );
        minPivot = this->m_diag.coeffRef( 0 );
        for ( auto &Di : this->m_diag ) {
            if ( Di < minPivot ) minPivot = Di;
            if ( Di < 0 ) {
                if ( update ) Di *= ( -1. );
                hasNegative = true;
                numNegEigs++;
            };
        };
        return hasNegative;
    }
};

template <typename _MatrixType, int _UpLo = Eigen ::Lower, typename _Ordering = Eigen ::AMDOrdering<typename _MatrixType::StorageIndex> >
class SimplicialLDLTderived;


/// Factorizations of EigenMtrx (the index into its cache of factorizations).
enum FactorizationType {
    FT_LLT,
    FT_LU,
    FT_QR,
    FT_LDLT
};

/**
 * Sparse matrix stored with the Eigen library (Eigen::SparseMatrix<double>, type SMT_EigenMtrx).
 * assemble only collects (row, column, value) triplets, skipping zero entries off the diagonal; the
 * compressed matrix is built from them on demand (createMatrixFromTriplets, giveMatrix) whenever the
 * version of the matrix has changed. The factorizations (LLT, LU, QR, LDL^T) are cached and recomputed
 * after a change of the version. getBlock, the static condensations and removeRowsAndColumns take
 * 1-based equation numbers.
 */
class OOFEM_EXPORT EigenMtrx : public SparseMtrx
{


protected:
    Eigen::SparseMatrix<double> EigMat;

    /// Whether each FactorizationType is computed for the version versionUpdate.
    std::vector<bool> areFactorized = { false, false, false, false };
    int versionUpdate;

    std::vector<Eigen::Triplet<double> > triplets; // Allocate vector of triplets
    /// Version of the matrix from which EigMat was last built.
    int BuiltFromTripletsAtVersion = 0;


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
    EigenMtrx( int n=0);

    EigenMtrx( Eigen::SparseMatrix<double>& EigMat );

    EigenMtrx( const EigenMtrx &EigenMtrxInput );

    //EigenMtrx &operator= ( Eigen::SparseMatrix<double> & EigMat );

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
    //void printYourself() const override;
    const char *giveClassName() const override { return "EigenMtrx"; }
    SparseMtrxType giveType() const override { return SMT_EigenMtrx; }
    bool isAsymmetric() const override { return true; }

    /// Returns the Eigen matrix, built from the triplets first if the matrix has changed.
    Eigen::SparseMatrix<double>& giveMatrix();

    //template <typename Derived>
    //std::shared_ptr<Eigen::SparseSolverBase<Derived> > giveFactorization( FactorizationType factorizationType );

    /// Returns the requested factorization, computed first if it does not belong to the current version.
    template <typename Derived>
    Eigen::SparseSolverBase<Derived>& giveFactorization( FactorizationType Factorization );


    //template <typename Derived>
    //void setFactorization( FactorizationType factorizationType, Eigen::SparseSolverBase<Derived> factorization );

    //void setLDLTFactorization(Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& factorization );

    //std::shared_ptr<SimplicialLDLTderived<Eigen::SparseMatrix<double> > > giveLDLTFactorization();
    //    SimplicialLDLTderived<Eigen::SparseMatrix<double> >& EigenMtrx::giveLDLTFactorization();
    /// Returns the LDL^T factorization (with access to D), computed first if needed.
    SimplicialLDLTderived<Eigen::SparseMatrix<double> >& giveLDLTFactorization();

    /// Computes the given factorization of the current matrix and marks it valid for the current version.
    void computeFactorization( FactorizationType factorizationType );

    /// Returns true if the given factorization belongs to the current version of the matrix.
    bool isFactorized( FactorizationType factorizationType );

    void times( const FloatArray &x, FloatArray &answer ) const override;

    /// Builds the compressed Eigen matrix from the triplets if the matrix has changed since the last build.
    void createMatrixFromTriplets();

    void printYourself() const override;

    /// Writes the matrix in Eigen's text format to a file.
    void printToFile( const std::string &filename );

    ///////////
    //void getBlock( int i, int j, int p, int q, EigenMtrx &BlockMatrix );
    /// Returns the block of p rows and q columns starting at row i, column j (taken from the stored Eigen matrix as it is).
    std::unique_ptr<EigenMtrx> getBlock( int i, int j, int p, int q );
    /**
     * Static condensation of the equations i, ..., n (block P) out of the matrix:
     * returns the Schur complement K_UU - K_UP K_PP^-1 K_PU of the equations 1, ..., i-1
     * (K_PP factorized by LDL^T, so it should be symmetric and nonsingular).
     */
    std::unique_ptr<EigenMtrx> doStaticCondensation( int i );
    /// The same Schur complement as doStaticCondensation (K_PP factorized by LU), stored in this matrix.
    void doStaticCondensationOwn( int i);
    /// Returns the block K_UU of the equations 1, ..., i-1 only, without condensation.
    std::unique_ptr<EigenMtrx> doStaticCondensationLagrange( int i );

    /// Returns the matrix without the rows and columns of the given equation numbers.
    std::unique_ptr<EigenMtrx> removeRowsAndColumns( const IntArray &indsToDelete );
    ///////////
};
} // end namespace oofem
#endif // eigenmtrx_h
