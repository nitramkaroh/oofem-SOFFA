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

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace oofem {
/**
 * Eigen's SimplicialLDLT (P A P^T = L D L^T with a fill-reducing permutation P) with read access to the
 * diagonal D, whose signs give the inertia of the symmetric matrix A, and with a solve that uses the
 * absolute values of D. Neither of them modifies the factorization, so that it keeps describing A.
 */
template <typename _MatrixType, int _UpLo = Eigen::Lower, typename _Ordering = Eigen::AMDOrdering<typename _MatrixType::StorageIndex> >
class SimplicialLDLTderived : public Eigen::SimplicialLDLT<_MatrixType, _UpLo, _Ordering>
{
public:
    typedef _MatrixType MatrixType;
    typedef Eigen::SimplicialLDLT<_MatrixType, _UpLo, _Ordering> Base;
    typedef typename Base::VectorType VectorType;


public:
    /** Default constructor */
    SimplicialLDLTderived() :
        Base() {}

    /** Constructs and performs the LDL^T factorization of \a matrix */
    explicit SimplicialLDLTderived( const MatrixType &matrix ) :
        Base( matrix ) {}

    /**
     * Reads the inertia of the factorized matrix from the diagonal D: by Sylvester's law of inertia the
     * number of negative entries of D is the number of negative eigenvalues of the matrix.
     * Eigen stops the factorization at the first zero pivot (info() != Success) and leaves the following
     * entries of D undefined, so only the pivots up to the breakdown are inspected then. Their negative
     * count is a lower bound of the number of negative eigenvalues, and the matrix is certainly not
     * positive definite (a positive definite matrix has no zero pivot).
     * @param minPivot Smallest inspected finite entry of D (a pivot: only its sign means the same as for
     * an eigenvalue); NaN if there is none.
     * @param numNegative Number of negative inspected entries of D.
     * @return True if the factorization succeeded and every entry of D is finite, i.e. numNegative is the
     * number of negative eigenvalues; false otherwise.
     */
    bool giveInertia( double &minPivot, int &numNegative ) const
    {
        minPivot    = std::numeric_limits<double>::quiet_NaN();
        numNegative = 0;
        if ( !this->m_factorizationIsOk ) {
            return false; // no factorization has been attempted
        }

        const VectorType D     = this->vectorD();
        const Eigen::Index n   = D.size();
        bool complete          = this->info() == Eigen::Success;
        Eigen::Index inspected = n;
        if ( !complete ) {
            // the pivots before the breakdown are nonzero, the breakdown pivot itself is stored as zero
            inspected = 0;
            while ( inspected < n && D[inspected] != 0. ) {
                inspected++;
            }
            if ( inspected < n ) {
                inspected++;
            }
        }

        double smallest = std::numeric_limits<double>::infinity();
        for ( Eigen::Index i = 0; i < inspected; i++ ) {
            const double Di = D[i];
            if ( !std::isfinite( Di ) ) {
                complete = false;
                continue;
            }
            if ( Di < smallest ) {
                smallest = Di;
            }
            if ( Di < 0. ) {
                numNegative++;
            }
        }
        if ( std::isfinite( smallest ) ) {
            minPivot = smallest;
        }
        return complete;
    }

    /**
     * Solves P^T L |D| L^T P x = b, i.e. with the absolute values of the pivots (the matrix made positive
     * definite by the Cholesky bifurcation update of EigenSolverStability), without modifying the
     * factorization. The steps are those of Eigen's own solve. Valid only after a successful factorization.
     */
    VectorType solveWithAbsolutePivots( const VectorType &b ) const
    {
        eigen_assert( this->m_factorizationIsOk && this->info() == Eigen::Success && "Simplicial LDLT not factorized" );
        VectorType x;
        if ( this->permutationP().size() > 0 ) {
            x = this->permutationP() * b;
        } else {
            x = b;
        }
        if ( this->m_matrix.nonZeros() > 0 ) { // otherwise L == I
            this->matrixL().solveInPlace( x );
        }
        const VectorType absD = this->vectorD().cwiseAbs();
        if ( absD.size() > 0 ) {
            x = absD.asDiagonal().inverse() * x;
        }
        if ( this->m_matrix.nonZeros() > 0 ) { // otherwise U == I
            this->matrixU().solveInPlace( x );
        }
        if ( this->permutationP().size() > 0 ) {
            x = this->permutationPinv() * x;
        }
        return x;
    }
};


/// Factorizations of EigenMtrx (the index into its cache of factorizations).
enum FactorizationType {
    FT_LLT,
    FT_LU,
    FT_QR,
    FT_LDLT
};

/**
 * Sparse matrix stored with the Eigen library (Eigen::SparseMatrix<double>, type SMT_EigenMtrx).
 * The compressed matrix EigMat holds the values. assemble only collects (row, column, value) triplets,
 * skipping zero entries off the diagonal; they are added to EigMat on demand (createMatrixFromTriplets,
 * called by every access to the values), so that the matrix is built once per assembly. Values changed
 * through at() stay in EigMat. The version counts the changes of the values; the factorizations
 * (LLT, LU, QR, LDL^T) are cached together with the version they were computed for. at() takes 1-based
 * indices like every SparseMtrx; getBlock, the static condensations and removeRowsAndColumns take
 * 1-based equation numbers as well.
 */
class OOFEM_EXPORT EigenMtrx : public SparseMtrx
{


protected:
    /// The matrix, kept compressed. Contributions still waiting in triplets are not included yet.
    mutable Eigen::SparseMatrix<double> EigMat;

    /// Contributions assembled since EigMat was last built, with 0-based indices.
    mutable std::vector<Eigen::Triplet<double> > triplets;
    /// True if triplets holds contributions that have not been added to EigMat yet.
    mutable bool tripletsPending = false;

    /// Version of the matrix each FactorizationType was last computed for (-1: never computed).
    std::array<SparseMtrxVersionType, 4> factorizedAtVersion = { { -1, -1, -1, -1 } };


    // factorizations
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > LLT_factorization;
    Eigen::SparseLU<Eigen::SparseMatrix<double> > LU_factorization;
    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > QR_factorization;
    SimplicialLDLTderived<Eigen::SparseMatrix<double> > LDLT_factorization;


public:
    /** Constructor. Before any operation an internal profile must be built.
     * @see buildInternalStructure
     */
    EigenMtrx( int n=0);

    /// Wraps a ready Eigen matrix: it is the value of the matrix, there are no pending contributions.
    EigenMtrx( Eigen::SparseMatrix<double> matrix );

    /// Copies the current values (pending contributions included); the factorizations are not copied.
    EigenMtrx( const EigenMtrx &EigenMtrxInput );

    /// Destructor
    virtual ~EigenMtrx() {}

    // Overloaded methods:
    int buildInternalStructure( EngngModel *, int, const UnknownNumberingScheme &s ) override;
    int assemble( const IntArray &loc, const FloatMatrix &mat ) override;
    int assemble( const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat ) override;
    bool canBeFactorized() const override { return false; }
    void zero() override;
    /**
     * Returns a reference to the entry (i, j), 1-based. An entry that is not stored is inserted first.
     * Counts as a change of the values (the version is increased), since the caller may write to it.
     */
    double &at( int i, int j ) override;
    /// Returns the entry (i, j), 1-based; zero if it is not stored.
    double at( int i, int j ) const override;
    const char *giveClassName() const override { return "EigenMtrx"; }
    SparseMtrxType giveType() const override { return SMT_EigenMtrx; }
    bool isAsymmetric() const override { return true; }

    /**
     * Returns the Eigen matrix (read only) with all assembled contributions added. Its values are changed
     * only through assemble, zero and at(), which keep the version and the cached factorizations consistent.
     */
    const Eigen::SparseMatrix<double> &giveMatrix() const;

    /// Returns the requested factorization, computed first if it does not belong to the current version.
    template <typename Derived>
    Eigen::SparseSolverBase<Derived>& giveFactorization( FactorizationType Factorization );

    /// Returns the LDL^T factorization (with access to D), computed first if needed; its info() tells whether it succeeded.
    SimplicialLDLTderived<Eigen::SparseMatrix<double> >& giveLDLTFactorization();

    /// Computes the given factorization of the current matrix and records the version it belongs to (also when it fails: check its info()).
    void computeFactorization( FactorizationType factorizationType );

    /// Returns true if the given factorization has been computed for the current version of the matrix.
    bool isFactorized( FactorizationType factorizationType ) const;

    void times( const FloatArray &x, FloatArray &answer ) const override;

    /// Adds the pending contributions (triplets) to the matrix. The version is not changed: the values were counted when assembled.
    void createMatrixFromTriplets() const;

    void printYourself() const override;

    /// Writes the matrix in Eigen's text format to a file.
    void printToFile( const std::string &filename ) const;

    ///////////
    /// Returns the block of p rows and q columns starting at row i, column j (1-based).
    std::unique_ptr<EigenMtrx> getBlock( int i, int j, int p, int q ) const;
    /**
     * Static condensation of the equations i, ..., n (block P) out of the matrix:
     * returns the Schur complement K_UU - K_UP K_PP^-1 K_PU of the equations 1, ..., i-1.
     * K_PP is factorized by LDL^T, or by LU if LDL^T breaks down on a zero pivot; a singular K_PP is an error.
     */
    std::unique_ptr<EigenMtrx> doStaticCondensation( int i ) const;
    /// The same Schur complement as doStaticCondensation, stored in this matrix (its size becomes i-1).
    void doStaticCondensationOwn( int i);
    /// Returns the block K_UU of the equations 1, ..., i-1 only, without condensation.
    std::unique_ptr<EigenMtrx> doStaticCondensationLagrange( int i ) const;

    /// Returns the matrix without the rows and columns of the given equation numbers.
    std::unique_ptr<EigenMtrx> removeRowsAndColumns( const IntArray &indsToDelete ) const;
    ///////////

protected:
    /// Stops with an error if (i, j) is not a 1-based position inside the matrix.
    void checkBounds( int i, int j ) const;
    /// Returns true if the entry (r, c), 0-based, is stored in EigMat (which must be compressed).
    bool hasStoredEntry( int r, int c ) const;
    /// Schur complement of the equations i, ..., n (1-based), shared by the static condensations.
    Eigen::SparseMatrix<double> computeSchurComplement( int i ) const;
};
} // end namespace oofem
#endif // eigenmtrx_h
