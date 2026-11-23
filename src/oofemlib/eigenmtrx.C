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


#include "eigenmtrx.h"
#include "floatarray.h"
#include "engngm.h"
#include "domain.h"
#include "element.h"
#include "sparsemtrxtype.h"
#include "activebc.h"
#include "classfactory.h"
#include <iostream>
#include <iomanip>

#include <algorithm>
#include <set>
#include <utility>


#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>
//#include <iostream>

namespace oofem {
REGISTER_SparseMtrx( EigenMtrx, SMT_EigenMtrx );


EigenMtrx ::EigenMtrx( int n ) :
    SparseMtrx( n, n )
{
    EigMat = Eigen::SparseMatrix<double>( n, n );
}

EigenMtrx ::EigenMtrx( Eigen::SparseMatrix<double> matrix ) :
    SparseMtrx( (int)matrix.rows(), (int)matrix.cols() ),
    EigMat( std::move( matrix ) )
{
    // A ready matrix is the value itself: there are no pending contributions to build it from,
    // so it must never be rebuilt from the (empty) triplets.
    EigMat.makeCompressed();
    this->version = 1;
}

EigenMtrx ::EigenMtrx( const EigenMtrx &EigenMtrxInput ) :
    SparseMtrx( EigenMtrxInput.giveNumberOfRows(), EigenMtrxInput.giveNumberOfColumns() )
{
    // The current values are copied: giveMatrix adds the pending contributions of the source first,
    // so values changed after the assembly (through at()) are kept, not rebuilt from old triplets.
    this->EigMat  = EigenMtrxInput.giveMatrix();
    this->version = EigenMtrxInput.version;
}

//EigenMtrx &EigenMtrx ::operator=( Eigen::SparseMatrix<double> &EigMat )
//{
//}

int EigenMtrx ::buildInternalStructure( EngngModel *eModel, int di, const UnknownNumberingScheme &s )
{
    IntArray loc;
    Domain *domain = eModel->giveDomain( di );
    int neq        = eModel->giveNumberOfDomainEquations( di, s );

    EigMat.resize( neq, neq ); // Resize the matrix (all values are dropped)
    triplets.clear(); // contributions to the previous structure would not fit the new one
    triplets.reserve( 500*neq );
    tripletsPending = false;

    nRows = nColumns = neq;


    //// allocation map
    //std ::vector<std ::set<int> > columns( neq );

    //// Eigen::Triplet<double> T;


    //for ( auto &elem : domain->giveElements() ) {
    //    elem->giveLocationArray( loc, s );

    //    for ( int ii : loc ) {
    //        if ( ii > 0 ) {
    //            for ( int jj : loc ) {
    //                if ( jj > 0 ) {
    //                    columns[jj - 1].insert( ii - 1 ); // for each column nonzero rows are stored

    //                    //Eigen::Triplet<double> T( ii - 1, jj - 1, 0. ); // Create triplet storing (row, col, val)
    //                    //this->triplets.push_back( T ); // Add to vector of triplets
    //                }
    //            }
    //        }
    //    }
    //}


    //// loop over active boundary conditions
    //std ::vector<IntArray> r_locs;
    //std ::vector<IntArray> c_locs;

    //for ( auto &gbc : domain->giveBcs() ) {
    //    ActiveBoundaryCondition *bc = dynamic_cast<ActiveBoundaryCondition *>( gbc.get() );
    //    if ( bc != NULL ) {
    //        bc->giveLocationArrays( r_locs, c_locs, UnknownCharType, s, s );
    //        for ( std ::size_t k = 0; k < r_locs.size(); k++ ) {
    //            IntArray &krloc = r_locs[k];
    //            IntArray &kcloc = c_locs[k];
    //            for ( int ii : krloc ) {
    //                if ( ii ) {
    //                    for ( int jj : kcloc ) {
    //                        if ( jj ) {
    //                            columns[jj - 1].insert( ii - 1 );

    //                            //Eigen::Triplet<double> T( ii - 0, jj - 0, 0. ); // Create triplet storing (row, col, val)
    //                            //this->triplets.push_back( T ); // Add to vector of triplets
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }
    //}

    //// std::vector<int> ColReserve( neq ); // Allocate vector of reserved number of nonzeros for each column
    //Eigen::VectorXi ColReserve( neq ); // Allocate vector of reserved number of nonzeros for each column

    //for ( int i = 0; i < neq; i++ ) {
    //    ColReserve[i] = columns[i].size() * 2; // overestimate two times
    //}


     //EigMat.reserve( ColReserve ); // Reserve the number of nonzeros
    //this->EigMat.setFromTriplets( this->triplets.begin(), this->triplets.end() );
    // this->EigMat.makeCompressed();

    this->version++;
    return true;
}


int EigenMtrx ::assemble( const IntArray &loc, const FloatMatrix &mat )
{
    int dim = mat.giveNumberOfRows();
    int ii, jj;
    for ( int j = 0; j < dim; j++ ) {
        jj = loc[j];
        if ( jj ) {
            for ( int i = 0; i < dim; i++ ) {
                ii = loc[i];
                if ( ii ) {
                    //EigMat.coeffRef( ii - 1, jj - 1 ) += mat( i, j );

                    // Structural zeros are not stored. An IGA patch element matrix is
                    // mostly zero, since each knot span couples only its own basis
                    // functions, so storing them inflates the triplet list by an order
                    // of magnitude. Diagonal entries are kept unconditionally, so that
                    // at() finds them without inserting into the compressed matrix.
                    if ( mat( i, j ) != 0. || ii == jj ) {
                        // when set from triplets
                        this->triplets.push_back( Eigen::Triplet<double>( ii - 1, jj - 1, mat( i, j ) ) ); // Add to vector of triplets
                    }
                }
            }
        }
    }

    this->tripletsPending = true;
    this->version++;
    return 1;
}

int EigenMtrx ::assemble( const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat )
{
    for ( int j = 0; j < mat.giveNumberOfColumns(); j++ ) {
        int jj = cloc[j];
        if ( jj ) {
            for ( int i = 0; i < mat.giveNumberOfRows(); i++ ) {
                int ii = rloc[i];
                if ( ii ) {
                    //EigMat.coeffRef( ii - 1, jj - 1 ) += mat( i, j );

                    // Structural zeros are not stored. An IGA patch element matrix is
                    // mostly zero, since each knot span couples only its own basis
                    // functions, so storing them inflates the triplet list by an order
                    // of magnitude. Diagonal entries are kept unconditionally, so that
                    // at() finds them without inserting into the compressed matrix.
                    if ( mat( i, j ) != 0. || ii == jj ) {
                        // when set from triplets
                        this->triplets.push_back( Eigen::Triplet<double>( ii - 1, jj - 1, mat( i, j ) ) ); // Add to vector of triplets
                    }
                }
            }
        }
    }

    this->tripletsPending = true;
    this->version++;
    return 1;
}

void EigenMtrx ::zero()
{
    EigMat.setZero();
    this->triplets.clear(); // When set from triplets
    this->tripletsPending = false;
    this->version++;
}

void EigenMtrx ::checkBounds( int i, int j ) const
{
    if ( i < 1 || i > this->giveNumberOfRows() || j < 1 || j > this->giveNumberOfColumns() ) {
        OOFEM_ERROR( "Array accessing exception -- (%d,%d) out of bounds (%d x %d)", i, j, this->giveNumberOfRows(), this->giveNumberOfColumns() );
    }
}

bool EigenMtrx ::hasStoredEntry( int r, int c ) const
{
    // in a compressed matrix the row indices of every column are sorted
    const int *outer = EigMat.outerIndexPtr();
    const int *inner = EigMat.innerIndexPtr();
    return std::binary_search( inner + outer[c], inner + outer[c + 1], r );
}

double &EigenMtrx ::at( int i, int j )
{
    this->checkBounds( i, j );
    this->createMatrixFromTriplets();
    if ( !EigMat.isCompressed() ) {
        EigMat.makeCompressed();
    }

    const int r = i - 1, c = j - 1; // at() is 1-based, Eigen is 0-based
    if ( !this->hasStoredEntry( r, c ) ) {
        // Inserting makes the matrix uncompressed; it is compressed again before the reference is taken,
        // so that the reference stays valid until the next change of the structure.
        EigMat.insert( r, c ) = 0.;
        EigMat.makeCompressed();
    }

    // The caller may change the value through the returned reference.
    this->version++;
    return EigMat.coeffRef( r, c );
}


double EigenMtrx ::at( int i, int j ) const
{
    this->checkBounds( i, j );
    return this->giveMatrix().coeff( i - 1, j - 1 ); // at() is 1-based, Eigen is 0-based
}

const Eigen::SparseMatrix<double> &EigenMtrx::giveMatrix() const
{
    this->createMatrixFromTriplets(); // add the pending contributions
    return EigMat;
}


template <typename Derived>
Eigen::SparseSolverBase<Derived> &EigenMtrx::giveFactorization( FactorizationType Factorization )
{
    if ( !this->isFactorized( Factorization ) ) {
      this->computeFactorization( Factorization );
    }

    switch ( Factorization ) {
    case FT_LLT:
        return this->LLT_factorization;
        break;
    case FT_LU:
        return this->LU_factorization;
        break;
    case FT_QR:
        return this->QR_factorization;
        break;
    case FT_LDLT:
        return this->LDLT_factorization;
        break;
    default:
        OOFEM_ERROR( "Unknown factorization type" );
    }
}


SimplicialLDLTderived<Eigen::SparseMatrix<double> > &EigenMtrx::giveLDLTFactorization()
{
    if ( !this->isFactorized( FT_LDLT ) ) {
        this->computeFactorization( FT_LDLT );
    }

    return this->LDLT_factorization;
}


bool EigenMtrx::isFactorized( FactorizationType factorizationType ) const
{
    if ( (int)factorizationType < (int)FT_LLT || (int)factorizationType > (int)FT_LDLT ) {
        OOFEM_ERROR( "Unknown factorization type" );
    }
    // each factorization remembers the version it was computed for, so a change of the values
    // (or the factorization of another type) can never leave a stale one marked as current
    return this->factorizedAtVersion[factorizationType] == this->version;
}


void EigenMtrx::computeFactorization( FactorizationType factorizationType )
{
    const Eigen::SparseMatrix<double> &A = this->giveMatrix();
    switch ( factorizationType ) {
    case FT_LLT:
        this->LLT_factorization.compute( A );
        break;
    case FT_LU:
        this->LU_factorization.compute( A );
        break;
    case FT_QR:
        this->QR_factorization.compute( A );
        break;
    case FT_LDLT:
        this->LDLT_factorization.compute( A );
        break;
    default:
        OOFEM_ERROR( "Unknown factorization type" );
    }
    // recorded also when the factorization failed (its info() says so), so it is not recomputed in vain
    this->factorizedAtVersion[factorizationType] = this->version;
}


//void EigenMtrx::printYourself() const
//{
//    // Create copy
//    Eigen::SparseMatrix<double> EigMatCopy( this->EigMat );
//    if ( this->BuiltFromTripletsAtVersion != this->version ) {
//        EigMatCopy.setFromTriplets( this->triplets.begin(), this->triplets.end() );
//    }
//
//    // create dense matrix
//    auto denseMat = Eigen::MatrixXd( EigMatCopy );
//    std::cout << denseMat << std::endl;
//}


void EigenMtrx::times( const FloatArray &x, FloatArray &answer ) const
{
    if ( x.giveSize() != this->giveNumberOfColumns() ) {
        OOFEM_ERROR( "dimension mismatch: matrix %d x %d, vector %d", this->giveNumberOfRows(), this->giveNumberOfColumns(), x.giveSize() );
    }

    // no copy of the matrix or of x
    Eigen::Map<const Eigen::VectorXd> xeig( x.givePointer(), x.giveSize() );
    Eigen::VectorXd answereig = this->giveMatrix() * xeig;

    answer                    = FloatArray( answereig.begin(), answereig.end() );
}

void EigenMtrx::createMatrixFromTriplets() const
{
    if ( !this->tripletsPending ) {
        return;
    }

    if ( this->EigMat.nonZeros() == 0 ) {
        // the usual case after zero(): the matrix is the sum of the assembled contributions
        this->EigMat.setFromTriplets( this->triplets.begin(), this->triplets.end() );
    } else {
        // contributions assembled into a matrix that already has values: they are added to them
        Eigen::SparseMatrix<double> contributions( this->EigMat.rows(), this->EigMat.cols() );
        contributions.setFromTriplets( this->triplets.begin(), this->triplets.end() );
        this->EigMat += contributions;
    }
    this->EigMat.makeCompressed();

    this->triplets.clear(); // the capacity is kept for the next assembly
    this->tripletsPending = false;
}

void EigenMtrx::printYourself() const
{
    std::cout << this->giveMatrix() << std::endl;
}

//void EigenMtrx::getBlock( int i, int j, int p, int q, EigenMtrx& BlockMatrix )
std::unique_ptr<EigenMtrx> EigenMtrx::getBlock( int i, int j, int p, int q ) const
{
    // i and j start from 1
    if ( i < 1 || j < 1 || p < 0 || q < 0 || i - 1 + p > this->giveNumberOfRows() || j - 1 + q > this->giveNumberOfColumns() ) {
        OOFEM_ERROR( "block at (%d,%d) of size %d x %d is outside the %d x %d matrix", i, j, p, q, this->giveNumberOfRows(), this->giveNumberOfColumns() );
    }
    Eigen::SparseMatrix<double> subMat = this->giveMatrix().block( i - 1, j - 1, p, q );
    return std::make_unique<EigenMtrx>( std::move( subMat ) );
}

Eigen::SparseMatrix<double> EigenMtrx::computeSchurComplement( int i ) const
{
    // i is position were the PP block starts
    int matSize = this->giveNumberOfRows();
    if ( i < 1 || i > matSize + 1 ) {
        OOFEM_ERROR( "the condensed block must start at an equation from 1 to %d, got %d", matSize + 1, i );
    }
    int uuSize = i - 1;
    int ppSize  = matSize - uuSize;

    const Eigen::SparseMatrix<double> &eigmat = this->giveMatrix();

    Eigen::SparseMatrix<double> UU = eigmat.block( 0, 0, uuSize, uuSize );
    if ( ppSize == 0 ) {
        return UU; // nothing to condense
    }
    Eigen::SparseMatrix<double> UP = eigmat.block( 0, i - 1, uuSize, ppSize );
    Eigen::SparseMatrix<double> PU = eigmat.block( i - 1, 0, ppSize, uuSize );
    Eigen::SparseMatrix<double> PP = eigmat.block( i - 1, i - 1, ppSize, ppSize );

    /////////////
    //std::cout << UU << std::endl;
    //std::cout << UP << std::endl;
    //std::cout << PU << std::endl;
    //std::cout << PP << std::endl;
    //std::cout << this->EigMat << std::endl;
    ////////////////

    // K_PP^-1 K_PU: LDL^T (no pivoting) breaks down on a zero pivot, which a nonsingular indefinite
    // K_PP may still have; LU with pivoting is used then. Without a check, a failed factorization
    // would silently leave K_UU (or garbage) as the result.
    Eigen::SparseMatrix<double> temp;
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > A_factorization( PP );
    if ( A_factorization.info() == Eigen::Success ) {
        temp = A_factorization.solve( PU ); // Solve the system
    } else {
        Eigen::SparseLU<Eigen::SparseMatrix<double> > LU_PP( PP );
        if ( LU_PP.info() != Eigen::Success ) {
            OOFEM_ERROR( "block K_PP of the equations %d to %d is singular: the Schur complement does not exist", i, matSize );
        }
        temp = LU_PP.solve( PU );
    }

    return UU - UP * temp;
}

std::unique_ptr<EigenMtrx> EigenMtrx::doStaticCondensation( int i ) const
{
    // the condensed matrix is a ready matrix: wrapped by the constructor that takes it as built
    return std::make_unique<EigenMtrx>( this->computeSchurComplement( i ) );
}

void EigenMtrx::doStaticCondensationOwn( int i)
{
    Eigen::SparseMatrix<double> Kcond = this->computeSchurComplement( i );

    // the whole state changes consistently: values, size, no pending contributions, new version
    this->EigMat = std::move( Kcond );
    this->EigMat.makeCompressed();
    this->triplets.clear();
    this->tripletsPending = false;
    this->nRows = this->nColumns = (int)this->EigMat.rows();
    this->version++;
}

std::unique_ptr<EigenMtrx> EigenMtrx::doStaticCondensationLagrange( int i ) const
{
    // i is position were the PP block starts
    int matSize = this->giveNumberOfRows();
    if ( i < 1 || i > matSize + 1 ) {
        OOFEM_ERROR( "the multiplier block must start at an equation from 1 to %d, got %d", matSize + 1, i );
    }
    int uuSize  = i - 1;

    Eigen::SparseMatrix<double> UU = this->giveMatrix().block( 0, 0, uuSize, uuSize );
    //Eigen::SparseMatrix<double> UP = this->EigMat.block( 0, i - 1, uuSize, ppSize );
    //Eigen::SparseMatrix<double> PU = this->EigMat.block( i - 1, 0, ppSize, uuSize );
    //Eigen::SparseMatrix<double> PP = this->EigMat.block( i - 1, i - 1, ppSize, ppSize );

    //// create factorization and solve for lambda
    //Eigen::SparseLU<Eigen::SparseMatrix<double> > A_factorization( PP );
    //Eigen::SparseMatrix<double> temp  = A_factorization.solve( PU ); // Solve the system
    //Eigen::SparseMatrix<double> Kcond = UU - UP * temp;

    // std::cout << PP << std::endl;
    return std::make_unique<EigenMtrx>( std::move( UU ) );

}

//#include "eigenmtrx.h" // Adjust header names as needed
//#include "intarray.h"
//#include <memory>
//#include <vector>
//#include <algorithm>

// Assuming EigenMtrx is a wrapper containing an Eigen::SparseMatrix<double>
// named mMatrix (or accessible via a getter). Adjust member names accordingly.

//std::unique_ptr<EigenMtrx> EigenMtrx::removeRowsAndColumns( const IntArray &indsToDelete )
//{
//
//    // Access the underlying Eigen sparse matrix.
//    // Replace 'this->mMatrix' with your actual internal member or getter.
//    //const Eigen::SparseMatrix<double> &A = this->mMatrix;
//    auto A = this->giveMatrix();
//
//    int old_rows = A.rows();
//    int old_cols = A.cols();
//
//    // 1. Create a boolean mask of indices to keep
//    std::vector<bool> keep_idx( std::max( old_rows, old_cols ), true );
//
//    // OOFEM IntArray uses 1-based indexing (1 to giveSize())
//    for ( int i = 1; i <= indsToDelete.giveSize(); ++i ) {
//        // NOTE: Assuming the values inside indsToDelete are 1-based OOFEM
//        // equation/DOF numbers. We subtract 1 to get 0-based C++ indices.
//        // If your array is already 0-based, remove the '- 1'.
//        int idx = indsToDelete.at( i ) - 1;
//
//        if ( idx >= 0 && idx < keep_idx.size() ) {
//            keep_idx[idx] = false;
//        }
//    }
//
//    // 2. Create mappings from old to new indices
//    std::vector<int> row_map( old_rows, -1 );
//    int new_rows = 0;
//    for ( int i = 0; i < old_rows; ++i ) {
//        if ( keep_idx[i] ) row_map[i] = new_rows++;
//    }
//
//    std::vector<int> col_map( old_cols, -1 );
//    int new_cols = 0;
//    for ( int j = 0; j < old_cols; ++j ) {
//        if ( keep_idx[j] ) col_map[j] = new_cols++;
//    }
//
//    // 3. Count non-zeros per new column for exact memory allocation
//    Eigen::VectorXi nnz_per_col = Eigen::VectorXi::Zero( new_cols );
//    for ( int k = 0; k < A.outerSize(); ++k ) {
//        if ( col_map[k] == -1 ) continue; // Skip removed columns
//
//        for ( Eigen::SparseMatrix<double>::InnerIterator it( A, k ); it; ++it ) {
//            if ( row_map[it.row()] != -1 ) { // Skip removed rows
//                nnz_per_col[col_map[k]]++;
//            }
//        }
//    }
//
//    // 4. Build the new Eigen matrix
//    Eigen::SparseMatrix<double> B( new_rows, new_cols );
//    B.reserve( nnz_per_col );
//
//    for ( int k = 0; k < A.outerSize(); ++k ) {
//        if ( col_map[k] == -1 ) continue;
//        int new_col = col_map[k];
//
//        for ( Eigen::SparseMatrix<double>::InnerIterator it( A, k ); it; ++it ) {
//            int new_row = row_map[it.row()];
//            if ( new_row != -1 ) {
//                // insert() is O(1) here because elements are inserted in order
//                B.insert( new_row, new_col ) = it.value();
//            }
//        }
//    }
//
//    B.makeCompressed();
//
//    // 5. Wrap it in a new OOFEM EigenMtrx and return
//    //auto condensedMtrx = std::make_unique<EigenMtrx>();
//
//    // Assign the new Eigen matrix to the OOFEM wrapper.
//    // Replace 'mMatrix' with your actual internal member or use a setter.
//    //condensedMtrx->mMatrix = std::move( B );
//    EigenMtrx condensedMtrx( B );
//
//    //return std::make_unique<EigenMtrx>( condensedMtrx );
//    return std::make_unique<EigenMtrx>( std::move( B ) );
//
//}

std::unique_ptr<EigenMtrx> EigenMtrx::removeRowsAndColumns( const IntArray &indsToDelete ) const
{
    const Eigen::SparseMatrix<double> &A = this->giveMatrix();

    const int n = (int)A.rows();

    std::vector<bool> remove( n, false );
    for ( int i = 1; i <= indsToDelete.giveSize(); ++i ) {
        int idx = indsToDelete.at( i );

        if ( idx < 1 || idx > n ) {
            OOFEM_ERROR( "Index %d out of range.", idx );
        }

        remove[idx - 1] = true;
    }

    std::vector<int> map( n, -1 );
    int newSize = 0;
    for ( int i = 0; i < n; ++i ) {
        if ( !remove[i] ) {
            map[i] = newSize++;
        }
    }

    Eigen::SparseMatrix<double> reduced( newSize, newSize );

    std::vector<Eigen::Triplet<double> > reducedTriplets;
    reducedTriplets.reserve( A.nonZeros() );

    for ( int k = 0; k < A.outerSize(); ++k ) {
        for ( Eigen::SparseMatrix<double>::InnerIterator it( A, k ); it; ++it ) {

            if ( remove[it.row()] || remove[it.col()] ) {
                continue;
            }

            reducedTriplets.emplace_back(
                map[it.row()],
                map[it.col()],
                it.value() );
        }
    }

    reduced.setFromTriplets( reducedTriplets.begin(), reducedTriplets.end() );

    // a ready matrix, wrapped by the constructor that takes it as built
    return std::make_unique<EigenMtrx>( std::move( reduced ) );
}

#include <fstream>

void EigenMtrx::printToFile( const std::string &filename ) const
{
    std::ofstream file( filename );

    if ( !file.is_open() ) {
        OOFEM_ERROR( "Cannot open file '%s' for writing.", filename.c_str() );
    }
    file << this->giveMatrix() << std::endl; // pending contributions are added first

    file.close();
}

//void EigenMtrx::printToFile( const std::string &filename ) const
//{
//    std::ofstream file( filename );
//
//    if ( !file.is_open() ) {
//        OOFEM_ERROR( "Cannot open file '%s' for writing.", filename.c_str() );
//    }
//
//    file << "% Rows " << EigMat.rows()
//         << " Cols " << EigMat.cols()
//         << " NNZ " << EigMat.nonZeros() << '\n';
//
//    for ( int k = 0; k < EigMat.outerSize(); ++k ) {
//        for ( Eigen::SparseMatrix<double>::InnerIterator it( EigMat, k ); it; ++it ) {
//            file << it.row() + 1 << " "
//                 << it.col() + 1 << " "
//                 << std::setprecision( 16 ) << it.value() << '\n';
//        }
//    }
//}

} // end namespace oofem
