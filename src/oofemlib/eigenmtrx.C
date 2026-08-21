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

#include <set>


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

EigenMtrx ::EigenMtrx( Eigen::SparseMatrix<double> &EigMatInput ) :
    SparseMtrx( EigMatInput.rows(), EigMatInput.cols() ),
    EigMat( EigMatInput )
{
}

EigenMtrx ::EigenMtrx( const EigenMtrx &EigenMtrxInput ) :
    SparseMtrx( EigenMtrxInput.giveNumberOfRows(), EigenMtrxInput.giveNumberOfColumns() ),
    EigMat( EigenMtrxInput.EigMat ),
    triplets( EigenMtrxInput.triplets )
{
    this->version = EigenMtrxInput.version;
    this->createMatrixFromTriplets();

}

//EigenMtrx &EigenMtrx ::operator=( Eigen::SparseMatrix<double> &EigMat )
//{
//}

int EigenMtrx ::buildInternalStructure( EngngModel *eModel, int di, const UnknownNumberingScheme &s )
{
    IntArray loc;
    Domain *domain = eModel->giveDomain( di );
    int neq        = eModel->giveNumberOfDomainEquations( di, s );

    EigMat.resize( neq, neq ); // Resize the matrix
    triplets.reserve( 500*neq );

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
                    // of magnitude. Diagonal entries are kept unconditionally: at() maps
                    // onto coeffRef, which would have to insert into a compressed matrix
                    // if the entry were missing.
                    if ( mat( i, j ) != 0. || ii == jj ) {
                        // when set from triplets
                        this->triplets.push_back( Eigen::Triplet<double>( ii - 1, jj - 1, mat( i, j ) ) ); // Add to vector of triplets
                    }
                }
            }
        }
    }

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
                    // of magnitude. Diagonal entries are kept unconditionally: at() maps
                    // onto coeffRef, which would have to insert into a compressed matrix
                    // if the entry were missing.
                    if ( mat( i, j ) != 0. || ii == jj ) {
                        // when set from triplets
                        this->triplets.push_back( Eigen::Triplet<double>( ii - 1, jj - 1, mat( i, j ) ) ); // Add to vector of triplets
                    }
                }
            }
        }
    }

    this->version++;
    return 1;
}

void EigenMtrx ::zero()
{
    EigMat.setZero();
    this->triplets.clear(); // When set from triplets
    this->version++;
}

double &EigenMtrx ::at( int i, int j )
{
    this->version++;
    if ( i > this->giveNumberOfRows() && j > this->giveNumberOfColumns() ) {
        OOFEM_ERROR( "Array accessing exception -- (%d,%d) out of bounds", i, j );
    } else {
        return this->giveMatrix().coeffRef( i, j );
    }
}


double EigenMtrx ::at( int i, int j ) const
{
    // Create copy
    Eigen::SparseMatrix<double> EigMatCopy( this->EigMat );
    if ( this->BuiltFromTripletsAtVersion != this->version ) {
        EigMatCopy.setFromTriplets( this->triplets.begin(), this->triplets.end() );
    }

    if ( i > this->giveNumberOfRows() && j > this->giveNumberOfColumns() ) {
        OOFEM_ERROR( "Array accessing exception -- (%d,%d) out of bounds", i, j );
    } else {
        //return EigMat.coeff( i, j );
        return EigMatCopy.coeff( i, j );
    }
}

Eigen::SparseMatrix<double>& EigenMtrx::giveMatrix()
{
    this->createMatrixFromTriplets(); // When the matrix is created from triplets
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


bool EigenMtrx::isFactorized( FactorizationType factorizationType )
{
    if ( this->versionUpdate != this->version ) {
        areFactorized = { false, false, false, false };
    }
    return this->areFactorized[factorizationType];
}


void EigenMtrx::computeFactorization( FactorizationType factorizationType )
{
    switch ( factorizationType ) {
    case FT_LLT:
        this->LLT_factorization.compute( this->giveMatrix() );
        break;
    case FT_LU:
        this->LU_factorization.compute( this->giveMatrix() );
        break;
    case FT_QR:
        this->QR_factorization.compute( this->giveMatrix() );
        break;
    case FT_LDLT:
        this->LDLT_factorization.compute( this->giveMatrix() );
        break;
    default:
        OOFEM_ERROR( "Unknown factorization type" );
    }
    this->versionUpdate                    = this->version;
    this->areFactorized[factorizationType] = true;
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
    // Create a copy
    Eigen::SparseMatrix<double> EigMatCopy(this->EigMat);
    if ( this->BuiltFromTripletsAtVersion != this->version ) {
        EigMatCopy.setFromTriplets( this->triplets.begin(), this->triplets.end() );
    }

    FloatArray xcopy          = x;
    Eigen::VectorXd xeig      = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>( xcopy.givePointer(), x.giveSize() );

    //Eigen::VectorXd answereig = this->EigMat * xeig;
    Eigen::VectorXd answereig = EigMatCopy * xeig; // If set from triplets

    answer                    = FloatArray( answereig.begin(), answereig.end() );
}

void EigenMtrx::createMatrixFromTriplets()
{
    if ( this->BuiltFromTripletsAtVersion != this->version || this->version  == 0) {
        this->EigMat.setFromTriplets( this->triplets.begin(), this->triplets.end() );
        this->EigMat.makeCompressed();
        BuiltFromTripletsAtVersion = this->version;
    }
}

void EigenMtrx::printYourself() const
{
    std::cout << this->EigMat << std::endl;
}

//void EigenMtrx::getBlock( int i, int j, int p, int q, EigenMtrx& BlockMatrix )
std::unique_ptr<EigenMtrx> EigenMtrx::getBlock( int i, int j, int p, int q )
{
    // i and j start from 1
    Eigen::SparseMatrix<double> subMat = this->EigMat.block(i-1,j-1, p, q );
    EigenMtrx BlockMatrix( subMat );
    return std::make_unique<EigenMtrx>( BlockMatrix );
}

std::unique_ptr<EigenMtrx> EigenMtrx::doStaticCondensation( int i )
{
    // i is position were the PP block starts
    int uuSize = i - 1;
    int matSize = this->giveNumberOfRows();
    int ppSize  = matSize - uuSize;

    auto eigmat = this->giveMatrix();

    Eigen::SparseMatrix<double> UU = eigmat.block( 0, 0, uuSize, uuSize );
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

    // create factorization
    //Eigen::SparseLU<Eigen::SparseMatrix<double> > A_factorization( PP );
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > A_factorization( PP );

    Eigen::SparseMatrix<double> temp = A_factorization.solve( PU ); // Solve the system
    Eigen::SparseMatrix<double> Kcond = UU - UP * temp;
    EigenMtrx KcondOOFEM( Kcond );

    //std::cout << PP << std::endl;
    return std::make_unique<EigenMtrx>( KcondOOFEM );
}

void EigenMtrx::doStaticCondensationOwn( int i)
{
    // i is position were the PP block starts
    int uuSize = i - 1;
    int matSize = this->giveNumberOfRows();
    int ppSize  = matSize - uuSize;

    auto eigmat = this->giveMatrix();

    Eigen::SparseMatrix<double> UU = eigmat.block( 0, 0, uuSize, uuSize );
    Eigen::SparseMatrix<double> UP = eigmat.block( 0, i - 1, uuSize, ppSize );   
    Eigen::SparseMatrix<double> PU = eigmat.block( i - 1, 0, ppSize, uuSize );
    Eigen::SparseMatrix<double> PP = eigmat.block( i - 1, i - 1, ppSize, ppSize );

    // create factorization
    Eigen::SparseLU<Eigen::SparseMatrix<double> > A_factorization( PP );
    Eigen::SparseMatrix<double> temp = A_factorization.solve( PU ); // Solve the system
    Eigen::SparseMatrix<double> Kcond = UU - UP * temp;

    this->EigMat = Kcond;
}

std::unique_ptr<EigenMtrx> EigenMtrx::doStaticCondensationLagrange( int i )
{
    // i is position were the PP block starts
    int uuSize  = i - 1;
    int matSize = this->giveNumberOfRows();
    int ppSize  = matSize - uuSize;
    auto eigmat                    = this->giveMatrix();

    Eigen::SparseMatrix<double> UU = eigmat.block( 0, 0, uuSize, uuSize );
    //Eigen::SparseMatrix<double> UP = this->EigMat.block( 0, i - 1, uuSize, ppSize );
    //Eigen::SparseMatrix<double> PU = this->EigMat.block( i - 1, 0, ppSize, uuSize );
    //Eigen::SparseMatrix<double> PP = this->EigMat.block( i - 1, i - 1, ppSize, ppSize );

    //// create factorization and solve for lambda
    //Eigen::SparseLU<Eigen::SparseMatrix<double> > A_factorization( PP );
    //Eigen::SparseMatrix<double> temp  = A_factorization.solve( PU ); // Solve the system
    //Eigen::SparseMatrix<double> Kcond = UU - UP * temp;
    EigenMtrx KcondOOFEM( UU );

    // std::cout << PP << std::endl;
    return std::make_unique<EigenMtrx>( KcondOOFEM );

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

std::unique_ptr<EigenMtrx> EigenMtrx::removeRowsAndColumns( const IntArray &indsToDelete )
{
    Eigen::SparseMatrix<double> &A = this->giveMatrix();

    const int n = A.rows();

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

    auto result = std::make_unique<EigenMtrx>( newSize );
    result->version = 1;

    Eigen::SparseMatrix<double> reduced( newSize, newSize );

    std::vector<Eigen::Triplet<double> > triplets;
    triplets.reserve( A.nonZeros() );

    for ( int k = 0; k < A.outerSize(); ++k ) {
        for ( Eigen::SparseMatrix<double>::InnerIterator it( A, k ); it; ++it ) {

            if ( remove[it.row()] || remove[it.col()] ) {
                continue;
            }

            triplets.emplace_back(
                map[it.row()],
                map[it.col()],
                it.value() );
        }
    }

    reduced.setFromTriplets( triplets.begin(), triplets.end() );
    reduced.makeCompressed();

    result->giveMatrix().swap( reduced );
    result->triplets = triplets;
    result->version  = 1;
    result->BuiltFromTripletsAtVersion = result->version;
    

    return result;
}

#include <fstream>

void EigenMtrx::printToFile( const std::string &filename ) 
{
    std::ofstream file( filename );

    if ( !file.is_open() ) {
        OOFEM_ERROR( "Cannot open file '%s' for writing.", filename.c_str() );
    }
    this->createMatrixFromTriplets(); // Ensure the matrix is up-to-date
    file << this->EigMat << std::endl;

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
