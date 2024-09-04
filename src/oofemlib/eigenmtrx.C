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

int EigenMtrx ::buildInternalStructure( EngngModel *eModel, int di, const UnknownNumberingScheme &s )
{
    IntArray loc;
    Domain *domain = eModel->giveDomain( di );
    int neq        = eModel->giveNumberOfDomainEquations( di, s );

    EigMat.resize( neq, neq ); // Resize the matrix
    triplets.reserve( 500*neq );


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

                    // when set from triplets
                    this->triplets.push_back( Eigen::Triplet<double>( ii - 1, jj - 1, mat( i, j ) ) ); // Add to vector of triplets
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

                    // when set from triplets
                    this->triplets.push_back( Eigen::Triplet<double>( ii - 1, jj - 1, mat( i, j ) ) ); // Add to vector of triplets
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
        this->computeFactorization( Factorization )
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


void EigenMtrx::printYourself() const
{
    // Create copy
    Eigen::SparseMatrix<double> EigMatCopy( this->EigMat );
    if ( this->BuiltFromTripletsAtVersion != this->version ) {
        EigMatCopy.setFromTriplets( this->triplets.begin(), this->triplets.end() );
    }

    // create dense matrix
    auto denseMat = Eigen::MatrixXd( EigMatCopy );
    std::cout << denseMat << std::endl;
}


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
    if ( this->BuiltFromTripletsAtVersion != this->version ) {
        this->EigMat.setFromTriplets( this->triplets.begin(), this->triplets.end() );
        this->EigMat.makeCompressed();
        BuiltFromTripletsAtVersion = this->version;
    }
}

} // end namespace oofem