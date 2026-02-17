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


#include "symildlmtrx.h"

#include "floatarray.h"
#include "engngm.h"
#include "domain.h"
#include "element.h"
#include "sparsemtrxtype.h"
#include "activebc.h"
#include "classfactory.h"
#include <iostream>

#include <set>


namespace oofem {
REGISTER_SparseMtrx( SymildlMtrx, SMT_SymildlMtrx );


SymildlMtrx ::SymildlMtrx( int n ) :
    SparseMtrx( n, n )
{
}

SymildlMtrx ::SymildlMtrx( const SymildlMtrx &symildlMtrxInput ) :
    SparseMtrx( symildlMtrxInput.giveNumberOfRows(), symildlMtrxInput.giveNumberOfColumns() )
{
}

int SymildlMtrx ::buildInternalStructure( EngngModel *eModel, int di, const UnknownNumberingScheme &s )
{
  IntArray loc;
  Domain *domain = eModel->giveDomain( di );
  int neq = eModel->giveNumberOfDomainEquations( di, s );

  triplets.reserve( 500 * neq );

  nRows = nColumns = neq;

  return true;
}


int SymildlMtrx ::assemble( const IntArray &loc, const FloatMatrix &mat )
{
  int dim = mat.giveNumberOfRows();
  int ii, jj;
  for ( int j = 0; j < dim; j++ ) {
    jj = loc[j];
    if ( jj ) {
      for ( int i = 0; i < dim; i++ ) {
        ii = loc[i];
        if ( ii ) {
          // set from triplets
          Triplet t;
          t.row = ii - 1;
          t.col = jj - 1;
          t.value = mat( i, j );
          this->triplets.push_back( t ); // Add to vector of triplets
        }
      }
    }
  }

  return 1;
}

int SymildlMtrx ::assemble( const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat )
{
  for ( int j = 0; j < mat.giveNumberOfColumns(); j++ ) {
    int jj = cloc[j];
    if ( jj ) {
      for ( int i = 0; i < mat.giveNumberOfRows(); i++ ) {
        int ii = rloc[i];
        if ( ii ) {
          // set from triplets
          Triplet t;
          t.row = ii - 1;
          t.col = jj - 1;
          t.value = mat( i, j );
          this->triplets.push_back( t ); // Add to vector of triplets
        }
      }
    }
  }

  return 1;
}

void SymildlMtrx ::zero()
{
  this->triplets.clear(); // When set from triplets
}

double &SymildlMtrx::at( int i, int j )
{
  if ( i > this->giveNumberOfRows() && j > this->giveNumberOfColumns() ) {
    OOFEM_ERROR( "Array accessing exception -- (%d,%d) out of bounds", i, j );
  } else {
    for ( int it = 0; it < triplets.size(); it++ ) {
      auto t = triplets[it];
      if ( t.row == i && t.col == j ) {
        return t.value;
      }
    }
    OOFEM_ERROR( "Array accessing exception -- (%d,%d) not allocated", i, j );
  }
  
}


double SymildlMtrx::at( int i, int j ) const
{
  if ( i > this->giveNumberOfRows() && j > this->giveNumberOfColumns() ) {
    OOFEM_ERROR( "Array accessing exception -- (%d,%d) out of bounds", i, j );
  } else {
    for ( auto &t : triplets ) {
      if ( t.row == i && t.col == j ) {
        return t.value;
      }
    }
    return 0.0;
  }
}


void SymildlMtrx::times( const FloatArray &x, FloatArray &answer ) const
{
  if ( x.giveSize() != nColumns ) {
    OOFEM_ERROR( "Incompatible sizes for matrix multiplication: A[%i,%i] * x[%i]", nRows, nColumns, x.giveSize() );
  }
  answer.resize(nColumns);

  for ( auto &t : triplets ) {
    answer.at(t.row+1) += t.value*x.at(t.col+1);
  }
}

void SymildlMtrx::times( double x )
{
  for ( auto &t : triplets ) {
    t.value *= x;
  }
}

void SymildlMtrx::printYourself() const
{
  std::vector<Triplet> unique_triplets = unifyTriplets( this->triplets );
  this->printTriplets(unique_triplets);
}

void SymildlMtrx::printTriplets( const std::vector<Triplet> triplets ) const
{
  for ( auto &t : triplets ) {
    std::cout << "a(" << t.row + 1 << ", " << t.col + 1 << ") = " << t.value << std::endl;
  }
}

CSCMatrix& SymildlMtrx::giveMatrix(){
  if (nColumns != nRows) OOFEM_ERROR("CSC building only supported for symmetric sparse matrices");

  this->csc = tripletToCSC(nColumns, this->triplets);
  return this->csc;
}

CSCMatrix SymildlMtrx::tripletToCSC( int n, const std::vector<Triplet> &triplets ) const 
{

  // Select only lower triangle
  std::vector<Triplet> lt_triplets;

  for ( size_t k = 0; k < triplets.size(); ++k ) {
    auto &t = triplets[k];
    if ( t.row >= t.col ) {
      lt_triplets.push_back(t);
    }
  }

  // unify
  lt_triplets = unifyTriplets(lt_triplets);

  // Sort by (col, row)
  std::sort( lt_triplets.begin(), lt_triplets.end(),
      []( const Triplet &a, const Triplet &b ) {
        if ( a.col != b.col ) return a.col < b.col;
        return a.row < b.row;
      } );

  std::vector<int> Ap( n + 1, 0 );
  std::vector<int> Ai;
  std::vector<double> Ax;

  Ai.reserve( lt_triplets.size() );
  Ax.reserve( lt_triplets.size() );

  int current_col = 0;
  int nnz = 0;
  Ap[0] = 0;

  for ( size_t k = 0; k < lt_triplets.size(); ++k ) {
    const auto &t = lt_triplets[k];

    while ( current_col < t.col )
      Ap[++current_col] = nnz;

    Ai.push_back( t.row );
    Ax.push_back( t.value );
    ++nnz;
  }

  while ( current_col < n )
    Ap[++current_col] = nnz;

  return { n, std::move( Ap ), std::move( Ai ), std::move( Ax ) };
}

std::vector<Triplet> SymildlMtrx::unifyTriplets( std::vector<Triplet> input ) const
{
  std::vector<Triplet> output(0);
  bool found = false;

  for ( int iit = 0; iit < input.size(); iit++){
    //search output for existing row, col
    found = false;
    for ( int iot = 0; iot < output.size(); iot++ ) {
      if (output[iot].row == input[iit].row && output[iot].col == input[iit].col){
        output[iot].value += input[iit].value;
        found=true;
        break;
      }
    }
    if (!found){
      output.push_back(input[iit]);
    }
  }

  return output;

}


} // end namespace oofem
