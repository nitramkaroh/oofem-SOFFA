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


#pragma once

#include "sparsemtrx.h"
#include "intarray.h"

#include "symildlwrapper.h"

#define _IFT_SymildlMtrx_Name "SymildlMtrx"

namespace oofem {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Triplet {
  int row;
  int col;
  double value;
};

struct CSCMatrix {
  int n;
  std::vector<int> Ap; // column pointers
  std::vector<int> Ai; // row indices
  std::vector<double> Ax; // values
};

class OOFEM_EXPORT SymildlMtrx : public SparseMtrx
{

  protected:

  std::vector<Triplet> triplets; // vector of triplets
  CSCMatrix csc;

  public:
  /** Constructor. Before any operation an internal profile must be built.
   * @see buildInternalStructure
   */
  SymildlMtrx( int n = 0 );

  SymildlMtrx( const SymildlMtrx &symildlMtrxInput );

  /// Destructor
  virtual ~SymildlMtrx() {}

  // Overloaded methods:
  int buildInternalStructure( EngngModel *, int, const UnknownNumberingScheme &s ) override;
  int assemble( const IntArray &loc, const FloatMatrix &mat ) override;
  int assemble( const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat ) override;
  bool canBeFactorized() const override { return false; }
  void zero() override;
  double &at( int i, int j ) override;
  double at( int i, int j ) const override;
  
  const char *giveClassName() const override { return "SymildlMtrx"; }
  SparseMtrxType giveType() const override { return SMT_SymildlMtrx; }
  bool isAsymmetric() const override { return false; }

  CSCMatrix &giveMatrix();

  //template <typename Derived>
  //Eigen::SparseSolverBase<Derived> &giveFactorization( FactorizationType Factorization );

  void times( const FloatArray &x, FloatArray &answer ) const override;
  void times( double factor ) override;

  void printYourself() const override;

  void printTriplets(const std::vector<Triplet> triplets) const;

  protected:
  // takes a vector of triplets and makes a CSC matrix representation
  CSCMatrix tripletToCSC( int n, const std::vector<Triplet> &triplets ) const;

  //checks given triplets for duplicates, adds them together, returns new triplet field
  std::vector<Triplet> unifyTriplets(std::vector<Triplet> input) const;

};
} // end namespace oofem
