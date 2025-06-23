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

#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "index.h"
#include "tensor2.h"
#pragma once

using namespace FTensor;



/**
 * Functions providing basic operations with tensors
 * based on FTensor library
 * It provides transformation from the second and fourth order tensors in 3d to FloatArrayF<9> and FloatArray<9,9>
 * It provides calculation of determinant, cofactor, inverse, and tensor cross product
 * More to be added gradually
 * @author Ondrej Faltus
 **/

namespace oofem {
class Tensor3_3d : public Tensor3<double, 3, 3, 3>
{
public:
    using Tensor3<double, 3, 3, 3>::Tensor3;

    /**
     * Creates second-order tensor in 3d
     * Initialized with zeros
     */

    Tensor3_3d()
    {
        this->operator()( 0, 0, 0 ) = 0;
        this->operator()( 0, 0, 1 ) = 0.;
        this->operator()( 0, 0, 2 ) = 0.;
        this->operator()( 0, 1, 0 ) = 0.;
        this->operator()( 0, 1, 1 ) = 0.;
        this->operator()( 0, 1, 2 ) = 0.;
        this->operator()( 0, 2, 0 ) = 0.;
        this->operator()( 0, 2, 1 ) = 0.;
        this->operator()( 0, 2, 2 ) = 0.;
        this->operator()( 1, 0, 0 ) = 0;
        this->operator()( 1, 0, 1 ) = 0.;
        this->operator()( 1, 0, 2 ) = 0.;
        this->operator()( 1, 1, 0 ) = 0.;
        this->operator()( 1, 1, 1 ) = 0.;
        this->operator()( 1, 1, 2 ) = 0.;
        this->operator()( 1, 2, 0 ) = 0.;
        this->operator()( 1, 2, 1 ) = 0.;
        this->operator()( 1, 2, 2 ) = 0.;
        this->operator()( 2, 0, 0 ) = 0;
        this->operator()( 2, 0, 1 ) = 0.;
        this->operator()( 2, 0, 2 ) = 0.;
        this->operator()( 2, 1, 0 ) = 0.;
        this->operator()( 2, 1, 1 ) = 0.;
        this->operator()( 2, 1, 2 ) = 0.;
        this->operator()( 2, 2, 0 ) = 0.;
        this->operator()( 2, 2, 1 ) = 0.;
        this->operator()( 2, 2, 2 ) = 0.;
    }

    /**
     * Creates a second-order order tensor in 3d from floatmatrixf<3, 9>
     */

    Tensor3_3d( const FloatMatrixF<3, 9> &matrix )
    {

        this->operator()( 0, 0, 0 ) = matrix.at( 1, 1 );
        this->operator()( 0, 0, 1 ) = matrix.at( 1, 6 );
        this->operator()( 0, 0, 2 ) = matrix.at( 1, 5 );
        this->operator()( 0, 1, 0 ) = matrix.at( 1, 9 );
        this->operator()( 0, 1, 1 ) = matrix.at( 1, 2 );
        this->operator()( 0, 1, 2 ) = matrix.at( 1, 4 );
        this->operator()( 0, 2, 0 ) = matrix.at( 1, 8 );
        this->operator()( 0, 2, 1 ) = matrix.at( 1, 7 );
        this->operator()( 0, 2, 2 ) = matrix.at( 1, 3 );

        this->operator()( 1, 0, 0 ) = matrix.at( 2, 1 );
        this->operator()( 1, 0, 1 ) = matrix.at( 2, 6 );
        this->operator()( 1, 0, 2 ) = matrix.at( 2, 5 );
        this->operator()( 1, 1, 0 ) = matrix.at( 2, 9 );
        this->operator()( 1, 1, 1 ) = matrix.at( 2, 2 );
        this->operator()( 1, 1, 2 ) = matrix.at( 2, 4 );
        this->operator()( 1, 2, 0 ) = matrix.at( 2, 8 );
        this->operator()( 1, 2, 1 ) = matrix.at( 2, 7 );
        this->operator()( 1, 2, 2 ) = matrix.at( 2, 3 );

        this->operator()( 2, 0, 0 ) = matrix.at( 3, 1 );
        this->operator()( 2, 0, 1 ) = matrix.at( 3, 6 );
        this->operator()( 2, 0, 2 ) = matrix.at( 3, 5 );
        this->operator()( 2, 1, 0 ) = matrix.at( 3, 9 );
        this->operator()( 2, 1, 1 ) = matrix.at( 3, 2 );
        this->operator()( 2, 1, 2 ) = matrix.at( 3, 4 );
        this->operator()( 2, 2, 0 ) = matrix.at( 3, 8 );
        this->operator()( 2, 2, 1 ) = matrix.at( 3, 7 );
        this->operator()( 2, 2, 2 ) = matrix.at( 3, 3 );
    }

    /**
     * Creates a second-order order tensor in 3d from floatmatrixf<9, 3>
     */

    Tensor3_3d( const FloatMatrixF<9, 3> &matrix )
    {

        this->operator()( 0, 0, 0 ) = matrix.at( 1, 1 );
        this->operator()( 0, 1, 0 ) = matrix.at( 6, 1 );
        this->operator()( 0, 2, 0 ) = matrix.at( 5, 1 );
        this->operator()( 1, 0, 0 ) = matrix.at( 9, 1 );
        this->operator()( 1, 1, 0 ) = matrix.at( 2, 1 );
        this->operator()( 1, 2, 0 ) = matrix.at( 4, 1 );
        this->operator()( 2, 0, 0 ) = matrix.at( 8, 1 );
        this->operator()( 2, 1, 0 ) = matrix.at( 7, 1 );
        this->operator()( 2, 2, 0 ) = matrix.at( 3, 1 );

        this->operator()( 0, 0, 1 ) = matrix.at( 1, 2 );
        this->operator()( 0, 1, 1 ) = matrix.at( 6, 2 );
        this->operator()( 0, 2, 1 ) = matrix.at( 5, 2 );
        this->operator()( 1, 0, 1 ) = matrix.at( 9, 2 );
        this->operator()( 1, 1, 1 ) = matrix.at( 2, 2 );
        this->operator()( 1, 2, 1 ) = matrix.at( 4, 2 );
        this->operator()( 2, 0, 1 ) = matrix.at( 8, 2 );
        this->operator()( 2, 1, 1 ) = matrix.at( 7, 2 );
        this->operator()( 2, 2, 1 ) = matrix.at( 3, 2 );

        this->operator()( 0, 0, 2 ) = matrix.at( 1, 3 );
        this->operator()( 0, 1, 2 ) = matrix.at( 6, 3 );
        this->operator()( 0, 2, 2 ) = matrix.at( 5, 3 );
        this->operator()( 1, 0, 2 ) = matrix.at( 9, 3 );
        this->operator()( 1, 1, 2 ) = matrix.at( 2, 3 );
        this->operator()( 1, 2, 2 ) = matrix.at( 4, 3 );
        this->operator()( 2, 0, 2 ) = matrix.at( 8, 3 );
        this->operator()( 2, 1, 2 ) = matrix.at( 7, 3 );
        this->operator()( 2, 2, 2 ) = matrix.at( 3, 3 );
    }

     /**
     * Creates a second-order order tensor in 3d from floatarrayf<27>
     */

    Tensor3_3d( const FloatArrayF<27> &array )
    {

      this->operator()( 0, 0, 0 ) = array.at( 1 );
      this->operator()( 0, 1, 0 ) = array.at( 6 );
      this->operator()( 0, 2, 0 ) = array.at( 5 );
      this->operator()( 1, 0, 0 ) = array.at( 9 );
      this->operator()( 1, 1, 0 ) = array.at( 2 );
      this->operator()( 1, 2, 0 ) = array.at( 4 );
      this->operator()( 2, 0, 0 ) = array.at( 8 );
      this->operator()( 2, 1, 0 ) = array.at( 7 );
      this->operator()( 2, 2, 0 ) = array.at( 3 );

      this->operator()( 0, 0, 1 ) = array.at( 10 );
      this->operator()( 0, 1, 1 ) = array.at( 15 );
      this->operator()( 0, 2, 1 ) = array.at( 14 );
      this->operator()( 1, 0, 1 ) = array.at( 18 );
      this->operator()( 1, 1, 1 ) = array.at( 11 );
      this->operator()( 1, 2, 1 ) = array.at( 13 );
      this->operator()( 2, 0, 1 ) = array.at( 17 );
      this->operator()( 2, 1, 1 ) = array.at( 16 );
      this->operator()( 2, 2, 1 ) = array.at( 12 );

      this->operator()( 0, 0, 2 ) = array.at( 19 );
      this->operator()( 0, 1, 2 ) = array.at( 24 );
      this->operator()( 0, 2, 2 ) = array.at( 23 );
      this->operator()( 1, 0, 2 ) = array.at( 27 );
      this->operator()( 1, 1, 2 ) = array.at( 20 );
      this->operator()( 1, 2, 2 ) = array.at( 22 );
      this->operator()( 2, 0, 2 ) = array.at( 26 );
      this->operator()( 2, 1, 2 ) = array.at( 25 );
      this->operator()( 2, 2, 2 ) = array.at( 21 );
    }


    /**
     * Transforms a second-order tensor into a floatmatrixf<3, 9>,  using the Voigt notation
     */

    const inline FloatMatrixF<3, 9> to_voigt_form_3x9()
    {
        FloatMatrixF<3, 9> answer;
        answer.at( 1, 1 ) = this->operator()( 0, 0, 0 );
        answer.at( 1, 2 ) = this->operator()( 0, 1, 1 );
        answer.at( 1, 3 ) = this->operator()( 0, 2, 2 );
        answer.at( 1, 4 ) = this->operator()( 0, 1, 2 );
        answer.at( 1, 5 ) = this->operator()( 0, 0, 2 );
        answer.at( 1, 6 ) = this->operator()( 0, 0, 1 );
        answer.at( 1, 7 ) = this->operator()( 0, 2, 1 );
        answer.at( 1, 8 ) = this->operator()( 0, 2, 0 );
        answer.at( 1, 9 ) = this->operator()( 0, 1, 0 );

        answer.at( 2, 1 ) = this->operator()( 1, 0, 0 );
        answer.at( 2, 2 ) = this->operator()( 1, 1, 1 );
        answer.at( 2, 3 ) = this->operator()( 1, 2, 2 );
        answer.at( 2, 4 ) = this->operator()( 1, 1, 2 );
        answer.at( 2, 5 ) = this->operator()( 1, 0, 2 );
        answer.at( 2, 6 ) = this->operator()( 1, 0, 1 );
        answer.at( 2, 7 ) = this->operator()( 1, 2, 1 );
        answer.at( 2, 8 ) = this->operator()( 1, 2, 0 );
        answer.at( 2, 9 ) = this->operator()( 1, 1, 0 );

        answer.at( 3, 1 ) = this->operator()( 2, 0, 0 );
        answer.at( 3, 2 ) = this->operator()( 2, 1, 1 );
        answer.at( 3, 3 ) = this->operator()( 2, 2, 2 );
        answer.at( 3, 4 ) = this->operator()( 2, 1, 2 );
        answer.at( 3, 5 ) = this->operator()( 2, 0, 2 );
        answer.at( 3, 6 ) = this->operator()( 2, 0, 1 );
        answer.at( 3, 7 ) = this->operator()( 2, 2, 1 );
        answer.at( 3, 8 ) = this->operator()( 2, 2, 0 );
        answer.at( 3, 9 ) = this->operator()( 2, 1, 0 );

        return answer;
    }


    /**
     * Transforms a second-order tensor into a floatmatrixf<9, 3>,  using the Voigt notation
     */

    const inline FloatMatrixF<9, 3> to_voigt_form_9x3()
    {
        FloatMatrixF<9, 3> answer;
        answer.at( 1, 1 ) = this->operator()( 0, 0, 0 );
        answer.at( 2, 1 ) = this->operator()( 1, 1, 0 );
        answer.at( 3, 1 ) = this->operator()( 2, 2, 0 );
        answer.at( 4, 1 ) = this->operator()( 1, 2, 0 );
        answer.at( 5, 1 ) = this->operator()( 0, 2, 0 );
        answer.at( 6, 1 ) = this->operator()( 0, 1, 0 );
        answer.at( 7, 1 ) = this->operator()( 2, 1, 0 );
        answer.at( 8, 1 ) = this->operator()( 2, 0, 0 );
        answer.at( 9, 1 ) = this->operator()( 1, 0, 0 );

        answer.at( 1, 2 ) = this->operator()( 0, 0, 1 );
        answer.at( 2, 2 ) = this->operator()( 1, 1, 1 );
        answer.at( 3, 2 ) = this->operator()( 2, 2, 1 );
        answer.at( 4, 2 ) = this->operator()( 1, 2, 1 );
        answer.at( 5, 2 ) = this->operator()( 0, 2, 1 );
        answer.at( 6, 2 ) = this->operator()( 0, 1, 1 );
        answer.at( 7, 2 ) = this->operator()( 2, 1, 1 );
        answer.at( 8, 2 ) = this->operator()( 2, 0, 1 );
        answer.at( 9, 2 ) = this->operator()( 1, 0, 1 );

        answer.at( 1, 3 ) = this->operator()( 0, 0, 2 );
        answer.at( 2, 3 ) = this->operator()( 1, 1, 2 );
        answer.at( 3, 3 ) = this->operator()( 2, 2, 2 );
        answer.at( 4, 3 ) = this->operator()( 1, 2, 2 );
        answer.at( 5, 3 ) = this->operator()( 0, 2, 2 );
        answer.at( 6, 3 ) = this->operator()( 0, 1, 2 );
        answer.at( 7, 3 ) = this->operator()( 2, 1, 2 );
        answer.at( 8, 3 ) = this->operator()( 2, 0, 2 );
        answer.at( 9, 3 ) = this->operator()( 1, 0, 2 );

        return answer;
    }

    /**
     * Transforms a second-order tensor into a floatarrayf<27>,  using the Voigt notation
     */

    const inline FloatArrayF<27> to_voigt_form_27()
    {
      FloatArrayF<27> answer;
      answer.at( 1 ) = this->operator()( 0, 0, 0 );
      answer.at( 2 ) = this->operator()( 1, 1, 0 );
      answer.at( 3 ) = this->operator()( 2, 2, 0 );
      answer.at( 4 ) = this->operator()( 1, 2, 0 );
      answer.at( 5 ) = this->operator()( 0, 2, 0 );
      answer.at( 6 ) = this->operator()( 0, 1, 0 );
      answer.at( 7 ) = this->operator()( 2, 1, 0 );
      answer.at( 8 ) = this->operator()( 2, 0, 0 );
      answer.at( 9 ) = this->operator()( 1, 0, 0 );

      answer.at( 10 ) = this->operator()( 0, 0, 1 );
      answer.at( 11 ) = this->operator()( 1, 1, 1 );
      answer.at( 12 ) = this->operator()( 2, 2, 1 );
      answer.at( 13 ) = this->operator()( 1, 2, 1 );
      answer.at( 14 ) = this->operator()( 0, 2, 1 );
      answer.at( 15 ) = this->operator()( 0, 1, 1 );
      answer.at( 16 ) = this->operator()( 2, 1, 1 );
      answer.at( 17 ) = this->operator()( 2, 0, 1 );
      answer.at( 18 ) = this->operator()( 1, 0, 1 );
                 
      answer.at( 19 ) = this->operator()( 0, 0, 2 );
      answer.at( 20 ) = this->operator()( 1, 1, 2 );
      answer.at( 21 ) = this->operator()( 2, 2, 2 );
      answer.at( 22 ) = this->operator()( 1, 2, 2 );
      answer.at( 23 ) = this->operator()( 0, 2, 2 );
      answer.at( 24 ) = this->operator()( 0, 1, 2 );
      answer.at( 25 ) = this->operator()( 2, 1, 2 );
      answer.at( 26 ) = this->operator()( 2, 0, 2 );
      answer.at( 27 ) = this->operator()( 1, 0, 2 );

      return answer;
    }

  void be_Levi_Civita()
  {
    for(int  i = 0; i < 3; i++) {
      for(int  j = 0; j < 3; j++) {
	for(int  k = 0; k < 3; k++) {
	  if ((i != j) &&(j != k) && (i != k)) {
	      if (((i == 1 ) && (j ==  2) && (k == 3 )) || ((i == 3 ) && (j ==  1) && (k == 2 )) || ((i == 2 ) && (j ==  3) && (k == 1 ))) {
		    this->operator()( i, j, k ) = 1.;
	      } else if (((i == 1 ) && (j ==  3) && (k == 2 )) || ((i == 3 ) && (j ==  2) && (k == 1 )) || ((i == 2 ) && (j ==  1) && (k == 3 ))) {
		this->operator()( i, j, k ) = -1.;
	      }
	  } else {
	    this->operator()( i, j, k ) = 0.;
	  }
	}
      }
    }

  }
  

  
};
} // namespace oofem
