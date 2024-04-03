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
 *               Copyright (C) 1993 - 2020   Borek Patzak
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

#include "basehyperelasticsurfacematerial.h"
#include "dynamicinputrecord.h"

namespace oofem {

double BaseHyperElasticSurfaceMaterial::compute_surface_determinant( const Tensor2_3d &F ) const {

    double p1 = std::pow( F( 0, 0 ) * F( 1, 1 ) - F( 0, 1 ) * F( 1, 0 ), 2 );
    double p2 = std::pow( F( 0, 0 ) * F( 2, 1 ) - F( 0, 1 ) * F( 2, 0 ), 2 );
    double p3 = std::pow( F( 1, 0 ) * F( 2, 1 ) - F( 1, 1 ) * F( 2, 0 ), 2 );

    return std::sqrt( p1 + p2 + p3 );
};

Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_cofactor( const Tensor2_3d &F ) const
{
    Tensor2_3d cofF;
    double F11, F12, F21, F22, F31, F32;
    F11 = F( 0, 0 );F12 = F( 0, 1 );F21 = F( 1, 0 );
    F22 = F( 1, 1 );F31 = F( 2, 0 );F32 = F( 2, 1 );

    cofF( 0, 0 ) = F11 * ( pow(F12, 2) + pow(F22 , 2) + pow(F32 , 2) ) - F12 * ( F11 * F12 + F21 * F22 + F31 * F32 );
    cofF( 0, 1 ) = F12 * ( pow(F11, 2) + pow(F21 , 2) + pow(F31 , 2) ) - F11 * ( F11 * F12 + F21 * F22 + F31 * F32 );
    cofF( 1, 0 ) = F21 * ( pow(F12, 2) + pow(F22 , 2) + pow(F32 , 2) ) - F22 * ( F11 * F12 + F21 * F22 + F31 * F32 );
    cofF( 1, 1 ) = F22 * ( pow(F11, 2) + pow(F21 , 2) + pow(F31 , 2) ) - F21 * ( F11 * F12 + F21 * F22 + F31 * F32 );
    cofF( 2, 0 ) = F31 * ( pow(F12, 2) + pow(F22 , 2) + pow(F32 , 2) ) - F32 * ( F11 * F12 + F21 * F22 + F31 * F32 );
    cofF( 2, 1 ) = F32 * ( pow(F11, 2) + pow(F21 , 2) + pow(F31 , 2) ) - F31 * ( F11 * F12 + F21 * F22 + F31 * F32 );

    cofF( i_3, j_3 ) = 1/this->compute_surface_determinant(F)*cofF( i_3, j_3 );

    return cofF;
};

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_dCof_dF( const Tensor2_3d &F ) const
{
    Tensor4_3d A1,A2,A;
    double F11, F12, F21, F22, F31, F32;
    F11 = F( 0, 0 ); F12 = F( 0, 1 );F21 = F( 1, 0 );
    F22 = F( 1, 1 ); F31 = F( 2, 0 );F32 = F( 2, 1 );

    A1( 0, 0, 0, 0 ) = pow( F22, 2 ) + pow( F32, 2 );
    A1( 0, 1, 0, 0 ) = -F21 * F22 - F31 * F32;
    A1( 1, 0, 0, 0 ) = -F12 * F22;
    A1( 1, 1, 0, 0 ) = 2 * F11 * F22 - F12 * F21;
    A1( 2, 0, 0, 0 ) = -F12 * F32;
    A1( 2, 1, 0, 0 ) = 2 * F11 * F32 - F12 * F31;

    A1( 0, 0, 1, 0 ) = -F12 * F22;
    A1( 0, 1, 1, 0 ) = 2 * F12 * F21 - F11 * F22;
    A1( 1, 0, 1, 0 ) = pow(F12 , 2) + pow(F32 , 2);
    A1( 1, 1, 1, 0 ) = -F11 * F12 - F31 * F32;
    A1( 2, 0, 1, 0 ) = -F22 * F32;
    A1( 2, 1, 1, 0 ) = 2 * F21 * F32 - F22 * F31;

    A1( 0, 0, 2, 0 ) = -F12 * F32;
    A1( 0, 1, 2, 0 ) = 2 * F12 * F31 - F11 * F32;
    A1( 1, 0, 2, 0 ) = -F22 * F32;
    A1( 1, 1, 2, 0 ) = 2 * F22 * F31 - F21 * F32;
    A1( 2, 0, 2, 0 ) = pow(F12 , 2) + pow(F22 , 2);
    A1( 2, 1, 2, 0 ) = -F11 * F12 - F21 * F22;

    A1( 0, 0, 0, 1 ) = -F21 * F22 - F31 * F32;
    A1( 0, 1, 0, 1 ) = pow( F21, 2 ) + pow( F31, 2 );
    A1( 1, 0, 0, 1 ) = 2 * F12 * F21 - F11 * F22;
    A1( 1, 1, 0, 1 ) = -F11 * F21;
    A1( 2, 0, 0, 1 ) = 2 * F12 * F31 - F11 * F32;
    A1( 2, 1, 0, 1 ) = -F11 * F31;

    A1( 0, 0, 1, 1 ) = 2 * F11 * F22 - F12 * F21;
    A1( 0, 1, 1, 1 ) = -F11 * F21;
    A1( 1, 0, 1, 1 ) = -F11 * F12 - F31 * F32;
    A1( 1, 1, 1, 1 ) = pow( F11, 2 ) + pow( F31, 2 );
    A1( 2, 0, 1, 1 ) = 2 * F22 * F31 - F21 * F32;
    A1( 2, 1, 1, 1 ) = -F21 * F31;

    A1( 0, 0, 2, 1 ) = 2 * F11 * F32 - F12 * F31;
    A1( 0, 1, 2, 1 ) = -F11 * F31;
    A1( 1, 0, 2, 1 ) = 2 * F21 * F32 - F22 * F31;
    A1( 1, 1, 2, 1 ) = -F21 * F31;
    A1( 2, 0, 2, 1 ) = -F11 * F12 - F21 * F22;
    A1( 2, 1, 2, 1 ) = pow( F11, 2 ) + pow( F21, 2 );

    Tensor2_3d cofF  = this->compute_surface_cofactor( F );
    A2( i_3, j_3, k_3, l_3 ) = cofF( i_3, j_3 ) * cofF( k_3, l_3 );

    A( i_3, j_3, k_3, l_3 ) = 1 / this->compute_surface_determinant( F ) * (A1( i_3, j_3, k_3, l_3 ) - A2( i_3, j_3, k_3, l_3 ));

    return A;
};


Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_d_normF_dF( const Tensor2_3d &F ) const
{ 
    Tensor2_3d answer;
    answer( i_3, j_3 ) = 1 / this->compute_surface_NormF( F ) * F( i_3, j_3 );
    return answer;
};

double BaseHyperElasticSurfaceMaterial::compute_surface_NormF( const Tensor2_3d &F ) const
{
    return sqrt( F( i_3, j_3 ) * F( i_3, j_3 ) );
}

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_d2_normF_dF2( const Tensor2_3d &F ) const
{
    Tensor4_3d answer, part1, part2;
    double normF = this->compute_surface_NormF( F );

    part1( 0, 0, 0, 0 ) = 1;part1( 1, 0, 1, 0 ) = 1;part1( 2, 0, 2, 0 ) = 1;
    part1( 0, 1, 0, 1 ) = 1;part1( 1, 1, 1, 1 ) = 1;part1( 2, 1, 2, 1 ) = 1;
    part1( 0, 2, 0, 2 ) = 1;part1( 1, 2, 1, 2 ) = 1;part1( 2, 2, 2, 2 ) = 1;

    part1( i_3, j_3, k_3, l_3 ) = 1 / normF * part1( i_3, j_3, k_3, l_3 );

    part2( i_3, j_3, k_3, l_3 ) = -1 / pow( normF, 3 ) * F( i_3, j_3 ) * F( k_3, l_3 );

    answer( i_3, j_3, k_3, l_3 ) = part1( i_3, j_3, k_3, l_3 ) + part2( i_3, j_3, k_3, l_3 );
    return answer;
};

} // end namespace oofem
