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
    //Tensor1_3d e1( { 1, 0, 0 } );
    //Tensor1_3d e2( { 0, 1, 0 } );

    //Tensor1_3d F1, F2;
    //F1( i_3 ) = F( i_3, j_3 ) * e1( j_3 );
    //F2( i_3 ) = F( i_3, j_3 ) * e2( j_3 );
    return 0.;
};

Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_cofactor( const Tensor2_3d &F ) const
{
    Tensor2_3d cofF;

    return cofF;
};

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_dCof_dF( const Tensor2_3d &F ) const
{
    Tensor4_3d A;

    return A;
};



} // end namespace oofem
