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

#include "hyperelasticsurfacematerial.h"
namespace oofem {
HyperElasticSurfaceMaterial::HyperElasticSurfaceMaterial( int n, Domain *d ) :
    BaseHyperElasticMaterial() ,
    StructuralMaterial( n, d )
{}

void HyperElasticSurfaceMaterial::initializeFrom( InputRecord &ir )
{
    BaseHyperElasticMaterial::initializeFrom( ir );
    StructuralMaterial::initializeFrom( ir );
}


FloatMatrixF<9, 9>
HyperElasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    // Default implementation used if this method is not overloaded by the particular material model.
    auto status    = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    const auto &vF = status->giveTempFVector();
    const auto &vS = status->giveTempStressVector();
    auto dSdE      = this->give3dSurfaceMaterialStiffnessMatrix( mode, normal, gp, tStep );
    return convert_dSdE_2_dPdF_3D( dSdE, vS, vF );

}


FloatArrayF<9>
HyperElasticSurfaceMaterial::giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    // Default implementation used if this method is not overloaded by the particular material model.
    // 1) Compute Green-Lagrange strain and call standard method for small strains.
    // 2) Treat stress as second Piola-Kirchhoff stress and convert to first Piola-Kirchhoff stress.
    // 3) Set state variables F, P

    auto F      = from_voigt_form( vF );
    auto E      = 0.5 * ( Tdot( F, F ) - eye<3>() );
    auto vE     = to_voigt_strain( E );
    //FloatArrayF<3> vE_red( vE[0],vE[1],vE[5] ); // Reduce to 2D

    auto vS = this->giveRealSurfaceStressVector_3d( vE, normal, gp, tStep );

    // Compute first PK stress from second PK stress
    auto status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

    //FloatArrayF<6> vS{ vS_red[0], vS_red[1], .0, .0, .0, vS_red[5] }; // Enlarge to 3d

    auto S      = from_voigt_stress( vS );
    auto P      = dot( F, S );
    auto vP     = to_voigt_form( P );
    status->letTempPVectorBe( vP );
    status->letTempFVectorBe( vF );

    return vP;
}

} // namespace oofem;
