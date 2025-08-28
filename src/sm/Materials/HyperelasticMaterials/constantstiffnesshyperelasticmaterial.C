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

#include "constantstiffnesshyperelasticmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"



namespace oofem {
REGISTER_Material( ConstantStiffnessHyperElasticMaterial );

ConstantStiffnessHyperElasticMaterial::ConstantStiffnessHyperElasticMaterial( int n, Domain *d ) :
    StructuralMaterial( n, d ), BaseHyperElasticMaterial()
{ }

FloatArrayF< 9 >
ConstantStiffnessHyperElasticMaterial::giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    Tensor2_3d F(vF), P;

    Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
    auto [J, cofF] = F.compute_determinant_and_cofactor();

    // compute the first Piola-Kirchhoff
    P( i_3, j_3 ) = this->compute_dVolumetricEnergy_dF( F )( i_3, j_3 )
      + lambda * delta( i_3, m_3 ) * ( F( k_3, k_3 ) - 3. ) * cofF( m_3, j_3 ) + mu * ( F( i_3, m_3 ) + F( m_3, i_3 ) - 2. * delta( i_3, m_3 ) ) * cofF( m_3, j_3 );

      
    auto vP = P.to_voigt_form();
    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(vP);
    //
    return vP;
}



FloatMatrixF< 9, 9 >
ConstantStiffnessHyperElasticMaterial::give3dMaterialStiffnessMatrix_dPdF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    FloatArrayF< 9 >vF(status->giveTempFVector() );
    Tensor2_3d F(vF);
    Tensor4_3d A;

    Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
    auto [J, cofF] = F.compute_determinant_and_cofactor();
    auto Fcross = F.compute_tensor_cross_product();

    A( i_3, j_3, k_3, l_3 ) = this->compute_d2VolumetricEnergy_dF2( F )( i_3, j_3, k_3, l_3 )
      + lambda * delta( i_3, m_3 ) * delta( k_3, l_3 ) * cofF( m_3, j_3 ) + mu * ( delta( i_3, k_3 ) * delta( m_3, l_3 ) + delta( i_3, l_3 ) * delta( m_3, k_3 ) ) * cofF( m_3, j_3 )
      + lambda * delta( i_3, m_3 ) * ( F( k_3, k_3 ) - 3. ) * Fcross( m_3, j_3, k_3, l_3 ) + mu * ( F( i_3, m_3 ) + F( m_3, i_3 ) - 2. * delta( i_3, m_3 ) ) * Fcross( m_3, j_3, k_3, l_3 );
    
    return A.to_voigt_form();
}




MaterialStatus *
ConstantStiffnessHyperElasticMaterial::CreateStatus( GaussPoint *gp ) const
{
    return new StructuralMaterialStatus(gp);
}


void ConstantStiffnessHyperElasticMaterial::initializeFrom( InputRecord &ir )
{
    StructuralMaterial::initializeFrom(ir);
    BaseHyperElasticMaterial::initializeFrom(ir);
    double E, nu;
    IR_GIVE_FIELD( ir, E, _IFT_ConstantStiffnessHyperElasticMaterial_E );
    IR_GIVE_FIELD( ir, nu, _IFT_ConstantStiffnessHyperElasticMaterial_nu );

    this->lambda = E * nu / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
    this->mu = E / ( 2. * ( 1. + nu ) );
}
} // end namespace oofem
