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

#include "mooneyrivlinhardmagnetic.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"
#include "domain.h"
#include "function.h"



namespace oofem {
REGISTER_Material(MooneyRivlinHardMagnetic);

FloatArrayF< 9 >
MooneyRivlinHardMagnetic::giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{

    //get stress from Mooney-Rivlin hyperelasticity
    FloatArrayF<9> vP_mr = MooneyRivlinCompressibleMaterial::giveFirstPKStressVector_3d( vF, gp, tStep );

    // calculate the magnetoelastic 1PK
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    FloatArrayF<3> B_app_at_time = load_level * B_app;

    Tensor2_3d F( vF ), P_me, delta(1., 0., 0., 0., 1., 0., 0., 0., 1. );
    
    Tensor1_3d Bapp( B_app_at_time ), Bres( B_res );
    auto [J, cofF] = F.compute_determinant_and_cofactor();


    P_me( k_3, l_3 ) = -(1/(2*J*J*mu_0))*(Bapp(i_3) - 2.0*Bres(i_3))*F(m_3,i_3)*F(m_3,j_3)*Bapp(j_3)*cofF(k_3, l_3)
        + (1/(2*J*mu_0)) * ( Bapp( i_3 ) - 2.0 * Bres( i_3 ) ) * Bapp(j_3) * (delta(i_3, l_3)*F(k_3,j_3) + delta(j_3, l_3)*F(k_3,i_3));

    auto vP_me = P_me.to_voigt_form();
    
    //add magnetoelastic and mooney-rivlin stress together
    auto vP    = vP_mr + vP_me;

    // update gp
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );
    
    return vP;
}



FloatMatrixF< 9, 9 >
MooneyRivlinHardMagnetic::give3dMaterialStiffnessMatrix_dPdF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    // get stiffness from Mooney-Rivlin hyperelasticity
    FloatMatrixF<9, 9> vD_mr = MooneyRivlinCompressibleMaterial::give3dMaterialStiffnessMatrix_dPdF( mode, gp, tStep );

    // calculate magnetoelastic stiffness
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );

    
    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    FloatArrayF<3> B_app_at_time = load_level * B_app;

    Tensor2_3d F( vF );
    Tensor4_3d D_me;

    Tensor2_3d delta(1., 0., 0., 0., 1., 0., 0., 0., 1. );

    Tensor1_3d Bapp ( B_app_at_time ), Bres( B_res );
    auto [J, cofF] = F.compute_determinant_and_cofactor();

    D_me( k_3, l_3, p_3, q_3 ) = (1/(J*J*J*mu_0)) * ( Bapp(i_3) - 2.0 * Bres(i_3))*F(m_3,i_3)*F(m_3,j_3)*Bapp(j_3)*cofF(k_3, l_3)*cofF(p_3, q_3)
        - (1/(2*J*J*mu_0)) * ( Bapp(i_3) - 2.0 * Bres(i_3)) * Bapp(j_3) * (delta(i_3, q_3)*F(p_3,j_3) + delta(j_3, q_3)*F(p_3,i_3)) * cofF(k_3, l_3) 
        - (1/(2*J*J*mu_0)) * ( Bapp(i_3) - 2.0 * Bres(i_3)) * F(m_3,i_3)*F(m_3,j_3) * Bapp(j_3) * F.compute_tensor_cross_product()(k_3, l_3, p_3, q_3)
        - (1/(2*J*J*mu_0)) * ( Bapp(i_3) - 2.0 * Bres(i_3)) * Bapp(j_3) * (delta(i_3, l_3)*F(k_3,j_3) + delta(j_3, l_3)*F(k_3,i_3)) * cofF(p_3, q_3)
        + (1/(2*J*mu_0)) * delta(k_3,p_3) * (( Bapp(l_3) - 2.0 * Bres(l_3))*Bapp(q_3) + ( Bapp(q_3) - 2.0 * Bres(q_3))*Bapp(l_3));
    
    auto vD_me = D_me.to_voigt_form();

    //add magnetoelastic and Mooney-Rivlin stiffness together
    auto vD = vD_me + vD_mr;

    return vD;
}


void MooneyRivlinHardMagnetic::initializeFrom( InputRecord &ir )
{   
    MooneyRivlinCompressibleMaterial::initializeFrom( ir );

    FloatArray B_app_temp, B_res_temp;

    IR_GIVE_FIELD( ir, B_app_temp, _IFT_MooneyRivlinHardMagnetic_B_app );
    IR_GIVE_FIELD( ir, B_res_temp, _IFT_MooneyRivlinHardMagnetic_B_res );
    IR_GIVE_FIELD( ir, mu_0, _IFT_MooneyRivlinHardMagnetic_mu_0 );
    IR_GIVE_FIELD( ir, ltf_index, _IFT_MooneyRivlinHardMagnetic_ltf );

    B_app = FloatArrayF<3>( B_app_temp );
    B_res = FloatArrayF<3>( B_res_temp );
}
} // end namespace oofem
