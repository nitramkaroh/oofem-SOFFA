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

#include "simopistermaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "domain.h"
#include "function.h"


namespace oofem {
REGISTER_Material(SimoPisterMaterial);

SimoPisterMaterial::SimoPisterMaterial(int n, Domain *d) : StructuralMaterial(n, d), BaseHyperElasticMaterial()
{ }

FloatMatrixF< 9, 9 >
SimoPisterMaterial::give3dMaterialStiffnessMatrix_dPdF(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    FloatArrayF< 9 >vF(status->giveTempFVector() );
    Tensor2_3d Fcurrent( vF ), F;

    //FloatArrayF<9> indent = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    //Tensor2_3d F0( indent );
    //FloatArray vF0_temp;
    //// Prestrain
    //if ( gp->giveElement()->giveIPValue( vF0_temp, gp, IST_PrestrainDeformationGradient, tStep ) ) {
    //    FloatArrayF<9> vF0( vF0_temp );
    //    F0 = Tensor2_3d::Tensor2_3d( vF0 );
    //}
    //F( i_3, j_3 ) = Fcurrent( i_3, k_3 ) * F0( k_3, j_3 );
    F( i_3, j_3 )  = Fcurrent( i_3, j_3 );

    auto [ J, cofF ] = F.compute_determinant_and_cofactor();
    
    Tensor4_3d A, Afull;
    //A(i_3, j_3, k_3, l_3) = G * ( 0.5 * this->compute_d2I1_C_dF2(F)(i_3, j_3, k_3, l_3) + 1 / J / J * cofF(i_3, j_3) * cofF(k_3, l_3) - 1. / J * F.compute_tensor_cross_product()(i_3, j_3, k_3, l_3) )  + this->compute_d2VolumetricEnergy_dF2(F)(i_3, j_3, k_3, l_3);
    double bulkMult = 1.;
    if ( this->bulk_ltf != 0 ) {
        bulkMult = domain->giveFunction( bulk_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    A( i_3, j_3, k_3, l_3 ) = G * ( 0.5 * this->compute_d2I1_C_dF2( F )( i_3, j_3, k_3, l_3 ) + 1 / J / J * cofF( i_3, j_3 ) * cofF( k_3, l_3 ) - 1. / J * F.compute_tensor_cross_product()( i_3, j_3, k_3, l_3 ) ) + bulkMult*this->compute_d2VolumetricEnergy_dF2( F )( i_3, j_3, k_3, l_3 );
    //double J0 = F0.compute_determinant();
    //Afull( i_3, j_3, k_3, l_3 ) = 1./J0*F0( j_3, m_3 ) * A( i_3, m_3, k_3, n_3 ) * F0( l_3, n_3 );
    Afull( i_3, j_3, k_3, l_3 ) = A( i_3, j_3, k_3, l_3 );

    return Afull.to_voigt_form();
}

FloatArrayF< 9 >
SimoPisterMaterial::giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    Tensor2_3d Fcurrent( vF ), F, P, Pfull;

    //FloatArray vF0_temp;
    //FloatArrayF<9> indent = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    //Tensor2_3d F0( indent );

    //// Prestrain
    //if ( gp->giveElement()->giveIPValue( vF0_temp, gp, IST_PrestrainDeformationGradient, tStep ) ) {
    //    FloatArrayF<9> vF0( vF0_temp );
    //    F0  = Tensor2_3d::Tensor2_3d( vF0 );
    //}
    //F( i_3, j_3 ) = Fcurrent( i_3, k_3 ) * F0( k_3, j_3 );
    F( i_3, j_3 ) = Fcurrent( i_3, j_3 );

    auto [ J, cofF ] = F.compute_determinant_and_cofactor();

    // compute the first Piola-Kirchhoff
    //P(i_3, j_3) =  G * ( 0.5 * this->compute_dI1_C_dF(F)(i_3, j_3) - 1. / J * cofF(i_3, j_3) ) + this->compute_dVolumetricEnergy_dF(F)(i_3, j_3);
    double bulkMult = 1.;
    if ( this->bulk_ltf != 0 ) {
        bulkMult = domain->giveFunction( bulk_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    P( i_3, j_3 ) = G * ( 0.5 * this->compute_dI1_C_dF( F )( i_3, j_3 ) - 1. / J * cofF( i_3, j_3 ) ) + bulkMult*this->compute_dVolumetricEnergy_dF( F )( i_3, j_3 );

    //double J0 = F0.compute_determinant();
    //Pfull( i_3, j_3 ) = 1. / J0 * P( i_3, k_3 ) * F0( j_3, k_3 ); 
    Pfull( i_3, j_3 ) = P( i_3, j_3 ); 

    auto vP = Pfull.to_voigt_form();
    // update gp
    status->letTempFVectorBe(vF);
    //status->letTempFVectorBe( F.to_voigt_form() );
    status->letTempPVectorBe(vP);

    return vP;
}


MaterialStatus *
SimoPisterMaterial::CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(gp);
}


void
SimoPisterMaterial::initializeFrom(InputRecord &ir)
{
    StructuralMaterial::initializeFrom(ir);
    BaseHyperElasticMaterial::initializeFrom(ir);
    IR_GIVE_FIELD(ir, G, _IFT_SimoPisterMaterial_g);
    IR_GIVE_OPTIONAL_FIELD( ir, bulk_ltf, _IFT_SimoPisterMaterial_bulkLTF );
}
} // end namespace oofem
