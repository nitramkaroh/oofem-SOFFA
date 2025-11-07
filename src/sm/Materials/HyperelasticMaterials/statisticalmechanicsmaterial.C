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

#include "statisticalmechanicsmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "domain.h"
#include "function.h"


namespace oofem {
REGISTER_Material( StatisticalMechanicsMaterial );

StatisticalMechanicsMaterial::StatisticalMechanicsMaterial( int n, Domain *d ) :
    StructuralMaterial( n, d ), BaseHyperElasticMaterial(), LargeStrainMixedPressureMaterialExtensionInterface( d )
{
}

FloatMatrixF<9, 9>
StatisticalMechanicsMaterial::give3dMaterialStiffnessMatrix_dPdF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    Tensor2_3d F( vF );

    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d I( indent3d );

    //auto [J, cofF] = F.compute_determinant_and_cofactor();

    double rho, beta, dbeta;
    Tensor2_3d integral1;
    Tensor4_3d integral2;

    FloatArrayF<3> zero = { 0., 0., 0. };
    Tensor1_3d test( zero );
    double t = 0., t2 = 0.;
    for ( int i = 0; i < nPoints; i++ ) {
        rho = this->rho0* sqrt( F( i_3, j_3 ) * this->integrationPositions[i]( j_3 ) * F( i_3, k_3 ) * this->integrationPositions[i]( k_3 ) );

        beta                 = rho * ( 3. - rho * rho ) / (1. - rho * rho );
        dbeta                = beta * beta * sinh( beta ) * sinh( beta ) / ( sinh( beta ) * sinh( beta ) - beta * beta );
        integral1( i_3, j_3 ) = integral1( i_3, j_3 ) + beta / rho * this->RoR[i]( i_3, j_3 ) * integrationWeights[i];
        integral2( i_3, j_3, k_3, l_3 ) = integral2( i_3, j_3, k_3, l_3 ) + ( dbeta * rho - beta ) / rho / rho / rho * this->RoRoRoR[i]( i_3, j_3, k_3, l_3 ) * integrationWeights[i];

        //////
        //test( i_3 ) = test( i_3 ) + this->integrationPositions[i]( i_3 ) * integrationWeights[i];
    };
    Tensor4_3d A1, A2, A3, A;

    //A1( i_3, k_3, m_3, n_3 ) = N0 * k * T / 4 / M_PI * 2 / 3 * pow( J, -1 / 3 ) * F( i_3, j_3 ) * integral1( j_3, k_3 ) * cofF( m_3, n_3 );
    //A2( i_3, k_3, o_3, n_3 ) = G / 4 / M_PI * pow( J, 4 / 3 ) * integral2( j_3, k_3, m_3, n_3 ) * F( i_3, j_3 ) * F( o_3, m_3 );
    //A3( i_3, k_3, j_3, l_3 ) = G / 4 / M_PI * pow( J, 2 / 3 ) * I( i_3, j_3 ) * integral1( k_3, l_3 );
    A2( i_3, k_3, o_3, n_3 ) = G / 4. / M_PI * integral2( j_3, k_3, m_3, n_3 ) * F( i_3, j_3 ) * F( o_3, m_3 );
    A3( i_3, k_3, j_3, l_3 ) = G / 4. / M_PI * I( i_3, j_3 ) * integral1( k_3, l_3 );

    A( i_3, k_3, j_3, l_3 ) = A2( i_3, k_3, j_3, l_3 ) + A3( i_3, k_3, j_3, l_3 );

    return A.to_voigt_form();
}

FloatArrayF<9>
StatisticalMechanicsMaterial::giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first Piola-Kirchhoff stress corresponding to the given deformation gradient
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    Tensor2_3d F( vF ), P;

    //auto [J, cofF] = F.compute_determinant_and_cofactor();

    double rho,beta;
    Tensor2_3d integral;
    for ( int i = 0; i < nPoints; i++ ) {
        rho = this->rho0 * sqrt( F( i_3, j_3 ) * this->integrationPositions[i]( j_3 ) * F( i_3, k_3 ) * this->integrationPositions[i]( k_3 ) );

        beta                 = rho * ( 3. - rho * rho ) / ( 1. - rho * rho );
        integral( i_3, j_3 ) = integral( i_3, j_3 ) + beta/rho* this->RoR[i]( i_3, j_3 )*integrationWeights[i];
    };

    //P( i_3, k_3 ) = N0 * k * T * pow( J, 2 / 3 ) / 4 / M_PI * F( i_3, j_3 ) * integral( j_3, k_3 );
    P( i_3, k_3 ) = G / 4 / M_PI * F( i_3, j_3 ) * integral( j_3, k_3 );


    auto vP = P.to_voigt_form();
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );

    return vP;
}

MaterialStatus *
StatisticalMechanicsMaterial::CreateStatus( GaussPoint *gp ) const
{
    return new StructuralMaterialStatus( gp );
}

void StatisticalMechanicsMaterial::initializeFrom( InputRecord &ir )
{
    StructuralMaterial::initializeFrom( ir );
    BaseHyperElasticMaterial::initializeFrom( ir );
    //IR_GIVE_OPTIONAL_FIELD( ir, N0, _IFT_StatisticalMechanicsMaterial_N0 );
    //IR_GIVE_OPTIONAL_FIELD( ir, l, _IFT_StatisticalMechanicsMaterial_l );
    //IR_GIVE_OPTIONAL_FIELD( ir, n, _IFT_StatisticalMechanicsMaterial_n );
    //IR_GIVE_OPTIONAL_FIELD( ir, T, _IFT_StatisticalMechanicsMaterial_T );
    //
    //IR_GIVE_OPTIONAL_FIELD( ir, R, _IFT_StatisticalMechanicsMaterial_R );
    IR_GIVE_OPTIONAL_FIELD( ir, rho0, _IFT_StatisticalMechanicsMaterial_rho0 );
    IR_GIVE_OPTIONAL_FIELD( ir, G, _IFT_StatisticalMechanicsMaterial_G );
    IR_GIVE_OPTIONAL_FIELD( ir, nPoints, _IFT_StatisticalMechanicsMaterial_nPoints );
    this->initializeData( nPoints);
    this->constructTensorProducts();
}

void StatisticalMechanicsMaterial ::initializeData( int numberOfPoints )
{

    //// definition of intergation points and weights
    // std::vector<FloatArrayF<3> > positionsHelp;
    // positionsHelp.resize( numberOfPoints );

    integrationPositions.resize( numberOfPoints );
    integrationWeights.resize( numberOfPoints );

    if ( numberOfPoints == 42 ) {
        for ( int i = 0; i < 3; i++ ) {
            integrationWeights[i] = 0.02652141274;
        }

        for ( int i = 3; i < 9; i++ ) {
            integrationWeights[i] = 0.01993014153;
        }

        for ( int i = 9; i < 21; i++ ) {
            integrationWeights[i] = 0.02507124272;
        }


        integrationPositions[0]( 0 ) = integrationPositions[1]( 1 ) = integrationPositions[2]( 2 ) = 1.;
        integrationPositions[0]( 1 ) = integrationPositions[0]( 2 ) = integrationPositions[1]( 0 ) = integrationPositions[1]( 2 ) = integrationPositions[2]( 0 ) = integrationPositions[2]( 1 ) = 0.;
        integrationPositions[3]( 0 ) = integrationPositions[3]( 1 ) = integrationPositions[4]( 0 ) = integrationPositions[5]( 0 ) = integrationPositions[5]( 2 ) = integrationPositions[6]( 0 ) = integrationPositions[7]( 1 ) = integrationPositions[7]( 2 ) = integrationPositions[8]( 1 ) = 0.7071067812;
        integrationPositions[4]( 1 ) = integrationPositions[6]( 2 ) = integrationPositions[8]( 2 ) = -0.7071067812;
        integrationPositions[3]( 2 ) = integrationPositions[4]( 2 ) = integrationPositions[5]( 1 ) = integrationPositions[6]( 1 ) = integrationPositions[7]( 0 ) = integrationPositions[8]( 0 ) = 0.;


        double help[3][3] = { { 0.3879072746, 0.3879072746, 0.8360956240 },
            { 0.3879072746, 0.8360956240, 0.3879072746 },
            { 0.8360956240, 0.3879072746, 0.3879072746 } };

        for ( int i = 0; i < 3; i++ ) {
            int i4                            = i * 4;
            integrationPositions[9 + i4]( 0 ) = help[i][0];
            integrationPositions[9 + i4]( 1 ) = help[i][1];
            integrationPositions[9 + i4]( 2 ) = help[i][2];

            integrationPositions[10 + i4]( 0 ) = help[i][0];
            integrationPositions[10 + i4]( 1 ) = help[i][1];
            integrationPositions[10 + i4]( 2 ) = -help[i][2];

            integrationPositions[11 + i4]( 0 ) = help[i][0];
            integrationPositions[11 + i4]( 1 ) = -help[i][1];
            integrationPositions[11 + i4]( 2 ) = help[i][2];

            integrationPositions[12 + i4]( 0 ) = help[i][0];
            integrationPositions[12 + i4]( 1 ) = -help[i][1];
            integrationPositions[12 + i4]( 2 ) = -help[i][2];
        }
        
        // For the other halfsphere mirror the values
        for ( int i = 0; i < numberOfPoints/2; i++ ) {
            integrationWeights[i + numberOfPoints / 2] = integrationWeights[i];

            integrationPositions[i + numberOfPoints / 2]( 0 ) = -integrationPositions[i]( 0 );
            integrationPositions[i + numberOfPoints / 2]( 1 ) = -integrationPositions[i]( 1 );
            integrationPositions[i + numberOfPoints / 2]( 2 ) = -integrationPositions[i]( 2 );
        }



    } else {
        OOFEM_ERROR( "Unsupported number of microplanes" );
    }
}


void StatisticalMechanicsMaterial::constructTensorProducts()
{
    this->RoR.resize( this->nPoints );
    this->RoRoRoR.resize( this->nPoints );

    for ( int i = 0; i < this->nPoints; i++ ) {
        this->RoR[i]( i_3, j_3 ) = this->integrationPositions[i]( i_3 ) * this->integrationPositions[i]( j_3 );
        this->RoRoRoR[i]( i_3, j_3, k_3, l_3 ) = this->integrationPositions[i]( i_3 ) * this->integrationPositions[i]( j_3 ) * 
            this->integrationPositions[i]( k_3 ) * this->integrationPositions[i]( l_3 );
    }
}


/////////////////////////// Functions for Mixed Pressure formulation

std::tuple<FloatArrayF<9>, double>
StatisticalMechanicsMaterial ::giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, double pressure, GaussPoint *gp, TimeStep *tStep )
{

    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    //
    Tensor2_3d F( vF ), P;
    //
    auto [detF, cofF] = F.compute_determinant_and_cofactor();

    double rho, beta;
    Tensor2_3d integral;
    for ( int i = 0; i < nPoints; i++ ) {
        rho = this->rho0 * sqrt( F( i_3, j_3 ) * this->integrationPositions[i]( j_3 ) * F( i_3, k_3 ) * this->integrationPositions[i]( k_3 ) );
        beta                 = rho * ( 3. - rho * rho ) / ( 1. - rho * rho );
        integral( i_3, j_3 ) = integral( i_3, j_3 ) + beta / rho * this->RoR[i]( i_3, j_3 ) * integrationWeights[i];
    };
    P( i_3, k_3 ) = G / 4 / M_PI * F( i_3, j_3 ) * integral( j_3, k_3 ) + pressure * cofF( i_3, k_3 );
    auto vP       = P.to_voigt_form();

    // compute cauchy stress
    Tensor2_3d sigma;
    sigma( i_3, j_3 ) = 1 / detF * P( i_3, k_3 ) * F( j_3, k_3 );
    auto sigmaP       = sigma.to_voigt_form();


    // update gp
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );
    status->letCVectorBe( sigmaP );
    //
    auto dJ = ( detF - 1. ) - pressure / this->K;
    return std::make_tuple( vP, dJ );
}

std::tuple<FloatArrayF<5>, double>
StatisticalMechanicsMaterial ::giveFirstPKStressVector_PlaneStrain( const FloatArrayF<5> &vF, double pressure, GaussPoint *gp, TimeStep *tStep )
{
    auto [vP, dJ] = this->giveFirstPKStressVector_3d( assemble<9>( vF, { 0, 1, 2, 5, 8 } ), pressure, gp, tStep );
    //
    return std::make_tuple( vP[{ 0, 1, 2, 5, 8 }], dJ );
}


std::tuple<FloatMatrixF<9, 9>, FloatArrayF<9>, double>
StatisticalMechanicsMaterial ::giveLargeStrainMixedPressureConstitutiveMatrices_3d( double pressure, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    Tensor2_3d F( vF );

    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d I( indent3d );

    double rho, beta, dbeta;
    Tensor2_3d integral1;
    Tensor4_3d integral2;

    for ( int i = 0; i < nPoints; i++ ) {
        rho = this->rho0 * sqrt( F( i_3, j_3 ) * this->integrationPositions[i]( j_3 ) * F( i_3, k_3 ) * this->integrationPositions[i]( k_3 ) );
        beta                            = rho * ( 3. - rho * rho ) / ( 1. - rho * rho );
        dbeta                           = beta * beta * sinh( beta ) * sinh( beta ) / ( sinh( beta ) * sinh( beta ) - beta * beta );
        integral1( i_3, j_3 )           = integral1( i_3, j_3 ) + beta / rho * this->RoR[i]( i_3, j_3 ) * integrationWeights[i];
        integral2( i_3, j_3, k_3, l_3 ) = integral2( i_3, j_3, k_3, l_3 ) + ( dbeta * rho - beta ) / rho / rho / rho * this->RoRoRoR[i]( i_3, j_3, k_3, l_3 ) * integrationWeights[i];
    };
    Tensor4_3d A2, A3, A;
    A2( i_3, k_3, o_3, n_3 ) = G / 4. / M_PI * integral2( j_3, k_3, m_3, n_3 ) * F( i_3, j_3 ) * F( o_3, m_3 );
    A3( i_3, k_3, j_3, l_3 ) = G / 4. / M_PI * I( i_3, j_3 ) * integral1( k_3, l_3 );
    A( i_3, k_3, j_3, l_3 )  = A2( i_3, k_3, j_3, l_3 ) + A3( i_3, k_3, j_3, l_3 ) + pressure * this->compute_d2J_dF2( F )( i_3, k_3, j_3, l_3 );

    //
    auto [detF, cofF] = F.compute_determinant_and_cofactor();
    //
    return std::make_tuple( A.to_voigt_form(), cofF.to_voigt_form(), -this->K );
}


std::tuple<FloatMatrixF<5, 5>, FloatArrayF<5>, double>
StatisticalMechanicsMaterial ::giveLargeStrainMixedPressureConstitutiveMatrices_PlaneStrain( double pressure, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    auto [A, cofF, minusK] = this->giveLargeStrainMixedPressureConstitutiveMatrices_3d( pressure, rMode, gp, tStep );
    //
    return std::make_tuple( A( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } ), cofF[{ 0, 1, 2, 5, 8 }], minusK );
}


} // end namespace oofem
