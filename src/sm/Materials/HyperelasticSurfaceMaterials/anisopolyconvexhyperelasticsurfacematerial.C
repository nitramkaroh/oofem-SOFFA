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

#include "anisopolyconvexhyperelasticsurfacematerial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"
#include "domain.h"
#include "function.h"
#include "sm/Elements/nlstructuralsurfaceelement.h"
#include "sm/Elements/SurfaceElements/quadsurface.h"
#include "structuralelement.h"
namespace oofem {
REGISTER_Material( AnisotropicPolyconvexHyperelasticSurfaceMaterial );

AnisotropicPolyconvexHyperelasticSurfaceMaterial::AnisotropicPolyconvexHyperelasticSurfaceMaterial( int n, Domain *d ) :
    IsotropicPolyconvexHyperelasticSurfaceMaterial( n, d )
    //beta(0.0)
{
}

FloatArrayF<9>
AnisotropicPolyconvexHyperelasticSurfaceMaterial::giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    //// Isotropic stress
    double alpha_t,eta_t;
    if ( this->gamma_ltf == 0 ) {
        //alpha_t = this->alpha1;
        eta_t = this->eta;
    } else {
        //alpha_t = this->alpha1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        eta_t   = this->eta * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }

    FloatArrayF<9> PisoMat = IsotropicPolyconvexHyperelasticSurfaceMaterial::giveFirstPKSurfaceStressVector_3d( vF, normal, gp, tStep );
    Tensor2_3d Piso( PisoMat ), F( vF ), P, S;

    //if ( abs( beta_t ) < 1e-6 ) {
    //    return PisoMat;
    //}
    // 
    //// Anisotropic part
    FloatMatrix matdirLoc = computeLocalMaterialDirection( gp );
    FloatMatrix Fmat, Fdot_a;
    F.toFloatMatrix( Fmat );
    Fdot_a.beProductOf( Fmat, matdirLoc );
    double normFa = sqrt( Fdot_a.at( 1, 1 ) * Fdot_a.at( 1, 1 ) + Fdot_a.at( 2, 1 ) * Fdot_a.at( 2, 1 ) + Fdot_a.at( 3, 1 ) * Fdot_a.at( 3, 1 ) );
    
    if ( normFa < 1e-6 )
    {
        return PisoMat;
    }

    Tensor2_3d Paniso, Paniso2;
    Tensor1_3d a,Fatens;
    a( 0 ) = matdirLoc.at( 1, 1 );
    a( 1 ) = matdirLoc.at( 2, 1 );
    a( 2 ) = matdirLoc.at( 3, 1 );

    Fatens( 0 ) = Fdot_a.at( 1, 1 );
    Fatens( 1 ) = Fdot_a.at( 2, 1 );
    Fatens( 2 ) = Fdot_a.at( 3, 1 );

    double a0 = a( 0 );
    double a1 = a( 1 );

    //Paniso( i_3, j_3 ) = alpha_t / normFa * F( i_3, k_3 ) * a( k_3 ) * a( j_3 );

    Paniso( i_3, j_3 ) = 1 / normFa * F( i_3, k_3 ) * a( k_3 ) * a( j_3 );
    //Paniso2( i_3, j_3 ) = 1 / normFa * Fatens( i_3 ) * a( j_3 );
    //Paniso.printYourself();
    //Paniso2.printYourself();


    // Sum isotropic and anisotropic part
    //P( i_3, j_3 ) = Piso( i_3, j_3 ) + this->eta* Paniso( i_3, j_3 );
    P( i_3, j_3 ) = Piso( i_3, j_3 ) + eta_t * Paniso( i_3, j_3 );


    //Tensor2_3d Paniso, a_dot_d_normal_dF;
    //FloatArray DefNormal  = this->compute_deformed_normal( F ); 
    //FloatArray matdirLoc = computeLocalMaterialDirection( gp );
    //double mult = DefNormal.dotProduct( matdirLoc );

    //// compute sign of the absolute valua of a dot n
    //double sign;
    //if ( abs( mult ) < 1e-6 ) {
    //    sign = 0.;
    //} else if ( mult > 0. ) {
    //    sign = 1.;
    //} else {
    //    sign = -1.;
    //}

    //a_dot_d_normal_dF( i_3, j_3 ) = this->compute_surface_a_dot_d_normal_dF( F, matdirLoc )( i_3, j_3 );

    //double energyIso = this->computeEnergy( F, gp, tStep );

    //Paniso( i_3, j_3 ) = beta_t * ( Piso( i_3, j_3 ) * abs( mult ) + sign * energyIso * a_dot_d_normal_dF( i_3, j_3 ) );

    //// Sum isotropic and anisotropic part
    //P( i_3, j_3 ) = Piso( i_3, j_3 ) + Paniso( i_3, j_3 );

    // compute Cauchy stress vector
    S( i_3, j_3 ) = ( 1 / this->compute_surface_determinant( F ) ) * P( i_3, k_3 ) * F( j_3, k_3 );

    auto vP = P.to_voigt_form();
    auto vS = S.to_voigt_form();

    // update gp
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );
    status->letTempCVectorBe( vS );

    return vP;
}


FloatMatrixF<9, 9>
AnisotropicPolyconvexHyperelasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );

    //// Isotropic part
    double alpha_t, eta_t;
    if ( this->gamma_ltf == 0 ) {
        //alpha_t = this->alpha1;
        eta_t   = this->eta;
    } else {
        //alpha_t = this->alpha1 * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        eta_t   = this->eta * domain->giveFunction( gamma_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    
    // Isotropic stiffness matrix
    FloatMatrixF<9, 9> AisoMat = IsotropicPolyconvexHyperelasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dPdF( mode, normal, gp, tStep );
    Tensor4_3d Aiso( AisoMat ), A, Aaniso;
    //if ( abs( beta_t ) < 1e-6 ) {
    //    return AisoMat;
    //}

    // Isotropic stress
    //FloatArrayF<9> PisoMat = IsotropicPolyconvexHyperelasticSurfaceMaterial::giveFirstPKSurfaceStressVector_3d( vF, normal, gp, tStep );
    //Tensor2_3d Piso( PisoMat ), F( vF );

    ////// Anisotropic part
    Tensor2_3d F( vF );
    FloatMatrix matdirLoc = computeLocalMaterialDirection( gp );
    FloatMatrix Fmat, Fdot_a;
    F.toFloatMatrix( Fmat );
    Fdot_a.beProductOf( Fmat, matdirLoc );
    double normFa = sqrt( Fdot_a.at( 1, 1 ) * Fdot_a.at( 1, 1 ) + Fdot_a.at( 2, 1 ) * Fdot_a.at( 2, 1 ) + Fdot_a.at( 3, 1 ) * Fdot_a.at( 3, 1 ) );
    
    if ( normFa < 1e-6 ) {
        return AisoMat;
    }

    Tensor1_3d a;
    a( 0 ) = matdirLoc.at( 1, 1 );
    a( 1 ) = matdirLoc.at( 2, 1 );
    a( 2 ) = matdirLoc.at( 3, 1 );

    double a0 = a( 0 );
    double a1 = a( 1 );


    FloatArrayF<9> indent = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Identity( indent ), aoa;
    aoa( i_3, j_3 )  = a( i_3 ) * a( j_3 );
    /*Aaniso( i_3, j_3, k_3, l_3 ) = alpha_t*( 1 / normFa * Identity( i_3, k_3 ) * aoa( j_3, l_3 ) 
                                    - 1 / normFa / normFa / normFa * F( i_3, m_3 ) * aoa( m_3, j_3 ) * F( k_3, n_3 ) * aoa( n_3, l_3 ) );*/
    Aaniso( i_3, j_3, k_3, l_3 ) = ( 1 / normFa * Identity( i_3, k_3 ) * aoa( j_3, l_3 ) - 
        1 / normFa / normFa / normFa * F( i_3, m_3 ) * aoa( m_3, j_3 ) * F( k_3, n_3 ) * aoa( n_3, l_3 ) );
    

    //A( i_3, j_3, k_3, l_3 ) = Aiso( i_3, j_3, k_3, l_3 ) + this->eta*Aaniso( i_3, j_3, k_3, l_3 );
    A( i_3, j_3, k_3, l_3 ) = Aiso( i_3, j_3, k_3, l_3 ) + eta_t * Aaniso( i_3, j_3, k_3, l_3 );



    ////// Anisotropic part
    //Tensor2_3d a_dot_d_normal_dF;
    //Tensor4_3d a_dot_d2_normal_d2F;

    //FloatArray DefNormal = this->compute_deformed_normal( F );
    //FloatArray matdirLoc = computeLocalMaterialDirection( gp );
    //double mult          = DefNormal.dotProduct( matdirLoc );

    //// compute sign of the absolute valua of a dot n
    //double sign;
    //if ( abs( mult ) < 1e-6 ) {
    //    sign = 0.;
    //} else if ( mult > 0. ) {
    //    sign = 1.;
    //} else {
    //    sign = -1.;
    //}

    //a_dot_d_normal_dF( i_3, j_3 ) = this->compute_surface_a_dot_d_normal_dF( F, matdirLoc )( i_3, j_3 );

    //a_dot_d2_normal_d2F( i_3, j_3, k_3, l_3 ) = this->compute_surface_a_dot_d2_normal_d2F( F, matdirLoc )( i_3, j_3, k_3, l_3 );

    //double energyIso = this->computeEnergy( F, gp, tStep );

    //Aaniso( i_3, j_3, k_3, l_3 ) = beta_t * ( Aiso( i_3, j_3, k_3, l_3 ) * abs( mult ) + sign * ( sign* Piso( i_3, j_3 ) * a_dot_d_normal_dF( k_3, l_3 ) + 
    //                                        Piso( k_3, l_3 ) * a_dot_d_normal_dF( i_3, j_3 ) + energyIso * a_dot_d2_normal_d2F( i_3, j_3, k_3, l_3 ) ) );

    //A( i_3, j_3, k_3, l_3 ) = Aiso( i_3, j_3, k_3, l_3 ) + Aaniso( i_3, j_3, k_3, l_3 );
    
    return A.to_voigt_form();
}


MaterialStatus *
AnisotropicPolyconvexHyperelasticSurfaceMaterial::CreateStatus( GaussPoint *gp ) const
{
    return new StructuralMaterialStatus( gp );
}


void AnisotropicPolyconvexHyperelasticSurfaceMaterial::initializeFrom( InputRecord &ir )
{
    IsotropicPolyconvexHyperelasticSurfaceMaterial::initializeFrom( ir );
    IR_GIVE_FIELD( ir, matdir, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_MatDir );
    IR_GIVE_OPTIONAL_FIELD( ir, beta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_beta );
    IR_GIVE_OPTIONAL_FIELD( ir, eta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_eta );

    matdir.normalize();
}

FloatArray AnisotropicPolyconvexHyperelasticSurfaceMaterial::computeLocalMaterialDirection( GaussPoint *gp ) const
{
    FloatArray matdirLoc( matdir );
    NLStructuralSurfaceElement *el = dynamic_cast<NLStructuralSurfaceElement *>( gp->giveElement() );
    FloatMatrix T; // transformation matrix
    int val = el->computeLoadGToLRotationMtrx( T );
    T.resizeWithData( 3, 3 );
    matdirLoc.rotatedWith( T, 'n' );

    return matdirLoc;
}

} // end namespace oofem
