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
//#include "sm/Elements/SurfaceElements/quadsurface.h"
#include "structuralelement.h"
#include "sm/Elements/SurfaceElements/triangularsurface.h"


namespace oofem {
REGISTER_Material( AnisotropicPolyconvexHyperelasticSurfaceMaterial );


AnisotropicPolyconvexHyperelasticSurfaceMaterial::AnisotropicPolyconvexHyperelasticSurfaceMaterial( int n, Domain *d ) :
    IsotropicPolyconvexHyperelasticSurfaceMaterial( n, d ),
    P1( 3 ),
    P2( 3 )
{
}

FloatArrayF<9>
AnisotropicPolyconvexHyperelasticSurfaceMaterial::giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{


    FloatArrayF<9> PisoMat = IsotropicPolyconvexHyperelasticSurfaceMaterial::giveFirstPKSurfaceStressVector_3d( vF, normal, gp, tStep );
    Tensor2_3d Piso( PisoMat ), Fcurrent( vF ), P, S, F;

    //// prestrain
    //FloatArray vF0_temp;
    //FloatArrayF<9> indent   = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    //Tensor2_3d F0( indent );
    //if ( gp->giveElement()->giveIPValue( vF0_temp, gp, IST_PrestrainDeformationGradient, tStep ) ) {
    //    FloatArrayF<9> vF0( vF0_temp );
    //    F0 = Tensor2_3d::Tensor2_3d( vF0 );
    //}
    //F( i_3, j_3 ) = Fcurrent( i_3, k_3 ) * F0( k_3, j_3 );
    F( i_3, j_3 ) = Fcurrent( i_3, j_3 );

    //// Anisotropic part
    Tensor1_3d a = computeLocalMaterialDirection( gp, tStep );

    Tensor2_3d Paniso1, Paniso2;
    //Tensor1_3d a;
    //a( 0 ) = matdirLoc.at( 1, 1 );
    //a( 1 ) = matdirLoc.at( 2, 1 );
    //a( 2 ) = matdirLoc.at( 3, 1 );

    Paniso1( i_3, j_3 ) = this->compute_surface_d_normFa_dF( F, a )( i_3, j_3 ); 
    Paniso2( i_3, j_3 ) = this->compute_surface_d_SqrtNormCa_dF( F, a )( i_3, j_3 ); 

    //// Check how it is in the undeformed original reference frame
    // Tensor2_3d F_rotated, P_rotated;
    // this->TransformTensor2( F, F_rotated, gp );
    // this->TransformTensor2( Paniso1, P_rotated, gp );
    // Tensor1_3d a_rotated;
    // this->TransformTensor1( a, a_rotated, gp );

    double eta_t, beta_t;
    if ( this->aniso_ltf == 0 ) {
        eta_t  = this->eta;
        beta_t = this->beta;
    } else {
        eta_t  = this->eta * domain->giveFunction( this->aniso_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        beta_t = this->beta * domain->giveFunction( this->aniso_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    P( i_3, j_3 ) = eta_t * Paniso1( i_3, j_3 ) + beta_t * Paniso2( i_3, j_3 );

    // Prestrain transform
    //double J0 = this->compute_surface_determinant( F0 );
    Tensor2_3d Pfull;
    //Pfull( i_3, j_3 ) = 1. / J0 * P( i_3, k_3 ) * F0( j_3, k_3 ) + Piso( i_3, j_3 ); 
    Pfull( i_3, j_3 ) = P( i_3, j_3 ) + Piso( i_3, j_3 ); 

    // compute Cauchy stress vector
    S( i_3, j_3 ) = ( 1 / this->compute_surface_determinant( F ) ) * Pfull( i_3, k_3 ) * F( j_3, k_3 );

    auto vP = Pfull.to_voigt_form();
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
  
    // Isotropic stiffness matrix
    FloatMatrixF<9, 9> AisoMat = IsotropicPolyconvexHyperelasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dPdF( mode, normal, gp, tStep );
    Tensor4_3d Aiso( AisoMat ), A, Aaniso1, Aaniso2, Aaniso2_a, Aaniso2_b, Aaniso2_c;

    // prestrain
    //FloatArray vF0_temp;
    Tensor2_3d Fcurrent( vF ), F;
    //FloatArrayF<9> indent = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    //Tensor2_3d F0( indent );
    //if ( gp->giveElement()->giveIPValue( vF0_temp, gp, IST_PrestrainDeformationGradient, tStep ) ) {
    //    FloatArrayF<9> vF0( vF0_temp );
    //    F0 = Tensor2_3d::Tensor2_3d( vF0 );
    //}
    //F( i_3, j_3 ) = Fcurrent( i_3, k_3 ) * F0( k_3, j_3 );
    F( i_3, j_3 ) = Fcurrent( i_3, j_3 ) ;

    Tensor1_3d a = computeLocalMaterialDirection( gp, tStep );

    Aaniso1( i_3, j_3, k_3, l_3 ) = this->compute_surface_d2_normFa_dF2( F, a )( i_3, j_3, k_3, l_3 );
    Aaniso2( i_3, j_3, k_3, l_3 ) = this->compute_surface_d2_SqrtNormCa_dF2( F, a )( i_3, j_3, k_3, l_3 ); 

    double eta_t, beta_t;
    if ( this->aniso_ltf == 0 ) {
        eta_t  = this->eta;
        beta_t = this->beta;
    } else {
        eta_t  = this->eta * domain->giveFunction( this->aniso_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
        beta_t = this->beta * domain->giveFunction( this->aniso_ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    A( i_3, j_3, k_3, l_3 ) =  eta_t * Aaniso1( i_3, j_3, k_3, l_3 ) + beta_t* Aaniso2( i_3, j_3, k_3, l_3 );
    
    // Transform in case of prestrain
    //double J0 = this->compute_surface_determinant( F0 );
    Tensor4_3d Afull;
    //Afull( i_3, j_3, k_3, l_3 ) = 1. / J0 * F0( j_3, m_3 ) * A( i_3, m_3, k_3, n_3 ) * F0( l_3, n_3 ) + Aiso( i_3, j_3, k_3, l_3 ) ;
    Afull( i_3, j_3, k_3, l_3 ) = A( i_3, j_3, k_3, l_3 ) + Aiso( i_3, j_3, k_3, l_3 );


    //// Transform to to original reference frame, just for testing
    // Tensor4_3d A_rotated;
    // Tensor2_3d F_rotated;
    // Tensor1_3d a_rotated;
    // this->TransformTensor1( a, a_rotated, gp );
    // this->TransformTensor2( F, F_rotated, gp );
    // this->TransformTensor4( A, A_rotated, gp );
    // FloatMatrix Ar = A_rotated.to_voigt_form();
    // FloatArray Fr = F_rotated.to_voigt_form();


    return Afull.to_voigt_form();
}


//MaterialStatus *
//AnisotropicPolyconvexHyperelasticSurfaceMaterial::CreateStatus( GaussPoint *gp ) const
//{
//    return new StructuralMaterialStatus( gp );
//}


void AnisotropicPolyconvexHyperelasticSurfaceMaterial::initializeFrom( InputRecord &ir )
{
    IsotropicPolyconvexHyperelasticSurfaceMaterial::initializeFrom( ir );
    //IR_GIVE_FIELD( ir, matdir, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_MatDir );
    IR_GIVE_FIELD( ir, P1, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P1 );
    IR_GIVE_FIELD( ir, P2, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P2 );

    if ( (P1 - P2).computeNorm() < 1e-6) {
        OOFEM_ERROR( "The P1 and P2 points specifying surface anisotropy are the same" );
    }

    IR_GIVE_OPTIONAL_FIELD( ir, beta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_beta );
    IR_GIVE_OPTIONAL_FIELD( ir, eta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_eta );
    IR_GIVE_OPTIONAL_FIELD( ir, zeta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_zeta );

    IR_GIVE_OPTIONAL_FIELD( ir, aniso_ltf, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_anisoLTF );

    if ( this->aniso_ltf == 0 ) {
        this->aniso_ltf = this->gamma_ltf;
    }


    //matdir.normalize();
}

Tensor1_3d AnisotropicPolyconvexHyperelasticSurfaceMaterial::computeLocalMaterialDirection( GaussPoint *gp, TimeStep *tStep ) const
{
    FloatArray matdir( 3 ), Naniso(3);
    NLStructuralSurfaceElement *el = dynamic_cast<NLStructuralSurfaceElement *>( gp->giveElement() );
    
    // If the direction of anisotropy is not directly specified, construct it 
    if ( !el->giveIPValue( matdir, gp, IST_AnisotropyDirection, tStep ) ) {
        FloatArray normal = el->giveNormal( gp );
        FloatArray gpCor = gp->giveGlobalCoordinates();
        Naniso.beVectorProductOf( gpCor - P1, gpCor - P2 );
        matdir.beVectorProductOf( Naniso, normal );
        matdir.normalize();
    }

    // Transform to local system
    FloatMatrix T; // transformation matrix
    int val = el->computeLoadGToLRotationMtrx( T );
    T.resizeWithData( 3, 3 );
    matdir.rotatedWith( T, 'n' );

    // get the vector rotated around the original reference normal
    FloatArray l( { 0., 0., 1. } );
    el->giveIPValue( l, gp, IST_ElementReferenceNormal, tStep );
    Tensor1_3d lt( FloatArrayF<3>::FloatArrayF<3>( l ) ); // Create tensor from it
    
    // Contruct the rotation matrix using levicivita tensor
    Tensor2_3d eps1( 0., 0., 0., 0., 0., 1., 0., -1., 0. ), eps2( 0., 0., -1., 0., 0., 0., 1., 0., 0. ), eps3( 0., 1., 0., -1., 0., 0., 0., 0., 0. );
    Tensor2_3d Qt; 
    Qt( i_3, j_3 ) = l.at( 1 ) * eps1( j_3, i_3 ) + l.at( 2 ) * eps2( j_3, i_3 ) + l.at( 3 ) * eps3( j_3, i_3 ) + lt( i_3 ) * lt( j_3 );
        
    Tensor1_3d matdirLocT, matdirT( FloatArrayF<3>::FloatArrayF<3>( matdir ) );
    matdirLocT( i_3 ) = Qt( i_3, j_3 ) * matdirT( j_3 );

    return matdirLocT;
}

} // end namespace oofem
