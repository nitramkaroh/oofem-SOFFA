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
{
    Tensor2_3d F( vF );
    Tensor2_3d Piso = this->computeMembraneStress( F, giveSurfaceIdentity( normal ), tStep );

    Tensor1_3d a = this->computeLocalMaterialDirection( gp, tStep, normal );
    Tensor2_3d Paniso1, Paniso2, P, Pfull;
    Paniso1( i_3, j_3 ) = this->compute_surface_d_normFa_dF( F, a )( i_3, j_3 );
    Paniso2( i_3, j_3 ) = this->compute_surface_d_SqrtNormCa_dF( F, a )( i_3, j_3 );

    double eta_t  = this->giveParameter( this->eta, calm_par_no, tStep );
    double beta_t = this->giveParameter( this->beta, calm_par_no, tStep );
    P( i_3, j_3 ) = eta_t * Paniso1( i_3, j_3 ) + beta_t * Paniso2( i_3, j_3 );

    Pfull( i_3, j_3 ) = P( i_3, j_3 ) + Piso( i_3, j_3 );
    this->storeStress( gp, vF, F, Pfull );
    return Pfull.to_voigt_form();
}


FloatMatrixF<9, 9>
AnisotropicPolyconvexHyperelasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    Tensor2_3d F( vF );
    Tensor4_3d Aiso = this->computeMembraneTangent( F, giveSurfaceIdentity( normal ), tStep );

    Tensor1_3d a = this->computeLocalMaterialDirection( gp, tStep, normal );
    Tensor4_3d Aaniso1, Aaniso2, A, Afull;
    Aaniso1( i_3, j_3, k_3, l_3 ) = this->compute_surface_d2_normFa_dF2( F, a )( i_3, j_3, k_3, l_3 );
    Aaniso2( i_3, j_3, k_3, l_3 ) = this->compute_surface_d2_SqrtNormCa_dF2( F, a )( i_3, j_3, k_3, l_3 );

    double eta_t  = this->giveParameter( this->eta, calm_par_no, tStep );
    double beta_t = this->giveParameter( this->beta, calm_par_no, tStep );
    A( i_3, j_3, k_3, l_3 ) = eta_t * Aaniso1( i_3, j_3, k_3, l_3 ) + beta_t * Aaniso2( i_3, j_3, k_3, l_3 );

    Afull( i_3, j_3, k_3, l_3 ) = A( i_3, j_3, k_3, l_3 ) + Aiso( i_3, j_3, k_3, l_3 );
    auto answer = Afull.to_voigt_form();
    if ( this->checkTangentTolerance != 0. ) {
        this->checkTangentsByFiniteDifferences( answer, normal, gp, tStep );
    }
    return answer;
}


void AnisotropicPolyconvexHyperelasticSurfaceMaterial::initializeFrom( InputRecord &ir )
{
    IsotropicPolyconvexHyperelasticSurfaceMaterial::initializeFrom( ir );
    IR_GIVE_FIELD( ir, P1, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P1 );
    IR_GIVE_FIELD( ir, P2, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P2 );
    if ( P1.giveSize() != 3 || P2.giveSize() != 3 ) {
        throw ValueInputException( ir, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P1, "p_1 and p_2 must have 3 components" );
    }
    if ( ( P1 - P2 ).computeNorm() < 1e-6 ) {
        OOFEM_ERROR( "The P1 and P2 points specifying surface anisotropy are the same" );
    }
    readParameter( ir, this->beta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_beta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_betaLTF );
    readParameter( ir, this->eta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_eta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_etaLTF );
}

Tensor1_3d AnisotropicPolyconvexHyperelasticSurfaceMaterial::computeLocalMaterialDirection( GaussPoint *gp, TimeStep *tStep, const FloatArray &normal ) const
{
    FloatArray matdir( 3 ), Naniso( 3 );
    // Only the 3D Lagrange surface elements are supported: the direction is built from the 3D position of the Gauss
    // point and rotated with the element's local frame. The IGA surfaces have no such frame and the 2D surfaces give
    // a one-component position and a nodal rotation, which would make the direction wrong (it ends up zero).
    auto el = dynamic_cast<Structural3DSurfaceElement *>( gp->giveElement() );
    if ( !el ) {
        OOFEM_ERROR( "element %d (%s) is not supported, the anisotropic surface material needs a 3D Lagrange surface element",
            gp->giveElement()->giveNumber(), gp->giveElement()->giveClassName() );
    }

    // Unless the element gives the direction (anisodir), (X - p_1) x (X - p_2) x n in global coordinates
    if ( !el->giveIPValue( matdir, gp, IST_AnisotropyDirection, tStep ) ) {
        // global geometry: normal of the element, Gauss point and p_1, p_2
        FloatArray geometricNormal = el->giveNormal( gp );
        FloatArray gpCor           = gp->giveGlobalCoordinates();
        Naniso.beVectorProductOf( gpCor - P1, gpCor - P2 );
        matdir.beVectorProductOf( Naniso, geometricNormal );
        matdir.normalize();
    }

    // Global to the frame of the element (first 3 x 3 block of its load rotation matrix)
    FloatMatrix T; // transformation matrix
    int val = el->computeLoadGToLRotationMtrx( T );
    T.resizeWithData( 3, 3 );
    matdir.rotatedWith( T, 'n' );

    // Rotate by 90 degrees around the reference normal in the frame of the element (e3, or the normal of the
    // original configuration with prestrain): a is then the tangential part of (X - p_1) x (X - p_2)
    const FloatArray &l = normal;
    FloatArrayF<3> l3( l );
    Tensor1_3d lt( l3 ); // Create tensor from it

    // Rotation matrix l x l + [l]_x, built with the Levi-Civita tensor
    Tensor2_3d eps1( 0., 0., 0., 0., 0., 1., 0., -1., 0. ), eps2( 0., 0., -1., 0., 0., 0., 1., 0., 0. ), eps3( 0., 1., 0., -1., 0., 0., 0., 0., 0. );
    Tensor2_3d Qt; 
    Qt( i_3, j_3 ) = l.at( 1 ) * eps1( j_3, i_3 ) + l.at( 2 ) * eps2( j_3, i_3 ) + l.at( 3 ) * eps3( j_3, i_3 ) + lt( i_3 ) * lt( j_3 );
        
    FloatArrayF<3> matdir3( matdir );
    Tensor1_3d matdirLocT, matdirT( matdir3 );
    matdirLocT( i_3 ) = Qt( i_3, j_3 ) * matdirT( j_3 );

    return matdirLocT;
}

} // end namespace oofem
