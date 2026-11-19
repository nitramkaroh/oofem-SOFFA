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

#include <cmath>


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

    Tensor1_3d a = this->computeLocalMaterialDirection( gp, tStep );
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

    Tensor1_3d a = this->computeLocalMaterialDirection( gp, tStep );
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
    // The axis p_1 p_2 is needed only by elements without anisodir; this is checked when the direction is computed.
    const bool hasP1 = ir.hasField( _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P1 );
    const bool hasP2 = ir.hasField( _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P2 );
    if ( hasP1 != hasP2 ) {
        throw ValueInputException( ir, hasP1 ? _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P2 : _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P1,
            "p_1 and p_2 must be given together" );
    }
    this->axisGiven = hasP1;
    if ( this->axisGiven ) {
        IR_GIVE_FIELD( ir, P1, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P1 );
        IR_GIVE_FIELD( ir, P2, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P2 );
        if ( P1.giveSize() != 3 || P2.giveSize() != 3 ) {
            throw ValueInputException( ir, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P1, "p_1 and p_2 must have 3 components" );
        }
        for ( double component : P1 ) {
            if ( !std::isfinite( component ) ) {
                throw ValueInputException( ir, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P1, "must have finite components" );
            }
        }
        for ( double component : P2 ) {
            if ( !std::isfinite( component ) ) {
                throw ValueInputException( ir, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P2, "must have finite components" );
            }
        }
        const double separation = ( P1 - P2 ).computeNorm();
        if ( !std::isfinite( separation ) || separation < 1e-6 ) {
            throw ValueInputException( ir, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P1,
                "p_1 and p_2 must be different points with a finite separation" );
        }
    }
    int meridionalFlag = 0;
    IR_GIVE_OPTIONAL_FIELD( ir, meridionalFlag, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_meridional );
    if ( meridionalFlag != 0 && meridionalFlag != 1 ) {
        throw ValueInputException( ir, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_meridional, "must be 0 or 1" );
    }
    if ( meridionalFlag == 1 && !this->axisGiven ) {
        throw ValueInputException( ir, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_meridional, "needs p_1 and p_2" );
    }
    this->meridional = meridionalFlag == 1;
    readParameter( ir, this->beta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_beta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_betaLTF );
    readParameter( ir, this->eta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_eta, _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_etaLTF );
}

Tensor1_3d AnisotropicPolyconvexHyperelasticSurfaceMaterial::computeLocalMaterialDirection( GaussPoint *gp, TimeStep *tStep ) const
{
    // Only the 3D Lagrange surface elements are supported: the direction is built from the 3D position of the Gauss
    // point and rotated with the element's local frame. The IGA surfaces have no such frame and the 2D surfaces give
    // a one-component position and a nodal rotation, which would make the direction wrong (it ends up zero).
    auto el = dynamic_cast<Structural3DSurfaceElement *>( gp->giveElement() );
    if ( !el ) {
        OOFEM_ERROR( "element %d (%s) is not supported, the anisotropic surface material needs a 3D Lagrange surface element",
            gp->giveElement()->giveNumber(), gp->giveElement()->giveClassName() );
    }

    // The direction is taken on the mesh, which is the reference configuration of the material only without prestrain.
    FloatArray direction;
    if ( el->giveIPValue( direction, gp, IST_PrestrainDeformationGradient, tStep ) ) {
        OOFEM_ERROR( "element %d: the anisotropic surface material does not support prestrain (the fibre direction is defined on the mesh)",
            el->giveNumber() );
    }

    // Global direction: anisodir of the element, or m = (X - p_1) x (X - p_2), normal to the plane through the
    // Gauss point X and the axis p_1 p_2, i.e. along the circle through X around the axis
    const bool anisodirGiven = el->giveIPValue( direction, gp, IST_AnisotropyDirection, tStep );
    if ( !anisodirGiven ) {
        if ( !this->axisGiven ) {
            OOFEM_ERROR( "element %d: no fibre direction, give anisodir on the element or p_1 and p_2 in material %d",
                el->giveNumber(), this->giveNumber() );
        }
        FloatArray gpCor = gp->giveGlobalCoordinates();
        direction.beVectorProductOf( gpCor - P1, gpCor - P2 );
    }
    if ( direction.giveSize() != 3 ) {
        OOFEM_ERROR( "element %d: the anisotropic material direction must have 3 components", el->giveNumber() );
    }
    for ( double component : direction ) {
        if ( !std::isfinite( component ) ) {
            OOFEM_ERROR( "element %d: the anisotropic material direction must have finite components", el->giveNumber() );
        }
    }

    // Global to the frame of the element (first 3 x 3 block of its load rotation matrix); its third axis is the normal
    FloatMatrix T;
    if ( !el->computeLoadGToLRotationMtrx( T ) || T.giveNumberOfRows() < 3 || T.giveNumberOfColumns() < 3 ) {
        OOFEM_ERROR( "element %d: the anisotropic material needs a 3D element frame", el->giveNumber() );
    }
    T.resizeWithData( 3, 3 );
    FloatArray local;
    local.beProductOf( T, direction );

    // The surface identity diag(1, 1, 0) keeps the part in the element plane, which is normalized. Without it,
    // F (with F e3 = 0) would ignore the normal part and silently shorten a, i.e. weaken eta and beta.
    const double length = std::hypot( local.at( 1 ), local.at( 2 ) );
    if ( !std::isfinite( length ) || !( length > 0. ) ) {
        if ( anisodirGiven ) {
            OOFEM_ERROR( "element %d: anisodir has no finite nonzero part in the element plane", el->giveNumber() );
        } else {
            OOFEM_ERROR( "element %d: no fibre direction at a Gauss point on the axis p_1 p_2 or where the circle around the axis is normal to the element",
                el->giveNumber() );
        }
    }
    const double a1 = local.at( 1 ) / length, a2 = local.at( 2 ) / length;
    if ( !anisodirGiven && this->meridional ) {
        // rotated by 90 degrees in the element plane: from the parallel to the meridian
        return Tensor1_3d( FloatArrayF<3> { -a2, a1, 0. } );
    }
    return Tensor1_3d( FloatArrayF<3> { a1, a2, 0. } );
}

} // end namespace oofem
