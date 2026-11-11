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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#include "sm/Elements/SurfaceElements/triangularsurface.h"
#include "fei3dquadlin.h"
#include "fei2dtrlin.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "floatmatrixf.h"
#include "floatarrayf.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "crosssection.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element( TriangularSurface );

FEI2dTrLin TriangularSurface ::interpolation( 1, 2 );

TriangularSurface ::TriangularSurface( int n, Domain *aDomain ) :
    Structural3DSurfaceElement( n, aDomain )
{
    this->GtoLRotationMatrix = NULL;
    numberOfDofMans          = 3;
    numberOfGaussPoints      = 1;
}


FEInterpolation *TriangularSurface ::giveInterpolation() const { return &interpolation; }


void TriangularSurface ::initializeFrom( InputRecord &ir )
{
    numberOfGaussPoints = 1;
    Structural3DSurfaceElement ::initializeFrom( ir );

    IR_GIVE_OPTIONAL_FIELD( ir, this->prestrain, _IFT_TriangularSurface_prestrain );
    if ( this->prestrain ) {
        FloatArray F0_temp;
        IR_GIVE_FIELD( ir, F0_temp, _IFT_TriangularSurface_f0 );
        if ( F0_temp.giveSize() != 9 ) {
            throw ValueInputException( ir, _IFT_TriangularSurface_f0, "must have 9 components" );
        }
        this->F0 = FloatArrayF<9>( F0_temp );
    }

    FloatArray normal;
    IR_GIVE_OPTIONAL_FIELD( ir, normal, _IFT_TriangularSurface_normal );
    if ( normal.giveSize() > 0 ) {
        if ( normal.giveSize() != 3 || normal.computeNorm() == 0. ) {
            throw ValueInputException( ir, _IFT_TriangularSurface_normal, "must be a nonzero vector with 3 components" );
        }
        if ( !this->prestrain ) {
            throw ValueInputException( ir, _IFT_TriangularSurface_normal, "is the normal of the original configuration of the prestrain, it needs prestrain and f0" );
        }
        normal.normalize();
        this->referenceNormal        = FloatArrayF<3>( normal );
        this->isReferenceNormalGiven = 1;
    }

    FloatArray adir;
    IR_GIVE_OPTIONAL_FIELD( ir, adir, _IFT_TriangularSurface_anisodir );
    if ( adir.giveSize() > 0 ) {
        if ( adir.giveSize() != 3 ) {
            throw ValueInputException( ir, _IFT_TriangularSurface_anisodir, "must have 3 components" );
        }
        this->anisoDir      = FloatArrayF<3>( adir );
        this->anisoDirGiven = 1;
    }
}


int TriangularSurface::checkConsistency()
{
    int result = Structural3DSurfaceElement::checkConsistency();
    if ( this->isReferenceNormalGiven ) {
        // the given normal must be the one the prestrain implies (up to the sign, which the energy does not see)
        GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint( 0 );
        FloatArray vF = { 1., 1., 0., 0., 0., 0., 0., 0., 0. }, vFmat, N;
        Tensor2_3d F0s;
        double J0;
        this->giveMaterialKinematics( vF, gp, nullptr, vFmat, N, F0s, J0 );
        FloatMatrix R;
        FloatArray given( this->referenceNormal );
        if ( this->giveRotationMatrix( R ) ) {
            IntArray pos( { 1, 2, 3 } );
            FloatMatrix R3x3;
            R3x3.beSubMatrixOf( R, pos, pos );
            given.beProductOf( R3x3, FloatArray( this->referenceNormal ) );
        }
        if ( 1. - std::fabs( given.dotProduct( N ) ) > 1.e-8 ) {
            OOFEM_ERROR( "element %d: the given normal does not belong to the prestrain f0, which maps the normal (%g %g %g) of the element frame to zero",
                this->giveNumber(), N.at( 1 ), N.at( 2 ), N.at( 3 ) );
        }
    }
    return result;
}


FloatArray TriangularSurface::giveNormal( GaussPoint *gp ) const
{
    FloatArray e1, e2, e3, help;
    e1.beDifferenceOf( this->giveNode( 2 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );
    help.beDifferenceOf( this->giveNode( 3 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );

    // let us normalize e1'
    e1.normalize();

    // compute e3' : vector product of e1' x help
    e3.beVectorProductOf( e1, help );
    e3.normalize();

    return e3;
}


FEICellGeometry *
TriangularSurface::giveCellGeometryWrapper()
{
    if ( cellGeometryWrapper ) {
        return cellGeometryWrapper;
    } else {
        this->computeLocalNodalCoordinates( lc );
        return ( cellGeometryWrapper = new FEIVertexListGeometryWrapper( lc ) );
    }
}


void TriangularSurface ::computeLocalNodalCoordinates( std::vector<FloatArray> &lxy )
// Returns the nodal coordinates rotated to the element frame, lxy[i] = R X_i (no translation).
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }


    lxy.resize( 3 );
    for ( int i = 0; i < 3; i++ ) {
        const auto &nc = this->giveNode( i + 1 )->giveCoordinates();
        lxy[i].beProductOf( *GtoLRotationMatrix, nc );
    }
}


const FloatMatrix *
TriangularSurface ::computeGtoLRotationMatrix()
// Returns the rotation matrix of the receiver of the size [3,3]
// coords(local) = T * coords(global)
//
// local coordinate (described by vector triplet e1',e2',e3') is defined as follows:
//
// e1'    : [N2-N1]    Ni - means i - th node
// help   : [N3-N1]
// e3'    : e1' x help
// e2'    : e3' x e1'
{
    if ( GtoLRotationMatrix == NULL ) {
        FloatArray e1, e2, e3, help;

        // compute e1' = [N2-N1]  and  help = [N3-N1]
        e1.beDifferenceOf( this->giveNode( 2 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );
        help.beDifferenceOf( this->giveNode( 3 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );

        // let us normalize e1'
        e1.normalize();

        // compute e3' : vector product of e1' x help
        e3.beVectorProductOf( e1, help );
        // let us normalize
        e3.normalize();

        // now from e3' x e1' compute e2'
        e2.beVectorProductOf( e3, e1 );

        //
        GtoLRotationMatrix = new FloatMatrix( 3, 3 );

        for ( int i = 1; i <= 3; i++ ) {
            GtoLRotationMatrix->at( 1, i ) = e1.at( i );
            GtoLRotationMatrix->at( 2, i ) = e2.at( i );
            GtoLRotationMatrix->at( 3, i ) = e3.at( i );
        }
    }
    return GtoLRotationMatrix;
}


bool TriangularSurface ::computeGtoLRotationMatrix( FloatMatrix &answer )
// Returns the rotation matrix of the receiver of the size [9,9]
// r(local) = T * r(global), for each node {u,v,w}(local) = R * {u,v,w}(global)
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize( 9, 9 );
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at( 1, i ) = answer.at( 4, i + 3 ) = answer.at( 7, i + 6 )  = GtoLRotationMatrix->at( 1, i );
        answer.at( 2, i ) = answer.at( 5, i + 3 ) = answer.at( 8, i + 6 )  = GtoLRotationMatrix->at( 2, i );
        answer.at( 3, i ) = answer.at( 6, i + 3 ) = answer.at( 9, i + 6 )  = GtoLRotationMatrix->at( 3, i );
    }


    return 1;
}
int TriangularSurface ::computeLoadGToLRotationMtrx( FloatMatrix &answer )
// Returns the rotation matrix of the receiver of the size [6,6], two blocks R
// f(local) = T * f(global); the anisotropic surface material rotates anisodir with its first block
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize( 6, 6 );
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at( 1, i ) = answer.at( 4, i + 3 ) = GtoLRotationMatrix->at( 1, i );
        answer.at( 2, i ) = answer.at( 5, i + 3 ) = GtoLRotationMatrix->at( 2, i );
        answer.at( 3, i ) = answer.at( 6, i + 3 ) = GtoLRotationMatrix->at( 3, i );
    }

    return 1;
}

int TriangularSurface ::computeGlobalCoordinates( FloatArray &answer, const FloatArray &lcoords )
{
    // Linear interpolation of the global nodal coordinates (area coordinates l1, l2, l3 of FEI2dTrLin).
    // Vertex coordinates
    const auto &p1 = this->giveNode( 1 )->giveCoordinates();
    const auto &p2 = this->giveNode( 2 )->giveCoordinates();
    const auto &p3 = this->giveNode( 3 )->giveCoordinates();

    // Local coordinates
    double l1 = lcoords.at( 1 );
    double l2 = lcoords.at( 2 );
    double l3 = 1.0 - l1 - l2;

    answer.resize( 3 );
    answer.zero();
    answer.at( 1 ) = l1 * p1.at( 1 ) + l2 * p2.at( 1 ) + l3 * p3.at( 1 );
    answer.at( 2 ) = l1 * p1.at( 2 ) + l2 * p2.at( 2 ) + l3 * p3.at( 2 );
    answer.at( 3 ) = l1 * p1.at( 3 ) + l2 * p2.at( 3 ) + l3 * p3.at( 3 );

    return true;
}

int TriangularSurface::giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep )
{
    if ( type == IST_PrestrainDeformationGradient ) {
        // f0 is given in global coordinates; the element works in its frame: R F0 R^T
        if ( this->prestrain ) {
            answer = FloatArray( this->F0 );
            FloatMatrix R;
            if ( this->giveRotationMatrix( R ) ) {
                IntArray pos( { 1, 2, 3 } );
                FloatMatrix R3x3, Fmat, FRt, RFRt;
                R3x3.beSubMatrixOf( R, pos, pos );

                Fmat.beMatrixForm( answer );
                FRt.beProductTOf( Fmat, R3x3 );
                RFRt.beProductOf( R3x3, FRt );

                answer.beVectorForm( RFRt );
            }

            return 1;
        } else {
            return 0;
        }  
    } else if ( type == IST_DisplacementVector ) {
        // the unknowns are in the element frame, the answer is the global displacement at gp
        FloatArray u;
        FloatMatrix N, G2L;
        this->computeVectorOf( VM_Total, tStep, u );

        if ( this->computeGtoLRotationMatrix( G2L ) ) {
            u.rotatedWith( G2L, 't' );
        }
        this->computeNmatrixAt( gp->giveSubPatchCoordinates(), N );
        answer.beProductOf( N, u );
        return 1;
    } else if ( type == IST_AnisotropyDirection ){
        if ( this->anisoDirGiven ) {
            answer = FloatArray( this->anisoDir );
            return 1;
        } else {
            return 0;
        }
    }
    else
    {
        return Structural3DSurfaceElement::giveIPValue( answer, gp, type, tStep );
    }
}

} // end namespace oofem
