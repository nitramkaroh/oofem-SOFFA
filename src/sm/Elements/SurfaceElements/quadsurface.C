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

#include "sm/Elements/SurfaceElements/quadsurface.h"
#include "fei3dquadlin.h"
#include "fei2dquadlin.h"
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
REGISTER_Element( QuadSurface );

FEI2dQuadLin QuadSurface ::interpolation( 1, 2 ); // x, y of the element frame

QuadSurface ::QuadSurface( int n, Domain *aDomain ) :
    Structural3DSurfaceElement( n, aDomain )
{
    this->GtoLRotationMatrix = NULL;
    numberOfDofMans     = 4;
    numberOfGaussPoints = 4;
}


FEInterpolation *QuadSurface ::giveInterpolation() const { return &interpolation; }


void QuadSurface ::initializeFrom( InputRecord &ir )
{
    numberOfGaussPoints = 4;
    Structural3DSurfaceElement ::initializeFrom( ir );
}


FloatArray QuadSurface::giveNormal( GaussPoint *gp ) const
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
QuadSurface::giveCellGeometryWrapper()
{
    if ( cellGeometryWrapper ) {
        return cellGeometryWrapper;
    } else {
        this->computeLocalNodalCoordinates( lc );
        return ( cellGeometryWrapper = new FEIVertexListGeometryWrapper( lc ) );
    }
}


void QuadSurface ::computeLocalNodalCoordinates( std::vector<FloatArray> &lxy )
// Returns the nodal coordinates rotated to the element frame, lxy[i] = R X_i (no translation).
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }


    lxy.resize( 4 );
    for ( int i = 0; i < 4; i++ ) {
        const auto &nc = this->giveNode( i + 1 )->giveCoordinates();
        lxy[i].beProductOf( *GtoLRotationMatrix, nc );
    }
}



const FloatMatrix *
QuadSurface ::computeGtoLRotationMatrix()
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


bool QuadSurface ::computeGtoLRotationMatrix( FloatMatrix &answer )
// Returns the rotation matrix of the receiver of the size [12,12]
// r(local) = T * r(global), for each node {u,v,w}(local) = R * {u,v,w}(global)
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize( 12, 12 );
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at( 1, i ) = answer.at( 4, i + 3 ) = answer.at( 7, i + 6 ) = answer.at( 10, i + 9 ) = GtoLRotationMatrix->at( 1, i );
        answer.at( 2, i ) = answer.at( 5, i + 3 ) = answer.at( 8, i + 6 ) = answer.at( 11, i + 9 ) = GtoLRotationMatrix->at( 2, i );
        answer.at( 3, i ) = answer.at( 6, i + 3 ) = answer.at( 9, i + 6 ) = answer.at( 12, i + 9 ) = GtoLRotationMatrix->at( 3, i );
    }
 

    return 1;
}
int QuadSurface ::computeLoadGToLRotationMtrx( FloatMatrix &answer )
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


int QuadSurface ::computeGlobalCoordinates( FloatArray &answer, const FloatArray &lcoords )
{
    // Bilinear interpolation of the global nodal coordinates (node order of FEI2dQuadLin).
    // Vertex coordinates
    const auto &p1 = this->giveNode( 1 )->giveCoordinates();
    const auto &p2 = this->giveNode( 2 )->giveCoordinates();
    const auto &p3 = this->giveNode( 3 )->giveCoordinates();
    const auto &p4 = this->giveNode( 4 )->giveCoordinates();

    // Local coordinates
    double ksi = lcoords.at( 1 );
    double eta = lcoords.at( 2 );

    // Shape functions
    double n1 = ( 1. + ksi ) * ( 1. + eta ) * 0.25;
    double n2 = ( 1. - ksi ) * ( 1. + eta ) * 0.25;
    double n3 = ( 1. - ksi ) * ( 1. - eta ) * 0.25;
    double n4 = ( 1. + ksi ) * ( 1. - eta ) * 0.25;


    answer = { n1 * p1.at( 1 ) + n2 * p2.at( 1 ) + n3 * p3.at( 1 ) + n4 * p4.at( 1 ),
               n1 * p1.at( 2 ) + n2 * p2.at( 2 ) + n3 * p3.at( 2 ) + n4 * p4.at( 2 ),
               n1 * p1.at( 3 ) + n2 * p2.at( 3 ) + n3 * p3.at( 3 ) + n4 * p4.at( 3 ) };

    

    return true;
}


} // end namespace oofem
