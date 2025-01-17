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

#include "sm/Elements/SurfaceElements/quadsurfacequadratic.h"
#include "fei2dquadquad.h"
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
REGISTER_Element( QuadSurfaceQuadratic );

FEI2dQuadQuad QuadSurfaceQuadratic ::interpolation( 1, 2 );
// FEI3dQuadLin QuadSurface ::interpolation;

QuadSurfaceQuadratic ::QuadSurfaceQuadratic( int n, Domain *aDomain ) :
    Structural3DSurfaceElement( n, aDomain )
{
    this->GtoLRotationMatrix = NULL;
    numberOfDofMans          = 8;
    numberOfGaussPoints      = 9;
}


FEInterpolation *QuadSurfaceQuadratic ::giveInterpolation() const { return &interpolation; }


void QuadSurfaceQuadratic ::initializeFrom( InputRecord &ir )
{
    numberOfGaussPoints = 16;
    Structural3DSurfaceElement ::initializeFrom( ir );
}


FloatArray QuadSurfaceQuadratic::giveNormal( GaussPoint *gp ) const
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
QuadSurfaceQuadratic::giveCellGeometryWrapper()
{
    if ( cellGeometryWrapper ) {
        return cellGeometryWrapper;
    } else {
        this->computeLocalNodalCoordinates( lc );
        return ( cellGeometryWrapper = new FEIVertexListGeometryWrapper( lc ) );
    }
}


void QuadSurfaceQuadratic ::computeLocalNodalCoordinates( std::vector<FloatArray> &lxy )
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }


    lxy.resize( 8 );
    for ( int i = 0; i < 8; i++ ) {
        const auto &nc = this->giveNode( i + 1 )->giveCoordinates();
        lxy[i].beProductOf( *GtoLRotationMatrix, nc );
    }
}


const FloatMatrix *
QuadSurfaceQuadratic ::computeGtoLRotationMatrix()
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


bool QuadSurfaceQuadratic ::computeGtoLRotationMatrix( FloatMatrix &answer )
// Returns the rotation matrix of the receiver of the size [8,12]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v} = T * {u,v,w} NOT NOW
// for one node (r written transposed): {u,v,w} = T * {u,v,w}
{
    // test if pereviously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    // answer.resize( 8, 12 );
    // answer.zero();

    // for ( int i = 1; i <= 3; i++ ) {
    //     answer.at( 1, i ) = answer.at( 3, i + 3 ) = answer.at( 5, i + 6 ) = answer.at( 7, i + 9 ) = GtoLRotationMatrix->at( 1, i );
    //     answer.at( 2, i ) = answer.at( 4, i + 3 ) = answer.at( 6, i + 6 ) = answer.at( 8, i + 9 ) = GtoLRotationMatrix->at( 2, i );
    // }

    answer.resize( 24, 24 );
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at( 1, i ) = answer.at( 4, i + 3 ) = answer.at( 7, i + 6 ) = answer.at( 10, i + 9 ) = GtoLRotationMatrix->at( 1, i );
        answer.at( 2, i ) = answer.at( 5, i + 3 ) = answer.at( 8, i + 6 ) = answer.at( 11, i + 9 ) = GtoLRotationMatrix->at( 2, i );
        answer.at( 3, i ) = answer.at( 6, i + 3 ) = answer.at( 9, i + 6 ) = answer.at( 12, i + 9 ) = GtoLRotationMatrix->at( 3, i );

        answer.at( 13, i + 12 ) = answer.at( 16, i + 15 ) = answer.at( 19, i + 18 ) = answer.at( 22, i + 21 ) = GtoLRotationMatrix->at( 1, i );
        answer.at( 14, i + 12 ) = answer.at( 17, i + 15 ) = answer.at( 20, i + 18 ) = answer.at( 23, i + 21 ) = GtoLRotationMatrix->at( 2, i );
        answer.at( 15, i + 12 ) = answer.at( 18, i + 15 ) = answer.at( 21, i + 18 ) = answer.at( 24, i + 21 ) = GtoLRotationMatrix->at( 3, i );
    }


    return 1;
}
int QuadSurfaceQuadratic ::computeLoadGToLRotationMtrx( FloatMatrix &answer )
// Returns the rotation matrix of the receiver of the size [6,6]
// f(local) = T * f(global)
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


int QuadSurfaceQuadratic ::computeGlobalCoordinates( FloatArray &answer, const FloatArray &lcoords )
{
    
    auto cellgeo = FEIElementGeometryWrapper( this );
    FloatArray n;

    this->giveInterpolation()->evalN( n, lcoords, cellgeo );

    answer.resize( 3 );
    answer.zero();
    for ( int i = 1; i <= n.giveSize(); i++ ) {
        answer.at( 1 ) += n.at( i ) * cellgeo.giveVertexCoordinates( i ).at( 1 );
        answer.at( 2 ) += n.at( i ) * cellgeo.giveVertexCoordinates( i ).at( 2 );
        answer.at( 3 ) += n.at( i ) * cellgeo.giveVertexCoordinates( i ).at( 3 );
    }

    return true;
}



} // end namespace oofem
