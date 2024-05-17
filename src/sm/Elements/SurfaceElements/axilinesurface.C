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

#include "sm/Elements/SurfaceElements/axilinesurface.h"
//#include "fei3dquadlin.h"
#include "fei1dlin.h"
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
REGISTER_Element( AxiLineSurface );

FEI1dLinAxi AxiLineSurface ::interpolation( 1 );

AxiLineSurface ::AxiLineSurface( int n, Domain *aDomain ) :
    AxisymSurfaceElement( n, aDomain )
{
    this->GtoLRotationMatrix = NULL;
    numberOfDofMans          = 2;
    numberOfGaussPoints      = 2;
}


FEInterpolation *AxiLineSurface ::giveInterpolation() const { return &interpolation; }


void AxiLineSurface ::initializeFrom( InputRecord &ir )
{
    numberOfGaussPoints = 2;
    AxisymSurfaceElement ::initializeFrom( ir );
}


FloatArray AxiLineSurface::giveNormal( GaussPoint *gp ) const
{
    FloatArray e1, normal;
    e1.beDifferenceOf( this->giveNode( 2 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );

    // let us normalize e1'
    e1.normalize();
    normal.resize( 3 );
    normal.at( 1 ) = -e1.at( 2 );
    normal.at( 2 ) = e1.at( 1 );

    return normal;
}


FEICellGeometry *
AxiLineSurface::giveCellGeometryWrapper()
{
    if ( cellGeometryWrapper ) {
        return cellGeometryWrapper;
    } else {
        this->computeLocalNodalCoordinates( lc );
        return ( cellGeometryWrapper = new FEIVertexListGeometryWrapper( lc ) );
    }
}


void AxiLineSurface ::computeLocalNodalCoordinates( std::vector<FloatArray> &lxy )
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }


    lxy.resize( 2 );
    for ( int i = 0; i < 2; i++ ) {
        const auto &nc = this->giveNode( i + 1 )->giveCoordinates();
        lxy[i].beProductOf( *GtoLRotationMatrix, nc );
    }
}


const FloatMatrix *
AxiLineSurface ::computeGtoLRotationMatrix()
// Returns the rotation matrix of the receiver of the size [2,2]
// coords(local) = T * coords(global)
//
{
    if ( GtoLRotationMatrix == NULL ) {
        FloatArray e1, e2;

        // compute e1' = [N2-N1]
        e1.beDifferenceOf( this->giveNode( 2 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );
        // let us normalize e1'
        e1.normalize();

        // Compute e2
        e2.resize( 2 );
        e2.at( 1 ) = -e1.at( 2 );
        e2.at( 2 ) = e1.at( 1 );

        //
        GtoLRotationMatrix = new FloatMatrix( 2, 2 );

        for ( int i = 1; i <= 2; i++ ) {
            GtoLRotationMatrix->at( 1, i ) = e1.at( i );
            GtoLRotationMatrix->at( 2, i ) = e2.at( i );
        }
    }

    return GtoLRotationMatrix;
}


bool AxiLineSurface ::computeGtoLRotationMatrix( FloatMatrix &answer )
// Returns the rotation matrix of the receiver of the size [4,4]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v} = T * {u,v}
{
    // test if pereviously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }


    answer.resize( 4, 4 );
    answer.zero();

    for ( int i = 1; i <= 2; i++ ) {
        answer.at( 1, i ) = answer.at( 3, i + 2 ) = GtoLRotationMatrix->at( 1, i );
        answer.at( 2, i ) = answer.at( 4, i + 2 ) = GtoLRotationMatrix->at( 2, i );
    }


    return 1;
}
int AxiLineSurface ::computeLoadGToLRotationMtrx( FloatMatrix &answer )
// Returns the rotation matrix of the receiver of the size [4,4]
// f(local) = T * f(global)
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize( 4, 4 );
    answer.zero();

    for ( int i = 1; i <= 2; i++ ) {
        answer.at( 1, i ) = answer.at( 3, i + 2 ) = GtoLRotationMatrix->at( 1, i );
        answer.at( 2, i ) = answer.at( 4, i + 2 ) = GtoLRotationMatrix->at( 2, i );
    }

    return 1;
}



} // end namespace oofem
