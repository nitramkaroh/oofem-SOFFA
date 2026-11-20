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

FEI1dLin AxiLineSurface ::interpolation( 1 );

AxiLineSurface ::AxiLineSurface( int n, Domain *aDomain ) :
    AxisymSurfaceElement( n, aDomain )
{
    this->GtoLRotationMatrix = NULL;
    numberOfDofMans          = 2;
    numberOfGaussPoints      = 2;
}


FEInterpolation *AxiLineSurface ::giveInterpolation() const { return &interpolation; }





FloatArray AxiLineSurface::giveNormal( GaussPoint *gp ) const
{
    // e1 = (X2 - X1) / |X2 - X1| in the xy plane; the normal (-e1_y, e1_x, 0) points to its left
    FloatArray e1, normal;
    e1 = { this->giveNode( 2 )->giveCoordinate( 1 ) - this->giveNode( 1 )->giveCoordinate( 1 ),
        this->giveNode( 2 )->giveCoordinate( 2 ) - this->giveNode( 1 )->giveCoordinate( 2 ) };

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
// Returns the nodal coordinates rotated to the element frame, lxy[i] = R X_i (no translation).
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }


    lxy.resize( 2 );
    for ( int i = 0; i < 2; i++ ) {
        const auto &nc = this->giveNode( i + 1 )->giveCoordinates();
        FloatArray xy = { nc.at( 1 ), nc.at( 2 ) }; // ignore an optional third coordinate for the planar line
        lxy[i].beProductOf( *GtoLRotationMatrix, xy );
    }
}


const FloatMatrix *
AxiLineSurface ::computeGtoLRotationMatrix()
// Returns the rotation matrix of the receiver of the size [2,2], rows e1 = (N2 - N1) / |N2 - N1| and
// the in-plane normal e2 = (-e1_y, e1_x) (the third axis of the material frame, see computeMaterialFrame)
// coords(local) = T * coords(global)
{
    if ( GtoLRotationMatrix == NULL ) {
        FloatArray e1, e2;

        // e1 = (X2 - X1) / |X2 - X1|, the unit tangent (x, y components)
        e1 = { this->giveNode( 2 )->giveCoordinate( 1 ) - this->giveNode( 1 )->giveCoordinate( 1 ),
        this->giveNode( 2 )->giveCoordinate( 2 ) - this->giveNode( 1 )->giveCoordinate( 2 ) };
        e1.normalize();

        // e2 = e1 rotated by +90 degrees, the in-plane normal to the left of N1 -> N2
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
// r(local) = T * r(global), for each node {u_t, u_n} = R * {u, v}: T = diag(R, R)
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
int AxiLineSurface::computeLoadGToLRotationMtrx( FloatMatrix &answer )
// Rotate pointwise load components before multiplication by N^T, not a whole element vector.
{
    answer = *this->computeGtoLRotationMatrix();
    return 1;
}



} // end namespace oofem
