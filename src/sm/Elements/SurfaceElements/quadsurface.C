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

FEI2dQuadLin QuadSurface ::interpolation( 1, 2 );
//FEI3dQuadLin QuadSurface ::interpolation;

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
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
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
// Returns the rotation matrix of the receiver of the size [8,12]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v} = T * {u,v,w} NOT NOW
// for one node (r written transposed): {u,v,w} = T * {u,v,w} 
{
    // test if pereviously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    //answer.resize( 8, 12 );
    //answer.zero();

    //for ( int i = 1; i <= 3; i++ ) {
    //    answer.at( 1, i ) = answer.at( 3, i + 3 ) = answer.at( 5, i + 6 ) = answer.at( 7, i + 9 ) = GtoLRotationMatrix->at( 1, i );
    //    answer.at( 2, i ) = answer.at( 4, i + 3 ) = answer.at( 6, i + 6 ) = answer.at( 8, i + 9 ) = GtoLRotationMatrix->at( 2, i );
    //}

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


int QuadSurface ::computeGlobalCoordinates( FloatArray &answer, const FloatArray &lcoords )
{
    //    FEInterpolation *fei = this->giveInterpolation();
    // #ifdef DEBUG
    //    if ( !fei ) {
    //        answer.clear();
    //        return false;
    //    }
    // #endif
    // fei->local2global3D( answer, lcoords, FEIElementGeometryWrapper( this ) );

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





//
//FloatMatrixF<3, 3> QuadSurface::giveJacobian( const FloatArrayF<3> &lcoords )
//{
//    // derivatives of interpolation functions
//    auto dn  = interpolation.evaldNdxi( lcoords[{ 0, 1 }] );
//    auto hk1 = dn.row<0>(); // dh(r1,r2)/dr1
//    auto hk2 = dn.row<1>(); // dh(r1,r2)/dr2
//
//    // interpolation functions - h(r1,r2)
//    //auto [h1, h2, h3, h4] = interp_lin.evalN(lcoords[{0, 1}]);
//    auto h    = interpolation.evalN( lcoords[{ 0, 1 }] );
//    double h1 = h[0];
//    double h2 = h[1];
//    double h3 = h[2];
//    double h4 = h[3];
//
//    //auto [a1, a2, a3, a4] = this->giveThickness();
//    //auto a    = this->giveThickness();
//    //double a1 = a[0];
//    //double a2 = a[1];
//    //double a3 = a[2];
//    //double a4 = a[3];
//    double a1 = 1;
//    double a2 = 1;
//    double a3 = 1;
//    double a4 = 1;
//
//    // get node coordinates
//    double r3 = lcoords.at( 3 );
//    auto c    = this->giveNodeCoordinates();
//    //auto [x1, y1, z1] = c[0];
//    //auto [x2, y2, z2] = c[1];
//    //auto [x3, y3, z3] = c[2];
//    //auto [x4, y4, z4] = c[3];
//    double x1 = c[0][0], y1 = c[0][1];
//    double x2 = c[1][0], y2 = c[1][1];
//    double x3 = c[2][0], y3 = c[2][1];
//    double x4 = c[3][0], y4 = c[3][1];
//
//    //auto [V1, V2, V3, V4] = this->giveLocalDirectorVectors();
//    auto V  = this->giveLocalDirectorVectors();
//    auto V1 = V[0], V2 = V[1], V3 = V[2], V4 = V[3];
//
//    // Jacobian Matrix
//    FloatMatrixF<3, 3> jacobianMatrix;
//    jacobianMatrix.at( 1, 1 ) = hk1.at( 1 ) * x1 + hk1.at( 2 ) * x2 + hk1.at( 3 ) * x3 + hk1.at( 4 ) * x4 + r3 / 2. * ( a1 * hk1.at( 1 ) * V1.at( 1 ) + a2 * hk1.at( 2 ) * V2.at( 1 ) + a3 * hk1.at( 3 ) * V3.at( 1 ) + a4 * hk1.at( 4 ) * V4.at( 1 ) );
//    jacobianMatrix.at( 1, 2 ) = hk1.at( 1 ) * y1 + hk1.at( 2 ) * y2 + hk1.at( 3 ) * y3 + hk1.at( 4 ) * y4 + r3 / 2. * ( a1 * hk1.at( 1 ) * V1.at( 2 ) + a2 * hk1.at( 2 ) * V2.at( 2 ) + a3 * hk1.at( 3 ) * V3.at( 2 ) + a4 * hk1.at( 4 ) * V4.at( 2 ) );
//    jacobianMatrix.at( 1, 3 ) = r3 / 2. * ( a1 * hk1.at( 1 ) * V1.at( 3 ) + a2 * hk1.at( 2 ) * V2.at( 3 ) + a3 * hk1.at( 3 ) * V3.at( 3 ) + a4 * hk1.at( 4 ) * V4.at( 3 ) );
//    jacobianMatrix.at( 2, 1 ) = hk2.at( 1 ) * x1 + hk2.at( 2 ) * x2 + hk2.at( 3 ) * x3 + hk2.at( 4 ) * x4 + r3 / 2. * ( a1 * hk2.at( 1 ) * V1.at( 1 ) + a2 * hk2.at( 2 ) * V2.at( 1 ) + a3 * hk2.at( 3 ) * V3.at( 1 ) + a4 * hk2.at( 4 ) * V4.at( 1 ) );
//    jacobianMatrix.at( 2, 2 ) = hk2.at( 1 ) * y1 + hk2.at( 2 ) * y2 + hk2.at( 3 ) * y3 + hk2.at( 4 ) * y4 + r3 / 2. * ( a1 * hk2.at( 1 ) * V1.at( 2 ) + a2 * hk2.at( 2 ) * V2.at( 2 ) + a3 * hk2.at( 3 ) * V3.at( 2 ) + a4 * hk2.at( 4 ) * V4.at( 2 ) );
//    jacobianMatrix.at( 2, 3 ) = r3 / 2. * ( a1 * hk2.at( 1 ) * V1.at( 3 ) + a2 * hk2.at( 2 ) * V2.at( 3 ) + a3 * hk2.at( 3 ) * V3.at( 3 ) + a4 * hk2.at( 4 ) * V4.at( 3 ) );
//    jacobianMatrix.at( 3, 1 ) = 1. / 2. * ( a1 * h1 * V1.at( 1 ) + a2 * h2 * V2.at( 1 ) + a3 * h3 * V3.at( 1 ) + a4 * h4 * V4.at( 1 ) );
//    jacobianMatrix.at( 3, 2 ) = 1. / 2. * ( a1 * h1 * V1.at( 2 ) + a2 * h2 * V2.at( 2 ) + a3 * h3 * V3.at( 2 ) + a4 * h4 * V4.at( 2 ) );
//    jacobianMatrix.at( 3, 3 ) = 1. / 2. * ( a1 * h1 * V1.at( 3 ) + a2 * h2 * V2.at( 3 ) + a3 * h3 * V3.at( 3 ) + a4 * h4 * V4.at( 3 ) );
//    return jacobianMatrix;
//}
//
//std::array<FloatArrayF<3>, 4>
//QuadSurface::giveNodeCoordinates()
//{
//    std::array<FloatArrayF<3>, 4> c;
//    for ( int i = 0; i < 4; ++i ) {
//        c[i] = this->giveLocalCoordinates( this->giveNode( i + 1 )->giveCoordinates() );
//    }
//    return c;
//}
//
//FloatArrayF<3>
//QuadSurface::giveLocalCoordinates( const FloatArrayF<3> &global )
//{
//    auto offset = global - FloatArrayF<3>( this->giveNode( 1 )->giveCoordinates() );
//    return dot( GtoLRotationMatrix, offset );
//}
//
//void QuadSurface::postInitialize()
//{
//    Structural3DSurfaceElement::postInitialize();
//
//    auto e = this->computeLocalBaseVectors();
//
//    for ( int i = 1; i <= 3; i++ ) {
//        GtoLRotationMatrix.at( 1, i ) = e[0].at( i );
//        GtoLRotationMatrix.at( 2, i ) = e[1].at( i );
//        GtoLRotationMatrix.at( 3, i ) = e[2].at( i );
//    }
//}
//
//std::array<FloatArrayF<3>, 3>
//QuadSurface::computeLocalBaseVectors()
//{
//    // compute A - (node2+node3)/2
//    auto coordA = 0.5 * ( FloatArrayF<3>( this->giveNode( 2 )->giveCoordinates() ) + FloatArrayF<3>( this->giveNode( 3 )->giveCoordinates() ) );
//    // compute B - (node1+node4)/2
//    auto coordB = 0.5 * ( FloatArrayF<3>( this->giveNode( 1 )->giveCoordinates() ) + FloatArrayF<3>( this->giveNode( 4 )->giveCoordinates() ) );
//    // compute e1' = [B-A]
//    auto e1 = normalize( coordB - coordA );
//
//    // compute A - (node3+node4)/2
//    auto coordC = 0.5 * ( FloatArrayF<3>( this->giveNode( 4 )->giveCoordinates() ) + FloatArrayF<3>( this->giveNode( 3 )->giveCoordinates() ) );
//    // compute B - (node2+node1)/2
//    auto coordD = 0.5 * ( FloatArrayF<3>( this->giveNode( 1 )->giveCoordinates() ) + FloatArrayF<3>( this->giveNode( 2 )->giveCoordinates() ) );
//
//    // compute help = [D-C]
//    auto help = coordD - coordC;
//    // compute e3' : vector product of e1' x help
//    auto e3 = normalize( cross( e1, help ) );
//    // now from e3' x e1' compute e2'
//    auto e2 = cross( e3, e1 );
//    return { e1, e2, e3 };
//}
//
//std::array<FloatArrayF<3>, 4>
//QuadSurface::giveLocalDirectorVectors()
//{
//    auto Vg = this->giveDirectorVectors();
//    return {
//        dot( GtoLRotationMatrix, Vg[0] ),
//        dot( GtoLRotationMatrix, Vg[1] ),
//        dot( GtoLRotationMatrix, Vg[2] ),
//        dot( GtoLRotationMatrix, Vg[3] ),
//    };
//}

//std::array<FloatArrayF<4>, 2>
//MITC4Shell::givedNdx( const FloatArrayF<3> &coords )
//{
//    auto dn    = interp_lin.evaldNdxi( coords[{ 0, 1 }] );
//    auto J     = this->giveJacobian( coords );
//    auto invJ  = inv( J );
//    auto invJ2 = invJ( { 0, 1 }, { 0, 1 } );
//    auto dndx  = dot( invJ2, dn );
//
//    auto hkx = dndx.row<0>();
//    auto hky = dndx.row<1>();
//    return { hkx, hky };
//}



} // end namespace oofem
