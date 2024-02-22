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

#include "sm/Elements/structural3dsurfaceelement.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "sm/CrossSections/simplesurfacecrosssection.h"
#include "gaussintegrationrule.h"
#include "mathfem.h"

namespace oofem {
Structural3DSurfaceElement::Structural3DSurfaceElement( int n, Domain *aDomain ) :
    NLStructuralSurfaceElement( n, aDomain ),
    matRotation( false )
{
    cellGeometryWrapper = NULL;
}


Structural3DSurfaceElement::~Structural3DSurfaceElement()
{
    if ( cellGeometryWrapper ) {
        delete cellGeometryWrapper;
    }
}

void Structural3DSurfaceElement::initializeFrom( InputRecord &ir )
{
    NLStructuralSurfaceElement::initializeFrom( ir );
    this->matRotation = ir.hasField( _IFT_Structural3DSurfaceElement_materialCoordinateSystem );
}


void Structural3DSurfaceElement ::computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int li, int ui )
//
// Returns the [6x (nno * 3) ] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// r = ( u1,v1,w1,u2,v2,w2,...)
// Note that three rows are zero, equivalent to plane strain
{
    FloatMatrix dnx;
    FEInterpolation *interp = this->giveInterpolation();
    interp->evaldNdx( dnx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() ); // nno x 2

    answer.resize( 6, dnx.giveNumberOfRows() * 3 );
    answer.zero();

    for ( int i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
        answer.at( 1, 3 * i - 2 ) = dnx.at( i, 1 );
        answer.at( 2, 3 * i - 1 ) = dnx.at( i, 2 );
        answer.at( 6, 3 * i - 2 ) = dnx.at( i, 2 );
        answer.at( 6, 3 * i - 1 ) = dnx.at( i, 1 );
    }
}

void Structural3DSurfaceElement ::computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer )
//
// Returns the [9x(nno * 3)] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// @todo not checked if correct
{
    FloatMatrix dNdx;
    FEInterpolation *interp = this->giveInterpolation();
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() ); // nno x 2
    answer.resize( 9, dNdx.giveNumberOfRows() * 3 );

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at( 1, 3 * i - 2 ) = dNdx.at( i, 1 ); // du/dx - 1
        answer.at( 2, 3 * i - 1 ) = dNdx.at( i, 2 ); // dv/dy - 2

        answer.at( 6, 3 * i - 2 ) = dNdx.at( i, 2 ); // du/dy - 6
        answer.at( 7, 3 * i - 0 ) = dNdx.at( i, 2 ); // dw/dy - 7
        answer.at( 8, 3 * i - 0 ) = dNdx.at( i, 1 ); // dw/dx - 8
        answer.at( 9, 3 * i - 1 ) = dNdx.at( i, 1 ); // dv/dx - 9
    }
}



MaterialMode
Structural3DSurfaceElement::giveMaterialMode()
{
    return _3dMat;
}


void Structural3DSurfaceElement::giveMaterialOrientationAt( FloatArray &x, FloatArray &y, const FloatArray &lcoords )
{

    if ( this->elemLocalCS.isNotEmpty() ) { // User specified orientation
        x = {
            elemLocalCS.at( 1, 1 ), elemLocalCS.at( 2, 1 )
        };
        y = {
            -x( 1 ), x( 0 )
        };
    } else {
        FloatMatrix jac;
        this->giveInterpolation()->giveJacobianMatrixAt( jac, lcoords, *this->giveCellGeometryWrapper() );
        x.beColumnOf( jac, 1 ); // This is {dx/dxi, dy/dxi, dz/dxi}
        x.normalize();
        y = {
            -x( 1 ), x( 0 )
        };
    }
}


void Structural3DSurfaceElement::computeStressVector( FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep )
{
    auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray normal = this->giveNormal( gp );

    if ( this->matRotation ) {
        // MATERIAL ROTATION NOT IMPLEMENTED
        /////@todo This won't work properly with "useUpdatedGpRecord" (!)
        //FloatArray x, y, z;

        //this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );
        //// Transform from global c.s. to material c.s.
        //FloatArrayF<3> rotStrain = {
        //    e[0] * x[0] * x[0] + e[2] * x[0] * x[1] + e[1] * x[1] * x[1],
        //    e[0] * y[0] * y[0] + e[2] * y[0] * y[1] + e[1] * y[1] * y[1],
        //    2 * e[0] * x[0] * y[0] + 2 * e[1] * x[1] * y[1] + e[2] * ( x[1] * y[0] + x[0] * y[1] )
        //};
        //auto s = cs->giveRealSurfaceStress_3d( rotStrain, normal, gp, tStep );
        //answer = {
        //    s[0] * x[0] * x[0] + 2 * s[2] * x[0] * y[0] + s[1] * y[0] * y[0],
        //    s[0] * x[1] * x[1] + 2 * s[2] * x[1] * y[1] + s[1] * y[1] * y[1],
        //    s[1] * y[0] * y[1] + s[0] * x[0] * x[1] + s[2] * ( x[1] * y[0] + x[0] * y[1] )
        //};
    } else {
        answer = cs->giveRealSurfaceStress_3d( e, normal, gp, tStep );
    }

}

void Structural3DSurfaceElement::computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{ 
    auto cs  = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray normal = this->giveNormal( gp );
    answer  = cs->giveSurfaceStiffnessMatrix_3d( rMode, normal, gp, tStep );

    // MATERIAL ROTATION NOT IMPLEMENTED
    //if ( this->matRotation ) {
    //    FloatArray x, y;
    //    FloatMatrix Q;

    //    this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );

    //    Q = {
    //        { x( 0 ) * x( 0 ), x( 1 ) * x( 1 ), x( 0 ) * x( 1 ) },
    //        { y( 0 ) * y( 0 ), y( 1 ) * y( 1 ), y( 0 ) * y( 1 ) },
    //        { 2 * x( 0 ) * y( 0 ), 2 * x( 1 ) * y( 1 ), x( 1 ) * y( 0 ) + x( 0 ) * y( 1 ) }
    //    };
    //    answer.rotatedWith( Q, 't' );
    //}
}


void Structural3DSurfaceElement::computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{

    auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray normal = this->giveNormal( gp );
    answer = cs->giveSurfaceStiffnessMatrix_dPdF_3d( rMode, normal, gp, tStep );

    // MATERIAL ROTATION NOT IMPLEMENTED
    //if ( this->matRotation ) {
    //    FloatArray x, y;
    //    FloatMatrix Q;

    //    this->giveMaterialOrientationAt( x, y, gp->giveNaturalCoordinates() );

    //    Q = {
    //        { x( 0 ) * x( 0 ), x( 1 ) * x( 1 ), x( 0 ) * x( 1 ), x( 1 ) * x( 0 ) },
    //        { y( 0 ) * y( 0 ), y( 1 ) * y( 1 ), y( 0 ) * y( 1 ), y( 1 ) * y( 0 ) },
    //        { x( 0 ) * y( 0 ), x( 1 ) * y( 1 ), x( 0 ) * y( 1 ), x( 1 ) * y( 0 ) },
    //        { y( 0 ) * x( 0 ), y( 1 ) * x( 1 ), y( 0 ) * x( 1 ), y( 1 ) * x( 0 ) }
    //    };
    //    answer.rotatedWith( Q, 't' );
    //}
}


void Structural3DSurfaceElement::giveDofManDofIDMask( int inode, IntArray &answer ) const
{
    answer = {
        D_u, D_v, D_w
    };
}


//int Structural3DSurfaceElement::computeNumberOfDofs()
//{
//    ///@todo move one hiearchy up and generalize
//    IntArray dofIdMask;
//    this->giveDofManDofIDMask( -1, dofIdMask ); // ok for standard elements
//    return this->giveInterpolation()->giveNumberOfNodes() * dofIdMask.giveSize();
//}


void Structural3DSurfaceElement::computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this, 1, 6 );
        this->giveCrossSection()->setupIntegrationPoints( *integrationRulesArray[0], numberOfGaussPoints, this );
    }
}


double
Structural3DSurfaceElement::computeVolumeAround( GaussPoint *gp )
// Returns the portion of the receiver which is attached to gp.
{
    double determinant = fabs( this->giveInterpolation()->giveTransformationJacobian( gp->giveNaturalCoordinates(),*this->giveCellGeometryWrapper() ) );
    double weight = gp->giveWeight();

    return determinant * weight;
}

FEICellGeometry *
Structural3DSurfaceElement::giveCellGeometryWrapper()
{
    if ( !cellGeometryWrapper ) {
        cellGeometryWrapper = new FEIElementGeometryWrapper( this );
    }

    return cellGeometryWrapper;
}

} // end namespace oofem
