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
// Returns the [ 6 x (nno * 3) ] small-strain matrix {B} of the receiver in the element frame,
// evaluated at gp: rows 1, 2, 6 = (epsilon_x, epsilon_y, gamma_xy) for r = ( u1, v1, w1, u2, v2, w2, ... );
// the other three rows are zero.
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
// Returns the [ 9 x (nno * 3) ] displacement gradient matrix {BH} of the receiver, evaluated at gp, in the
// element frame (x, y in the plane, z = normal) for the rotated nodal unknowns (u, v, w). Voigt rows
// 11, 22, 33, 23, 13, 12, 32, 31, 21; the derivatives with respect to z (rows 3, 4, 5) are zero.
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
    FloatArray normal = { 0., 0., 1. }; // the reference normal in the frame of the element

    if ( this->matRotation ) {
        ///@todo Material rotation (matcs) is not implemented; the stress is not computed then.
    } else {
        answer = cs->giveRealSurfaceStress_3d( e, normal, gp, tStep );
    }

}

void Structural3DSurfaceElement::computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{ 
    auto cs  = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray normal = { 0., 0., 1. }; // the reference normal in the frame of the element
    answer  = cs->giveSurfaceStiffnessMatrix_3d( rMode, normal, gp, tStep );
    ///@todo Material rotation (matcs) is not implemented and ignored here.
}


void Structural3DSurfaceElement::giveDofManDofIDMask( int inode, IntArray &answer ) const
{
    answer = {
        D_u, D_v, D_w
    };
}


void Structural3DSurfaceElement::computeGaussPoints()
// Sets up one Gauss integration rule with numberOfGaussPoints points (the element default, or nip).
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this, 1, 6 );
        this->giveCrossSection()->setupIntegrationPoints( *integrationRulesArray[0], numberOfGaussPoints, this );
    }
}


double
Structural3DSurfaceElement::computeVolumeAround( GaussPoint *gp )
// Returns the area of the receiver attached to gp (no thickness).
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
