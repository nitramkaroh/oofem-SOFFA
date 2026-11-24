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

#include "sm/Elements/igasurfaceelements.h"
#include "sm/CrossSections/simplesurfacecrosssection.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "crosssection.h"
#include "mathfem.h"
#include "iga/iga.h"
#include "iga/feinurbs.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element( BsplinePlaneStrainSurfaceElement );


BsplinePlaneStrainSurfaceElement ::BsplinePlaneStrainSurfaceElement( int n, Domain *aDomain ) :
    IGAElement( n, aDomain ), Space2dStructuralSurfaceElementEvaluator() {}


void BsplinePlaneStrainSurfaceElement ::initializeFrom( InputRecord &ir )
{
    Space2dStructuralSurfaceElementEvaluator::initializeFrom( ir );
    IGAElement ::initializeFrom( ir );
}


int BsplinePlaneStrainSurfaceElement ::checkConsistency()
{
    BSplineInterpolation *interpol = static_cast<BSplineInterpolation *>( this->giveInterpolation() );
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints( 1 ) ) {
        OOFEM_WARNING( "number of control points mismatch" );
        return 0;
    }
    // the surface evaluator passes the normal to the surface material through it
    if ( !dynamic_cast<SimpleSurfaceCrossSection *>( this->giveCrossSection() ) ) {
        OOFEM_WARNING( "element %d (%s): a SimpleSurfaceCrossSection (simplesurfacecs) is required", this->giveNumber(), this->giveClassName() );
        return 0;
    }
    return 1;
}

FloatArray BsplinePlaneStrainSurfaceElement::giveNormal( GaussPoint *gp ) const
{
    // (t_y, -t_x, 0) / |t| in the knot span of the rule of gp (an interior knot takes the tangent of that span);
    // the boundary number of boundaryEvalNormal is ignored by BSpline2dLineInterpolation
    FloatArray normal;
    this->interpolation.boundaryEvalNormal( normal, 1, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    return normal;
};


//////////////////
// axisym

REGISTER_Element( BsplineAxisymSurfaceElement );


BsplineAxisymSurfaceElement ::BsplineAxisymSurfaceElement( int n, Domain *aDomain ) :
    IGAElement( n, aDomain ), AxisymSurfaceElementEvaluator() {}


void BsplineAxisymSurfaceElement ::initializeFrom( InputRecord &ir )
{
    AxisymSurfaceElementEvaluator::initializeFrom( ir );
    IGAElement ::initializeFrom( ir );
}


int BsplineAxisymSurfaceElement ::checkConsistency()
{
    BSplineInterpolation *interpol = static_cast<BSplineInterpolation *>( this->giveInterpolation() );
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints( 1 ) ) {
        OOFEM_WARNING( "number of control points mismatch" );
        return 0;
    }
    // the surface evaluator passes the normal to the surface material through it
    if ( !dynamic_cast<SimpleSurfaceCrossSection *>( this->giveCrossSection() ) ) {
        OOFEM_WARNING( "element %d (%s): a SimpleSurfaceCrossSection (simplesurfacecs) is required", this->giveNumber(), this->giveClassName() );
        return 0;
    }
    return 1;
}

FloatArray BsplineAxisymSurfaceElement::giveNormal( GaussPoint *gp ) const
{
    // as for the plane strain curve: (t_y, -t_x, 0) / |t| in the knot span of the rule of gp
    FloatArray normal;
    this->interpolation.boundaryEvalNormal( normal, 1, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    return normal;
};


//////////////////
// surface in space

REGISTER_Element( Bspline3dSurfaceElement );


Bspline3dSurfaceElement ::Bspline3dSurfaceElement( int n, Domain *aDomain ) :
    Bspline3dSurfaceElement( n, aDomain, std::make_unique<BSpline3dSurfaceInterpolation>() ) {}


Bspline3dSurfaceElement ::Bspline3dSurfaceElement( int n, Domain *aDomain, std::unique_ptr<BSpline3dSurfaceInterpolation> interp ) :
    IGAElement( n, aDomain ), Space3dStructuralSurfaceElementEvaluator(), interpolation( std::move( interp ) ) {}


void Bspline3dSurfaceElement ::initializeFrom( InputRecord &ir )
{
    Space3dStructuralSurfaceElementEvaluator::initializeFrom( ir );
    IGAElement ::initializeFrom( ir );
}


int Bspline3dSurfaceElement ::checkConsistency()
{
    // nu x nv control points (u index fastest), one node each
    int nu = this->interpolation->giveNumberOfControlPoints( 1 );
    int nv = this->interpolation->giveNumberOfControlPoints( 2 );
    if ( giveNumberOfDofManagers() != nu * nv ) {
        OOFEM_WARNING( "element %d (%s): %d nodes given, the knot vectors and degrees need %d x %d control points",
            this->giveNumber(), this->giveClassName(), giveNumberOfDofManagers(), nu, nv );
        return 0;
    }
    // the interpolation reads all three coordinates of the control points
    for ( int i = 1; i <= giveNumberOfDofManagers(); i++ ) {
        if ( this->giveNode( i )->giveCoordinates().giveSize() != 3 ) {
            OOFEM_WARNING( "element %d (%s): node %d needs three coordinates (coords 3 x y z)",
                this->giveNumber(), this->giveClassName(), this->giveNode( i )->giveNumber() );
            return 0;
        }
    }
    // the surface evaluator passes the normal to the surface material through it
    if ( !dynamic_cast<SimpleSurfaceCrossSection *>( this->giveCrossSection() ) ) {
        OOFEM_WARNING( "element %d (%s): a SimpleSurfaceCrossSection (simplesurfacecs) is required", this->giveNumber(), this->giveClassName() );
        return 0;
    }
    return 1;
}


FloatArray Bspline3dSurfaceElement::giveNormal( GaussPoint *gp ) const
{
    // A_1 x A_2 / |A_1 x A_2| in the knot spans of the rule of gp (one-sided tangents on a knot line)
    FloatArray normal;
    this->interpolation->evalSurfaceNormal( normal, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    return normal;
}


//////////////////
// NURBS surface in space

REGISTER_Element( Nurbs3dSurfaceElement );


Nurbs3dSurfaceElement ::Nurbs3dSurfaceElement( int n, Domain *aDomain ) :
    Bspline3dSurfaceElement( n, aDomain, std::make_unique<NURBS3dSurfaceInterpolation>() ) {}


} // end namespace oofem
