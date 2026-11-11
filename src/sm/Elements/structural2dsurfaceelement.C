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

#include "sm/Elements/structural2dsurfaceelement.h"
#include "feinterpol1d.h"
#include "gausspoint.h"
//#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/CrossSections/simplesurfacecrosssection.h"
#include "gaussintegrationrule.h"
#include "mathfem.h"

namespace oofem {
Structural2DSurfaceElement::Structural2DSurfaceElement(int n, Domain *aDomain) :
    NLStructuralSurfaceElement(n, aDomain),
    matRotation(false)
{
    cellGeometryWrapper = NULL;
}

Structural2DSurfaceElement::~Structural2DSurfaceElement()
{
    if ( cellGeometryWrapper ) {
        delete cellGeometryWrapper;
    }
}


void
Structural2DSurfaceElement::initializeFrom(InputRecord &ir)
{
    NLStructuralSurfaceElement::initializeFrom(ir);
    matRotation = ir.hasField(_IFT_Structural2DSurfaceElement_materialCoordinateSystem); //|| this->elemLocalCS.isNotEmpty();
}


void
Structural2DSurfaceElement::postInitialize()
{
    // Element must be created before giveNumberOfNodes can be called
    StructuralElement::postInitialize();
    this->numberOfDofMans = this->giveNumberOfNodes();
}


int
Structural2DSurfaceElement::giveNumberOfNodes() const
{
    return this->giveInterpolation()->giveNumberOfNodes();
}


FEICellGeometry *
Structural2DSurfaceElement::giveCellGeometryWrapper()
{
    if ( !cellGeometryWrapper ) {
        cellGeometryWrapper = new FEIElementGeometryWrapper(this);
    }

    return cellGeometryWrapper;
}


int
Structural2DSurfaceElement::computeNumberOfDofs()
{
    ///@todo move one hiearchy up and generalize
    IntArray dofIdMask;
    this->giveDofManDofIDMask(-1, dofIdMask); // ok for standard elements
    return this->giveInterpolation()->giveNumberOfNodes() * dofIdMask.giveSize();
}


void
Structural2DSurfaceElement::giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_v
    };
}


double
Structural2DSurfaceElement::computeVolumeAround(GaussPoint *gp)
{
    // Computes the volume element dV associated with the given gp.

    double weight = gp->giveWeight();
    const FloatArray &lCoords = gp->giveNaturalCoordinates(); // local/natural coords of the gp (parent domain)
    double detJ = fabs(this->giveInterpolation()->giveTransformationJacobian(lCoords, * this->giveCellGeometryWrapper() ) );
    return detJ * weight; // dV
}


void
Structural2DSurfaceElement::computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points
    // Default: create one integration rule
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ] = std::make_unique< GaussIntegrationRule >(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}




// Axisymmetry

AxisymSurfaceElement::AxisymSurfaceElement(int n, Domain *aDomain) :
    Structural2DSurfaceElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{
}


double
AxisymSurfaceElement::computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    // note: radius is accounted by interpolation (of Fei2d*Axi type)
    /*double determinant = fabs(static_cast< FEInterpolation1d * >( this->giveInterpolation() )->
                              giveTransformationJacobian(gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() ) );*/

    FloatArray N;
    FloatArray lcoords = gp->giveNaturalCoordinates();
    FEInterpolation1d* interpolation = static_cast<FEInterpolation1d *>( this->giveInterpolation() );

    interpolation->evalN( N, lcoords, *this->giveCellGeometryWrapper() );
    

    double r = 0., x;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        x = this->giveNode( i )->giveCoordinate( 1 );
        r += x * N.at( i );
    }

    double determinant = fabs( r * interpolation->giveTransformationJacobian( lcoords, *this->giveCellGeometryWrapper() ) );


    double weight = gp->giveWeight();
    return determinant * weight;
}


void
AxisymSurfaceElement::computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    OOFEM_ERROR( "small-strain formulation is not implemented, use nlgeo 1" );
}


void
AxisymSurfaceElement::computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [ 9 x (nno*2) ] displacement gradient matrix {BH} of the receiver, evaluated at gp,
// in the frame of the surface material (1 = tangent, 2 = circumferential, 3 = in-plane normal)
// for the local nodal unknowns (u_t, u_n). Voigt rows: 1 (11) = du_t/ds, 8 (31) = du_n/ds,
// 2 (22) = u_R/R, where the circumferential stretch r/R = 1 + u_R/R holds for large deformations.
{
    FloatArray n;
    FloatMatrix dnx, Q;
    FEInterpolation1d *interp = static_cast< FEInterpolation1d * >( this->giveInterpolation() );

    interp->evalN(n, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() );
    interp->evaldNdx(dnx, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() );

    // u_local = Q u_global, so the global radial displacement of node i is
    // u_R = Q(2i-1, 2i-1) u_t + Q(2i, 2i-1) u_n, whatever the orientation of the element is
    if ( !this->computeGtoLRotationMatrix(Q) ) {
        OOFEM_ERROR( "the unknowns are expected in the element local frame" );
    }

    // reference radius, interpolated from the global radial coordinates of the nodes
    double r = 0.;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        r += this->giveNode(i)->giveCoordinate(1) * n.at(i);
    }

    int nRows = dnx.giveNumberOfRows();
    answer.resize(9, nRows * 2);
    answer.zero();

    for ( int i = 1; i <= nRows; i++ ) {
        answer.at( 1, 2 * i - 1 ) = dnx.at( i, 1 ); // du_t/ds
        answer.at( 8, 2 * i )     = dnx.at( i, 1 ); // du_n/ds
        answer.at( 2, 2 * i - 1 ) = Q.at( 2 * i - 1, 2 * i - 1 ) * n.at( i ) / r; // u_R/R
        answer.at( 2, 2 * i )     = Q.at( 2 * i, 2 * i - 1 ) * n.at( i ) / r;
    }
}


void
AxisymSurfaceElement::computeGaussPoints()
{
    // Sets up the integration rule array which contains all the Gauss points
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ] = std::make_unique< GaussIntegrationRule >(1, this, 1, 6);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}

void
AxisymSurfaceElement::computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
    //auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    //FloatArray normal = this->giveNormal( gp );
    //answer = cs->giveRealSurfaceStress_3d( e, normal, gp, tStep );
    OOFEM_ERROR( "not implemented " );

}

void
AxisymSurfaceElement::computeConstitutiveMatrixAt(FloatMatrix &answer,
                                           MatResponseMode rMode, GaussPoint *gp,
                                           TimeStep *tStep)
{
    //auto cs           = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    //FloatArray normal = this->giveNormal( gp );
    //answer            = cs->giveSurfaceStiffnessMatrix_3d( rMode, normal, gp, tStep );
    OOFEM_ERROR( "not implemented " );

}

} // end namespace oofem
