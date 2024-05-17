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
//
// Returns the [ 6 x (nno*2) ] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,...,Gamma_xy) = B . r
// r = ( u1,v1,u2,v2)
// NOT TESTED
{
    FEInterpolation *interp = this->giveInterpolation();

    FloatArray N;
    interp->evalN(N, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() );
    double r = 0.0;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        double x = this->giveNode(i)->giveCoordinate(1);
        r += x * N.at(i);
    }

    FloatMatrix dNdx;
    interp->evaldNdx(dNdx, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() );
    answer.resize(6, dNdx.giveNumberOfRows() * 1);
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = dNdx.at(i, 1);
        answer.at(3, 2 * i - 0) = -N.at(i) / r;
        answer.at(6, 2 * i - 0) = dNdx.at(i, 1);

    }

    OOFEM_ERROR( "THIS FUNCTION IS NOT TESTED " );
}


void
AxisymSurfaceElement::computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [ 9 x (nno) ] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz = u/r, 0, 0, du/dy,  0, 0, dv/dx
///@todo not checked if correct, is dw/dz = u/r for large deformations? /JB
{
    FloatArray n;
    FloatMatrix dnx;
    FEInterpolation1d *interp = static_cast< FEInterpolation1d * >( this->giveInterpolation() );

    interp->evalN(n, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() );
    interp->evaldNdx(dnx, gp->giveNaturalCoordinates(), * this->giveCellGeometryWrapper() );


    int nRows = dnx.giveNumberOfRows();
    answer.resize(9, nRows * 2);
    answer.zero();

    double r = 0., x;
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        x  = this->giveNode(i)->giveCoordinate(1);
        r += x * n.at(i);
    }
    // mode is _3dMat !!!!!! z corresponds to outward normal here!, y is out of plane 
    for ( int i = 1; i <= nRows; i++ ) {
        answer.at( 1, 2 * i - 1 ) = dnx.at( i, 1 ); // du/dx
        answer.at( 9, 2 * i ) = dnx.at( i, 1 ); // dv/dx
    }

    for ( int i = 0; i < this->giveNumberOfDofManagers(); i++ ) {
        answer.at( 3, 2 * i + 2 ) = -n.at( i + 1 ) / r; // dw/dz
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

void
AxisymSurfaceElement::computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    auto cs  = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray normalTemp = this->giveNormal( gp );

    // surface material law is formulated for coordinate system with z as outward axis, therefore the 2d x-y formulation needs to be tranformed
    FloatArray normal;
    normal.resize( 3 );
    // transform normal
    normal.at( 1 ) = normalTemp.at( 1 );
    normal.at( 2 ) = normalTemp.at( 3 );
    normal.at( 3 ) = normalTemp.at( 2 );


    // transformation matrix
    FloatMatrix T( 9, 9 );
    T.at( 1, 1 ) = 1;
    T.at( 2, 3 ) = 1;
    T.at( 3, 2 ) = 1;
    T.at( 4, 7 ) = 1;
    T.at( 5, 6 ) = 1;
    T.at( 6, 5 ) = 1;
    T.at( 7, 4 ) = 1;
    T.at( 8, 9 ) = 1;
    T.at( 9, 8 ) = 1;

    // tranform stiffness matrix
    FloatMatrix Ktr = cs->giveSurfaceStiffnessMatrix_dPdF_3d( rMode, normal, gp, tStep );
    FloatMatrix TK;
    TK.beProductOf( T, Ktr );
    answer.beProductOf( TK, T );
}


void AxisymSurfaceElement::computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );

    FloatArray normal = this->giveNormal( gp );
    FloatArray vF;
    this->computeDeformationGradientVector( vF, gp, tStep );

    // surface material law is formulated for coordinate system with z as outward axis, therefore the 2d x-y formulation needs to be tranformed
    // transformation matrix
    FloatMatrix T( 9, 9 );
    T.at( 1, 1 )      = 1;
    T.at( 2, 3 )      = 1;
    T.at( 3, 2 )      = 1;
    T.at( 4, 7 )      = 1;
    T.at( 5, 6 )      = 1;
    T.at( 6, 5 )      = 1;
    T.at( 7, 4 )      = 1;
    T.at( 8, 9 )      = 1;
    T.at( 9, 8 )      = 1;

    FloatArray vFr;
    vFr.beProductOf( T, vF );

    FloatArray streesTr = cs->giveFirstPKSurfaceStresses( vFr, normal, gp, tStep );

    answer.beProductOf( T, streesTr );
}
void AxisymSurfaceElement::computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    this->computeVectorOf( { D_u, D_v }, VM_Total, tStep, u ); // solution vector
    if ( initialDisplacements ) {
        u.subtract( *initialDisplacements );
    }

    // Displacement gradient H = du/dX
    FloatMatrix B;
    this->computeBHmatrixAt( gp, B );
    answer.beProductOf( B, u );

    // Deformation gradient F = H + I
    answer.at( 1 ) += 1.0;
    answer.at( 3 ) += 1.0;
}

} // end namespace oofem
