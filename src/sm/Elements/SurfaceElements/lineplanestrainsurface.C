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

#include "sm/Elements/SurfaceElements/lineplanestrainsurface.h"
// #include "fei3dquadlin.h"
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
#include "sm/CrossSections/simplesurfacecrosssection.h"

namespace oofem {
REGISTER_Element( LinePlaneStrainSurface );

FEI1dLin LinePlaneStrainSurface ::interpolation( 1 );

LinePlaneStrainSurface ::LinePlaneStrainSurface( int n, Domain *aDomain ) :
    Structural2DSurfaceElement( n, aDomain )
{
    this->GtoLRotationMatrix = NULL;
    numberOfDofMans          = 2;
    numberOfGaussPoints      = 2;
}


FEInterpolation *LinePlaneStrainSurface ::giveInterpolation() const { return &interpolation; }


void LinePlaneStrainSurface ::initializeFrom( InputRecord &ir )
{
    numberOfGaussPoints = 2;
    Structural2DSurfaceElement ::initializeFrom( ir );
}


FloatArray LinePlaneStrainSurface::giveNormal( GaussPoint *gp ) const
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
LinePlaneStrainSurface::giveCellGeometryWrapper()
{
    if ( cellGeometryWrapper ) {
        return cellGeometryWrapper;
    } else {
        this->computeLocalNodalCoordinates( lc );
        return ( cellGeometryWrapper = new FEIVertexListGeometryWrapper( lc ) );
    }
}


void LinePlaneStrainSurface ::computeLocalNodalCoordinates( std::vector<FloatArray> &lxy )
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
LinePlaneStrainSurface ::computeGtoLRotationMatrix()
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


bool LinePlaneStrainSurface ::computeGtoLRotationMatrix( FloatMatrix &answer )
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
int LinePlaneStrainSurface ::computeLoadGToLRotationMtrx( FloatMatrix &answer )
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


void LinePlaneStrainSurface::computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int li, int ui )
//
// Returns the [ 6 x (nno*2) ] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,...,Gamma_xy) = B . r
// r = ( u1,v1,u2,v2)
// NOT TESTED
{
    OOFEM_ERROR( "THIS FUNCTION IS NOT IMPLEMENTED " );
}


void LinePlaneStrainSurface::computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer )
// Returns the [ 9 x (nno*2) ] displacement gradient matrix {BH} of the receiver, evaluated at gp,
// in the frame of the surface material (1 = tangent, 2 = out-of-plane, 3 = in-plane normal)
// for the local nodal unknowns (u_t, u_n). Voigt rows: 1 (11) = du_t/ds, 8 (31) = du_n/ds;
// the out-of-plane stretch is 1 (plane strain) and comes from the surface identity only.
{
    FloatMatrix dnx;
    FEInterpolation1d *interp = static_cast<FEInterpolation1d *>( this->giveInterpolation() );

    interp->evaldNdx( dnx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );

    int nRows = dnx.giveNumberOfRows();
    answer.resize( 9, nRows * 2 );
    answer.zero();

    for ( int i = 1; i <= nRows; i++ ) {
        answer.at( 1, 2 * i - 1 ) = dnx.at( i, 1 ); // du_t/ds
        answer.at( 8, 2 * i )     = dnx.at( i, 1 ); // du_n/ds
    }
}


void LinePlaneStrainSurface::computeStressVector( FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep )
{
    // auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    // FloatArray normal = this->giveNormal( gp );
    // answer = cs->giveRealSurfaceStress_3d( e, normal, gp, tStep );
    OOFEM_ERROR( "not implemented " );
}

void LinePlaneStrainSurface::computeConstitutiveMatrixAt( FloatMatrix &answer,
    MatResponseMode rMode, GaussPoint *gp,
    TimeStep *tStep )
{
    // auto cs           = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    // FloatArray normal = this->giveNormal( gp );
    // answer            = cs->giveSurfaceStiffnessMatrix_3d( rMode, normal, gp, tStep );
    OOFEM_ERROR( "not implemented " );
}

} // end namespace oofem
