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

#include "sm/Elements/SurfaceElements/space3dsurfaceelementevaluator.h"
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
#include "tensor/tensor1.h"
#include "tensor/tensor2.h"
#include "tensor/tensor3.h"

namespace oofem {
/* 3D Space Elements */
void Space3dStructuralSurfaceElementEvaluator::computeNMatrixAt( FloatMatrix &answer, GaussPoint *gp )
{
    FloatArray N;
    Element *element        = this->giveElement();
    FEInterpolation *interp = element->giveInterpolation();

    interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( element, gp->giveIntegrationRule()->giveKnotSpan() ) );
    answer.beNMatrixOf( N, 3 );
}


void Space3dStructuralSurfaceElementEvaluator::computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer )
// Returns the [ 9 x (nno * 3) ] displacement gradient matrix {BH} of the receiver,
// evaluated at gp (dNdx is the surface gradient in x, y, z).
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz, dv/dz, du/dz, du/dy, dw/dy, dw/dx, dv/dx
{
    Element *element        = this->giveElement();
    FEInterpolation *interp = element->giveInterpolation();
    FloatMatrix dNdx;
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( element, gp->giveIntegrationRule()->giveKnotSpan() ) );

    answer.resize( 9, dNdx.giveNumberOfRows() * 3 );
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at( 1, 3 * i - 2 ) = dNdx.at( i, 1 ); // du/dx - 1
        answer.at( 2, 3 * i - 1 ) = dNdx.at( i, 2 ); // dv/dy - 2
        answer.at( 3, 3 * i - 0 ) = dNdx.at( i, 3 ); // dw/dz - 3
        answer.at( 4, 3 * i - 1 ) = dNdx.at( i, 3 ); // dv/dz - 4
        answer.at( 5, 3 * i - 2 ) = dNdx.at( i, 3 ); // du/dz - 5
        answer.at( 6, 3 * i - 2 ) = dNdx.at( i, 2 ); // du/dy - 6
        answer.at( 7, 3 * i - 0 ) = dNdx.at( i, 2 ); // dw/dy - 7
        answer.at( 8, 3 * i - 0 ) = dNdx.at( i, 1 ); // dw/dx - 8
        answer.at( 9, 3 * i - 1 ) = dNdx.at( i, 1 ); // dv/dx - 9
    }
}

void Space3dStructuralSurfaceElementEvaluator::computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer )
// Returns the [ 27 x (nno * 3) ] second gradient matrix of the receiver, evaluated at gp.
{
    // d2Ndx2 has the columns xx, yy, zz, yz, xz, xy of the tangential Hessian of the basis functions
    Element *element = this->giveElement();
    FloatMatrix d2Ndx2;
    element->giveInterpolation()->evald2Ndx2( d2Ndx2, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( element, gp->giveIntegrationRule()->giveKnotSpan() ) );

    // v[i][j]: position (1-based) of the component ij in the 9 component order 11, 22, 33, 23, 13, 12, 32, 31, 21
    static const int v[3][3] = { { 1, 6, 5 }, { 9, 2, 4 }, { 8, 7, 3 } };
    // c[j][k]: column of d2Ndx2 with the component jk of the symmetric Hessian
    static const int c[3][3] = { { 0, 5, 4 }, { 5, 1, 3 }, { 4, 3, 2 } };

    int nno = d2Ndx2.giveNumberOfRows();
    answer.resize( 27, 3 * nno );
    answer.zero();
    for ( int a = 0; a < nno; a++ ) {
        for ( int i = 0; i < 3; i++ ) { // dof (component of x) of node a
            for ( int j = 0; j < 3; j++ ) {
                for ( int k = 0; k < 3; k++ ) { // derivative direction: rows 9 k + 1 ... 9 k + 9
                    answer.at( 9 * k + v[i][j], 3 * a + i + 1 ) = d2Ndx2( a, c[j][k] );
                }
            }
        }
    }
}


double Space3dStructuralSurfaceElementEvaluator::computeVolumeAround( GaussPoint *gp )
{
    double determinant = fabs( this->giveElement()->giveInterpolation()->giveTransformationJacobian( gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( this->giveElement(),
            gp->giveIntegrationRule()->giveKnotSpan() ) ) );
    return determinant * gp->giveWeight();
}


bool Space3dStructuralSurfaceElementEvaluator::giveCurvatureExportValue( FloatArray &answer, InternalStateType type, GaussPoint *gp, TimeStep *tStep )
{
    if ( type != IST_CurvatureTensor ) {
        return false;
    }
    // only the bending (second gradient) materials have a curvature, and only then G is available
    if ( !this->hasSecondGradient ) {
        OOFEM_ERROR( "element %d: the curvature is computed only by the bending surface materials (surfacetensionmat2)",
            this->giveElement()->giveNumber() );
    }
    FloatArray vF, vG;
    this->computeDeformationGradientVector( vF, gp, tStep );
    this->computeGradientOfDeformationGradientVector( vG, gp, tStep );
    FloatArrayF<3> vN( static_cast<IGAElement *>( this->giveElement() )->giveNormal( gp ) );
    FloatArrayF<9> vF9( vF );
    FloatArrayF<27> vG27( vG );

    // the current normal and curvature as in the bending material (surfacetensionmat2)
    Tensor2_3d F( vF9 ), K;
    Tensor2_3d cofF = F.compute_cofactor();
    Tensor3_3d G( vG27 );
    Tensor1_3d N( vN ), n;
    n( i_3 ) = cofF( i_3, j_3 ) * N( j_3 );
    n.normalize();
    K( j_3, k_3 ) = -n( i_3 ) * G( i_3, j_3, k_3 );

    // VTK tensors are written row by row
    answer = { K( 0, 0 ), K( 0, 1 ), K( 0, 2 ), K( 1, 0 ), K( 1, 1 ), K( 1, 2 ), K( 2, 0 ), K( 2, 1 ), K( 2, 2 ) };
    return true;
}

} // end namespace oofem
