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

#include "sm/Elements/SurfaceElements/space2dsurfaceelementevaluator.h"
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
#include "iga/feibspline.h"
#include "sm/Materials/structuralsurfacems.h"
#include "tensor/tensor1.h"
#include "tensor/tensor2.h"
#include "tensor/tensor3.h"

namespace oofem {
void Space2dStructuralSurfaceElementEvaluator ::computeNMatrixAt( FloatMatrix &answer, GaussPoint *gp )
{
    FloatArray N;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    answer.beNMatrixOf( N, 2 );
}


void Space2dStructuralSurfaceElementEvaluator::computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    // Returns the [ 9 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
    // evaluated at gp (dNdx is already the surface gradient in x, y).

    FloatMatrix dNdx;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    // one row of dNdx per nonzero basis function of the span; columns u, v of each control point
    int nno = dNdx.giveNumberOfRows();
    answer.resize( 9, nno * 2 );
    answer.zero();

    for ( int i = 1; i <= nno; i++ ) {
        answer.at( 1, 2 * i - 1 ) = dNdx.at( i, 1 ); // du/dx -1
        answer.at( 2, 2 * i - 0 ) = dNdx.at( i, 2 ); // dv/dy -2
        answer.at( 6, 2 * i - 1 ) = dNdx.at( i, 2 ); // du/dy -6
        answer.at( 9, 2 * i - 0 ) = dNdx.at( i, 1 ); // dv/dx -9
    }
}

void Space2dStructuralSurfaceElementEvaluator::computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    // Returns the [ 27 x (nno*2) ] second gradient matrix of the receiver, evaluated at gp.
    // d2Ndx2 has the columns xx, yy, xy of d^2 N / dS^2 (t x t).
    /// @todo No patch test of the G matrices yet (surface and bulk).

    FloatMatrix d2Ndx2;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    interp->evald2Ndx2( d2Ndx2, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    // row 9 (k - 1) + v(i, j) holds d^2 x_i / dX_j dX_k; nothing depends on z, so rows 19-27 (k = z) and the
    // components with i or j = z stay zero
    int nno = d2Ndx2.giveNumberOfRows();
    answer.resize( 27, nno * 2 );
    answer.zero();

    for ( int i = 1; i <= nno; i++ ) {
        // derivative with respect to x (rows 1-9)
        answer.at( 1, 2 * i - 1 ) = d2Ndx2.at( i, 1 ); // d2u/dxdx
        answer.at( 2, 2 * i - 0 ) = d2Ndx2.at( i, 3 ); // d2v/dydx
        answer.at( 6, 2 * i - 1 ) = d2Ndx2.at( i, 3 ); // d2u/dydx
        answer.at( 9, 2 * i - 0 ) = d2Ndx2.at( i, 1 ); // d2v/dxdx

        // derivative with respect to y (rows 10-18)
        answer.at( 1 + 9, 2 * i - 1 ) = d2Ndx2.at( i, 3 ); // d2u/dxdy
        answer.at( 2 + 9, 2 * i - 0 ) = d2Ndx2.at( i, 2 ); // d2v/dydy
        answer.at( 6 + 9, 2 * i - 1 ) = d2Ndx2.at( i, 2 ); // d2u/dydy
        answer.at( 9 + 9, 2 * i - 0 ) = d2Ndx2.at( i, 3 ); // d2v/dxdy
    }
}


double Space2dStructuralSurfaceElementEvaluator ::computeVolumeAround( GaussPoint *gp )
{
    // dL = |X_,xi| dxi; the Gauss weight already contains the length of the knot span in xi
    double determinant, weight, volume;
    determinant = fabs( this->giveElement()->giveInterpolation()->giveTransformationJacobian( gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( this->giveElement(),
            gp->giveIntegrationRule()->giveKnotSpan() ) ) );
    weight      = gp->giveWeight();
    volume      = determinant * weight;

    return volume;
}


bool Space2dStructuralSurfaceElementEvaluator::giveCurvatureExportValue( FloatArray &answer, InternalStateType type, GaussPoint *gp, int RuleIndex, TimeStep *tStep )
{
    if ( type != IST_CurvatureTensor ) {
        return false;
    }
    // only the bending (second gradient) materials have a curvature, and only then G is available
    if ( !this->hasSecondGradient ) {
        OOFEM_ERROR( "element %d: the curvature is computed only by the bending surface materials (surfacetensionmat2)",
            this->giveElement()->giveNumber() );
    }
    // F and grad F at the export point, in the knot span of its export rule
    FloatArray vF, vG;
    this->computeDeformationGradientVector( vF, gp, tStep );
    this->computeGradientOfDeformationGradientVector( vG, gp, tStep );
    FloatArrayF<3> vN( static_cast<IGAElement *>( this->giveElement() )->giveNormal( gp ) );
    FloatArrayF<9> vF9( vF );
    FloatArrayF<27> vG27( vG );

    // the current normal and curvature as in the bending material (surfacetensionmat2):
    // n = cof(F) N / |cof(F) N|, K_jk = -n_i (grad F)_ijk
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


/////////////////////
// Axisymmetric formulation

void AxisymSurfaceElementEvaluator::computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    // [ 9 x (nno*2) ] displacement gradient matrix: the plane one with the hoop part, F_33 = 1 + u / r
    // (the hoop circumference grows from 2 pi r to 2 pi (r + u), r the reference radius at gp)
    Space2dStructuralSurfaceElementEvaluator::computeBHmatrixAt( gp, answer );

    FloatArray N; // nonzero basis functions of the span
    gp->giveElement()->giveInterpolation()->evalN( N, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    double r = this->computeRadiusAt( gp );

    for ( int i = 0; i < N.giveSize(); i++ ) {
        answer.at( 3, 2 * i + 1 ) = N.at( i + 1 ) / r; // u / r
    }
}

void AxisymSurfaceElementEvaluator::computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    FloatMatrix ders;
    auto interp = dynamic_cast< const BSpline2dLineInterpolation * >( gp->giveElement()->giveInterpolation() );
    if ( !interp ) {
        OOFEM_ERROR( "element %d: the second gradient needs a B-spline curve (BSpline2dLineInterpolation)", gp->giveElement()->giveGlobalNumber() );
    }
    /// @todo No patch test of the G matrices yet (surface and bulk).

    // ders( d, k ) = d-th derivative of the k-th nonzero basis function with respect to xi
    auto knotSpan = gp->giveIntegrationRule()->giveKnotSpan();
    FEIIGAElementGeometryWrapper cellgeo( gp->giveElement(), knotSpan );
    interp->evalDerivatives( 2, ders, gp->giveNaturalCoordinates(), cellgeo );

    // reference curve X(xi) = (R, Y) at gp: R_,xi, Y_,xi, R_,xixi, Y_,xixi and the radius R
    double Rxi = 0., Yxi = 0., Rxixi = 0., Yxixi = 0.;
    double R = 0.; // undeformed radius (X1)

    IntArray mask;
    interp->giveKnotSpanBasisFuncMask( *knotSpan, mask );

    // only the p+1 control points of the span (not the coordinates of the whole patch)
    for ( int k = 0; k < ders.giveNumberOfColumns(); k++ ) {
        const auto &vertexCoords = cellgeo.giveVertexCoordinates( mask.at( k + 1 ) );
        double vertexCoord_1     = vertexCoords[0];
        double vertexCoord_2     = vertexCoords[1];

        Rxi += ders( 1, k ) * vertexCoord_1;
        Yxi += ders( 1, k ) * vertexCoord_2;

        Rxixi += ders( 2, k ) * vertexCoord_1;
        Yxixi += ders( 2, k ) * vertexCoord_2;

        R += vertexCoord_1 * ders( 0, k ); // compute radial coordinate at the gauss point
    }
    checkRadius( R, gp );

    // The reference surface X(xi, theta) = R(xi) e_r(theta) + Y(xi) e_y has the tangents X_,xi and
    // X_,theta = R e_theta, the dual vectors A^xi = g^11 X_,xi and A^theta = e_theta / R, and the nonzero
    // Christoffel symbols Gamma^xi_xixi, Gamma^xi_thetatheta and Gamma^theta_xitheta = R_,xi / R. For the
    // current position p = r e_r + y e_y (e_r,theta = e_theta, e_r,thetatheta = -e_r), grad F = G p is the
    // covariant Hessian (p_,ab - Gamma^c_ab p_,c) (x) A^a (x) A^b, (x) the tensor product. At theta = 0, e_r, e_y
    // and e_theta are the global x, y and z.

    // Inverse (contravariant) metric coefficient g^{11} = 1 / |X_,xi|^2
    double gInv = 1. / ( Rxi * Rxi + Yxi * Yxi ); // Same as 1/J^2

    // Christoffel symbols Gamma^xi_xixi = g^11 X_,xi . X_,xixi and Gamma^xi_thetatheta = -g^11 R R_,xi
    double Gamma_XiXiXi       = gInv * ( Rxi * Rxixi + Yxi * Yxixi );
    double Gamma_XiThetaTheta = -gInv * ( R * Rxi );

    // Base vectors e_r (Ir), e_y (Iy), e_theta (Is) and the dual vector A^xi (Exi_contra)
    FloatArrayF<3> IrFA = { 1., 0., 0. }, IyFA = { 0., 1., 0. }, IsFA = { 0., 0., 1. };
    FloatArrayF<3> Exi_contraFA = gInv * ( Rxi * IrFA + Yxi * IyFA );

    Tensor1_3d Ir( IrFA ), Iy( IyFA ), Is( IsFA ), Exi_contra( Exi_contraFA );

    // Create part proportional to r-coordinate (x_1)
    // xi xi: e_r (x) A^xi (x) A^xi; xi theta: e_theta (x) e_theta (x) A^xi + e_theta (x) A^xi (x) e_theta;
    // theta theta: e_r (x) e_theta (x) e_theta
    Tensor3_3d dF_r1, dF_r2, dF_r3;
    dF_r1( i_3, j_3, k_3 ) = Ir( i_3 ) * Exi_contra( j_3 ) * Exi_contra( k_3 );
    dF_r2( i_3, j_3, k_3 ) = Is( i_3 ) * Is( j_3 ) * Exi_contra( k_3 ) + Is( i_3 ) * Exi_contra( j_3 )* Is( k_3 );
    dF_r3( i_3, j_3, k_3 ) = Ir( i_3 ) * Is( j_3 ) * Is( k_3 );
    // change to voigt
    FloatArray dF_r1v = dF_r1.to_voigt_form_27();
    FloatArray dF_r2v = dF_r2.to_voigt_form_27();
    FloatArray dF_r3v = dF_r3.to_voigt_form_27();


    // Create part proportional to y-coordinate (x_2)
    // xi xi: e_y (x) A^xi (x) A^xi; theta theta: e_y (x) e_theta (x) e_theta
    Tensor3_3d dF_y1, dF_y2;
    dF_y1( i_3, j_3, k_3 ) = Iy( i_3 ) * Exi_contra( j_3 ) * Exi_contra( k_3 );
    dF_y2( i_3, j_3, k_3 ) = Iy( i_3 ) * Is( j_3 ) * Is( k_3 );
    // change to voigt
    FloatArray dF_y1v = dF_y1.to_voigt_form_27();
    FloatArray dF_y2v = dF_y2.to_voigt_form_27();

    int nno = ders.giveNumberOfColumns(); // nonzero basis functions of the span
    answer.resize( 27, nno * 2 );
    answer.zero();

    for ( int i = 1; i <= nno; i++ ) {
        // coefficients of the radial dof (column 2i-1), with N, N', N'' = ders( 0..2, i - 1 ) of basis function
        // i; each A^theta = e_theta / R gives a factor 1 / R:
        // cR1 = N'' - Gamma^xi_xixi N', cR2 = (N' - Gamma^theta_xitheta N) / R = (N' R - R_,xi N) / R^2,
        // cR3 = (-N - Gamma^xi_thetatheta N') / R^2 (the -N from e_r,thetatheta = -e_r)
        double cR1 = ders( 2, i - 1 ) - ders( 1, i - 1 ) * Gamma_XiXiXi;
        double cR2 = (ders( 1, i - 1 )*R - Rxi*ders( 0, i - 1 ))/R/R ;
        double cR3 = -( ders( 1, i - 1 ) * Gamma_XiThetaTheta + ders(0, i - 1 ) ) / R / R;

        // coefficients of the axial dof (column 2i): e_y does not turn with theta
        double cY1 = cR1;
        double cY2 = -( ders( 1, i - 1 ) * Gamma_XiThetaTheta ) / R / R;

        FloatArray col1 = cR1 * dF_r1v + cR2 * dF_r2v + cR3 * dF_r3v;
        FloatArray col2 = cY1 * dF_y1v + cY2 * dF_y2v;

        answer.setColumn( col1,2*i-1);
        answer.setColumn( col2, 2 * i);
    }
}


double AxisymSurfaceElementEvaluator ::computeVolumeAround( GaussPoint *gp )
{
    // area per radian: the length of the reference curve times the reference radius R
    return Space2dStructuralSurfaceElementEvaluator::computeVolumeAround( gp ) * this->computeRadiusAt( gp );
}


} // end namespace oofem
