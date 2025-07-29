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
//#include "sm/CrossSections/structuralcrosssection.h"
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

namespace oofem {
void Space2dStructuralSurfaceElementEvaluator ::computeNMatrixAt( FloatMatrix &answer, GaussPoint *gp )
{
    FloatArray N;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    answer.beNMatrixOf( N, 2 );

    // Transform to global coordinates!!
}


void Space2dStructuralSurfaceElementEvaluator::computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    // Returns the [ 4 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
    // evaluated at gp.
    /// @todo not checked if correct

    FloatMatrix dNdx;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );



    // Transform to global coordinates!!
    int nRows = dNdx.giveNumberOfRows();
    answer.resize( 9, nRows * 2 );
    answer.zero();

    //for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
    //    answer.at( 1, 2 * i - 1 ) = dNdx.at( i, 1 ); // du/dx -1
    //    answer.at( 2, 2 * i - 0 ) = dNdx.at( i, 2 ); // dv/dy -2
    //    answer.at( 4, 2 * i - 1 ) = dNdx.at( i, 2 ); // du/dy -6
    //    answer.at( 5, 2 * i - 0 ) = dNdx.at( i, 1 ); // dv/dx -9
    //}

    // assume 
    for ( int i = 1; i <= nRows; i++ ) {
        answer.at( 1, 2 * i - 1 ) = dNdx.at( i, 1 ); // du/dx -1
        answer.at( 2, 2 * i - 0 ) = dNdx.at( i, 2 ); // dv/dy -2
        answer.at( 6, 2 * i - 1 ) = dNdx.at( i, 2 ); // du/dy -6
        answer.at( 9, 2 * i - 0 ) = dNdx.at( i, 1 ); // dv/dx -9
    }
}

void Space2dStructuralSurfaceElementEvaluator::computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    // Returns the [ 4 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
    // evaluated at gp.
    /// @todo not checked if correct

    FloatMatrix d2Ndx2;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evald2Ndx2( d2Ndx2, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );


    // Transform to global coordinates!!
    int nRows = d2Ndx2.giveNumberOfRows();
    answer.resize( 27, nRows * 2 );
    answer.zero();

    // assume
    for ( int i = 1; i <= nRows; i++ ) {
        answer.at( 1, 2 * i - 1 ) = d2Ndx2.at( i, 1 ); // du/dx -1
        answer.at( 2, 2 * i - 0 ) = d2Ndx2.at( i, 3 ); // dv/dy -2
        answer.at( 6, 2 * i - 1 ) = d2Ndx2.at( i, 3 ); // du/dy -6
        answer.at( 9, 2 * i - 0 ) = d2Ndx2.at( i, 1 ); // dv/dx -9

        answer.at( 1 + 9, 2 * i - 1 ) = d2Ndx2.at( i, 3 ); // du/dx -1
        answer.at( 2 + 9, 2 * i - 0 ) = d2Ndx2.at( i, 2 ); // dv/dy -2
        answer.at( 6 + 9, 2 * i - 1 ) = d2Ndx2.at( i, 2 ); // du/dy -6
        answer.at( 9 + 9, 2 * i - 0 ) = d2Ndx2.at( i, 3 ); // dv/dx -9
    }
}


double Space2dStructuralSurfaceElementEvaluator ::computeVolumeAround( GaussPoint *gp )
{
    double determinant, weight, thickness, volume;
    determinant = fabs( this->giveElement()->giveInterpolation()->giveTransformationJacobian( gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( this->giveElement(),
            gp->giveIntegrationRule()->giveKnotSpan() ) ) );
    weight      = gp->giveWeight();
    volume      = determinant * weight;

    return volume;
}


void Space2dStructuralSurfaceElementEvaluator ::initializeFrom( InputRecord &ir )
{
    StructuralSurfaceElementEvaluator::initializeFrom( ir );
}

int Space2dStructuralSurfaceElementEvaluator ::giveIntegrationElementLocalCodeNumbers( IntArray &answer, Element *elem,
    IntegrationRule *ie )
{
    int nsd;
    IntArray mask, nodeDofIDMask, nodalArray;
    int dofmandof;

    // get number of dofs in node
    elem->giveDofManDofIDMask( 1, nodeDofIDMask );
    dofmandof = nodeDofIDMask.giveSize();

    nsd = elem->giveInterpolation()->giveNsd();

    if ( nsd == 1 ) { // Change compared to the general funciton!!
        nsd = 2;
    }

    // first evaluate nonzero basis function mask
    if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
        IGAIntegrationElement *ee = static_cast<IGAIntegrationElement *>( ie );
        elem->giveInterpolation()->giveKnotSpanBasisFuncMask( *ee->giveKnotSpan(), mask );
        // loop over nonzero shape functions and assemble localization array
        answer.clear();
        for ( int i = 1; i <= mask.giveSize(); i++ ) {
            nodalArray.resize( nodeDofIDMask.giveSize() );
            for ( int j = 1; j <= nsd; j++ ) { 
                nodalArray.at( j ) = dofmandof * ( mask.at( i ) - 1 ) + j;
            }

            answer.followedBy( nodalArray );
        }

        return 1;
    } else {
        return 0;
    }
}


//int Space2dStructuralSurfaceElementEvaluator ::giveIntegrationElementLocalCodeNumbers( IntArray &answer, Element *elem,
//    FloatArray &lcoords )
//{
//    int nsd;
//    IntArray mask, nodeDofIDMask, nodalArray;
//    int dofmandof;
//
//    // get number of dofs in node
//    elem->giveDofManDofIDMask( 1, nodeDofIDMask );
//    dofmandof = nodeDofIDMask.giveSize();
//
//    nsd = elem->giveInterpolation()->giveNsd();
//
//    if ( nsd == 1 ) { // Change compared to the general funciton!!
//        nsd = 2;
//    }
//
//    // first evaluate nonzero basis function mask
//    if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
//        IGAIntegrationElement *ee = static_cast<IGAIntegrationElement *>( ie );
//
//
//        BSpline2dLineInterpolation *ee = static_cast<BSpline2dLineInterpolation *>( elem->giveInterpolation() );
//
//        elem->giveInterpolation().findSpan( numberOfControlPoints[i], degree[i], lcoords[i], knotVector[i] );
//        //{
//        //    span = *gw.knotSpan;
//        //}
//        //else
//        //{
//        //    for ( int i = 0; i < nsd; i++ ) {
//        //        span[i] = this->findSpan( numberOfControlPoints[i], degree[i], lcoords[i], knotVector[i] );
//        //    }
//        //}
//
//        elem->giveInterpolation()->giveKnotSpanBasisFuncMask( *ee->giveKnotSpan(), mask );
//        // loop over nonzero shape functions and assemble localization array
//        answer.clear();
//        for ( int i = 1; i <= mask.giveSize(); i++ ) {
//            nodalArray.resize( nodeDofIDMask.giveSize() );
//            for ( int j = 1; j <= nsd; j++ ) {
//                nodalArray.at( j ) = dofmandof * ( mask.at( i ) - 1 ) + j;
//            }
//
//            answer.followedBy( nodalArray );
//        }
//
//        return 1;
//    } else {
//        return 0;
//    }
//}

// void Space2dStructuralSurfaceElementEvaluator ::computeBMatrixAt( FloatMatrix &answer, GaussPoint *gp )
//{
//     FloatMatrix d;
//
//     FEInterpolation *interp = gp->giveElement()->giveInterpolation();
//     // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
//     interp->evaldNdx( d, gp->giveNaturalCoordinates(),
//         FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
//
//     answer.resize( 4, d.giveNumberOfRows() * 2 );
//     answer.zero();
//
//     // Transform to global coordinates!!
//
//     for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
//         answer.at( 1, i * 2 - 1 ) = d.at( i, 1 );
//         answer.at( 2, i * 2 - 0 ) = d.at( i, 2 );
//
//         answer.at( 4, 2 * i - 1 ) = d.at( i, 2 );
//         answer.at( 4, 2 * i - 0 ) = d.at( i, 1 );
//     }
// }


/////////////////////
// Axisym formulation

void AxisymSurfaceElementEvaluator::computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    // Returns the [ 4 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
    // evaluated at gp.
    /// @todo not checked if correct

    FloatMatrix dNdx;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    //BSplineInterpolation *interp = static_cast<BSplineInterpolation *>( gp->giveElement()->giveInterpolation() );

    // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    // Transform to global coordinates!!
    int nRows = dNdx.giveNumberOfRows();
    answer.resize( 9, nRows * 2 );
    answer.zero();

    // assume
    for ( int i = 1; i <= nRows; i++ ) {
        answer.at( 1, 2 * i - 1 ) = dNdx.at( i, 1 ); // du/dx -1
        answer.at( 2, 2 * i - 0 ) = dNdx.at( i, 2 ); // dv/dy -2
        answer.at( 6, 2 * i - 1 ) = dNdx.at( i, 2 ); // du/dy -6
        answer.at( 9, 2 * i - 0 ) = dNdx.at( i, 1 ); // dv/dx -9
    }


    // compute the axisymmetrix part 
    FloatArray n;
    //interp->evalN( n, gp->giveNaturalCoordinates(), 
    //    FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    // Check if these two options are the same

    FloatMatrix nmat;
    // return index where nonzero basis functions start, 1-based
    int indexStart = interp->evalDerivatives( 0, nmat, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    FloatMatrix nmatT;
    nmatT.beTranspositionOf( nmat );
    n.beColumnOf( nmatT, 1 );


    //// Get nodal coordinates
    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );
    FloatArray nodesCoords;
    elem->giveNodalCoordinates( nodesCoords );

    //// also works
    //double r = 0., x;
    //for ( int i = 1; i <= n.giveSize(); i++ ) {
    //    int indCords = indexStart + i - 1;
    //    x = nodesCoords.at( 2 * indCords - 1 ); // get radial (x) coordinate of the nodes
    //    r += x * n.at( i ); // compute radial coordinate at the gauss point
    //}

    FloatArray gcoords;
    interp->local2global( gcoords, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    double r = gcoords.at( 1 ); 

    for ( int i = 0; i < n.giveSize(); i++ ) {
        answer.at( 3, 2 * i + 1 ) = n.at( i + 1 ) / r; // dUs/ds
    }
}

void AxisymSurfaceElementEvaluator::computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    FloatMatrix ders;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    //BSplineInterpolation *interpBS = static_cast<BSplineInterpolation *>( gp->giveElement()->giveInterpolation() );

    //// this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    //interp->evald2Ndx2( d2Ndx2, gp->giveNaturalCoordinates(),
    //    FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    // returns index where nonzero basis functions start (1-based)
    int indexStart = interp->evalDerivatives( 2, ders, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    ///////////////
    //// Get nodal coordinates
    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );
    FloatArray nodesCoords;
    elem->giveNodalCoordinates( nodesCoords );

    double Rxi = 0., Yxi = 0., Rxixi = 0., Yxixi = 0.;
    double R = 0.; // undeformed radius (X1)

    for ( int k = 0; k < ders.giveNumberOfColumns(); k++ ) {
        int indCords = indexStart + k;
        double vertexCoord_1 = nodesCoords.at( 2 * indCords - 1 );
        double vertexCoord_2  = nodesCoords.at( 2 * indCords );

        Rxi += ders( 1, k ) * vertexCoord_1;
        Yxi += ders( 1, k ) * vertexCoord_2;

        Rxixi += ders( 2, k ) * vertexCoord_1;
        Yxixi += ders( 2, k ) * vertexCoord_2;

        R += vertexCoord_1 * ders( 0, k ); // compute radial coordinate at the gauss point
    }

    //FloatArray gcoords;
    //interp->local2global( gcoords, gp->giveNaturalCoordinates(),
    //FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    //double r = gcoords.at( 1 ); 

    // Compute covariant metric coefficient
    double G_11 = 1. / ( Rxi * Rxi + Yxi * Yxi ); // Same as 1/J^2

    // Christoffel symbol
    double Gamma_XiXiXi       = G_11 * ( Rxi * Rxixi + Yxi * Yxixi );
    double Gamma_XiThetaTheta = -G_11 * ( R * Rxi );

    // Create base vectors and covectors
    FloatArrayF<3> IrFA = { 1., 0., 0. }, IyFA = { 0., 1., 0. }, IsFA = { 0., 0., 1. };
    FloatArrayF<3> Exi_contraFA = G_11 * ( Rxi * IrFA + Yxi * IyFA );

    Tensor1_3d Ir( IrFA ), Iy( IyFA ), Is( IsFA ), Exi_contra( Exi_contraFA );

    // Create part proportional to r-coordinate (x_1)
    Tensor3_3d dF_r1, dF_r2, dF_r3;
    dF_r1( i_3, j_3, k_3 ) = Ir( i_3 ) * Exi_contra( j_3 ) * Exi_contra( k_3 );
    dF_r2( i_3, j_3, k_3 ) = Is( i_3 ) * Is( j_3 ) * Exi_contra( k_3 ) + Is( i_3 ) * Exi_contra( j_3 )* Is( k_3 );
    dF_r3( i_3, j_3, k_3 ) = Ir( i_3 ) * Is( j_3 ) * Is( k_3 );
    // change to voigt
    FloatArray dF_r1v = dF_r1.to_voigt_form_27();
    FloatArray dF_r2v = dF_r2.to_voigt_form_27();
    FloatArray dF_r3v = dF_r3.to_voigt_form_27();


    // Create part proportional to y-coordinate (x_2)
    Tensor3_3d dF_y1, dF_y2;
    dF_y1( i_3, j_3, k_3 ) = Iy( i_3 ) * Exi_contra( j_3 ) * Exi_contra( k_3 );
    dF_y2( i_3, j_3, k_3 ) = Iy( i_3 ) * Is( j_3 ) * Is( k_3 );
    // change to voigt
    FloatArray dF_y1v = dF_y1.to_voigt_form_27();
    FloatArray dF_y2v = dF_y2.to_voigt_form_27();

    // Transform to global coordinates!!
    int nRows = ders.giveNumberOfColumns();
    answer.resize( 27, nRows * 2 );
    answer.zero();

    // assume
    for ( int i = 1; i <= nRows; i++ ) {
        double cF1 = ders( 2, i - 1 ) - ders( 1, i - 1 ) * Gamma_XiXiXi;
        double cF2 = (ders( 1, i - 1 )*R - Rxi*ders( 0, i - 1 ))/R/R ;
        double cF3 = -( ders( 1, i - 1 ) * Gamma_XiThetaTheta + ders(0, i - 1 ) ) / R / R;   

        double dF1 = cF1;
        double dF2 = -( ders( 1, i - 1 ) * Gamma_XiThetaTheta ) / R / R;   

        FloatArray col1 = cF1 * dF_r1v + cF2 * dF_r2v + cF3 * dF_r3v;
        FloatArray col2 = dF1 * dF_y1v + dF2 * dF_y2v;

        answer.setColumn( col1,2*i-1);
        answer.setColumn( col2, 2 * i);
    }
}


double AxisymSurfaceElementEvaluator ::computeVolumeAround( GaussPoint *gp )
{
    double determinant, weight, thickness, volume;
    determinant = fabs( this->giveElement()->giveInterpolation()->giveTransformationJacobian( gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( this->giveElement(),
            gp->giveIntegrationRule()->giveKnotSpan() ) ) );
    weight      = gp->giveWeight();
    volume      = determinant * weight;

    // Compute R
    FloatMatrix ders;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    // returns index where nonzero basis functions start (1-based)
    int indexStart = interp->evalDerivatives( 0, ders, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    //// Get nodal coordinates
    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );
    FloatArray nodesCoords;
    elem->giveNodalCoordinates( nodesCoords );

    //// also valid approach
    //double R = 0.; // undeformed radius (X1)
    //for ( int k = 0; k < ders.giveNumberOfColumns(); k++ ) {
    //    int indCords              = indexStart + k;
    //    double vertexCoord_1 = nodesCoords.at( 2 * indCords - 1 );
    //    R += vertexCoord_1 * ders( 0, k ); // compute radial coordinate at the gauss point
    //}

    FloatArray gcoords;
    interp->local2global( gcoords, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    double R = gcoords.at( 1 ); 

    return volume*R;
}


void AxisymSurfaceElementEvaluator ::initializeFrom( InputRecord &ir )
{
    Space2dStructuralSurfaceElementEvaluator::initializeFrom( ir );
}



} // end namespace oofem
