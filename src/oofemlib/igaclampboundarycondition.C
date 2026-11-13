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


#include "igaclampboundarycondition.h"
#include "set.h"
#include "domain.h"
#include "node.h"
#include "masterdof.h"
#include "floatmatrix.h"
#include "unknownnumberingscheme.h"
#include "sparsemtrx.h"
#include "classfactory.h"
#include "iga/feibspline.h"
#include "iga/iga.h"

#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace oofem {
REGISTER_BoundaryCondition( IGAClampBoundaryCondition );


IGAClampBoundaryCondition ::IGAClampBoundaryCondition( int n, Domain *d ) :
    ActiveBoundaryCondition( n, d )
{
}


void IGAClampBoundaryCondition ::initializeFrom( InputRecord &ir )
{
    ActiveBoundaryCondition ::initializeFrom( ir );
    IR_GIVE_FIELD( ir, this->elementNumber, _IFT_IGAClampBoundaryCondition_ElementNumber );
    double lp;
    IR_GIVE_FIELD( ir, lp, _IFT_IGAClampBoundaryCondition_LocalPosition );
    this->localPosition.resize( 1 );
    this->localPosition.at(1) = lp;

    int ks;
    IR_GIVE_FIELD( ir, ks, _IFT_IGAClampBoundaryCondition_BcKnotSpan );
    this->BcKnotSpan.resize( 1 );
    this->BcKnotSpan.at( 1 ) = ks;

    IR_GIVE_OPTIONAL_FIELD( ir, this->penalty, _IFT_IGAClampBoundaryCondition_p );
    if ( !( this->penalty > 0. ) || !std::isfinite( this->penalty ) ) {
        throw ValueInputException( ir, _IFT_IGAClampBoundaryCondition_p, "must be positive and finite" );
    }
}


void IGAClampBoundaryCondition ::checkClampLocation( Domain *d, int elementNumber, double localPosition, int knotSpan )
{
    if ( elementNumber < 1 || elementNumber > d->giveNumberOfElements() ) {
        OOFEM_ERROR( "numel %d is outside the element range", elementNumber );
    }
    auto elem   = dynamic_cast<IGAElement *>( d->giveElement( elementNumber ) );
    auto interp = elem ? dynamic_cast<BSpline2dLineInterpolation *>( elem->giveInterpolation() ) : nullptr;
    if ( !interp ) {
        OOFEM_ERROR( "numel %d must reference a planar IGA curve element", elementNumber );
    }

    const auto &knots = interp->giveKnotVector()[0];
    const int degree  = interp->giveInterpolationOrder();
    const int ncp     = interp->giveNumberOfControlPoints( 1 );
    if ( degree < 1 ) {
        OOFEM_ERROR( "numel %d must be at least of degree 1 to be clamped", elementNumber );
    }
    if ( knotSpan < degree || knotSpan >= ncp ) {
        OOFEM_ERROR( "bcknotspan %d is outside the zero-based knot-index range [%d, %d]", knotSpan, degree, ncp - 1 );
    }
    if ( !( knots[knotSpan] < knots[knotSpan + 1] ) ) {
        OOFEM_ERROR( "bcknotspan %d must identify a nonzero knot interval", knotSpan );
    }
    // Include both ends so an internal knot can be approached from either selected span (NaN fails both comparisons).
    if ( !std::isfinite( localPosition ) || localPosition < knots[knotSpan] || localPosition > knots[knotSpan + 1] ) {
        OOFEM_ERROR( "locpos %g is outside bcknotspan %d [%g, %g]", localPosition, knotSpan, knots[knotSpan], knots[knotSpan + 1] );
    }
}


void IGAClampBoundaryCondition ::postInitialize()
{
    ActiveBoundaryCondition ::postInitialize();
    checkClampLocation( this->giveDomain(), this->elementNumber, this->localPosition.at( 1 ), this->BcKnotSpan.at( 1 ) );
}


void IGAClampBoundaryCondition ::computeClampVector( FloatArray &answer, Element *elem, const FloatArray &localPosition, const IntArray &knotSpan )
{
    auto interp = static_cast<const BSpline2dLineInterpolation *>( elem->giveInterpolation() ); // checked in checkClampLocation
    FEIIGAElementGeometryWrapper gw( elem, &knotSpan );

    // normal of the reference curve in the knot span; rows of ders: basis functions and their first
    // derivatives w.r.t. xi (nonzero functions of the span)
    FloatArray normal;
    FloatMatrix ders;
    interp->evalCurveNormal( normal, localPosition, gw );
    interp->evalDerivatives( 1, ders, localPosition, gw );

    // du/dxi = dNdxi_mat * u_span (2 x 2*nnz), c = dNdxi_mat^T n
    FloatMatrix dNdxi_mat;
    dNdxi_mat.resize( 2, 2 * ders.giveNumberOfColumns() );
    for ( int i = 1; i <= ders.giveNumberOfColumns(); i++ ) {
        dNdxi_mat.at( 1, 2 * i - 1 ) = ders.at( 2, i );
        dNdxi_mat.at( 2, 2 * i )     = ders.at( 2, i );
    }
    answer.beTProductOf( dNdxi_mat, normal );
}


void IGAClampBoundaryCondition ::giveSpanCodeNumbers( IntArray &answer, Element *elem, const IntArray &knotSpan )
{
    // all nodes have the dofs of node 1, u and v first (the IGA surface elements have just these two)
    IntArray mask, nodeDofIDMask;
    elem->giveDofManDofIDMask( 1, nodeDofIDMask );
    int dofmandof = nodeDofIDMask.giveSize();

    elem->giveInterpolation()->giveKnotSpanBasisFuncMask( knotSpan, mask );
    answer.resize( 2 * mask.giveSize() );
    for ( int i = 1; i <= mask.giveSize(); i++ ) {
        answer.at( 2 * i - 1 ) = dofmandof * ( mask.at( i ) - 1 ) + 1;
        answer.at( 2 * i )     = dofmandof * ( mask.at( i ) - 1 ) + 2;
    }
}


void IGAClampBoundaryCondition ::giveSpanLocationArray( IntArray &answer, Element *elem, const IntArray &knotSpan, const UnknownNumberingScheme &s )
{
    IntArray lc, elemLoc;
    giveSpanCodeNumbers( lc, elem, knotSpan );
    elem->giveLocationArray( elemLoc, s );
    answer.resize( lc.giveSize() );
    for ( int i = 1; i <= lc.giveSize(); i++ ) {
        answer.at( i ) = elemLoc.at( lc.at( i ) );
    }
}


void IGAClampBoundaryCondition ::giveSpanDisplacements( FloatArray &answer, Element *elem, const IntArray &knotSpan, TimeStep *tStep )
{
    IntArray lc;
    FloatArray u;
    giveSpanCodeNumbers( lc, elem, knotSpan );
    elem->computeVectorOf( VM_Total, tStep, u );
    answer.resize( lc.giveSize() );
    for ( int i = 1; i <= lc.giveSize(); i++ ) {
        answer.at( i ) = u.at( lc.at( i ) );
    }
}


void IGAClampBoundaryCondition ::assemble( SparseMtrx &answer, TimeStep *tStep,
    CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s,
    double scale,
    void *lock )
{
    // the clamp term is linear in the unknowns, so all stiffness types get the same matrix
    if ( type != TangentStiffnessMatrix && type != SecantStiffnessMatrix && type != ElasticStiffnessMatrix ) {
        return;
    }

    auto el = this->giveDomain()->giveElement( this->elementNumber );
    FloatArray c;
    IntArray loc_r, loc_c;
    computeClampVector( c, el, this->localPosition, this->BcKnotSpan );
    giveSpanLocationArray( loc_r, el, this->BcKnotSpan, r_s );
    giveSpanLocationArray( loc_c, el, this->BcKnotSpan, c_s );

    // penalty matrix p c c^T
    FloatMatrix ccT;
    ccT.beDyadicProductOf( c, c );
    ccT.times( this->penalty * scale );

#ifdef _OPENMP
    if ( lock ) {
        omp_set_lock( static_cast<omp_lock_t *>( lock ) );
    }
#endif
    answer.assemble( loc_r, loc_c, ccT );
#ifdef _OPENMP
    if ( lock ) {
        omp_unset_lock( static_cast<omp_lock_t *>( lock ) );
    }
#endif
}


void IGAClampBoundaryCondition ::assembleVector( FloatArray &answer, TimeStep *tStep,
    CharType type, ValueModeType mode,
    const UnknownNumberingScheme &s,
    FloatArray *eNorms,
    void *lock )
{
    if ( type != InternalForcesVector ) {
        return;
    }

    auto el = this->giveDomain()->giveElement( this->elementNumber );
    FloatArray c, u;
    IntArray loc;
    computeClampVector( c, el, this->localPosition, this->BcKnotSpan );
    giveSpanDisplacements( u, el, this->BcKnotSpan, tStep );
    giveSpanLocationArray( loc, el, this->BcKnotSpan, s );

    // internal forces p c (c . u_span)
    FloatArray f = c.dotProduct( u ) * c;
    f.times( this->penalty );

#ifdef _OPENMP
    if ( lock ) {
        omp_set_lock( static_cast<omp_lock_t *>( lock ) );
    }
#endif
    answer.assemble( f, loc );
#ifdef _OPENMP
    if ( lock ) {
        omp_unset_lock( static_cast<omp_lock_t *>( lock ) );
    }
#endif
}

} // namespace oofem
