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


void IGAClampBoundaryCondition ::assemble( SparseMtrx &answer, TimeStep *tStep,
    CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s,
    double scale,
    void *lock )
{
    // the clamp term is linear in the unknowns, so all stiffness types get the same matrix
    if ( type != TangentStiffnessMatrix && type != SecantStiffnessMatrix && type != ElasticStiffnessMatrix ) {
        return;
    }

    auto el                 = this->giveDomain()->giveElement( this->elementNumber );
    FEInterpolation *interp = el->giveInterpolation();
    // normal of the reference curve in the selected knot span
    FloatArray normal;
    auto gw = FEIIGAElementGeometryWrapper( el, &this->BcKnotSpan );
    interp->boundaryEvalNormal( normal, 1, this->localPosition, gw );
    normal.resize( 2 ); // drop the out-of-plane component

    // rows: basis functions and their first derivatives w.r.t. xi (nonzero functions of the span)
    FloatMatrix ders, dersT;
    interp->evalDerivatives( 1, ders, this->localPosition, gw );

    // du/dxi = dNdxi_mat * u_span (2 x 2*nnz)
    FloatMatrix dNdxi_mat;
    dNdxi_mat.resize( 2, 2 * ders.giveNumberOfColumns() );
    for ( int i = 1; i <= ders.giveNumberOfColumns(); i++ ) {
        dNdxi_mat.at( 1, 2 * i - 1 ) = ders.at( 2, i );
        dNdxi_mat.at( 2, 2 * i )     = ders.at( 2, i );
    }

    // c = dNdxi_mat^T n, so that n . du/dxi = c . u_span
    FloatArray dNdxi;
    dNdxi.beTProductOf( dNdxi_mat, normal );

    // penalty matrix p c c^T
    FloatMatrix dNdNt;
    dNdNt.beDyadicProductOf( dNdxi, dNdxi );
    dNdNt.times( this->penalty * scale );


    IGAElement *elem = static_cast<IGAElement *>( el );
    // element-local code numbers of the dofs of the current knot span
    IntArray lc;
    this->giveLocalCodeNumbers( lc, elem, this->BcKnotSpan );


    // global equation numbers of the dofs of the clamped knot span only
    IntArray lcGlob_r, lcGlob_c, loc_r( lc.giveSize() ), loc_c( lc.giveSize() );
    elem->giveLocationArray( lcGlob_r, r_s );
    elem->giveLocationArray( lcGlob_c, c_s );
    for ( int i = 1; i <= lc.giveSize(); i++ ) {
        loc_r.at( i ) = lc.at( i ) ? lcGlob_r.at( lc.at( i ) ) : 0;
        loc_c.at( i ) = lc.at( i ) ? lcGlob_c.at( lc.at( i ) ) : 0;
    }
#ifdef _OPENMP
    if ( lock ) {
        omp_set_lock( static_cast<omp_lock_t *>( lock ) );
    }
#endif
    answer.assemble( loc_r, loc_c, dNdNt );
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

    if ( type == InternalForcesVector ) {

        auto el                 = this->giveDomain()->giveElement( this->elementNumber );
        FEInterpolation *interp = el->giveInterpolation();

        // normal of the reference curve in the selected knot span
        FloatArray normal;
        auto gw = FEIIGAElementGeometryWrapper( el, &this->BcKnotSpan );
        interp->boundaryEvalNormal( normal, 1, this->localPosition, gw );
        normal.resize( 2 ); // drop the out-of-plane component

        FloatMatrix ders, dersT;
        // rows: basis functions and their first derivatives w.r.t. xi (nonzero functions of the span)
        interp->evalDerivatives( 1, ders, this->localPosition, gw );

        // du/dxi = dNdxi_mat * u_span (2 x 2*nnz)
        FloatMatrix dNdxi_mat;
        dNdxi_mat.resize( 2, 2 * ders.giveNumberOfColumns() );
        for ( int i = 1; i <= ders.giveNumberOfColumns(); i++ ) {
            dNdxi_mat.at( 1, 2 * i - 1 ) = ders.at( 2, i );
            dNdxi_mat.at( 2, 2 * i )     = ders.at( 2, i );
        }

        // c = dNdxi_mat^T n, so that n . du/dxi = c . u_span
        FloatArray dNdxi;
        dNdxi.beTProductOf( dNdxi_mat, normal );

        // total displacements of all control points of the patch
        FloatArray u;
        IGAElement *elem = static_cast<IGAElement *>( el );
        elem->computeVectorOf( VM_Total, tStep, u );

        // element-local code numbers of the dofs of the current knot span
        IntArray lc,lcGlob;
        FloatArray ur, xr;
        this->giveLocalCodeNumbers( lc, elem, this->BcKnotSpan );

        // displacements of the current knot span
        ur.resize( dNdxi.giveSize() );
        for ( int i = 1; i <= lc.giveSize(); i++ ) {
            if ( lc.at( i ) ) {
                ur.at( i ) = u.at( lc.at( i ) );
            }
        }

        // internal forces p c (c . u_span)
        FloatArray prod2 = dNdxi.dotProduct( ur ) * dNdxi;
        prod2.times( this->penalty );

        // global equation numbers of the dofs of the clamped knot span only
        IntArray loc( lc.giveSize() );
        elem->giveLocationArray( lcGlob, s );
        for ( int i = 1; i <= lc.giveSize(); i++ ) {
            loc.at( i ) = lc.at( i ) ? lcGlob.at( lc.at( i ) ) : 0;
        }

#ifdef _OPENMP
        if ( lock ) {
            omp_set_lock( static_cast<omp_lock_t *>( lock ) );
        }
#endif
        answer.assemble( prod2, loc );
#ifdef _OPENMP
        if ( lock ) {
            omp_unset_lock( static_cast<omp_lock_t *>( lock ) );
        }
#endif
    }
}

int IGAClampBoundaryCondition ::giveLocalCodeNumbers( IntArray &answer, Element *elem, IntArray &span )
{
    int nsd;
    IntArray mask, nodeDofIDMask, nodalArray;
    int dofmandof;

    // get number of dofs in node
    elem->giveDofManDofIDMask( 1, nodeDofIDMask );
    dofmandof = nodeDofIDMask.giveSize();

    nsd = elem->giveInterpolation()->giveNsd();

    // unlike StructuralElementEvaluator::giveIntegrationElementLocalCodeNumbers: the planar curve has
    // one parametric dimension but two displacement dofs (u, v) per control point
    if ( nsd == 1 ) {
        nsd = 2;
    }

    // first evaluate nonzero basis function mask
    if ( elem->giveInterpolation()->hasSubPatchFormulation() ) {
        elem->giveInterpolation()->giveKnotSpanBasisFuncMask( span, mask );
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

} // namespace oofem
