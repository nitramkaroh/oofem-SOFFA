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


#include "igaclampboundaryconditionlagrange.h"
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
#include "igaclampboundarycondition.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace oofem {
REGISTER_BoundaryCondition( IGAClampBoundaryConditionLagrange );


IGAClampBoundaryConditionLagrange ::IGAClampBoundaryConditionLagrange( int n, Domain *d ) :
    ActiveBoundaryCondition( n, d )
{
}


void IGAClampBoundaryConditionLagrange ::initializeFrom( InputRecord &ir )
{
    ActiveBoundaryCondition ::initializeFrom( ir );
    IR_GIVE_FIELD( ir, this->elementNumber, _IFT_IGAClampBoundaryConditionLagrange_ElementNumber );
    double lp;
    IR_GIVE_FIELD( ir, lp, _IFT_IGAClampBoundaryConditionLagrange_LocalPosition );
    this->localPosition.resize( 1 );
    this->localPosition.at( 1 ) = lp;

    int ks;
    IR_GIVE_FIELD( ir, ks, _IFT_IGAClampBoundaryConditionLagrange_BcKnotSpan );
    this->BcKnotSpan.resize( 1 );
    this->BcKnotSpan.at( 1 ) = ks;


    // internal dof manager with one dof (new dof id) for the Lagrange multiplier
    for ( int pos = 0; pos < 1; ++pos ) {
        lmdm =  std ::unique_ptr<Node>( new Node( 0, domain ) ) ;
        lmdm->appendDof( new MasterDof( this->lmdm.get(), (DofIDItem)( this->giveDomain()->giveNextFreeDofID() ) ) );
    }
}


void IGAClampBoundaryConditionLagrange ::postInitialize()
{
    ActiveBoundaryCondition ::postInitialize();
    IGAClampBoundaryCondition ::checkClampLocation( this->giveDomain(), this->elementNumber, this->localPosition.at( 1 ), this->BcKnotSpan.at( 1 ) );
}


void IGAClampBoundaryConditionLagrange ::assemble( SparseMtrx &answer, TimeStep *tStep,
    CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s,
    double scale,
    void *lock )
{
    // the constraint is linear in the unknowns, so all stiffness types get the same matrix
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


    IGAElement *elem = static_cast<IGAElement *>( el );
    // element-local code numbers of the dofs of the current knot span
    IntArray lc;
    this->giveLocalCodeNumbers( lc, elem, this->BcKnotSpan );

    // the constraint row (and its transpose) for the dofs of the clamped knot span only
    FloatMatrix dNdxi_m( dNdxi.giveSize(), 1 );
    for ( int i = 1; i <= dNdxi.giveSize(); i++ ) {
        dNdxi_m.at( i, 1 ) = scale * dNdxi.at( i );
    }

    // global equation numbers of the dofs of the clamped knot span
    IntArray lcGlob_r, lcGlob_c, loc_r( lc.giveSize() ), loc_c( lc.giveSize() );
    elem->giveLocationArray( lcGlob_r, r_s );
    elem->giveLocationArray( lcGlob_c, c_s );
    for ( int i = 1; i <= lc.giveSize(); i++ ) {
        loc_r.at( i ) = lc.at( i ) ? lcGlob_r.at( lc.at( i ) ) : 0;
        loc_c.at( i ) = lc.at( i ) ? lcGlob_c.at( lc.at( i ) ) : 0;
    }

    // equation number of the multiplier
    IntArray lcGlob_lam_r(1), lcGlob_lam_c(1);
    lcGlob_lam_r.at( 1 ) = r_s.giveDofEquationNumber( *lmdm->begin() );
    lcGlob_lam_c.at( 1 ) = c_s.giveDofEquationNumber( *lmdm->begin() );

    FloatMatrix dNdxi_mT;
    dNdxi_mT.beTranspositionOf( dNdxi_m );

#ifdef _OPENMP
    if ( lock ) {
        omp_set_lock( static_cast<omp_lock_t *>( lock ) );
    }
#endif
    answer.assemble( loc_r, lcGlob_lam_c, dNdxi_m );
    answer.assemble( lcGlob_lam_r, loc_c, dNdxi_mT );
#ifdef _OPENMP
    if ( lock ) {
        omp_unset_lock( static_cast<omp_lock_t *>( lock ) );
    }
#endif


}

void IGAClampBoundaryConditionLagrange ::assembleVector( FloatArray &answer, TimeStep *tStep,
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
        IntArray lc, lcGlob;
        FloatArray ur, xr;
        this->giveLocalCodeNumbers( lc, elem, this->BcKnotSpan );

        // displacements of the current knot span
        ur.resize( dNdxi.giveSize() );
        for ( int i = 1; i <= lc.giveSize(); i++ ) {
            if ( lc.at( i ) ) {
                ur.at( i ) = u.at( lc.at( i ) );
            }
        }

        // current value of the multiplier
        double lambda;
        Dof *mdof = *( lmdm->begin() );
        lambda =  mdof->giveUnknown( mode, tStep ) ;


        // global equation numbers of the dofs of the clamped knot span only
        IntArray loc( lc.giveSize() );
        elem->giveLocationArray( lcGlob, s );
        for ( int i = 1; i <= lc.giveSize(); i++ ) {
            loc.at( i ) = lc.at( i ) ? lcGlob.at( lc.at( i ) ) : 0;
        }
        // forces of the multiplier on the span dofs
        FloatArray fe = lambda * dNdxi;

        // residual of the constraint, n . du/dxi = c . u_span
        FloatArray der(1);
        der.at(1) = dNdxi.dotProduct( ur );
        IntArray lcGlob_lam(1);
        lcGlob_lam.at(1) = s.giveDofEquationNumber( *lmdm->begin() );

#ifdef _OPENMP
        if ( lock ) {
            omp_set_lock( static_cast<omp_lock_t *>( lock ) );
        }
#endif
        answer.assemble( fe, loc );
        answer.assemble( der, lcGlob_lam );
#ifdef _OPENMP
        if ( lock ) {
            omp_unset_lock( static_cast<omp_lock_t *>( lock ) );
        }
#endif
    }
}

int IGAClampBoundaryConditionLagrange ::giveLocalCodeNumbers( IntArray &answer, Element *elem, IntArray &span )
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

void IGAClampBoundaryConditionLagrange ::giveLagrangianMultiplierLocationArray( const UnknownNumberingScheme &r_s, IntArray &answer )
{
    answer.resize( 1 );
    int l  = r_s.giveDofEquationNumber( *this->lmdm->begin() );
    answer.at( 1 ) = l;
}

void IGAClampBoundaryConditionLagrange ::giveLocationArrays( std ::vector<IntArray> &rows, std ::vector<IntArray> &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s )
{
    IntArray r_loc, c_loc;
    rows.resize( 3 );
    cols.resize( 3  );
    IntArray dofIdArray = {
        D_u, D_v
    };
    IntArray lambdaeq( 1 );
    lambdaeq.at( 1 ) = r_s.giveDofEquationNumber( *lmdm->begin() );

    // all dofs of the element (the clamped span is a subset)
    auto elem = this->giveDomain()->giveElement( this->elementNumber );
    elem->giveLocationArray( r_loc, r_s );
    elem->giveLocationArray( c_loc, c_s );


    // column block
    rows[0 ] = r_loc;
    cols[0 ] = lambdaeq;
    // row block
    cols[1 ] = c_loc;
    rows[1 ] = lambdaeq;
    // diagonal entry (some sparse matrix implementations require it)
    rows[2 ] = lambdaeq;
    cols[2 ] = lambdaeq;

}
} // namespace oofem
