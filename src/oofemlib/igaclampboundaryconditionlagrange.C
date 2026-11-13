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
    this->lmdm = std ::make_unique<Node>( 0, this->giveDomain() );
    this->lmdm->appendDof( new MasterDof( this->lmdm.get(), (DofIDItem)( this->giveDomain()->giveNextFreeDofID() ) ) );
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

    auto el = this->giveDomain()->giveElement( this->elementNumber );
    FloatArray c;
    IntArray loc_r, loc_c;
    IGAClampBoundaryCondition ::computeClampVector( c, el, this->localPosition, this->BcKnotSpan );
    IGAClampBoundaryCondition ::giveSpanLocationArray( loc_r, el, this->BcKnotSpan, r_s );
    IGAClampBoundaryCondition ::giveSpanLocationArray( loc_c, el, this->BcKnotSpan, c_s );

    // the constraint column c and row c^T coupling the dofs of the clamped knot span and the multiplier
    FloatMatrix cm( c.giveSize(), 1 ), cmT;
    for ( int i = 1; i <= c.giveSize(); i++ ) {
        cm.at( i, 1 ) = scale * c.at( i );
    }
    cmT.beTranspositionOf( cm );

    // equation number of the multiplier
    IntArray lambda_r = { r_s.giveDofEquationNumber( *this->lmdm->begin() ) };
    IntArray lambda_c = { c_s.giveDofEquationNumber( *this->lmdm->begin() ) };

#ifdef _OPENMP
    if ( lock ) {
        omp_set_lock( static_cast<omp_lock_t *>( lock ) );
    }
#endif
    answer.assemble( loc_r, lambda_c, cm );
    answer.assemble( lambda_r, loc_c, cmT );
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
    if ( type != InternalForcesVector ) {
        return;
    }

    auto el = this->giveDomain()->giveElement( this->elementNumber );
    FloatArray c, u;
    IntArray loc;
    IGAClampBoundaryCondition ::computeClampVector( c, el, this->localPosition, this->BcKnotSpan );
    IGAClampBoundaryCondition ::giveSpanDisplacements( u, el, this->BcKnotSpan, tStep );
    IGAClampBoundaryCondition ::giveSpanLocationArray( loc, el, this->BcKnotSpan, s );

    // forces lambda c of the multiplier on the span dofs and the residual c . u_span of the constraint
    double lambda = ( *this->lmdm->begin() )->giveUnknown( mode, tStep );
    FloatArray fe = lambda * c;
    FloatArray der = { c.dotProduct( u ) };
    IntArray lambdaLoc = { s.giveDofEquationNumber( *this->lmdm->begin() ) };

#ifdef _OPENMP
    if ( lock ) {
        omp_set_lock( static_cast<omp_lock_t *>( lock ) );
    }
#endif
    answer.assemble( fe, loc );
    answer.assemble( der, lambdaLoc );
#ifdef _OPENMP
    if ( lock ) {
        omp_unset_lock( static_cast<omp_lock_t *>( lock ) );
    }
#endif
}

void IGAClampBoundaryConditionLagrange ::giveLagrangianMultiplierLocationArray( const UnknownNumberingScheme &r_s, IntArray &answer )
{
    answer.resize( 1 );
    int l  = r_s.giveDofEquationNumber( *this->lmdm->begin() );
    answer.at( 1 ) = l;
}

void IGAClampBoundaryConditionLagrange ::giveLocationArrays( std ::vector<IntArray> &rows, std ::vector<IntArray> &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s )
{
    // all dofs of the element (the clamped span is a subset) and the multiplier
    IntArray r_loc, c_loc;
    auto elem = this->giveDomain()->giveElement( this->elementNumber );
    elem->giveLocationArray( r_loc, r_s );
    elem->giveLocationArray( c_loc, c_s );
    IntArray lambda_r = { r_s.giveDofEquationNumber( *this->lmdm->begin() ) };
    IntArray lambda_c = { c_s.giveDofEquationNumber( *this->lmdm->begin() ) };

    rows.resize( 3 );
    cols.resize( 3 );
    // column block
    rows[0] = r_loc;
    cols[0] = lambda_c;
    // row block
    rows[1] = lambda_r;
    cols[1] = c_loc;
    // diagonal entry (some sparse matrix implementations require it)
    rows[2] = lambda_r;
    cols[2] = lambda_c;
}
} // namespace oofem
