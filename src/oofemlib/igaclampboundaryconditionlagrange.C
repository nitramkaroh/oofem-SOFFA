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
// #include "feitspline.h"

// #ifdef _OPENMP
// #include <omp.h>
// #endif

namespace oofem {
REGISTER_BoundaryCondition( IGAClampBoundaryConditionLagrange );


IGAClampBoundaryConditionLagrange ::IGAClampBoundaryConditionLagrange( int n, Domain *d ) :
    ActiveBoundaryCondition( n, d )
{
}


void IGAClampBoundaryConditionLagrange ::initializeFrom( InputRecord &ir )
{
    ActiveBoundaryCondition ::initializeFrom( ir );
    // BoundaryCondition ::initializeFrom( ir );
    IR_GIVE_FIELD( ir, this->elementNumber, _IFT_IGAClampBoundaryConditionLagrange_ElementNumber );
    double lp;
    IR_GIVE_FIELD( ir, lp, _IFT_IGAClampBoundaryConditionLagrange_LocalPosition );
    this->localPosition.resize( 1 );
    this->localPosition.at( 1 ) = lp;

    int ks;
    IR_GIVE_FIELD( ir, ks, _IFT_IGAClampBoundaryConditionLagrange_BcKnotSpan );
    this->BcKnotSpan.resize( 1 );
    this->BcKnotSpan.at( 1 ) = ks;


     //these are internal lagrange multipliers
     for ( int pos = 0; pos < 1; ++pos ) {
        //lmdm.push_back( std ::unique_ptr<DofManager>( nullptr ) );
        lmdm =  std ::unique_ptr<Node>( new Node( 0, domain ) ) ;
        lmdm->appendDof( new MasterDof( this->lmdm.get(), (DofIDItem)( this->giveDomain()->giveNextFreeDofID() ) ) );
         //lmdm = std ::make_unique<DofManager>( 1, domain );
         //lmdm->appendDof( new MasterDof( this->lmdm.get(), (DofIDItem)( this->giveDomain()->giveNextFreeDofID() ) ) );

        //this->pressureDof = std::make_unique<ElementDofManager>( 1, aDomain, this );
        //this->pressureDof->appendDof( new MasterDof( this->pressureDof.get(), P_f ) );
    }
}


void IGAClampBoundaryConditionLagrange ::assemble( SparseMtrx &answer, TimeStep *tStep,
    CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s,
    double scale,
    void *lock )
{
    auto el                 = this->giveDomain()->giveElement( this->elementNumber );
    FEInterpolation *interp = el->giveInterpolation();
    // compute normal
    FloatArray normal;
    auto gw = FEIIGAElementGeometryWrapper( el, &this->BcKnotSpan );
    interp->boundaryEvalNormal( normal, 1, this->localPosition, gw );
    normal.resize( 2 ); // delete the third outofplane direction

    FloatMatrix ders, dersT;
    //interp->evalDerivatives( 1, ders, this->localPosition );
    //interp->evalDerivatives( 1, ders, this->localPosition, this->BcKnotSpan );
    interp->evalDerivatives( 1, ders, this->localPosition, gw );

    FloatMatrix dNdxi_mat;
    dNdxi_mat.resize( 2, 2 * ders.giveNumberOfColumns() );
    for ( int i = 1; i <= ders.giveNumberOfColumns(); i++ ) {
        dNdxi_mat.at( 1, 2 * i - 1 ) = ders.at( 2, i );
        dNdxi_mat.at( 2, 2 * i )     = ders.at( 2, i );
    }

    // dNdXi projected to normal
    FloatArray dNdxi;
    dNdxi.beTProductOf( dNdxi_mat, normal );


    //// Create the matrix
    //FloatMatrix dNdNt;
    //dNdNt.beDyadicProductOf( dNdxi, dNdxi );


    IGAElement *elem = static_cast<IGAElement *>( el );
    // get local code numbers corresponding to currect knotspan
    IntArray lc;
    this->giveLocalCodeNumbers( lc, elem, this->BcKnotSpan );


    //// assemble to the element vector
    //FloatMatrix answerTemp;
    //int ndofs = elem->computeNumberOfDofs();
    //answerTemp.resize( ndofs, ndofs );
    //answerTemp.assemble( dNdNt, lc );

    // assemble to the element vector
    FloatMatrix dNdxi_m;
    int ndofs = elem->computeNumberOfDofs();
    dNdxi_m.resize( ndofs, 1 );
    IntArray lc1(1);
    lc1.at( 1 ) = 1;
    dNdxi_m.assemble( dNdxi, lc, lc1 );

    // assemble to the global vector
    IntArray lcGlob_r, lcGlob_c;
    elem->giveLocationArray( lcGlob_r, r_s );
    elem->giveLocationArray( lcGlob_c, c_s );
    //answer.assemble( lcGlob_r, lcGlob_c, this->penalty * answerTemp );


    ////
    //IntArray dofIdArray  = domain->giveDefaultNodeDofIDArry();
    //IntArray dofIdArray = {
    //    D_u, D_v
    //};

    IntArray lcGlob_lam_r(1), lcGlob_lam_c(1);
    lcGlob_lam_r.at( 1 ) = r_s.giveDofEquationNumber( *lmdm->begin() );
    lcGlob_lam_c.at( 1 ) = c_s.giveDofEquationNumber( *lmdm->begin() );

    //lmdm->giveLocationArray( dofIdArray, lcGlob_lam_r, r_s );
    //lmdm->giveLocationArray( dofIdArray, lcGlob_lam_c, c_s );

    //FloatMatrix dNdxi_mat( dNdxi ), dNdxi_matT;
    FloatMatrix dNdxi_mT;
    dNdxi_mT.beTranspositionOf( dNdxi_m );

    answer.assemble( lcGlob_r, lcGlob_lam_c, dNdxi_m );
    answer.assemble( lcGlob_lam_r, lcGlob_c, dNdxi_mT );


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

        // compute normal
        FloatArray normal;
        auto gw = FEIIGAElementGeometryWrapper( el, &this->BcKnotSpan );
        interp->boundaryEvalNormal( normal, 1, this->localPosition, gw );
        normal.resize( 2 ); // delete the third outofplane direction

        FloatMatrix ders, dersT;
        //FloatMatrix ders2;
        // obtain 0 and 1st derivative of nonzero shapefunctions
        //interp->evalDerivatives( 1, ders, this->localPosition );
        //interp->evalDerivatives( 1, ders, this->localPosition, this->BcKnotSpan );
        interp->evalDerivatives( 1, ders, this->localPosition, gw );

        //( ders2 - ders ).printYourself();
        // dersT.beTranspositionOf( ders );

        FloatMatrix dNdxi_mat;
        dNdxi_mat.resize( 2, 2 * ders.giveNumberOfColumns() );
        for ( int i = 1; i <= ders.giveNumberOfColumns(); i++ ) {
            dNdxi_mat.at( 1, 2 * i - 1 ) = ders.at( 2, i );
            dNdxi_mat.at( 2, 2 * i )     = ders.at( 2, i );
        }


        // dNdXi projected to normal
        FloatArray dNdxi;
        dNdxi.beTProductOf( dNdxi_mat, normal );


        // get displacements of all element nodes
        FloatArray u;
        IGAElement *elem = static_cast<IGAElement *>( el );
        elem->computeVectorOf( VM_Total, tStep, u );

        // get local code numbers corresponding to currect knotspan
        IntArray lc, lcGlob;
        FloatArray ur, xr;
        this->giveLocalCodeNumbers( lc, elem, this->BcKnotSpan );

        // extract displacements related to the current knotspan
        ur.resize( dNdxi.giveSize() );
        for ( int i = 1; i <= lc.giveSize(); i++ ) {
            if ( lc.at( i ) ) {
                ur.at( i ) = u.at( lc.at( i ) );
            }
        }


        // get unknown lambda ... 
        double lambda;
        Dof *mdof = *( lmdm->begin() );
        lambda =  mdof->giveUnknown( mode, tStep ) ;


        // assemble to the element vector, 
        FloatArray answerTemp;
        int ndofs = elem->computeNumberOfDofs();
        answerTemp.resize( ndofs );
        answerTemp.assemble( lambda * dNdxi, lc );


        // assemble to the global vector at positions of classical displacement
        elem->giveLocationArray( lcGlob, s );
        answer.assemble(  answerTemp, lcGlob );

        ////////////////////////////////////
        // Assemble the lambda part
        //IntArray dofIdArray = {
        //    D_u, D_v
        //};
        //IntArray dofIdArray = domain->giveDefaultNodeDofIDArry();

        FloatArray der(1);
        der.at(1) = dNdxi.dotProduct( ur );
        // 
        IntArray lcGlob_lam(1);
        //lmdm->giveLocationArray( dofIdArray, lcGlob_lam, s );
        lcGlob_lam.at(1) = s.giveDofEquationNumber( *lmdm->begin() );

        answer.assemble( der, lcGlob_lam );
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

    if ( nsd == 1 ) { // Change compared to the general funciton!!
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


// double
// IGAClampBoundaryCondition ::computeTangentFromContact( FloatMatrix &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep )
//{
//     double gap;
//     FloatArray Nv;
//     this->computeGap( gap, masterNode, slaveNode, tStep );
//     this->computeNormalMatrixAt( Nv, masterNode, slaveNode, tStep );
//     answer.initFromVector( Nv, false );
//
//     return gap;
//     // answer.times(this->penalty);
//     /*    if ( gap > 0.0 ) {
//      *  int size = this->masterSet.giveSize();
//      *  FloatMatrix help(2 * size, 2 * size);
//      *  answer.beUnitMatrix();
//      *  }*/
// }
//
// void IGAClampBoundaryCondition ::computeGap( double &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep )
//{
//     FloatArray uS, uM;
//     auto xs           = slaveNode->giveCoordinates();
//     auto xm           = masterNode->giveCoordinates();
//     FloatArray normal = xs - xm;
//     double norm       = normal.computeNorm();
//     if ( norm < 1.0e-8 ) {
//         OOFEM_ERROR( "Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.", masterNode->giveGlobalNumber(), slaveNode->giveGlobalNumber() );
//     } else {
//         normal.times( 1.0 / norm );
//     }
//
//     slaveNode->giveUnknownVector( uS, { D_u, D_v }, VM_Total, tStep, true );
//     masterNode->giveUnknownVector( uM, { D_u, D_v }, VM_Total, tStep, true );
//     xs.add( uS );
//     xm.add( uM );
//     FloatArray dx = xs - xm;
//     answer        = dx.dotProduct( normal );
// }
//
//
// void IGAClampBoundaryCondition ::computeNormalMatrixAt( FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *TimeStep )
//{
//     const auto &xs = slaveNode->giveCoordinates();
//     const auto &xm = masterNode->giveCoordinates();
//     auto normal    = xs - xm;
//     double norm    = normal.computeNorm();
//     if ( norm < 1.0e-8 ) {
//         OOFEM_ERROR( "Couldn't compute normal between master node (num %d) and slave node (num %d), nodes are too close to each other.", masterNode->giveGlobalNumber(), slaveNode->giveGlobalNumber() );
//     } else {
//         normal.times( 1.0 / norm );
//     }
//     // The normal is not updated for node2node which is for small deformations only
//     // C = {n -n}
//     answer = {
//         normal.at( 1 ), normal.at( 2 ),
//         -normal.at( 1 ), -normal.at( 2 )
//     };
// }
//
//
// void IGAClampBoundaryCondition ::computeExternalForcesFromContact( FloatArray &answer, Node *masterNode, Node *slaveNode, TimeStep *tStep )
//{
//     answer.resize( 1 );
//     this->computeGap( answer.at( 1 ), masterNode, slaveNode, tStep );
//     if ( answer.at( 1 ) > 0.0 ) {
//         answer.at( 1 ) = 0.0;
//     }
// }
//
//
void IGAClampBoundaryConditionLagrange ::giveLagrangianMultiplierLocationArray( const UnknownNumberingScheme &r_s, IntArray &answer )
{
    // assemble location array
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
     //std ::vector<IntArray> lambdaeq;
     //this->giveLagrangianMultiplierLocationArray( r_s, lambdaeq );
     IntArray lambdaeq( 1 );
     lambdaeq.at( 1 ) = r_s.giveDofEquationNumber( *lmdm->begin() );



    auto elem = this->giveDomain()->giveElement( this->elementNumber );
    elem->giveLocationArray( r_loc, r_s );
    elem->giveLocationArray( c_loc, c_s );


    // column block
    rows[0 ] = r_loc;
    cols[0 ] = lambdaeq;
    // row block
    cols[1 ] = c_loc;
    rows[1 ] = lambdaeq;
    // diagonal enry (some sparse mtrx implementation requaire this)
    rows[2 ] = lambdaeq;
    cols[2 ] = lambdaeq;

 }
} // namespace oofem
