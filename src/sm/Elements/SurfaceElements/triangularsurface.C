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

#include "sm/Elements/SurfaceElements/triangularsurface.h"
#include "fei3dquadlin.h"
#include "fei2dtrlin.h"
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

namespace oofem {
REGISTER_Element( TriangularSurface );

FEI2dTrLin TriangularSurface ::interpolation( 1, 2 );

TriangularSurface ::TriangularSurface( int n, Domain *aDomain ) :
    Structural3DSurfaceElement( n, aDomain )
{
    this->GtoLRotationMatrix = NULL;
    numberOfDofMans          = 3;
    numberOfGaussPoints      = 1;
}


FEInterpolation *TriangularSurface ::giveInterpolation() const { return &interpolation; }


void TriangularSurface ::initializeFrom( InputRecord &ir )
{
    numberOfGaussPoints = 1;
    Structural3DSurfaceElement ::initializeFrom( ir );

    IR_GIVE_OPTIONAL_FIELD( ir, this->prestrain, _IFT_TriangularSurface_prestrain );
    if ( this->prestrain ) {
        FloatArray F0_temp;
        IR_GIVE_OPTIONAL_FIELD( ir, F0_temp, _IFT_TriangularSurface_f0 );
        this->F0 = FloatArrayF<9>( F0_temp );
    }  

    FloatArray normal;
    IR_GIVE_OPTIONAL_FIELD( ir, normal, _IFT_TriangularSurface_normal );
    if ( normal.giveSize()>0) {
        this->referenceNormal = FloatArrayF<3>( normal );
        this->isReferenceNormalGiven = 1;
    }

    FloatArray G2L;
    IR_GIVE_OPTIONAL_FIELD( ir, G2L, _IFT_TriangularSurface_g2l );
    if ( G2L.giveSize() > 0 ) {
        this->G2L_prescribed        = FloatArrayF<9>( G2L );
        this->prescribedOrientation = 1;
    }

    FloatArray adir;
    IR_GIVE_OPTIONAL_FIELD( ir, adir, _IFT_TriangularSurface_anisodir );
    if ( adir.giveSize() > 0 ) {
        this->anisoDir = FloatArrayF<3>( adir );
        this->anisoDirGiven  = 1;
    }
    
}


FloatArray TriangularSurface::giveNormal( GaussPoint *gp ) const
{
    FloatArray e1, e2, e3, help;
    e1.beDifferenceOf( this->giveNode( 2 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );
    help.beDifferenceOf( this->giveNode( 3 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );

    // let us normalize e1'
    e1.normalize();

    // compute e3' : vector product of e1' x help
    e3.beVectorProductOf( e1, help );
    e3.normalize();

    return e3;
}


FEICellGeometry *
TriangularSurface::giveCellGeometryWrapper()
{
    if ( cellGeometryWrapper ) {
        return cellGeometryWrapper;
    } else {
        this->computeLocalNodalCoordinates( lc );
        return ( cellGeometryWrapper = new FEIVertexListGeometryWrapper( lc ) );
    }
}


void TriangularSurface ::computeLocalNodalCoordinates( std::vector<FloatArray> &lxy )
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }


    lxy.resize( 3 );
    for ( int i = 0; i < 3; i++ ) {
        const auto &nc = this->giveNode( i + 1 )->giveCoordinates();
        lxy[i].beProductOf( *GtoLRotationMatrix, nc );
    }
}


const FloatMatrix *
TriangularSurface ::computeGtoLRotationMatrix()
// Returns the rotation matrix of the receiver of the size [3,3]
// coords(local) = T * coords(global)
//
// local coordinate (described by vector triplet e1',e2',e3') is defined as follows:
//
// e1'    : [N2-N1]    Ni - means i - th node
// help   : [N3-N1]
// e3'    : e1' x help
// e2'    : e3' x e1'
{
    if ( GtoLRotationMatrix == NULL ) {
        //if ( this->prescribedOrientation ) {
        //    ////GtoLRotationMatrix = new FloatMatrix( 3, 3 );
        //    //FloatArray g2l( this->G2L_prescribed );
        //    //FloatMatrix g2l_test;
        //    //GtoLRotationMatrix = new FloatMatrix( 3, 3 );
        //    //GtoLRotationMatrix->beMatrixForm( g2l );

        //    ////////////////////////////////////
        //    //g2l_test = *GtoLRotationMatrix;
        //    FloatArray e1, e2, e3, help;

        //    // compute e1' = [N2-N1]  and  help = [N3-N1]
        //    e1.beDifferenceOf( this->giveNode( 2 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );
        //    help.beDifferenceOf( this->giveNode( 3 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );

        //    // let us normalize e1'
        //    e1.normalize();

        //    // compute e3' : vector product of e1' x help
        //    e3.beVectorProductOf( e1, help );
        //    // let us normalize
        //    e3.normalize();

        //    // now from e3' x e1' compute e2'
        //    e2.beVectorProductOf( e3, e1 );

        //    //
        //    GtoLRotationMatrix = new FloatMatrix( 3, 3 );

        //    for ( int i = 1; i <= 3; i++ ) {
        //        GtoLRotationMatrix->at( 1, i ) = e1.at( i );
        //        GtoLRotationMatrix->at( 2, i ) = e2.at( i );
        //        GtoLRotationMatrix->at( 3, i ) = e3.at( i );
        //    }

        //    /////////////
        //    //FloatMatrix err = g2l_test - *GtoLRotationMatrix;

        //} else {
            FloatArray e1, e2, e3, help;

            // compute e1' = [N2-N1]  and  help = [N3-N1]
            e1.beDifferenceOf( this->giveNode( 2 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );
            help.beDifferenceOf( this->giveNode( 3 )->giveCoordinates(), this->giveNode( 1 )->giveCoordinates() );

            // let us normalize e1'
            e1.normalize();

            // compute e3' : vector product of e1' x help
            e3.beVectorProductOf( e1, help );
            // let us normalize
            e3.normalize();

            // now from e3' x e1' compute e2'
            e2.beVectorProductOf( e3, e1 );

            //
            GtoLRotationMatrix = new FloatMatrix( 3, 3 );

            for ( int i = 1; i <= 3; i++ ) {
                GtoLRotationMatrix->at( 1, i ) = e1.at( i );
                GtoLRotationMatrix->at( 2, i ) = e2.at( i );
                GtoLRotationMatrix->at( 3, i ) = e3.at( i );
            }
        //}
    }
    //GtoLRotationMatrix->printYourself();
    return GtoLRotationMatrix;
}


bool TriangularSurface ::computeGtoLRotationMatrix( FloatMatrix &answer )
// Returns the rotation matrix of the receiver of the size [8,12]
// r(local) = T * r(global)
// for one node (r written transposed): {u,v} = T * {u,v,w} NOT NOW
// for one node (r written transposed): {u,v,w} = T * {u,v,w}
{
    // test if pereviously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize( 9, 9 );
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at( 1, i ) = answer.at( 4, i + 3 ) = answer.at( 7, i + 6 )  = GtoLRotationMatrix->at( 1, i );
        answer.at( 2, i ) = answer.at( 5, i + 3 ) = answer.at( 8, i + 6 )  = GtoLRotationMatrix->at( 2, i );
        answer.at( 3, i ) = answer.at( 6, i + 3 ) = answer.at( 9, i + 6 )  = GtoLRotationMatrix->at( 3, i );
    }


    return 1;
}
int TriangularSurface ::computeLoadGToLRotationMtrx( FloatMatrix &answer )
// Returns the rotation matrix of the receiver of the size [6,6]
// f(local) = T * f(global)
{
    // test if previously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }

    answer.resize( 6, 6 );
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at( 1, i ) = answer.at( 4, i + 3 ) = GtoLRotationMatrix->at( 1, i );
        answer.at( 2, i ) = answer.at( 5, i + 3 ) = GtoLRotationMatrix->at( 2, i );
        answer.at( 3, i ) = answer.at( 6, i + 3 ) = GtoLRotationMatrix->at( 3, i );
    }

    return 1;
}

int TriangularSurface ::computeGlobalCoordinates( FloatArray &answer, const FloatArray &lcoords )
{
//    FEInterpolation *fei = this->giveInterpolation();
//#ifdef DEBUG
//    if ( !fei ) {
//        answer.clear();
//        return false;
//    }
//#endif
    //fei->local2global3D( answer, lcoords, FEIElementGeometryWrapper( this ) );
   
    // Vertex coordinates
    const auto &p1 = this->giveNode( 1 )->giveCoordinates();
    const auto &p2 = this->giveNode( 2 )->giveCoordinates();
    const auto &p3 = this->giveNode( 3 )->giveCoordinates();

    // Local coordinates
    double l1 = lcoords.at( 1 );
    double l2 = lcoords.at( 2 );
    double l3 = 1.0 - l1 - l2;

    answer.resize( 3 );
    answer.zero();
    answer.at( 1 ) = l1 * p1.at( 1 ) + l2 * p2.at( 1 ) + l3 * p3.at( 1 );
    answer.at( 2 ) = l1 * p1.at( 2 ) + l2 * p2.at( 2 ) + l3 * p3.at( 2 );
    answer.at( 3 ) = l1 * p1.at( 3 ) + l2 * p2.at( 3 ) + l3 * p3.at( 3 );

    return true;
}

int TriangularSurface::giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep )
{
    if ( type == IST_PrestrainDeformationGradient ) {
        if ( this->prestrain ) {
            answer = FloatArray( this->F0 );
            FloatMatrix R;
            if ( this->giveRotationMatrix( R ) ) {
                IntArray pos( { 1, 2, 3 } );
                FloatMatrix R3x3, Fmat, FRt, RFRt;
                R3x3.beSubMatrixOf( R, pos, pos );

                Fmat.beMatrixForm( answer );
                FRt.beProductTOf( Fmat, R3x3 );
                RFRt.beProductOf( R3x3, FRt );
                ///////////
                //FloatMatrix FRt2, RFRt2;
                //FRt2.beProductOf( Fmat, R3x3 );
                //RFRt2.beTProductOf( R3x3, FRt2 );

                //FloatMatrix err = RFRt2 - RFRt;
                ///////////

                answer.beVectorForm( RFRt );
                //answer.beVectorForm( RFRt2 );
            }

            return 1;
        } else {
            return 0;
        }  
    } else if ( type == IST_ElementReferenceNormal ) {
        if ( this->isReferenceNormalGiven ) {
            answer = FloatArray( this->referenceNormal );
            
            FloatMatrix R;
            if ( this->giveRotationMatrix( R ) ) {
                IntArray pos( { 1, 2, 3 } );
                FloatMatrix R3x3;
                R3x3.beSubMatrixOf( R, pos, pos );

                FloatArray RN;
                RN.beProductOf( R3x3, answer );
                answer = RN;
            }
            return 1;
        } else {
            return 0;
        }

    }else if ( type == IST_DisplacementVector ) {
        FloatArray u;
        FloatMatrix N, G2L;
        this->computeVectorOf( VM_Total, tStep, u );

        if ( this->computeGtoLRotationMatrix( G2L ) ) {
            u.rotatedWith( G2L, 't' );
        }
        this->computeNmatrixAt( gp->giveSubPatchCoordinates(), N );
        answer.beProductOf( N, u );
        return 1;


        
    } else if ( type == IST_AnisotropyDirection ){
        if ( this->anisoDirGiven ) {
            answer = FloatArray( this->anisoDir );
            return 1;
        } else {
            return 0;
        }
    }
    else
    {
        return Structural3DSurfaceElement::giveIPValue( answer, gp, type, tStep );
    }
}

//void TriangularSurface::computeBHmatrixAt(GaussPoint* gp, FloatMatrix& answer) {
//    
//    //Structural3DSurfaceElement::computeBHmatrixAt( gp, answer );
//
//    ///////////////////////////////////////// Just some checks
//    // Check transformations
//    FloatMatrix dNdx;
//    FEInterpolation *interp = this->giveInterpolation();
//    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() ); // nno x 2
//    // enlarge it
//    FloatMatrix dNdx_L( 3, 3 );
//    dNdx_L.setSubMatrix( dNdx, 1, 1 );
//
//
//    // Check how it transforms
//    FloatMatrix R, R3x3, Rotation_Ref_2_Ref_Old; 
//
//    this->giveRotationMatrix( R );
//    IntArray pos( { 1, 2, 3 } );
//    R3x3.beSubMatrixOf( R, pos, pos );
//    FloatMatrix RefRotation( R3x3 );
//
//    //this->givePrescribedRotation( RefRotation ); 
//
//    //if ( this->givePrescribedRotation( RefRotation ) ) {
//    //    Rotation_Ref_2_Ref_Old.beProductTOf( RefRotation, R3x3 );
//    //} else {
//    //    Rotation_Ref_2_Ref_Old.beTranspositionOf( R3x3 );
//    //}
//
//    //Rotation_Ref_2_Ref_Old.beProductTOf( RefRotation, R3x3 );
//
//    // Get F0
//    FloatMatrix F0m;
//    FloatArray F0a( { 1., 1., 0., 0., 0., 0., 0., 0., 0. } );
//    if ( this->prestrain ) {
//        F0a = FloatArray( this->F0 );
//        FloatMatrix Fmat, FRt;
//        Fmat.beMatrixForm( F0a );
//        FRt.beProductTOf( Fmat, R3x3 );
//        F0m.beProductOf( R3x3, FRt );
//    } else {
//        F0m.beMatrixForm( F0a );
//    }
//    
//
//    // Try for 1 SF
//    FloatArray dN1dX;
//    FloatMatrix dNdx_Lt;
//    dNdx_Lt.beTranspositionOf( dNdx_L );
//
//
//    FloatMatrix dNdX_F0;
//    dNdX_F0.beTProductOf( F0m, dNdx_Lt );
//    FloatMatrix dNdX_F0_t;
//    dNdX_F0_t.beTranspositionOf( dNdX_F0 );
//
//    answer.resize( 9, dNdX_F0_t.giveNumberOfRows() * 3 );
//    answer.zero();
//
//    for ( int i = 1; i <= dNdX_F0_t.giveNumberOfRows(); i++ ) {
//        answer.at( 1, 3 * i - 2 ) = dNdX_F0_t.at( i, 1 ); // du/dx
//        answer.at( 2, 3 * i - 1 ) = dNdX_F0_t.at( i, 2 ); // dv/dy
//        answer.at( 3, 3 * i - 0 ) = dNdX_F0_t.at( i, 3 ); // dw/dz
//        answer.at( 4, 3 * i - 1 ) = dNdX_F0_t.at( i, 3 ); // dv/dz
//        answer.at( 7, 3 * i - 0 ) = dNdX_F0_t.at( i, 2 ); // dw/dy
//        answer.at( 5, 3 * i - 2 ) = dNdX_F0_t.at( i, 3 ); // du/dz
//        answer.at( 8, 3 * i - 0 ) = dNdX_F0_t.at( i, 1 ); // dw/dx
//        answer.at( 6, 3 * i - 2 ) = dNdX_F0_t.at( i, 2 ); // du/dy
//        answer.at( 9, 3 * i - 1 ) = dNdX_F0_t.at( i, 1 ); // dv/dx
//    }
//
//    //FloatMatrix answerTest;
//    //Structural3DSurfaceElement::computeBHmatrixAt( gp, answerTest );
//    //FloatMatrix Err = answerTest - answer;
//
//
//    //// transformation to the original reference frame
//    //FloatMatrix dNdX_F0_ref;
//    //dNdX_F0_ref.beProductTOf( dNdX_F0_t, Rotation_Ref_2_Ref_Old);
//
//
//}
//
//void TriangularSurface::computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) 
//{
//    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
//    // the receiver at time step tStep.
//    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.
//
//    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
//    FloatArray u;
//    this->computeVectorOf( { D_u, D_v, D_w }, VM_Total, tStep, u ); // solution vector
//    if ( initialDisplacements ) {
//        u.subtract( *initialDisplacements );
//    }
//
//    // Displacement gradient H = du/dX
//    FloatMatrix B;
//    Structural3DSurfaceElement::computeBHmatrixAt( gp, B );
//
//    answer.beProductOf( B, u );
//
//    // Deformation gradient F = H + I
//    answer.at( 1 ) += 1.0;
//    answer.at( 2 ) += 1.0;
//}


} // end namespace oofem
