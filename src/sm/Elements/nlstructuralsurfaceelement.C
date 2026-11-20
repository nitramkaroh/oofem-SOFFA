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

#include "sm/Elements/nlstructuralsurfaceelement.h"

#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "feinterpol.h"
#include "domain.h"
#include "material.h"
#include "sm/CrossSections/simplesurfacecrosssection.h"
#include "integrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dynamicinputrecord.h"
#include "gausspoint.h"
#include "engngm.h"
#include "mathfem.h"
#include "pointload.h"
#include "node.h"

#include <cmath>

namespace oofem {
NLStructuralSurfaceElement::NLStructuralSurfaceElement( int n, Domain *aDomain ) :
    NLStructuralElement( n, aDomain )
// Constructor. Creates an element with number n, belonging to aDomain.
{}

void NLStructuralSurfaceElement::computeNmatrixAt( const FloatArray &lcoords, FloatMatrix &answer )
{
    // ndofs x (ndofs nno): each shape function multiplies every translational DOF of its node (2 or 3 per node)
    FloatArray N;
    this->giveInterpolation()->evalN( N, lcoords, FEIElementGeometryWrapper( this ) );
    IntArray dofs;
    this->giveDofManDofIDMask( 1, dofs );
    answer.beNMatrixOf( N, dofs.giveSize() );
}

int NLStructuralSurfaceElement::computeGlobalCoordinates( FloatArray &answer, const FloatArray &lcoords )
{
    FloatArray N;
    this->giveInterpolation()->evalN( N, lcoords, FEIElementGeometryWrapper( this ) );
    answer.resize( 3 );
    answer.zero();
    for ( int i = 0; i < N.giveSize(); ++i ) {
        const auto &coords = this->giveNode( i + 1 )->giveCoordinates();
        for ( int j = 0; j < coords.giveSize() && j < 3; ++j ) {
            answer[j] += N[i] * coords[j];
        }
    }
    return 1;
}

void NLStructuralSurfaceElement::computeLoadVector( FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep )
{
    // Reject before the base wrapper tries its spatial localization, which is not implemented in the surface frame.
    if ( type == ExternalForcesVector && dynamic_cast<PointLoad *>( load ) ) {
        OOFEM_ERROR( "element %d (%s): PointLoad is not supported by surface elements; use nodal or distributed loads", this->giveNumber(), this->giveClassName() );
    }
    StructuralElement::computeLoadVector( answer, load, type, mode, tStep );
}

void NLStructuralSurfaceElement::computeMaterialFrame( FloatMatrix &answer )
{
    FloatMatrix R;
    if ( !this->computeLoadGToLRotationMtrx( R ) ) {
        OOFEM_ERROR( "element %d: the surface material frame is not available", this->giveNumber() );
    }
    if ( R.giveNumberOfRows() == 3 && R.giveNumberOfColumns() == 3 ) {
        answer = R;
    } else if ( R.giveNumberOfRows() == 2 && R.giveNumberOfColumns() == 2 ) {
        // The nodal frame is (tangent t, left normal n). The material axes are (tangent, out of plane, normal):
        // rows e1 = (t_x, t_y, 0), e2 = (0, 0, -1), e3 = (n_x, n_y, 0). Choosing the second axis opposite to
        // global z makes this a right-handed spatial frame (e1 x e2 = e3).
        answer.resize( 3, 3 );
        answer.at( 1, 1 ) = R.at( 1, 1 );
        answer.at( 1, 2 ) = R.at( 1, 2 );
        answer.at( 2, 3 ) = -1.;
        answer.at( 3, 1 ) = R.at( 2, 1 );
        answer.at( 3, 2 ) = R.at( 2, 2 );
    } else {
        OOFEM_ERROR( "element %d: invalid surface material frame dimensions", this->giveNumber() );
    }
}

int NLStructuralSurfaceElement::giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep )
{
    if ( type == IST_DisplacementVector ) {
        // u = N * u_global at gp: the unknowns come in the element frame, R^T rotates them back
        FloatArray u;
        FloatMatrix N, R;
        this->computeVectorOf( VM_Total, tStep, u );
        if ( this->computeGtoLRotationMatrix( R ) ) {
            u.rotatedWith( R, 't' );
        }
        this->computeNmatrixAt( gp->giveSubPatchCoordinates(), N );
        answer.beProductOf( N, u );
        answer.resizeWithValues( 3 ); // planar elements: (u, v, 0)
        return 1;
    } else if ( type == IST_GlobalDeformationGradientTensor ) {
        auto *status = static_cast<StructuralMaterialStatus *>( this->giveCrossSection()->giveMaterial( gp )->giveStatus( gp ) );
        FloatArray vFmat, normal;
        Tensor2_3d F0s;
        double J0;
        // The status holds the kinematic F (the element stores it after the constitutive call), so the
        // prestrain is applied exactly once here: vFmat = F * F0 in the element frame.
        this->giveMaterialKinematics( status->giveFVector(), gp, vFmat, normal, F0s, J0 );
        // to global coordinates: F_global = R^T (F F0) R, R = global-to-element rotation (rows e1, e2, e3)
        FloatMatrix R, F, FR, globalF;
        this->computeMaterialFrame( R );
        F.beMatrixForm( vFmat );
        FR.beProductOf( F, R );
        globalF.beTProductOf( R, FR );
        answer.beVectorForm( globalF );
        return 1;
    }
    return NLStructuralElement::giveIPValue( answer, gp, type, tStep );
}

bool NLStructuralSurfaceElement::giveMaterialKinematics( const FloatArray &vF, GaussPoint *gp, FloatArray &vFmat, FloatArray &normal, Tensor2_3d &F0s, double &J0 )
{
    Tensor2_3d F0;
    if ( !this->givePrestrainDeformationGradient( F0, gp ) ) {
        // No prestrain, F0 = I_s: F * F0 = F (F e3 = 0), the normal is e3 and the area ratio 1.
        vFmat  = vF;
        normal = { 0., 0., 1. };
        F0s    = Tensor2_3d( 1., 0., 0., 0., 1., 0., 0., 0., 0. );
        J0     = 1.;
        return false;
    }
    vFmat = applyPrestrain( vF, F0 );

    // F maps e3 to zero (F = F I_s), so F F0 = F F0s with F0s = I_s F0, the rows 1, 2 of F0: only the part of
    // F0 in the element plane counts. F0s maps the plane of the original configuration onto the element plane
    // and its normal N to zero. With the third row zero, the only nonzero row of cof(F0s) is the third one,
    // row1(F0) x row2(F0) = J0 N, i.e. cof(F0s) = J0 e3 N^T: N = cof(F0s)^T e3 / J0 and J0 its length.
    F0s = Tensor2_3d( F0( 0, 0 ), F0( 0, 1 ), F0( 0, 2 ), F0( 1, 0 ), F0( 1, 1 ), F0( 1, 2 ), 0., 0., 0. );
    Tensor2_3d cofF0 = F0s.compute_cofactor();
    normal           = { cofF0( 2, 0 ), cofF0( 2, 1 ), cofF0( 2, 2 ) };
    J0               = normal.computeNorm();
    if ( !std::isfinite( J0 ) || !( J0 > 0. ) ) {
        OOFEM_ERROR( "element %d: the prestrain deformation gradient is singular in the element plane or has a non-finite area ratio", this->giveNumber() );
    }
    normal.times( 1. / J0 );
    return true;
}


void NLStructuralSurfaceElement::computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray vF, vFmat, normal;
    this->computeDeformationGradientVector( vF, gp, tStep );

    // The material sees vFmat = F * F0 with the normal of its reference configuration; the returned
    // stress is pulled back to the kinematic F.
    Tensor2_3d F0s;
    double J0;
    bool prestrained = this->giveMaterialKinematics( vF, gp, vFmat, normal, F0s, J0 );

    FloatArray vPmat = cs->giveFirstPKSurfaceStresses( vFmat, normal, gp, tStep );

    if ( prestrained ) {
        // The energy per unit area of the mesh is W(F F0) / J0 and F F0 = F F0s, so by the chain rule
        // P = P_material(F F0) F0s^T / J0.
        FloatArrayF<9> vP( vPmat );
        Tensor2_3d P( vP ), Pmod;
        Pmod( i_3, j_3 ) = 1. / J0 * P( i_3, k_3 ) * F0s( j_3, k_3 );
        answer = Pmod.to_voigt_form();
    } else {
        // P I_s: drop the components with the reference index 3 (Voigt 33, 23, 13), which BH does not use.
        answer = vPmat;
        answer.at( 3 ) = 0.;
        answer.at( 4 ) = 0.;
        answer.at( 5 ) = 0.;
    }

    // Store the work-conjugate F and P in the frame of BH: the internal forces and reactions computed
    // from the status (cached stresses) use this P.
    auto mat = cs->giveMaterial( gp );
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( mat->giveStatus( gp ) );
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( answer );
}


void NLStructuralSurfaceElement::computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    auto cs  = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    auto mat = cs->giveMaterial( gp );
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( mat->giveStatus( gp ) );
    // the kinematic F, stored by computeFirstPKStressVector
    FloatArray vF( status->giveTempFVector() ), vFmat, normal;

    Tensor2_3d F0s;
    double J0;
    if ( !this->giveMaterialKinematics( vF, gp, vFmat, normal, F0s, J0 ) ) {
        // No prestrain: the status already holds F; F0 = I_s on both reference indices only drops the rows
        // and columns with the reference index 3 (Voigt 33, 23, 13).
        answer = cs->giveSurfaceStiffnessMatrix_dPdF_3d( rMode, normal, gp, tStep );
        for ( int i = 1; i <= 9; ++i ) {
            for ( int k = 3; k <= 5; ++k ) {
                answer.at( k, i ) = 0.;
                answer.at( i, k ) = 0.;
            }
        }
        return;
    }

    status->letTempFVectorBe( vFmat ); // the material takes F from the status: give it F * F0
    FloatMatrix Amat = cs->giveSurfaceStiffnessMatrix_dPdF_3d( rMode, normal, gp, tStep );
    status->letTempFVectorBe( vF ); // Restore kinematic F after the constitutive evaluation.

    // Chain rule for P = P_material(F F0s) F0s^T / J0: F0s on both reference indices,
    // A_ijkl = F0s_jm A_material_imkn F0s_ln / J0.
    FloatMatrixF<9, 9> mA( Amat );
    Tensor4_3d A( mA ), Amod;
    Amod( i_3, j_3, k_3, l_3 ) = 1. / J0 * F0s( j_3, m_3 ) * A( i_3, m_3, k_3, n_3 ) * F0s( l_3, n_3 );
    answer = Amod.to_voigt_form();
}


void NLStructuralSurfaceElement ::computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep )
{
    // small-strain path (not reached: checkConsistency requires nlgeo 1)
    auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray normal = { 0., 0., 1. }; // the reference normal in the frame of the element
    answer = cs->giveRealSurfaceStresses( strain, normal, gp, tStep );
}


int NLStructuralSurfaceElement::checkConsistency()
{
    int result = NLStructuralElement::checkConsistency();
    auto *cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveCrossSection() );
    if ( !cs ) {
        OOFEM_ERROR( "element %d (%s): a SimpleSurfaceCrossSection (simplesurfacecs) is required", this->giveNumber(), this->giveClassName() );
    }
    // the small strain path is not implemented by any surface material (isolesurf stops at input)
    if ( this->nlGeometry != 1 ) {
        OOFEM_ERROR( "element %d (%s): only nlgeo 1 is implemented for surface elements", this->giveNumber(), this->giveClassName() );
    }
    // C0 interpolation gives no second gradient, which a material with bending stiffness needs
    for ( auto &iRule : this->integrationRulesArray ) {
        for ( GaussPoint *gp : *iRule ) {
            auto *mat = cs->giveSurfaceMaterial( gp );
            if ( mat->supportsSecondGradient() ) {
                OOFEM_ERROR( "element %d (%s): material %d (%s) needs the second gradient (bending stiffness), which only C1 IGA surface elements provide (degree >= 2, interior knot multiplicities <= degree - 1)",
                    this->giveNumber(), this->giveClassName(), mat->giveNumber(), mat->giveClassName() );
            }
        }
    }
    return result;
}

void NLStructuralSurfaceElement::computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present).
    // The element's own dof mask is used, so planar (u, v) and spatial (u, v, w) surfaces share this code.
    FloatArray u;
    this->computeVectorOf( VM_Total, tStep, u ); // Element unknowns, rotated to the local frame when provided.
    if ( initialDisplacements ) {
        u.subtract( *initialDisplacements );
    }

    // Surface displacement gradient in the element's reference surface frame, H = BH u.
    FloatMatrix BH;
    this->computeBHmatrixAt( gp, BH );

    answer.beProductOf( BH, u );

    // Surface deformation gradient F = H + I_s, where I_s = diag(1, 1, 0).
    answer.at( 1 ) += 1.0;
    answer.at( 2 ) += 1.0;

}



} // end namespace oofem
