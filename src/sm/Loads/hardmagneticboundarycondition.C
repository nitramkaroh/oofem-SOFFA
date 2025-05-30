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
 *               Copyright (C) 1993 - 2020   Borek Patzak
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

#include "hardmagneticboundarycondition.h"
#include "classfactory.h"
#include "MagnetoelasticMaterials/mooneyrivlinhardmagnetic.h"
#include "sparsemtrx.h"
#include "feinterpol.h"
#include "magneticloadinterface.h"
#include "interfacetype.h"
#include "mathfem.h"
#include "function.h"
#include "tensor/tensor3.h"
#include "nlstructuralelement.h"
#include "feinterpol3d.h"
#include "feinterpol2d.h"
namespace oofem {
REGISTER_BoundaryCondition( HardMagneticBoundaryCondition );

void HardMagneticBoundaryCondition::initializeFrom( InputRecord &ir )
{
    ActiveBoundaryCondition::initializeFrom( ir );

    IR_GIVE_FIELD( ir, bcMode, _IFT_HardMagneticBoundaryCondition_mode );

    FloatArray b_temp, m_temp;
    IR_GIVE_FIELD( ir, b_temp, _IFT_HardMagneticBoundaryCondition_b_ext );
    b_app = FloatArrayF<3>( b_temp );

    if ( bcMode == 2 || bcMode == 3 ) {
        IR_GIVE_FIELD( ir, m_temp, _IFT_HardMagneticBoundaryCondition_m );
        m_app = FloatArrayF<3>( m_temp );
    }

    IR_GIVE_FIELD( ir, ltf_index, _IFT_HardMagneticBoundaryCondition_ltf );
    IR_GIVE_FIELD( ir, mltf_index, _IFT_HardMagneticBoundaryCondition_mltf );

    mu0 = BASE_VACUUM_PERMEABILITY_MU_0;
    IR_GIVE_OPTIONAL_FIELD( ir, mu0, _IFT_HardMagneticBoundaryCondition_mu_0 );

    if ( bcMode == 1 ) {
        this->evaluateFreeSpaceStress();
    }
    int pullBackType_temp;
    IR_GIVE_FIELD( ir, pullBackType_temp, _IFT_HardMagneticBoundaryCondition_pullBackType );

    this->pullBackType = (MooneyRivlinHardMagnetic::PullBackType)pullBackType_temp;
}

void HardMagneticBoundaryCondition::assemble( SparseMtrx &answer, TimeStep *tStep,
    CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale, void *lock )
{

    if ( type != TangentStiffnessMatrix ) {
        return;
    }

    FloatMatrix Ke;
    IntArray r_loc, c_loc, bNodes, dofIDs;

    dofIDs = domain->giveDefaultNodeDofIDArry();

    Set *set                   = this->giveDomain()->giveSet( this->set );
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e   = this->giveDomain()->giveElement( boundaries.at( pos * 2 - 1 ) );
        int boundary = boundaries.at( pos * 2 );
        if ( e->giveSpatialDimension() == 2 ) {
            bNodes = e->giveInterpolation()->boundaryEdgeGiveNodes( boundary );
        } else {
            bNodes = e->giveInterpolation()->boundarySurfaceGiveNodes( boundary );
        }

        e->giveBoundaryLocationArray( r_loc, bNodes, dofIDs /* this->dofs*/, r_s );
        e->giveBoundaryLocationArray( c_loc, bNodes, dofIDs /* this->dofs*/, c_s );
        this->computeTangentFromElement( Ke, e, boundary, tStep );

        // debug
        // compute numerical tangent
        /*
        FloatMatrix Ke_num;
        double pert = 1.e-6;
        this->computeNumericalTangentFromElement( Ke_num, e, boundary, tStep, pert );

        //symmetry check
        FloatMatrix KeT, Ke_numT, errorKe, errorKe_num, Ke_dif;
        KeT.beTranspositionOf( Ke );
        Ke_numT.beTranspositionOf( Ke_num );
        errorKe = KeT;
        errorKe_num = Ke_numT;
        errorKe.add( -1. * Ke );
        errorKe_num.add( -1. * Ke_num );

        Ke_dif = Ke - Ke_num;

        if ( pos == 0 ) {

            OOFEM_LOG_INFO( "Debugging surface %i\n", boundary );
            OOFEM_LOG_INFO( "-------Analytical stiffness-----------\n" );
            Ke.printYourself();
            OOFEM_LOG_INFO( "-------Numerical stiffness------------\n" );
            Ke_num.printYourself();
            OOFEM_LOG_INFO( "-------Difference---------------------\n" );
            Ke_dif.printYourself();
        }

        answer.assemble( r_loc, c_loc, Ke);
        */
    }
}

void HardMagneticBoundaryCondition::assembleVector( FloatArray &answer, TimeStep *tStep,
    CharType type, ValueModeType mode,
    const UnknownNumberingScheme &s, FloatArray *eNorms, void *lock )
{

    if ( type != InternalForcesVector ) {
        return;
    }

    FloatArray fe;
    IntArray loc, masterdofids, bNodes, dofIDs;
    dofIDs = domain->giveDefaultNodeDofIDArry();

    Set *set                   = this->giveDomain()->giveSet( this->set );
    const IntArray &boundaries = set->giveBoundaryList();
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e   = this->giveDomain()->giveElement( boundaries.at( pos * 2 - 1 ) );
        int boundary = boundaries.at( pos * 2 );
        if ( e->giveSpatialDimension() == 2 ) {
            bNodes = e->giveInterpolation()->boundaryEdgeGiveNodes( boundary );
        } else {
            bNodes = e->giveInterpolation()->boundarySurfaceGiveNodes( boundary );
        }
        e->giveBoundaryLocationArray( loc, bNodes, dofIDs /* this->dofs*/, s, &masterdofids );
        this->computeLoadVectorFromElement( fe, e, boundary, tStep, mode );
        answer.assemble( -1. * fe, loc );
        if ( eNorms ) {
            eNorms->assembleSquared( fe, masterdofids );
        }
    }
}

void HardMagneticBoundaryCondition::giveLocationArrays( std::vector<IntArray> &rows, std::vector<IntArray> &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s )
{
    if ( type != TangentStiffnessMatrix ) {
        return;
    }

    Set *set                   = this->giveDomain()->giveSet( this->set );
    const IntArray &boundaries = set->giveBoundaryList();
    IntArray bNodes, dofIDs;
    dofIDs = domain->giveDefaultNodeDofIDArry();

    rows.resize( boundaries.giveSize() / 2 );
    cols.resize( boundaries.giveSize() / 2 );

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e   = this->giveDomain()->giveElement( boundaries.at( pos * 2 - 1 ) );
        int boundary = boundaries.at( pos * 2 );

        bNodes = e->giveInterpolation()->boundaryGiveNodes( boundary );

        e->giveBoundaryLocationArray( rows[pos], bNodes, dofIDs /* this->dofs*/, r_s );
        e->giveBoundaryLocationArray( cols[pos], bNodes, dofIDs /* this->dofs*/, c_s );
    }
}

void HardMagneticBoundaryCondition::computeLoadVectorFromElement( FloatArray &answer, Element *e, int iSurf, TimeStep *tStep, ValueModeType mode )
{
    if ( iSurf == -1 ) {
        iSurf = 1;
    }

    switch ( bcMode ) {
    case 1:
        this->computeElementLoadVectorContribution_maxwell( answer, e, iSurf, tStep );
        break;
    case 2:
        this->computeElementLoadVectorContribution_bla( answer, e, iSurf, tStep );
        break;
    case 3:
        this->computeElementLoadVectorContribution_maxwell2( answer, e, iSurf, tStep );
        break;
    default:
        OOFEM_ERROR( "Unknown mode for HardMagneticBoundaryCondition (%i)", bcMode );
        break;
    }
}

void HardMagneticBoundaryCondition::computeTangentFromElement( FloatMatrix &answer, Element *e, int iSurf, TimeStep *tStep )
{

    if ( iSurf == -1 ) {
        iSurf = 1;
    }

    switch ( bcMode ) {
    case 1:
        this->computeElementTangentContribution_maxwell( answer, e, iSurf, tStep );
        break;
    case 2:
        this->computeElementTangentContribution_bla( answer, e, iSurf, tStep );
        break;
    case 3:
        this->computeElementTangentContribution_maxwell2( answer, e, iSurf, tStep );
        break;
    default:
        OOFEM_ERROR( "Unknown mode for HardMagneticBoundaryCondition (%i)", bcMode );
        break;
    }
}

void HardMagneticBoundaryCondition::evaluateFreeSpaceStress()
{
    auto b = Tensor1_3d( b_app );
    Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
    maxwell_stress( i_3, j_3 ) = 1 / mu0 * ( b( i_3 ) * b( j_3 ) - 0.5 * b( p_3 ) * b( p_3 ) * delta( i_3, j_3 ) );
}


//-------------------------------------------------------------------Maxwell approach-----------------------------------------------------------------------------

void HardMagneticBoundaryCondition::computeElementLoadVectorContribution_maxwell( FloatArray &answer, Element *e, int iSurf, TimeStep *tStep )
{

    NLStructuralElement *nle         = dynamic_cast<NLStructuralElement *>( e );
    FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>( e->giveInterpolation() );
    if ( nle == nullptr || interpolation == nullptr ) {
        OOFEM_ERROR( "Nonlinear elements required for magnetic BCs" );
    }

    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );

    answer.clear();

    int order  = 4;
    auto iRule = nle->giveBoundarySurfaceIntegrationRule( order, iSurf );

    Tensor1_3d maxwellStressCofFN;
    Tensor2_3d F;
    FloatArray vF;

    for ( GaussPoint *gp : *iRule ) {
        // compute normal vector
        auto [dA, n] = interpolation->surfaceEvalUnitNormal( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( e ) );
        // compute surface N and B matrix
        FloatMatrix N;
        nle->computeSurfaceNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
        // compute deformation gradient
        nle->computeDeformationGradientVectorAtBoundary( vF, gp, iSurf, tStep );

        // compute the force vector
        auto F                    = Tensor2_3d( FloatArrayF<9>( vF ) );
        auto Normal               = Tensor1_3d( FloatArrayF<3>( n ) );
        auto [J, cofF]            = F.compute_determinant_and_cofactor();
        maxwellStressCofFN( i_3 ) = load_level * load_level * maxwell_stress( i_3, j_3 ) * cofF( j_3, k_3 ) * Normal( k_3 );
        //
        answer.plusProduct( N, maxwellStressCofFN.to_voigt_form(), -gp->giveWeight() * dA );
    }
}

void HardMagneticBoundaryCondition::computeElementTangentContribution_maxwell( FloatMatrix &answer, Element *e, int iSurf, TimeStep *tStep )
{
    NLStructuralElement *nle         = dynamic_cast<NLStructuralElement *>( e );
    FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>( e->giveInterpolation() );
    if ( nle == nullptr || interpolation == nullptr ) {
        OOFEM_ERROR( "Nonlinear elements required for magnetic BCs" );
    }

    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );

    answer.clear();
    int order  = 4;
    auto iRule = nle->giveBoundarySurfaceIntegrationRule( order, iSurf );

    FloatMatrix N, B, Nt;

    if ( e->giveSpatialDimension() == 3 ) {
        Tensor3_3d maxwellStressFcrossN;
        FloatArray vF;

        for ( GaussPoint *gp : *iRule ) {
            // compute normal vector
            auto [dA, n] = interpolation->surfaceEvalUnitNormal( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( e ) );
            // compute surface N and B matrix
            nle->computeSurfaceNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
            nle->computeBHmatrixAtBoundary( gp, B, iSurf );
            // compute deformation gradient
            nle->computeDeformationGradientVectorAtBoundary( vF, gp, iSurf, tStep );
            // compute the force vector
            auto F         = Tensor2_3d( FloatArrayF<9>( vF ) );
            auto Normal    = Tensor1_3d( FloatArrayF<3>( n ) );
            auto [J, cofF] = F.compute_determinant_and_cofactor();
            //
            maxwellStressFcrossN( i_3, m_3, n_3 ) = load_level * load_level * maxwell_stress( i_3, j_3 ) * F.compute_tensor_cross_product()( j_3, k_3, m_3, n_3 ) * Normal( k_3 );
            //
            Nt.beTranspositionOf( N );
            answer += -gp->giveWeight() * dA * ( Nt * maxwellStressFcrossN.to_voigt_form_3x9() * B );
        }
    } else if ( e->giveSpatialDimension() == 2 ) {
        OOFEM_ERROR( "Magnetic boundary condition not implemented for 2D domains." );
    }
}

//-------------------------------------------------------------------Body load approach-----------------------------------------------------------------------------

void HardMagneticBoundaryCondition::computeElementLoadVectorContribution_bla( FloatArray &answer, Element *e, int iSurf, TimeStep *tStep )
{
    NLStructuralElement *nle         = dynamic_cast<NLStructuralElement *>( e );
    FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>( e->giveInterpolation() );
    if ( nle == nullptr || interpolation == nullptr ) {
        OOFEM_ERROR( "Nonlinear elements required for magnetic BCs" );
    }

    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );

    answer.clear();

    int order  = 4;
    auto iRule = nle->giveBoundarySurfaceIntegrationRule( order, iSurf );

    Tensor1_3d contribution;
    auto B_app = Tensor1_3d( load_level * b_app );
    auto M     = Tensor1_3d( m_app );
    FloatArray vF;
    FloatMatrix N;

    for ( GaussPoint *gp : *iRule ) {
        // compute normal vector
        auto [dA, n] = interpolation->surfaceEvalUnitNormal( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( e ) );
        // compute surface N matrix
        nle->computeSurfaceNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
        // compute deformation gradient
        nle->computeDeformationGradientVectorAtBoundary( vF, gp, iSurf, tStep );

        // compute the force vector
        auto F         = Tensor2_3d( FloatArrayF<9>( vF ) );
        auto Normal    = Tensor1_3d( FloatArrayF<3>( n ) );
        auto [J, cofF] = F.compute_determinant_and_cofactor();

        contribution( i_3 ) = M( m_3 ) * Normal( m_3 ) * B_app( i_3 )
            + ( mu0 / 2. ) * Normal( k_3 ) * M( k_3 ) * Normal( l_3 ) * M( l_3 ) * cofF( i_3, m_3 ) * Normal( m_3 )
            - ( mu0 / ( 2. * J * J ) ) * F( k_3, p_3 ) * M( p_3 ) * F( k_3, q_3 ) * M( q_3 ) * cofF( i_3, m_3 ) * Normal( m_3 );
        //
        answer.plusProduct( N, contribution.to_voigt_form(), gp->giveWeight() * dA );
    }
}

void HardMagneticBoundaryCondition::computeElementTangentContribution_bla( FloatMatrix &answer, Element *e, int iSurf, TimeStep *tStep )
{
    NLStructuralElement *nle         = dynamic_cast<NLStructuralElement *>( e );
    FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>( e->giveInterpolation() );
    if ( nle == nullptr || interpolation == nullptr ) {
        OOFEM_ERROR( "Nonlinear elements required for magnetic BCs" );
    }

    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );

    answer.clear();
    int order  = 4;
    auto iRule = nle->giveBoundarySurfaceIntegrationRule( order, iSurf );

    Tensor3_3d contribution;
    auto B_app = Tensor1_3d( load_level * b_app );
    auto M     = Tensor1_3d( m_app );
    FloatArray vF;
    FloatMatrix N, B, Nt;

    if ( e->giveSpatialDimension() == 3 ) {

        for ( GaussPoint *gp : *iRule ) {
            // compute normal vector
            auto [dA, n] = interpolation->surfaceEvalUnitNormal( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( e ) );
            // compute surface N and B matrix
            nle->computeSurfaceNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
            nle->computeBHmatrixAtBoundary( gp, B, iSurf );
            // compute deformation gradient
            nle->computeDeformationGradientVectorAtBoundary( vF, gp, iSurf, tStep );
            // compute the force vector
            auto F         = Tensor2_3d( FloatArrayF<9>( vF ) );
            auto Normal    = Tensor1_3d( FloatArrayF<3>( n ) );
            auto [J, cofF] = F.compute_determinant_and_cofactor();
            auto Fcross    = F.compute_tensor_cross_product();
            //
            contribution( i_3, k_3, l_3 ) = ( mu0 / 2. ) * Normal( p_3 ) * M( p_3 ) * Normal( q_3 ) * Normal( q_3 ) * Fcross( i_3, m_3, k_3, l_3 ) * Normal( m_3 )
                + ( mu0 / ( J * J * J ) ) * F( j_3, p_3 ) * M( p_3 ) * F( j_3, q_3 ) * M( q_3 ) * cofF( i_3, m_3 ) * cofF( k_3, l_3 ) * Normal( m_3 )
                - ( mu0 / ( J * J ) ) * M( l_3 ) * F( k_3, q_3 ) * M( q_3 ) * cofF( i_3, m_3 ) * Normal( m_3 )
                - ( mu0 / ( 2. * J * J ) ) * F( j_3, p_3 ) * M( p_3 ) * F( j_3, q_3 ) * M( q_3 ) * Fcross( i_3, m_3, k_3, l_3 ) * Normal( m_3 );
            //
            Nt.beTranspositionOf( N );
            answer += gp->giveWeight() * dA * ( Nt * contribution.to_voigt_form_3x9() * B );
        }
    } else if ( e->giveSpatialDimension() == 2 ) {
        OOFEM_ERROR( "Magnetic boundary condition not implemented for 2D domains." );
    }
}

//-------------------------------------------------------------------Maxwell 2 approach-----------------------------------------------------------------------------

void HardMagneticBoundaryCondition::computeElementLoadVectorContribution_maxwell2( FloatArray &answer, Element *e, int iSurf, TimeStep *tStep )
{
    NLStructuralElement *nle = dynamic_cast<NLStructuralElement *>( e );
    // FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>( e->giveInterpolation() );
    FEInterpolation2d *interpolation = dynamic_cast<FEInterpolation2d *>( e->giveInterpolation() );
    if ( nle == nullptr || interpolation == nullptr ) {
        OOFEM_ERROR( "Nonlinear elements required for magnetic BCs" );
    }

    double load_level  = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    double mload_level = this->giveDomain()->giveFunction( mltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );

    answer.clear();

    int order = 1;
    // auto iRule = nle->giveBoundarySurfaceIntegrationRule( order, iSurf );
    auto iRule = nle->giveBoundaryEdgeIntegrationRule( order, iSurf );
    Tensor1_3d contribution, B_ref;
    auto ba = Tensor1_3d( load_level * b_app );
    auto M = Tensor1_3d( mload_level * m_app );
    FloatArray vF, vFred;
    FloatMatrix N;
    Tensor1_3d t;
    for ( GaussPoint *gp : *iRule ) {
        // compute normal vector
        auto [dA, nred]  = interpolation->surfaceEvalUnitNormal( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( e ) );
        FloatArrayF<3> n = { nred.at( 1 ), nred.at( 2 ), 0.0 };
        // compute surface N matrix
        // nle->computeSurfaceNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
        nle->computeEdgeNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
        // compute deformation gradient
        FloatArray coords;
        this->giveGaussPointCoordinates( coords, gp, iSurf );
        // correct this
        auto gp0 = std::make_unique<GaussPoint>( nullptr, 1, coords, gp->giveWeight(), _PlaneStrain );
        nle->computeDeformationGradientVector( vFred, gp0.get(), tStep );
        // nle->computeDeformationGradientVectorAtBoundary( vF, gp, iSurf, tStep );
        StructuralMaterial::giveFullVectorForm( vF, vFred, _PlaneStrain );
        // compute the force vector auxilliaries
        auto F         = Tensor2_3d( FloatArrayF<9>( vF ) );
        auto Normal    = Tensor1_3d( FloatArrayF<3>( n ) );
        auto [J, cofF] = F.compute_determinant_and_cofactor();
        // compute the actual contribution
        Tensor1_3d normal, m;
        normal( i_3 ) = cofF( i_3, j_3 ) * Normal( j_3 );
        auto norm2    = ( normal( i_3 ) * normal( i_3 ) );
        m             = this->giveActualMagnetization( F, M );
	//
	auto b = ba;
	b(i_3) += mu0 * m(i_3);
        //
        t( i_3 ) = 1. / mu0 * normal( k_3 ) * b( k_3 ) * b( i_3 ) - normal( k_3 ) * b( k_3 ) * m( i_3 );
	//+ 0.5 * mu0 * pow((m(k_3) * normal(k_3)),2)/norm2 * normal(i_3) - 0.5 / mu0 * b( k_3 ) * b( k_3 ) * normal( i_3 );
		     //+ b( k_3 ) * m( k_3 ) * normal( i_3 ) - 0.5 / mu0 * b( k_3 ) * b( k_3 ) * normal( i_3 ) )  + 0.5* mu0 * (pow((m(k_3) * normal(k_3)),2)/norm2- m(k_3)*m(k_3)) * normal(i_3);
	FloatArray t2d;
	IntArray indx = {1,2};
	t2d.beSubArrayOf(t.to_voigt_form(),indx);
	//
        answer.plusProduct( N, t2d, gp->giveWeight() * dA );
    }
}

Tensor1_3d HardMagneticBoundaryCondition ::giveActualMagnetization( const Tensor2_3d &F, const Tensor1_3d &M )
{
    Tensor1_3d m;
    Tensor2_3d iF;
    Tensor2_3d C, U, invU, R;
    auto J = F.compute_determinant();
    switch ( this->pullBackType ) {
    case MooneyRivlinHardMagnetic::PullBackType::PBT_F:
        m( i_3 ) = 1. / J * F( i_3, k_3 ) * M( k_3 );
        break;
    case MooneyRivlinHardMagnetic::PullBackType::PBT_R:
        C( i_3, j_3 ) = F( k_3, i_3 ) * F( k_3, j_3 );
        //
        U( k_3, l_3 ) = C.computeTensorPower( 0.5 )( k_3, l_3 );
        invU          = U.compute_inverse();
        R( i_3, j_3 ) = F( i_3, k_3 ) * invU( k_3, j_3 );
        //
        m( i_3 ) = 1. / J * R( i_3, k_3 ) * M( k_3 );
        break;
    case MooneyRivlinHardMagnetic::PullBackType::PBT_iFt: // theoretically wrong -MA
        iF       = F.compute_inverse();
        m( i_3 ) = 1. / J * iF( k_3, i_3 ) * M( k_3 );
        break;
    case MooneyRivlinHardMagnetic::PullBackType::PBT_iFtnoJ:
        iF       = F.compute_inverse();
        m( i_3 ) = iF( k_3, i_3 ) * M( k_3 );
        break;
    default:
        OOFEM_ERROR( "Unknown pullback type." );
        break;
    }

    return m;
}

void HardMagneticBoundaryCondition ::giveGaussPointCoordinates( FloatArray &coords, GaussPoint *gp, int iSurf )
{
    coords.resize( 2 );
    //
    if ( iSurf == 1 ) {
        coords.at( 1 ) = gp->giveNaturalCoordinates().at( 1 );
        coords.at( 2 ) = -1;
    } else if ( iSurf == 2 ) {
        coords.at( 1 ) = 1;
        coords.at( 2 ) = gp->giveNaturalCoordinates().at( 1 );
    } else if ( iSurf == 3 ) {
        coords.at( 1 ) = gp->giveNaturalCoordinates().at( 1 );
        coords.at( 2 ) = 1;
    } else if ( iSurf == 4 ) {
        coords.at( 1 ) = -1;
        coords.at( 2 ) = gp->giveNaturalCoordinates().at( 1 );
    }
}

void HardMagneticBoundaryCondition::computeElementTangentContribution_maxwell2( FloatMatrix &answer, Element *e, int iSurf, TimeStep *tStep )
{
  
    NLStructuralElement *nle = dynamic_cast<NLStructuralElement *>( e );
    // FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>( e->giveInterpolation() );
    FEInterpolation2d *interpolation = dynamic_cast<FEInterpolation2d *>( e->giveInterpolation() );
    if ( nle == nullptr || interpolation == nullptr ) {
        OOFEM_ERROR( "Nonlinear elements required for magnetic BCs" );
    }

    double load_level  = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
    double mload_level = this->giveDomain()->giveFunction( mltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );

    answer.clear();

    int order = 1;
    // auto iRule = nle->giveBoundarySurfaceIntegrationRule( order, iSurf );
    auto iRule = nle->giveBoundaryEdgeIntegrationRule( order, iSurf );
    Tensor1_3d contribution, B_ref;
    auto b = Tensor1_3d( load_level * b_app );
    auto M = Tensor1_3d( mload_level * m_app );
    FloatArray vF, vFred;
    FloatMatrix N, B, Bp, Ndt;
    Tensor1_3d t;
    for ( GaussPoint *gp : *iRule ) {
        // compute normal vector
        auto [dA, nred]  = interpolation->surfaceEvalUnitNormal( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( e ) );
        FloatArrayF<3> n = { nred.at( 1 ), nred.at( 2 ), 0.0 };
        // compute surface N matrix
        // nle->computeSurfaceNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
        nle->computeEdgeNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
        // compute deformation gradient
        FloatArray coords;
        this->giveGaussPointCoordinates( coords, gp, iSurf );
        // correct this
        auto gp0 = std::make_unique<GaussPoint>( nullptr, 1, coords, gp->giveWeight(), _PlaneStrain );
        nle->computeDeformationGradientVector( vFred, gp0.get(), tStep );
        // nle->computeDeformationGradientVectorAtBoundary( vF, gp, iSurf, tStep );
        StructuralMaterial::giveFullVectorForm( vF, vFred, _PlaneStrain );
	//
	nle->computeBHmatrixAtBoundary(gp0.get(), B, iSurf);
        // compute the force vector auxilliaries
        auto F         = Tensor2_3d( FloatArrayF<9>( vF ) );
        auto Normal    = Tensor1_3d( FloatArrayF<3>( n ) );
        auto [J, cofF] = F.compute_determinant_and_cofactor();
        // compute the actual contribution
        Tensor1_3d normal, m;
        normal( i_3 ) = cofF( i_3, j_3 ) * Normal( j_3 );
        m             = this->giveActualMagnetization( F, M );
	Tensor1_3d bapp;
	bapp = b;
	b(i_3) += mu0 * m(i_3);
	//
	Tensor2_3d iC;
	iC(i_3,j_3) = 1./J/J * cofF(k_3,i_3) * cofF(k_3,j_3);
	Tensor4_3d Fcross;
	Fcross(i_3, j_3, k_3, l_3) = F.compute_tensor_cross_product()(i_3, j_3, k_3, l_3);
	Tensor4_3d d_iC_dF;
	d_iC_dF(i_3, j_3, k_3, l_3) = -F(i_3,p_3) *  (iC(p_3,k_3) * iC(l_3, j_3) + iC(p_3,l_3) * iC(k_3, j_3));
	Tensor3_3d dtdF;
	dtdF(i_3,k_3,l_3) = 1. / mu0 * bapp(m_3) * Fcross(m_3,n_3,k_3,l_3) * Normal(n_3) * bapp(i_3) + (Normal(m_3) * iC(m_3,n_3) * M(n_3)) * bapp(i_3)*cofF(k_3,l_3) + J * Normal(m_3) * d_iC_dF(m_3,n_3,k_3,l_3) * M(n_3) * bapp(i_3);
	//
	FloatMatrix dtdF_3d(dtdF.to_voigt_form_3x9());
	FloatMatrix dtdF_2d;
	dtdF_2d.beSubMatrixOf(dtdF_3d, 1, 2, 1, 5);
	//
	Ndt.beTProductOf( N, dtdF_2d);
	FloatMatrix NtB;
	NtB.beProductOf(Ndt,B);
	answer.add( -1. * gp->giveWeight() * dA, NtB);
	//answer.plusProductUnsym( Ndt	, B,  -1. * gp->giveWeight() * dA);
	
    }
  
}


//----------------------------------------------------------------------DEBUG-------------------------------------------------------------------------------------

void HardMagneticBoundaryCondition::computeNumericalTangentFromElement( FloatMatrix &answer, Element *e, int side, TimeStep *tStep, double perturb )
{
    // debugging, numerically computing stiffness
    answer.resize( 12, 12 );

    FloatArray perturbedForceFront, perturbedForceBack;
    for ( int i = 1; i <= 12; i++ ) {

        this->computePerturbedLoadVectorFromElement( perturbedForceFront, e, side, tStep, perturb, i );
        this->computePerturbedLoadVectorFromElement( perturbedForceBack, e, side, tStep, -perturb, i );

        for ( int j = 1; j <= 12; j++ ) {

            answer.at( i, j ) = ( perturbedForceFront.at( j ) - perturbedForceBack.at( j ) ) / ( 2 * perturb );
        }
    }
}

void HardMagneticBoundaryCondition::computePerturbedLoadVectorFromElement( FloatArray &answer, Element *e, int iSurf, TimeStep *tStep, double perturb, int index )
{
    if ( iSurf == -1 ) {
        iSurf = 1;
    }

    switch ( bcMode ) {
    case 1:
        this->computePerturbedLoadVectorFromElement_maxwell( answer, e, iSurf, tStep, perturb, index );
        break;
    case 2:
        this->computePerturbedLoadVectorFromElement_bla( answer, e, iSurf, tStep, perturb, index );
        break;
    case 3:
        this->computePerturbedLoadVectorFromElement_maxwell2( answer, e, iSurf, tStep, perturb, index );
        break;
    default:
        OOFEM_ERROR( "Unknown mode for HardMagneticBoundaryCondition (%i)", bcMode );
        break;
    }
}

void HardMagneticBoundaryCondition::computePerturbedLoadVectorFromElement_bla( FloatArray &answer, Element *e, int iSurf, TimeStep *tStep, double perturb, int index )
{
    // debugging, numerically perturbed load vector
    IntArray bNodes = e->giveInterpolation()->boundarySurfaceGiveNodes( iSurf );

    int ndof = 3 * bNodes.giveSize();
    answer.resize( ndof );

    // standard load vector computation... almost:
    NLStructuralElement *nle         = dynamic_cast<NLStructuralElement *>( e );
    FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>( e->giveInterpolation() );

    if ( nle == nullptr || interpolation == nullptr ) {
        OOFEM_ERROR( "Nonlinear elements required for magnetic BCs" );
    }

    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );


    answer.clear();

    int order  = 4;
    auto iRule = nle->giveBoundarySurfaceIntegrationRule( order, iSurf );

    Tensor1_3d contribution;
    auto B_app = Tensor1_3d( load_level * b_app );
    auto M     = Tensor1_3d( m_app );
    Tensor2_3d F;
    FloatArray vF, uPert;

    for ( GaussPoint *gp : *iRule ) {
        // compute normal vector
        auto [dA, n] = interpolation->surfaceEvalUnitNormal( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( e ) );
        // compute surface N matirx
        FloatMatrix N, B;
        nle->computeSurfaceNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
        nle->computeBHmatrixAtBoundary( gp, B, iSurf );
        // compute deformation gradient
        // do it from perturbed displacements
        nle->computeBoundaryVectorOf( bNodes, { D_u, D_v, D_w }, VM_Total, tStep, uPert );
        uPert.at( index ) += perturb;
        vF.beProductOf( B, uPert );
        vF.at( 1 ) += 1.0;
        vF.at( 2 ) += 1.0;
        vF.at( 3 ) += 1.0;

        // compute the force vector
        auto F              = Tensor2_3d( FloatArrayF<9>( vF ) );
        auto Normal         = Tensor1_3d( FloatArrayF<3>( n ) );
        auto [J, cofF]      = F.compute_determinant_and_cofactor();
        contribution( i_3 ) = M( m_3 ) * Normal( m_3 ) * B_app( i_3 )
            + ( mu0 / 2. ) * Normal( k_3 ) * M( k_3 ) * Normal( l_3 ) * M( l_3 ) * cofF( i_3, m_3 ) * Normal( m_3 )
            - ( mu0 / ( 2. * J * J ) ) * F( k_3, p_3 ) * M( p_3 ) * F( k_3, q_3 ) * M( q_3 ) * cofF( i_3, m_3 ) * Normal( m_3 );
        //
        answer.plusProduct( N, contribution.to_voigt_form(), gp->giveWeight() * dA );
    }
}

void HardMagneticBoundaryCondition::computePerturbedLoadVectorFromElement_maxwell2( FloatArray &answer, Element *e, int iSurf, TimeStep *tStep, double perturb, int index )
{
    // debugging, numerically perturbed load vector
    IntArray bNodes = e->giveInterpolation()->boundarySurfaceGiveNodes( iSurf );

    int ndof = 3 * bNodes.giveSize();
    answer.resize( ndof );

    // standard load vector computation... almost:
    NLStructuralElement *nle         = dynamic_cast<NLStructuralElement *>( e );
    FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>( e->giveInterpolation() );

    if ( nle == nullptr || interpolation == nullptr ) {
        OOFEM_ERROR( "Nonlinear elements required for magnetic BCs" );
    }

    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );


    answer.clear();

    int order  = 4;
    auto iRule = nle->giveBoundarySurfaceIntegrationRule( order, iSurf );

    Tensor1_3d contribution, B_ref;
    auto B_app = Tensor1_3d( load_level * b_app );
    auto M     = Tensor1_3d( m_app );
    Tensor2_3d F;
    FloatArray vF, uPert;

    for ( GaussPoint *gp : *iRule ) {
        // compute normal vector
        auto [dA, n] = interpolation->surfaceEvalUnitNormal( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( e ) );
        // compute surface N matirx
        FloatMatrix N, B;
        nle->computeSurfaceNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
        nle->computeBHmatrixAtBoundary( gp, B, iSurf );
        // compute deformation gradient
        // do it from perturbed displacements
        nle->computeBoundaryVectorOf( bNodes, { D_u, D_v, D_w }, VM_Total, tStep, uPert );
        uPert.at( index ) += perturb;
        vF.beProductOf( B, uPert );
        vF.at( 1 ) += 1.0;
        vF.at( 2 ) += 1.0;
        vF.at( 3 ) += 1.0;

        // compute the force vector
        auto F         = Tensor2_3d( FloatArrayF<9>( vF ) );
        auto Normal    = Tensor1_3d( FloatArrayF<3>( n ) );
        auto [J, cofF] = F.compute_determinant_and_cofactor();

        // convert B to reference value
        B_ref( i_3 ) = J * F.compute_inverse()( i_3, j_3 ) * B_app( j_3 );

        // compute the actual contribution
        contribution( i_3 ) = 1. / ( mu0 * J ) * Normal( k_3 ) * B_ref( k_3 ) * F( i_3, m_3 ) * B_ref( m_3 )
            - ( 1. / J ) * Normal( k_3 ) * B_ref( k_3 ) * F( i_3, m_3 ) * M( m_3 )
            - ( 1. / ( 2. * mu0 * J * J ) ) * F( k_3, p_3 ) * B_ref( p_3 ) * F( k_3, q_3 ) * B_ref( q_3 ) * cofF( i_3, m_3 ) * Normal( m_3 )
            + ( 1. / J * J ) * F( k_3, p_3 ) * B_ref( p_3 ) * F( k_3, q_3 ) * M( q_3 ) * cofF( i_3, m_3 ) * Normal( m_3 )
            - ( mu0 / ( 2. * J * J ) ) * F( k_3, p_3 ) * M( p_3 ) * F( k_3, q_3 ) * M( q_3 ) * cofF( i_3, m_3 ) * Normal( m_3 )
            + ( mu0 / 2. ) * M( k_3 ) * Normal( k_3 ) * M( l_3 ) * Normal( l_3 ) * cofF( i_3, m_3 ) * Normal( m_3 );
        //
        answer.plusProduct( N, contribution.to_voigt_form(), -gp->giveWeight() * dA );
    }
}

void HardMagneticBoundaryCondition::computePerturbedLoadVectorFromElement_maxwell( FloatArray &answer, Element *e, int iSurf, TimeStep *tStep, double perturb, int index )
{
    // debugging, numerically perturbed load vector
    IntArray bNodes = e->giveInterpolation()->boundarySurfaceGiveNodes( iSurf );

    int ndof = 3 * bNodes.giveSize();
    answer.resize( ndof );

    // standard load vector computation... almost:
    NLStructuralElement *nle         = dynamic_cast<NLStructuralElement *>( e );
    FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>( e->giveInterpolation() );

    if ( nle == nullptr || interpolation == nullptr ) {
        OOFEM_ERROR( "Nonlinear elements required for magnetic BCs" );
    }

    double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );


    answer.clear();

    int order  = 4;
    auto iRule = nle->giveBoundarySurfaceIntegrationRule( order, iSurf );

    Tensor1_3d maxwellStressCofFN;

    Tensor2_3d F;
    FloatArray vF, uPert;

    for ( GaussPoint *gp : *iRule ) {
        // compute normal vector
        auto [dA, n] = interpolation->surfaceEvalUnitNormal( iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( e ) );
        // compute surface N matirx
        FloatMatrix N, B;
        nle->computeSurfaceNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
        nle->computeBHmatrixAtBoundary( gp, B, iSurf );
        // compute deformation gradient
        // do it from perturbed displacements
        nle->computeBoundaryVectorOf( bNodes, { D_u, D_v, D_w }, VM_Total, tStep, uPert );
        uPert.at( index ) += perturb;
        vF.beProductOf( B, uPert );
        vF.at( 1 ) += 1.0;
        vF.at( 2 ) += 1.0;
        vF.at( 3 ) += 1.0;

        // compute the force vector
        auto F                    = Tensor2_3d( FloatArrayF<9>( vF ) );
        auto Normal               = Tensor1_3d( FloatArrayF<3>( n ) );
        auto [J, cofF]            = F.compute_determinant_and_cofactor();
        maxwellStressCofFN( i_3 ) = load_level * load_level * maxwell_stress( i_3, j_3 ) * cofF( j_3, k_3 ) * Normal( k_3 );
        //
        answer.plusProduct( N, maxwellStressCofFN.to_voigt_form(), -gp->giveWeight() * dA );
    }
}

} // namespace oofem
