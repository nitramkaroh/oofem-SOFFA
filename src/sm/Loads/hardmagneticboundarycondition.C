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
#include "sparsemtrx.h"
#include "feinterpol.h"
#include "mathfem.h"
#include "function.h"
#include "tensor/tensor3.h"
#include "nlstructuralelement.h"
#include "feinterpol3d.h"
#include "feinterpol2d.h"
#include "floatarrayf.h"
namespace oofem {
REGISTER_BoundaryCondition( HardMagneticEdgeBoundaryCondition );
REGISTER_BoundaryCondition( HardMagneticSurfaceBoundaryCondition );

void HardMagneticBoundaryCondition::initializeFrom( InputRecord &ir )
{
    ActiveBoundaryCondition::initializeFrom( ir );
    int matNum;
    IR_GIVE_FIELD( ir, matNum, _IFT_HardMagneticBoundaryCondition_reducedHardMagneticMaterialNumber );
    auto m = domain->giveMaterial(matNum);
    // Attempt to cast to ReducedHardMagneticMaterial
    hmMat = dynamic_cast<ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent*>(m);

    if (hmMat == nullptr) {
      OOFEM_ERROR("Material number %d is not a ReducedHardMagneticMaterial.", matNum);
    }
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
	bNodes = e->giveInterpolation()->boundaryGiveNodes( boundary, e->giveGeometryType());
	//
        e->giveBoundaryLocationArray( r_loc, bNodes, dofIDs /* this->dofs*/, r_s );
        e->giveBoundaryLocationArray( c_loc, bNodes, dofIDs /* this->dofs*/, c_s );
        this->computeTangentFromElement( Ke, e, boundary, tStep );
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
	bNodes = e->giveInterpolation()->boundaryGiveNodes( boundary, e->giveGeometryType() );
	//        
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

        bNodes = e->giveInterpolation()->boundaryGiveNodes( boundary, e->giveGeometryType() );

        e->giveBoundaryLocationArray( rows[pos], bNodes, dofIDs /* this->dofs*/, r_s );
        e->giveBoundaryLocationArray( cols[pos], bNodes, dofIDs /* this->dofs*/, c_s );
    }
}



//-------------------------------------------------------------------Edge----------------------------------------------------------------------------
void HardMagneticEdgeBoundaryCondition ::giveGaussPointCoordinates( FloatArray &coords, GaussPoint *gp, int iSurf )
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



  
  void HardMagneticEdgeBoundaryCondition::computeLoadVectorFromElement( FloatArray &answer, Element *e, int iSurf, TimeStep *tStep, ValueModeType mode )
{
    NLStructuralElement *nle = dynamic_cast<NLStructuralElement *>( e );
    FEInterpolation *interpolation =  e->giveInterpolation();
    if ( nle == nullptr || interpolation == nullptr ) {
        OOFEM_ERROR( "Wrong element or interpolation for edge magnetic BCs" );
    }
    answer.clear();

    int order = 1;
    auto iRule = nle->giveBoundaryEdgeIntegrationRule( order, iSurf );
    FloatMatrix N;
    Tensor1_3d t;
    for ( GaussPoint *gp : *iRule ) {
        // compute normal vector
        FloatArray nred;
        auto dA =  interpolation->boundaryEdgeEvalNormal(nred, iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( e ));
	FloatArrayF<3> n = { nred.at( 1 ), nred.at( 2 ), 0.0 };
	//
	nle->computeEdgeNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
	// surface deformation gradient
	 // compute deformation gradient
        FloatArray coords;
        this->giveGaussPointCoordinates( coords, gp, iSurf );
        // correct this
        auto gp0 = std::make_unique<GaussPoint>( nullptr, 1, coords, gp->giveWeight(), _PlaneStrain );
	
	FloatArray vFred, vF;
        nle->computeDeformationGradientVector( vFred, gp0.get(), tStep );
        // nle->computeDeformationGradientVectorAtBoundary( vF, gp, iSurf, tStep );
        StructuralMaterial::giveFullVectorForm( vF, vFred, _PlaneStrain );
	auto F = Tensor2_3d(FloatArrayF<9> (vF));
	// nle->computeSurfaceDeformationGradientVector( vFred, gp, tStep, iSurf );
	//auto F = Tensor2_3d(oofem::assemble< 9 >(static_cast<const FloatArrayF<5>&>(vFred), { 0, 1, 2, 5, 8 }));
        auto [J, cofF] = F.compute_determinant_and_cofactor();
	// compute the actual unit normal
	auto Normal    = Tensor1_3d( FloatArrayF<3>( n ) );
        Tensor1_3d normalb, normal;
        normalb( i_3 ) = cofF( i_3, j_3 ) * Normal( j_3 );
	auto norm = normalb( k_3 ) * normalb( k_3 );
	normal( i_3 ) = 1./norm * normalb( i_3 );  
	// ask material for magnetization and magnetic induction
	auto [m, b] = hmMat->give_EulerianMagnetization_MagneticInduction_FromF(F, tStep );
        // compute the actual contribution
	t( i_3 ) = 1. / mu_0 * normal( k_3 ) * b( k_3 )* (b(i_3) - mu_0 * m( i_3 ) );
	//+ 0.5 * mu0 * pow((m(k_3) * normal(k_3)),2)/norm2 * normal(i_3) - 0.5 / mu0 * b( k_3 ) * b( k_3 ) * normal( i_3 );
		     //+ b( k_3 ) * m( k_3 ) * normal( i_3 ) - 0.5 / mu0 * b( k_3 ) * b( k_3 ) * normal( i_3 ) )  + 0.5* mu0 * (pow((m(k_3) * normal(k_3)),2)/norm2- m(k_3)*m(k_3)) * normal(i_3);

	FloatArray vt;
	vt.beSubArrayOf(t.to_voigt_form(),{1,2});
	answer.plusProduct( N, vt, gp->giveWeight() * dA );
    }
}



void HardMagneticEdgeBoundaryCondition::computeTangentFromElement( FloatMatrix &answer, Element *e, int iSurf, TimeStep *tStep )
{
  
    NLStructuralElement *nle = dynamic_cast<NLStructuralElement *>( e );
    // FEInterpolation3d *interpolation = dynamic_cast<FEInterpolation3d *>( e->giveInterpolation() );
    FEInterpolation2d *interpolation = dynamic_cast<FEInterpolation2d *>( e->giveInterpolation() );
    if ( nle == nullptr || interpolation == nullptr ) {
        OOFEM_ERROR( "Nonlinear elements required for magnetic BCs" );
    }

    answer.clear();

    int order = 1;
    // auto iRule = nle->giveBoundarySurfaceIntegrationRule( order, iSurf );
    auto iRule = nle->giveBoundaryEdgeIntegrationRule( order, iSurf );
    Tensor1_3d contribution, B_ref;

    FloatMatrix N, B, Bp, Ndt;
    Tensor1_3d t;
    for ( GaussPoint *gp : *iRule ) {
        // compute normal vector
        FloatArray nred;
        auto dA =  interpolation->boundaryEdgeEvalNormal(nred, iSurf, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( e ));
	FloatArrayF<3> n = { nred.at( 1 ), nred.at( 2 ), 0.0 };
	//
	nle->computeEdgeNMatrix( N, iSurf, gp->giveNaturalCoordinates() );
	// surface deformation gradient
	FloatArray vFred;
        nle->computeSurfaceDeformationGradientVector( vFred, gp, tStep, iSurf );
	//
	auto F = Tensor2_3d(oofem::assemble< 9 >(static_cast<const FloatArrayF<5>&>(vFred), { 0, 1, 2, 5, 8 }));
	//
	auto [J, cofF] = F.compute_determinant_and_cofactor();
	// compute auxiliary variables
	Tensor2_3d iC;
	iC(i_3,j_3) = 1./J/J * cofF(k_3,i_3) * cofF(k_3,j_3);
	Tensor4_3d Fcross;
	Fcross(i_3, j_3, k_3, l_3) = F.compute_tensor_cross_product()(i_3, j_3, k_3, l_3);
	Tensor4_3d d_iC_dF;
	d_iC_dF(i_3, j_3, k_3, l_3) = -F(i_3,p_3) *  (iC(p_3,k_3) * iC(l_3, j_3) + iC(p_3,l_3) * iC(k_3, j_3));
	// compute the actual unit normal
	auto Normal    = Tensor1_3d( FloatArrayF<3>( n ) );
        Tensor1_3d normal;
        normal( i_3 ) = cofF( i_3, j_3 ) * Normal( j_3 );
	// ask material for magnetization and magnetic induction
	auto [m, b] = hmMat->give_EulerianMagnetization_MagneticInduction_FromF(F, tStep );
	auto M = hmMat->giveLagrangianMagnetization(tStep);
	//
	Tensor1_3d bapp;
	bapp(m_3) = (b(m_3)-mu_0 * m(m_3));
	//
	Tensor3_3d dtdF;
	dtdF(i_3,k_3,l_3) = 1. / mu_0 * bapp(m_3) * Fcross(m_3,n_3,k_3,l_3) * Normal(n_3) * bapp(i_3) + (Normal(m_3) * iC(m_3,n_3) * M(n_3)) * bapp(i_3)*cofF(k_3,l_3) + J * Normal(m_3) * d_iC_dF(m_3,n_3,k_3,l_3) * M(n_3) * bapp(i_3);
	//	dtdF(i_3,k_3,l_3) = 1. / mu_00 * b(m_3) * Fcross(m_3,n_3,k_3,l_3) * Normal(n_3) * b(i_3) + (Normal(m_3) * iC(m_3,n_3) * M(n_3)) * b(i_3)*cofF(k_3,l_3) + J * Normal(m_3) * d_iC_dF(m_3,n_3,k_3,l_3) * M(n_3) * b(i_3);
	//
	FloatMatrix dtdF_3d(dtdF.to_voigt_form_3x9());
	FloatMatrix dtdF_2d;
	dtdF_2d.beSubMatrixOf(dtdF_3d, 1, 2, 1, 5);
	//
	nle->computeSurfaceBHmatrixAt(gp, B, iSurf);
	//
	Ndt.beTProductOf( N, dtdF_2d);
	FloatMatrix NtB;
	NtB.beProductOf(Ndt,B);
	//answer.add( -gp->giveWeight() * dA, NtB);
	
    }
  
}




} // namespace oofem
