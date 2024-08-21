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


namespace oofem {
REGISTER_BoundaryCondition( HardMagneticBoundaryCondition );

    void HardMagneticBoundaryCondition::initializeFrom( InputRecord &ir )
    {
        ActiveBoundaryCondition::initializeFrom( ir );

        IR_GIVE_FIELD( ir, b_ext, _IFT_HardMagneticBoundaryCondition_b_ext );
        IR_GIVE_FIELD( ir, ltf_index, _IFT_HardMagneticBoundaryCondition_ltf );

        mu0 = BASE_VACUUM_PERMEABILITY_MU_0;
        IR_GIVE_OPTIONAL_FIELD( ir, mu0, _IFT_HardMagneticBoundaryCondition_mu_0 );

        this->evaluateFreeSpaceStress();
    }

    void HardMagneticBoundaryCondition::assemble( SparseMtrx &answer, TimeStep *tStep,
        CharType type , const UnknownNumberingScheme& r_s , const UnknownNumberingScheme& c_s , double scale , void* lock) {
        
        if ( type != TangentStiffnessMatrix ) {
            return;
        }

        FloatMatrix Ke;
        IntArray r_loc, c_loc, bNodes;

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

            e->giveBoundaryLocationArray( r_loc, bNodes, this->dofs, r_s );
            e->giveBoundaryLocationArray( c_loc, bNodes, this->dofs, c_s );
            this->computeTangentFromElement( Ke, e, boundary, tStep );
            answer.assemble( r_loc, c_loc, Ke );
        }
    }
    
    void HardMagneticBoundaryCondition::assembleVector(FloatArray& answer , TimeStep* tStep ,
        CharType type , ValueModeType mode ,
        const UnknownNumberingScheme& s , FloatArray* eNorms, void* lock ) {

        if ( type != ExternalForcesVector ) {
            return;
        }

        FloatArray fe;
        IntArray loc, masterdofids, bNodes;

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
            e->giveBoundaryLocationArray( loc, bNodes, this->dofs, s, &masterdofids );
            this->computeLoadVectorFromElement( fe, e, boundary, tStep, mode );
            answer.assemble( fe, loc );
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
        IntArray bNodes;

        rows.resize( boundaries.giveSize() / 2 );
        cols.resize( boundaries.giveSize() / 2 );

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e   = this->giveDomain()->giveElement( boundaries.at( pos * 2 - 1 ) );
            int boundary = boundaries.at( pos * 2 );

            bNodes = e->giveInterpolation()->boundaryGiveNodes( boundary );

            e->giveBoundaryLocationArray( rows[pos], bNodes, this->dofs, r_s );
            e->giveBoundaryLocationArray( cols[pos], bNodes, this->dofs, c_s );
        }
    }

    void HardMagneticBoundaryCondition::computeLoadVectorFromElement( FloatArray &answer, Element *e, int iSurf, TimeStep *tStep, ValueModeType mode )
    {
        MagneticLoadElementInterface *mfli = static_cast<MagneticLoadElementInterface *>( e->giveInterface( MagneticLoadElementInterfaceType ) );

        if ( !mfli ) {
            OOFEM_ERROR( "Element doesn't implement the required magnetic load interface!" );
        }

        if ( iSurf == -1 ) {
            iSurf = 1;
        }

        double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );


        answer.clear();
	/*
        int order = 4;
        IntegrationRule *iRule = mfli->surfaceGiveIntegrationRule( order, iSurf );
       
        Tensor1_3d SigStarTimesCofN, Normal;
	Tensor2_3d SigStar( load_level*load_level*sigma_star );
        Tensor2_3d F;
        FloatArray vF;

        for ( GaussPoint *gp : *iRule ) {
            // compute normal vector
            FloatArray n, dxdksi, dxdeta;
            mfli->surfaceEvalNormalAt( n, dxdksi, dxdeta, iSurf, gp, tStep );
            // compute surface N matirx
            FloatMatrix N;
            mfli->surfaceEvalNmatrixAt( N, iSurf, gp );

            // compute deformation gradient
            mfli->surfaceEvalDeformationGradientAt( vF, iSurf, gp, tStep );

            //compute the force vector
            F = Tensor2_3d( vF );
            Normal = Tensor2_3d( n );
            auto [J, cofF] = F.compute_determinant_and_cofactor();

            SigStarTimesCofN( i_3 ) = SigStar( i_3, j_3 ) * cofF( j_3, k_3 ) * Normal(k_3);

            answer.plusProduct( N, SigStarTimesCofN.to_voigt_form(), gp->giveWeight() );
        }
	*/
        
    }

    void HardMagneticBoundaryCondition::computeTangentFromElement( FloatMatrix &answer, Element *e, int iSurf, TimeStep *tStep )
    {
       
        MagneticLoadElementInterface *mfli = static_cast<MagneticLoadElementInterface *>( e->giveInterface( MagneticLoadElementInterfaceType ) );

        if ( !mfli ) {
            OOFEM_ERROR( "Element doesn't implement the required magnetic load interface!" );
        }

        if ( iSurf == -1 ) {
            iSurf = 1;
        }

        double load_level = this->giveDomain()->giveFunction( ltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );

        answer.clear();
        int order = 4;
        IntegrationRule *iRule = mfli->surfaceGiveIntegrationRule( order, iSurf );

        IntArray bNodes;
        bNodes = e->giveBoundarySurfaceNodes( iSurf );
        double nNodes = bNodes.giveSize();
        FloatMatrix K;
        FloatMatrix testAnswer( 12, 12 );
	
        if ( e->giveSpatialDimension() == 3 ) {

            Tensor1_3d Normal;
	    //@todo: this is not possible
	    Tensor2_3d SigStar( load_level * load_level * sigma_star );         
            Tensor2_3d F;
            Tensor3_3d SigStarFcrossN;
            FloatArray vF;
	   
            for ( GaussPoint *gp : *iRule ) {
                FloatMatrix dN;

                // compute normal vector
                FloatArray n;
                mfli->surfaceEvalNormalAt( n, iSurf, gp, tStep );

                // compute surface N and B matrix
                FloatMatrix N, B;
                mfli->surfaceEvalNmatrixAt( N, iSurf, gp );
                mfli->surfaceEvalBmatrixAt( B, iSurf, gp );

                // compute deformation gradient
                mfli->surfaceEvalDeformationGradientAt( vF, iSurf, gp, tStep );

                //compute the force vector
                F = Tensor2_3d( FloatArrayF< 9 >(vF) );
		Normal = Tensor1_3d( FloatArrayF< 3 >(n) );
                auto [J, cofF] = F.compute_determinant_and_cofactor();

                SigStarFcrossN( i_3, m_3, n_3 ) = SigStar( i_3, j_3) * F.compute_tensor_cross_product()( j_3, k_3, m_3, n_3) * Normal(k_3);

                N.beTranspositionOf(N);
                FloatMatrix NtimesTensor, NtimesTensortimesB;
                NtimesTensor.beProductOf(N, SigStarFcrossN.to_voigt_form_3x9());
                NtimesTensortimesB.beProductOf(NtimesTensor, B);
                NtimesTensortimesB.times(gp->giveWeight());
                answer.add(NtimesTensortimesB);

            }
        } else if ( e->giveSpatialDimension() == 2 ) {
            OOFEM_ERROR( "Magnetic boundary condition not implemented for 2D domains." );
        }
	    
    }

    void HardMagneticBoundaryCondition::evaluateFreeSpaceStress()
    {
      
        Tensor2_3d SigStar;
	FloatArrayF<3> bb( b_ext );
        Tensor1_3d Bext(bb);
        Tensor2_3d delta(1., 0., 0., 0., 1., 0., 0., 0., 1. );

        SigStar( i_3, j_3 ) = 1 / mu0 * (Bext(i_3) * Bext(j_3) - 0.5 * Bext(p_3) * Bext(p_3) * delta(i_3,j_3));

        sigma_star = SigStar.to_voigt_form();
      
    }

    
}//end namespace oofem;
