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

#include "linearelasticsurfacematerial.h"
#include "gausspoint.h"
#include "sm/CrossSections/simplecrosssection.h"
#include "sm/Materials/structuralms.h"
#include "dynamicinputrecord.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

namespace oofem {
void LinearElasticSurfaceMaterial ::initializeFrom( InputRecord &ir )
{
    StructuralMaterial ::initializeFrom( ir );

    preCastStiffnessReduction = 0.99999999;
    IR_GIVE_OPTIONAL_FIELD( ir, preCastStiffnessReduction, _IFT_LinearElasticSurfaceMaterial_preCastStiffRed );
}


void LinearElasticSurfaceMaterial ::giveInputRecord( DynamicInputRecord &input )
{
    StructuralMaterial ::giveInputRecord( input );
    input.setField( this->preCastStiffnessReduction, _IFT_LinearElasticSurfaceMaterial_preCastStiffRed );
}

void LinearElasticSurfaceMaterial ::computesSubTangents()
{
    tangentPlaneStrain = {
        tangent( 0, 0 ),
        tangent( 1, 0 ),
        tangent( 2, 0 ),
        tangent( 3, 0 ),
        tangent( 0, 1 ),
        tangent( 1, 1 ),
        tangent( 2, 1 ),
        tangent( 3, 1 ),
        tangent( 0, 2 ),
        tangent( 1, 2 ),
        tangent( 2, 2 ),
        tangent( 3, 2 ),
        tangent( 0, 3 ),
        tangent( 1, 3 ),
        tangent( 2, 3 ),
        tangent( 3, 3 ),
    };

    auto c                     = inv( tangent );
    FloatMatrixF<3, 3> reduced = {
        c( 0, 0 ),
        c( 0, 1 ),
        c( 0, 5 ),
        c( 1, 0 ),
        c( 1, 1 ),
        c( 1, 5 ),
        c( 5, 0 ),
        c( 5, 1 ),
        c( 5, 5 ),
    };
    tangentPlaneStress = inv( reduced );
}


FloatMatrixF<6, 6>
LinearElasticSurfaceMaterial ::give3dSurfaceMaterialStiffnessMatrix( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    if ( tStep->giveIntrinsicTime() < this->castingTime ) {
        return tangent * ( 1. - this->preCastStiffnessReduction );
    } else {
        return tangent;
    }
}


FloatArrayF<6>
LinearElasticSurfaceMaterial ::giveThermalDilatationVector( GaussPoint *gp, TimeStep *tStep ) const
{
    return alpha;
}


FloatArrayF<6>
LinearElasticSurfaceMaterial ::giveRealSurfaceStressVector_3d( const FloatArrayF<6> &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

    auto d = this->give3dSurfaceMaterialStiffnessMatrix( TangentStiffness, normal, gp, tStep );

    FloatArrayF<6> stress;
    if ( this->castingTime < 0. ) { // no changes in material stiffness ->> total formulation
        auto thermalStrain = this->computeStressIndependentStrainVector_3d( gp, tStep, VM_Total );
        auto strainVector  = strain - thermalStrain;
        stress             = dot( d, strainVector );
    } else { // changes in material stiffness ->> incremental formulation
        auto thermalStrain   = this->computeStressIndependentStrainVector_3d( gp, tStep, VM_Incremental );
        auto strainIncrement = strain - thermalStrain - FloatArrayF<6>( status->giveStrainVector() );

        stress = dot( d, strainIncrement ) + status->giveStressVector();
    }

    // update gp
    status->letTempStrainVectorBe( strain );
    status->letTempStressVectorBe( stress );
    return stress;
}


double
LinearElasticSurfaceMaterial ::giveEnergyDensity( GaussPoint *gp, TimeStep *tStep )
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    const FloatArray &strain         = status->giveTempStrainVector();
    const FloatArray &stress         = status->giveTempStressVector();

    return 0.5 * stress.dotProduct( strain );
}

int LinearElasticSurfaceMaterial ::giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep )
{
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    if ( type == IST_ElasticStrainTensor ) {
        this->giveStressDependentPartOfStrainVector( answer, gp, status->giveStrainVector(), tStep, VM_Total );
        return 1;
    } else if ( type == IST_ThermalStrainTensor ) {
        answer = this->computeStressIndependentStrainVector_3d( gp, tStep, VM_Total );
        return 1;
    } else if ( type == IST_CreepStrainTensor ) {
        answer.resize( 6 );
        answer.zero();
        return 1;
    } else if ( type == IST_FreeEnergyDensity ) {
        answer.resize( 1 );
        answer.at( 1 ) = this->giveEnergyDensity( gp, tStep );
        return 1;
    } else {
        return StructuralMaterial ::giveIPValue( answer, gp, type, tStep );
    }
}


MaterialStatus *
LinearElasticSurfaceMaterial ::CreateStatus( GaussPoint *gp ) const
{
    return new StructuralMaterialStatus( gp );
}
} // end namespace oofem
