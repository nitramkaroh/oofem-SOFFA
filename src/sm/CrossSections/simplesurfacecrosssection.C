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

#include "sm/CrossSections/simplesurfacecrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Elements/structuralelement.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "datastream.h"
#include "contextioerr.h"
#include "engngm.h"
#include "sm/Elements/nlstructuralsurfaceelement.h"

namespace oofem {
REGISTER_CrossSection( SimpleSurfaceCrossSection );

FloatArray
SimpleSurfaceCrossSection::giveRealSurfaceStresses( const FloatArray &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( mode == _3dMat ) {
        return this->giveRealSurfaceStress_3d( strain, normal, gp, tStep );
    //} else if ( mode == _PlaneStrain ) {
    //    return this->giveRealSurfaceStress_PlaneStrain( strain, normal, gp, tStep );
    //} else if ( mode == _PlaneStress ) {
    //    return this->giveRealSurfaceStress_PlaneStress( strain, normal, gp, tStep );
    } else {
        OOFEM_ERROR( "unsupported mode" );
    }
}

FloatArrayF<6>
SimpleSurfaceCrossSection::giveRealSurfaceStress_3d( const FloatArrayF<6> &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
    return mat->giveRealSurfaceStressVector_3d( strain, normal, gp, tStep );
}





FloatArray
SimpleSurfaceCrossSection::giveFirstPKSurfaceStresses( const FloatArray &fullF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        return this->giveFirstPKSurfaceStress_3d( fullF, normal, gp, tStep );
    //} else if ( mode == _PlaneStrain ) {
    //    return this->giveFirstPKSurfaceStress_PlaneStrain( reducedF, normal, gp, tStep );
    //} else if ( mode == _PlaneStress ) {
    //    return this->giveFirstPKSurfaceStress_PlaneStress( reducedF, normal, gp, tStep );
    } else {
        OOFEM_ERROR( "unsupported mode" );
    }
}

FloatArray
SimpleSurfaceCrossSection::giveSecondOrderSurfaceStresses( const FloatArray &fullF, const FloatArray &fullG, const FloatArray &gradI, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        return this->giveSecondOrderSurfaceStress_3d( fullF, fullG, gradI, normal, gp, tStep );
    } else {
        OOFEM_ERROR( "unsupported mode" );
    }
}

FloatArrayF<9>
SimpleSurfaceCrossSection::giveFirstPKSurfaceStress_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
    return mat->giveFirstPKSurfaceStressVector_3d( vF, normal, gp, tStep );
}

FloatArrayF<27>
SimpleSurfaceCrossSection::giveSecondOrderSurfaceStress_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
    return mat->giveSecondOrderSurfaceStressVector_3d( vF, vG, gradI, normal, gp, tStep );
}

FloatMatrixF<9, 9>
SimpleSurfaceCrossSection::giveSurfaceStiffnessMatrix_dPdF_3d( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
    return mat->give3dSurfaceMaterialStiffnessMatrix_dPdF( rMode, normal, gp, tStep );
}

FloatMatrixF<27, 27>
SimpleSurfaceCrossSection::giveSurfaceStiffnessMatrix_dAddF_3d( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
    return mat->give3dSurfaceMaterialStiffnessMatrix_dAddF( rMode, normal, gp, tStep );
}

FloatMatrixF<27, 9>
SimpleSurfaceCrossSection::giveSurfaceStiffnessMatrix_dAdF_3d( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
    return mat->give3dSurfaceMaterialStiffnessMatrix_dAdF( rMode, normal, gp, tStep );
}

FloatMatrixF<6, 6>
SimpleSurfaceCrossSection::giveSurfaceStiffnessMatrix_3d( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
    return mat->give3dSurfaceMaterialStiffnessMatrix( rMode, normal, gp, tStep );
}

//
//FloatArrayF<4>
//SimpleSurfaceCrossSection::giveRealSurfaceStress_PlaneStrain( const FloatArrayF<4> &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
//{
//    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//    return mat->giveRealSurfaceStressVector_PlaneStrain( strain, normal, gp, tStep );
//}

//FloatArrayF<3>
//SimpleSurfaceCrossSection::giveRealSurfaceStress_PlaneStress( const FloatArrayF<3> &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
//{
//    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//    return mat->giveRealSurfaceStressVector_PlaneStress( strain, normal, gp, tStep );
//}
// 
//FloatArrayF<5>
//SimpleSurfaceCrossSection::giveFirstPKSurfaceStress_PlaneStrain( const FloatArrayF<5> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
//{
//    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//    return mat->giveFirstPKSurfaceStressVector_PlaneStrain( vF, normal, gp, tStep );
//}

//FloatArrayF<4>
//SimpleSurfaceCrossSection::giveFirstPKSurfaceStress_PlaneStress( const FloatArrayF<4> &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
//{
//    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//    return mat->giveFirstPKSurfaceStressVector_PlaneStress( strain, normal, gp, tStep );
//}




//void SimpleSurfaceCrossSection::giveCauchyStresses( FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep )
//{
//    // This function returns the Cauchy stress in vector format
//    // corresponding to a given deformation gradient according to the stress-deformation
//    // mode stored in the each gp.
//
//    auto el = dynamic_cast<NLStructuralSurfaceElement *>( gp->giveElement() );
//    FloatArray normal = el->giveNormal( gp );
//
//    MaterialMode mode = gp->giveMaterialMode();
//    auto mat          = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//
//    if ( mode == _3dMat ) {
//        mat->giveCauchyStressVector_3d( answer, normal, gp, reducedvF, tStep );
//    } else if ( mode == _PlaneStrain ) {
//        mat->giveCauchyStressVector_PlaneStrain( answer, normal, gp, reducedvF, tStep );
//    } else if ( mode == _PlaneStress ) {
//        mat->giveCauchyStressVector_PlaneStress( answer, normal, gp, reducedvF, tStep );
//    }
//}


//void SimpleSurfaceCrossSection::giveStiffnessMatrix_dCde( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp,TimeStep *tStep )
//{
//    auto el = dynamic_cast<NLStructuralSurfaceElement *>( gp->giveElement() );
//    FloatArray normal = el->giveNormal( gp );
//
//    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//
//    MaterialMode mode = gp->giveMaterialMode();
//    if ( mode == _3dMat ) {
//        mat->give3dMaterialStiffnessMatrix_dCde( answer, normal, rMode, gp, tStep );
//    } else if ( mode == _PlaneStress ) {
//        mat->givePlaneStressStiffMtrx_dCde( answer, normal, rMode, gp, tStep );
//    } else if ( mode == _PlaneStrain ) {
//        mat->givePlaneStrainStiffMtrx_dCde( answer, normal, rMode, gp, tStep );
//    } else {
//        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString( mode ) );
//    }
//}


//
// 
//void SimpleSurfaceCrossSection::giveCharSurfaceMaterialStiffnessMatrix_dPdF( FloatMatrix &answer, const FloatArray &normal, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
//{
//    MaterialMode mode = gp->giveMaterialMode();
//    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//    if ( mode == _3dMat ) {
//        answer = mat->give3dSurfaceMaterialStiffnessMatrix_dPdF( rMode, normal, gp, tStep );
//    } else if ( mode == _PlaneStress ) {
//        answer = mat->givePlaneStressSurfaceStiffnessMatrix_dPdF( rMode, normal, gp, tStep );
//    } else if ( mode == _PlaneStrain ) {
//        answer = mat->givePlaneStrainSurfaceStiffnessMatrix_dPdF( rMode, normal, gp, tStep );
//    } else {
//        OOFEM_ERROR( "Unsupported material mode for large strain analysis" );
//    }
//}




//FloatMatrixF<5, 5>
//SimpleSurfaceCrossSection::giveSurfaceStiffnessMatrix_dPdF_PlaneStrain( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
//{
//    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//    return mat->givePlaneStrainSurfaceStiffnessMatrix_dPdF( rMode, normal, gp, tStep );
//}

//FloatMatrixF<4, 4>
//SimpleSurfaceCrossSection::giveSurfaceStiffnessMatrix_dPdF_PlaneStress( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
//{
//    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//    return mat->givePlaneStressSurfaceStiffnessMatrix_dPdF( rMode, normal, gp, tStep );
//}





//FloatMatrixF<3, 3>
//SimpleSurfaceCrossSection::giveSurfaceStiffnessMatrix_PlaneStress( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
//{
//    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//    return mat->givePlaneStressSurfaceStiffMtrx( rMode, normal, gp, tStep );
//}
//
//FloatMatrixF<4, 4>
//SimpleSurfaceCrossSection::giveSurfaceStiffnessMatrix_PlaneStrain( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
//{
//    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//    return mat->givePlaneStrainSurfaceStiffMtrx( rMode, normal, gp, tStep );
//}

//void SimpleSurfaceCrossSection::giveCharSurfaceMaterialStiffnessMatrix( FloatMatrix &answer, const FloatArray &normal, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
//{
//    MaterialMode mode = gp->giveMaterialMode();
//
//    auto mat = dynamic_cast<HyperElasticSurfaceMaterial *>( this->giveMaterial( gp ) );
//
//    if ( mode == _3dMat ) {
//        answer = mat->give3dSurfaceMaterialStiffnessMatrix( rMode, normal, gp, tStep );
//    } else if ( mode == _PlaneStress ) {
//        answer = mat->givePlaneStressSurfaceStiffMtrx( rMode, normal, gp, tStep );
//    } else if ( mode == _PlaneStrain ) {
//        answer = mat->givePlaneStrainSurfaceStiffMtrx( rMode, normal, gp, tStep );
//    }
//
//}








} // end namespace oofem
