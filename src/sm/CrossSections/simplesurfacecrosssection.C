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

} // end namespace oofem
