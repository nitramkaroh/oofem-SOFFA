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

#include "mooneyrivlinhardmagnetic.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"



namespace oofem {
REGISTER_Material(MooneyRivlinHardMagnetic);

FloatArrayF< 9 >
MooneyRivlinHardMagnetic::giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{

    //get stress from Mooney-Rivlin hyperelasticity
    FloatArrayF<9> mr_pkstress = MooneyRivlinCompressibleMaterial::giveFirstPKStressVector_3d( vF, gp, tStep );

    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    return mr_pkstress;
}



FloatMatrixF< 9, 9 >
MooneyRivlinHardMagnetic::give3dMaterialStiffnessMatrix_dPdF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    // get stiffness from Mooney-Rivlin hyperelasticity
    FloatMatrixF<9, 9> mr_stiffness = MooneyRivlinCompressibleMaterial::give3dMaterialStiffnessMatrix_dPdF( mode, gp, tStep );

    return mr_stiffness;
}


void MooneyRivlinHardMagnetic::initializeFrom( InputRecord &ir )
{   
    MooneyRivlinCompressibleMaterial::initializeFrom( ir );

    IR_GIVE_FIELD( ir, B_app, _IFT_MooneyRivlinHardMagnetic_B_app );
    IR_GIVE_FIELD( ir, B_res, _IFT_MooneyRivlinHardMagnetic_B_res );
}
} // end namespace oofem
