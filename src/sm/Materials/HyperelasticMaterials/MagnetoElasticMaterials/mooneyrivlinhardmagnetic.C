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
#include "domain.h"
#include "function.h"



namespace oofem {
REGISTER_Material(ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent);

  ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent::ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent( int n, Domain *d ) : : public MooneyRivlinCompressibleMaterial(n, d)
{
  hmMat  = new HardMagneticMooneyRivlinCompressibleMaterial( n, d );
}


  
void
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent::initializeFrom( InputRecord &ir )
{   
     hmMat->initializeFrom( ir );
}

FloatArrayF<9>
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent::giveFirstPKStressVector_3d_consistent( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
{
    auto [P, B] = hmMat = give_FirstPKStressVector_MagneticInduction_3d( vF, vH, gp, tStep );
    //
    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    //
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );
    //
    return P.to_voigt_form();
}

FloatMatrixF<9, 9>
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent::give3dMaterialStiffnessMatrix_dPdF_consistent( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{

    StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    //
    FloatArrayF<3> vH( status->giveTempHVector() );
    Tensor1_3d H( vH );
    Tensor2_3d F( vF );
    //
    auto [dPdF, dPdH, dBdF, dBdH] = hmMat->giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d(mode, gp, tStep );
    //
    return dPdF.to_voigt_form();
}
  

int
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
  return hmMat->giveIPValue(answer, gp, type, tStep);
 
}

  
} // end namespace oofem
