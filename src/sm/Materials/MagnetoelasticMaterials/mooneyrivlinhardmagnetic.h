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

#pragma once

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/HyperelasticMaterials/mooneyrivlincompressiblematerial.h"
#include "sm/Materials/structuralms.h"
#include "floatarray.h"


///@name Input fields for MooneyRivlinHardMagnetic
//@{
#define _IFT_MooneyRivlinHardMagnetic_Name "mooneyrivlinhardmagnetic"
#define _IFT_MooneyRivlinHardMagnetic_B_app "b_app" //applied magnetic field
#define _IFT_MooneyRivlinHardMagnetic_B_res "b_res" //residual magnetic field in material
#define _IFT_MooneyRivlinHardMagnetic_H "h" //magnetic field in material
#define _IFT_MooneyRivlinHardMagnetic_mu_0 "mu_0" //magnetic permeability of free space (optional, default value is 4*pi*10^7)
#define _IFT_MooneyRivlinHardMagnetic_ltf "ltf" //load time function for applied field
#define _IFT_MooneyRivlinHardMagnetic_mode "mode" //mode of magnetic energy 1 = consistent (default), 2 = MIT, 3 - fix, 4 - with Ogden pullback, 5 - multiphysics
#define _IFT_MooneyRivlinHardMagnetic_referenceB "referenceb" //whether b_app is given in reference configuration, default false
//@}

// base vacuum permeability constant definition for OOFEM, in basic SI units
#define BASE_VACUUM_PERMEABILITY_MU_0 1.25663706143e-6

namespace oofem {
/**
 * This class implements a hard magnetic material based on Mooney-Rivlin.
 *
 * @author Ondrej Faltus
 *
 * @note
 */
class MooneyRivlinHardMagnetic : public MooneyRivlinCompressibleMaterial//, public MagnetoelasticMaterial
{
protected:
    // Material parameters
    FloatArrayF<3> B_app;
    FloatArrayF<3> B_res;
    FloatArrayF<3> H;
    double mu_0;
    int ltf_index, materialMode;
    bool referenceB;

public:
    MooneyRivlinHardMagnetic( int n, Domain *d ) :
        MooneyRivlinCompressibleMaterial( n, d ){};

    void initializeFrom(InputRecord &ir) override;

    //hyperelasticity
    FloatMatrixF< 9, 9 >give3dMaterialStiffnessMatrix_dPdF(MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const override;

    ////magnetoelastic multiphysics
    //FloatArrayF<3> giveMagneticInduction_3d( const FloatArrayF<3> &vH, GaussPoint *gp, TimeStep *tStep ) const override;
    //FloatMatrixF<3, 9> give3dMagnetoelasticCouplingTensor_dBdF(MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    //FloatMatrixF<3, 3> give3dPermeabilityMatrix_dBdH( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const override;

    const char *giveInputRecordName() const override { return _IFT_MooneyRivlinHardMagnetic_Name; }
    const char *giveClassName() const override { return "MooneyRivlinHardMagnetic"; }

private:
      FloatArrayF< 9 > computeFirstPKStressVector_3d_consistent(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const;
      FloatMatrixF< 9, 9 > compute3dMaterialStiffnessMatrix_dPdF_consistent(MatResponseMode matResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const;

      FloatArrayF< 9 > computeFirstPKStressVector_3d_mit(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const;
      FloatMatrixF< 9, 9 > compute3dMaterialStiffnessMatrix_dPdF_mit(MatResponseMode matResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const;

      FloatArrayF< 9 > computeFirstPKStressVector_3d_consistentfix(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const;
      FloatMatrixF< 9, 9 > compute3dMaterialStiffnessMatrix_dPdF_consistentfix(MatResponseMode matResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const;

      FloatArrayF< 9 > computeFirstPKStressVector_3d_ogdenpullback(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const;
      FloatMatrixF< 9, 9 > compute3dMaterialStiffnessMatrix_dPdF_ogdenpullback(MatResponseMode matResponseMode,
                                                           GaussPoint *gp, TimeStep *tStep ) const;

      FloatArrayF< 9 > computeFirstPKStressVector_3d_universal(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const;
      FloatMatrixF< 9, 9 > compute3dMaterialStiffnessMatrix_dPdF_universal(MatResponseMode matResponseMode,
                                                           GaussPoint *gp, TimeStep *tStep ) const;

      FloatArrayF< 9 > computeFirstPKStressVector_3d_universal_nobb(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const;
      FloatMatrixF< 9, 9 > compute3dMaterialStiffnessMatrix_dPdF_universal_nobb(MatResponseMode matResponseMode,
                                                           GaussPoint *gp, TimeStep *tStep ) const;

      //FloatArrayF< 9 > computeFirstPKStressVector_3d_multiphysics(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const;
      //FloatMatrixF<9, 9> compute3dMaterialStiffnessMatrix_dPdF_multiphysics( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const;
      //FloatArrayF<3> computeMagneticInduction_3d_multiphysics( const FloatArrayF<3> &vH, GaussPoint *gp, TimeStep *tStep ) const;
      //FloatMatrixF<3, 9> compute3dMagnetoelasticCouplingTensor_dBdF_multiphysics( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const;
      //FloatMatrixF<3, 3> compute3dPermeabilityMatrix_dBdH_multiphysics( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const;

      //todo remove this \/
      //int giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep ) override;
};
} // end namespace oofem
