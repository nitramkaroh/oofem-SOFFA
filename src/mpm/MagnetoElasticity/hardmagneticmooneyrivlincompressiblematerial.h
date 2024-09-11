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

#ifndef hardmagneticmooneyrivlinmaterial_h
#define hardmagneticmooneyrivlinmaterial_h

#include "sm/Materials/Hyperelastic/mooneyrivlincompressiblematerial.h"


///@name Input fields for MooneyRivlinMaterial
//@{
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_Name "hardmagneticmooneyrivlincompressiblemat"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_magnetization "magnetization"
//@}

namespace oofem {
/**
 * This class implements Hard Magnetic Compressible Mooney - Rivlin material.
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 *
 */
class MooneyRivlinCompressibleMaterial : public StructuralMaterial, public BaseHyperElasticMaterial
{
protected:
    // Material parameters
    Tensor1_3d magnetization;


public:
    HardMagneticMooneyRivlinCompressibleMaterial(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;


  std::tuple<FloatArrayF< 9 >,FloatArrayF< 3 >> give_1PKStressVector_MagneticInduction_3d(FloatMatrixF<9> vF, FloatMatrixF<3> vH, GuassPoint *gp, TimeStep *tStep) const override;
  std::tuple<FloatArrayF< 5 >,FloatArrayF< 2 >> give_1PKStressVector_MagneticInduction_PlaneStrain(FloatMatrixF<5> vF, FloatMatrixF<2> vH, GuassPoint *gp, TimeStep *tStep) const override;

  std::tuple<FloatMatrixF< 9, 9 >, FloatMatrixF< 3, 3 >, FloatMatrixF< 9, 3 >>  giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
   std::tuple<FloatMatrixF< 5, 5 >, FloatMatrixF< 5, 2 >, FloatMatrixF< 5, 2 >> giveConstitutiveMatrices_dPdF_dBdH_dPdH_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;


    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    const char *giveInputRecordName() const override { return _IFT_HardMagneticMooneyRivlinCompressibleMaterial_Name; }
    const char *giveClassName() const override { return "HardMagneticMooneyRivlinCompressibleMaterial"; }
};
} // end namespace oofem
#endif
