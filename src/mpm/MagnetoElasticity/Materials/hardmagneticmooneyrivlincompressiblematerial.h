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

#include "sm/Materials/HyperelasticMaterials/mooneyrivlincompressiblematerial.h"
#include "sm/Materials/HyperelasticMaterials/basehyperelasticmaterial.h"
#include "magnetoelasticmaterial.h"


///@name Input fields for MooneyRivlinMaterial
//@{
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_Name "hardmagneticmooneyrivlincompressiblemat"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_magnetization "magnetization"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_c1 "c1"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_c2 "c2"

//@}

namespace oofem {
/**
 * This class implements Hard Magnetic Compressible Mooney - Rivlin material.
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 *
 */
class HardMagneticMooneyRivlinCompressibleMaterial : public MagnetoElasticMaterial, public BaseHyperElasticMaterial
{
protected:
    // Material parameters
    Tensor1_3d M;
  // Material parameters
    double C1;
    double C2;
    double mu_0 = 1.25663706143e-6;



public:
    HardMagneticMooneyRivlinCompressibleMaterial(int n, Domain *d);

  void initializeFrom(InputRecord &ir) override;


  std::tuple<FloatArrayF< 9 >,FloatArrayF< 3 >> give_FirstPKStressVector_MagneticInduction_3d(const FloatArrayF<9> &vF, const FloatArrayF<3> &vH, GaussPoint *gp, TimeStep *tStep) const override;
   std::tuple<FloatArrayF< 5 >,FloatArrayF< 2 >> give_FirstPKStressVector_MagneticInduction_PlaneStrain(const FloatArrayF<5> &vF, const FloatArrayF<2> &vH, GaussPoint *gp, TimeStep *tStep) const override;

  std::tuple<FloatMatrixF< 9, 9 >, FloatMatrixF< 9, 3 >, FloatMatrixF< 3, 9 >, FloatMatrixF< 3, 3 > >  giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
  std::tuple<FloatMatrixF< 5, 5 >, FloatMatrixF< 5, 2 >, FloatMatrixF< 2, 5 >,FloatMatrixF< 2, 2 >> giveConstitutiveMatrices_dPdF_dBdH_dPdH_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
  //
  std::tuple<FloatMatrixF< 9, 9 >, FloatMatrixF< 9, 3 >, FloatMatrixF< 3, 9 >, FloatMatrixF< 3, 3 >> compute_stiff_num(const FloatArray &vF, const FloatArray &vH, GaussPoint *gp, TimeStep *tStep);
  //
    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    const char *giveInputRecordName() const override { return _IFT_HardMagneticMooneyRivlinCompressibleMaterial_Name; }
    const char *giveClassName() const override { return "HardMagneticMooneyRivlinCompressibleMaterial"; }
};
} // end namespace oofem
#endif
