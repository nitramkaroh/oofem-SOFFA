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

#ifndef reducedmooneyrivlinmaterial_consistent_h
#define reducedmooneyrivlinmaterial_consistent_h

#include "sm/Materials/HyperelasticMaterials/mooneyrivlincompressiblematerial.h"
#include "mpm/MagnetoElasticity/Materials/hardmagneticmooneyrivlincompressiblematerial.h"
#include "mpm/MagnetoElasticity/Materials/magnetoelasticmaterialstatus.h"

///@name Input fields for MooneyRivlinMaterial
//@{
#define _IFT_ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent_Name "reducedhardmagneticmooneyrivlincompressiblemat"
#define _IFT_ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent_Demagnetization "demag"
#define _IFT_ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent_AppliedInduction "bapp"
#define _IFT_ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent_bltf "bltf" //load time function for applied field
#define _IFT_ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent_hmMatNumber "hmmatnumber"

//@}

namespace oofem {
  //Class MagnetoElasticMaterialStatus;
/**
 * This class implements reduced hard magnetic Compressible Mooney - Rivlin material.
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 *
 * @note our super cool paper that we have to write
 * Free energy is considered as: I think this model has no energy
 * Compressible Neo-Hookean model is obtained by setting @f$C_2 = 0@f$

  
 */
class ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent : public MooneyRivlinCompressibleMaterial
{
protected:
    // Material parameters
    double mu_0 = 1.25663706143e-6;
    FloatArray demagFactor;
    HardMagneticMooneyRivlinCompressibleMaterial *hmMat;
    Tensor1_3d b_app;
    int bltf_index;
    int hmMatNumber;
    
public:
    ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent (int n, Domain *d);
  ~ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent (){delete hmMat;}

    void initializeFrom(InputRecord &ir) override;

  std::tuple<Tensor1_3d,Tensor1_3d> give_EulerianMagnetization_MagneticInduction_FromF(const Tensor2_3d &F, TimeStep *tStep);
  // 
    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    FloatMatrixF< 9, 9 >give3dMaterialStiffnessMatrix_dPdF(MatResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const override;

    FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const override;

  MaterialStatus *CreateStatus(GaussPoint *gp) const override;
  Tensor1_3d giveLagrangianMagnetization(TimeStep *tStep) const;
  
  int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
  void postInitialize() override;

    const char *giveInputRecordName() const override { return _IFT_ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent_Name; }
    const char *giveClassName() const override { return "ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent"; }
};





class ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent : public StructuralMaterialStatus
{
protected:
    int slaveMat;
    /// Gauss Point for slave material status
    GaussPoint *slaveGp;
    /// Slave material status
  MagnetoElasticMaterialStatus *slaveStatus;
  //
  Domain *domain;
    
public:
    ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent(GaussPoint *g, Domain *d, int s);
    ~ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent();

    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void initTempStatus() override;
    void updateYourself(TimeStep *) override;

    GaussPoint* giveSlaveGaussPoint()
    { return this->slaveGp;}
    MaterialStatus *giveSlaveStatus() const
    {return this->slaveStatus;}
    void createSlaveStatus();

  
    const char *giveClassName() const override { return "ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent"; }


};

  
} // end namespace oofem
#endif
