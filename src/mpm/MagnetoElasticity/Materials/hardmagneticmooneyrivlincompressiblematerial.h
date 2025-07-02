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
#include "tensor/tensor1.h"
#include "tensor/tensor2.h"
#include "tensor/tensor3.h"
#include "tensor/tensor4.h"
#include "tensor/tensor5.h"
#include "tensor/tensor6.h"
#include "tensor/tensor7.h"


///@name Input fields for MooneyRivlinMaterial
//@{
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_Name "hardmagneticmooneyrivlincompressiblemat"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_magnetization "magnetization"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_happ "happ"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_m_ltf "m_ltf"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_hload_ltf "hload_ltf"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_c1 "c1"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_c2 "c2"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_PullBackType "pullbacktype"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_kappaGradFGradF "kappagradfgradf"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterial_kappaGradJGradJ "kappagradjgradj"

//@}

namespace oofem {
/**
 * This class implements Hard Magnetic Compressible Mooney - Rivlin material.
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 *
 */
  enum PullBackType {
        PBT_F       = 0,
        PBT_R       = 1,
        PBT_iFt     = 2,
    };
  
class HardMagneticMooneyRivlinCompressibleMaterial : public MagnetoElasticMaterial, public BaseHyperElasticMaterial
{
protected:
  // Material parameters
  Tensor1_3d M, h_app;
  // Material parameters
  double C1;
  double C2;
  double mu_0 = 1.25663706143e-6;
  int m_ltf = 0, hload_ltf = 0;
  int pb = 1;
  // second gradient material parameter
  double kappaGradFGradF = 0, kappaGradJGradJ = 0;

  PullBackType pullBackType;

public:
  HardMagneticMooneyRivlinCompressibleMaterial( int n, Domain *d );

  void initializeFrom( InputRecord &ir ) override;


  std::tuple<FloatArrayF<9>, FloatArrayF<3> > give_FirstPKStressVector_MagneticInduction_3d( const FloatArrayF<9> &vF, const FloatArrayF<3> &vH, GaussPoint *gp, TimeStep *tStep ) const override;
  std::tuple<FloatArrayF<5>, FloatArrayF<2> > give_FirstPKStressVector_MagneticInduction_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<2> &vH, GaussPoint *gp, TimeStep *tStep ) const override;

  std::tuple<FloatMatrixF< 9, 9 >, FloatMatrixF< 9, 3 >, FloatMatrixF< 3, 9 >, FloatMatrixF< 3, 3 > >  giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
  std::tuple<FloatMatrixF< 5, 5 >, FloatMatrixF< 5, 2 >, FloatMatrixF< 2, 5 >,FloatMatrixF< 2, 2 >> giveConstitutiveMatrices_dPdF_dBdH_dPdH_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
  //
  std::tuple<FloatMatrixF< 9, 9 >, FloatMatrixF< 9, 3 >, FloatMatrixF< 3, 9 >, FloatMatrixF< 3, 3 >> compute_stiff_num(const FloatArray &vF, const FloatArray &vH, GaussPoint *gp, TimeStep *tStep);

  //
    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

  const char *giveInputRecordName() const override { return _IFT_HardMagneticMooneyRivlinCompressibleMaterial_Name; }
  const char *giveClassName() const override { return "HardMagneticMooneyRivlinCompressibleMaterial"; }

  int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
  ////////////////////////////////////////////////////
  FloatArrayF< 8 > give_SecondGradient_FirstPKStressVector_PlaneStrain(const FloatArrayF<8> &vG, GaussPoint *gp, TimeStep *tStep) override;
  FloatArrayF< 27 > give_SecondGradient_FirstPKStressVector_3d(const FloatArrayF<27> &vG, GaussPoint *gp, TimeStep *tStep) override;
  FloatMatrixF<1,1> give_SecondGradient_ConstitutiveMatrix_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
  FloatMatrixF<1,1> give_SecondGradient_ConstitutiveMatrix_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;

  std::tuple<FloatArrayF<5>, FloatArrayF<8> > give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep ) override;
  std::tuple<FloatArrayF<9>, FloatArrayF<27> > give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vGradF, GaussPoint *gp, TimeStep *tStep ) override;
  std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 27>, FloatMatrixF<27, 9>, FloatMatrixF<27, 27> > give_JacobianGradient_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) override;
  std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> > give_JacobianGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) override;


  /************** interface to reduced hard magnetic material ****************************/
  Tensor1_3d giveLagrangianMagnetization(TimeStep *tStep);
  std::tuple<Tensor1_3d,Tensor1_3d> give_EulerianMagnetization_MagneticInduction_FromF(const Tensor2_3d &F, TimeStep *tStep);
  
  PullBackType givePullBackType() const {return pullBackType;}
  /*****************************************************************************************/


  
protected:
  std::tuple<Tensor2_3d, Tensor1_3d> computeFirstPKStressMagneticInductionTensors_3d( const Tensor2_3d &F, const Tensor1_3d &H, const Tensor1_3d &M ) const;
  std::tuple<Tensor4_3d, Tensor3_3d, Tensor3_3d, Tensor2_3d> computeStiffnessTensors_dPdF_dBdH_dPdH_3d( const Tensor2_3d &F, const Tensor1_3d &H, const Tensor1_3d &M ) const;

private:
  Tensor1_3d compute_gradJ_3d( const Tensor2_3d &F, const Tensor3_3d &G ) const;
  std::tuple<Tensor3_3d, Tensor4_3d> compute_gradJ_derivatives_3d( const Tensor2_3d &F, const Tensor3_3d &G ) const;
  std::tuple<Tensor5_3d, Tensor6_3d, Tensor6_3d, Tensor7_3d> compute_gradJ_secondDerivatives_3d( const Tensor2_3d &F, const Tensor3_3d &G ) const;
};





  
} // end namespace oofem
#endif
