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

#pragma once

#include "material.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "mathfem.h"
#include "matconst.h"
#include "matstatus.h"
#include "valuemodetype.h"
#include <vector>
#include "tensor/tensor1.h"
#include "tensor/tensor2.h"
#include "tensor/tensor3.h"
#include "tensor/tensor4.h"
#include "tensor/tensor5.h"
#include "tensor/tensor6.h"
#include "tensor/tensor7.h"
#include "thirdmediummaterialstatus.h"

///@name Input fields for Third Medium Material
//@{
//@}
#define _IFT_ThirdMediumMaterial_kappaGradFGradF "kappagradfgradf"
#define _IFT_ThirdMediumMaterial_kappaGradJGradJ "kappagradjgradj"
#define _IFT_ThirdMediumMaterial_kappaGradRGradR "kappagradrgradr"
#define _IFT_ThirdMediumMaterial_kappaJbar "kappajbar"
#define _IFT_ThirdMediumMaterial_kappaFbar "kappafbar"
#define _IFT_ThirdMediumMaterial_kappaFskewBar "kappafskewbar"
#define _IFT_ThirdMediumMaterial_eLinear "e_lin"
#define _IFT_ThirdMediumMaterial_nuLinear "nu_lin"

namespace oofem {

class GaussPoint;
///@todo Update documentation
class OOFEM_EXPORT ThirdMediumMaterial
{
protected:
  // material parameters
  double kappaGradFGradF = 0, kappaGradJGradJ = 0, kappaGradRGradR = 0, kappaJbar = 0, kappaFbar = 0, kappaFskewBar = 0;

  //elastic tangents
  double E = 0., nu = 0.;
  FloatMatrixF<6, 6> linearElasticTangent;
  FloatMatrixF<4, 4> linearElasticTangentPlaneStrain;

  public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    ThirdMediumMaterial( )  { ; }

    const char *giveClassName() const { return "ThirdMediumMaterial"; }

    void initializeFrom(InputRecord &ir);    

    virtual ThirdMediumMaterialStatus* giveStatus(GaussPoint* gp);

    ////////////////////////SECOND_GRADIENT////////////////////////////
    virtual FloatArrayF<8> give_SecondGradient_FirstPKStressVector_PlaneStrain( const FloatArrayF<8> &vG, GaussPoint *gp, TimeStep *tStep );
    virtual FloatArrayF<27> give_SecondGradient_FirstPKStressVector_3d( const FloatArrayF<27> &vG, GaussPoint *gp, TimeStep *tStep );
    virtual FloatMatrixF<1, 1> give_SecondGradient_ConstitutiveMatrix_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep );
    virtual FloatMatrixF<1, 1> give_SecondGradient_ConstitutiveMatrix_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep );

    ////////////////////////FIRST AND SECOND GRADIENT////////////////////////////
    virtual std::tuple<FloatArrayF<5>, FloatArrayF<8> > give_FirstSecondGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep );
    virtual std::tuple<FloatArrayF<9>, FloatArrayF<27> > give_FirstSecondGradient_FirstPKStressVector_SecondOrderStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vGradF, GaussPoint *gp, TimeStep *tStep );
    virtual std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 27>, FloatMatrixF<27, 9>, FloatMatrixF<27, 27> > give_FirstSecondGradient_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep );
    virtual std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> > give_FirstSecondGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep );

    ////////////////////////FBAR////////////////////////////
    virtual std::tuple<FloatArrayF<5>, FloatArrayF<5> > give_Fbar_FirstPKStressVector_FPKbarStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<5> &vFbar, GaussPoint *gp, TimeStep *tStep );
    virtual std::tuple<FloatArrayF<9>, FloatArrayF<9> > give_Fbar_FirstPKStressVector_FPKbarStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<9> &vFbar, GaussPoint *gp, TimeStep *tStep );
    virtual std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9> > give_Fbar_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep );
    virtual std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 5>, FloatMatrixF<5, 5>, FloatMatrixF<5, 5> > give_Fbar_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep );

    ////////////////////////LINEAR ELASTICITY////////////////////////////
    virtual FloatArrayF<4> give_LinearElasticity_RealStressVector_PlaneStrain( const FloatArrayF<4> &eps, GaussPoint *gp, TimeStep *tStep );
    virtual FloatArrayF<6> give_LinearElasticity_RealStressVector_3d( const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep );
    virtual FloatMatrixF<6, 6> give_LinearElasticity_ConstitutiveMatrix_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep );
    virtual FloatMatrixF<4, 4> give_LinearElasticity_ConstitutiveMatrix_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep );

 private:
     
    // Jacobian
    std::tuple<FloatArrayF<5>, FloatArrayF<8> > give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep );
    std::tuple<FloatArrayF<9>, FloatArrayF<27> > give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vGradF, GaussPoint *gp, TimeStep *tStep );
    std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 27>, FloatMatrixF<27, 9>, FloatMatrixF<27, 27> > give_JacobianGradient_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep );
    std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> > give_JacobianGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep );

    // for Jacobian gradients
    Tensor1_3d compute_gradJ_3d( const Tensor2_3d &F, const Tensor3_3d &G ) const;
    std::tuple<Tensor3_3d, Tensor4_3d> compute_gradJ_derivatives_3d( const Tensor2_3d &F, const Tensor3_3d &G ) const;
    std::tuple<Tensor5_3d, Tensor6_3d, Tensor6_3d, Tensor7_3d> compute_gradJ_secondDerivatives_3d( const Tensor2_3d &F, const Tensor3_3d &G ) const;

    // Rotation by Wriggers
    std::tuple<FloatArrayF<5>, FloatArrayF<8> > give_RotationGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep );
    std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> > give_RotationGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep );

    // for rotation gradients
    double compute_Z_fraction_PlaneStrain( const Tensor2_3d &F ) const;
    Tensor2_3d compute_Z_fraction_derivative_PlaneStrain( const Tensor2_3d &F ) const;
    Tensor4_3d compute_Z_fraction_secondDerivative_PlaneStrain( const Tensor2_3d &F ) const;
    Tensor6_3d compute_Z_fraction_thirdDerivative_PlaneStrain( const Tensor2_3d &F ) const;
    Tensor1_3d compute_gradPhi_PlaneStrain( const Tensor2_3d &F, const Tensor3_3d &gradF ) const;
    std::tuple<Tensor3_3d, Tensor4_3d> compute_gradPhi_derivatives_PlaneStrain( const Tensor2_3d &F, const Tensor3_3d &gradF ) const;
    std::tuple<Tensor5_3d, Tensor6_3d, Tensor7_3d> compute_gradPhi_secondDerivatives_PlaneStrain( const Tensor2_3d &F, const Tensor3_3d &gradF ) const;

    //for linear elasticity
    void initTangents();

};
} // end namespace oofem
