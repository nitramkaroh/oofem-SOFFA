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
#include "sm/Materials/structuralms.h"
#include "basehyperelasticmaterial.h"
#include "tensor/tensor2.h"


///@name Input fields for QuadraticBiotHyperelasticMaterial
//@{
#define _IFT_QuadraticBiotHyperelasticMaterial_Name "QuadraticBiotHyperelasticMat"
#define _IFT_QuadraticBiotHyperelasticMaterial_mu "mu"
#define _IFT_QuadraticBiotHyperelasticMaterial_lambda "lambda"
//@}

namespace oofem {
/**
 * Quadratic Biot hyperelastic material:
 *
 * W = mu * ||U - I||^2 + (lambda/2) * (tr(U - I))^2,   F = R U,  U = sqrt(F^T F)
 *
 * The model is the classical Hooke law between the Biot strain (U - I) and the
 * Biot-type stress t = dW/dU = 2 mu (U - I) + lambda tr(U - I) I; the first
 * Piola stress is P = dW/dF = R t.
 *
 * The class is defined in plane strain only, using the algebraic 2x2 closed
 * forms (no polar decomposition or eigensolver), embedded into the 3d
 * Tensor2_3d / Tensor4_3d representation with F33 = 1:
 *
 *   j = det(F_2d),  s = tr(U_2d) = sqrt(tr(C_2d) + 2 j),
 *   U_2d     = (C_2d + j delta_2d) / s,
 *   U_2d^{-1}= ((tr(C_2d) + j) delta_2d - C_2d) / (s j),
 *   R_2d     = F_2d U_2d^{-1},
 *   t_2d     = 2 mu (U_2d - delta_2d) + lambda (s - 2) delta_2d,
 *   P = R t  (out-of-plane: P33 = lambda (s - 2))
 *
 * The tangent A = dP/dF follows by the product rule
 * A = (dR/dF) t + R (dt/dF) from the closed-form derivatives dU/dF, dR/dF,
 * dt/dF (the 2x2 plane-strain constitutive equations of the model).
 *
 * @author Ondrej Faltus
 */
class QuadraticBiotHyperelasticMaterial : public StructuralMaterial, public BaseHyperElasticMaterial
{
protected:
    // Material parameters (Lame parameters of the in-plane block)
    double mu;
    double lambda;


public:
    QuadraticBiotHyperelasticMaterial( int n, Domain *d );

    void initializeFrom(InputRecord &ir) override;

    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }

    FloatMatrixF<9, 9> give3dMaterialStiffnessMatrix_dPdF( MatResponseMode, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for plane strain only" ); }
    FloatArrayF<9> giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for plane strain only" ); }

    FloatMatrixF<5, 5> givePlaneStrainStiffnessMatrix_dPdF( MatResponseMode, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatArrayF<5> giveFirstPKStressVector_PlaneStrain( const FloatArrayF<5> &vF, GaussPoint *gp, TimeStep *tStep ) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    const char *giveInputRecordName() const override { return _IFT_QuadraticBiotHyperelasticMaterial_Name; }
    const char *giveClassName() const override { return "QuadraticBiotHyperelasticMaterial"; }

private:

    double compute_s_from_F( const Tensor2_3d &F ) const;

    Tensor2_3d compute_U_from_F( const Tensor2_3d &F ) const;
    Tensor2_3d compute_Uinv_from_F( const Tensor2_3d &F ) const;
    Tensor2_3d compute_R_from_F( const Tensor2_3d &F ) const;
    Tensor2_3d compute_t_from_F( const Tensor2_3d &F ) const;
    Tensor2_3d compute_dsdF_from_F( const Tensor2_3d &F ) const;

    Tensor4_3d compute_dUdF_from_F( const Tensor2_3d &F ) const;
    Tensor4_3d compute_dRdF_from_F( const Tensor2_3d &F ) const;
    Tensor4_3d compute_dtdF_from_F( const Tensor2_3d &F ) const;

};
} // end namespace oofem
