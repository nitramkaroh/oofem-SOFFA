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


///@name Input fields for HarmonicHyperelasticMaterial
//@{
#define _IFT_HarmonicHyperelasticMaterial_Name "harmonichyperelasticmat"
#define _IFT_HarmonicHyperelasticMaterial_alpha "alpha"
#define _IFT_HarmonicHyperelasticMaterial_beta "beta"
#define _IFT_HarmonicHyperelasticMaterial_mu "mu"
//@}

namespace oofem {
/**
 * This class implements a harmonic hyperelastic material
 *
 * @author Ondrej Faltus
 *
 * @note References: C.Q. Ru: On complex-variable formulation for finite plane elastostatics of harmonic materials, Acta Mechanica, 2002
 *
 * Free energy is considered as:
 *
 * W = 2 mu (f(I) - J)
 * 
 * with df/dI = 1/(4 alpha) * (I + sqrt(I^2 - 16*alpha*beta))
 * 
 * where mu, alpha, beta are material parameters, and
 * I1 is a trace of a 2D stretch tensor, i.e. I = lambda_1 + lambda_2
 * 
 * The material is thus defined in plane strain only.
 * 
 */
class HarmonicHyperelasticMaterial : public StructuralMaterial, public BaseHyperElasticMaterial
{
protected:
    // Material parameters
    double alpha;
    double beta;
    double mu;


public:
    HarmonicHyperelasticMaterial( int n, Domain *d );

    void initializeFrom(InputRecord &ir) override;

    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    
    FloatMatrixF<9, 9> give3dMaterialStiffnessMatrix_dPdF( MatResponseMode, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for plane strain only" ); }
    FloatArrayF<9> giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for plane strain only" ); }

    FloatMatrixF<5, 5> givePlaneStrainStiffnessMatrix_dPdF( MatResponseMode, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatArrayF<5> giveFirstPKStressVector_PlaneStrain( const FloatArrayF<5> &vF, GaussPoint *gp, TimeStep *tStep ) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    const char *giveInputRecordName() const override { return _IFT_HarmonicHyperelasticMaterial_Name; }
    const char *giveClassName() const override { return "HarmonicHyperelasticMaterial"; }

private:

    double compute_dFunctiondI_from_F(const Tensor2_3d &F) const;
    double compute_d2FunctiondIdI_from_F( const Tensor2_3d &F ) const;

    double compute_I_of_U_2d_from_F( const Tensor2_3d &F) const;
    Tensor2_3d compute_dI_of_U_2d_dF_from_F( const Tensor2_3d &F ) const;
    Tensor4_3d compute_d2I_of_U_2d_dFdF_from_F( const Tensor2_3d &F ) const;

};
} // end namespace oofem

