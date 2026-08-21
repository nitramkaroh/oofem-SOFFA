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


#ifndef simopistermaterial_h
#define simopistermaterial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "basehyperelasticmaterial.h"
#include "tensor/tensor1.h"
#include "tensor/tensor3.h"
#include "tensor/tensor5.h"
#include "tensor/tensor6.h"


///@name Input fields for SimoPisterMaterial
//@{
#define _IFT_SimoPisterMaterial_Name "simopistermat"
#define _IFT_SimoPisterMaterial_g "g"
#define _IFT_SimoPisterMaterial_bulkLTF "bulk_ltf"
//@}

// For the gradient version of the material
#define _IFT_GradientSimoPisterMaterial_Name "gradientsimopistermat"
#define _IFT_GradientSimoPisterMaterial_omega "omega"

namespace oofem {
/**
 * Free energy is considered as:
 * \f$[
 * \rho_0 \psi = U(J) + G( 0.5 * I_1 - ln J) ]\f$
 *  This form of energy corresponds to a neo-Hookean material which is extended to the compressible * range by adding an extra function depending on J.
 * @author Martin Horak, nitramkaroh@seznam.cz
 * @note Reference: article{simo1984remarks,
 * title={Remarks on rate constitutive equations for finite deformation problems: computational implications},
 * author={Simo, Juan C and Pister, Karl S},
 * journal={Computer Methods in Applied Mechanics and Engineering},
 * volume={46},
 * number={2},
 * pages={201--215},
 * year={1984},
 * publisher={Elsevier}
 * }
 *
 *
 */
class SimoPisterMaterial : public StructuralMaterial, public BaseHyperElasticMaterial
{
protected:
    double G = 0.; ///< Shear modulus.
    int bulk_ltf = 0;

public:
    SimoPisterMaterial(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }


    FloatMatrixF< 9, 9 >give3dMaterialStiffnessMatrix_dPdF(MatResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const override;

    FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    const char *giveInputRecordName() const override { return _IFT_SimoPisterMaterial_Name; }
    const char *giveClassName() const override { return "SimoPisterMaterial"; }
};


// Second gradient formulation of Simo-Pister material, see article simo1984remarks for details.
class GradientSimoPisterMaterial : public SimoPisterMaterial
{
protected:
    // Gradient length-scale / stiffness parameter
    double omega;

public:
    GradientSimoPisterMaterial( int n, Domain *d );

    void initializeFrom( InputRecord &ir ) override;

    const char *giveInputRecordName() const override { return _IFT_GradientSimoPisterMaterial_Name; }
    const char *giveClassName() const override { return "GradientSimoPisterMaterial"; }

    //// Standard PK1 and its derivative (Inherited from base, but overridden to sync status if needed)
    //FloatArrayF<9> giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const override;
    //FloatMatrixF<9, 9> give3dMaterialStiffnessMatrix_dPdF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const override;

    virtual FloatArrayF<27> giveSecondOrderStressVector_3d( const FloatArrayF<9> &reducedF, const FloatArrayF<27> &reducedG, GaussPoint *gp, TimeStep *tStep ) const override;
    virtual FloatMatrixF<27, 27> give3dMaterialStiffnessMatrix_dAddF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const override;
    virtual FloatMatrixF<27, 9> give3dMaterialStiffnessMatrix_dAdF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const override;

    MaterialStatus *CreateStatus( GaussPoint *gp ) const override;

    bool supportsSecondGradient() const override { return true; }


};


} // end namespace oofem
#endif
