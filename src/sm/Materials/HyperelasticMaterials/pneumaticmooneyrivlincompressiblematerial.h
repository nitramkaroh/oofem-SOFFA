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

#include "sm/Materials/HyperelasticMaterials/mooneyrivlincompressiblematerial.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "basehyperelasticmaterial.h"


///@name Input fields for PneumaticMooneyRivlinMaterial
//@{
#define _IFT_PneumaticMooneyRivlinCompressibleMaterial_Name "pneumaticmooneyrivlincompressiblemat"
#define _IFT_PneumaticMooneyRivlinCompressibleMaterial_pressure "dp"
#define _IFT_PneumaticMooneyRivlinCompressibleMaterial_pressure_ltf "pressure_ltf"
//@}

namespace oofem {
/**
 * This class implements a pneumatic Compressible Mooney - Rivlin material.
 *
 * @author Ondrej Faltus
 *
 * @note References: Faltus, Horak, Doskar, Rokos: Third Medium Finite Element Contact Formulation for Pneumatically Actuated Systems
 *
 * An extension of the Mooney - Rivlin material to introduced a prescribed pneumatic load given in the input file
 */
class PneumaticMooneyRivlinCompressibleMaterial : public MooneyRivlinCompressibleMaterial
{
protected:
    // Material parameters
    double pressure = 0.0;
    int pressure_ltf;


public:
    PneumaticMooneyRivlinCompressibleMaterial(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    FloatMatrixF< 9, 9 >give3dMaterialStiffnessMatrix_dPdF(MatResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const override;

    FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const override;

    const char *giveInputRecordName() const override { return _IFT_PneumaticMooneyRivlinCompressibleMaterial_Name; }
    const char *giveClassName() const override { return "PneumaticMooneyRivlinCompressibleMaterial"; }
};
} // end namespace oofem
