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
#define _IFT_MooneyRivlinHardMagnetic_Name "mooneyrivlincompressiblemat"
#define _IFT_MooneyRivlinHardMagnetic_B_app "B_app"
#define _IFT_MooneyRivlinHardMagnetic_B_res "B_res"
#define _IFT_MooneyRivlinHardMagnetic_mu_0 "mu_0"
//@}

namespace oofem {
/**
 * This class implements a hard magnetic material based on Mooney-Rivlin.
 *
 * @author Ondrej Faltus
 *
 * @note
 */
class MooneyRivlinHardMagnetic : public MooneyRivlinCompressibleMaterial
{
protected:
    // Material parameters
    FloatArrayF<3> B_app;
    FloatArrayF<3> B_res;
    double mu_0;

public:
    MooneyRivlinHardMagnetic( int n, Domain *d ) :
        MooneyRivlinCompressibleMaterial( n, d ){};

    void initializeFrom(InputRecord &ir) override;

    FloatMatrixF< 9, 9 >give3dMaterialStiffnessMatrix_dPdF(MatResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const override;

    FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const override;

    const char *giveInputRecordName() const override { return _IFT_MooneyRivlinHardMagnetic_Name; }
    const char *giveClassName() const override { return "MooneyRivlinHardMagnetic"; }
};
} // end namespace oofem
