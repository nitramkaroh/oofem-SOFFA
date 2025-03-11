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

#include "hardmagneticmooneyrivlincompressiblematerial.h"

#define _IFT_HardMagneticMooneyRivlinCompressibleMaterialPrestrain_Name "hardmagneticmooneyrivlincompressiblematprestrain"
#define _IFT_HardMagneticMooneyRivlinCompressibleMaterialPrestrain_magEndTime "magendtime"

namespace oofem {
/**
 * This class implements Hard Magnetic Compressible Mooney - Rivlin material with prestrain considerations.
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 *
 */
class HardMagneticMooneyRivlinCompressibleMaterialPrestrain : public HardMagneticMooneyRivlinCompressibleMaterial
{
  protected:
    double magEndTime;

  public:
    HardMagneticMooneyRivlinCompressibleMaterialPrestrain( int n, Domain *d );

    void initializeFrom( InputRecord &ir ) override{
      HardMagneticMooneyRivlinCompressibleMaterial::initializeFrom( ir );
      magEndTime = 0.;
      IR_GIVE_OPTIONAL_FIELD( ir, magEndTime, _IFT_HardMagneticMooneyRivlinCompressibleMaterialPrestrain_magEndTime );
    }

    std::tuple<FloatArrayF<9>, FloatArrayF<3> > give_FirstPKStressVector_MagneticInduction_3d( const FloatArrayF<9> &vF, const FloatArrayF<3> &vH, GaussPoint *gp, TimeStep *tStep ) const override;
    std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 3>, FloatMatrixF<3, 9>, FloatMatrixF<3, 3> > giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) override;

    const char *giveInputRecordName() const override { return _IFT_HardMagneticMooneyRivlinCompressibleMaterialPrestrain_Name; }
    const char *giveClassName() const override { return "HardMagneticMooneyRivlinCompressibleMaterialPrestrain"; }
};
} // end namespace oofem
