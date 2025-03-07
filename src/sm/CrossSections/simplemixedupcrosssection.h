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

#ifndef simplemixedupcrosssection_h
#define simplemixedupcrosssection_h

#include "sm/CrossSections/simplecrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"

///@name Input fields for SimpleCrossSection
//@{
#define _IFT_SimpleMixedUPCrossSection_Name "simplemixedupcs"
#define _IFT_SimpleMixedUPCrossSection_thick "thick"
#define _IFT_SimpleMixedUPCrossSection_MaterialNumber "material" ///< Material number for the bulk//@}

namespace oofem {
/**
 * Class implementing "simple" cross section model in finite element problem.
 */
class OOFEM_EXPORT SimpleMixedUPCrossSection : public SimpleCrossSection
{
protected:
    int materialNumber = 0;   ///< Material number

public:
    /**
     * Constructor.
     * @param n Cross section number.
     * @param d Associated domain.
     */
    SimpleMixedUPCrossSection(int n, Domain *d) : SimpleCrossSection(n, d) { }

    std::tuple<FloatArray, double> giveRealStresses(const FloatArray &strain, double p, GaussPoint *gp, TimeStep *tStep) const;
    std::tuple<FloatArrayF< 6 >, double> giveRealStress_3d(const FloatArrayF< 6 > &strain, double p, GaussPoint *gp, TimeStep *tStep) const;

    std::tuple<FloatArrayF< 4 >, double> giveRealStress_PlaneStrain(const FloatArrayF< 4 > &strain, double p, GaussPoint *gp, TimeStep *tStep) const;

    std::tuple<FloatMatrixF<6,6>, FloatArrayF<6>, double> giveSmallStrainMixedPressureConstitutiveMatrices_3d(MatResponseMode, GaussPoint *gp, TimeStep *tStep);

    std::tuple<FloatMatrixF<4,4>, FloatArrayF<4>, double> giveSmallStrainMixedPressureConstitutiveMatrices_PlaneStrain(MatResponseMode, GaussPoint *gp, TimeStep *tStep);


    std::tuple<FloatArray, double> giveFirstPKStresses(const FloatArray &reducedF, double p, GaussPoint *gp, TimeStep *tStep) const;


    std::tuple<FloatArrayF< 9 >, double > giveFirstPKStress_3d(const FloatArrayF< 9 > &vF, double p, GaussPoint *gp, TimeStep *tStep) const;

    std::tuple< FloatArrayF< 5 >, double > giveFirstPKStress_PlaneStrain(const FloatArrayF< 5 > &vF, double p, GaussPoint *gp, TimeStep *tStep) const;

    
    
    void giveSmallStrainStiffnessMatrix(FloatMatrix &answer_uu, FloatArray &answer_up, double& K, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);


    void giveLargeStrainStiffnessMatrix(FloatMatrix &answer_uu, FloatArray &answer_up, double& K, double pressure, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    

    std::tuple<FloatMatrixF<9,9>, FloatArrayF<9>, double> giveLargeStrainMixedPressureConstitutiveMatrices_3d(double pressure, MatResponseMode, GaussPoint *gp, TimeStep *tStep);


    std::tuple<FloatMatrixF<5,5>, FloatArrayF<5>, double> giveLargeStrainMixedPressureConstitutiveMatrices_PlaneStrain(double pressure, MatResponseMode, GaussPoint *gp, TimeStep *tStep);


    const char *giveClassName() const override { return "SimpleMixedUPCrossSection"; }
    const char *giveInputRecordName() const override { return _IFT_SimpleMixedUPCrossSection_Name; }


};
} // end namespace oofem
#endif // simplecrosssection_h
