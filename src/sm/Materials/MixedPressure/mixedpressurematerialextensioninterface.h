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

#ifndef mixedpressurematerialextensioninterface_h
#define mixedpressurematerialextensioninterface_h

#include "interface.h"
#include "matresponsemode.h"
#include "domain.h"
///@name micromorphicmaterialextensioninterface
//@{

//@}

namespace oofem {
class FloatMatrix;
class FloatArray;
template <std::size_t N> class FloatArrayF;
template <std::size_t N, std::size_t M> class FloatMatrixF;
class GaussPoint;
class TimeStep;



/**
 * Material interface for gradient material models.
 */
class MixedPressureMaterialExtensionInterface : public Interface
{
protected:
    Domain *dom = nullptr;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    MixedPressureMaterialExtensionInterface(Domain *d);
    /// Destructor.
    virtual ~MixedPressureMaterialExtensionInterface() { }
};


class SmallStrainMixedPressureMaterialExtensionInterface : public MixedPressureMaterialExtensionInterface
{
protected:

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    SmallStrainMixedPressureMaterialExtensionInterface(Domain *d);
    /// Destructor.
    virtual ~SmallStrainMixedPressureMaterialExtensionInterface() { }
  // move to cross-section 
  //void giveRealStressVector(FloatArray &answer, double &volChange, GaussPoint *gp, const FloatArray &reducedStrain, double pressure, TimeStep *tStep);
  //void giveConstitutiveMatrices(FloatMatrix &answer, FloatArray &delta, double K, GaussPoint *gp,  TimeStep *tStep); 
  
  virtual std::tuple<FloatMatrixF<6,6>, FloatArrayF<6>, double> giveSmallStrainMixedPressureConstitutiveMatrices_3d(MatResponseMode, GaussPoint *gp, TimeStep *tStep) = 0;
  virtual std::tuple<FloatMatrixF<4,4>, FloatArrayF<4>, double> giveSmallStrainMixedPressureConstitutiveMatrices_PlaneStrain(MatResponseMode, GaussPoint *gp, TimeStep *tStep);
  //
  virtual std::tuple<FloatArrayF<6>, double> giveRealStressVector_3d(const FloatArrayF<6> &reducedStrain, double pressure, GaussPoint *gp,TimeStep *tStep)  = 0;
  virtual std::tuple<FloatArrayF<4>, double> giveRealStressVector_PlaneStrain(const FloatArrayF<4> &reducedStrain, double pressure, GaussPoint *gp,TimeStep *tStep);

};


class LargeStrainMixedPressureMaterialExtensionInterface : public MixedPressureMaterialExtensionInterface
{
protected:

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    LargeStrainMixedPressureMaterialExtensionInterface(Domain *d);
    /// Destructor.
    virtual ~LargeStrainMixedPressureMaterialExtensionInterface() { }

  // move to cross-section 
  //void giveFirstPKStressVector(FloatArray &answer, double &volChange, GaussPoint *gp, const FloatArray &reducedStrain, double pressure, TimeStep *tStep);
  //void giveConstitutiveMatrices_dPdF(FloatMatrix &answer, FloatArray &delta, double K, GaussPoint *gp,  TimeStep *tStep); 
  
  
  virtual std::tuple<FloatArrayF<9>, double> giveFirstPKStressVector_3d(const FloatArrayF<9> &F, double pressure, GaussPoint *gp, TimeStep *tStep) = 0;

  virtual std::tuple<FloatArrayF<5>, double> giveFirstPKStressVector_PlaneStrain(const FloatArrayF<5> &F, double pressure, GaussPoint *gp,TimeStep *tStep);
  virtual std::tuple<FloatMatrixF<9,9>, FloatArrayF<9>, double> giveLargeStrainMixedPressureConstitutiveMatrices_3d(double pressure, MatResponseMode, GaussPoint *gp, TimeStep *tStep) = 0;
  virtual std::tuple<FloatMatrixF<5,5>, FloatArrayF<5>, double> giveLargeStrainMixedPressureConstitutiveMatrices_PlaneStrain(double pressure, MatResponseMode, GaussPoint *gp, TimeStep *tStep);  
};



  
}
#endif
