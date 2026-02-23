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

#include "matstatus.h"
#include "structuralms.h"
#include "floatarray.h"
#include "tensor/tensor1.h"

namespace oofem {
class GaussPoint;
class Dictionary;
class Domain;

/**
 * This class implements a structural material status information. It is attribute of
 * gaussPoint. This is only an abstract class, for every instance of material class
 * there should be specialized derived class, which handles are history variables.
 *
 * This is a base class for all material statuses corresponding to materials derived from
 * structural material class.
 * It defines stress and strain vectors and their increments.
 * Functions for accessing these components are defined.
 *
 * Tasks:
 * This is abstract class - only basic functionality is supported like:
 * - maintaining and providing access to stress and strain vectors
 *   (including their increments)
 * - storing and restoring status on tape
 * - printingYourself()
 * - updating Yourself after a new equilibrium state has been reached.
 */
class OOFEM_EXPORT ThirdMediumMaterialStatus : public StructuralMaterialStatus
{
protected:

    /// Equilibrated second deformation gradient in reduced form
    FloatArray gradFVector;
    /// Temporary second deformation gradient in reduced form (to find balanced state)
    FloatArray tempGradFVector;
    /// Equilibrated second order stress in reduced form
    FloatArray TVector;
    /// Temporary second order stress (to find balanced state)
    FloatArray tempTVector;

    /// Equilibrated deformation gradient vector in the centroid in reduced form
    FloatArray FbarVector;
    /// Temporary deformation gradient vector in the centroid in reduced form (to find balanced state)
    FloatArray tempFbarVector;
    /// Equilibrated Pbar stress in reduced form
    FloatArray PbarVector;
    /// Temporary Pbar stress in reduced form
    FloatArray tempPbarVector;
    
public:
    /// Constructor. Creates new StructuralMaterialStatus with IntegrationPoint g.
    ThirdMediumMaterialStatus( GaussPoint *g );

    void printOutputAt(FILE *file, TimeStep *tStep) const override{;}

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override{;}
    void restoreContext(DataStream &stream, ContextMode mode) override{;}
 
    /// Returns the const pointer to receiver's second deformation gradient vector.
    const FloatArray &giveGradFVector() const { return gradFVector; }    
    const FloatArray &giveTempGradFVector() const { return tempGradFVector; }

    /// Returns the const pointer to receiver's second order stress vector
    const FloatArray &giveTVector() const { return TVector; }
    const FloatArray &giveTempTVector() const { return tempTVector; }

    /// Returns the const pointer to receiver's deformation gradient vector in the centroid.
    const FloatArray &giveFbarVector() const { return FbarVector; }
    const FloatArray &giveTempFbarVector() const { return tempFbarVector; }

    /// Returns the const pointer to receiver's Pbar stress vector in the centroid.
    const FloatArray &givePbarVector() const { return PbarVector; }
    const FloatArray &giveTempPbarVector() const { return tempPbarVector; }

    /// Assigns tempGradFVector to given vector v
    void letTempGradFVectorBe( const FloatArray &v ) { tempGradFVector = v; }
    /// Assigns tempTVector to given vector v
    void letTempTVectorBe(const FloatArray &v) { tempTVector = v; }
    /// Assigns tempFbarVector to given vector v
    void letTempFbarVectorBe( const FloatArray &v ) { tempFbarVector = v; }
    /// Assigns tempPbarVector to given vector v
    void letTempPbarVectorBe( const FloatArray &v ) { tempPbarVector = v; }

    /// Functions for MaterialStatusMapperInterface
    void copyStateVariables(const MaterialStatus &iStatus) override;

    const char *giveClassName() const override { return "ThirdMediumMaterialStatus"; }
};
} // end namespace oofem
