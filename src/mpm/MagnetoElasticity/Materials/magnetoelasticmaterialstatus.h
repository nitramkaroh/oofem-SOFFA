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

#ifndef magnetoelasticmaterialstatus_h
#define magnetoelasticmaterialstatus_h

#include "matstatus.h"
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
class OOFEM_EXPORT MagnetoElasticMaterialStatus : public MaterialStatus
{
protected:
    /// Equilibrated first Piola-Kirchhoff stress vector
    FloatArray PVector;
    /// Temporary first Piola-Kirchhoff stress vector (to find balanced state)
    FloatArray tempPVector;
    /// Equilibrated deformation gradient in reduced form
    FloatArray FVector;
    /// Temporary deformation gradient in reduced form (to find balanced state)
    FloatArray tempFVector;

    /// Equilibrated second deformation gradient in reduced form
    FloatArray gradFVector;
    /// Temporary second deformation gradient in reduced form (to find balanced state)
    FloatArray tempGradFVector;

    /// Equilibrated first Piola-Kirchhoff stress vector
    FloatArray HVector;
    /// Temporary first Piola-Kirchhoff stress vector (to find balanced state)
    FloatArray tempHVector;
    /// Equilibrated deformation gradient in reduced form
    FloatArray BVector;
    /// Temporary deformation gradient in reduced form (to find balanced state)
    FloatArray tempBVector;

    //Magnetization-induced H tensor
    Tensor1_3d HmagTensor;


    
public:
    /// Constructor. Creates new StructuralMaterialStatus with IntegrationPoint g.
    MagnetoElasticMaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override{;}

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override{;}
    void restoreContext(DataStream &stream, ContextMode mode) override{;}
    /// Returns the const pointer to receiver's first Piola-Kirchhoff stress vector.
    const FloatArray &givePVector() const { return PVector; }
    /// Returns the const pointer to receiver's deformation gradient vector.
    const FloatArray &giveFVector() const { return FVector; }
    /// Returns the const pointer to receiver's second deformation gradient vector.
    const FloatArray &giveGradFVector() const { return gradFVector; }
    /// Returns the const pointer to receiver's first Piola-Kirchhoff stress vector.
    const FloatArray &giveHVector() const { return HVector; }
    /// Returns the const pointer to receiver's deformation gradient vector.
    const FloatArray &giveBVector() const { return BVector; }
    
    const FloatArray &giveTempPVector() const { return tempPVector; }
    /// Returns the const pointer to receiver's temporary deformation gradient vector.
    const FloatArray &giveTempFVector() const { return tempFVector; }
    const FloatArray &giveTempGradFVector() const { return tempGradFVector; }
    const FloatArray &giveTempHVector() const { return tempHVector; }
    /// Returns the const pointer to receiver's temporary deformation gradient vector.
    const FloatArray &giveTempBVector() const { return tempBVector; }

    /// Assigns tempPVector to given vector v
    void letTempPVectorBe(const FloatArray &v) { tempPVector = v; }
    /// Assigns tempFVector to given vector v
    void letTempFVectorBe(const FloatArray &v) { tempFVector = v; }
    /// Assigns tempGradFVector to given vector v
    void letTempGradFVectorBe( const FloatArray &v ) { tempGradFVector = v; }
    /// Assigns tempPVector to given vector v
    void letTempHVectorBe(const FloatArray &v) { tempHVector = v; }
    /// Assigns tempFVector to given vector v
    void letTempBVectorBe(const FloatArray &v) { tempBVector = v; }

    /// Assigns Htensor to given tensor t
    void letHmagTensorBe( const Tensor1_3d &t ) { HmagTensor = t; }
    /// Returns the const pointer to receiver's Hmag tensor
    const Tensor1_3d giveHmagTensor() const { return HmagTensor; }

    
    const char *giveClassName() const override { return "MagnetoElasticMaterialStatus"; }
};
} // end namespace oofem
#endif // structuralms_h
