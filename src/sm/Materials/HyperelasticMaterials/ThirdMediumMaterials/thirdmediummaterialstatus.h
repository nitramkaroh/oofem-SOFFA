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

    /// Assigns tempGradFVector to given vector v
    void letTempGradFVectorBe( const FloatArray &v ) { tempGradFVector = v; }

    const char *giveClassName() const override { return "ThirdMediumMaterialStatus"; }
};
} // end namespace oofem
