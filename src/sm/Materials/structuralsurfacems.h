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

#ifndef structuralsurfacems_h
#define structuralsurfacems_h

#include "structuralms.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matstatmapperint.h"

namespace oofem {
class GaussPoint;
class Dictionary;
class Domain;

class OOFEM_EXPORT StructuralSurfaceMaterialStatus : public StructuralMaterialStatus
{
protected:

    
    FloatArray GVector; // Second gradient vector
    FloatArray tempGVector;

    FloatArray MVector; // Second order stress conjugate to Curvature tensor viewed from ref. cofnig (see Silhous)
    FloatArray tempMVector;

    FloatArray dIVector; // gradient of reference identity
    FloatArray tempdIVector;

    FloatMatrix D2eDK2_Matrix; // second order stiffness
    FloatMatrix tempD2eDK2_Matrix;


public:
    /// Constructor. Creates new StructuralMaterialStatus with IntegrationPoint g.
    StructuralSurfaceMaterialStatus( GaussPoint *g );

    void initTempStatus() override;
    void updateYourself( TimeStep *tStep ) override;

    const FloatArray &giveMVector() const { return MVector; }
    const FloatArray &giveTempMVector() const { return tempMVector; }
    const FloatArray &giveGVector() const { return GVector; }
    const FloatArray &giveTempGVector() const { return tempGVector; }
    const FloatArray &givedIVector() const { return dIVector; }
    const FloatArray &giveTempdIVector() const { return tempdIVector; }
    const FloatMatrix &giveD2eDK2_Matrix() const { return D2eDK2_Matrix; }
    const FloatMatrix &giveTempD2eDK2_Matrix() const { return tempD2eDK2_Matrix; }

    void letMVectorBe( const FloatArray &v ) { MVector = v; }
    void letTempMVectorBe( const FloatArray &v ) { tempMVector = v; }
    void letGVectorBe( const FloatArray &v ) { GVector = v; }
    void letTempGVectorBe( const FloatArray &v ) { tempGVector = v; }
    void letdIVectorBe( const FloatArray &v ) { dIVector = v; }
    void letTempdIVectorBe( const FloatArray &v ) { tempdIVector = v; }
    void letD2eDK2_MatrixBe( const FloatMatrix &v ) { D2eDK2_Matrix = v; }
    void lettempD2eDK2_MatrixBe( const FloatMatrix &v ) { tempD2eDK2_Matrix = v; }

    const char *giveClassName() const override { return "StructuralSurfaceMaterialStatus"; }

    /// Functions for MaterialStatusMapperInterface
    void copyStateVariables( const MaterialStatus &iStatus ) override;
};
} // end namespace oofem
#endif // structuralms_h
