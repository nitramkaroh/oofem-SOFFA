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

/**
 * Status of the surface materials. Besides F, P and the Cauchy stress of StructuralMaterialStatus it holds, for the
 * second gradient, grad F and the gradient of the reference surface identity (27 components each, as left by the
 * last evaluation of the second order stress; the stress and the tangents of the bending material read them), and
 * for postprocessing the curvature K = -n . grad F and the moment M = de/dK (9 components each, stored by the
 * bending material; the 2D surface evaluator interpolates K between Gauss points). The initial F is the
 * surface identity of the element frame: diag(1, 1, 0) for the Lagrange elements, I - N x N for the IGA elements.
 */
class OOFEM_EXPORT StructuralSurfaceMaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatArray GVector; ///< grad F
    FloatArray tempGVector;

    FloatArray MVector; ///< Second order stress M = de/dK conjugate to the curvature K
    FloatArray tempMVector;

    FloatArray KVector; ///< Current curvature K = -n . grad F, in the reference frame
    FloatArray tempKVector;

    FloatArray dIVector; ///< Gradient of the reference surface identity
    FloatArray tempdIVector;

    FloatMatrix D2eDK2_Matrix; ///< d2e/dK2 (9 x 9); no longer written by the materials, stays zero
    FloatMatrix tempD2eDK2_Matrix;

public:
    /// Constructor. Creates new StructuralSurfaceMaterialStatus with IntegrationPoint g.
    StructuralSurfaceMaterialStatus( GaussPoint *g );

    void initTempStatus() override;
    void updateYourself( TimeStep *tStep ) override;

    const FloatArray &giveMVector() const { return MVector; }
    const FloatArray &giveTempMVector() const { return tempMVector; }
    const FloatArray &giveKVector() const { return KVector; }
    const FloatArray &giveTempKVector() const { return tempKVector; }
    const FloatArray &giveGVector() const { return GVector; }
    const FloatArray &giveTempGVector() const { return tempGVector; }
    const FloatArray &givedIVector() const { return dIVector; }
    const FloatArray &giveTempdIVector() const { return tempdIVector; }
    const FloatMatrix &giveD2eDK2_Matrix() const { return D2eDK2_Matrix; }
    const FloatMatrix &giveTempD2eDK2_Matrix() const { return tempD2eDK2_Matrix; }

    void letMVectorBe( const FloatArray &v ) { MVector = v; }
    void letTempMVectorBe( const FloatArray &v ) { tempMVector = v; }
    void letKVectorBe( const FloatArray &v ) { KVector = v; }
    void letTempKVectorBe( const FloatArray &v ) { tempKVector = v; }
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
#endif // structuralsurfacems_h
