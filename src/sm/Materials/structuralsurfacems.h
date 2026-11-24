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

namespace oofem {
class GaussPoint;

/**
 * Status of the surface materials. Besides F, P and the Cauchy stress of StructuralMaterialStatus it holds, for the
 * second gradient, grad F and the gradient of the reference surface identity (27 components each, as left by the
 * last evaluation of the second order stress; the stress and the tangents of the bending material read them; empty
 * in a fresh status until the IGA surface evaluator sets both to the reference grad I before the first tangent), and
 * for postprocessing the curvature K = -n . grad F (9 components, stored by the bending material;
 * VTK evaluates K directly at the export point). The initial F is the
 * surface identity of the element frame: diag(1, 1, 0) for the Lagrange elements, I - N x N for the IGA elements.
 */
class OOFEM_EXPORT StructuralSurfaceMaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatArray GVector; ///< grad F, dF_ij/dX_k (27 components, see HyperElasticSurfaceMaterial for the order)
    FloatArray tempGVector; ///< grad F of the current iteration

    FloatArray KVector; ///< Current curvature K = -n . grad F, in the reference frame (9 components)
    FloatArray tempKVector; ///< K of the current iteration

    FloatArray dIVector; ///< Gradient of the reference surface identity (27 components)
    FloatArray tempdIVector; ///< Gradient of the reference surface identity of the current iteration

public:
    /// Constructor. Creates new StructuralSurfaceMaterialStatus with IntegrationPoint g.
    StructuralSurfaceMaterialStatus( GaussPoint *g );

    /// Also copies the committed K, grad F, grad I and the surface Cauchy stress to the temporary ones.
    void initTempStatus() override;
    /// Also commits K, grad F, grad I and the surface Cauchy stress.
    void updateYourself( TimeStep *tStep ) override;
    /// Prints all nine components of F and P, in the frame of the material evaluation (see above).
    void printOutputAt( FILE *file, TimeStep *tStep ) const override;

    /// Also stores the committed second gradients and curvature.
    void saveContext( DataStream &stream, ContextMode mode ) override;
    /// Also restores the committed second gradients and curvature and their temporary values.
    void restoreContext( DataStream &stream, ContextMode mode ) override;

    const FloatArray &giveKVector() const { return KVector; }
    const FloatArray &giveTempKVector() const { return tempKVector; }
    const FloatArray &giveGVector() const { return GVector; }
    const FloatArray &giveTempGVector() const { return tempGVector; }
    const FloatArray &givedIVector() const { return dIVector; }
    const FloatArray &giveTempdIVector() const { return tempdIVector; }

    void letTempKVectorBe( const FloatArray &v ) { tempKVector = v; }
    void letTempGVectorBe( const FloatArray &v ) { tempGVector = v; }
    void letTempdIVectorBe( const FloatArray &v ) { tempdIVector = v; }

    const char *giveClassName() const override { return "StructuralSurfaceMaterialStatus"; }

    /// Functions for MaterialStatusMapperInterface
    void copyStateVariables( const MaterialStatus &iStatus ) override;
};
} // end namespace oofem
#endif // structuralsurfacems_h
