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

#ifndef structuralsurfaceelementevaluator_h
#define structuralsurfaceelementevaluator_h

#include "iga/iga.h"
#include "matresponsemode.h"
#include "sm/Elements/structuralelementevaluator.h"

namespace oofem {
/**
 * Base of the evaluators of IGA surface elements, which carry a surface energy (see StructuralElementEvaluator
 * for the concept). Only the large strain formulation (nlgeo 1) is implemented.
 * The surface deformation gradient is F = I_s + H with the displacement gradient H = BH u and the surface
 * identity I_s = I - N x N of the reference unit normal N (IGAElement::giveNormal), 9 components in the
 * order 11, 22, 33, 23, 13, 12, 32, 31, 21. The cross section has to be a SimpleSurfaceCrossSection; it
 * passes F and N (and for second gradient materials grad F and grad I) to the surface material, which
 * returns the first Piola-Kirchhoff stress P (9 components) and the third order stress A (27 components).
 * The element must provide the normal and the derived evaluator the BH matrix (9 rows) and, for the
 * second gradient, the G matrix (27 rows).
 *
 * The bending material (second gradient) stores F, grad F and grad I in the temporary status when it computes
 * A; P takes grad F and grad I from there (and checks that F is the same) and the tangents take all three, so
 * StructuralElementEvaluator evaluates A before P.
 */
class StructuralSurfaceElementEvaluator : public StructuralElementEvaluator
{
protected:

    StructuralSurfaceElementEvaluator();
    virtual ~StructuralSurfaceElementEvaluator(){}

public:
    /// Reads nlgeo like StructuralElementEvaluator and rejects everything but 1 (only large strains are implemented).
    void initializeFrom( InputRecord &ir ) override;

protected:
    /**
     * Computes the first Piola-Kirchhoff stress (9 components) from the surface deformation gradient F
     * and the reference normal of the element, through the SimpleSurfaceCrossSection.
     *
     * @param answer Computed stress vector in Voigt form.
     * @param gp Gauss point at which the stress is evaluated.
     * @param tStep Time step.
     */
    virtual void computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    /**
     * Third order stress A (27 components) conjugate to grad F, from F, grad F = G x (current positions),
     * grad I = G X (reference positions) and the reference normal. The bending material also stores F, grad F,
     * grad I and the curvature K in the temporary status.
     */
    void computeSecondOrderStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    // small strains (nlgeo 0) are not implemented for the surface elements
    void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override{
        OOFEM_ERROR( "Small strains not implemented" );
        return;
    };
    void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override{
        OOFEM_ERROR( "Small strains not implemented" );
        return;
    };

    void computeBMatrixAt( FloatMatrix &answer, GaussPoint *gp ) override{
        OOFEM_ERROR( "Small strains not implemented" );
        return;
    };

    // The tangents below come from the material, which takes F (and for bending grad F and grad I) from the
    // temporary status; initializeReferenceGradients supplies grad F and grad I of a fresh status first.

    /// dP/dF, 9 x 9.
    void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    /// dA/d(grad F), 27 x 27.
    void computeConstitutiveMatrix_dAddF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    /// dA/dF, 27 x 9 (dP/d(grad F) is its transpose).
    void computeConstitutiveMatrix_dAdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    /**
     * For a second gradient material: if the status is fresh (no grad I stored yet), stores grad F = grad I =
     * G X, the values of the reference configuration F = I_s, in the temporary status, so that a tangent
     * requested before the first stress evaluation finds them. Does nothing otherwise.
     */
    void initializeReferenceGradients( GaussPoint *gp, TimeStep *tStep );

    /**
     * Computes the surface deformation gradient F = I_s + BH u (9 components) at the Gauss point,
     * with I_s = I - N x N of the reference normal N.
     *
     * @param answer Deformation gradient vector
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    void computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    /**
     * grad I = G X (27 components), as grad F with the reference positions X, i.e. grad F of the reference
     * configuration. It carries the reference curvature: -N . grad I is the gradient of the reference normal
     * (the bending material takes the reference mean curvature from it, useh0 1).
     */
    void computeSurfaceIdentityGradient( FloatArray &answer, GaussPoint *gp, TimeStep *tStep );

};
} // end namespace oofem
#endif // structuralsurfaceelementevaluator_h
