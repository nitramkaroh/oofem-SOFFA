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
 * This class represent a new concept on how to define elements.
 * Traditionally, Elements are derived from problem specific class (structural element, for example)
 * define their interpolation and implement methods to evaluate shape function matrix, geometrical matrix, etc.
 * This new concept, here represented by  StructuralElementEvaluator and derived classes allows to
 * define all problem specific methods (including shape function matrix, geometrical matrix evaluation) to be defined
 * once for all elements of the same type (plane stress elements, space3d elements, etc) just relying on
 * services of finite element interpolation classes.
 * Definition of particular element is then done simply by deriving if from evaluator and providing interpolation.
 *
 * StructuralElementEvaluator - base class of all structural elements
 * Individual elements supposed to be derived from StructuralElementEvaluator and IGAElement
 * @todo{Class is missing much documentation.}
 */
class StructuralSurfaceElementEvaluator : public StructuralElementEvaluator
{
protected:

    StructuralSurfaceElementEvaluator();
    virtual ~StructuralSurfaceElementEvaluator(){}

protected:
    /**
     * Computes the first Piola-Kirchhoff stress tensor on Voigt format. This method will
     * be called if nlGeo = 1 and mode = TL. This method computes the deformation gradient F and passes
     * it on to the crossection which then asks for the stress from the material.
     * @note P is related to S through F*S.
     *
     * @param answer Computed stress vector in Voigt form.
     * @param gp Gauss point at which the stress is evaluated.
     * @param tStep Time step.
     */
    virtual void computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    void computeSecondOrderStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

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

    void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    void computeConstitutiveMatrix_dAddF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    void computeConstitutiveMatrix_dAdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    /**
     * Computes the deformation gradient in Voigt form at integration point ip and at time
     * step tStep. Computes the displacement gradient and adds an identitiy tensor.
     *
     * @param answer Deformation gradient vector
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    void computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    void computeGradientOfDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    void computeSurfaceIdentityGradient( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

public:
    virtual void initializeFrom( InputRecord &ir );
};
} // end namespace oofem
#endif // structuralelementevaluator_h
