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

#ifndef nlstructuralsurfaceelement_h
#define nlstructuralsurfaceelement_h

#include "sm/Elements/nlstructuralelement.h"

///@name Input fields for NLStructuralSurfaceElement
//@{
#define _IFT_NLStructuralSurfaceElement_nlgeoflag "nlgeo"
//@}

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Abstract base class for structural surface finite elements with geometrical nonlinearities.
 *
 * The deformation gradient is F = BH * u + I_s with I_s = diag(1, 1, 0): computeBHmatrixAt must
 * return the displacement gradient in the element's reference surface frame: the first two axes
 * are tangent to the reference surface and the third is its reference normal. The kinematic F
 * and the work-conjugate P are stored in this frame, so the stored P matches BH (e.g. for reactions).
 * A constitutive call may use F*F0 to include prestrain; the returned P is then pulled back to F.
 * F and P are 9-component Voigt vectors (11, 22, 33, 23, 13, 12, 32, 31, 21) and BH has 9 rows.
 *
 * The material receives F*F0 together with the unit normal N of the configuration this deformation
 * gradient starts from, in the same frame and with (F*F0) N = 0 (see giveMaterialKinematics): e3 without
 * prestrain, the normal of the original configuration with it. These elements interpolate C0 only, so a
 * material which needs the second gradient (bending stiffness) is rejected in checkConsistency.
 */
class NLStructuralSurfaceElement : public NLStructuralElement
{
protected:

public:

    NLStructuralSurfaceElement( int n, Domain *d );

    virtual ~NLStructuralSurfaceElement() {}

    void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;

    void computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    void computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    /**
     * Kinematics handed to the surface material, in the frame of the element.
     * @param vF Kinematic deformation gradient F = BH u + I_s, F e3 = 0.
     * @param vFmat F * F0 with the prestrain F0 of the element (none: F0 = diag(1, 1, 0) and vFmat = F).
     * @param normal Unit normal N of the configuration F * F0 starts from, (F * F0) N = 0: the normal of
     * F0s = I_s F0, the part of F0 in the element plane, N = cof(F0s)^T e3 / J0 (e3 without prestrain).
     * @param F0s The part of F0 in the element plane, used to pull the stress and the tangent back to F.
     * @param J0 Area ratio of F0s, J0 = |cof(F0s)^T e3| (1 without prestrain).
     */
    void giveMaterialKinematics( const FloatArray &vF, GaussPoint *gp, TimeStep *tStep, FloatArray &vFmat, FloatArray &normal, Tensor2_3d &F0s, double &J0 );

    int checkConsistency() override;

    /// Unit normal to the reference geometry, expressed in global coordinates (the geometry of the element,
    /// used e.g. for the anisotropy direction). The material gets the normal of giveMaterialKinematics instead.
    virtual FloatArray giveNormal( GaussPoint *gp ) const = 0;

    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override
    {
        return false;
    };
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override
    {
        return 0;
    };
};
} // end namespace oofem
#endif // nlstructuralsurfaceelement_h
