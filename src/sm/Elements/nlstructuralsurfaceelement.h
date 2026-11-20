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

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Abstract base class of the Lagrange surface elements (geometrically nonlinear, nlgeo 1 only), which carry
 * a surface energy given by a SimpleSurfaceCrossSection with a surface material.
 *
 * Each element works in an orthonormal frame (e1, e2, e3) whose third axis e3 is the unit normal of the
 * reference surface; computeMaterialFrame gives it as a global-to-element rotation. The derived elements
 * rotate the nodal unknowns to their frame (computeGtoLRotationMatrix), and computeBHmatrixAt returns the
 * displacement gradient in this frame, so the deformation gradient is F = BH * u + I_s with
 * I_s = diag(1, 1, 0), and F e3 = 0. F and P are 9-component Voigt vectors (11, 22, 33, 23, 13, 12, 32, 31, 21)
 * and BH has 9 rows. The kinematic F and the work-conjugate P are stored in this frame, so the stored P
 * matches BH (internal forces and reactions computed from the status).
 *
 * The material receives F together with the unit normal N of the configuration F starts from, in the same
 * frame and with F N = 0 (see giveMaterialKinematics): N = e3 without prestrain. An element with a prestrain
 * F0 (givePrestrainDeformationGradient, only TriangularSurface) hands F * F0 and the normal of the original
 * configuration to the material and pulls the returned stress and tangent back to F. These elements
 * interpolate C0 only, so a material which needs the second gradient (bending stiffness) is rejected in
 * checkConsistency.
 */
class NLStructuralSurfaceElement : public NLStructuralElement
{
protected:
    /**
     * Frame of the surface material as a 3 x 3 global-to-element rotation: the rows are the axes e1, e2, e3
     * in global coordinates, e3 = the reference normal. It is the physical frame of computeLoadGToLRotationMtrx
     * (without nodal coordinate systems or master DOFs). The 2 x 2 rotation of a planar element (rows tangent
     * and in-plane normal) is extended to the material order: e1 = tangent, e2 = -z, e3 = in-plane normal.
     */
    void computeMaterialFrame( FloatMatrix &answer );
public:

    NLStructuralSurfaceElement( int n, Domain *d );

    virtual ~NLStructuralSurfaceElement() {}

    /// Small-strain stress, reference normal e3 (not used: checkConsistency requires nlgeo 1).
    void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;

    /**
     * First Piola-Kirchhoff stress P at gp, in the element frame (9 Voigt components). The material gets the
     * kinematics of giveMaterialKinematics. Without prestrain the components with the reference index 3
     * (Voigt 33, 23, 13) are zeroed, P I_s; with prestrain the material stress is pulled back, P F0s^T / J0.
     * The kinematic F and this P are stored in the status (temporary values).
     */
    void computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    /**
     * Tangent dP/dF (9 x 9, Voigt) at gp, in the element frame, at the kinematic F of the status (stored by
     * computeFirstPKStressVector). Without prestrain the rows and columns with the reference index 3
     * (Voigt 33, 23, 13) are zeroed; with prestrain the material tangent A at F * F0 is pulled back,
     * F0s_jm A_imkn F0s_ln / J0.
     */
    void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    /// F = BH * u + I_s at gp in the element frame (9 Voigt components); u = the unknowns in the element frame minus the initial displacements.
    void computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    /// Interpolate every translational DOF of a node, independently of the dimension of the parameter space.
    void computeNmatrixAt( const FloatArray &lcoords, FloatMatrix &answer ) override;
    /// Global coordinates of lcoords, interpolated from the global nodal coordinates (up to 3 components).
    int computeGlobalCoordinates( FloatArray &answer, const FloatArray &lcoords ) override;
    /**
     * Besides the values of NLStructuralElement: IST_DisplacementVector (global displacement at gp, 3 components)
     * and IST_GlobalDeformationGradientTensor (F * F0, or F without prestrain, rotated to global coordinates:
     * R^T (F F0) R with the rotation R of computeMaterialFrame).
     */
    int giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep ) override;

    /// Body forces use a rotation of their pointwise components. Spatial PointLoad localization is unsupported.
    void computeLoadVector( FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep ) override;

    /**
     * Kinematics handed to the surface material, in the frame of the element.
     * @param vF Kinematic deformation gradient F = BH u + I_s, F e3 = 0.
     * @param gp Integration point (for the prestrain of the element).
     * @param vFmat F * F0 with the prestrain F0 of the element (none: vFmat = F).
     * @param normal Unit normal N of the configuration F * F0 starts from, (F * F0) N = 0. With the part of F0
     * in the element plane, F0s = I_s F0: N = cof(F0s)^T e3 / J0, so that F0s N = 0 (e3 without prestrain).
     * @param F0s The part of F0 in the element plane, I_s F0, used to pull the stress and the tangent back to F
     * (I_s without prestrain).
     * @param J0 Area ratio of F0s, J0 = |cof(F0s)^T e3| (1 without prestrain).
     * @return Whether the element has a prestrain. Without it, F0s = I_s only drops the stress and tangent
     * components with the reference index 3 (Voigt 33, 23, 13), which BH does not use.
     */
    bool giveMaterialKinematics( const FloatArray &vF, GaussPoint *gp, FloatArray &vFmat, FloatArray &normal, Tensor2_3d &F0s, double &J0 );

    /// Requires a simplesurfacecs cross section and nlgeo 1, and rejects a material with bending stiffness (second gradient).
    int checkConsistency() override;

    /// Unit normal of the reference geometry in global coordinates, i.e. the axis e3 of the element frame (the
    /// in-plane normal of the planar elements). The material gets the normal of giveMaterialKinematics instead.
    virtual FloatArray giveNormal( GaussPoint *gp ) const = 0;

    /// No rotation by default; the derived elements return the rotation of their unknowns to the element frame.
    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override
    {
        return false;
    };
    /// No frame by default; the derived elements return the global-to-element rotation of a point (2 x 2 or 3 x 3).
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override
    {
        return 0;
    };
};
} // end namespace oofem
#endif // nlstructuralsurfaceelement_h
