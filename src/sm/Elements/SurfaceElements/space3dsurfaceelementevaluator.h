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

#ifndef space3dsurfaceelementevaluator_h
#define space3dsurfaceelementevaluator_h

#include "sm/Elements/structuralsurfaceelementevaluator.h"
#include "internalstatetype.h"

namespace oofem {
/**
 * Evaluator of a surface energy on a surface in space (interpolation with two parametric directions and
 * control points with three coordinates, BSpline3dSurfaceInterpolation or NURBS3dSurfaceInterpolation; dofs u, v,
 * w), per unit area of the
 * reference surface. F and P are in the global basis with 9 components (11, 22, 33, 23, 13, 12, 32, 31, 21).
 * The derivatives dN/dx of the interpolation are surface gradients, so H = BH u is tangential in its second
 * index and F = I_s + H satisfies F N = 0. The second gradient (bending) is built from the tangential Hessian
 * of the interpolation, see computeGMatrixAt.
 */
class Space3dStructuralSurfaceElementEvaluator : public StructuralSurfaceElementEvaluator
{
public:
    Space3dStructuralSurfaceElementEvaluator() :
        StructuralSurfaceElementEvaluator() {}
    virtual ~Space3dStructuralSurfaceElementEvaluator() {}

    /**
     * IST_CurvatureTensor at the VTK export point gp, evaluated there from the displacements:
     * K = -n . grad F with the current unit normal n = cof(F) N / |cof(F) N| (as in the bending material),
     * 9 components row by row. An export point on a knot line uses the knot span of its export rule, so K
     * may jump between neighbouring spans (the second derivatives of a C1 basis do). Stops with an error
     * unless the element has a bending material.
     * @return false for the other types.
     */
    bool giveCurvatureExportValue( FloatArray &answer, InternalStateType type, GaussPoint *gp, TimeStep *tStep );

protected:
    /**
     * Interpolation matrix 3 x (3 nno), [N1 0 0 N2 ...; 0 N1 0 0 N2 ...; 0 0 N1 0 0 N2 ...], of the nno nonzero
     * basis functions of the knot span of gp.
     */
    void computeNMatrixAt( FloatMatrix &answer, GaussPoint *gp ) override;
    /**
     * Displacement gradient matrix, 9 x (3 nno) for the nno nonzero basis functions of the knot span,
     * the rows of the 9 component order filled with the surface gradient as in the 3D bulk evaluator.
     */
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /**
     * Matrix of the second gradient, 27 x (3 nno): grad F = G x with the nodal positions x. Component
     * (i, j, k) = d^2 x_i / dX_j dX_k has the index 9 (k - 1) + v(i, j), v the 9 component order of F, and
     * for the dof i of a node it is the tangential Hessian (j, k) of its basis function (evald2Ndx2).
     * The normal part of the full derivative (N in the slot j) is left out: the surface materials contract
     * grad F in j and k only with tangential tensors, so it does not contribute.
     */
    void computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /// Area of the reference surface around gp (|X_,xi x X_,eta| times the weight).
    double computeVolumeAround( GaussPoint *gp ) override;

    /// Dofs u, v, w of every control point.
    void giveDofManDofIDMask( int inode, IntArray &answer ) const
    {
        answer = { D_u, D_v, D_w };
    }

};
} // end namespace oofem
#endif // space3dsurfaceelementevaluator_h
