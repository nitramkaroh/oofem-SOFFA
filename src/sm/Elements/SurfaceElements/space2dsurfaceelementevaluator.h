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

#ifndef space2dstructuraleurfaceelementevaluator_h
#define space2dstructuraleurfaceelementevaluator_h

#include "sm/Elements/structuralsurfaceelementevaluator.h"
#include "internalstatetype.h"

namespace oofem {

/**
 * Evaluator of a surface energy on a planar curve (interpolation with one parametric direction,
 * control points in the xy plane, dofs u, v), for the boundary of a plane strain body: per unit length
 * of the curve and unit thickness.
 * F, P and grad F are in the global basis (x, y in the plane, z out of plane), F has 9 components
 * (11, 22, 33, 23, 13, 12, 32, 31, 21); F_33 = 1 comes from the surface identity.
 * The B-spline derivatives dN/dx are surface gradients, dN/dS times the unit tangent.
 * The reference normal N = (t_y, -t_x, 0) of the element enters through the surface identity I_s = I - N x N.
 */
class Space2dStructuralSurfaceElementEvaluator : public StructuralSurfaceElementEvaluator
{
public:
    Space2dStructuralSurfaceElementEvaluator() :
        StructuralSurfaceElementEvaluator() {}

protected:
    /**
     * Interpolation matrix 2 x (2 nno), [N1 0 N2 0 ...; 0 N1 0 N2 ...], of the nno nonzero basis functions
     * of the knot span of gp.
     */
    void computeNMatrixAt( FloatMatrix &answer, GaussPoint *gp ) override;
    /**
     * Displacement gradient matrix, 9 x (2 nno) for the nno nonzero basis functions of the knot span:
     * rows 1 (du/dx), 2 (dv/dy), 6 (du/dy) and 9 (dv/dx) of the 9 component order.
     */
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /**
     * Matrix of the second gradient, 27 x (2 nno): grad F = G x with the nodal positions x.
     * Component (i, j, k) = d^2 x_i / dX_j dX_k has the index 9 (k - 1) + v(i, j), v the 9 component
     * order of F; only the tangential part d^2 N / dS^2 (t x t) is assembled. The part with the normal in
     * the slot j (from the derivative of the unit tangent along the curve) is left out: the surface materials
     * contract grad F in j and k only with tangential tensors, so it does not contribute.
     */
    virtual void computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /// Length of the reference curve around gp (Jacobian times weight).
    double computeVolumeAround( GaussPoint *gp ) override;
    /// Dofs u, v of every control point.
    void giveDofManDofIDMask( int inode, IntArray &answer ) const
    {
        answer = { D_u, D_v };
    }

public:
    /**
     * IST_CurvatureTensor evaluated directly at the VTK export point gp of the knot span RuleIndex, as
     * 9 components row by row; stops with an error unless the element has a bending material (only those
     * compute the curvature). @return false for the other types.
     */
    bool giveCurvatureExportValue( FloatArray &answer, InternalStateType type, GaussPoint *gp, int RuleIndex, TimeStep *tStep );
};

/**
 * Evaluator of a surface energy on the surface of revolution of a curve in the xy plane: x is the radius,
 * y the axis of rotation and z the hoop direction, measured per radian. As the plane version, with the
 * hoop stretch F_33 = 1 + u / r and the hoop terms of the second gradient. r is the reference radius
 * (x at the Gauss point), which must be positive.
 */
class AxisymSurfaceElementEvaluator : public Space2dStructuralSurfaceElementEvaluator
{
public:
    AxisymSurfaceElementEvaluator() :
        Space2dStructuralSurfaceElementEvaluator() {}

protected:
    /// The plane version with row 3 (hoop, u / r) added.
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /**
     * Second gradient matrix 27 x (2 nno): grad F = G x is the covariant (tangential) Hessian of the position x
     * on the surface of revolution (coordinates xi and the hoop angle theta, evaluated at theta = 0), with the
     * hoop terms given by the Christoffel symbols of the reference surface. Index order as in the plane
     * version; as there, the normal part is left out.
     */
    virtual void computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /// Area per radian around gp, the length of the plane version times the radius r.
    double computeVolumeAround( GaussPoint *gp ) override;

};


} // end namespace oofem
#endif 
