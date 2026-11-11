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

namespace oofem {

/**
 * Evaluator of a surface energy on a planar curve (interpolation with one parametric direction,
 * control points in the xy plane, dofs u, v), for the boundary of a plane strain body: per unit length
 * of the curve and unit thickness.
 * F, P and grad F are in the global basis (x, y in the plane, z out of plane), F has 9 components
 * (11, 22, 33, 23, 13, 12, 32, 31, 21); F_33 = 1 comes from the surface identity.
 * The B-spline derivatives dN/dx are surface gradients, dN/dS times the unit tangent.
 */
class Space2dStructuralSurfaceElementEvaluator : public StructuralSurfaceElementEvaluator
{
public:
    Space2dStructuralSurfaceElementEvaluator() :
        StructuralSurfaceElementEvaluator() {}

protected:
    /**
     * Assemble interpolation matrix at given IP.
     * In case of IGAElements, N is assumed to contain only nonzero interpolation functions.
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
     * order of F; only the tangential part d^2 N / dS^2 (t x t) is assembled.
     */
    virtual void computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /// Length of the reference curve around gp (Jacobian times weight).
    double computeVolumeAround( GaussPoint *gp ) override;
    void giveDofManDofIDMask( int inode, IntArray &answer ) const
    {
        answer = { D_u, D_v };
    }

public:
    void initializeFrom( InputRecord &ir ) override;

    /// Code numbers of the knot span ie; a curve has nsd 1 but two dofs per control point.
    int giveIntegrationElementLocalCodeNumbers( IntArray &answer, Element *elem, IntegrationRule *ie ) override;

    /**
     * Curvature K (temporary value stored by the material at the Gauss points of the integration rule
     * RuleIndex) at lcoords, interpolated linearly between the two nearest Gauss points; used for the VTK
     * export. A rule with one point gives its value.
     */
    void computeReferenceCurvature( FloatArray &RefCurvature, TimeStep *tStep, FloatArray &lcoords, int RuleIndex ) override;
};

/**
 * Evaluator of a surface energy on the surface of revolution of a curve in the xy plane: x is the radius,
 * y the axis of rotation and z the hoop direction, measured per radian. As the plane version, with the
 * hoop stretch F_33 = 1 + u / r and the hoop terms of the second gradient.
 */
class AxisymSurfaceElementEvaluator : public Space2dStructuralSurfaceElementEvaluator
{
public:
    AxisymSurfaceElementEvaluator() :
        Space2dStructuralSurfaceElementEvaluator() {}

protected:
    /// As the plane version, with row 3 (hoop, u / r) added.
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /// Second gradient matrix 27 x (2 nno) including the hoop terms (Christoffel symbols of the surface of revolution).
    virtual void computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /// Area per radian around gp, r times the length.
    double computeVolumeAround( GaussPoint *gp ) override;

public:
    void initializeFrom( InputRecord &ir ) override;
};


} // end namespace oofem
#endif 
