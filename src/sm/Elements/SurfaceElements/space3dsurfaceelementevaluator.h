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

namespace oofem {
/**
 * Evaluator of a surface energy on a surface in space (interpolation with two parametric directions and
 * control points with three coordinates, BSpline3dSurfaceInterpolation; dofs u, v, w), per unit area of the
 * reference surface. F and P are in the global basis with 9 components (11, 22, 33, 23, 13, 12, 32, 31, 21).
 * The derivatives dN/dx of the interpolation are surface gradients, so H = BH u is tangential in its second
 * index and F = I_s + H satisfies F N = 0. Membranes only: the second gradient (bending) is not implemented.
 */
class Space3dStructuralSurfaceElementEvaluator : public StructuralSurfaceElementEvaluator
{
public:
    Space3dStructuralSurfaceElementEvaluator() :
        StructuralSurfaceElementEvaluator() {}
    virtual ~Space3dStructuralSurfaceElementEvaluator() {}

protected:
    /**
     * Assemble interpolation matrix at given IP
     * In case of IGAElements, N is assumed to contain only nonzero interpolation functions
     */
    void computeNMatrixAt( FloatMatrix &answer, GaussPoint *gp ) override;
    /**
     * Displacement gradient matrix, 9 x (3 nno) for the nno nonzero basis functions of the knot span,
     * the rows of the 9 component order filled with the surface gradient as in the 3D bulk evaluator.
     */
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /// Area of the reference surface around gp (|X_,xi x X_,eta| times the weight).
    double computeVolumeAround( GaussPoint *gp ) override;

    void giveDofManDofIDMask( int inode, IntArray &answer ) const
    {
        answer = { D_u, D_v, D_w };
    }

};
} // end namespace oofem
#endif // space3dsurfaceelementevaluator_h
