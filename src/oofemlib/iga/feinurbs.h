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

#ifndef feinurbs_h
#define feinurbs_h

#include "feibspline.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "mathfem.h"

///@name Input fields for NURBSInterpolation
//@{
#define _IFT_NURBSInterpolation_weights "weights"
//@}

namespace oofem {
/**
 * Interpolation class for NURBS.
 */
class OOFEM_EXPORT NURBSInterpolation : public BSplineInterpolation
{
protected:
    /// Weight of each control point (finite and positive).
    FloatArray weights;

public:
    NURBSInterpolation(int nsd) : BSplineInterpolation(nsd) { }

    /// Reads the B-spline data and one weight per control point.
    void initializeFrom(InputRecord &ir) override;

    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
        OOFEM_ERROR("Not yet implemented.");
        return 0;
    }
    void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    /// Not implemented (the B-spline version would ignore the weights).
    void evalSecondDerivatives(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("second derivatives of NURBS basis functions are not implemented"); }

    const char *giveClassName() const { return "NURBSInterpolation"; }
};


/**
 * NURBS surface in space: the surface of BSpline3dSurfaceInterpolation with one weight w_A per control point
 * (keyword weights, u index fastest, finite and positive), i.e. with the rational basis functions
 * R_A = w_A N_A / W, W = sum_B w_B N_B, which represent conic sections such as the circles of a surface of
 * revolution exactly. All services of the base class (tangents, surface gradient, Hessian, normal, area element)
 * use the rational functions and their parametric derivatives up to the second order; the displacements are
 * interpolated with them as well.
 */
class OOFEM_EXPORT NURBS3dSurfaceInterpolation : public BSpline3dSurfaceInterpolation
{
protected:
    /// Weight of each control point, u index fastest, divided by the largest input weight (in (0, 1]).
    FloatArray weights;

    /**
     * The rational functions of evalBasis: with W, W_,a and W_,ab the sums of w_A times the B-spline values,
     * R = w N / W, R_,a = (w N_,a - R W_,a) / W and R_,ab = (w N_,ab - R_,a W_,b - R_,b W_,a - R W_,ab) / W.
     */
    void evalBasis( int nder, const FloatArray &lcoords, const FEICellGeometry &cellgeo, FloatMatrix &answer, IntArray &cps ) const override;

public:
    NURBS3dSurfaceInterpolation() : BSpline3dSurfaceInterpolation() {}

    /**
     * Reads the B-spline data (degree >= 1 in both directions) and one weight per control point (finite and
     * positive). The weights are divided by the largest one, which leaves the surface unchanged and keeps the
     * computation independent of a common factor; a ratio of the smallest to the largest weight below the normal
     * range of double is rejected.
     */
    void initializeFrom( InputRecord &ir ) override;
    /// The rational basis functions of the knot span of lcoords, u index fastest.
    void evalN( FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;
    /// Position on the surface (three coordinates) from the rational basis.
    void local2global( FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;
    /// The rational functions with their derivatives up to the second order, in the layout of BSplineInterpolation.
    void evalSecondDerivatives( FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    const char *giveClassName() const { return "NURBS3dSurfaceInterpolation"; }
};
} // end namespace oofem
#endif // feinurbs_h
