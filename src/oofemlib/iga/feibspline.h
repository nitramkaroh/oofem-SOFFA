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

#ifndef feibspline_h
#define feibspline_h

#include "feinterpol.h"
#include "floatarray.h"
#include <array>

///@name Input fields for BSplineInterpolation
//@{
#define _IFT_BSplineInterpolation_degree "degree"
#define _IFT_BSplineInterpolation_knotVectorU "knotvectoru"
#define _IFT_BSplineInterpolation_knotVectorV "knotvectorv"
#define _IFT_BSplineInterpolation_knotVectorW "knotvectorw"
#define _IFT_BSplineInterpolation_knotMultiplicityU "knotmultiplicityu"
#define _IFT_BSplineInterpolation_knotMultiplicityV "knotmultiplicityv"
#define _IFT_BSplineInterpolation_knotMultiplicityW "knotmultiplicityw"
//@}

namespace oofem {
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Interpolation for B-splines.
 */
class OOFEM_EXPORT BSplineInterpolation : public FEInterpolation
{
protected:
    /// Number of parametric directions (one for a curve, even when embedded in the xy plane).
    int nsd;
    /// Degree in each direction.
    std::array<int, 3> degree;                                   // eg. 2
    /// Knot values [nsd]
    std::array<FloatArray, 3> knotValues;                        // eg. 0 1 2 3 4 5
    /// Knot multiplicity [nsd]
    std::array<IntArray, 3> knotMultiplicity;                    // eg. 3 1 1 1 2 3
    /// numberOfControlPoints[nsd]
    /**
     * For TSpline this is filled by values corresponding to case when there are no T-junctions
     * (i.e. TSpline = BSpline)
     */
    std::array<int, 3> numberOfControlPoints;
    /// Knot vectors [nsd][knot_vector_size]
    std::array<FloatArray, 3> knotVector;                           // eg. 0 0 0 1 2 3 4 4 5 5 5
    /// Nonzero spans in each directions [nsd]
    std::array<int, 3> numberOfKnotSpans;                        // eg. 5 (0-1,1-2,2-3,3-4,4-5)
public:
    BSplineInterpolation(int nsd) : FEInterpolation(0),
        nsd(nsd)
    {}

    integrationDomain giveIntegrationDomain() const override {
        if ( nsd == 3 ) {
            return _Cube;
        } else if ( nsd == 2 ) {
            return _Square;
        } else if ( nsd == 1 ) {
            return _Line;
        } else {
            return _UnknownIntegrationDomain;
        }
    }
    Element_Geometry_Type giveGeometryType() const override { return EGT_unknown; }

    integrationDomain giveBoundaryIntegrationDomain(int ib) const override {
        if ( nsd == 3 ) {
            return _Square;
        } else if ( nsd == 2 ) {
            return _Line;
        } else if ( nsd == 1 ) {
            return _Point;
        } else {
            return _UnknownIntegrationDomain;
        }
    }
    integrationDomain giveBoundarySurfaceIntegrationDomain(int isurf) const override {
        if ( nsd == 3 ) {
            return _Square;
        } else if ( nsd == 2 ) {
            return _Square;
        } else {
            return _UnknownIntegrationDomain;
        }
    }
    integrationDomain giveBoundaryEdgeIntegrationDomain(int iedge) const override {
        if ( nsd >= 2 ) {
            return _Line;
        } else {
            return _UnknownIntegrationDomain;
        }
    }

    int giveNsd()  const override { return nsd; }
    void initializeFrom(InputRecord &ir) override;

    IntArray boundaryEdgeGiveNodes(int boundary) const override
    { OOFEM_ERROR("Functions not supported for this interpolator."); return {}; }
    void boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("Functions not supported for this interpolator."); }
    double boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("Functions not supported for this interpolator."); return 0.; }
    void boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("Functions not supported for this interpolator."); }


    /**@name Surface interpolation services */
    //@{
    void boundarySurfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("Functions not supported for this interpolator."); }
    void boundarySurfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("Functions not supported for this interpolator."); }
    double boundarySurfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("Functions not supported for this interpolator."); }
    void boundarySurfaceLocal2global(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("Functions not supported for this interpolator."); }
    double boundarySurfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("Functions not supported for this interpolator."); }
    IntArray boundarySurfaceGiveNodes(int boundary) const override
    { OOFEM_ERROR("Functions not supported for this interpolator."); return {}; }
    //@}


    IntArray boundaryGiveNodes(int boundary) const override
    { OOFEM_ERROR("Not implemented"); return {}; }
    void boundaryEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("Not implemented"); }
    double boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("Not implemented"); return 0.; }
    double boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("boundaryGiveTransformationJacobian - Not implemented"); return 0.; }
    void boundaryLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override
    { OOFEM_ERROR("boundaryLocal2Global - Not implemented"); }


    int giveNumberOfKnotSpans(int dim) const override { return numberOfKnotSpans [ dim - 1 ]; }
    virtual int giveNumberOfControlPoints(int dim) const { return numberOfControlPoints [ dim - 1 ]; }
    const FloatArray *giveKnotVector()  const override {
        return this->knotVector.data();
    }
    const IntArray *giveKnotMultiplicity(int dim) const override { return & this->knotMultiplicity [ dim - 1 ]; }
    const FloatArray *giveKnotValues(int dim) const override { return & this->knotValues [ dim - 1 ]; }
    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
        OOFEM_ERROR("Not yet implemented.");
        return 0;
    }
    void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int giveKnotSpanBasisFuncMask(const IntArray &knotSpan, IntArray &mask) const override;
    int giveNumberOfKnotSpanBasisFunctions(const IntArray &knotSpan) const override;

    const char *giveClassName() const { return "BSplineInterpolation"; }
    bool hasSubPatchFormulation() const override { return true; }

    std::unique_ptr<IntegrationRule> giveIntegrationRule(int order) const override
    { OOFEM_ERROR("Not supported."); return nullptr; }
    std::unique_ptr<IntegrationRule> giveBoundaryIntegrationRule(int order, int boundary) const override
    { OOFEM_ERROR("Not supported."); return nullptr; }
    std::unique_ptr<IntegrationRule> giveBoundaryEdgeIntegrationRule(int order, int boundary) const override
    { OOFEM_ERROR("Not supported."); return nullptr; }

    //////////////////////
    // Evaluate up to second derivatives of basis functions
    int evalSecondDerivatives(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;


protected:
    /**
     * Evaluates the nonvanishing basis functions of 1d BSpline (algorithm A2.2 from NURBS book)
     * @param span Knot span index (zero based).
     * @param u Value at which to evaluate.
     * @param p Degree.
     * @param U Knot vector.
     * @param N Computed p+1 nonvanishing functions (N_{span-p,p}-N_{span,p})
     * @warning Parameter u and span must be in a valid range.
     */
    void basisFuns(FloatArray &N, int span, double u, int p, const FloatArray &U) const;
    /**
     * Computes nonzero basis functions and their derivatives at u.
     *
     * For information on the algorithm, see A2.3 on p72 of
     * the NURBS book. The result is a matrix of size (n+1,p+1), using zero-based
     * derivative rows and local basis columns:
     * @f[
     * N^{(r)}_{span-p+j,p}(u) = \mathit{ders}(r,j),
     * \qquad r=0,\ldots,n,\quad j=0,\ldots,p.
     * @f]
     * Thus column j refers to the global basis function with index span-p+j.
     *
     * @param n Highest derivative order to compute.
     * @param u Parametric value.
     * @param span Knot span index (zero based).
     * @param p Degree.
     * @param U Knot vector.
     * @param ders Matrix containing the derivatives of the basis functions.
     *
     * @warning Parameters n, u and span must be in a valid range.
     */
    void dersBasisFuns(int n, double u, int span, int p, const FloatArray &U, FloatMatrix &ders) const;
    /**
     * Determines the knot span index (Algorithm A2.1 from the NURBS book)
     *
     * Returns the index k such that U[k] <= u < U[k+1]. At an internal knot this
     * selects the nonzero span on its right. The right endpoint u = U[n+1] is
     * assigned to the last span, k = n.
     *
     * @param n Number of control points - 1.
     * @param u Parametric value.
     * @param p Degree.
     * @param U Knot vector.
     * @return Span index at u (zero based).
     * @warning Parameters outside [U[p], U[n+1]] are rejected; no endpoint tolerance is applied.
     */
    int findSpan(int n, int p, double u, const FloatArray &U) const;

    

    /**
     * Returns the range of nonzero basis functions for given knot span and given degree.
     */
    void giveNonzeroBasisFuncInterval(int span, int deg, int &s, int &e) const {
        s = span - deg;
        e = span;
    }
};


class OOFEM_EXPORT BSpline2dLineInterpolation : public BSplineInterpolation
{
public:
    BSpline2dLineInterpolation( int nsd ) :
        BSplineInterpolation( nsd ){}
    void initializeFrom( InputRecord &ir ) override {BSplineInterpolation::initializeFrom( ir );}
    double evaldNdx( FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    void evald2Ndx2( FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    int evalDerivatives( int maxOrder, FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    double evaldNdxGeneral( FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const;

    void evald2Ndx2General( FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const;

    /**
     * Unit normal of the reference curve, answer = (dY/dxi, -dX/dxi, 0) / |X_,xi|: the tangent
     * rotated clockwise in the xy plane. The boundary number is ignored; the knot span is taken
     * from the geometry wrapper when it is given, otherwise it is found from lcoords.
     * @return Always 1, not the boundary transformation Jacobian |X_,xi| promised by FEInterpolation.
     * Current callers use only the normal vector.
     */
    double boundaryEvalNormal( FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    /// Position on the curve in the xy plane (both coordinates; BSplineInterpolation gives one per parametric dimension).
    void local2global( FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    void giveJacobianMatrixAt( FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

};

} // end namespace oofem
#endif // feibspline_h
