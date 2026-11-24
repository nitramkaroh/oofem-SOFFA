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
    /// Distinct knot values [nsd], divided by the length of the knot vector (not shifted)
    std::array<FloatArray, 3> knotValues;                        // eg. 0 0.2 0.4 0.6 0.8 1
    /// Knot multiplicity [nsd]
    std::array<IntArray, 3> knotMultiplicity;                    // eg. 3 1 1 1 2 3
    /// numberOfControlPoints[nsd]
    /**
     * For TSpline this is filled by values corresponding to case when there are no T-junctions
     * (i.e. TSpline = BSpline)
     */
    std::array<int, 3> numberOfControlPoints;
    /// Knot vectors with the multiplicities [nsd][knot_vector_size]
    std::array<FloatArray, 3> knotVector;                        // eg. 0 0 0 0.2 0.4 0.6 0.8 0.8 1 1 1
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

    /**
     * Nonzero basis functions of the knot span with their first and second derivatives with respect
     * to the parameters (nsd = 2 only): answer is count x 6 with the columns N, N_,u, N_,v, N_,uu,
     * N_,vv, N_,uv. The knot span is taken from the geometry wrapper when given, otherwise found from lcoords.
     * The control points of the rows are given by giveKnotSpanBasisFuncMask.
     * NURBS and T-splines stop with an error (the rational functions are not implemented).
     * Needs degree >= 2 in both directions (dersBasisFuns gives derivatives only up to the degree). It is not
     * checked here; the elements require it through supportsSecondGradient before they use the second gradient.
     */
    virtual void evalSecondDerivatives(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const;

    /**
     * Lowest continuity C^k of the basis inside the patch in parametric direction dim (one based):
     * the degree minus the largest multiplicity of an interior knot (the degree when there is none).
     * The second gradient needs C^1, i.e. interior knot multiplicities at most degree - 1.
     */
    int giveMinimalContinuity(int dim) const;


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
     * @param n Highest derivative order to compute (at most p).
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
};


/**
 * B-spline curve in the xy plane: one parametric direction xi (nsd = 1) of degree p >= 1 with the control
 * points (element nodes) in the xy plane. Used by the IGA surface elements and the IGA clamp.
 * Derivatives with respect to x, y are the surface (tangential) derivatives along the reference curve,
 * with the tangent X_,xi, the metric g_11 = |X_,xi|^2, g^11 = 1/g_11 and the unit tangent t.
 * The cell geometry must be an FEIIGAElementGeometryWrapper. Its knot span is used when given (at an
 * interior knot it selects the left or right derivative), otherwise the span is found from lcoords.
 * The rows of the answers are the p+1 nonzero basis functions of the span (giveKnotSpanBasisFuncMask).
 * The methods that need the tangent stop with an error for a degenerate curve (|X_,xi|^2 zero, not finite
 * or outside the normal range of doubles, e.g. coinciding control points).
 */
class OOFEM_EXPORT BSpline2dLineInterpolation : public BSplineInterpolation
{
public:
    BSpline2dLineInterpolation() : BSplineInterpolation( 1 ) {}

    /// Reads the B-spline data and requires degree >= 1 (a curve with a tangent).
    void initializeFrom( InputRecord &ir ) override;

    /**
     * Basis functions and their derivatives with respect to xi up to maxOrder >= 0: answer is
     * (maxOrder+1) x (p+1), row r holds the r-th derivatives of the nonzero functions of the span.
     */
    void evalDerivatives( int maxOrder, FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const;

    /**
     * Surface gradient of the basis functions, dN/dX = N' g^11 X_,xi = (N'/|X_,xi|) t: answer is (p+1) x 2 (x, y).
     * @return |X_,xi|.
     */
    double evaldNdx( FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    /**
     * Intrinsic (covariant) Hessian of the basis functions along the curve, d2N/ds2 t x t with
     * d2N/ds2 = ( N'' - Gamma N' ) g^11 and Gamma = g^11 X_,xi . X_,xixi: answer is (p+1) x 3 (xx, yy, xy,
     * no engineering factor 2). Both slots are tangential. The full derivative of the vector field grad N
     * along the curve has the additional term (N'/|X_,xi|) dt/ds x t, which is not included; the IGA surface
     * elements build their second gradient from this projected Hessian.
     */
    void evald2Ndx2( FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    /**
     * Unit normal of the whole reference curve at xi, answer = (Y_,xi, -X_,xi) / |X_,xi| (two components):
     * the tangent rotated clockwise in the xy plane, i.e. opposite to LinePlaneStrainSurface::giveNormal.
     * @return |X_,xi|.
     */
    double evalCurveNormal( FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const;

    /**
     * Adapter of FEInterpolation: the normal of the whole curve (evalCurveNormal) with a third, zero
     * component, not the normal at an end point of the curve. The boundary number is ignored.
     * @return |X_,xi|.
     */
    double boundaryEvalNormal( FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    /// Position on the curve in the xy plane (both coordinates; BSplineInterpolation gives one per parametric dimension).
    void local2global( FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    /// 1 x 1 matrix with |X_,xi|.
    void giveJacobianMatrixAt( FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

private:
    /// Knot span of lcoords: the one of the geometry wrapper when given, otherwise found from lcoords.
    int giveSpan( const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const;
    /**
     * Derivatives of the nonzero basis functions up to order nder (1 or 2) into ders, the tangent a = X_,xi and,
     * for nder = 2, b = X_,xixi of the reference curve. Stops for a degenerate curve (see the class description).
     */
    void evalTangent( int nder, const FloatArray &lcoords, const FEICellGeometry &cellgeo, FloatMatrix &ders,
        std::array<double, 2> &a, std::array<double, 2> &b ) const;
};


/**
 * B-spline surface in space: two parametric directions xi, eta (nsd = 2) of degrees p, q >= 1 with the control
 * points (element nodes, u index fastest) given by three coordinates. Used by the IGA surface element in 3D.
 * With the tangents A_1 = X_,xi and A_2 = X_,eta, the metric A_ab = A_a . A_b, its inverse A^ab and the dual
 * basis A^a = A^ab A_b, the derivatives with respect to x, y, z are the surface (tangential) gradients
 * N_,a A^a, the unit normal is A_1 x A_2 / |A_1 x A_2| and the area element is |A_1 x A_2|.
 * BSplineInterpolation with nsd = 2 is the planar patch: it uses only x and y of the control points.
 * The cell geometry must be an FEIIGAElementGeometryWrapper. Its knot span is used when given (at an interior
 * knot it selects the one sided derivatives), otherwise the span is found from lcoords. The rows of the answers
 * are the (p+1)(q+1) nonzero basis functions of the span, u index fastest (giveKnotSpanBasisFuncMask).
 * The methods that need the tangents stop with an error when the area |A_1 x A_2| is zero or not finite
 * (e.g. coinciding control points). The surface gradient uses reciprocal cross products rather than
 * inverse metric combinations to avoid cancellation on skew patches and unnecessary squared norms.
 * All services take the basis functions and their parametric derivatives from evalBasis, which
 * NURBS3dSurfaceInterpolation replaces by the rational functions.
 */
class OOFEM_EXPORT BSpline3dSurfaceInterpolation : public BSplineInterpolation
{
public:
    BSpline3dSurfaceInterpolation() : BSplineInterpolation( 2 ) {}

    /// Reads the B-spline data and requires degree >= 1 in both directions (a surface with two tangents).
    void initializeFrom( InputRecord &ir ) override;

    /**
     * Surface gradient of the basis functions, N_,xi A^1 + N_,eta A^2: answer is (p+1)(q+1) x 3 (x, y, z).
     * @return |A_1 x A_2|.
     */
    double evaldNdx( FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    /**
     * Tangential (covariant) Hessian of the basis functions on the surface,
     * (N_,ab - Gamma^c_ab N_,c) A^a x A^b with the Christoffel symbols Gamma^c_ab = A^c . X_,ab, i.e.
     * (N_,ab - X_,ab . grad_s N) A^a x A^b: answer is (p+1)(q+1) x 6 (xx, yy, zz, yz, xz, xy, no engineering
     * factor 2). Both slots are tangential. The full surface derivative of the vector field grad_s N has the
     * additional term N_,a B^a_b N x A^b with the curvature B of the reference surface, which is not included;
     * the IGA surface element builds its second gradient from this projected Hessian, like the curves.
     * Needs degree >= 2 in both directions.
     */
    void evald2Ndx2( FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    /**
     * Unit normal of the reference surface, A_1 x A_2 / |A_1 x A_2| (oriented by xi, eta with the right hand rule).
     * @return |A_1 x A_2|.
     */
    double evalSurfaceNormal( FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const;

    /// Position on the surface (three coordinates; BSplineInterpolation gives one per parametric dimension).
    void local2global( FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    /// Area element |A_1 x A_2|.
    double giveTransformationJacobian( const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

    /// Stops with an error: the tangents of a surface in space form no square matrix (see giveTransformationJacobian).
    void giveJacobianMatrixAt( FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const override;

protected:
    /// Knot spans of lcoords: those of the geometry wrapper when given, otherwise found from lcoords.
    std::array<int, 2> giveSpans( const FloatArray &lcoords, const FEICellGeometry &cellgeo ) const;
    /**
     * The nonzero basis functions of the knot span of lcoords and their derivatives with respect to xi and eta up
     * to the order nder (0, 1 or 2): answer is (p+1)(q+1) x 6 with the columns N, N_,xi, N_,eta, N_,xixi,
     * N_,etaeta, N_,xieta (the columns above the order nder are zero), rows u index fastest; cps gets the one
     * based numbers of their control points. Here the tensor products of the B-spline functions.
     */
    virtual void evalBasis( int nder, const FloatArray &lcoords, const FEICellGeometry &cellgeo, FloatMatrix &answer, IntArray &cps ) const;
    /**
     * The basis of evalBasis up to the order nder (1 or 2) and the tangents a1 = X_,xi, a2 = X_,eta of the
     * reference surface; for nder = 2 also b = { X_,xixi, X_,etaeta, X_,xieta }. Stops for a degenerate surface
     * (see the class description).
     * @return |a1 x a2|.
     */
    double evalTangents( int nder, const FloatArray &lcoords, const FEICellGeometry &cellgeo, FloatMatrix &basis, IntArray &cps,
        std::array<double, 3> &a1, std::array<double, 3> &a2, std::array<std::array<double, 3>, 3> *b = nullptr ) const;
    /// Reciprocal tangents d1 = A^1 = (a2 x N) / J and d2 = A^2 = (N x a1) / J with N = (a1 x a2) / J and J = |a1 x a2|.
    static void reciprocalTangents( const std::array<double, 3> &a1, const std::array<double, 3> &a2, double J,
        const FloatArray &lcoords, std::array<double, 3> &d1, std::array<double, 3> &d2 );
    /// Coordinates of control point i of the cell; a surface in space needs three (further ones are not used).
    static const FloatArray &giveControlPointCoordinates( const FEICellGeometry &cellgeo, int i );
};

} // end namespace oofem
#endif // feibspline_h
