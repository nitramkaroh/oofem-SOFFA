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

#ifndef iga_h
#define iga_h

#include "inputrecord.h"
#include "intarray.h"
#include "feinterpol.h"
#include "gaussintegrationrule.h"
#include "vtkxmlexportmodule.h"
#include "internalstatetype.h"

///@name Input fields for IGAElement
//@{
#define _IFT_IGAElement_KnotSpanParallelMode "knotspanparmode"
#define _IFT_IGAElement_numberOfExportPoints "numexp" ///< VTK export points per knot span and direction (optional, default 2, at least 2).
//@}

namespace oofem {
#ifdef __OOFEG
class StructuralElementEvaluator;
void drawIGAPatchDeformedGeometry(Element *elem, StructuralElementEvaluator *se, oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
#endif

/**
 * Geometry wrapper for IGA elements: the control points are the element nodes.
 */
class OOFEM_EXPORT FEIIGAElementGeometryWrapper : public FEICellGeometry
{
public:
    /**
     * Knot span in each parametric direction, as the zero based index into the knot vector with
     * multiplicities (see BSplineInterpolation::findSpan); nullptr = the interpolation finds the span
     * from the local coordinates (at an interior knot the span on its right).
     */
    const IntArray *knotSpan;
    Element *elem;
public:
    FEIIGAElementGeometryWrapper(Element *elem, const IntArray *knotSpan=nullptr) : FEICellGeometry(), knotSpan(knotSpan), elem(elem) { }

    int giveNumberOfVertices() const override { return elem->giveNumberOfNodes(); }
    const FloatArray &giveVertexCoordinates(int i) const override { return elem->giveNode(i)->giveCoordinates(); }
};


/**
 * IntegrationElement represent nonzero knot span, derived from Integration Rule.
 */
class OOFEM_EXPORT IGAIntegrationElement : public GaussIntegrationRule
{
protected:
    IntArray knotSpan;     // knot_span(nsd)
public:
    IGAIntegrationElement(int _n, Element *e, IntArray knotSpan) :
        GaussIntegrationRule(_n, e, 0, 0, false),
        knotSpan(std::move(knotSpan)) { }
    const IntArray *giveKnotSpan() override { return & this->knotSpan; }

    /**
     * Creates nPoints points (on a line), floor(sqrt(nPoints))^2 points (on a square) or n^3 points with the
     * largest n such that n^3 <= nPoints (in a cube) with zero coordinates and weights, to be placed by the
     * caller (the VTK export points of IGAElement).
     * @return Number of created points.
     */
    int SetUpDummyPointsOnLine( int nPoints, MaterialMode mode );
    /// @copydoc SetUpDummyPointsOnLine
    int SetUpDummyPointsOnSquare( int nPoints, MaterialMode mode );
    /// @copydoc SetUpDummyPointsOnLine
    int SetUpDummyPointsOnCube( int nPoints, MaterialMode mode );
};


/**
 * Implements base IGAElement, supposed to be a parent class of all elements with B-spline or NURBS based interpolation.
 */
class OOFEM_EXPORT IGAElement : public Element, public VTKXMLExportModuleElementInterface
{
protected:
    /// VTK export points, one rule per nonzero knot span (B-spline and NURBS elements, not T-splines).
    std::vector<std ::unique_ptr<IGAIntegrationElement> > ExportIntegrationRulesArray;
    /// Number of export points per knot span and parametric direction (numexp).
    int numExportPoints = 2;
    /// Set when the warning about VTK hexahedra that do not follow the geometry was given (once per element).
    bool exportCellWarningGiven = false;

#ifdef __PARALLEL_MODE
    IntArray knotSpanParallelMode;
#endif
public:
    IGAElement( int n, Domain *aDomain ) :
        Element( n, aDomain ), VTKXMLExportModuleElementInterface()
    {}
    /**
     * Reads the nodes, material, cross section and the interpolation, creates one integration rule
     * with nip points per nonzero knot span (a x a points for nsd = 2, nip a square number) and the
     * VTK export points (numexp). nlgeo is read by the element evaluators.
     */
    void initializeFrom(InputRecord &ir) override;

    /// Element and the VTK export interface both declare it, so the IGAElement code needs this override; the elements give their own names.
    const char *giveClassName() const override { return "IGAElement"; }

    /// Unit normal of the reference geometry at the Gauss point (surface elements).
    virtual FloatArray giveNormal( GaussPoint *gp ) const {
        OOFEM_ERROR( "method not implemented for this element" );
        FloatArray normal;
        return normal;
    };

    /// The IGA interpolations need FEIIGAElementGeometryWrapper (Element's version passes FEIElementGeometryWrapper).
    int computeGlobalCoordinates( FloatArray &answer, const FloatArray &lcoords ) override;

    /// Composite geometry, exported by giveCompositeExportData.
    Element_Geometry_Type giveGeometryType() const override { return EGT_Composite; }
    /**
     * VTK export of a curve (nsd 1), a patch with two parametric directions (nsd 2, planar or a surface in
     * space) or a patch with three (nsd 3) sampled at the export points: numexp points per knot span and
     * parametric direction, u index fastest, giving numexp - 1 lines, (numexp - 1)^2 quads or (numexp - 1)^3
     * hexahedra per knot span. A hexahedron with a negative volume of its own (a left-handed parametrization)
     * is written top face first, so that its volume is positive in the VTK convention. If the Jacobian of a
     * hexahedron, sampled at its corners, Gauss points and centre, changes sign or vanishes inside, the export
     * points are too coarse for the patch: a warning (once per element) asks for a larger numexp. The
     * positions and DisplacementVector (padded to 3 components) are interpolated from the control points of
     * the span in the directions of the dofs of the element (u, v or u, v, w), the internal variables come from
     * giveExportPointValue. Other primary variables, unsupported internal variables and cell variables stop
     * with an error.
     */
    void giveCompositeExportData( std::vector<VTKPiece> &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;
    /**
     * Value of the internal variable type at the export point gp of the knot span ir (index into the export
     * rules), as written to the VTK file (tensors as 9 components, row by row).
     * @return false if the element does not export the variable (the default for all).
     */
    virtual bool giveExportPointValue( FloatArray &answer, InternalStateType type, GaussPoint *gp, int ir, TimeStep *tStep ) { return false; }
    Interface *giveInterface( InterfaceType it ) override;

#ifdef __PARALLEL_MODE
    elementParallelMode giveKnotSpanParallelMode(int) const override;
#endif

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

protected:
    /// Number of parametric directions (the nsd of the interpolation).
    virtual int giveNsd() = 0;
};


/**
 * IGATSplineElement setups integration rules differently from IGAElement.
 */
class OOFEM_EXPORT IGATSplineElement : public IGAElement
{
public:
    IGATSplineElement(int n, Domain * aDomain) : IGAElement(n, aDomain) { }
    void initializeFrom(InputRecord &ir) override;
};
} // end namespace oofem
#endif //iga_h
