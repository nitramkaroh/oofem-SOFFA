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

///@name Input fields for IGAElement
//@{
#define _IFT_IGAElement_KnotSpanParallelMode "knotspanparmode"
#define _IFT_IGAElement_numberOfExportPoints "numexp"
//@}

namespace oofem {
#ifdef __OOFEG
class StructuralElementEvaluator;
void drawIGAPatchDeformedGeometry(Element *elem, StructuralElementEvaluator *se, oofegGraphicContext &gc, TimeStep *tStep, UnknownType);
#endif

/**
 * Geometry wrapper for IGA elements.
 */
class OOFEM_EXPORT FEIIGAElementGeometryWrapper : public FEICellGeometry
{
public:
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
    void setKnotSpan1(const IntArray &src) { this->knotSpan = src; }

    int SetUpDummyPointsOnLine( int nPoints, MaterialMode mode );
    int SetUpDummyPointsOnSquare( int nPoints, MaterialMode mode );
};


/**
 * Implements base IGAElement, supposed to be a parent class of all elements with B-spline or NURBS based interpolation.
 */
class OOFEM_EXPORT IGAElement : public Element, public VTKXMLExportModuleElementInterface
{
protected:
    std::vector<std ::unique_ptr<IGAIntegrationElement> > ExportIntegrationRulesArray;
    int numExportPoints = 2;

#ifdef __PARALLEL_MODE
    IntArray knotSpanParallelMode;
#endif
public:
    IGAElement( int n, Domain *aDomain ) :
        Element( n, aDomain ), VTKXMLExportModuleElementInterface()
    {}
    void initializeFrom(InputRecord &ir) override;

    virtual FloatArray giveNormal( GaussPoint *gp ) const { 
        OOFEM_ERROR( "method not implemented for this element" ); 
        FloatArray normal;
        return normal;
    };

    virtual Element *giveElement() { return this; };

    virtual void giveNodalCoordinates( FloatArray &NodalCoordinates ) const;


    // composite type - so we can do the postprocessing in the next function
    Element_Geometry_Type giveGeometryType() const override { return EGT_Composite; }
    // divide the beam into some small subelements
    void giveCompositeExportData( std::vector<VTKPiece> &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;
    Interface *giveInterface( InterfaceType it ) override;

#ifdef __PARALLEL_MODE
    elementParallelMode giveKnotSpanParallelMode(int) const override;
#endif

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

protected:
    virtual int giveNsd() = 0; // this info is available also from interpolation. Do we need it here ???
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
