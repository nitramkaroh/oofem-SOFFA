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

#ifndef igasurfaceelements_h
#define igasurfaceelements_h

#include "iga/iga.h"
#include "iga/feibspline.h"
#include "iga/feinurbs.h"
#include "iga/feitspline.h"

#include "sm/Elements/SurfaceElements/space2dsurfaceelementevaluator.h"

#include "sm/Elements/3D/space3delementevaluator.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matresponsemode.h"
#include "mathfem.h"

#define _IFT_BsplinePlaneStrainSurfaceElement_Name "bsplineplanestrainsurfaceelement"
#define _IFT_BsplineAxisymSurfaceElement_Name "bsplineaxisymsurfaceelement"
namespace oofem {

/**
 * Surface energy on a planar B-spline curve (one parametric direction, control points in the xy plane),
 * typically the boundary of a plane strain B-spline patch sharing its control points. Per unit length
 * of the curve and unit thickness. Large strains only (nlgeo 1); the cross section is a
 * SimpleSurfaceCrossSection with a surface material.
 */
class BsplinePlaneStrainSurfaceElement : public IGAElement, public Space2dStructuralSurfaceElementEvaluator
{
protected:
    BSpline2dLineInterpolation interpolation;

public:
    BsplinePlaneStrainSurfaceElement( int n, Domain *aDomain );

    void initializeFrom( InputRecord &ir ) override;
    void postInitialize() override { IGAElement :: postInitialize(); this->setupEvaluator(); }
    int checkConsistency() override;

    void giveCharacteristicMatrix( FloatMatrix &answer, CharType mtrx, TimeStep *tStep ) override
    {
        Space2dStructuralSurfaceElementEvaluator ::giveCharacteristicMatrix( answer, mtrx, tStep );
    }
    void giveCharacteristicVector( FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t ) override
    {
        Space2dStructuralSurfaceElementEvaluator ::giveCharacteristicVector( answer, type, mode, t );
    }

    FEInterpolation *giveInterpolation() const override { return const_cast<BSpline2dLineInterpolation *>( &this->interpolation ); }
    Element *giveElement() override { return this; }
    void giveDofManDofIDMask( int inode, IntArray &answer ) const override
    {
        Space2dStructuralSurfaceElementEvaluator ::giveDofManDofIDMask( inode, answer );
    }
    int computeNumberOfDofs() override { return numberOfDofMans * 2; }
    void updateInternalState( TimeStep *tStep ) override { Space2dStructuralSurfaceElementEvaluator ::updateInternalState( tStep ); }
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_BsplinePlaneStrainSurfaceElement_Name; }
    const char *giveClassName() const override { return "BsplinePlaneStrainSurfaceElement"; }

    MaterialMode giveMaterialMode() override { return _3dMat; }

    /// Unit normal of the reference curve, (t_y, -t_x, 0) / |t| with the tangent t = dX/dxi.
    FloatArray giveNormal( GaussPoint *gp ) const override;

    /// VTK export of the curve sampled in numexp points per knot span (displacements only).
    void giveCompositeExportData( std::vector<VTKPiece> &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;

    /// The second gradient (bending) needs C1 continuity: degree at least 2 (repeated knots are not checked).
    bool supportsSecondGradient() const override { return this->interpolation.giveInterpolationOrder() >= 2; }

protected:
    int giveNsd() override { return 1; }
};


/**
 * Surface energy on the surface of revolution of a B-spline curve in the xy plane (x the radius,
 * y the axis of rotation), per radian. Large strains only (nlgeo 1); the cross section is a
 * SimpleSurfaceCrossSection with a surface material.
 */
class BsplineAxisymSurfaceElement : public IGAElement, public AxisymSurfaceElementEvaluator
{
protected:
    BSpline2dLineInterpolation interpolation;

public:
    BsplineAxisymSurfaceElement( int n, Domain *aDomain );

    void initializeFrom( InputRecord &ir ) override;
    void postInitialize() override { IGAElement :: postInitialize(); this->setupEvaluator(); }
    int checkConsistency() override;

    void giveCharacteristicMatrix( FloatMatrix &answer, CharType mtrx, TimeStep *tStep ) override
    {
        AxisymSurfaceElementEvaluator ::giveCharacteristicMatrix( answer, mtrx, tStep );
    }
    void giveCharacteristicVector( FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t ) override
    {
        AxisymSurfaceElementEvaluator ::giveCharacteristicVector( answer, type, mode, t );
    }

    FEInterpolation *giveInterpolation() const override { return const_cast<BSpline2dLineInterpolation *>( &this->interpolation ); }
    Element *giveElement() override { return this; }
    void giveDofManDofIDMask( int inode, IntArray &answer ) const override
    {
        AxisymSurfaceElementEvaluator ::giveDofManDofIDMask( inode, answer );
    }
    int computeNumberOfDofs() override { return numberOfDofMans * 2; }
    void updateInternalState( TimeStep *tStep ) override { AxisymSurfaceElementEvaluator ::updateInternalState( tStep ); }
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_BsplineAxisymSurfaceElement_Name; }
    const char *giveClassName() const override { return "BsplineAxisymSurfaceElement"; }

    MaterialMode giveMaterialMode() override { return _3dMat; }

    /// Unit normal of the reference curve, (t_y, -t_x, 0) / |t| with the tangent t = dX/dxi.
    FloatArray giveNormal( GaussPoint *gp ) const override;

    /**
     * VTK export of the curve sampled in numexp points per knot span: displacements and, as the only
     * internal variable, IST_CurvatureTensor (see computeReferenceCurvature).
     */
    void giveCompositeExportData( std::vector<VTKPiece> &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;

    /// The second gradient (bending) needs C1 continuity: degree at least 2 (repeated knots are not checked).
    bool supportsSecondGradient() const override { return this->interpolation.giveInterpolationOrder() >= 2; }

protected:
    int giveNsd() override { return 1; }
};

} // end namespace oofem
#endif // igasurfaceelements_h
