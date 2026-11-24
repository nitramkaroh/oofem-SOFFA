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

#include "sm/Elements/SurfaceElements/space2dsurfaceelementevaluator.h"
#include "sm/Elements/SurfaceElements/space3dsurfaceelementevaluator.h"

#include "floatarray.h"
#include "floatmatrix.h"
#include "matresponsemode.h"
#include "mathfem.h"

#define _IFT_BsplinePlaneStrainSurfaceElement_Name "bsplineplanestrainsurfaceelement"
#define _IFT_BsplineAxisymSurfaceElement_Name "bsplineaxisymsurfaceelement"
#define _IFT_Bspline3dSurfaceElement_Name "bspline3dsurfaceelement"
namespace oofem {

/**
 * Surface energy on a planar B-spline curve (one parametric direction, control points in the xy plane),
 * typically the boundary of a plane strain B-spline patch sharing its control points. Per unit length
 * of the curve and unit thickness. Large strains only (nlgeo 1); the cross section is a
 * SimpleSurfaceCrossSection with a surface material.
 *
 * One element is the whole curve and its integration rules are the nonzero knot spans. The dofs are u, v
 * of the control points; the kinematics, stresses and stiffness come from
 * Space2dStructuralSurfaceElementEvaluator. The reference normal (giveNormal) points to the right of the
 * direction of increasing xi, i.e. of the node order (outwards for a boundary traversed counterclockwise);
 * it fixes the sign of the curvature of a bending material.
 */
class BsplinePlaneStrainSurfaceElement : public IGAElement, public Space2dStructuralSurfaceElementEvaluator
{
protected:
    BSpline2dLineInterpolation interpolation;

public:
    BsplinePlaneStrainSurfaceElement( int n, Domain *aDomain );

    /// Reads nlgeo (the evaluator accepts only 1) and the curve (IGAElement::initializeFrom).
    void initializeFrom( InputRecord &ir ) override;
    /// Base post-initialization, then setupEvaluator: decides from the material whether the second gradient is used.
    void postInitialize() override { IGAElement :: postInitialize(); this->setupEvaluator(); }
    /// The number of nodes must equal the number of control points and the cross section be a SimpleSurfaceCrossSection.
    int checkConsistency() override;

    // the element matrices and vectors, the dof masks and the state update are those of the evaluator
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

    /// The surface materials work with the full 3D deformation gradient (9 components).
    MaterialMode giveMaterialMode() override { return _3dMat; }

    /// Unit normal of the reference curve, (t_y, -t_x, 0) / |t| with the tangent t = dX/dxi.
    FloatArray giveNormal( GaussPoint *gp ) const override;

    /**
     * IST_CurvatureTensor for the VTK export of IGAElement (bending materials only), evaluated at the export
     * point gp from the displacements (ir is not needed).
     */
    bool giveExportPointValue( FloatArray &answer, InternalStateType type, GaussPoint *gp, int ir, TimeStep *tStep ) override
    {
        return this->giveCurvatureExportValue( answer, type, gp, ir, tStep );
    }

    /// The second gradient (bending) needs C1 continuity: degree >= 2 and interior knot multiplicities <= degree - 1.
    bool supportsSecondGradient() const override
    {
        return this->interpolation.giveInterpolationOrder() >= 2 && this->interpolation.giveMinimalContinuity( 1 ) >= 1;
    }

protected:
    int giveNsd() override { return 1; }
};


/**
 * Surface energy on the surface of revolution of a B-spline curve in the xy plane (x the radius,
 * y the axis of rotation), per radian. Large strains only (nlgeo 1); the cross section is a
 * SimpleSurfaceCrossSection with a surface material.
 *
 * One element is the whole curve and its integration rules are the nonzero knot spans. The dofs are u
 * (radial) and v (axial) of the control points; the kinematics, stresses and stiffness come from
 * AxisymSurfaceElementEvaluator, which needs a positive reference radius at the Gauss points. The normal
 * is oriented as in BsplinePlaneStrainSurfaceElement.
 */
class BsplineAxisymSurfaceElement : public IGAElement, public AxisymSurfaceElementEvaluator
{
protected:
    BSpline2dLineInterpolation interpolation;

public:
    BsplineAxisymSurfaceElement( int n, Domain *aDomain );

    /// Reads nlgeo (the evaluator accepts only 1) and the curve (IGAElement::initializeFrom).
    void initializeFrom( InputRecord &ir ) override;
    /// Base post-initialization, then setupEvaluator: decides from the material whether the second gradient is used.
    void postInitialize() override { IGAElement :: postInitialize(); this->setupEvaluator(); }
    /// The number of nodes must equal the number of control points and the cross section be a SimpleSurfaceCrossSection.
    int checkConsistency() override;

    // the element matrices and vectors, the dof masks and the state update are those of the evaluator
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

    /// The surface materials work with the full 3D deformation gradient (9 components, z the hoop direction).
    MaterialMode giveMaterialMode() override { return _3dMat; }

    /// Unit normal of the reference curve, (t_y, -t_x, 0) / |t| with the tangent t = dX/dxi.
    FloatArray giveNormal( GaussPoint *gp ) const override;

    /**
     * IST_CurvatureTensor for the VTK export of IGAElement (bending materials only), evaluated at the export
     * point gp from the displacements (ir is not needed).
     */
    bool giveExportPointValue( FloatArray &answer, InternalStateType type, GaussPoint *gp, int ir, TimeStep *tStep ) override
    {
        return this->giveCurvatureExportValue( answer, type, gp, ir, tStep );
    }

    /// The second gradient (bending) needs C1 continuity: degree >= 2 and interior knot multiplicities <= degree - 1.
    bool supportsSecondGradient() const override
    {
        return this->interpolation.giveInterpolationOrder() >= 2 && this->interpolation.giveMinimalContinuity( 1 ) >= 1;
    }

protected:
    int giveNsd() override { return 1; }
};


/**
 * Surface energy on a B-spline surface in space (two parametric directions, control points with three
 * coordinates and the dofs u, v, w), per unit area of the reference surface. Bending materials (second
 * gradient) need C1 continuity in both directions. Large strains only (nlgeo 1); the cross section is a
 * SimpleSurfaceCrossSection with a surface material.
 *
 * One element is the whole patch (nodes = control points, u index fastest) and its integration rules are
 * the nonzero knot spans; the kinematics, stresses and stiffness come from
 * Space3dStructuralSurfaceElementEvaluator. The reference normal N = A_1 x A_2 / |A_1 x A_2| with
 * A_1 = X_,xi and A_2 = X_,eta (right-hand rule in xi, eta) fixes the sign of the curvature of a bending
 * material.
 */
class Bspline3dSurfaceElement : public IGAElement, public Space3dStructuralSurfaceElementEvaluator
{
protected:
    BSpline3dSurfaceInterpolation interpolation;

public:
    Bspline3dSurfaceElement( int n, Domain *aDomain );

    /// Reads nlgeo (the evaluator accepts only 1) and the patch (IGAElement::initializeFrom).
    void initializeFrom( InputRecord &ir ) override;
    /// Base post-initialization, then setupEvaluator: decides from the material whether the second gradient is used.
    void postInitialize() override { IGAElement :: postInitialize(); this->setupEvaluator(); }
    /// The number of nodes and their coordinates (three each) and the cross section (SimpleSurfaceCrossSection).
    int checkConsistency() override;

    // the element matrices and vectors, the dof masks and the state update are those of the evaluator
    void giveCharacteristicMatrix( FloatMatrix &answer, CharType mtrx, TimeStep *tStep ) override
    {
        Space3dStructuralSurfaceElementEvaluator ::giveCharacteristicMatrix( answer, mtrx, tStep );
    }
    void giveCharacteristicVector( FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t ) override
    {
        Space3dStructuralSurfaceElementEvaluator ::giveCharacteristicVector( answer, type, mode, t );
    }

    FEInterpolation *giveInterpolation() const override { return const_cast<BSpline3dSurfaceInterpolation *>( &this->interpolation ); }
    Element *giveElement() override { return this; }
    void giveDofManDofIDMask( int inode, IntArray &answer ) const override
    {
        Space3dStructuralSurfaceElementEvaluator ::giveDofManDofIDMask( inode, answer );
    }
    int computeNumberOfDofs() override { return numberOfDofMans * 3; }
    void updateInternalState( TimeStep *tStep ) override { Space3dStructuralSurfaceElementEvaluator ::updateInternalState( tStep ); }
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_Bspline3dSurfaceElement_Name; }
    const char *giveClassName() const override { return "Bspline3dSurfaceElement"; }

    /// The surface materials work with the full 3D deformation gradient (9 components).
    MaterialMode giveMaterialMode() override { return _3dMat; }

    /// Unit normal of the reference surface, X_,xi x X_,eta / |X_,xi x X_,eta|.
    FloatArray giveNormal( GaussPoint *gp ) const override;

    /**
     * IST_CurvatureTensor for the VTK export of IGAElement (bending materials only), evaluated at the export
     * point gp from the displacements (ir is not needed).
     */
    bool giveExportPointValue( FloatArray &answer, InternalStateType type, GaussPoint *gp, int ir, TimeStep *tStep ) override
    {
        return this->giveCurvatureExportValue( answer, type, gp, tStep );
    }

    /**
     * The second gradient (bending) needs C1 continuity in both directions: degree >= 2 and interior knot
     * multiplicities <= degree - 1.
     */
    bool supportsSecondGradient() const override
    {
        return this->interpolation.giveInterpolationOrder() >= 2 && this->interpolation.giveMinimalContinuity( 1 ) >= 1 &&
               this->interpolation.giveMinimalContinuity( 2 ) >= 1;
    }

protected:
    int giveNsd() override { return 2; }
};

} // end namespace oofem
#endif // igasurfaceelements_h
