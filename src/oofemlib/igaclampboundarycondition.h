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

#ifndef igaclampboundarycondition_h
#define igaclampboundarycondition_h

#include "activebc.h"
#include "iga/iga.h"
#include "boundarycondition.h"


///@name Input fields for IGAClampBoundaryCondition
//@{
#define _IFT_IGAClampBoundaryCondition_Name "igaclampboundarycondition"
#define _IFT_IGAClampBoundaryCondition_ElementNumber "numel" ///< Number of the clamped IGA curve element.
#define _IFT_IGAClampBoundaryCondition_LocalPosition "locpos" ///< Curve parameter xi of the clamp.
#define _IFT_IGAClampBoundaryCondition_BcKnotSpan "bcknotspan" ///< Zero-based knot index of the span containing locpos.
#define _IFT_IGAClampBoundaryCondition_p "p" ///< Penalty factor (optional, default 1e6).
//@}

namespace oofem {
/**
 * Clamp of a planar IGA curve (BSpline2dLineInterpolation, i.e. the IGA surface elements) at one value
 * xi_c of the curve parameter, enforced by a penalty: n . du/dxi (xi_c) = 0 with the normal n of the
 * reference curve, so the tangent at xi_c may stretch but not rotate. The condition is linear in the
 * displacements: with c = dN/dxi projected to n (over the u, v dofs of the control points of the
 * selected span) the matrix p c c^T is the same for all stiffness types and the internal forces are
 * p c (c . u).
 *
 * The normal and the basis derivatives are evaluated in the explicitly given knot span, so a clamp at
 * an internal knot uses the side chosen in the input. The location is checked in postInitialize.
 *
 * Input: numel (element number), locpos (curve parameter, not a physical coordinate), bcknotspan
 * (zero-based index into the full knot vector including the repeated end knots; the span
 * [knot(k), knot(k+1)] must have nonzero length), optional p (penalty, default 1e6).
 * @see IGAClampBoundaryConditionLagrange for the same condition with a Lagrange multiplier.
 */
class OOFEM_EXPORT IGAClampBoundaryCondition : public ActiveBoundaryCondition
{
protected:
    int elementNumber; ///< IGA curve element on which the condition is imposed.
    FloatArray localPosition; ///< One curve parameter xi, not a physical coordinate.
    IntArray BcKnotSpan; ///< Zero-based knot-vector index, not the ordinal of the nonzero interval.
    double penalty = 1e6; ///< Penalty factor.

public:
    /// Constructor.
    IGAClampBoundaryCondition( int n, Domain *d );
    /// Destructor.
    virtual ~IGAClampBoundaryCondition() {};

    void initializeFrom( InputRecord &ir ) override;
    void postInitialize() override;

    /**
     * Checks the clamp location; shared by IGAClampBoundaryCondition and IGAClampBoundaryConditionLagrange.
     * Element elementNumber must be a planar IGA curve element of degree at least 1, and knotSpan
     * (a zero-based knot index) must be a nonzero knot interval containing localPosition.
     * Both endpoints are included, so an internal knot may be evaluated from either selected span.
     */
    static void checkClampLocation( Domain *d, int elementNumber, double localPosition, int knotSpan );

    /// Assembles the penalty matrix for the tangent, secant and elastic stiffness; other types are ignored.
    void assemble( SparseMtrx &answer, TimeStep *tStep,
        CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0,
        void *lock = nullptr ) override;

    /// Assembles the penalty internal forces (InternalForcesVector only).
    void assembleVector( FloatArray &answer, TimeStep *tStep,
        CharType type, ValueModeType mode,
        const UnknownNumberingScheme &s, FloatArray *eNorms = NULL,
        void *lock = nullptr ) override;


    const char *giveClassName() const override { return "IGAClampBoundaryCondition"; }
    const char *giveInputRecordName() const override { return _IFT_IGAClampBoundaryCondition_Name; }

    /**
     * Gives the element-local code numbers of the u, v dofs of the control points whose basis functions
     * are nonzero in the knot span (two per control point; further dofs of a node give 0).
     * @return 1 if the interpolation has a sub-patch (knot span) formulation, 0 otherwise (answer untouched).
     */
    int giveLocalCodeNumbers( IntArray &answer, Element *elem, IntArray &span );
};
} // end namespace oofem
#endif // igaclampboundarycondition_h
