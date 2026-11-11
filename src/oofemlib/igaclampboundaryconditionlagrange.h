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

#ifndef igaclampboundaryconditionlagrange_h
#define igaclampboundaryconditionlagrange_h

#include "activebc.h"
#include "iga/iga.h"
#include "boundarycondition.h"


///@name Input fields for IGAClampBoundaryConditionLagrange
//@{
#define _IFT_IGAClampBoundaryConditionLagrange_Name "igaclampboundaryconditionlagrange"
#define _IFT_IGAClampBoundaryConditionLagrange_ElementNumber "numel" ///< Number of the clamped IGA curve element.
#define _IFT_IGAClampBoundaryConditionLagrange_LocalPosition "locpos" ///< Curve parameter xi of the clamp.
#define _IFT_IGAClampBoundaryConditionLagrange_BcKnotSpan "bcknotspan" ///< Zero-based knot index of the span containing locpos.
//@}

namespace oofem {
/**
 * The clamp n . du/dxi (xi_c) = 0 of IGAClampBoundaryCondition (same input except the penalty, same
 * location check) enforced exactly by a Lagrange multiplier lambda, an additional dof of an internal
 * dof manager of this boundary condition. With c = dN/dxi projected to the reference normal n (over the
 * u, v dofs of the control points of the selected knot span) it adds the blocks c and c^T coupling the
 * span dofs and lambda to the stiffness matrix, and lambda c and c . u to the internal forces.
 */
class OOFEM_EXPORT IGAClampBoundaryConditionLagrange : public ActiveBoundaryCondition
{
protected:
    int elementNumber; ///< IGA curve element on which the condition is imposed.
    FloatArray localPosition; ///< One curve parameter xi, not a physical coordinate.
    IntArray BcKnotSpan; ///< Zero-based knot-vector index, not the ordinal of the nonzero interval.
    std ::unique_ptr<DofManager> lmdm; ///< Internal dof manager with the Lagrange multiplier.

public:
    /// Constructor.
    IGAClampBoundaryConditionLagrange( int n, Domain *d );
    /// Destructor.
    virtual ~IGAClampBoundaryConditionLagrange() {};

    void initializeFrom( InputRecord &ir ) override;
    void postInitialize() override;

    /// Assembles the multiplier blocks for the tangent, secant and elastic stiffness; other types are ignored.
    void assemble( SparseMtrx &answer, TimeStep *tStep,
        CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale = 1.0,
        void *lock = nullptr ) override;

    /// Assembles lambda c to the span dofs and the constraint value c . u to the multiplier (InternalForcesVector only).
    void assembleVector( FloatArray &answer, TimeStep *tStep,
        CharType type, ValueModeType mode,
        const UnknownNumberingScheme &s, FloatArray *eNorms = NULL,
        void *lock = nullptr ) override;


    const char *giveClassName() const override { return "IGAClampBoundaryConditionLagrange"; }
    const char *giveInputRecordName() const override { return _IFT_IGAClampBoundaryConditionLagrange_Name; }

    /// Same as IGAClampBoundaryCondition::giveLocalCodeNumbers.
    int giveLocalCodeNumbers( IntArray &answer, Element *elem, IntArray &span );

    int giveNumberOfInternalDofManagers() override { return 1; }
    DofManager *giveInternalDofManager( int i ) override { return this->lmdm.get(); }

    /// Location arrays of the coupling blocks (all dofs of the element and the multiplier) and of the multiplier diagonal.
    void giveLocationArrays( std ::vector<IntArray> &rows, std ::vector<IntArray> &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s ) override;
    /// Equation number of the multiplier.
    void giveLagrangianMultiplierLocationArray( const UnknownNumberingScheme &r_s, IntArray &answer );
};
} // end namespace oofem
#endif // igaclampboundaryconditionlagrange_h
