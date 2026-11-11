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

#ifndef planestresselementevaluator_h
#define planestresselementevaluator_h

#include "sm/Elements/structuralelementevaluator.h"
#include "oofemlib/tensor/tensor1.h"
#include "oofemlib/tensor/tensor2.h"
#include "oofemlib/tensor/tensor3.h"

namespace oofem {
/**
 * General purpose Plane stress structural element evaluator.
 */
class PlaneStressStructuralElementEvaluator : public StructuralElementEvaluator
{
public:
    PlaneStressStructuralElementEvaluator() : StructuralElementEvaluator() { }

protected:
    /**
     * Assemble interpolation matrix at given IP.
     * In case of IGAElements, N is assumed to contain only nonzero interpolation functions.
     */
    void computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp) override;
    /**
     * Assembles the strain-displacement matrix of the receiver at given integration point.
     * In case of IGAElements, B is assumed to contain only contribution from nonzero interpolation functions.
     */
    void computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp) override;
    /// Displacement gradient matrix for nlgeo 1, 4 x (2 nno): rows 11, 22, 12, 21.
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    double computeVolumeAround(GaussPoint *gp) override;
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const {
        answer = {D_u, D_v};
    }
    void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

public:
    void initializeFrom( InputRecord &ir ) override;


}; // end of PlaneStressStructuralElementEvaluator definition


/**
 * Plane strain evaluator (dofs u, v), small strains (nlgeo 0) or the first Piola-Kirchhoff formulation
 * (nlgeo 1), the latter also with second gradient materials. The volume includes the cross section thickness.
 */
class PlaneStrainStructuralElementEvaluator : public StructuralElementEvaluator
{
public:
    PlaneStrainStructuralElementEvaluator() :
        StructuralElementEvaluator() {}

protected:
    /**
     * Assemble interpolation matrix at given IP.
     * In case of IGAElements, N is assumed to contain only nonzero interpolation functions.
     */
    void computeNMatrixAt( FloatMatrix &answer, GaussPoint *gp ) override;
    /**
     * Assembles the strain-displacement matrix of the receiver at given integration point.
     * In case of IGAElements, B is assumed to contain only contribution from nonzero interpolation functions.
     */
    void computeBMatrixAt( FloatMatrix &answer, GaussPoint *gp ) override;
    /// Displacement gradient matrix for nlgeo 1, 5 x (2 nno): rows 11, 22, 33 (zero), 12, 21.
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    double computeVolumeAround( GaussPoint *gp ) override;
    void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;
    void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;
    void giveDofManDofIDMask( int inode, IntArray &answer ) const
    {
        answer = { D_u, D_v };
    }
    void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    /// Second gradient matrix, 27 x (2 nno) (full 27 component grad F, see StructuralElementEvaluator).
    void computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /// grad F = G x with the current positions x of the control points. @todo The same in all evaluators.
    void computeGradientOfDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

    /// dA/dF reduced to the plane strain components of F, 27 x 5
    void computeConstitutiveMatrix_dAdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

public:
    void initializeFrom( InputRecord &ir ) override;

}; // end of PlaneStrainStructuralElementEvaluator definition

/**
 * Axisymmetric evaluator (dofs u, v; x is the radius, y the axis, z the hoop direction), first
 * Piola-Kirchhoff formulation only (nlgeo 1), also with second gradient materials. The material is
 * evaluated in the full 3d mode with the hoop stretch F_33 = 1 + u / r; the volume is per radian
 * (r dr dz), the cross section thickness is not used.
 */
class AxisymStructuralElementEvaluator : public StructuralElementEvaluator
{
public:
    AxisymStructuralElementEvaluator() :
        StructuralElementEvaluator() {}

protected:
    /**
     * Assemble interpolation matrix at given IP.
     * In case of IGAElements, N is assumed to contain only nonzero interpolation functions.
     */
    void computeNMatrixAt( FloatMatrix &answer, GaussPoint *gp ) override;
    /// Displacement gradient matrix, 9 x (2 nno), with the hoop row 3 (u / r).
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;
    double computeVolumeAround( GaussPoint *gp ) override;
    
    void giveDofManDofIDMask( int inode, IntArray &answer ) const
    {
        answer = { D_u, D_v };
    }
    void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    // small strains (nlgeo 0) are not implemented
    void computeBMatrixAt( FloatMatrix &answer, GaussPoint *gp ) override
    {
        OOFEM_ERROR( "This element is only for large strains so far" );
        return;
    }
    void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override
    {
        OOFEM_ERROR( "This element is only for large strains so far" );
        return;
    }

    void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override
    {
        OOFEM_ERROR( "This element is only for large strains so far" );
        return;
    }

    /// Second gradient matrix, 27 x (2 nno), including the hoop terms.
    void computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    /// grad F = G x with the current positions x of the control points. @todo The same in all evaluators.
    void computeGradientOfDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

public:
    void initializeFrom( InputRecord &ir ) override;

}; // end of AxisymStructuralElementEvaluator definition



} // end namespace oofem
#endif //planestresselementevaluator_h
