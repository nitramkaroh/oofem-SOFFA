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

#ifndef space2dstructuraleurfaceelementevaluator_h
#define space2dstructuraleurfaceelementevaluator_h

#include "sm/Elements/structuralsurfaceelementevaluator.h"

namespace oofem {

// Plane strain formulation for surface 
class Space2dStructuralSurfaceElementEvaluator : public StructuralSurfaceElementEvaluator
{
public:
    Space2dStructuralSurfaceElementEvaluator() :
        StructuralSurfaceElementEvaluator() {}

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
    //void computeBMatrixAt( FloatMatrix &answer, GaussPoint *gp ) override;
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    virtual void computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    double computeVolumeAround( GaussPoint *gp ) override;
    //void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;
    //void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;
    void giveDofManDofIDMask( int inode, IntArray &answer ) const
    {
        answer = { D_u, D_v };
    }
    //void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

public:
    void initializeFrom( InputRecord &ir ) override;

    int giveIntegrationElementLocalCodeNumbers( IntArray &answer, Element *elem, IntegrationRule *ie ) override;

    //int giveIntegrationElementLocalCodeNumbers( IntArray &answer, Element *elem, FloatArray &lcoords );


}; // end 

////////////////
// Axisymmetrix formulation
class AxisymSurfaceElementEvaluator : public Space2dStructuralSurfaceElementEvaluator
{
public:
    AxisymSurfaceElementEvaluator() :
        Space2dStructuralSurfaceElementEvaluator() {}

protected:

    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    virtual void computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    double computeVolumeAround( GaussPoint *gp ) override;

public:
    void initializeFrom( InputRecord &ir ) override;
}; // end 


} // end namespace oofem
#endif 
