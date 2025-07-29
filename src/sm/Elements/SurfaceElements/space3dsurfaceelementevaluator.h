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

#ifndef space3dsurfaceelementevaluator_h
#define space3dsurfaceelementevaluator_h

#include "sm/Elements/structuralsurfaceelementevaluator.h"

namespace oofem {
/**
 * General purpose 3d structural element evaluator.
 */
class Space3dStructuralSurfaceElementEvaluator : public StructuralSurfaceElementEvaluator
{
public:
    Space3dStructuralSurfaceElementEvaluator() :
        StructuralSurfaceElementEvaluator() {}
    virtual ~Space3dStructuralSurfaceElementEvaluator() {}

protected:
    /**
     * Assemble interpolation matrix at given IP
     * In case of IGAElements, N is assumed to contain only nonzero interpolation functions
     */
    void computeNMatrixAt( FloatMatrix &answer, GaussPoint *gp ) override;
    /**
     * Assembles the strain-displacement matrix of the receiver at given integration point
     * In case of IGAElements, B is assumed to contain only contribution from nonzero interpolation functions
     */

    //void computeBMatrixAt( FloatMatrix &answer, GaussPoint *gp ) override;
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;

    double computeVolumeAround( GaussPoint *gp ) override;

    void giveDofManDofIDMask( int inode, IntArray &answer ) const
    {
        answer = { D_u, D_v, D_w };
    }

public:
    void initializeFrom( InputRecord &ir ) override;
};
} // end namespace oofem
#endif // space3delementevaluator_h
