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

#pragma once

#include "interface.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "element.h"
#include "feinterpol.h"
#include "gausspoint.h"


namespace oofem {
class GaussPoint;
class IntegrationRule;
class TimeStep;
/**
 * Provides Pressure Follower Load for an element.
 * @author Martin Horak
 */
class OOFEM_EXPORT MagneticLoadElementInterface : public Interface
{
public:
    Element *element;


    /// Constructor.
    MagneticLoadElementInterface(Element *e);
    virtual ~MagneticLoadElementInterface();

    virtual const char *giveClassName() const override { return "MagneticLoadElementInterface"; }


    
    // private:
    virtual void surfaceEvalNmatrixAt( FloatMatrix &answer, int iSurf, GaussPoint *gp ) = 0;
    virtual void surfaceEvalBmatrixAt( FloatMatrix &answer, int iSurf, GaussPoint *gp ) = 0;
    virtual void surfaceEvalNormalAt(FloatArray &answer, int iSurf, GaussPoint *gp, TimeStep *tStep){;}
    virtual void surfaceEvalDeformationGradientAt( FloatArray &answer, int isurf, GaussPoint *gp, TimeStep *tStep ) { ; }
    
    virtual IntegrationRule* surfaceGiveIntegrationRule(int order, int iSurf) = 0;

};
} // end namespace oofem
