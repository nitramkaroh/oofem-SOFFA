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

#include "sm/Elements/3D/lspacefb.h"
#include "fei3dhexalin.h"
#include "gausspoint.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(LSpaceFB);

  LSpaceFB :: LSpaceFB(int n, Domain *aDomain) : LSpace(n, aDomain), FbarElementExtensionInterface()
{ }

void
LSpaceFB :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)

{
  FbarElementExtensionInterface :: computeFbarStiffnessMatrix_3d(answer, rMode, tStep, this);
}



void
LSpaceFB ::computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *ts)
{
  FbarElementExtensionInterface :: computeFirstPKStressVector(answer, gp, ts, this);

}

  

  

} // end namespace oofem
