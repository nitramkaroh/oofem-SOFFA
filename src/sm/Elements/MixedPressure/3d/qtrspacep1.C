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

#include "Elements/MixedPressure/3d/qtrspacep1.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "classfactory.h"
#include "elementinternaldofman.h"
#include "masterdof.h"
#include "fei3dtetlin.h"


namespace oofem {
REGISTER_Element(QTRSpaceP1);


FEI3dTetLin QTRSpaceP1 :: interpolation_lin;
  
  
QTRSpaceP1 :: QTRSpaceP1(int n, Domain *aDomain) : QTRSpace(n, aDomain), BaseMixedPressureElement()
{
    displacementDofsOrdering = {
        1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25,26, 27, 28, 29,30,31,32,33,34
    };
    pressureDofsOrdering = {
      4,8,12,16
    };

}

void QTRSpaceP1 ::computePressureNMatrixAt( GaussPoint *gp, FloatArray &answer )
{
    this->interpolation_lin.evalN( answer, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( this ) );
}



void
QTRSpaceP1 :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode <= 4 ) {
        answer = {
            D_u, D_v, D_w,P_f
        };
    } else {
        answer = {
            D_u, D_v, D_w
        };
    }
}


void
QTRSpaceP1 :: giveDofManDofIDMask_u(IntArray &answer)
{
    answer = {
      D_u, D_v, D_w
    };
}


void
QTRSpaceP1 :: giveDofManDofIDMask_p(IntArray &answer)
{
    answer = {
        P_f
    };
}
  
  
void
QTRSpaceP1 ::  postInitialize()
{
    BaseMixedPressureElement :: postInitialize();
    QTRSpace :: postInitialize();
}
} // end namespace oofem
