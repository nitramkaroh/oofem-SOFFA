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

#include "Elements/MixedPressure/AxiSymmetric/l4axisymmp0.h"
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


namespace oofem {
REGISTER_Element( L4AxisymmP0 );


L4AxisymmP0 ::L4AxisymmP0( int n, Domain *aDomain ) :
    L4Axisymm( n, aDomain ), BaseMixedPressureElement()
{
    displacementDofsOrdering = {
        1, 2, 3, 4, 5, 6, 7, 8
    };
    pressureDofsOrdering = {
        9
    };
    this->internalDof = std::make_unique<ElementDofManager>( 1, aDomain, this );
    this->internalDof->appendDof( new MasterDof( this->internalDof.get(), P_f ) );
}

void L4AxisymmP0 ::computePressureNMatrixAt( GaussPoint *gp, FloatArray &answer )
{
    answer.resize( 1 );
    answer.at( 1 ) = 1.0;
}


void L4AxisymmP0 ::giveDofManDofIDMask( int inode, IntArray &answer ) const
{
    answer = {
        D_u, D_v
    };
}


void L4AxisymmP0 ::giveDofManDofIDMask_u( IntArray &answer )
{
    answer = {
        D_u, D_v
    };
}


void L4AxisymmP0 ::giveDofManDofIDMask_p( IntArray &answer )
{
    answer = {
        P_f
    };
}


DofManager *
L4AxisymmP0 ::giveInternalDofManager( int i ) const
{
    return this->internalDof.get();
}

void L4AxisymmP0 ::giveInternalDofManDofIDMask( int i, IntArray &answer ) const
{
    answer = { P_f };
}


void L4AxisymmP0 ::postInitialize()
{
    BaseMixedPressureElement ::postInitialize();
    L4Axisymm ::postInitialize();
}
} // end namespace oofem
