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

#include "mpm/MagnetoElasticity/Materials/magnetoelasticmaterialstatus.h"
#include "contextioerr.h"
#include "gausspoint.h"

namespace oofem {
MagnetoElasticMaterialStatus :: MagnetoElasticMaterialStatus(GaussPoint *g) :
  MaterialStatus(g), FVector(), tempFVector(), HVector(), tempHVector(), gradFVector(), tempGradFVector(), FbarVector(), tempFbarVector()
{
    PVector.resize(9);
    FVector.resize(9);
    FVector.at(1) = FVector.at(2) = FVector.at(3) = 1.;
    gradFVector.resize(27);
    FbarVector.resize( 9 );
    FbarVector.at( 1 ) = FbarVector.at( 2 ) = FbarVector.at( 3 ) = 1.;
    tempPVector = PVector;
    tempFVector = FVector;
    tempGradFVector = gradFVector;
    tempFbarVector = FbarVector;
    
    HVector.resize(3);
    BVector.resize(3);
    tempHVector = HVector;
    tempBVector = BVector;
    

    

}

void MagnetoElasticMaterialStatus :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    MaterialStatus :: updateYourself(tStep);

    PVector = tempPVector;
    FVector = tempFVector;
    gradFVector = tempGradFVector;
    FbarVector = tempFbarVector;
    HVector = tempHVector;
    BVector = tempBVector;
}


void MagnetoElasticMaterialStatus :: initTempStatus()
//
// initialize record at the begining of new load step
//
{
    MaterialStatus :: initTempStatus();

    tempPVector = PVector;
    tempFVector = FVector;
    tempGradFVector = gradFVector;
    tempFbarVector = FbarVector;

    tempHVector = HVector;
    tempBVector = BVector;

}

  /*
void
MagnetoElasticMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    MaterialStatus :: saveContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = strainVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stressVector.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}


void
MagnetoElasticMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    MaterialStatus :: restoreContext(stream, mode);

    contextIOResultType iores;
    if ( ( iores = strainVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = stressVector.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }
}
  */

} // end namespace oofem
