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
  MaterialStatus(g), FVector(), tempFVector(), HVector(), tempHVector()
{
    PVector.resize(9);
    FVector.resize(9);
    FVector.at(1) = FVector.at(2) = FVector.at(3) = 1.;
    tempPVector = PVector;
    tempFVector = FVector;
    tempGradFVector = gradFVector;
    
    HVector.resize(3);
    BVector.resize(3);
    tempHVector = HVector;
    tempBVector = BVector;
    

    

}

  /*
void StructuralMaterialStatus :: printOutputAt(FILE *File, TimeStep *tStep) const
// Prints the strains and stresses on the data file.
{
    FloatArray helpVec;

    MaterialStatus :: printOutputAt(File, tStep);
    NLStructuralElement * el = static_cast< NLStructuralElement * >( gp->giveElement());
    if ( el->giveGeometryMode() == 1) {
      fprintf(File, "  F ");
      StructuralMaterial :: giveFullVectorFormF( helpVec, FVector, gp->giveMaterialMode() );
      for ( auto &var : helpVec ) {
        fprintf( File, " %+.4e", var );
      }

      fprintf(File, "\n  P");
      StructuralMaterial :: giveFullVectorForm( helpVec, PVector, gp->giveMaterialMode() );
      for ( auto &var : helpVec ) {
        fprintf( File, " %+.4e", var );
      }
    } else {
      fprintf(File, "  strains ");
      StructuralMaterial :: giveFullSymVectorForm( helpVec, strainVector, gp->giveMaterialMode() );
      for ( auto &var : helpVec ) {
        fprintf( File, " %+.4e", var );
      }
      
      fprintf(File, "\n              stresses");
      StructuralMaterial :: giveFullSymVectorForm( helpVec, stressVector, gp->giveMaterialMode() );
      
      for ( auto &var : helpVec ) {
        fprintf( File, " %+.4e", var );
      }
    }
    fprintf(File, "\n");
}
  */

void MagnetoElasticMaterialStatus :: updateYourself(TimeStep *tStep)
// Performs end-of-step updates.
{
    MaterialStatus :: updateYourself(tStep);

    PVector = tempPVector;
    FVector = tempFVector;
    gradFVector = tempGradFVector;
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
