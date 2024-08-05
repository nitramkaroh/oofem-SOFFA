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

#include <algorithm>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <iterator>

#include "elementmatlabexportmodule.h"
#include "elementmatlabexportmoduleinterface.h"
#include "timestep.h"
#include "classfactory.h"
#include "element.h"
#include "engngm.h"


namespace oofem {

REGISTER_ExportModule( ElementMatlabExportModule )

ElementMatlabExportModule :: ElementMatlabExportModule(int n, EngngModel *e) : ExportModule(n, e)
{}


ElementMatlabExportModule :: ~ElementMatlabExportModule()
{ }


void
ElementMatlabExportModule :: doOutput(TimeStep *tStep, bool forcedOutput)
{
    if ( !( testTimeStepOutput(tStep) || forcedOutput ) ) {
        return;
    }


    //int nelem = this->elList.giveSize();
    /*    if ( nelem == 0 ) { // no list given, export all elements
        this->elList.enumerate(this->emodel->giveDomain(1)->giveNumberOfElements());
    }
    */
    Domain *domain  = emodel->giveDomain(1);
    for ( int ielem = 1; ielem <= domain->giveNumberOfElements(); ielem++ ) {
      //Element *el = domain->giveElement( this->elList.at(ielem) );
      Element *el = domain->giveElement( ielem );
      ElementMatlabExportModuleInterface* ei = dynamic_cast<ElementMatlabExportModuleInterface * >(el);
      if(!ei) {
	OOFEM_WARNING("Elemend %d is not shooting method based beam element, no matlab output generated", el->giveNumber());
      } else {
	ei->doMatlabOutput(this->giveOutputString(tStep));
      }

    }
}





std::string 
ElementMatlabExportModule :: giveOutputString(TimeStep *tStep)
{
    char fext[100];
    sprintf( fext, "_m%d_%d", this->number, tStep->giveNumber() );

    if ( this->testSubStepOutput() ) {
        // include tStep version in output file name
        if ( this->emodel->isParallel() && this->emodel->giveNumberOfProcesses() > 1 ) {
            sprintf( fext, "_%03d_m%d_%d_%d", emodel->giveRank(), this->number, tStep->giveNumber(), tStep->giveSubStepNumber() );
        } else {
            sprintf( fext, "_m%d_%d_%d", this->number, tStep->giveNumber(), tStep->giveSubStepNumber() );
        }
    } else   {
        if ( this->emodel->isParallel() && this->emodel->giveNumberOfProcesses() > 1 ) {
            sprintf( fext, "_%03d_m%d_%d", emodel->giveRank(), this->number, tStep->giveNumber() );
        } else {
            sprintf( fext, "_m%d_%d", this->number, tStep->giveNumber() );
        }
    }

    std :: string fileName;

    fileName = this->emodel->giveOutputBaseFileName();

    size_t foundDot;
    foundDot = fileName.rfind(".");

    // CARL
    while (foundDot != std :: string :: npos) {
        fileName.replace(foundDot, 1, "_");
        foundDot = fileName.rfind(".");
    }

    fileName += fext;
    return fileName;
}


} // end namespace oofem
