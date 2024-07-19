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

#ifndef elementmatlabexportmodule_h_
#define elementmatlabexportmodule_h_

#include <vector>

#include "exportmodule.h"

///@name Input fields for MatlabExportModule
//@{
#define _IFT_ElementMatlabExportModule_Name "elementmatlab"
//@}

namespace oofem {
/**
 * The Element Matlab export module enables oofem to export the results on an element to a file in a matlab format. 
 * The format is given by an element, which has to inherit from the appropriate element matlab export modue interface

 */
class OOFEM_EXPORT ElementMatlabExportModule : public ExportModule
{
 private:
  
  std::string giveOutputString(TimeStep *tStep);

public:
  ElementMatlabExportModule(int n, EngngModel * e);
  virtual ~ElementMatlabExportModule();

  //void initializeFrom(InputRecord &ir) override;
    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;
    const char *giveClassName() const override { return "ElementMatlabExportModule"; }
    const char *giveInputRecordName() const { return _IFT_ElementMatlabExportModule_Name; }
};
} // end namespace oofem
#endif // matlabexportmodule_h_
