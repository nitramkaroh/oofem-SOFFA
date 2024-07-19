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

#ifndef stiffnessexportmodule_h_
#define stiffnessexportmodule_h_

#include <vector>

#include "exportmodule.h"

///@name Input fields for MatlabExportModule
//@{
#define _IFT_StiffnessExportModule_Name "stiffnessexportmodule"
//@}

namespace oofem {
/**
 * 
 */
class OOFEM_EXPORT StiffnessExportModule : public ExportModule
{
protected:
    /// list of InternalStateType values, identifying the selected vars for export


public:
    StiffnessExportModule(int n, EngngModel * e);
    virtual ~StiffnessExportModule();

    //void initializeFrom(InputRecord &ir) override;
    void doOutput(TimeStep *tStep, bool forcedOutput = false) override;
    std::string giveOutputFileName(TimeStep *tStep);
    const char *giveClassName() const override { return "StiffnessExportModule"; }
    const char *giveInputRecordName() const { return _IFT_StiffnessExportModule_Name; }
};
} // end namespace oofem
#endif // stiffnessexportmodule_h_
