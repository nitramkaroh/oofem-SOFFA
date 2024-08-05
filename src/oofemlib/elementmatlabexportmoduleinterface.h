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

#ifndef elementmatlabexportmoduleinterface_h
#define elementmatlabexportmoduleinterface_h

#include "interface.h"

namespace oofem {
/**
 * Element Matlab Export Module Interface allowing to export each element separately its variables to matlab. 
 * Currently used for beam elements
 */
class ElementMatlabExportModuleInterface : public Interface
{
protected:


public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    ElementMatlabExportModuleInterface(){;}
    /// Destructor.
    virtual ~ElementMatlabExportModuleInterface() { }
  //virtual void initializeFrom(InputRecord &ir) = 0;
    virtual void doMatlabOutput(std::string s) = 0;
  const char *giveClassName() const override {return "ElementMatlabExportModulInterface";}

};

}
#endif
