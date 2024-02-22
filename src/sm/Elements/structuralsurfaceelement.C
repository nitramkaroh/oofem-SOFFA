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

#include "sm/Elements/structuralsurfaceelement.h"
#include "sm/Elements/structuralelement.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"
#include "sm/Materials/LatticeMaterials/latticematstatus.h"
#include "Loads/structtemperatureload.h"
#include "sm/Materials/structuralnonlocalmaterialext.h"
#include "Loads/structeigenstrainload.h"
#include "feinterpol.h"
#include "domain.h"
#include "material.h"
#include "nonlocalmaterialext.h"
#include "load.h"
#include "boundaryload.h"
#include "pointload.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "nonlocmatstiffinterface.h"
#include "mathfem.h"
#include "materialmapperinterface.h"
#include "unknownnumberingscheme.h"
#include "set.h"

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "connectivitytable.h"
#endif


namespace oofem {
StructuralSurfaceElement ::StructuralSurfaceElement( int n, Domain *aDomain ) :
    StructuralElement( n, aDomain )
{
}


StructuralSurfaceElement ::~StructuralSurfaceElement()
{
}


} // end namespace oofem