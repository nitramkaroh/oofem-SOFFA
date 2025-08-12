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
#pragma once

#include "element.h"
#include "dofmanager.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "intarray.h"
#include "classfactory.h"


namespace oofem {

/*
 * Note: someone should be able to return for given cell and variable vector of unknowns.
 * this depends on interpolation (constant, linear, etc), cell type and variable (defined the physical meaning of unknown(s))
 * Interpolation should identify (or even introduce) cell nodes needed (quadratic element, linear interpolation), variable should assign to these nodes DOFs.
 * interpolation.getCellNodes(cell)
 */

/**
 * @brief Class representing unknown field (or test feld) in a weak psolution.
 * The variable has its interpolation, type (scalar, vector), size.
 * When test (dual) field, it keeps reference to its primary (unknown) variable.
 * @todo The history parameter dermines how many time steps to remember.
 */
class Variable
{
  public:
  enum VariableType {
    scalar,
    vector
  };

  enum VariableQuantity {
    Displacement,
    Velocity,
    Temperature,
    Pressure,
    VolumeFraction,
    MagneticPotential
  };

  const FEInterpolation &interpolation;
  Variable *dualVar; //? or just bool?
  VariableType type;
  VariableQuantity q;
  int size;
  IntArray dofIDs;

  Variable( const FEInterpolation &i, Variable::VariableQuantity q, Variable::VariableType t, int size, Variable *dual = NULL, std ::initializer_list<int> dofIDs = {} ) :
      interpolation( i ),
      dualVar( dual ),
      q( q ),
      dofIDs( dofIDs )
  {
    this->type = t;
    this->size = size;
  }
  Variable( const FEInterpolation &i, Variable::VariableQuantity q, Variable::VariableType t, int size, IntArray &dofIDs, Variable *dual = NULL ) :
      interpolation( i ),
      dualVar( dual ),
      q( q ),
      dofIDs( dofIDs )
  {
    this->type = t;
    this->size = size;
  }


  /// Returns DodIF mask in node; need generalization (which dofMan)
  const IntArray &getDofManDofIDs() const { return this->dofIDs; }
};

} //end namespace oofem
