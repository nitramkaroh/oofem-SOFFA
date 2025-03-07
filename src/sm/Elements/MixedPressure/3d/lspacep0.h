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
 *               Copyright (C) 1993 - 2015   Borek Patzak
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

#ifndef lspacep0_h
#define lspacep0_h

#include "../sm/Elements/3D/lspace.h"
#include "../sm/Elements/MixedPressure/basemixedpressureelement.h"


#define _IFT_LSpaceP0_Name "lspacep0"

namespace oofem {
class ElementDofManager;

  
class LSpaceP0 : public LSpace, public BaseMixedPressureElement
{
protected:
    /// The extra dofs from the bubble
    std :: unique_ptr< ElementDofManager > internalDof;

public:
    LSpaceP0(int n, Domain *d);
    virtual ~LSpaceP0() { }

protected:
    void computePressureNMatrixAt(GaussPoint *gp, FloatArray &Np) override;
    NLStructuralElement *giveElement() override { return this; }

public:
    const char *giveInputRecordName() const override { return _IFT_LSpaceP0_Name; }
    const char *giveClassName() const override { return "LSpaceP0"; }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    void giveDofManDofIDMask_u(IntArray &answer) override;
    void giveDofManDofIDMask_p(IntArray &answer) override;

    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep) override { BaseMixedPressureElement :: computeStiffnessMatrix(answer, mode, tStep); }
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override { BaseMixedPressureElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }
  ///////////
  void giveInternalDofManDofIDMask(int inode, IntArray &answer) const override;
  int giveNumberOfInternalDofManagers() const override { return 1; }
  DofManager *giveInternalDofManager(int i) const override;
  //////
  int giveNumberOfPressureDofs() override { return 1; }
  int giveNumberOfDisplacementDofs() override { return 24; }
  int giveNumberOfDofs() override { return 25; }
  void postInitialize() override;
};
} // end namespace oofem
#endif // quad1planestrainp1_h
