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

#ifndef qtrspacep1_h
#define qtrspacep1_h

#include "../sm/Elements/3D/qtrspace.h"
#include "../sm/Elements/MixedPressure/basemixedpressureelement.h"


#define _IFT_QTRSpaceP1_Name "qtrspacep1"

namespace oofem {
class ElementDofManager;
class FEI3dTetLin;
  
class QTRSpaceP1 : public QTRSpace, public BaseMixedPressureElement
{
protected:
    static FEI3dTetLin interpolation_lin;

public:
    QTRSpaceP1(int n, Domain *d);
    virtual ~QTRSpaceP1() {}

protected:
    void computePressureNMatrixAt(GaussPoint *gp, FloatArray &Np) override;
    NLStructuralElement *giveElement() override { return this; }

public:
    const char *giveInputRecordName() const override { return _IFT_QTRSpaceP1_Name; }
    const char *giveClassName() const override { return "QTRSpaceP1"; }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    void giveDofManDofIDMask_u(IntArray &answer) override;
    void giveDofManDofIDMask_p(IntArray &answer) override;

    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep) override { BaseMixedPressureElement :: computeStiffnessMatrix(answer, mode, tStep); }
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override { BaseMixedPressureElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord); }

  int giveNumberOfPressureDofs() override { return 4; }
  int giveNumberOfDisplacementDofs() override { return 30; }
  int giveNumberOfDofs() override { return 34; }
  void postInitialize() override;
};
} // end namespace oofem
#endif // ltrspacep0_h
