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

#ifndef qspacefb_h
#define qspacefb_h

#include "sm/Elements/3D/qspace.h"
#include "sm/Elements/Interfaces/fbarelementinterface.h"

#define _IFT_QSpaceFB_Name "qspacefb"

namespace oofem {
/**
 * Three dimensional tetrahedron with linear approximation, suitable for incompressible settings in large strain setting
 * This is achieved by F-bar formulation
 */
class QSpaceFB : public QSpace, public FbarElementExtensionInterface
{
public:
    QSpaceFB( int n, Domain *d );
    virtual ~QSpaceFB() {}

    const char *giveInputRecordName() const override { return _IFT_QSpaceFB_Name; }
    const char *giveClassName() const override { return "QSpaceFB"; }
    void giveElementParametricCentroid( FloatArray &centroidCoords ) override { centroidCoords = { 0., 0., 0. }; }
    void computeStiffnessMatrix( FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep ) override;
    void computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *ts ) override;
};
} // end namespace oofem
#endif // lspacefb_h
