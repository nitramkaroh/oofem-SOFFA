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

#ifndef structural3surfacedelement_h
#define structural3dsurfaceelement_h

#include "sm/Elements/nlstructuralsurfaceelement.h"
#include "feinterpol2d.h"


#define _IFT_Structural3DSurfaceElement_materialCoordinateSystem "matcs" ///< [optional] Support for material directions based on element orientation.

namespace oofem {
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Base class 3D surface elements.
 */
class Structural3DSurfaceElement : public NLStructuralSurfaceElement
{
protected:
    bool matRotation;
    FEICellGeometry *cellGeometryWrapper;

public:
    /**
     * Constructor. Creates element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    Structural3DSurfaceElement( int n, Domain *d );
    /// Destructor.
    virtual ~Structural3DSurfaceElement();

    void initializeFrom( InputRecord &ir ) override;

    MaterialMode giveMaterialMode() override;
    //int computeNumberOfDofs() override;
    void giveDofManDofIDMask( int inode, IntArray &answer ) const override;
    double computeVolumeAround( GaussPoint *gp ) override;

    //double giveCharacteristicLength( const FloatArray &normalToCrackPlane ) override;

    void giveMaterialOrientationAt( FloatArray &x, FloatArray &y, const FloatArray &lcoords );
    void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;
    void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;
    void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    virtual FEICellGeometry *giveCellGeometryWrapper();

protected:
    void computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS ) override;
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;
    void computeGaussPoints() override;
};
} // end namespace oofem
#endif // structural3delement_h
