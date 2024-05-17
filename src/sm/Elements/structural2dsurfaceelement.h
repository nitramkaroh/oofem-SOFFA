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

#ifndef structural2dsurfaceelement_h
#define structural2dsurfaceelement_h

#include "sm/Elements/nlstructuralsurfaceelement.h"
#include "feinterpol2d.h"

#define _IFT_Structural2DSurfaceElement_materialCoordinateSystem "matcs" ///< [optional] Support for material directions based on element orientation.

namespace oofem {
/**
 * Base class for planar 2D surface elements.
 */
class Structural2DSurfaceElement : public NLStructuralSurfaceElement
{
protected:
    /**

     */
    FEICellGeometry *cellGeometryWrapper;

    bool matRotation;

public:
    /**
     * Constructor. Creates element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    Structural2DSurfaceElement( int n, Domain *d );
    /// Destructor.
    virtual ~Structural2DSurfaceElement();
    void postInitialize() override;
    int giveNumberOfNodes() const override;
    /**
     * Returns the Cell Geometry Wrapper. Default inplementation creates FEIElementGeometryWrapper.
     */
    virtual FEICellGeometry *giveCellGeometryWrapper();

    int computeNumberOfDofs() override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    double computeVolumeAround(GaussPoint *gp) override;

    void initializeFrom(InputRecord &ir) override;

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override = 0;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override = 0;
    void computeGaussPoints() override;
};



class AxisymSurfaceElement : public Structural2DSurfaceElement
{
public:
    AxisymSurfaceElement( int n, Domain *d );
    virtual ~AxisymSurfaceElement() {}
    MaterialMode giveMaterialMode() override { return _3dMat; }
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    double computeVolumeAround(GaussPoint *gp) override;

    void computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;
    void computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) override;

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void computeGaussPoints() override;
};
} // end namespace oofem
#endif // structural2delement_h
