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

#ifndef structural3dsurfaceelement_h
#define structural3dsurfaceelement_h

#include "sm/Elements/nlstructuralsurfaceelement.h"
#include "feinterpol2d.h"


#define _IFT_Structural3DSurfaceElement_materialCoordinateSystem "matcs" ///< [optional] Material directions from the element orientation (not implemented).

namespace oofem {
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Base class for the planar Lagrange surface elements in 3D (triangles and quadrilaterals).
 *
 * The derived elements define a frame e1, e2 (in the element plane) and e3 (normal) and rotate the
 * unknowns D_u, D_v, D_w of each node to it; the 2D interpolation works with the nodal coordinates in
 * this frame (giveCellGeometryWrapper), so BH is the displacement gradient in the frame of the surface
 * material (1, 2 tangent, 3 normal), 9 x 3 nno. The material mode is _3dMat; the surface measure is the
 * area, without thickness.
 */
class Structural3DSurfaceElement : public NLStructuralSurfaceElement
{
protected:
    /// Set by matcs; material rotation is not implemented (see computeStressVector).
    bool matRotation;
    /// Geometry of the element for the interpolation, created on demand (see giveCellGeometryWrapper).
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
    void giveDofManDofIDMask( int inode, IntArray &answer ) const override;
    double computeVolumeAround( GaussPoint *gp ) override;

    /// In-plane material directions from elemLocalCS or the element geometry (meant for matcs, not used).
    void giveMaterialOrientationAt( FloatArray &x, FloatArray &y, const FloatArray &lcoords );
    void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;
    void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    virtual FEICellGeometry *giveCellGeometryWrapper();

protected:
    void computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS ) override;
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;
    void computeGaussPoints() override;
};
} // end namespace oofem
#endif // structural3dsurfaceelement_h
