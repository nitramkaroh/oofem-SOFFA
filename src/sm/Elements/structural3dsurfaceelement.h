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


    /// _3dMat: the material gets the full 9-component F.
    MaterialMode giveMaterialMode() override;
    /// D_u, D_v, D_w in every node.
    void giveDofManDofIDMask( int inode, IntArray &answer ) const override;
    /// Area attached to gp: the Jacobian of the 2D interpolation in the element plane times the weight (no thickness).
    double computeVolumeAround( GaussPoint *gp ) override;

    /// In-plane material directions from elemLocalCS or the element geometry (reserved for the small-strain path, currently unused).
    void giveMaterialOrientationAt( FloatArray &x, FloatArray &y, const FloatArray &lcoords );
    /// Small-strain stress, reference normal e3 (not used: checkConsistency requires nlgeo 1).
    void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;
    /// Small-strain stiffness, reference normal e3 (not used: checkConsistency requires nlgeo 1).
    void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

    /// Geometry for the interpolation; the derived elements give the nodal coordinates in the element frame.
    virtual FEICellGeometry *giveCellGeometryWrapper();

protected:
    /// Small-strain B (6 x 3 nno, in-plane rows 1, 2, 6 only), reserved for the small-strain path.
    void computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS ) override;
    /// Displacement gradient matrix BH (9 x 3 nno) in the element frame, for the rotated unknowns (u, v, w).
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;
    /// One Gauss integration rule with numberOfGaussPoints points (the element default or nip).
    void computeGaussPoints() override;
};
} // end namespace oofem
#endif // structural3dsurfaceelement_h
