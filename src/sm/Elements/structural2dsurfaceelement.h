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


namespace oofem {
/**
 * Base class for the planar surface elements: segments in the xy plane on the boundary of a plane strain
 * body (LinePlaneStrainSurface) or, through AxisymSurfaceElement, of an axisymmetric body.
 *
 * Unknowns D_u, D_v per node, rotated to the nodal frame of the element: u_t along the unit tangent t and
 * u_n along the in-plane normal n = (-t_y, t_x), to the left of the direction from node 1 to node 2. The
 * displacement gradient BH is assembled directly in the frame of the surface material (see
 * NLStructuralSurfaceElement::computeMaterialFrame): axis 1 = t, axis 2 = out of the plane (-z, the
 * circumferential direction for axisymmetry), axis 3 = n; these are the e1, e2, e3 of the element manual.
 * F, P and BH have 9 Voigt components (BH is 9 x 2 nno). The surface measure has no thickness (per unit
 * length out of the plane, per radian for axisymmetry).
 */
class Structural2DSurfaceElement : public NLStructuralSurfaceElement
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
    Structural2DSurfaceElement( int n, Domain *d );
    /// Destructor.
    virtual ~Structural2DSurfaceElement();
    /// Sets the number of nodes from the interpolation.
    void postInitialize() override;
    /// Number of nodes of the interpolation.
    int giveNumberOfNodes() const override;
    /**
     * Returns the Cell Geometry Wrapper. Default implementation creates FEIElementGeometryWrapper.
     */
    virtual FEICellGeometry *giveCellGeometryWrapper();

    /// Two DOFs (D_u, D_v) per node of the interpolation.
    int computeNumberOfDofs() override;
    /// D_u, D_v in every node.
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    /// Length of the segment attached to gp, |dX/dxi| times the weight (no thickness).
    double computeVolumeAround(GaussPoint *gp) override;


protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override = 0;
    /// Displacement gradient matrix BH (9 x 2 nno) in the material frame, for the unknowns (u_t, u_n) of the nodes.
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override = 0;
    /// One Gauss integration rule with numberOfGaussPoints points (the element default or nip).
    void computeGaussPoints() override;
};


/**
 * Base class for axisymmetric surface elements: a line in the (r, z) plane, x = radius, y = axis of
 * rotation, describing a surface of revolution. The surface material frame is 1 = tangent,
 * 2 = circumferential, 3 = in-plane normal; the circumferential stretch is r / R = 1 + u_R / R (exact
 * for large deformations), with the radial displacement u_R and the reference radius R. The material
 * mode is _3dMat and the surface measure is per radian, R |dX/ds| ds.
 * Only the large strain formulation (nlgeo 1) is implemented.
 */
class AxisymSurfaceElement : public Structural2DSurfaceElement
{
public:
    AxisymSurfaceElement( int n, Domain *d );
    virtual ~AxisymSurfaceElement() {}
    MaterialMode giveMaterialMode() override { return _3dMat; }
    /// Small-strain path, not implemented (stops with an error).
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    /// Small-strain path, not implemented (stops with an error).
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    /// Reference area per radian attached to gp, R |dX/dxi| times the weight; R must be positive.
    double computeVolumeAround(GaussPoint *gp) override;

protected:
    /// Small-strain B, not implemented (stops with an error).
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override;
    /// BH (9 x 2 nno) with the rows 11 = du_t/ds, 31 = du_n/ds and the hoop row 22 = u_R / R.
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    /// One Gauss integration rule with numberOfGaussPoints points (the element default or nip).
    void computeGaussPoints() override;
};
} // end namespace oofem
#endif // structural2dsurfaceelement_h
