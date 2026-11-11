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

#ifndef lineplanestrainsurface_h
#define lineplanestrainsurface_h

#include "sm/Elements/structural2dsurfaceelement.h"
#include "floatmatrixf.h"

#define _IFT_LinePlaneStrainSurface_Name "lineplanestrainsurface"

namespace oofem {
class FEI1dLin;

/**
 * Two-node linear line element (FEI1dLin) carrying a surface energy on the boundary of a plane strain body,
 * given by a SimpleSurfaceCrossSection with a surface material; per unit length out of the plane.
 *
 * The frame is the tangent e1 = (N2 - N1) / |N2 - N1| and the in-plane normal e2 = (-e1_y, e1_x), to the left of
 * the direction N1 -> N2 (the IGA curves use the normal to the right); the unknowns D_u, D_v of each node are
 * rotated to (u_t, u_n). BH is 9 x 4 in the frame of the surface material (see Structural2DSurfaceElement),
 * 2 Gauss points by default (nip). Only the large strain formulation (nlgeo 1) is implemented.
 * No element specific input.
 */
class LinePlaneStrainSurface : public Structural2DSurfaceElement
{
protected:
    static FEI1dLin interpolation;

    /// Nodal coordinates in the element frame (for the geometry wrapper).
    std ::vector<FloatArray> lc;
    /**
     * Cached 2x2 global-to-local rotation of the nodal components in the plane.
     * The element vector uses this block once for each node.
     */
    FloatMatrix *GtoLRotationMatrix;


public:
    LinePlaneStrainSurface( int n, Domain *d );
    virtual ~LinePlaneStrainSurface()
    {
        delete this->GtoLRotationMatrix;
    };
    FEInterpolation *giveInterpolation() const override;


    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_LinePlaneStrainSurface_Name; }
    const char *giveClassName() const override { return "LinePlaneStrainSurface"; }
    void initializeFrom( InputRecord &ir ) override;

    FloatArray giveNormal( GaussPoint *gp ) const override;

    FEICellGeometry *giveCellGeometryWrapper() override;
    void computeLocalNodalCoordinates( std::vector<FloatArray> &lxy );

    int computeNumberOfDofs() override { return 4; }
    int computeNumberOfGlobalDofs() override { return 4; }

    const FloatMatrix *computeGtoLRotationMatrix();
    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override;
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override;


    MaterialMode giveMaterialMode() override { return _3dMat; }
    void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override;
    void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override;

protected:
    void computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS ) override;
    void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override;


protected:
    int giveNumberOfIPForMassMtrxIntegration() override { return 2; }


};
} // end namespace oofem
#endif // lineplanestrainsurface_h
