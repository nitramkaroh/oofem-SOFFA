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

#ifndef axilinesurface_h
#define axilinesurface_h

#include "sm/Elements/structural2dsurfaceelement.h"
#include "floatmatrixf.h"

#define _IFT_AxiLineSurface_Name "axilinesurface"

namespace oofem {
class FEI1dLinAxi;

/**
 * Two-node linear line element (FEI1dLinAxi) in the (r, z) plane, x = radius, y = axis, carrying a surface energy
 * on the surface of revolution, given by a SimpleSurfaceCrossSection with a surface material; per radian.
 *
 * The frame is the tangent e1 = (N2 - N1) / |N2 - N1| and the in-plane normal e2 = (-e1_y, e1_x), to the left of
 * the direction N1 -> N2 (the IGA curves use the normal to the right); the unknowns D_u, D_v of each node are
 * rotated to (u_t, u_n). BH is 9 x 4 with the hoop stretch 1 + u_R / R (see AxisymSurfaceElement),
 * 2 Gauss points by default (nip). No element specific input.
 */
class AxiLineSurface : public AxisymSurfaceElement
{
protected:
    static FEI1dLinAxi interpolation;

    /// Nodal coordinates in the element frame (for the geometry wrapper).
    std ::vector<FloatArray> lc;
    /// Cached 2x2 rotation from the global to the element frame, rows e1 (tangent) and e2 (in-plane normal).
    FloatMatrix *GtoLRotationMatrix;


public:
    AxiLineSurface( int n, Domain *d );
    virtual ~AxiLineSurface()
    {
        delete this->GtoLRotationMatrix;
    };
    FEInterpolation *giveInterpolation() const override;


    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_AxiLineSurface_Name; }
    const char *giveClassName() const override { return "AxiLineSurface"; }
    void initializeFrom( InputRecord &ir ) override;

    FloatArray giveNormal( GaussPoint *gp ) const override;

    FEICellGeometry *giveCellGeometryWrapper() override;
    void computeLocalNodalCoordinates( std::vector<FloatArray> &lxy );

    int computeNumberOfDofs() override { return 4; }
    int computeNumberOfGlobalDofs() override { return 4; }

    const FloatMatrix *computeGtoLRotationMatrix();
    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override;
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override;

protected:
    int giveNumberOfIPForMassMtrxIntegration() override { return 2; }


};
} // end namespace oofem
#endif // axilinesurface_h