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
class FEI1dLin;

/**
 * Two-node linear line element (FEI1dLin) in the (r, z) plane, x = radius, y = axis, carrying a surface energy
 * on the surface of revolution, given by a SimpleSurfaceCrossSection with a surface material; per radian.
 *
 * The frame is the tangent e1 = (N2 - N1) / |N2 - N1| and the in-plane normal e2 = (-e1_y, e1_x), to the left of
 * the direction N1 -> N2 (the IGA curves use the normal to the right); the unknowns D_u, D_v of each node are
 * rotated to (u_t, u_n). The in-plane normal is the third axis of the material frame (see
 * Structural2DSurfaceElement). BH is 9 x 4 with the hoop stretch 1 + u_R / R (see AxisymSurfaceElement),
 * 2 Gauss points by default (nip). No element specific input.
 */
class AxiLineSurface : public AxisymSurfaceElement
{
protected:
    static FEI1dLin interpolation;

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

    /// In-plane normal (-e1_y, e1_x, 0) in global coordinates.
    FloatArray giveNormal( GaussPoint *gp ) const override;

    /// Vertex list geometry with the nodal coordinates in the element frame (see computeLocalNodalCoordinates).
    FEICellGeometry *giveCellGeometryWrapper() override;
    /// Nodal coordinates rotated to the element frame, lxy[i] = R X_i (2 components, no translation).
    void computeLocalNodalCoordinates( std::vector<FloatArray> &lxy );

    int computeNumberOfDofs() override { return 4; }
    int computeNumberOfGlobalDofs() override { return 4; }

    /// The cached 2 x 2 rotation R (rows e1 and the in-plane normal), computed on the first call.
    const FloatMatrix *computeGtoLRotationMatrix();
    /// Rotation of the element unknowns, 4 x 4 block diagonal with R: (u_t, u_n) = R (u, v) in each node.
    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override;
    /// R itself (2 x 2): rotates the pointwise components of a body load, and gives the material frame.
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override;

protected:
    int giveNumberOfIPForMassMtrxIntegration() override { return 2; }


};
} // end namespace oofem
#endif // axilinesurface_h