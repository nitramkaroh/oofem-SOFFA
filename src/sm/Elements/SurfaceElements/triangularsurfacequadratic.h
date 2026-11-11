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

#ifndef triangularsurfacequadratic_h
#define triangularsurfacequadratic_h

#include "sm/Elements/structural3dsurfaceelement.h"
#include "floatmatrixf.h"


#define _IFT_TriangularSurfaceQuadratic_Name "triangularsurfacequadratic"

namespace oofem {
class FEI2dTrQuad;

/**
 * Six-node quadratic triangular surface element in 3D (FEI2dTrQuad, vertex nodes 1-3, midside nodes 4-6),
 * for a surface energy given by a SimpleSurfaceCrossSection with a surface material.
 *
 * The element is planar. Its frame is e1 = (N2 - N1) / |N2 - N1|, the normal e3 = e1 x (N3 - N1) normalized and
 * e2 = e3 x e1; the interpolation works with the nodal coordinates in this frame (the coordinate along e3 is
 * not used), and the unknowns D_u, D_v, D_w of each node are rotated to it.
 * BH is 9 x 18 (see Structural3DSurfaceElement), 4 Gauss points by default (nip). No element specific input.
 */
class TriangularSurfaceQuadratic : public Structural3DSurfaceElement
{
protected:
    static FEI2dTrQuad interpolation;

    /// Nodal coordinates in the element frame (for the geometry wrapper).
    std ::vector<FloatArray> lc;
    /// Cached 3x3 rotation from the global to the element frame, rows e1, e2, e3.
    FloatMatrix *GtoLRotationMatrix;


public:
    TriangularSurfaceQuadratic( int n, Domain *d );
    virtual ~TriangularSurfaceQuadratic()
    {
        delete this->GtoLRotationMatrix;
    };
    FEInterpolation *giveInterpolation() const override;


    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_TriangularSurfaceQuadratic_Name; }
    const char *giveClassName() const override { return "TriangularSurfaceQuadratic"; }
    void initializeFrom( InputRecord &ir ) override;

    FloatArray giveNormal( GaussPoint *gp ) const override;

    FEICellGeometry *giveCellGeometryWrapper() override;
    void computeLocalNodalCoordinates( std::vector<FloatArray> &lxy );

    int computeNumberOfDofs() override { return 18; }
    int computeNumberOfGlobalDofs() override { return 18; }

    const FloatMatrix *computeGtoLRotationMatrix();
    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override;
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override;

    int computeGlobalCoordinates( FloatArray &answer, const FloatArray &lcoords ) override;

protected:
    int giveNumberOfIPForMassMtrxIntegration() override { return 4; }


};
} // end namespace oofem
#endif // triangularsurfacequadratic_h