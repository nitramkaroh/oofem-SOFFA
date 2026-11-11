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

#ifndef triangularsurface_h
#define triangularsurface_h

#include "sm/Elements/structural3dsurfaceelement.h"
#include "floatmatrixf.h"


#define _IFT_TriangularSurface_Name "triangularsurface"
#define _IFT_TriangularSurface_f0 "f0"
#define _IFT_TriangularSurface_prestrain "prestrain"
#define _IFT_TriangularSurface_normal "normal"
#define _IFT_TriangularSurface_anisodir "anisodir"

namespace oofem {
class FEI2dTrLin;

/**
 * Three-node linear triangular surface element in 3D (FEI2dTrLin), for a surface energy given by a
 * SimpleSurfaceCrossSection with a surface material.
 *
 * The element is planar. Its frame is e1 = (N2 - N1) / |N2 - N1|, the normal e3 = e1 x (N3 - N1) normalized and
 * e2 = e3 x e1; the interpolation works with the nodal coordinates in this frame (the coordinate along e3 is
 * not used), and the unknowns D_u, D_v, D_w of each node are rotated to it.
 * BH is 9 x 9 (see Structural3DSurfaceElement), 1 Gauss point by default (nip).
 *
 * Optional input:
 * - prestrain 1 f0 9 ...: prestrain deformation gradient F0 (global coordinates, one for the element) from the
 *   original, stress-free configuration to the mesh; the material works with F * F0 and the normal of the
 *   original configuration, which follows from F0 (see NLStructuralSurfaceElement::giveMaterialKinematics).
 * - normal 3 ...: that normal (global coordinates), only checked against the one computed from F0.
 * - anisodir 3 ...: fibre direction for the anisotropic material (global coordinates).
 *
 * giveIPValue reports IST_PrestrainDeformationGradient (F0 in the element frame, R F0 R^T),
 * IST_AnisotropyDirection (anisodir as given) and IST_DisplacementVector (global, at the Gauss point).
 */
class TriangularSurface : public Structural3DSurfaceElement
{
protected:
    static FEI2dTrLin interpolation;

    /// Nodal coordinates in the element frame (for the geometry wrapper).
    std ::vector<FloatArray> lc;
    /// Cached 3x3 rotation from the global to the element frame, rows e1, e2, e3.
    FloatMatrix *GtoLRotationMatrix;
    FloatArrayF<9> F0; // prestrain deformation gradient f0 (global coordinates, one for the element)
    int prestrain = 0;
    int isReferenceNormalGiven = 0;
    FloatArrayF<3> referenceNormal; // normal of the original configuration, checked against F0
    int anisoDirGiven = 0;
    FloatArrayF<3> anisoDir; // direction of anisotropy (global coordinates, one for the element)


public:
    TriangularSurface( int n, Domain *d );
    virtual ~TriangularSurface()
    {
        delete this->GtoLRotationMatrix;
    };
    FEInterpolation *giveInterpolation() const override;


    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_TriangularSurface_Name; }
    const char *giveClassName() const override { return "TriangularSurface"; }
    void initializeFrom( InputRecord &ir ) override;

    FloatArray giveNormal( GaussPoint *gp ) const override;

    FEICellGeometry *giveCellGeometryWrapper() override;
    void computeLocalNodalCoordinates( std::vector<FloatArray> &lxy );

    int computeNumberOfDofs() override { return 9; }
    int computeNumberOfGlobalDofs() override { return 9; }

    const FloatMatrix *computeGtoLRotationMatrix();
    bool computeGtoLRotationMatrix( FloatMatrix &answer ) override;
    int computeLoadGToLRotationMtrx( FloatMatrix &answer ) override;

    int computeGlobalCoordinates( FloatArray &answer, const FloatArray &lcoords ) override;

    int giveSpatialDimension() override { return 3; };

    int giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep ) override;

    int checkConsistency() override;

protected:
    int giveNumberOfIPForMassMtrxIntegration() override { return 3; }


};
} // end namespace oofem
#endif // triangularsurface_h