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
 *               Copyright (C) 1993 - 2020   Borek Patzak
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

#ifndef anisopolyconvexhyperelasticsurfacematerial_h
#define anisopolyconvexhyperelasticsurfacematerial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/HyperelasticSurfaceMaterials/isopolyconvexhyperelasticsurfacematerial.h"

///@name Input fields for AnisotropicPolyconvexHyperelasticSurfaceMaterial
//@{
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_Name "anisosurfacetensionmat"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P1 "p_1"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P2 "p_2"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_meridional "meridional"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_beta "beta"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_eta "eta"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_betaLTF "beta_ltf"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_etaLTF "eta_ltf"
//@}

namespace oofem {
/**
 * Anisotropic membrane: the energy of IsotropicPolyconvexHyperelasticSurfaceMaterial plus
 * eta |F a| + beta sqrt(|C a|) (Horak, Smejkal, Kruzik, JMPS 2025, eq. 97) with one family of fibres along the
 * unit vector a in the plane of the element (the energy is even in a, so its sign does not matter).
 *
 * The fibre direction is given in global coordinates in one of two ways. It is rotated into the element frame,
 * where the surface identity diag(1, 1, 0) keeps its part in the element plane, which is then normalized:
 * - anisodir 3 dx dy dz on the element (TriangularSurface only) is a general direction for any fibre layout,
 *   e.g. computed for each element by a script.
 * - Otherwise p_1 and p_2 of the material define an axis, a construction made for fibres wound around it (the
 *   sphere in the paper). m = (X - p_1) x (X - p_2) is normal to the plane through the Gauss point X and the axis,
 *   so it points along the circle through X around the axis. By default a is the part of m in the element plane,
 *   fibres along the parallels as in the paper; with meridional 1 it is rotated by 90 degrees in that plane,
 *   fibres along the meridians.
 *   On a surface of revolution about the axis (a sphere with its centre on the axis, a cylinder, a cone) m is
 *   tangent, so the two options are exactly the parallels and the meridians. On a mesh of flat elements this
 *   holds only approximately, and on other surfaces m has a normal part: the meridional fibres then lie in the
 *   plane through the axis and X, the default ones are perpendicular to them in the surface. On the axis, or
 *   where the circle is normal to the element, there is no direction and the analysis stops.
 *
 * Made for elements without prestrain: the direction is taken on the mesh, which is the reference configuration
 * only without it, so prestrain stops with an error.
 * Parameters: the membrane ones (no bending: omega1 and arclengthpar 2 are rejected); p_1 and p_2 (3 components
 * each, different points, given together; required unless every element gives anisodir); meridional (0 or 1,
 * needs p_1 and p_2); eta and beta (default 0) with their time functions eta_ltf and beta_ltf, otherwise constant.
 * Only the 3D Lagrange surface elements are supported (the direction is rotated with their frame).
 */
class AnisotropicPolyconvexHyperelasticSurfaceMaterial : public IsotropicPolyconvexHyperelasticSurfaceMaterial
{
protected:
    FloatArray P1, P2; ///< Points p_1, p_2 on the axis of the fibres (global coordinates).
    bool axisGiven  = false; ///< p_1 and p_2 are given
    bool meridional = false; ///< fibres along the meridians of the axis p_1 p_2 instead of the parallels
    SurfaceMaterialParameter beta, eta; ///< Fibre coefficients beta, eta with their time functions.

public:
    AnisotropicPolyconvexHyperelasticSurfaceMaterial( int n, Domain *d );

    void initializeFrom( InputRecord &ir ) override;

    /// dP/dF of the membrane plus eta d2|F a|/dF2 + beta d2 sqrt(|C a|)/dF2.
    FloatMatrixF<9, 9> give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;
    /// P of the membrane plus eta d|F a|/dF + beta d sqrt(|C a|)/dF.
    FloatArrayF<9> giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;

    /// Unit fibre direction a = (a1, a2, 0) in the frame of the element, in its plane (see the class description).
    Tensor1_3d computeLocalMaterialDirection( GaussPoint *gp, TimeStep *tStep ) const;

    const char *giveInputRecordName() const override { return _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_Name; }
    const char *giveClassName() const override { return "AnisotropicPolyconvexHyperelasticSurfaceMaterial"; }
};
} // end namespace oofem
#endif
