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
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_beta "beta"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_eta "eta"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_betaLTF "beta_ltf"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_etaLTF "eta_ltf"
//@}

namespace oofem {
/**
 * Anisotropic membrane: the energy of IsotropicPolyconvexHyperelasticSurfaceMaterial plus
 * eta |F a| + beta sqrt(|C a|) (Horak, Smejkal, Kruzik, JMPS 2025, eq. 97), with the reference fibre vector a.
 * In the element frame a = N x d + (N . d) N, a rotation of d by 90 degrees about the material reference normal N.
 * Without anisodir, d is the normalized vector ((X - p_1) x (X - p_2)) x n_e, rotated into the element frame.
 * Without prestrain this gives a unit tangential fibre direction, circumferential for p_1 and p_2 on an axis
 * of revolution. With prestrain, N belongs to the original configuration and this construction need not make
 * a tangent to that configuration.
 * Parameters: the membrane ones (no bending: omega1 and arclengthpar 2 are rejected), p_1 and p_2 (required,
 * 3 components each, different points), eta and beta (default 0) with their time functions eta_ltf and beta_ltf,
 * otherwise constant.
 * Only the 3D Lagrange surface elements are supported (the direction is rotated with their frame).
 * @note anisodir of TriangularSurface supplies d in global coordinates before rotation into the element frame.
 * Its magnitude is retained; it must have a finite nonzero tangential part. The vector a is unit and tangent
 * only when d is unit and perpendicular to N; p_1 and p_2 are still required.
 */
class AnisotropicPolyconvexHyperelasticSurfaceMaterial : public IsotropicPolyconvexHyperelasticSurfaceMaterial
{
protected:
    FloatArray P1, P2;
    SurfaceMaterialParameter beta, eta;

public:
    AnisotropicPolyconvexHyperelasticSurfaceMaterial( int n, Domain *d );

    void initializeFrom( InputRecord &ir ) override;

    FloatMatrixF<9, 9> give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatArrayF<9> giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;

    /// Fibre direction in the frame of the element; normal is the reference normal in that frame.
    Tensor1_3d computeLocalMaterialDirection( GaussPoint *gp, TimeStep *tStep, const FloatArray &normal ) const;

    const char *giveInputRecordName() const override { return _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_Name; }
    const char *giveClassName() const override { return "AnisotropicPolyconvexHyperelasticSurfaceMaterial"; }
};
} // end namespace oofem
#endif
