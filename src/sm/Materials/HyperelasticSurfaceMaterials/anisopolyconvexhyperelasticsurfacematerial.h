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
//#include "sm/Elements/SurfaceElements/triangularsurface.h"



///@name Input fields for SurfaceTensionMaterial
//@{
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_Name "anisosurfacetensionmat"
//#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_MatDir "matdir"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P1 "p_1"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_P2 "p_2"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_beta "beta"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_eta "eta"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_betaLTF "beta_ltf"
#define _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_etaLTF "eta_ltf"

//@

namespace oofem {
/**
 * Anisotropic membrane: the energy of IsotropicPolyconvexHyperelasticSurfaceMaterial plus
 * eta |F a| + beta sqrt(|C a|) (Horak, Smejkal, Kruzik, JMPS 2025, eq. 97), with the unit fibre direction a
 * tangent to the reference surface. a is the direction anisodir given by the element, or the tangential part of
 * (X - p_1) x (X - p_2): for p_1 and p_2 on an axis of revolution the circumferential direction around it.
 * eta and beta can follow their time functions eta_ltf and beta_ltf, otherwise they are constant.
 * Only the 3D Lagrange surface elements are supported (the direction is rotated with their frame).
 */

class AnisotropicPolyconvexHyperelasticSurfaceMaterial : public IsotropicPolyconvexHyperelasticSurfaceMaterial
{
protected:
    // Material parameters
    //FloatArray matdir;
    FloatArray P1, P2;
    SurfaceMaterialParameter beta, eta;


public:
    AnisotropicPolyconvexHyperelasticSurfaceMaterial( int n, Domain *d );

    void initializeFrom( InputRecord &ir ) override;

    FloatMatrixF<9, 9> give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatArrayF<9> giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;

    /// Fibre direction in the frame of the element; normal is the reference normal in that frame.
    Tensor1_3d computeLocalMaterialDirection( GaussPoint *gp, TimeStep *tStep, const FloatArray &normal ) const;

    //MaterialStatus *CreateStatus( GaussPoint *gp ) const override;

    const char *giveInputRecordName() const override { return _IFT_AnisotropicPolyconvexHyperelasticSurfaceMaterial_Name; }
    const char *giveClassName() const override { return "AnisotropicPolyconvexHyperelasticSurfaceMaterial"; }
};
} // end namespace oofem
#endif
