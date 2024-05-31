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

#ifndef isopolyconvexhyperelasticsurfacematerial_h
#define isopolyconvexhyperelasticsurfacematerial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/HyperelasticSurfaceMaterials/hyperelasticsurfacematerial.h"


///@name Input fields for SurfaceTensionMaterial
//@{
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_Name "surfacetensionmat"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gamma1 "gamma1"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gamma2 "gamma2"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha1 "alpha1"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha2 "alpha2"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_delta "delta"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gammaLTF "gamma_ltf"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gaminit "gaminit"

//@}

namespace oofem {
/**
 * This class implements basic surface tension material for fluids
 */
class IsotropicPolyconvexHyperelasticSurfaceMaterial : public HyperElasticSurfaceMaterial
{
protected:
    // Material parameters
    /*double gamma;*/
    double gamma1, gamma2, gaminit;
    int gamma_ltf = 0;
    double alpha1, alpha2; 
    double delta; 


public:
    IsotropicPolyconvexHyperelasticSurfaceMaterial( int n, Domain *d );

    void initializeFrom( InputRecord &ir ) override;

    FloatMatrixF<6, 6> give3dSurfaceMaterialStiffnessMatrix( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for large strains only" ); }
    FloatArrayF<6> giveRealSurfaceStressVector_3d( const FloatArrayF<6> &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for large strains only" ); }
    
    FloatMatrixF<9, 9> give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatArrayF<9> giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;

    MaterialStatus *CreateStatus( GaussPoint *gp ) const override;

    double computeEnergy( const Tensor2_3d &F, GaussPoint *gp, TimeStep *tStep ) const;

    const char *giveInputRecordName() const override { return _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_Name; }
    const char *giveClassName() const override { return "IsotropicPolyconvexHyperelasticSurfaceMaterial"; }
};
} // end namespace oofem
#endif