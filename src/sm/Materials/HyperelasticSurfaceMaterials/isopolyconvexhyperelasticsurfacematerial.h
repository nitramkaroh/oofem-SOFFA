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
//#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralsurfacems.h"
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
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_omega1 "omega1"

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
    double omega1; // For second order stress tensor


public:
    IsotropicPolyconvexHyperelasticSurfaceMaterial( int n, Domain *d );

    void initializeFrom( InputRecord &ir ) override;

    FloatMatrixF<6, 6> give3dSurfaceMaterialStiffnessMatrix( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for large strains only" ); }
    FloatArrayF<6> giveRealSurfaceStressVector_3d( const FloatArrayF<6> &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for large strains only" ); }
    
    FloatMatrixF<9, 9> give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatArrayF<9> giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;
    Tensor2_3d giveFirstPKSurfaceStressBending( const Tensor2_3d &F, Tensor2_3d &IdentityUD, const Tensor1_3d &normal, const Tensor3_3d &G, const Tensor2_3d &MK ) const;


    FloatArrayF<27> giveSecondOrderSurfaceStressVector_3d( const FloatArrayF<9> &reducedF, const FloatArrayF<27> &reducedG, const FloatArrayF<27> &gradI, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;

    FloatMatrixF<27, 27> give3dSurfaceMaterialStiffnessMatrix_dAddF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatMatrixF<27, 9> give3dSurfaceMaterialStiffnessMatrix_dAdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;


    MaterialStatus *CreateStatus( GaussPoint *gp ) const override;

    double computeEnergy( const Tensor2_3d &F, GaussPoint *gp, TimeStep *tStep ) const;

    const char *giveInputRecordName() const override { return _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_Name; }
    const char *giveClassName() const override { return "IsotropicPolyconvexHyperelasticSurfaceMaterial"; }

    Tensor2_3d DoRTRtTransform( const Tensor2_3d &T, const FloatMatrix &R ) const;

    void TransformTensor1( const Tensor1_3d &T, Tensor1_3d &Ttrans, GaussPoint *gp ) const;
    void TransformTensor2( const Tensor2_3d &T, Tensor2_3d &Ttrans, GaussPoint *gp ) const;
    void TransformTensor4( const Tensor4_3d &T, Tensor4_3d &Ttrans, GaussPoint *gp ) const;

    Tensor4_3d computeTangentPK1Num( const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, double EPSILON, const Tensor2_3d &P, const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    Tensor4_3d computeTangentPK1NumB( double EPSILON, const Tensor2_3d &F, Tensor2_3d &IdentityUD, const Tensor1_3d &normal, const Tensor3_3d &G, const Tensor2_3d &MK ) const;
    Tensor6_3d computeTangentSecondOrderNum( const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, double EPSILON, const Tensor2_3d &P, const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    Tensor5_3d computeTangentMixedNum( const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, double EPSILON, const Tensor2_3d &P, const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    Tensor5_3d computeTangentMixedNum2( const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, double EPSILON, const Tensor2_3d &P, const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;

    Tensor3_3d dndFnum( double EPSILON, const Tensor2_3d &F, const Tensor1_3d &normal ) const;
    Tensor4_3d dFinvdFnum( double EPSILON, const Tensor2_3d &F, const Tensor2_3d &IdentityUD ) const;

};
} // end namespace oofem
#endif
