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
#pragma once
#include "./HyperelasticSurfaceMaterials/basehyperelasticsurfacematerial.h"
#include "./structuralmaterial.h"
#include "./structuralms.h"

namespace oofem {

/**
 * Base class of the surface materials, evaluated by the cross section simplesurfacecs.
 *
 * All quantities are three-dimensional and given in one orthonormal frame chosen by the element: the surface
 * deformation gradient F (9 components, F N = 0), the unit reference normal N (argument normal) and, for the
 * second gradient, grad F with the components dF_ij/dX_k (27) and the gradient of the reference surface identity
 * (27). The 9 components are in the order 11 22 33 23 13 12 32 31 21, the 27 components are three such blocks
 * for k = 1, 2, 3. The energy is per unit reference area.
 *
 * Large strains: the first Piola-Kirchhoff stress P (9) with dP/dF (9 x 9), and for second gradient materials
 * the second order stress A conjugate to grad F (27) with dA/dgradF (27 x 27) and dA/dF (27 x 9); the element
 * uses dP/dgradF = (dA/dF)^T. The tangents are evaluated at the F (and grad F) stored in the temporary status by
 * the last stress evaluation. The surface Cauchy stress (output IST_CauchyStressTensor) is the one stored in the
 * status by the material. Small strains (6 components) are only declared; the current materials do not implement them.
 */
class OOFEM_EXPORT HyperElasticSurfaceMaterial : public BaseHyperElasticSurfaceMaterial, public StructuralMaterial
{
public:
    HyperElasticSurfaceMaterial( int n, Domain *d );
    ~HyperElasticSurfaceMaterial() { }
    void initializeFrom( InputRecord &ir ) override;

    /// Small strains: stress (6) for the strain (6).
    virtual FloatArrayF<6> giveRealSurfaceStressVector_3d( const FloatArrayF<6> &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const { OOFEM_ERROR( "not implemented " ); };
    /// First Piola-Kirchhoff stress P (9) for F (9). The default converts the small strain stress, used as the second Piola-Kirchhoff stress.
    virtual FloatArrayF<9> giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    /// dP/dF (9 x 9) at the F stored in the temporary status.
    virtual FloatMatrixF<9, 9> give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    /// Small strains: stiffness (6 x 6).
    virtual FloatMatrixF<6, 6> give3dSurfaceMaterialStiffnessMatrix( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const { OOFEM_ERROR( "not implemented " ); }
    int giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep ) override;

    /// Second order stress A (27) for F (9), grad F (27) and the gradient of the reference surface identity gradI (27).
    virtual FloatArrayF<27> giveSecondOrderSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
    {
        OOFEM_ERROR( "not implemented " );
    }
    /// dA/dgradF (27 x 27).
    virtual FloatMatrixF<27, 27> give3dSurfaceMaterialStiffnessMatrix_dAddF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
    {
        OOFEM_ERROR( "not implemented " );
    }
    /// dA/dF (27 x 9), its transpose is dP/dgradF.
    virtual FloatMatrixF<27, 9> give3dSurfaceMaterialStiffnessMatrix_dAdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
    {
        OOFEM_ERROR( "not implemented " );
    }
};
} // end namespace oofem
