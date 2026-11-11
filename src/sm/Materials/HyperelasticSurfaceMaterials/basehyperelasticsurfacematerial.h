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

#include "dynamicinputrecord.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "tensor/tensor1.h"
#include "tensor/tensor2.h"
#include "tensor/tensor3.h"
#include "tensor/tensor4.h"
#include "tensor/tensor5.h"
#include "tensor/tensor6.h"

#pragma once


namespace oofem {
class InputRecord;
/**
 * Kinematic helpers of the hyperelastic surface materials. F is the surface deformation gradient (3 x 3, rank 2,
 * F N = 0 for the unit reference normal N), C = F^T F, I (or Is) the surface identity I - N x N and
 * J = |cof F N| the area ratio. The functions with the argument I hold for any N. The one-argument versions of
 * the cofactor, its derivative and d2|F|/dF2, and compute_deformed_normal, are not used by the current materials;
 * the cofactor ones and compute_deformed_normal assume N = e3 (F with a zero third column).
 */
class OOFEM_EXPORT BaseHyperElasticSurfaceMaterial
{
public:
    BaseHyperElasticSurfaceMaterial() { ; }
    ~BaseHyperElasticSurfaceMaterial() { ; }

    /// Area ratio J = sqrt( ((tr C)^2 - tr C^2) / 2 ), i.e. the product of the two nonzero stretches.
    double compute_surface_determinant( const Tensor2_3d &F ) const;

    /// dJ/dF for N = e3 (explicit formula).
    Tensor2_3d compute_surface_cofactor( const Tensor2_3d &F ) const;

    /// d2J/dF2 for N = e3 (explicit formula).
    Tensor4_3d compute_surface_dCof_dF( const Tensor2_3d &F ) const;

    /// d|F|/dF = F / |F|.
    Tensor2_3d compute_surface_d_normF_dF( const Tensor2_3d &F ) const;

    /// Frobenius norm |F|.
    double compute_surface_NormF( const Tensor2_3d &F ) const;

    /// d2|F|/dF2 over all 9 components of F.
    Tensor4_3d compute_surface_d2_normF_dF2( const Tensor2_3d &F ) const;

    /// Current normal as the cross product of the first two columns of F divided by J (N = e3).
    FloatArray compute_deformed_normal( const Tensor2_3d &F ) const;

    /// Current unit normal n = cof F N / |cof F N|.
    Tensor1_3d compute_surface_deformed_normal( const Tensor2_3d &F, const Tensor1_3d &undeformedNormal ) const;

    /// Pseudoinverse F^+ = (dJ/dF)^T / J, with F^+ F = I - N x N.
    Tensor2_3d compute_surface_pseudoinverse( const Tensor2_3d &F, const Tensor2_3d &I ) const;

    /// dJ/dF = F (tr C I - C) / J = J F^{+T}.
    Tensor2_3d compute_surface_cofactor( const Tensor2_3d &F, const Tensor2_3d &I ) const;

    /// d2J/dF2.
    Tensor4_3d compute_surface_dCof_dF( const Tensor2_3d &F, const Tensor2_3d &I ) const;

    /// d2|F|/dF2 = delta_ik I_jl / |F| - F_ij F_kl / |F|^3 (derivatives in the tangent directions only).
    Tensor4_3d compute_surface_d2_normF_dF2( const Tensor2_3d &F, const Tensor2_3d &I ) const;

    /// d|F a|/dF for the fibre direction a.
    Tensor2_3d compute_surface_d_normFa_dF( const Tensor2_3d &F, const Tensor1_3d &a ) const;

    /// d2|F a|/dF2.
    Tensor4_3d compute_surface_d2_normFa_dF2( const Tensor2_3d &F, const Tensor1_3d &a ) const;

    /// d sqrt(|C a|)/dF.
    Tensor2_3d compute_surface_d_SqrtNormCa_dF( const Tensor2_3d &F, const Tensor1_3d &a ) const;

    /// d2 sqrt(|C a|)/dF2.
    Tensor4_3d compute_surface_d2_SqrtNormCa_dF2( const Tensor2_3d &F, const Tensor1_3d &a ) const;
};
} // end namespace oofem
