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
 * Abstract base class for hyperelastic surface materials
 **/

class OOFEM_EXPORT BaseHyperElasticSurfaceMaterial
{
protected:



public:
    BaseHyperElasticSurfaceMaterial() { ; }
    ~BaseHyperElasticSurfaceMaterial() { ; }


    double compute_surface_determinant( const Tensor2_3d &F ) const;

    Tensor2_3d compute_surface_cofactor( const Tensor2_3d &F ) const;

    Tensor4_3d compute_surface_dCof_dF( const Tensor2_3d &F ) const;

    Tensor2_3d compute_surface_d_normF_dF( const Tensor2_3d &F ) const;

    double compute_surface_NormF( const Tensor2_3d &F ) const;

    Tensor4_3d compute_surface_d2_normF_dF2( const Tensor2_3d &F ) const;

    Tensor2_3d compute_surface_d_traceU_dF( const Tensor2_3d &F ) const;

    Tensor4_3d compute_surface_d2_traceU_dF2( const Tensor2_3d &F ) const;

    FloatArray compute_deformed_normal( const Tensor2_3d &F ) const;

    Tensor2_3d compute_surface_a_dot_d_normal_dF( const Tensor2_3d &F, const FloatArray &a ) const;

    Tensor4_3d compute_surface_a_dot_d2_normal_d2F( const Tensor2_3d &F, const FloatArray &a ) const;

    Tensor1_3d compute_surface_deformed_normal( const Tensor2_3d &F, const Tensor1_3d &undeformedNormal ) const;

    Tensor2_3d compute_surface_pseudoinverse( const Tensor2_3d &F, const Tensor2_3d &I ) const;

    // Implementations valid also for prestrained reference configuration
    Tensor2_3d compute_surface_cofactor( const Tensor2_3d &F, const Tensor2_3d &I ) const;

    Tensor4_3d compute_surface_dCof_dF( const Tensor2_3d &F, const Tensor2_3d &I ) const;

    Tensor4_3d compute_surface_d2_normF_dF2( const Tensor2_3d &F, const Tensor2_3d &I ) const;

    Tensor2_3d compute_surface_d_normFa_dF( const Tensor2_3d &F, const Tensor1_3d &a ) const;

    Tensor4_3d compute_surface_d2_normFa_dF2( const Tensor2_3d &F, const Tensor1_3d &a ) const;

    Tensor2_3d compute_surface_d_SqrtNormCa_dF( const Tensor2_3d &F, const Tensor1_3d &a ) const;

    Tensor4_3d compute_surface_d2_SqrtNormCa_dF2( const Tensor2_3d &F, const Tensor1_3d &a ) const;
};
} // end namespace oofem
