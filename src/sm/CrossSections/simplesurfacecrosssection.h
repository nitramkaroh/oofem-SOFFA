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

#ifndef SimpleSurfaceCrossSection_h
#define SimpleSurfaceCrossSection_h

#include "sm/CrossSections/simplecrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/HyperelasticSurfaceMaterials/hyperelasticsurfacematerial.h"

#define _IFT_SimpleSurfaceCrossSection_Name "simplesurfacecs"

namespace oofem {

/**
 * Cross section of the surface elements: passes F (9 components), grad F (27), the gradient of the reference
 * surface identity (27) and the unit reference normal to the surface material (HyperElasticSurfaceMaterial, see
 * there for the frame and the component order) and returns its stresses and tangents. Only the 3D material mode
 * is supported. Input as SimpleCrossSection; the thickness is not used by the surface elements.
 */
class OOFEM_EXPORT SimpleSurfaceCrossSection : public SimpleCrossSection
{
public:
    SimpleSurfaceCrossSection( int n, Domain *d ) : SimpleCrossSection( n, d ){};

    /// Small strains: stress (6) for the strain (6).
    FloatArray giveRealSurfaceStresses( const FloatArray &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    FloatArrayF<6> giveRealSurfaceStress_3d( const FloatArrayF<6> &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    /// First Piola-Kirchhoff stress (9) for F (9).
    FloatArray giveFirstPKSurfaceStresses( const FloatArray &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    /// Second order stress (27) conjugate to grad F.
    FloatArray giveSecondOrderSurfaceStresses( const FloatArray &vF, const FloatArray &vG, const FloatArray &gradI, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;

    FloatArrayF<9> giveFirstPKSurfaceStress_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    FloatArrayF<27> giveSecondOrderSurfaceStress_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;

    /// Small strains: stiffness (6 x 6), not implemented by the current surface materials.
    FloatMatrixF<6, 6> giveSurfaceStiffnessMatrix_3d( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;

    /// dP/dF (9 x 9).
    FloatMatrixF<9, 9> giveSurfaceStiffnessMatrix_dPdF_3d( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    /// dA/dgradF (27 x 27).
    FloatMatrixF<27, 27> giveSurfaceStiffnessMatrix_dAddF_3d( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    /// dA/dF (27 x 9).
    FloatMatrixF<27, 9> giveSurfaceStiffnessMatrix_dAdF_3d( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;

    const char *giveClassName() const override { return "SimpleSurfaceCrossSection"; }
    const char *giveInputRecordName() const override { return _IFT_SimpleSurfaceCrossSection_Name; }
};
} // end namespace oofem
#endif // SimpleSurfaceCrossSection_h
