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

class OOFEM_EXPORT SimpleSurfaceCrossSection : public SimpleCrossSection
{
protected:


public:

    SimpleSurfaceCrossSection( int n, Domain *d ) : SimpleCrossSection( n, d ){};

    FloatArray giveRealSurfaceStresses( const FloatArray &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    FloatArrayF<6> giveRealSurfaceStress_3d( const FloatArrayF<6> &reducedStrain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    FloatArray giveFirstPKSurfaceStresses( const FloatArray &reducedF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
    FloatArrayF<9> giveFirstPKSurfaceStress_3d( const FloatArrayF<9> &reducedF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
 
    //void giveCharSurfaceMaterialStiffnessMatrix( FloatMatrix &answer, const FloatArray &normal, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ); // NOT SURE
    FloatMatrixF<6, 6> giveSurfaceStiffnessMatrix_3d( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const; // NOT SURE

    //void giveCharSurfaceMaterialStiffnessMatrix_dPdF( FloatMatrix &answer, const FloatArray &normal, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep );
    FloatMatrixF<9, 9> giveSurfaceStiffnessMatrix_dPdF_3d( MatResponseMode rMode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;


    const char *giveClassName() const override { return "SimpleSurfaceCrossSection"; }
    const char *giveInputRecordName() const override { return _IFT_SimpleSurfaceCrossSection_Name; }

    //void giveTemperatureVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep ) const; // NOT SURE

  
};
} // end namespace oofem
#endif // SimpleSurfaceCrossSection_h
