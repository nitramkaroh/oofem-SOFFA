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

#include "sm/CrossSections/simplecrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"

///@name Input fields for SimpleCrossSection
//@{
#define _IFT_ThirdMediumCrossSection_Name "thirdmediumcs"
//@}

namespace oofem {
/**
 * Class implementing a cross section model for third medium, supporting integration of local second gradient terms
 */
class OOFEM_EXPORT ThirdMediumCrossSection : public SimpleCrossSection
{

public:
    /**
     * Constructor.
     * @param n Cross section number.
     * @param d Associated domain.
     */
    ThirdMediumCrossSection( int n, Domain *d ) :
        SimpleCrossSection( n, d ) {}

    const char *giveClassName() const override { return "ThirdMediumCrossSection"; }
    const char *giveInputRecordName() const override { return _IFT_ThirdMediumCrossSection_Name; }

    // terms with second gradient handling
    FloatArray give_SecondGradient_FluxVector( const FloatArray &vGrad, GaussPoint *gp, TimeStep *tStep ) const;
    FloatArrayF<27> give_SecondGradient_FluxVector_3d( const FloatArrayF<27> &vG, GaussPoint *gp, TimeStep *tStep ) const;
    FloatArrayF<8> give_SecondGradient_FluxVector_PlaneStrain( const FloatArrayF<8> &vG, GaussPoint *gp, TimeStep *tStep ) const;

    void give_SecondGradient_dFlux_dGrad( FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep ) const;
    FloatMatrixF<1, 1> give_SecondGradient_ConstitutiveMatrix_3d( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const;
    FloatMatrixF<1, 1> give_SecondGradient_ConstitutiveMatrix_PlaneStrain( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const;

    //terms with first and second gradient handling
    std::tuple<FloatArray, FloatArray> give_FirstSecondGradient_FluxVectors( const FloatArray &vF, const FloatArray &vGradF, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatArrayF<9>, FloatArrayF<27> > give_FirstSecondGradient_FluxVectors_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vGradF, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatArrayF<5>, FloatArrayF<8> > give_FirstSecondGradient_FluxVectors_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep ) const;

    std::tuple<FloatMatrix, FloatMatrix, FloatMatrix, FloatMatrix> give_FirstSecondGradient_dFluxes_dGrads( MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 27>, FloatMatrixF<27, 9>, FloatMatrixF<27, 27> > give_FirstSecondGradient_ConstitutiveMatrices_3d( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> > give_FirstSecondGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const;

    //terms with F and Fbar handling
    std::tuple<FloatArray, FloatArray> give_Fbar_FluxVectors( const FloatArray &vF, const FloatArray &vFbar, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatArrayF<9>, FloatArrayF<9> > give_Fbar_FluxVectors_3d( const FloatArrayF<9> &vF, const FloatArrayF<9> &vFbar, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatArrayF<5>, FloatArrayF<5> > give_Fbar_FluxVectors_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<5> &vFbar, GaussPoint *gp, TimeStep *tStep ) const;

    std::tuple<FloatMatrix, FloatMatrix, FloatMatrix, FloatMatrix> give_Fbar_dFluxes_dGrads( MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9> > give_Fbar_ConstitutiveMatrices_3d( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 5>, FloatMatrixF<5, 5>, FloatMatrixF<5, 5> > give_Fbar_ConstitutiveMatrices_PlaneStrain( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const;

    //terms with linear elasticity handling
    FloatArray give_LinearElasticity_FluxVector( const FloatArray &eps, GaussPoint *gp, TimeStep *tStep ) const;
    FloatArrayF<6> give_LinearElasticity_FluxVector_3d( const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep ) const;
    FloatArrayF<4> give_LinearElasticity_FluxVector_PlaneStrain( const FloatArrayF<4> &eps, GaussPoint *gp, TimeStep *tStep ) const;

    FloatMatrix give_LinearElasticity_dFlux_dGrad( MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep ) const;
    FloatMatrixF<6, 6> give_LinearElasticity_ConstitutiveMatrix_3d( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const;
    FloatMatrixF<4, 4> give_LinearElasticity_ConstitutiveMatrix_PlaneStrain( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const;

};
} // end namespace oofem