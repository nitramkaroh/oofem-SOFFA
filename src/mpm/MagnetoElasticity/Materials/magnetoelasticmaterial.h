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

#ifndef magnetoelasticmaterial_h
#define magnetoelasticmaterial_h

#include "material.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "mathfem.h"
#include "matconst.h"
#include "matstatus.h"
#include "valuemodetype.h"
#include <vector>

///@name Input fields for StructuralMaterial
//@{
//@}

namespace oofem {

class GaussPoint;
///@todo Update documentation
class OOFEM_EXPORT MagnetoElasticMaterial : public Material
{
protected:

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
  MagnetoElasticMaterial(int n, Domain *d):Material(n,d){;}

    const char *giveClassName() const override { return "MagnetoElasticMaterial"; }

    void initializeFrom(InputRecord &ir) override = 0;
    void giveInputRecord(DynamicInputRecord &input) override{;}
    //@{
    /// Default implementation relies on giveRealStressVector for second Piola-Kirchoff stress
    virtual std::tuple<FloatArrayF< 9 >,FloatArrayF< 3 >> give_FirstPKStressVector_MagneticInduction_3d(const FloatArrayF<9> &vF, const FloatArrayF<3> &vH, GaussPoint *gp, TimeStep *tStep) const = 0;
     virtual std::tuple<FloatArrayF< 5 >,FloatArrayF< 2 >> give_FirstPKStressVector_MagneticInduction_PlaneStrain(const FloatArrayF<5> &vF, const FloatArrayF<2> &vH, GaussPoint *gp, TimeStep *tStep) const = 0;
     //
     virtual std::tuple<FloatMatrixF< 9, 9 >,  FloatMatrixF< 9, 3 >,  FloatMatrixF< 3, 9 >, FloatMatrixF< 3, 3 >>  giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
     virtual std::tuple<FloatMatrixF< 5, 5 >, FloatMatrixF< 5, 2 >, FloatMatrixF< 2, 5 >, FloatMatrixF< 2, 2 >> giveConstitutiveMatrices_dPdF_dBdH_dPdH_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;

  ////////////////////////////////////////////////////
     virtual FloatArrayF< 8 > give_SecondGradient_FirstPKStressVector_PlaneStrain(const FloatArrayF<8> &vG, GaussPoint *gp, TimeStep *tStep) {    return zeros< 8 >(); }
     virtual FloatArrayF< 27 > give_SecondGradient_FirstPKStressVector_3d(const FloatArrayF<27> &vG, GaussPoint *gp, TimeStep *tStep)  {    return zeros< 27 >(); }
     virtual FloatMatrixF<1,1> give_SecondGradient_ConstitutiveMatrix_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) {    return FloatMatrixF< 1, 1 >();}
     virtual FloatMatrixF<1,1> give_SecondGradient_ConstitutiveMatrix_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) {    return FloatMatrixF< 1, 1 >();}

     virtual std::tuple<FloatArrayF<5>, FloatArrayF<8> > give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep ) { return std::make_tuple( zeros<5>(), zeros<8>() ); }
     virtual std::tuple<FloatArrayF<9>, FloatArrayF<27> > give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vGradF, GaussPoint *gp, TimeStep *tStep ) { return std::make_tuple( zeros<9>(), zeros<27>() ); }
     virtual std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 27>, FloatMatrixF<27, 9>, FloatMatrixF<27, 27> > give_JacobianGradient_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) { return std::make_tuple( FloatMatrixF<9, 9>(), FloatMatrixF<9, 27>(), FloatMatrixF<27, 9>(), FloatMatrixF<27, 27>() ); }
     virtual std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> > give_JacobianGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) { return std::make_tuple( FloatMatrixF<5, 5>(), FloatMatrixF<5, 8>(), FloatMatrixF<8, 5>(), FloatMatrixF<8, 8>() ); }

     
};
} // end namespace oofem
#endif // structuralmaterial_h
