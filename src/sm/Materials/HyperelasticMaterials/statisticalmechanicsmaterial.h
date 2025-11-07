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


#ifndef statisticalmechanicsmaterial_h
#define statisticalmechanicsmaterial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "basehyperelasticmaterial.h"
#include "sm/Materials/MixedPressure/mixedpressurematerialextensioninterface.h"



///@name Input fields for StatisticalMechanicsMaterial
//@{
#define _IFT_StatisticalMechanicsMaterial_Name "statisticalmechmat"
#define _IFT_StatisticalMechanicsMaterial_nPoints "npoints"
//#define _IFT_StatisticalMechanicsMaterial_N0 "N0"
//#define _IFT_StatisticalMechanicsMaterial_n "n"
//#define _IFT_StatisticalMechanicsMaterial_l "l"
//#define _IFT_StatisticalMechanicsMaterial_T "T"
//#define _IFT_StatisticalMechanicsMaterial_nPoints "nPoints"
//#define _IFT_StatisticalMechanicsMaterial_R "R"
#define _IFT_StatisticalMechanicsMaterial_rho0 "rho0"
#define _IFT_StatisticalMechanicsMaterial_G "g"
    //@}

namespace oofem {

class StatisticalMechanicsMaterial : public StructuralMaterial, public BaseHyperElasticMaterial, public LargeStrainMixedPressureMaterialExtensionInterface
{
protected:
    //double N0 = 1e10; // chain density
    //int n  = 1000; // number of links per chain
    //double l  = 1e-8; // chain link length
    //double R    = 1e-7; // end-to-end length
    //double T  = 293.15; // Thermodynamic temperature
    //double k  = 1.380649e-23; // Boltzman constant
    double G = 30.; // k*T*N0
    double rho0 = 0.29;  // R/n/l

    int nPoints = 42; // number of integration points, foer the whole sphere

    FloatArray integrationWeights;
    std::vector<Tensor1_3d> integrationPositions;

    std::vector<Tensor2_3d> RoR;
    std::vector<Tensor4_3d> RoRoRoR;


public:
    StatisticalMechanicsMaterial( int n, Domain *d );

    void initializeFrom( InputRecord &ir ) override;

    FloatMatrixF<6, 6> give3dMaterialStiffnessMatrix( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for large strains only" ); }

    FloatArrayF<6> giveRealStressVector_3d( const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for large strains only" ); }


    FloatMatrixF<9, 9> give3dMaterialStiffnessMatrix_dPdF( MatResponseMode,
        GaussPoint *gp,
        TimeStep *tStep ) const override;

    FloatArrayF<9> giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const override;

    MaterialStatus *CreateStatus( GaussPoint *gp ) const override;

    const char *giveInputRecordName() const override { return _IFT_StatisticalMechanicsMaterial_Name; }
    const char *giveClassName() const override { return "StatisticalMechanicsMaterial"; }

    void initializeData( int nPoints);

    void constructTensorProducts();

    ////////////////////////////////// Functions for Mixed Pressure formulation
    std::tuple<FloatArrayF<9>, double> giveFirstPKStressVector_3d( const FloatArrayF<9> &F, double pressure, GaussPoint *gp, TimeStep *tStep ) override;
    std::tuple<FloatArrayF<5>, double> giveFirstPKStressVector_PlaneStrain( const FloatArrayF<5> &F, double pressure, GaussPoint *gp, TimeStep *tStep ) override;
    ///////////////////////////////////////////////////////////////////////
    std::tuple<FloatMatrixF<9, 9>, FloatArrayF<9>, double> giveLargeStrainMixedPressureConstitutiveMatrices_3d( double pressure, MatResponseMode, GaussPoint *gp, TimeStep *tStep ) override;
    std::tuple<FloatMatrixF<5, 5>, FloatArrayF<5>, double> giveLargeStrainMixedPressureConstitutiveMatrices_PlaneStrain( double pressure, MatResponseMode, GaussPoint *gp, TimeStep *tStep ) override;

    Interface *giveInterface( InterfaceType t ) override
    {
        if ( t == LargeStrainMixedPressureMaterialExtensionInterfaceType ) {
            return static_cast<LargeStrainMixedPressureMaterialExtensionInterface *>( this );
        } else {
            return nullptr;
        }
    }
};
} // end namespace oofem
#endif
