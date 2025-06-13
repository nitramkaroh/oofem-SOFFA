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

#include "reducedhardmagneticmooneyrivlincompressiblematerial_consistent.h"
#include "mpm/MagnetoElasticity/Materials/magnetoelasticmaterialstatus.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"
#include "domain.h"
#include "function.h"



namespace oofem {
REGISTER_Material(ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent);

  ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent::ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent( int n, Domain *d ) : MooneyRivlinCompressibleMaterial(n, d)
{
  
}


 
  
void
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent::initializeFrom( InputRecord &ir )
{   

     demagFactor.resize(3);
     IR_GIVE_OPTIONAL_FIELD( ir, demagFactor, _IFT_ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent_Demagnetization );
     if(demagFactor.giveSize() != 3) {
       OOFEM_ERROR("Wrong size of the demagnetization factor");
     }

     FloatArray b;
     IR_GIVE_OPTIONAL_FIELD( ir, b, _IFT_ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent_AppliedInduction);
     if(b.giveSize() != 3) {
       OOFEM_ERROR("Wrong size of the applied induction");
     }
     b_app = Tensor1_3d(FloatArrayF<3>( b ));
     //
     IR_GIVE_FIELD( ir, bltf_index, _IFT_ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent_bltf );
     //
     IR_GIVE_FIELD( ir, hmMatNumber, _IFT_ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent_hmMatNumber );
}


void
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent::postInitialize()
{
  hmMat = dynamic_cast< HardMagneticMooneyRivlinCompressibleMaterial * >( domain->giveMaterial(hmMatNumber) );
  if (!hmMat) {
    OOFEM_ERROR("Material ID %d is not a HardMagneticMooneyRivlinCompressibleMaterial", hmMatNumber);
  }
}


MaterialStatus *
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent :: CreateStatus(GaussPoint *gp) const
{
    return new ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent(gp, domain, hmMatNumber);
}

  
FloatArrayF<9>
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent::giveFirstPKStressVector_3d( const FloatArrayF<9> &vF, GaussPoint *gp, TimeStep *tStep ) const
{
  ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent *status = static_cast< ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent * >( this->giveStatus(gp) );
  Tensor2_3d F( vF );
  //
  auto [J, G] = F.compute_determinant_and_cofactor();
  auto M = hmMat->giveLagrangianMagnetization(tStep);
  //
  double bload_level = this->giveDomain()->giveFunction( bltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
  //
  Tensor1_3d Happ;
  Happ(i_3)  = 1/mu_0 * bload_level * b_app( j_3 ) * G(j_3,i_3);
  //
  FloatArrayF<3> vHm(M.to_voigt_form()), vHapp(Happ.to_voigt_form()), vH;
  //
  vH.at(1) = -demagFactor.at(1) * vHm.at(1) + vHapp.at(1);
  vH.at(2) = -demagFactor.at(2) * vHm.at(2) + vHapp.at(2);
  vH.at(3) = -demagFactor.at(3) * vHm.at(3) + vHapp.at(3);
  //
  auto [vP, B] = hmMat->give_FirstPKStressVector_MagneticInduction_3d( vF, vH, status->giveSlaveGaussPoint(), tStep );
  //
  return vP;
}


  Tensor1_3d
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent :: giveLagrangianMagnetization(TimeStep *tStep) const
  {

    return hmMat->giveLagrangianMagnetization(tStep);

  }


  std::tuple<Tensor1_3d,Tensor1_3d>
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent :: give_EulerianMagnetization_MagneticInduction_FromF(const Tensor2_3d &F, TimeStep *tStep)
{

  
    auto M = hmMat->giveLagrangianMagnetization(tStep);
    Tensor1_3d m;
    auto [J, G] = F.compute_determinant_and_cofactor();
    if(hmMat->givePullBackType() ==   PBT_F) { //
       m(i_3) = 1./J * F(i_3,j_3) * M(j_3);
    } else if (hmMat->givePullBackType() ==  PBT_iFt) {
      m(i_3) = 1./J * G(i_3, j_3) * M(j_3);
    } else if( hmMat->givePullBackType() == PBT_R) {
      OOFEM_ERROR("Not Implemented");
    }
    //
     double bload_level = this->giveDomain()->giveFunction( bltf_index )->evaluateAtTime( tStep->giveIntrinsicTime() );
 
    Tensor1_3d Happ, H, h, b;
    Happ(i_3)  = 1/mu_0 * bload_level * b_app( j_3 ) * G(j_3,i_3);
    //
    H(0) = -demagFactor.at(1) * M(0) + Happ(0);
    H(1) = -demagFactor.at(2) * M(1) + Happ(1);
    H(2) = -demagFactor.at(3) * M(2) + Happ(2);    
    //
    h(i_3) = 1./J * G(i_3, j_3) * H(j_3);
    //
    b(i_3) = mu_0 * (h(i_3) + m(i_3));
    //
    return std::make_tuple( m, b );

}




FloatMatrixF<9, 9>
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent::give3dMaterialStiffnessMatrix_dPdF( MatResponseMode matResponseMode, GaussPoint *gp, TimeStep *tStep ) const
{
  ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent *status = static_cast< ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent * >( this->giveStatus(gp) );
  auto [dPdF, dPdH, dBdF, dBdH] = hmMat->giveConstitutiveMatrices_dPdF_dBdH_dPdH_3d(matResponseMode, status->giveSlaveGaussPoint(), tStep );
  //
  return dPdF;
}
  

int
ReducedHardMagneticMooneyRivlinCompressibleMaterial_Consistent :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
  return hmMat->giveIPValue(answer, gp, type, tStep);
 
}

  /////status
ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent  :: ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent (GaussPoint *g, Domain *d, int s) :
    StructuralMaterialStatus(g),
    domain(d),
    slaveMat(s)
{

  FloatArray coords  = {0.,0.,0.};
  slaveGp = new GaussPoint( NULL, 1, coords, 1.0, _3dMat);
  this->createSlaveStatus();

}


  ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent :: ~ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent()
{
  delete slaveGp;
  delete slaveStatus;
}


void
ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent :: createSlaveStatus()
{
    slaveStatus = new MagnetoElasticMaterialStatus(slaveGp);
    
}

  

void
ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    auto sMat = static_cast< HardMagneticMooneyRivlinCompressibleMaterial * >( domain->giveMaterial(slaveMat) );
    auto mS = sMat->giveStatus(slaveGp);
    mS->printOutputAt(file, tStep);
}


  

  // initializes temporary variables based on their values at the previous equlibrium state
void
ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent :: initTempStatus()
{
    auto sMat = static_cast< HardMagneticMooneyRivlinCompressibleMaterial * >( domain->giveMaterial(slaveMat) );
    auto mS = sMat->giveStatus(slaveGp);
    mS->initTempStatus();
    StructuralMaterialStatus :: initTempStatus();
}


// updates internal variables when equilibrium is reached
void
ReducedHardMagneticMooneyRivlinCompressibleMaterialStatus_Consistent :: updateYourself(TimeStep *tStep)
{
    auto sMat = static_cast< HardMagneticMooneyRivlinCompressibleMaterial * >( domain->giveMaterial(slaveMat) );
    auto mS = sMat->giveStatus(slaveGp);
    mS->updateYourself(tStep);
    StructuralMaterialStatus :: updateYourself(tStep);
}

  

  
} // end namespace oofem
