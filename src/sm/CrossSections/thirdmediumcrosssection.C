

#include "thirdmediumcrosssection.h"
#include "sm/Materials/HyperelasticMaterials/ThirdMediumMaterials/thirdmediummaterial.h"
#include "classfactory.h"

namespace oofem {

REGISTER_CrossSection( ThirdMediumCrossSection );

FloatArray
ThirdMediumCrossSection::give_SecondGradient_FluxVector( const FloatArray &vGrad, GaussPoint *gp, TimeStep *tStep ) const
{
  MaterialMode mode = gp->giveMaterialMode();

  if ( mode == _3dMat ) {
    return this->give_SecondGradient_FluxVector_3d( vGrad, gp, tStep );
  } else if ( mode == _PlaneStrain ) {
    return this->give_SecondGradient_FluxVector_PlaneStrain( vGrad, gp, tStep );
  } else {
    OOFEM_ERROR( "unsupported mode" );
  }
}


FloatArrayF<27>
ThirdMediumCrossSection ::give_SecondGradient_FluxVector_3d( const FloatArrayF<27> &vG, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_SecondGradient_FirstPKStressVector_3d( vG, gp, tStep );
}


FloatArrayF<8>
ThirdMediumCrossSection ::give_SecondGradient_FluxVector_PlaneStrain( const FloatArrayF<8> &vG, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_SecondGradient_FirstPKStressVector_PlaneStrain( vG, gp, tStep );
}


void ThirdMediumCrossSection ::give_SecondGradient_dFlux_dGrad( FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep )
{
  MaterialMode mode = gp->giveMaterialMode();
  if ( mode == _3dMat ) {
    answer = this->give_SecondGradient_ConstitutiveMatrix_3d( rmode, gp, tStep );
  } else if ( mode == _PlaneStrain ) {
    answer = this->give_SecondGradient_ConstitutiveMatrix_PlaneStrain( rmode, gp, tStep );
  } else {
    OOFEM_ERROR( "unsupported mode" );
  }
}


FloatMatrixF<1, 1>
ThirdMediumCrossSection ::give_SecondGradient_ConstitutiveMatrix_3d( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_SecondGradient_ConstitutiveMatrix_3d( rMode, gp, tStep );
}


FloatMatrixF<1, 1>
ThirdMediumCrossSection ::give_SecondGradient_ConstitutiveMatrix_PlaneStrain( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_SecondGradient_ConstitutiveMatrix_PlaneStrain( rMode, gp, tStep );
}

//////////////////////////////////////JACOBIAN_GRADIENT//////////////////////////////////////////////////////////////
std::tuple<FloatArray, FloatArray>
ThirdMediumCrossSection::give_JacobianGradient_FluxVectors( const FloatArray &vF, const FloatArray &vGradF, GaussPoint *gp, TimeStep *tStep ) const
{
  MaterialMode mode = gp->giveMaterialMode();

  if ( mode == _3dMat ) {
    return this->give_JacobianGradient_FluxVectors_3d( vF, vGradF, gp, tStep );
  } else if ( mode == _PlaneStrain ) {
    return this->give_JacobianGradient_FluxVectors_PlaneStrain( vF, vGradF, gp, tStep );
  } else {
    OOFEM_ERROR( "unsupported mode" );
  }
}


std::tuple<FloatArrayF<9>, FloatArrayF<27> >
ThirdMediumCrossSection ::give_JacobianGradient_FluxVectors_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vGradF, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_3d( vF, vGradF, gp, tStep );
}


std::tuple<FloatArrayF<5>, FloatArrayF<8> >
ThirdMediumCrossSection ::give_JacobianGradient_FluxVectors_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( vF, vGradF, gp, tStep );
}


std::tuple<FloatMatrix, FloatMatrix, FloatMatrix, FloatMatrix> ThirdMediumCrossSection ::give_JacobianGradient_dFluxes_dGrads( MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep )
{
  MaterialMode mode = gp->giveMaterialMode();
  if ( mode == _3dMat ) {
    return this->give_JacobianGradient_ConstitutiveMatrices_3d( rmode, gp, tStep );
  } else if ( mode == _PlaneStrain ) {
    return this->give_JacobianGradient_ConstitutiveMatrices_PlaneStrain( rmode, gp, tStep );
  } else {
    OOFEM_ERROR( "unsupported mode" );
  }
}


std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 27>, FloatMatrixF<27, 9>, FloatMatrixF<27, 27> >
ThirdMediumCrossSection ::give_JacobianGradient_ConstitutiveMatrices_3d( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_JacobianGradient_ConstitutiveMatrices_3d( rMode, gp, tStep );
}


std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> >
ThirdMediumCrossSection ::give_JacobianGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_JacobianGradient_ConstitutiveMatrices_PlaneStrain( rMode, gp, tStep );
}

//////////////////////////////////////FBAR//////////////////////////////////////////////////////////////
std::tuple<FloatArray, FloatArray>
ThirdMediumCrossSection::give_Fbar_FluxVectors( const FloatArray &vF, const FloatArray &vFbar, GaussPoint *gp, TimeStep *tStep ) const
{
  MaterialMode mode = gp->giveMaterialMode();

  if ( mode == _3dMat ) {
    return this->give_JacobianGradient_FluxVectors_3d( vF, vFbar, gp, tStep );
  } else if ( mode == _PlaneStrain ) {
    return this->give_Fbar_FluxVectors_PlaneStrain( vF, vFbar, gp, tStep );
  } else {
    OOFEM_ERROR( "unsupported mode" );
  }
}


std::tuple<FloatArrayF<9>, FloatArrayF<9> >
ThirdMediumCrossSection ::give_Fbar_FluxVectors_3d( const FloatArrayF<9> &vF, const FloatArrayF<9> &vFbar, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_Fbar_FirstPKStressVector_FPKbarStressVector_3d( vF, vFbar, gp, tStep );
}


std::tuple<FloatArrayF<5>, FloatArrayF<5> >
ThirdMediumCrossSection ::give_Fbar_FluxVectors_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<5> &vFbar, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_Fbar_FirstPKStressVector_FPKbarStressVector_PlaneStrain( vF, vFbar, gp, tStep );
}


std::tuple<FloatMatrix, FloatMatrix, FloatMatrix, FloatMatrix> ThirdMediumCrossSection ::give_Fbar_dFluxes_dGrads( MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep )
{
  MaterialMode mode = gp->giveMaterialMode();
  if ( mode == _3dMat ) {
    return this->give_Fbar_ConstitutiveMatrices_3d( rmode, gp, tStep );
  } else if ( mode == _PlaneStrain ) {
    return this->give_Fbar_ConstitutiveMatrices_PlaneStrain( rmode, gp, tStep );
  } else {
    OOFEM_ERROR( "unsupported mode" );
  }
}


std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9> >
ThirdMediumCrossSection ::give_Fbar_ConstitutiveMatrices_3d( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_Fbar_ConstitutiveMatrices_3d( rMode, gp, tStep );
}


std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 5>, FloatMatrixF<5, 5>, FloatMatrixF<5, 5> >
ThirdMediumCrossSection ::give_Fbar_ConstitutiveMatrices_PlaneStrain( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const
{
  auto mat = dynamic_cast<ThirdMediumMaterial *>( this->giveMaterial( gp ) );
  return mat->give_Fbar_ConstitutiveMatrices_PlaneStrain( rMode, gp, tStep );
}


} // namespace oofem