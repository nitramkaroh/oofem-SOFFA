
#pragma once

#include "thirdmediummaterial.h"
#include "tensor/tensor1.h"
#include "tensor/tensor2.h"
#include "tensor/tensor3.h"
#include "tensor/tensor4.h"
#include "tensor/tensor5.h"
#include "tensor/tensor6.h"
#include "tensor/tensor7.h"
#include "thirdmediummaterialstatus.h"

namespace oofem {

  void ThirdMediumMaterial::initializeFrom(InputRecord& ir)
  {
    IR_GIVE_OPTIONAL_FIELD( ir, this->kappaGradFGradF, _IFT_ThirdMediumMaterial_kappaGradFGradF );
    IR_GIVE_OPTIONAL_FIELD( ir, this->kappaGradJGradJ, _IFT_ThirdMediumMaterial_kappaGradJGradJ );
    IR_GIVE_OPTIONAL_FIELD( ir, this->kappaJbar, _IFT_ThirdMediumMaterial_kappaJbar );
  }

  ThirdMediumMaterialStatus *ThirdMediumMaterial::giveStatus( GaussPoint *gp )
  {
    //this is a terrible hack but every ThirdMediumMaterial will always be also a descendant of Material, so it should work
    Material* material = dynamic_cast<Material *> (this);
    return dynamic_cast<ThirdMediumMaterialStatus*> (material->giveStatus(gp));
  }

  ////////////////////////////////////////////////////////////////////////////

 FloatArrayF<8>
ThirdMediumMaterial ::give_SecondGradient_FirstPKStressVector_PlaneStrain( const FloatArrayF<8> &vG, GaussPoint *gp, TimeStep *tStep )
{
  auto vT = kappaGradFGradF * vG;
  return vT;
}


FloatArrayF<27>
ThirdMediumMaterial ::give_SecondGradient_FirstPKStressVector_3d( const FloatArrayF<27> &vG, GaussPoint *gp, TimeStep *tStep )
{
  auto vT = kappaGradFGradF * vG;
  return vT;
}


FloatMatrixF<1, 1>
ThirdMediumMaterial ::give_SecondGradient_ConstitutiveMatrix_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  return kappaGradFGradF;
}

FloatMatrixF<1, 1>
ThirdMediumMaterial ::give_SecondGradient_ConstitutiveMatrix_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  return kappaGradFGradF;
}

 ////////////////////////JACOBIAN////////////////////////////

std::tuple<FloatArrayF<5>, FloatArrayF<8> >
ThirdMediumMaterial ::give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep )
{
  // these elements should give all Tensor3 components without the third dimension in them
  //  according to the established Voigt notation, see Tensor3_3d::to_voigt_form_27()
  auto [vP, vT] = give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_3d( assemble<9>( vF, { 0, 1, 2, 5, 8 } ), assemble<27>( vGradF, { 0, 1, 5, 8, 9, 10, 14, 17 } ), gp, tStep );
  //
  return std::make_tuple( vP[{ 0, 1, 2, 5, 8 }], vT[{ 0, 1, 5, 8, 9, 10, 14, 17 }] );
}


std::tuple<FloatArrayF<9>, FloatArrayF<27> >
ThirdMediumMaterial ::give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vGradF, GaussPoint *gp, TimeStep *tStep )
{
  Tensor2_3d P, F( vF );
  Tensor3_3d T, gradF( vGradF );

  auto gradJ = compute_gradJ_3d( F, gradF );
  auto [dGradJdF, dGradJdGradF] = compute_gradJ_derivatives_3d( F, gradF );

  P( r_3, s_3 ) = kappaGradJGradJ * gradJ( m_3 ) * dGradJdF( m_3, r_3, s_3 );
  T( r_3, s_3, t_3 ) = kappaGradJGradJ * gradJ( m_3 ) * dGradJdGradF( m_3, r_3, s_3, t_3 );

  // save gradients
  ThirdMediumMaterialStatus *status = this->giveStatus( gp );
  status->letTempFVectorBe( vF );
  status->letTempGradFVectorBe( vGradF );

  return std::make_tuple( P.to_voigt_form(), T.to_voigt_form_27() );
}


std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 27>, FloatMatrixF<27, 9>, FloatMatrixF<27, 27> >
ThirdMediumMaterial ::give_JacobianGradient_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  // load F and GradF
  ThirdMediumMaterialStatus *status = this->giveStatus( gp );
  FloatArrayF<9> vF = status->giveTempFVector();
  FloatArrayF<27> vGradF = status->giveTempGradFVector();

  Tensor2_3d F( vF );
  Tensor3_3d gradF( vGradF );

  Tensor4_3d D_FF;
  Tensor5_3d C_FG;
  Tensor6_3d C_GG;

  auto gradJ = compute_gradJ_3d( F, gradF );
  auto [dGradJdF, dGradJdGradF] = compute_gradJ_derivatives_3d( F, gradF );
  auto [d2GradJdFdF, d2GradJdFdGradF, d2GradJdGradFdF, d2GradJdGradFdGradF] = compute_gradJ_secondDerivatives_3d( F, gradF );

  D_FF( r_3, s_3, u_3, v_3 ) = kappaGradJGradJ * dGradJdF( m_3, u_3, v_3 ) * dGradJdF( m_3, r_3, s_3 ) + kappaGradJGradJ * gradJ( m_3 ) * d2GradJdFdF( m_3, r_3, s_3, u_3, v_3 );
  C_FG( r_3, s_3, u_3, v_3, w_3 ) = kappaGradJGradJ * dGradJdGradF( m_3, u_3, v_3, w_3 ) * dGradJdF( m_3, r_3, s_3 ) + kappaGradJGradJ * gradJ( m_3 ) * d2GradJdFdGradF( m_3, r_3, s_3, u_3, v_3, w_3 );
  C_GG( r_3, s_3, t_3, u_3, v_3, w_3 ) = kappaGradJGradJ * dGradJdGradF( m_3, u_3, v_3, w_3 ) * dGradJdGradF( m_3, r_3, s_3, t_3 ) + kappaGradJGradJ * gradJ( m_3 ) * d2GradJdGradFdGradF( m_3, r_3, s_3, t_3, u_3, v_3, w_3 );

  FloatMatrix vC_FG( C_FG.to_voigt_form_9x27() ), vC_GF;
  vC_GF.beTranspositionOf( vC_FG );

  return std::make_tuple( D_FF.to_voigt_form(), vC_FG, vC_GF, C_GG.to_voigt_form() );
}

std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> >
ThirdMediumMaterial ::give_JacobianGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  //these elements should give all Tensor3 components without the third dimension in them
  // according to the established Voigt notation, see Tensor3_3d::to_voigt_form_27()
  auto [vD_FF_3d, vC_FG_3d, vC_GF_3d, vC_GG_3d] = this->give_JacobianGradient_ConstitutiveMatrices_3d(mode, gp, tStep);
  auto vD_FF_red = vD_FF_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  auto vC_FG_red = vC_FG_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 5, 8, 9, 10, 14, 17 } );
  auto vC_GF_red = vC_GF_3d( { 0, 1, 5, 8, 9, 10, 14, 17 }, { 0, 1, 2, 5, 8 } );
  auto vC_GG_red = vC_GG_3d( { 0, 1, 5, 8, 9, 10, 14, 17 }, { 0, 1, 5, 8, 9, 10, 14, 17 } );
  return std::make_tuple( vD_FF_red, vC_FG_red, vC_GF_red, vC_GG_red );
}

Tensor1_3d ThirdMediumMaterial::compute_gradJ_3d( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor1_3d gradJ;
  Tensor3_3d eps;
  eps.be_Levi_Civita();

  gradJ( m_3 ) = 0.5 * eps( i_3, k_3, p_3 ) * ( eps( j_3, l_3, q_3 ) * F( k_3, l_3 ) ) * F( p_3, q_3 ) * gradF( i_3, j_3, m_3 );

  return gradJ;
}

std::tuple<Tensor3_3d, Tensor4_3d> ThirdMediumMaterial::compute_gradJ_derivatives_3d( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor3_3d eps, dGradJdF;
  Tensor4_3d dGradJdGradF;
  eps.be_Levi_Civita();

  dGradJdF( m_3, r_3, s_3 ) = eps( i_3, r_3, p_3 ) * ( eps( j_3, s_3, q_3 ) * F( p_3, q_3 ) ) * gradF( i_3, j_3, m_3 );
  dGradJdGradF( m_3, r_3, s_3, t_3 ) = 0.5 * eps( r_3, k_3, p_3 ) * ( eps( s_3, l_3, q_3 ) * F( k_3, l_3 ) ) * ( F( p_3, q_3 ) * delta( m_3, t_3 ) );

  return std::make_tuple( dGradJdF, dGradJdGradF );
}

std::tuple<Tensor5_3d, Tensor6_3d, Tensor6_3d, Tensor7_3d> ThirdMediumMaterial::compute_gradJ_secondDerivatives_3d( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor3_3d eps;
  Tensor5_3d d2GradJdFdF;
  Tensor6_3d d2GradJdFdGradF, d2GradJdGradFdF;
  Tensor7_3d d2GradJdGradFdGradF;
  eps.be_Levi_Civita();

  d2GradJdFdF( m_3, r_3, s_3, u_3, v_3 ) = eps( i_3, r_3, u_3 ) * ( eps( j_3, s_3, v_3 ) * gradF( i_3, j_3, m_3 ) );
  d2GradJdFdGradF( m_3, r_3, s_3, u_3, v_3, w_3 ) = eps( u_3, r_3, p_3 ) * ( eps( v_3, s_3, q_3 ) * ( F( p_3, q_3 ) * delta( m_3, w_3 ) ) );
  d2GradJdGradFdF( m_3, r_3, s_3, t_3, u_3, v_3 ) = eps( r_3, u_3, p_3 ) * ( eps( s_3, v_3, q_3 ) * ( F( p_3, q_3 ) * delta( m_3, t_3 ) ) );
  d2GradJdGradFdGradF( m_3, r_3, s_3, t_3, u_3, v_3, w_3 ) = 0.;

  return std::make_tuple( d2GradJdFdF, d2GradJdFdGradF, d2GradJdGradFdF, d2GradJdGradFdGradF );
}

 ////////////////////////FBAR////////////////////////////
 
std::tuple<FloatArrayF<5>, FloatArrayF<5> >
ThirdMediumMaterial ::give_Fbar_FirstPKStressVector_FPKbarStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<5> &vFbar, GaussPoint *gp, TimeStep *tStep )
{
  // these elements should give all Tensor3 components without the third dimension in them
  //  according to the established Voigt notation, see Tensor3_3d::to_voigt_form_27()
  auto [vP, vPbar] = give_Fbar_FirstPKStressVector_FPKbarStressVector_3d( assemble<9>( vF, { 0, 1, 2, 5, 8 } ), assemble<9>( vFbar, { 0, 1, 2, 5, 8 } ), gp, tStep );
  //
  return std::make_tuple( vP[{ 0, 1, 2, 5, 8 }], vPbar[{ 0, 1, 2, 5, 8 }] );
}


std::tuple<FloatArrayF<9>, FloatArrayF<9> >
ThirdMediumMaterial ::give_Fbar_FirstPKStressVector_FPKbarStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<9> &vFbar, GaussPoint *gp, TimeStep *tStep )
{
  Tensor2_3d P, Pbar, F( vF ), Fbar( vFbar );

  auto [ J, cofF ] = F.compute_determinant_and_cofactor();
  auto [ Jbar, cofFbar ] = Fbar.compute_determinant_and_cofactor();

  P( r_3, s_3 ) = kappaJbar * (J-Jbar)*cofF(r_3, s_3);
  Pbar( r_3, s_3 ) = -kappaJbar * (J-Jbar)*cofFbar(r_3, s_3);

  // save gradients
  ThirdMediumMaterialStatus *status = this->giveStatus( gp );
  status->letTempFVectorBe( vF );
  status->letTempFbarVectorBe( vFbar );

  return std::make_tuple( P.to_voigt_form(), Pbar.to_voigt_form() );
}


std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9> >
ThirdMediumMaterial ::give_Fbar_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  // load F and GradF
  ThirdMediumMaterialStatus *status = this->giveStatus( gp );
  FloatArrayF<9> vF = status->giveTempFVector();
  FloatArrayF<9> vFbar = status->giveTempFbarVector();

  Tensor2_3d F( vF ), Fbar( vFbar );

  Tensor4_3d dPdF, dPdFbar, dPbardF, dPbardFbar;

  auto [ J, cofF ] = F.compute_determinant_and_cofactor();
  auto [ Jbar, cofFbar ] = Fbar.compute_determinant_and_cofactor();
  auto Fcross = F.compute_tensor_cross_product();
  auto Fbarcross = Fbar.compute_tensor_cross_product();

  dPdF(r_3, s_3, u_3, v_3) = kappaJbar * cofF(u_3, v_3) * cofF(r_3, s_3) + kappaJbar * J * Fcross(r_3, s_3, u_3, v_3) - kappaJbar * Jbar * Fcross(r_3, s_3, u_3, v_3);
  dPdFbar(r_3, s_3, u_3, v_3) = - kappaJbar * cofFbar(u_3, v_3) * cofF(r_3, s_3);
  dPbardF(r_3, s_3, u_3, v_3) = - kappaJbar * cofF(u_3, v_3) * cofFbar(r_3, s_3);
  dPbardFbar(r_3, s_3, u_3, v_3) = - kappaJbar * J * Fbarcross(r_3, s_3, u_3, v_3) + kappaJbar * cofFbar(u_3, v_3) * cofFbar(r_3, s_3) + kappaJbar * Jbar * Fbarcross(r_3, s_3, u_3, v_3);

  return std::make_tuple( dPdF.to_voigt_form(), dPdFbar.to_voigt_form(), dPbardF.to_voigt_form(), dPbardFbar.to_voigt_form() );
}

std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 5>, FloatMatrixF<5, 5>, FloatMatrixF<5, 5> >
ThirdMediumMaterial ::give_Fbar_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  //these elements should give all Tensor3 components without the third dimension in them
  // according to the established Voigt notation, see Tensor3_3d::to_voigt_form_27()
  auto [vdPdF_3d, vdPdFbar_3d, vdPbardF_3d, vdPbardFbar_3d] = this->give_Fbar_ConstitutiveMatrices_3d(mode, gp, tStep);
  auto vdPdF_red = vdPdF_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  auto vdPdFbar_red = vdPdFbar_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  auto vdPbardF_red = vdPbardF_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  auto vdPbardFbar_red = vdPbardFbar_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  return std::make_tuple( vdPdF_red, vdPdFbar_red, vdPbardF_red, vdPbardFbar_red );
}



} // end namespace oofem

