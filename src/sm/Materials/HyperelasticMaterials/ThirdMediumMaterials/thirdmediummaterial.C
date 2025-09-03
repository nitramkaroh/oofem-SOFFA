
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
#include "gausspoint.h"
#include "element.h"

namespace oofem {

  void ThirdMediumMaterial::initializeFrom(InputRecord& ir)
  {
    IR_GIVE_OPTIONAL_FIELD( ir, this->kappaGradFGradF, _IFT_ThirdMediumMaterial_kappaGradFGradF );
    IR_GIVE_OPTIONAL_FIELD( ir, this->kappaGradJGradJ, _IFT_ThirdMediumMaterial_kappaGradJGradJ );
    IR_GIVE_OPTIONAL_FIELD( ir, this->kappaGradRGradR, _IFT_ThirdMediumMaterial_kappaGradRGradR );
    IR_GIVE_OPTIONAL_FIELD( ir, this->kappaJbar, _IFT_ThirdMediumMaterial_kappaJbar );
    IR_GIVE_OPTIONAL_FIELD( ir, this->kappaJbarVol, _IFT_ThirdMediumMaterial_kappaJbarVol );
    IR_GIVE_OPTIONAL_FIELD( ir, this->kappaFbar, _IFT_ThirdMediumMaterial_kappaFbar );
    
    IR_GIVE_OPTIONAL_FIELD( ir, this->E, _IFT_ThirdMediumMaterial_eLinear );
    IR_GIVE_OPTIONAL_FIELD( ir, this->nu, _IFT_ThirdMediumMaterial_nuLinear );
    this->initTangents();

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

 ////////////////////////FIRST AND SECOND GRADIENTS TOGETHER////////////////////////////

std::tuple<FloatArrayF<5>, FloatArrayF<8> >
ThirdMediumMaterial::give_FirstSecondGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep )
{
  FloatArrayF<5> vP_ans;
  FloatArrayF<8> vT_ans;
  
  if (kappaGradJGradJ != 0.0){
    auto [vP_gradJ, vT_gradJ] = this->give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain(vF,vGradF,gp,tStep);
    vP_ans += vP_gradJ;
    vT_ans += vT_gradJ;
  }
  if (kappaGradRGradR != 0.0){
    auto [vP_gradR, vT_gradR] = this->give_RotationGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( vF, vGradF, gp, tStep );
    vP_ans += vP_gradR;
    vT_ans += vT_gradR;
  }

  return std::make_tuple(vP_ans, vT_ans);
  
}

std::tuple<FloatArrayF<9>, FloatArrayF<27> >
ThirdMediumMaterial::give_FirstSecondGradient_FirstPKStressVector_SecondOrderStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vGradF, GaussPoint *gp, TimeStep *tStep )
{
  FloatArrayF<9> vP_ans;
  FloatArrayF<27> vT_ans;

  if ( kappaGradJGradJ != 0.0 ) {
    auto [vP_gradJ, vT_gradJ] = this->give_JacobianGradient_FirstPKStressVector_SecondOrderStressVector_3d( vF, vGradF, gp, tStep );
    vP_ans += vP_gradJ;
    vT_ans += vT_gradJ;
  }
  if (kappaGradRGradR != 0.0 ){
    OOFEM_ERROR("gradRgradR regularization unsupported outside of Plane Strain formulations.");
  }

  return std::make_tuple( vP_ans, vT_ans );
}

std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 27>, FloatMatrixF<27, 9>, FloatMatrixF<27, 27> >
ThirdMediumMaterial::give_FirstSecondGradient_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  FloatMatrixF<9, 9> vdPdF_ans;
  FloatMatrixF<9, 27> vdPdGradF_ans;
  FloatMatrixF<27, 9> vdTdF_ans;
  FloatMatrixF<27, 27> vdTdGradF_ans;

  if ( kappaGradJGradJ != 0.0 ) {
    auto [vdPdF_gradJ, vdPdGradF_gradJ, vdTdF_gradJ, vdTdGradF_gradJ] = this->give_JacobianGradient_ConstitutiveMatrices_3d( mode, gp, tStep );
    vdPdF_ans += vdPdF_gradJ;
    vdPdGradF_ans += vdPdGradF_gradJ;
    vdTdF_ans += vdTdF_gradJ;
    vdTdGradF_ans += vdTdGradF_gradJ;
  }
  if ( kappaGradRGradR != 0.0 ) {
    OOFEM_ERROR( "gradRgradR regularization unsupported outside of Plane Strain formulations." );
  }

  return std::make_tuple( vdPdF_ans, vdPdGradF_ans, vdTdF_ans, vdTdGradF_ans );
}

std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> >
ThirdMediumMaterial::give_FirstSecondGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  FloatMatrixF<5, 5> vdPdF_ans;
  FloatMatrixF<5, 8> vdPdGradF_ans;
  FloatMatrixF<8, 5> vdTdF_ans;
  FloatMatrixF<8, 8> vdTdGradF_ans;

  if ( kappaGradJGradJ != 0.0 ) {
    auto [vdPdF_gradJ, vdPdGradF_gradJ, vdTdF_gradJ, vdTdGradF_gradJ] = this->give_JacobianGradient_ConstitutiveMatrices_PlaneStrain( mode, gp, tStep );
    vdPdF_ans += vdPdF_gradJ;
    vdPdGradF_ans += vdPdGradF_gradJ;
    vdTdF_ans += vdTdF_gradJ;
    vdTdGradF_ans += vdTdGradF_gradJ;
  }
  if ( kappaGradRGradR != 0.0 ) {
    auto [vdPdF_gradR, vdPdGradF_gradR, vdTdF_gradR, vdTdGradF_gradR] = this->give_RotationGradient_ConstitutiveMatrices_PlaneStrain( mode, gp, tStep );
    vdPdF_ans += vdPdF_gradR;
    vdPdGradF_ans += vdPdGradF_gradR;
    vdTdF_ans += vdTdF_gradR;
    vdTdGradF_ans += vdTdGradF_gradR;
  }

  return std::make_tuple( vdPdF_ans, vdPdGradF_ans, vdTdF_ans, vdTdGradF_ans );
}



////////////////////////JACOBIAN GRADIENT///////////////////////////////////////////////

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
  Tensor2_3d F( vF ), Fbar( vFbar );

  Tensor2_3d P_Jbar, Pbar_Jbar;
  Tensor2_3d Pbar_JbarVol;
  Tensor2_3d P_Fbar, Pbar_Fbar;

  auto [ J, cofF ] = F.compute_determinant_and_cofactor();
  auto [ Jbar, cofFbar ] = Fbar.compute_determinant_and_cofactor();

  /*if ( J < 0 || std::isnan(J)) {
    OOFEM_ERROR("J =%e in the element %i and in the Gauss point %i.", J, gp->giveElement()->giveNumber(), gp->giveNumber());
  }

  if ( Jbar < 0 || std::isnan(Jbar)) {
    OOFEM_ERROR("Jbar =%e in the element %i and in the Gauss point %i.", Jbar, gp->giveElement()->giveNumber(), gp->giveNumber());
  }*/

  P_Jbar(r_3, s_3) = kappaJbar * (J - Jbar) * cofF(r_3, s_3);
  Pbar_Jbar(r_3, s_3) = - kappaJbar * (J - Jbar) * cofFbar(r_3, s_3);

  Pbar_JbarVol(r_3, s_3) = (( kappaJbarVol * log(Jbar)) / Jbar )* cofFbar(r_3, s_3);

  P_Fbar(r_3, s_3) = kappaFbar * (F(r_3, s_3) - Fbar(r_3, s_3));
  Pbar_Fbar(r_3, s_3) = - kappaFbar * (F(r_3, s_3) - Fbar(r_3, s_3));

  // save gradients
  ThirdMediumMaterialStatus *status = this->giveStatus( gp );
  status->letTempFVectorBe( vF );
  status->letTempFbarVectorBe( vFbar );

  return std::make_tuple( P_Jbar.to_voigt_form() + P_Fbar.to_voigt_form(), Pbar_Jbar.to_voigt_form() + Pbar_JbarVol.to_voigt_form() + Pbar_Fbar.to_voigt_form() );
}

std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9>, FloatMatrixF<9, 9> >
ThirdMediumMaterial ::give_Fbar_ConstitutiveMatrices_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  // load F and GradF
  ThirdMediumMaterialStatus *status = this->giveStatus( gp );
  FloatArrayF<9> vF = status->giveTempFVector();
  FloatArrayF<9> vFbar = status->giveTempFbarVector();

  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. ), F( vF ), Fbar( vFbar );

  Tensor4_3d dPdF_Jbar, dPdFbar_Jbar, dPbardF_Jbar, dPbardFbar_Jbar;
  Tensor4_3d dPbardFbar_JbarVol;
  Tensor4_3d dPdF_Fbar, dPdFbar_Fbar, dPbardF_Fbar, dPbardFbar_Fbar;

  auto [ J, cofF ] = F.compute_determinant_and_cofactor();
  auto [ Jbar, cofFbar ] = Fbar.compute_determinant_and_cofactor();
  auto Fcross = F.compute_tensor_cross_product();
  auto Fbarcross = Fbar.compute_tensor_cross_product();

  dPdF_Jbar(r_3, s_3, u_3, v_3) = kappaJbar * (cofF(u_3, v_3) * cofF(r_3, s_3) + (J - Jbar) * Fcross(r_3, s_3, u_3, v_3));
  dPdFbar_Jbar(r_3, s_3, u_3, v_3) = - kappaJbar * cofFbar(u_3, v_3) * cofF(r_3, s_3);
  dPbardF_Jbar(r_3, s_3, u_3, v_3) = - kappaJbar * cofF(u_3, v_3) * cofFbar(r_3, s_3);
  dPbardFbar_Jbar(r_3, s_3, u_3, v_3) = kappaJbar * (cofFbar(u_3, v_3) * cofFbar(r_3, s_3) - (J - Jbar) * Fbarcross(r_3, s_3, u_3, v_3));

  dPbardFbar_JbarVol(r_3, s_3, u_3, v_3) = (kappaJbarVol/Jbar) * (log(Jbar) * Fbarcross(r_3, s_3, u_3, v_3) + ((1 - log(Jbar)) / Jbar) * cofFbar(r_3, s_3) * cofFbar(u_3, v_3));

  dPdF_Fbar(r_3, s_3, u_3, v_3) = kappaFbar * delta(r_3, u_3) * delta(s_3, v_3);
  dPdFbar_Fbar(r_3, s_3, u_3, v_3) = - kappaFbar * delta(r_3, u_3) * delta(s_3, v_3);
  dPbardF_Fbar(r_3, s_3, u_3, v_3) = - kappaFbar * delta(r_3, u_3) * delta(s_3, v_3);
  dPbardFbar_Fbar(r_3, s_3, u_3, v_3) = kappaFbar * delta(r_3, u_3) * delta(s_3, v_3);

  return std::make_tuple( dPdF_Jbar.to_voigt_form() + dPdF_Fbar.to_voigt_form(), dPdFbar_Jbar.to_voigt_form() + dPdFbar_Fbar.to_voigt_form(), dPbardF_Jbar.to_voigt_form() + dPbardF_Fbar.to_voigt_form(), dPbardFbar_Jbar.to_voigt_form() + dPbardFbar_JbarVol.to_voigt_form() + dPbardFbar_Fbar.to_voigt_form() );
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

////////////////////////////////// ROTATION GRADIENTS ///////////////////////////////////////////////////////////////////////////////////////

std::tuple<FloatArrayF<5>, FloatArrayF<8> >
ThirdMediumMaterial ::give_RotationGradient_FirstPKStressVector_SecondOrderStressVector_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep )
{
  // get tensors
  FloatArrayF<9> vF_3d = assemble<9>( vF, { 0, 1, 2, 5, 8 } );
  FloatArrayF<27> vGradF_3d = assemble<27>( vGradF, { 0, 1, 5, 8, 9, 10, 14, 17 } );
  Tensor2_3d F( vF_3d ), P;
  Tensor3_3d gradF( vGradF_3d ), T;

  auto gradPhi = compute_gradPhi_PlaneStrain( F, gradF );
  auto [dGradPhidF, dGradPhidGradF] = compute_gradPhi_derivatives_PlaneStrain( F, gradF );

  P( i_3, j_3 ) = kappaGradRGradR * 2 * gradPhi( m_3 ) * dGradPhidF( m_3, i_3, j_3 );
  T( i_3, j_3, k_3 ) = kappaGradRGradR * 2 * gradPhi( m_3 ) * dGradPhidGradF( m_3, i_3, j_3, k_3 );

  // save gradients
  ThirdMediumMaterialStatus *status = static_cast<ThirdMediumMaterialStatus *>( this->giveStatus( gp ) );
  status->letTempFVectorBe( vF_3d );
  status->letTempGradFVectorBe( vGradF_3d );

  // output
  auto vP_3d = P.to_voigt_form();
  auto vT_3d = T.to_voigt_form_27();
  return std::make_tuple( vP_3d[{ 0, 1, 2, 5, 8 }], vT_3d[{ 0, 1, 5, 8, 9, 10, 14, 17 }] );
}


std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> >
ThirdMediumMaterial ::give_RotationGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  // load F and GradF
  ThirdMediumMaterialStatus *status = static_cast<ThirdMediumMaterialStatus *>( this->giveStatus( gp ) );
  FloatArrayF<9> vF_3d = status->giveTempFVector();
  FloatArrayF<27> vGradF_3d = status->giveTempGradFVector();

  // get tensors
  Tensor2_3d F( vF_3d );
  Tensor3_3d gradF( vGradF_3d );
  Tensor4_3d D_FF;
  Tensor5_3d C_FG;
  Tensor6_3d C_GG;

  auto gradPhi = compute_gradPhi_PlaneStrain( F, gradF );
  auto [dGradPhidF, dGradPhidGradF] = compute_gradPhi_derivatives_PlaneStrain( F, gradF );
  auto [d2GradPhidFdF, d2GradPhidFdGradF, d2GradPhidGradFdGradF] = compute_gradPhi_secondDerivatives_PlaneStrain( F, gradF );

  D_FF( i_3, j_3, r_3, s_3 ) = kappaGradRGradR * 2 * dGradPhidF( m_3, r_3, s_3 ) * dGradPhidF( m_3, i_3, j_3 ) + kappaGradRGradR * 2 * gradPhi( m_3 ) * d2GradPhidFdF( m_3, i_3, j_3, r_3, s_3 );
  C_FG( i_3, j_3, r_3, s_3, t_3 ) = kappaGradRGradR * 2 * dGradPhidGradF( m_3, r_3, s_3, t_3 ) * dGradPhidF( m_3, i_3, j_3 ) + kappaGradRGradR * 2 * gradPhi( m_3 ) * d2GradPhidFdGradF( m_3, i_3, j_3, r_3, s_3, t_3 );
  C_GG( i_3, j_3, k_3, r_3, s_3, t_3 ) = kappaGradRGradR * 2 * dGradPhidGradF( m_3, r_3, s_3, t_3 ) * dGradPhidGradF( m_3, i_3, j_3, k_3 ) + kappaGradRGradR * 2 * gradPhi( m_3 ) * d2GradPhidGradFdGradF( m_3, i_3, j_3, k_3, r_3, s_3, t_3 );

  // output
  auto vD_FF_3d = D_FF.to_voigt_form();
  FloatMatrix vC_FG_3d( C_FG.to_voigt_form_9x27() ), vC_GF_3d;
  vC_GF_3d.beTranspositionOf( vC_FG_3d );
  auto vC_GG_3d = C_GG.to_voigt_form();

  auto vD_FF = vD_FF_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  auto vC_FG = ( (FloatMatrixF<9, 27>)vC_FG_3d ).operator()( { 0, 1, 2, 5, 8 }, { 0, 1, 5, 8, 9, 10, 14, 17 } );
  auto vC_GF = ( (FloatMatrixF<27, 9>)vC_GF_3d ).operator()( { 0, 1, 5, 8, 9, 10, 14, 17 }, { 0, 1, 2, 5, 8 } );
  auto vC_GG = vC_GG_3d( { 0, 1, 5, 8, 9, 10, 14, 17 }, { 0, 1, 5, 8, 9, 10, 14, 17 } );
  return std::make_tuple( vD_FF, vC_FG, vC_GF, vC_GG );
}


Tensor1_3d
ThirdMediumMaterial::compute_gradPhi_PlaneStrain( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor1_3d gradPhi;

  auto Z = compute_Z_fraction_PlaneStrain( F );
  auto dZdF = compute_Z_fraction_derivative_PlaneStrain( F );

  gradPhi( m_3 ) = ( 1. / ( Z * Z + 1 ) ) * dZdF( p_3, q_3 ) * gradF( p_3, q_3, m_3 );

  return gradPhi;
}

std::tuple<Tensor3_3d, Tensor4_3d>
ThirdMediumMaterial::compute_gradPhi_derivatives_PlaneStrain( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor3_3d dGradPhidF;
  Tensor4_3d dGradPhidGradF;

  auto Z = compute_Z_fraction_PlaneStrain( F );
  auto dZdF = compute_Z_fraction_derivative_PlaneStrain( F );
  auto d2ZdFdF = compute_Z_fraction_secondDerivative_PlaneStrain( F );

  dGradPhidF( m_3, i_3, j_3 ) = ( -2. * Z / ( Z * Z + 1. ) * ( Z * Z + 1. ) ) * dZdF( i_3, j_3 ) * dZdF( p_3, q_3 ) * gradF( p_3, q_3, m_3 ) + ( 1. / ( Z * Z + 1 ) ) * d2ZdFdF( p_3, q_3, i_3, j_3 ) * gradF( p_3, q_3, m_3 );
  dGradPhidGradF( m_3, i_3, j_3, k_3 ) = ( 1. / ( Z * Z + 1 ) ) * dZdF( i_3, j_3 ) * delta( m_3, k_3 );

  return std::make_tuple( dGradPhidF, dGradPhidGradF );
}

std::tuple<Tensor5_3d, Tensor6_3d, Tensor7_3d>
ThirdMediumMaterial::compute_gradPhi_secondDerivatives_PlaneStrain( const Tensor2_3d &F, const Tensor3_3d &gradF ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor5_3d d2GradPhidFdF;
  Tensor6_3d d2GradPhidFdGradF;
  Tensor7_3d d2GradPhidGradFdGradF;

  auto Z = compute_Z_fraction_PlaneStrain( F );
  auto dZdF = compute_Z_fraction_derivative_PlaneStrain( F );
  auto d2ZdFdF = compute_Z_fraction_secondDerivative_PlaneStrain( F );
  auto d3ZdFdFdF = compute_Z_fraction_thirdDerivative_PlaneStrain( F );

  d2GradPhidFdF( m_3, i_3, j_3, r_3, s_3 ) = ( 8. * Z * Z / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * dZdF( r_3, s_3 ) * dZdF( i_3, j_3 ) * dZdF( p_3, q_3 ) * gradF( p_3, q_3, m_3 )
      - ( 2. / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * dZdF( r_3, s_3 ) * dZdF( i_3, j_3 ) * dZdF( p_3, q_3 ) * gradF( p_3, q_3, m_3 )
      - ( 2. * Z / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * ( d2ZdFdF( i_3, j_3, r_3, s_3 ) * dZdF( p_3, q_3 ) * gradF( p_3, q_3, m_3 ) )
      - ( 2. * Z / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * dZdF( i_3, j_3 ) * d2ZdFdF( p_3, q_3, r_3, s_3 ) * gradF( p_3, q_3, m_3 )
      - ( 2. * Z / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * dZdF( r_3, s_3 ) * d2ZdFdF( p_3, q_3, i_3, j_3 ) * gradF( p_3, q_3, m_3 )
      + ( 1. / ( Z * Z + 1 ) ) * d3ZdFdFdF( p_3, q_3, i_3, j_3, r_3, s_3 ) * gradF( p_3, q_3, m_3 );

  d2GradPhidFdGradF( m_3, i_3, j_3, r_3, s_3, t_3 ) = -( 2. * Z / ( ( Z * Z + 1. ) * ( Z * Z + 1. ) ) ) * dZdF( i_3, j_3 ) * ( dZdF( r_3, s_3 ) * delta( m_3, t_3 ) )
      + ( 1. / ( Z * Z + 1 ) ) * d2ZdFdF( r_3, s_3, i_3, j_3 ) * delta( m_3, t_3 );


  d2GradPhidGradFdGradF( m_3, i_3, j_3, k_3, r_3, s_3, t_3 ) = 0.;

  return std::make_tuple( d2GradPhidFdF, d2GradPhidFdGradF, d2GradPhidGradFdGradF );
}

double
ThirdMediumMaterial::compute_Z_fraction_PlaneStrain( const Tensor2_3d &F ) const
{
  return ( F( 0, 1 ) - F( 1, 0 ) ) / ( F( 0, 0 ) + F( 1, 1 ) );
}

Tensor2_3d
ThirdMediumMaterial::compute_Z_fraction_derivative_PlaneStrain( const Tensor2_3d &F ) const
{
  Tensor2_3d dZdF, delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );

  dZdF( p_3, q_3 ) = -( F( 0, 1 ) - F( 1, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( p_3, 0 ) * delta( q_3, 0 ) + delta( p_3, 1 ) * delta( q_3, 1 ) )
      + ( delta( p_3, 0 ) * delta( q_3, 1 ) - delta( p_3, 1 ) * delta( q_3, 0 ) ) / ( F( 0, 0 ) + F( 1, 1 ) );

  return dZdF;
}

Tensor4_3d
ThirdMediumMaterial::compute_Z_fraction_secondDerivative_PlaneStrain( const Tensor2_3d &F ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor4_3d d2ZdFdF;

  d2ZdFdF( p_3, q_3, i_3, j_3 ) = -( delta( i_3, 0 ) * delta( j_3, 1 ) - delta( i_3, 1 ) * delta( j_3, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( p_3, 0 ) * delta( q_3, 0 ) + delta( p_3, 1 ) * delta( q_3, 1 ) )
      + 2. * ( F( 0, 1 ) - F( 1, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( p_3, 0 ) * delta( q_3, 0 ) + delta( p_3, 1 ) * delta( q_3, 1 ) ) * ( delta( i_3, 0 ) * delta( j_3, 0 ) + delta( i_3, 1 ) * delta( j_3, 1 ) )
      - ( delta( p_3, 0 ) * delta( q_3, 1 ) - delta( p_3, 1 ) * delta( q_3, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( i_3, 0 ) * delta( j_3, 0 ) + delta( i_3, 1 ) * delta( j_3, 1 ) );

  return d2ZdFdF;
}

Tensor6_3d
ThirdMediumMaterial::compute_Z_fraction_thirdDerivative_PlaneStrain( const Tensor2_3d &F ) const
{
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor6_3d d3ZdFdFdF;

  d3ZdFdFdF( p_3, q_3, i_3, j_3, r_3, s_3 ) = 2. * ( delta( i_3, 0 ) * delta( j_3, 1 ) - delta( i_3, 1 ) * delta( j_3, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( p_3, 0 ) * delta( q_3, 0 ) + delta( p_3, 1 ) * delta( q_3, 1 ) ) * ( delta( r_3, 0 ) * delta( s_3, 0 ) + delta( r_3, 1 ) * delta( s_3, 1 ) )
      + 2. * ( delta( r_3, 0 ) * delta( s_3, 1 ) - delta( r_3, 1 ) * delta( s_3, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( p_3, 0 ) * delta( q_3, 0 ) + delta( p_3, 1 ) * delta( q_3, 1 ) ) * ( delta( i_3, 0 ) * delta( j_3, 0 ) + delta( i_3, 1 ) * delta( j_3, 1 ) )
      - 6. * ( F( 0, 1 ) - F( 1, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( p_3, 0 ) * delta( q_3, 0 ) + delta( p_3, 1 ) * delta( q_3, 1 ) ) * ( delta( i_3, 0 ) * delta( j_3, 0 ) + delta( i_3, 1 ) * delta( j_3, 1 ) ) * ( delta( r_3, 0 ) * delta( s_3, 0 ) + delta( r_3, 1 ) * delta( s_3, 1 ) )
      + 2. * ( delta( p_3, 0 ) * delta( q_3, 1 ) - delta( p_3, 1 ) * delta( q_3, 0 ) ) / ( ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) * ( F( 0, 0 ) + F( 1, 1 ) ) ) * ( delta( r_3, 0 ) * delta( s_3, 0 ) + delta( r_3, 1 ) * delta( s_3, 1 ) ) * ( delta( i_3, 0 ) * delta( j_3, 0 ) + delta( i_3, 1 ) * delta( j_3, 1 ) );

  return d3ZdFdFdF;
}

////////////////////////////////// LINEAR ELASTICITY ///////////////////////////////////////////////////////////////////////////////////////

FloatArrayF<4>
ThirdMediumMaterial::give_LinearElasticity_RealStressVector_PlaneStrain( const FloatArrayF<4> &eps, GaussPoint *gp, TimeStep *tStep )
{
  auto vS = this->give_LinearElasticity_RealStressVector_3d( assemble<6>( eps, { 0, 1, 2, 5 } ), gp, tStep );
  return vS[{ 0, 1, 2, 5 }];
}

FloatArrayF<6>
ThirdMediumMaterial::give_LinearElasticity_RealStressVector_3d( const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep )
{
  auto status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

  auto d = this->give_LinearElasticity_ConstitutiveMatrix_3d( TangentStiffness, gp, tStep );

  FloatArrayF<6> sig;
  sig = dot( d, eps );
  // update gp
  status->letTempStrainVectorBe( eps );
  status->letTempStressVectorBe( sig );
  return sig;
}

FloatMatrixF<4,4>
ThirdMediumMaterial::give_LinearElasticity_ConstitutiveMatrix_PlaneStrain( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  return linearElasticTangentPlaneStrain;
}

FloatMatrixF<6, 6>
ThirdMediumMaterial::give_LinearElasticity_ConstitutiveMatrix_3d( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep )
{
  return linearElasticTangent;
}

void
ThirdMediumMaterial ::initTangents()
{
  double G = E / ( 2.0 * ( 1. + nu ) );
  double K = E / ( 3.0 * ( 1. - 2. * nu ) );
  this->linearElasticTangent = 2 * G * I_dev6 + K * I6_I6;
  
  linearElasticTangentPlaneStrain = {
    linearElasticTangent( 0, 0 ),
    linearElasticTangent( 1, 0 ),
    linearElasticTangent( 2, 0 ),
    linearElasticTangent( 3, 0 ),
    linearElasticTangent( 0, 1 ),
    linearElasticTangent( 1, 1 ),
    linearElasticTangent( 2, 1 ),
    linearElasticTangent( 3, 1 ),
    linearElasticTangent( 0, 2 ),
    linearElasticTangent( 1, 2 ),
    linearElasticTangent( 2, 2 ),
    linearElasticTangent( 3, 2 ),
    linearElasticTangent( 0, 3 ),
    linearElasticTangent( 1, 3 ),
    linearElasticTangent( 2, 3 ),
    linearElasticTangent( 3, 3 ),
  };
}



} // end namespace oofem

