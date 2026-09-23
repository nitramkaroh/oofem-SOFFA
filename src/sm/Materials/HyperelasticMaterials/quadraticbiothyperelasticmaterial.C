#include "quadraticbiothyperelasticmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"


namespace oofem {
REGISTER_Material( QuadraticBiotHyperelasticMaterial );

QuadraticBiotHyperelasticMaterial::QuadraticBiotHyperelasticMaterial( int n, Domain *d ) :
    StructuralMaterial( n, d ), BaseHyperElasticMaterial()
{
}

FloatArrayF<5>
QuadraticBiotHyperelasticMaterial::giveFirstPKStressVector_PlaneStrain( const FloatArrayF<5> &vF_ps, GaussPoint *gp, TimeStep *tStep ) const
// returns 5 components of the first piola kirchhoff stress corresponding to the given deformation gradient
{
  StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

  // from plane strain to 3d
  auto vF_3d = assemble<9>( vF_ps, { 0, 1, 2, 5, 8 } );

  Tensor2_3d F( vF_3d );
  auto t = this->compute_t_from_F( F );
  auto R = this->compute_R_from_F( F );

  // first Piola stress P = R t
  Tensor2_3d P;
  P( i_3, j_3 ) = R( i_3, k_3 ) * t( k_3, j_3 );

  auto vP_3d = P.to_voigt_form();
  // update gp
  status->letTempFVectorBe( vF_3d );
  status->letTempPVectorBe( vP_3d );
  //
  auto vP_red = vP_3d[ { 0, 1, 2, 5, 8 } ];
  return vP_red;
}


FloatMatrixF<5, 5>
QuadraticBiotHyperelasticMaterial::givePlaneStrainStiffnessMatrix_dPdF( MatResponseMode mode, GaussPoint *gp, TimeStep *tStep ) const
// returns the 5x5 tangent stiffness matrix - dP/dF
{
  StructuralMaterialStatus *status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );

  FloatArrayF<9> vF_3d( status->giveTempFVector() );

  Tensor2_3d F( vF_3d );
  Tensor4_3d A;
  auto t = this->compute_t_from_F( F );
  auto R = this->compute_R_from_F( F );
  auto dRdF = this->compute_dRdF_from_F( F );
  auto dtFd = this->compute_dtdF_from_F( F );

  // product rule for P = R t, A = (dR/dF) t + R (dt/dF)
  A( i_3, j_3, k_3, l_3 ) = dRdF( i_3, p_3, k_3, l_3 ) * t( p_3, j_3 ) + R( i_3, p_3 ) * dtFd( p_3, j_3, k_3, l_3 );

  //back to plane strain
  auto vA_3d = A.to_voigt_form();
  auto vA_red = vA_3d( { 0, 1, 2, 5, 8 }, { 0, 1, 2, 5, 8 } );
  return vA_red;
}


MaterialStatus *
QuadraticBiotHyperelasticMaterial::CreateStatus( GaussPoint *gp ) const
{
  return new StructuralMaterialStatus( gp );
}


void QuadraticBiotHyperelasticMaterial::initializeFrom( InputRecord &ir )
{
  StructuralMaterial::initializeFrom( ir );
  BaseHyperElasticMaterial::initializeFrom( ir );
  IR_GIVE_FIELD( ir, mu, _IFT_QuadraticBiotHyperelasticMaterial_mu );
  IR_GIVE_FIELD( ir, lambda, _IFT_QuadraticBiotHyperelasticMaterial_lambda );
}


double QuadraticBiotHyperelasticMaterial::compute_s_from_F( const Tensor2_3d &F ) const
{
  // s = tr(U_2d) = sqrt( tr(C_2d) + 2 j );
  // in plane strain the out-of-plane stretch is 1, so tr(C_2d) = tr(C) - 1
  auto J = F.compute_determinant();
  Tensor2_3d C;
  C( i_3, j_3 ) = F( k_3, i_3 ) * F( k_3, j_3 );
  double traceC = C( i_3, i_3 );
  return std::sqrt( traceC - 1. + 2. * J );
}


Tensor2_3d QuadraticBiotHyperelasticMaterial::compute_U_from_F( const Tensor2_3d &F ) const
{
  // U_2d = (C_2d + j delta_2d) / s, embedded with U_33 = 1
  auto J = F.compute_determinant();
  double s = this->compute_s_from_F( F );
  Tensor2_3d C, U;
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  C( i_3, j_3 ) = F( k_3, i_3 ) * F( k_3, j_3 );
  U( i_3, j_3 ) = ( C( i_3, j_3 ) + J * delta( i_3, j_3 ) ) / s;
  U( 2, 2 ) = 1.;
  return U;
}


Tensor2_3d QuadraticBiotHyperelasticMaterial::compute_Uinv_from_F( const Tensor2_3d &F ) const
{
  // U_2d^{-1} = ((tr(C_2d) + j) delta_2d - C_2d) / (s j), embedded with (U^{-1})_33 = 1
  auto J = F.compute_determinant();
  double s = this->compute_s_from_F( F );
  Tensor2_3d C, Uinv;
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  C( i_3, j_3 ) = F( k_3, i_3 ) * F( k_3, j_3 );
  double traceC = C( i_3, i_3 );
  Uinv( i_3, j_3 ) = ( ( traceC - 1. + J ) * delta( i_3, j_3 ) - C( i_3, j_3 ) ) / ( s * J );
  Uinv( 2, 2 ) = 1.;
  return Uinv;
}


Tensor2_3d QuadraticBiotHyperelasticMaterial::compute_R_from_F( const Tensor2_3d &F ) const
{
  // R = F U^{-1}; the plane-strain embedding yields diag(R_2d, 1)
  auto Uinv = this->compute_Uinv_from_F( F );
  Tensor2_3d R;
  R( i_3, j_3 ) = F( i_3, k_3 ) * Uinv( k_3, j_3 );
  return R;
}


Tensor2_3d QuadraticBiotHyperelasticMaterial::compute_t_from_F( const Tensor2_3d &F ) const
{
  // t = 2 mu (U - I) + lambda tr(U - I) I,  tr(U - I) = s - 2
  // (with U_33 = 1 this is exact in all 9 components: t_33 = lambda (s - 2))
  double s = this->compute_s_from_F( F );
  auto U = this->compute_U_from_F( F );
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor2_3d t;
  t( i_3, j_3 ) = 2. * mu * ( U( i_3, j_3 ) - delta( i_3, j_3 ) ) + lambda * ( s - 2. ) * delta( i_3, j_3 );
  return t;
}


Tensor2_3d QuadraticBiotHyperelasticMaterial::compute_dsdF_from_F( const Tensor2_3d &F ) const
{
  // ds/dF_kl = (F_kl + cof_kl) / s, in-plane block
  double s = this->compute_s_from_F( F );
  auto cofF = F.compute_cofactor();
  Tensor2_3d dsdF;
  dsdF( i_3, j_3 ) = ( F( i_3, j_3 ) + cofF( i_3, j_3 ) ) / s;
  // s = tr(U_2d) depends only on the in-plane block of F
  for ( int ii = 0; ii < 3; ++ii )
    for ( int jj = 0; jj < 3; ++jj )
      if ( ii == 2 || jj == 2 )
        dsdF( ii, jj ) = 0.;
  return dsdF;
}


Tensor4_3d QuadraticBiotHyperelasticMaterial::compute_dUdF_from_F( const Tensor2_3d &F ) const
{
  // dU_ij/dF_kl = (delta_il F_kj + delta_jl F_ki + delta_ij cof_kl) / s
  //                 - (F_kl + cof_kl) U_ij / s^2
  double s = this->compute_s_from_F( F );
  auto U = this->compute_U_from_F( F );
  auto cofF = F.compute_cofactor();
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor4_3d dUdF;
  dUdF( i_3, j_3, k_3, l_3 ) = ( delta( i_3, l_3 ) * F( k_3, j_3 ) + delta( j_3, l_3 ) * F( k_3, i_3 ) + delta( i_3, j_3 ) * cofF( k_3, l_3 ) ) * ( 1. / s )
      - ( F( k_3, l_3 ) + cofF( k_3, l_3 ) ) * U( i_3, j_3 ) * ( 1. / ( s * s ) );
  // plane strain: U(F) = diag(U_2d(F_2d), 1) depends only on the in-plane
  // block of F, hence all derivatives with an out-of-plane index vanish
  for ( int ii = 0; ii < 3; ++ii )
    for ( int jj = 0; jj < 3; ++jj )
      for ( int kk = 0; kk < 3; ++kk )
        for ( int ll = 0; ll < 3; ++ll )
          if ( ii == 2 || jj == 2 || kk == 2 || ll == 2 )
            dUdF( ii, jj, kk, ll ) = 0.;
  return dUdF;
}


Tensor4_3d QuadraticBiotHyperelasticMaterial::compute_dRdF_from_F( const Tensor2_3d &F ) const
{
  // dR_ia/dF_kl = delta_ik Uinv_la - R_ic Uinv_da dU_cd/dF_kl
  auto R = this->compute_R_from_F( F );
  auto Uinv = this->compute_Uinv_from_F( F );
  auto dUdF = this->compute_dUdF_from_F( F );
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor4_3d M, term2, dRdF;
  // M(p,m,k,l) = Uinv(n,p) dUdF(m,n,k,l), summed over n
  M( p_3, m_3, k_3, l_3 ) = Uinv( n_3, p_3 ) * dUdF( m_3, n_3, k_3, l_3 );
  // term2(i,p,k,l) = R(i,m) M(p,m,k,l), summed over m
  term2( i_3, p_3, k_3, l_3 ) = M( p_3, m_3, k_3, l_3 ) * R( i_3, m_3 );
  dRdF( i_3, p_3, k_3, l_3 ) = delta( i_3, k_3 ) * Uinv( l_3, p_3 ) - term2( i_3, p_3, k_3, l_3 );
  // plane strain: R(F) = diag(R_2d(F_2d), 1) depends only on the in-plane block of F
  for ( int ii = 0; ii < 3; ++ii )
    for ( int aa = 0; aa < 3; ++aa )
      for ( int kk = 0; kk < 3; ++kk )
        for ( int ll = 0; ll < 3; ++ll )
          if ( kk == 2 || ll == 2 )
            dRdF( ii, aa, kk, ll ) = 0.;
  return dRdF;
}


Tensor4_3d QuadraticBiotHyperelasticMaterial::compute_dtdF_from_F( const Tensor2_3d &F ) const
{
  // dt_ij/dF_kl = 2 mu dU_ij/dF_kl + lambda delta_ij ds/dF_kl
  auto dUdF = this->compute_dUdF_from_F( F );
  auto dsdF = this->compute_dsdF_from_F( F );
  Tensor2_3d delta( 1., 0., 0., 0., 1., 0., 0., 0., 1. );
  Tensor4_3d dtFd;
  dtFd( i_3, j_3, k_3, l_3 ) = 2. * mu * dUdF( i_3, j_3, k_3, l_3 ) + lambda * delta( i_3, j_3 ) * dsdF( k_3, l_3 );
  return dtFd;
}


} // end namespace oofem
