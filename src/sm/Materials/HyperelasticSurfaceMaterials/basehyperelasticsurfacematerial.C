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

#include "basehyperelasticsurfacematerial.h"
#include "dynamicinputrecord.h"

namespace oofem {

double BaseHyperElasticSurfaceMaterial::compute_surface_determinant( const Tensor2_3d &F ) const
{

    // double p1 = std::pow( F( 0, 0 ) * F( 1, 1 ) - F( 0, 1 ) * F( 1, 0 ), 2 );
    // double p2 = std::pow( F( 0, 0 ) * F( 2, 1 ) - F( 0, 1 ) * F( 2, 0 ), 2 );
    // double p3 = std::pow( F( 1, 0 ) * F( 2, 1 ) - F( 1, 1 ) * F( 2, 0 ), 2 );

    // double J = std::sqrt( p1 + p2 + p3 );

    // Alternatively (holds also for prestrained reference configuration)
    Tensor2_3d C, CtC;
    C( i_3, j_3 )   = F( k_3, i_3 ) * F( k_3, j_3 );
    CtC( i_3, j_3 ) = C( k_3, i_3 ) * C( k_3, j_3 );
    double trC      = C( i_3, i_3 );
    double trCtC    = CtC( i_3, i_3 );
    double J2       = std::sqrt( 0.5 * ( trC * trC - trCtC ) );

    return J2;
};

Tensor1_3d BaseHyperElasticSurfaceMaterial::compute_surface_deformed_normal( const Tensor2_3d &F, const Tensor1_3d &Nud ) const
{
    Tensor2_3d cofF = F.compute_cofactor();
    Tensor1_3d n;
    n( i_3 ) = cofF( i_3, j_3 ) * Nud( j_3 );
    n.normalize();
    return n;
};

Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_pseudoinverse( const Tensor2_3d &F, const Tensor2_3d &I ) const
{
    Tensor2_3d JFinv = this->compute_surface_cofactor( F, I );
    Tensor2_3d Finv;
    Finv( i_3, k_3 ) = JFinv( k_3, i_3 ) / this->compute_surface_determinant( F);
    return Finv;
};


Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_cofactor( const Tensor2_3d &F ) const
{
    // This implementation holds only if the reference configuration is not strained
    Tensor2_3d cofF;
    double F11, F12, F21, F22, F31, F32;
    F11 = F( 0, 0 );
    F12 = F( 0, 1 );
    F21 = F( 1, 0 );
    F22 = F( 1, 1 );
    F31 = F( 2, 0 );
    F32 = F( 2, 1 );

    cofF( 0, 0 ) = F11 * ( pow( F12, 2 ) + pow( F22, 2 ) + pow( F32, 2 ) ) - F12 * ( F11 * F12 + F21 * F22 + F31 * F32 );
    cofF( 0, 1 ) = F12 * ( pow( F11, 2 ) + pow( F21, 2 ) + pow( F31, 2 ) ) - F11 * ( F11 * F12 + F21 * F22 + F31 * F32 );
    cofF( 1, 0 ) = F21 * ( pow( F12, 2 ) + pow( F22, 2 ) + pow( F32, 2 ) ) - F22 * ( F11 * F12 + F21 * F22 + F31 * F32 );
    cofF( 1, 1 ) = F22 * ( pow( F11, 2 ) + pow( F21, 2 ) + pow( F31, 2 ) ) - F21 * ( F11 * F12 + F21 * F22 + F31 * F32 );
    cofF( 2, 0 ) = F31 * ( pow( F12, 2 ) + pow( F22, 2 ) + pow( F32, 2 ) ) - F32 * ( F11 * F12 + F21 * F22 + F31 * F32 );
    cofF( 2, 1 ) = F32 * ( pow( F11, 2 ) + pow( F21, 2 ) + pow( F31, 2 ) ) - F31 * ( F11 * F12 + F21 * F22 + F31 * F32 );

    /////
    // Alternatively
    Tensor2_3d C, cofF2;
    C( i_3, j_3 )         = F( k_3, i_3 ) * F( k_3, j_3 );
    FloatArrayF<9> indent = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d I( indent );
    auto trC          = C( i_3, i_3 );
    cofF2( i_3, k_3 ) = F( i_3, j_3 ) * ( trC * I( j_3, k_3 ) - C( j_3, k_3 ) );
    /////

    cofF( i_3, j_3 ) = 1. / this->compute_surface_determinant( F ) * cofF( i_3, j_3 );
    // cofF( i_3, j_3 ) = 1. / this->compute_surface_determinant( F ) * cofF2( i_3, j_3 );

    return cofF;
};

Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_cofactor( const Tensor2_3d &F, const Tensor2_3d &I ) const
{
    Tensor2_3d cofF, C, cofF2;
    C( i_3, j_3 )     = F( k_3, i_3 ) * F( k_3, j_3 );
    auto trC          = C( i_3, i_3 );
    cofF2( i_3, k_3 ) = F( i_3, j_3 ) * ( trC * I( j_3, k_3 ) - C( j_3, k_3 ) );
    cofF( i_3, j_3 )  = 1. / this->compute_surface_determinant( F ) * cofF2( i_3, j_3 );

    return cofF;
};

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_dCof_dF( const Tensor2_3d &F ) const
{
    // This implementation holds only if the reference configuration is not strained
    Tensor4_3d A1, A2, A;
    double F11, F12, F21, F22, F31, F32;
    F11 = F( 0, 0 );
    F12 = F( 0, 1 );
    F21 = F( 1, 0 );
    F22 = F( 1, 1 );
    F31 = F( 2, 0 );
    F32 = F( 2, 1 );

    A1( 0, 0, 0, 0 ) = pow( F22, 2 ) + pow( F32, 2 );
    A1( 0, 1, 0, 0 ) = -F21 * F22 - F31 * F32;
    A1( 1, 0, 0, 0 ) = -F12 * F22;
    A1( 1, 1, 0, 0 ) = 2 * F11 * F22 - F12 * F21;
    A1( 2, 0, 0, 0 ) = -F12 * F32;
    A1( 2, 1, 0, 0 ) = 2 * F11 * F32 - F12 * F31;

    A1( 0, 0, 1, 0 ) = -F12 * F22;
    A1( 0, 1, 1, 0 ) = 2 * F12 * F21 - F11 * F22;
    A1( 1, 0, 1, 0 ) = pow( F12, 2 ) + pow( F32, 2 );
    A1( 1, 1, 1, 0 ) = -F11 * F12 - F31 * F32;
    A1( 2, 0, 1, 0 ) = -F22 * F32;
    A1( 2, 1, 1, 0 ) = 2 * F21 * F32 - F22 * F31;

    A1( 0, 0, 2, 0 ) = -F12 * F32;
    A1( 0, 1, 2, 0 ) = 2 * F12 * F31 - F11 * F32;
    A1( 1, 0, 2, 0 ) = -F22 * F32;
    A1( 1, 1, 2, 0 ) = 2 * F22 * F31 - F21 * F32;
    A1( 2, 0, 2, 0 ) = pow( F12, 2 ) + pow( F22, 2 );
    A1( 2, 1, 2, 0 ) = -F11 * F12 - F21 * F22;

    A1( 0, 0, 0, 1 ) = -F21 * F22 - F31 * F32;
    A1( 0, 1, 0, 1 ) = pow( F21, 2 ) + pow( F31, 2 );
    A1( 1, 0, 0, 1 ) = 2 * F12 * F21 - F11 * F22;
    A1( 1, 1, 0, 1 ) = -F11 * F21;
    A1( 2, 0, 0, 1 ) = 2 * F12 * F31 - F11 * F32;
    A1( 2, 1, 0, 1 ) = -F11 * F31;

    A1( 0, 0, 1, 1 ) = 2 * F11 * F22 - F12 * F21;
    A1( 0, 1, 1, 1 ) = -F11 * F21;
    A1( 1, 0, 1, 1 ) = -F11 * F12 - F31 * F32;
    A1( 1, 1, 1, 1 ) = pow( F11, 2 ) + pow( F31, 2 );
    A1( 2, 0, 1, 1 ) = 2 * F22 * F31 - F21 * F32;
    A1( 2, 1, 1, 1 ) = -F21 * F31;

    A1( 0, 0, 2, 1 ) = 2 * F11 * F32 - F12 * F31;
    A1( 0, 1, 2, 1 ) = -F11 * F31;
    A1( 1, 0, 2, 1 ) = 2 * F21 * F32 - F22 * F31;
    A1( 1, 1, 2, 1 ) = -F21 * F31;
    A1( 2, 0, 2, 1 ) = -F11 * F12 - F21 * F22;
    A1( 2, 1, 2, 1 ) = pow( F11, 2 ) + pow( F21, 2 );

    Tensor2_3d cofF          = this->compute_surface_cofactor( F );
    A2( i_3, j_3, k_3, l_3 ) = -cofF( i_3, j_3 ) * cofF( k_3, l_3 );
    A( i_3, j_3, k_3, l_3 )  = 1. / this->compute_surface_determinant( F ) * ( A1( i_3, j_3, k_3, l_3 ) + A2( i_3, j_3, k_3, l_3 ) );

    //// Alternatively
    // Tensor4_3d A3, A4, A5, A6, A_alt;
    // Tensor2_3d C, help, B;
    // C( i_3, j_3 )  = F( k_3, i_3 ) * F( k_3, j_3 );
    // B( i_3, j_3 )  = F( i_3, k_3 ) * F( j_3, k_3 );
    // FloatArrayF<9> indent_s = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    // FloatArrayF<9> indent = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    // Tensor2_3d Is( indent_s ), I( indent );
    // auto trC = C( i_3, i_3 );
    ////help( j_3, k_3 ) = trC * Is( j_3, k_3 ) - C( j_3, k_3 );

    // A3( i_3, j_3, k_3, l_3 ) = trC*I( i_3, k_3 ) * Is( j_3, l_3 ) - I( i_3, k_3 ) * C( j_3, l_3 );

    // A4( i_3, j_3, k_3, l_3 ) = 2 * F( i_3, j_3 ) * F( k_3, l_3 );

    // A5( i_3, j_3, k_3, l_3 ) = -F( i_3, l_3 ) * F( k_3, j_3 );

    // A6( i_3, j_3, k_3, l_3 ) = -B( i_3, k_3 ) * Is( j_3, l_3 );

    // A_alt( i_3, j_3, k_3, l_3 ) = 1. / this->compute_surface_determinant( F ) * ( A2( i_3, j_3, k_3, l_3 ) + A3( i_3, j_3, k_3, l_3 )
    //                                     + A4( i_3, j_3, k_3, l_3 ) + A5( i_3, j_3, k_3, l_3 ) + A6( i_3, j_3, k_3, l_3 ) );


    return A;
};

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_dCof_dF( const Tensor2_3d &F, const Tensor2_3d &I ) const
{
    Tensor4_3d A2;
    Tensor2_3d cofF          = this->compute_surface_cofactor( F, I );
    A2( i_3, j_3, k_3, l_3 ) = -cofF( i_3, j_3 ) * cofF( k_3, l_3 );

    Tensor4_3d A3, A4, A5, A6, A_alt;
    Tensor2_3d C, B;
    C( i_3, j_3 )         = F( k_3, i_3 ) * F( k_3, j_3 );
    B( i_3, j_3 )         = F( i_3, k_3 ) * F( j_3, k_3 );
    FloatArrayF<9> indent = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Ib( indent );
    Tensor2_3d Is( I );
    auto trC = C( i_3, i_3 );

    A3( i_3, j_3, k_3, l_3 ) = trC * Ib( i_3, k_3 ) * Is( j_3, l_3 ) - Ib( i_3, k_3 ) * C( j_3, l_3 );
    A4( i_3, j_3, k_3, l_3 ) = 2 * F( i_3, j_3 ) * F( k_3, l_3 );
    A5( i_3, j_3, k_3, l_3 ) = -F( i_3, l_3 ) * F( k_3, j_3 );
    A6( i_3, j_3, k_3, l_3 ) = -B( i_3, k_3 ) * Is( j_3, l_3 );

    A_alt( i_3, j_3, k_3, l_3 ) = 1. / this->compute_surface_determinant( F ) * ( A2( i_3, j_3, k_3, l_3 ) + A3( i_3, j_3, k_3, l_3 ) + A4( i_3, j_3, k_3, l_3 ) + A5( i_3, j_3, k_3, l_3 ) + A6( i_3, j_3, k_3, l_3 ) );

    return A_alt;
};


Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_d_normF_dF( const Tensor2_3d &F ) const
{
    Tensor2_3d answer;
    answer( i_3, j_3 ) = 1. / this->compute_surface_NormF( F ) * F( i_3, j_3 );
    return answer;
};

double BaseHyperElasticSurfaceMaterial::compute_surface_NormF( const Tensor2_3d &F ) const
{
    return sqrt( F( i_3, j_3 ) * F( i_3, j_3 ) );
}

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_d2_normF_dF2( const Tensor2_3d &F ) const
{
    Tensor4_3d answer, part1, part2;
    double normF = this->compute_surface_NormF( F );

    part1( 0, 0, 0, 0 ) = 1;
    part1( 1, 0, 1, 0 ) = 1;
    part1( 2, 0, 2, 0 ) = 1;
    part1( 0, 1, 0, 1 ) = 1;
    part1( 1, 1, 1, 1 ) = 1;
    part1( 2, 1, 2, 1 ) = 1;
    part1( 0, 2, 0, 2 ) = 1;
    part1( 1, 2, 1, 2 ) = 1;
    part1( 2, 2, 2, 2 ) = 1;

    part1( i_3, j_3, k_3, l_3 ) = 1 / normF * part1( i_3, j_3, k_3, l_3 );

    part2( i_3, j_3, k_3, l_3 ) = -1 / pow( normF, 3 ) * F( i_3, j_3 ) * F( k_3, l_3 );

    answer( i_3, j_3, k_3, l_3 ) = part1( i_3, j_3, k_3, l_3 ) + part2( i_3, j_3, k_3, l_3 );
    return answer;
};

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_d2_normF_dF2( const Tensor2_3d &F, const Tensor2_3d &I ) const
{
    Tensor4_3d answer, part1, part2;
    double normF = this->compute_surface_NormF( F );

    FloatArrayF<9> indent = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Ib( indent );

    part1( i_3, j_3, k_3, l_3 ) = 1. / normF * Ib( i_3, k_3 ) * I( j_3, l_3 );
    part2( i_3, j_3, k_3, l_3 ) = -1. / pow( normF, 3 ) * F( i_3, j_3 ) * F( k_3, l_3 );

    answer( i_3, j_3, k_3, l_3 ) = part1( i_3, j_3, k_3, l_3 ) + part2( i_3, j_3, k_3, l_3 );
    //// compare
    // Tensor4_3d Acomp = compute_surface_d2_normF_dF2( F );
    // Tensor4_3d dA;
    // dA( i_3, j_3, k_3, l_3 ) = answer( i_3, j_3, k_3, l_3 ) - Acomp( i_3, j_3, k_3, l_3 );

    return answer;
};

Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_d_traceU_dF( const Tensor2_3d &F ) const
{
    FloatMatrix Fm;
    F.toFloatMatrix( Fm );


    Tensor2_3d C;
    // Compute C
    C( j_3, k_3 ) = F( i_3, j_3 ) * F( i_3, k_3 );
    FloatMatrix Cmat;
    C.toFloatMatrix( Cmat );
    Cmat.resizeWithData( 2, 2 );

    // Compute svd of C
    FloatMatrix Uc, Sc, Vc;
    Cmat.computeSVD2x2( Uc, Sc, Vc );

    // Get svd of F
    double lam1 = sqrt( Sc.at( 1, 1 ) );
    double lam2 = sqrt( Sc.at( 2, 2 ) );

    FloatMatrix Uf( Uc ), Vf( Vc ), Sinv, VSinv;
    Vf.resizeWithData( 3, 3 );
    Vf.at( 3, 3 ) = 1;

    Sinv.resize( 3, 3 );
    Sinv.at( 1, 1 ) = 1. / lam1;
    Sinv.at( 2, 2 ) = 1. / lam2;

    VSinv.beProductOf( Vf, Sinv );
    Uf.beProductOf( Fm, VSinv );

    // Compute third column as cross product of the first two
    Uf.at( 1, 3 ) = Uf.at( 2, 1 ) * Uf.at( 3, 2 ) - Uf.at( 3, 1 ) * Uf.at( 2, 2 );
    Uf.at( 2, 3 ) = Uf.at( 3, 1 ) * Uf.at( 1, 2 ) - Uf.at( 1, 1 ) * Uf.at( 3, 2 );
    Uf.at( 3, 3 ) = Uf.at( 1, 1 ) * Uf.at( 2, 2 ) - Uf.at( 2, 1 ) * Uf.at( 1, 2 );

    // Compute rotation R
    FloatMatrix R;
    R.beProductTOf( Uf, Vf );

    Tensor2_3d answer( R );
    return answer;
};

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_d2_traceU_dF2( const Tensor2_3d &F ) const
{
    FloatMatrix Fm;
    F.toFloatMatrix( Fm );


    Tensor2_3d C;
    // Compute C
    C( j_3, k_3 ) = F( i_3, j_3 ) * F( i_3, k_3 );
    FloatMatrix Cmat;
    C.toFloatMatrix( Cmat );
    Cmat.resizeWithData( 2, 2 );

    // Compute svd of C
    FloatMatrix Uc, Sc, Vc;
    Cmat.computeSVD2x2( Uc, Sc, Vc );

    // Get svd of F
    double lam1 = sqrt( Sc.at( 1, 1 ) );
    double lam2 = sqrt( Sc.at( 2, 2 ) );

    FloatMatrix Uf( Uc ), Vf( Vc ), Sinv, VSinv;
    Vf.resizeWithData( 3, 3 );
    Vf.at( 3, 3 ) = 1;

    Sinv.resize( 3, 3 );
    Sinv.at( 1, 1 ) = 1 / lam1;
    Sinv.at( 2, 2 ) = 1 / lam2;

    VSinv.beProductOf( Vf, Sinv );
    Uf.beProductOf( Fm, VSinv );

    // Compute third column as cross product of the first two
    Uf.at( 1, 3 ) = Uf.at( 2, 1 ) * Uf.at( 3, 2 ) - Uf.at( 3, 1 ) * Uf.at( 2, 2 );
    Uf.at( 2, 3 ) = Uf.at( 3, 1 ) * Uf.at( 1, 2 ) - Uf.at( 1, 1 ) * Uf.at( 3, 2 );
    Uf.at( 3, 3 ) = Uf.at( 1, 1 ) * Uf.at( 2, 2 ) - Uf.at( 2, 1 ) * Uf.at( 1, 2 );

    Vf.resizeWithData( 3, 3 );
    Vf.at( 3, 3 ) = 1;


    // Compute T matrices
    FloatMatrix I1, I2, I3, T1, T2, T3, temp1, temp2, temp3;

    I1.resize( 3, 3 );
    I2.resize( 3, 3 );
    I3.resize( 3, 3 );

    I1.at( 2, 3 ) = -1;
    I1.at( 3, 2 ) = 1;
    I2.at( 1, 3 ) = -1;
    I2.at( 3, 1 ) = 1;
    I3.at( 1, 2 ) = -1;
    I3.at( 2, 1 ) = 1;

    temp1.beProductTOf( I1, Vf );
    T1.beProductOf( Uf, temp1 );

    temp2.beProductTOf( I2, Vf );
    T2.beProductOf( Uf, temp2 );

    temp3.beProductTOf( I3, Vf );
    T3.beProductOf( Uf, temp3 );

    Tensor2_3d T1t( T1 ), T2t( T2 ), T3t( T3 );

    Tensor4_3d answer;
    answer( i_3, j_3, k_3, l_3 ) = 1 / lam2 * T1t( i_3, j_3 ) * T1t( k_3, l_3 ) + 1 / lam1 * T2t( i_3, j_3 ) * T2t( k_3, l_3 ) + 1 / ( lam1 + lam2 ) * T3t( i_3, j_3 ) * T3t( k_3, l_3 );

    return answer;
};

FloatArray BaseHyperElasticSurfaceMaterial::compute_deformed_normal( const Tensor2_3d &F ) const
{
    FloatArray normal( 3 );
    FloatMatrix Fm;
    F.toFloatMatrix( Fm );
    // Compute cross product of first two columns
    normal.at( 1 ) = Fm.at( 2, 1 ) * Fm.at( 3, 2 ) - Fm.at( 3, 1 ) * Fm.at( 2, 2 );
    normal.at( 2 ) = Fm.at( 3, 1 ) * Fm.at( 1, 2 ) - Fm.at( 1, 1 ) * Fm.at( 3, 2 );
    normal.at( 3 ) = Fm.at( 1, 1 ) * Fm.at( 2, 2 ) - Fm.at( 2, 1 ) * Fm.at( 1, 2 );

    normal = 1 / this->compute_surface_determinant( F ) * normal;

    return normal;
}

Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_a_dot_d_normal_dF( const Tensor2_3d &F, const FloatArray &a ) const
{
    double J = this->compute_surface_determinant( F );
    Tensor2_3d part1, part2, result;
    // part 1
    part1( 0, 0 )     = F( 1, 1 ) * a( 2 ) - F( 2, 1 ) * a( 1 );
    part1( 1, 0 )     = F( 2, 1 ) * a( 0 ) - F( 0, 1 ) * a( 2 );
    part1( 2, 0 )     = F( 0, 1 ) * a( 1 ) - F( 1, 1 ) * a( 0 );
    part1( 0, 1 )     = F( 2, 0 ) * a( 1 ) - F( 1, 0 ) * a( 2 );
    part1( 1, 1 )     = F( 0, 0 ) * a( 2 ) - F( 2, 0 ) * a( 0 );
    part1( 2, 1 )     = F( 1, 0 ) * a( 0 ) - F( 0, 0 ) * a( 1 );
    part1( i_3, j_3 ) = part1( i_3, j_3 ) / J;

    // part 02
    FloatArray normal = compute_deformed_normal( F );
    double mult       = normal.dotProduct( a );
    part2( i_3, j_3 ) = -mult / J * compute_surface_cofactor( F )( i_3, j_3 );

    // result
    result( i_3, j_3 ) = part1( i_3, j_3 ) + part2( i_3, j_3 );
    return result;
}


Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_a_dot_d2_normal_d2F( const Tensor2_3d &F, const FloatArray &a ) const
{
    double J          = this->compute_surface_determinant( F );
    FloatArray normal = compute_deformed_normal( F );
    double mult       = normal.dotProduct( a );

    Tensor4_3d result, result1, result2, result3, result4;
    Tensor2_3d cofF, a_dot_dndF, res2, res4;
    cofF = this->compute_surface_cofactor( F );

    // First part
    result1( 1, 1, 0, 0 ) = a( 2 );
    result1( 2, 1, 0, 0 ) = -a( 1 );
    result1( 0, 1, 1, 0 ) = -a( 2 );
    result1( 2, 1, 1, 0 ) = a( 0 );
    result1( 0, 1, 2, 0 ) = a( 1 );
    result1( 1, 1, 2, 0 ) = -a( 0 );

    result1( 1, 0, 0, 1 )         = -a( 2 );
    result1( 2, 0, 0, 1 )         = a( 1 );
    result1( 0, 0, 1, 1 )         = a( 2 );
    result1( 2, 0, 1, 1 )         = -a( 0 );
    result1( 0, 0, 2, 1 )         = -a( 1 );
    result1( 1, 0, 2, 1 )         = a( 0 );
    result1( i_3, j_3, k_3, l_3 ) = 1 / J * result1( i_3, j_3, k_3, l_3 );
    // second part
    res2( 0, 0 )                  = F( 1, 1 ) * a( 2 ) - F( 2, 1 ) * a( 1 );
    res2( 1, 0 )                  = F( 2, 1 ) * a( 0 ) - F( 0, 1 ) * a( 2 );
    res2( 2, 0 )                  = F( 0, 1 ) * a( 1 ) - F( 1, 1 ) * a( 0 );
    res2( 0, 1 )                  = F( 2, 0 ) * a( 1 ) - F( 1, 0 ) * a( 2 );
    res2( 1, 1 )                  = F( 0, 0 ) * a( 2 ) - F( 2, 0 ) * a( 0 );
    res2( 2, 1 )                  = F( 1, 0 ) * a( 0 ) - F( 0, 0 ) * a( 1 );
    result2( i_3, j_3, k_3, l_3 ) = -1 / J / J * res2( i_3, j_3 ) * cofF( k_3, l_3 ) - 1 / J / J * res2( k_3, l_3 ) * cofF( i_3, j_3 );
    // third part
    result3                       = this->compute_surface_dCof_dF( F );
    result3( i_3, j_3, k_3, l_3 ) = -1 / J * mult * result3( i_3, j_3, k_3, l_3 );
    // fourth part
    result4( i_3, j_3, k_3, l_3 ) = 2 * mult / J / J * cofF( i_3, j_3 ) * cofF( k_3, l_3 );
    // sum
    result( i_3, j_3, k_3, l_3 ) = result1( i_3, j_3, k_3, l_3 ) + result2( i_3, j_3, k_3, l_3 ) + result3( i_3, j_3, k_3, l_3 ) + result4( i_3, j_3, k_3, l_3 );
    return result;
}

Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_d_normFa_dF( const Tensor2_3d &F, const Tensor1_3d &a ) const
{
    Tensor1_3d Fa;
    Tensor2_3d answer;
    Fa( i_3 ) = F( i_3, j_3 ) * a( j_3 );
    double normFa = sqrt( Fa( i_3 ) * Fa( i_3 ) );
    answer( i_3, j_3 ) = 1 / normFa * F( i_3, k_3 ) * a( k_3 ) * a( j_3 );
    return answer;
}

 Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_d2_normFa_dF2( const Tensor2_3d &F, const Tensor1_3d &a ) const
{
     FloatArrayF<9> indent3D = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
     Tensor2_3d Identity( indent3D ), aoa;
     Tensor1_3d Fa;
     Tensor4_3d answer;
     Fa( i_3 )     = F( i_3, j_3 ) * a( j_3 );
     double normFa = sqrt( Fa( i_3 ) * Fa( i_3 ) );

     aoa( i_3, j_3 ) = a( i_3 ) * a( j_3 );
     answer( i_3, j_3, k_3, l_3 ) = ( 1 / normFa * Identity( i_3, k_3 ) * aoa( j_3, l_3 ) - 1 / normFa / normFa / normFa * F( i_3, m_3 ) * aoa( m_3, j_3 ) * F( k_3, n_3 ) * aoa( n_3, l_3 ) );
     return answer;
 }

Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_d_SqrtNormCa_dF( const Tensor2_3d &F, const Tensor1_3d &a ) const
{
    Tensor2_3d answer, C;
    Tensor1_3d Ca;
    C( i_3, j_3 ) = F( k_3, i_3 ) * F( k_3, j_3 );
    Ca( i_3 )     = C( i_3, j_3 ) * a( j_3 );
    double CaNorm = sqrt( Ca( i_3 ) * Ca( i_3 ) );

    answer( i_3, l_3 ) = 0.5 / CaNorm / sqrt( CaNorm ) * F( i_3, j_3 ) * ( Ca( j_3 ) * a( l_3 ) + a( j_3 ) * Ca( l_3 ) );
    return answer;
}

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_d2_SqrtNormCa_dF2( const Tensor2_3d &F, const Tensor1_3d &a ) const
{
    Tensor4_3d answer, Aaniso2_a, Aaniso2_b, Aaniso2_c;
    Tensor1_3d Ca, Fa;
    FloatArrayF<9> indent3D = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Identity( indent3D ), aoa, C, aFa, Faa;

    aoa( i_3, j_3 ) = a( i_3 ) * a( j_3 );
    C( i_3, j_3 )   = F( k_3, i_3 ) * F( k_3, j_3 );
    Ca( i_3 )       = C( i_3, j_3 ) * a( j_3 );
    double CaNorm   = sqrt( Ca( i_3 ) * Ca( i_3 ) );
    
    Fa( i_3 )       = F( i_3, j_3 ) * a( j_3 );
    aFa( i_3, k_3 ) = a( i_3 ) * Fa( k_3 );
    Faa( i_3, k_3 ) = Fa( i_3 ) * a( k_3 );
    

    Aaniso2_a( i_3, l_3, o_3, p_3 ) = -3. / 4. / pow( CaNorm, 3 ) / sqrt( CaNorm ) * F( i_3, j_3 ) * ( Ca( j_3 ) * a( l_3 ) + a( j_3 ) * Ca( l_3 ) ) * F( o_3, k_3 ) * ( Ca( k_3 ) * a( p_3 ) + a( k_3 ) * Ca( p_3 ) );
    Aaniso2_b( i_3, l_3, o_3, p_3 ) = 0.5 / CaNorm / sqrt( CaNorm ) * ( F( i_3, p_3 ) * aFa( l_3, o_3 ) + F( i_3, j_3 ) * F( o_3, j_3 ) * aoa( l_3, p_3 ) + Fa( i_3 ) * Fa( o_3 ) * Identity( l_3, p_3 ) + Faa( i_3, p_3 ) * F( o_3, l_3 ) );
    Aaniso2_c( i_3, l_3, o_3, p_3 ) = 0.5 / CaNorm / sqrt( CaNorm ) * Identity( i_3, o_3 ) * ( Ca( l_3 ) * a( p_3 ) + a( l_3 ) * Ca( p_3 ) );
    answer( i_3, l_3, o_3, p_3 ) = Aaniso2_a( i_3, l_3, o_3, p_3 ) + Aaniso2_b( i_3, l_3, o_3, p_3 ) + Aaniso2_c( i_3, l_3, o_3, p_3 );

    return answer;
}

} // end namespace oofem
