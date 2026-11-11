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
    // second invariant of C; C has the eigenvalues lambda1^2, lambda2^2 and 0, so this is lambda1 lambda2 for any N
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
    // Explicit formula for N = e3 (F with a zero third column)
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

    // cofF2 (the general formula of the version with I, for I = diag(1, 1, 0)) is computed but not used
    Tensor2_3d C, cofF2;
    C( i_3, j_3 )         = F( k_3, i_3 ) * F( k_3, j_3 );
    FloatArrayF<9> indent = { 1., 1., 0., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d I( indent );
    auto trC          = C( i_3, i_3 );
    cofF2( i_3, k_3 ) = F( i_3, j_3 ) * ( trC * I( j_3, k_3 ) - C( j_3, k_3 ) );

    cofF( i_3, j_3 ) = 1. / this->compute_surface_determinant( F ) * cofF( i_3, j_3 );

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
    // Explicit formula for N = e3 (F with a zero third column); the version with I below is the general form
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

    return A;
};

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_dCof_dF( const Tensor2_3d &F, const Tensor2_3d &I ) const
{
    // d2J/dF2 = ( d(F (tr C I - C))/dF - dJ/dF x dJ/dF ) / J
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
