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
    double J        = std::sqrt( 0.5 * ( trC * trC - trCtC ) );

    return J;
};

Tensor1_3d BaseHyperElasticSurfaceMaterial::compute_surface_deformed_normal( const Tensor2_3d &F, const Tensor1_3d &N ) const
{
    // n = cof(F) N / |cof(F) N|: cof(F) maps the reference area vector N dA to the current one n da
    Tensor2_3d cofF = F.compute_cofactor();
    Tensor1_3d n;
    n( i_3 ) = cofF( i_3, j_3 ) * N( j_3 );
    n.normalize();
    return n;
};

Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_pseudoinverse( const Tensor2_3d &F, const Tensor2_3d &I ) const
{
    // dJ/dF = J F^+T, so F^+ = (dJ/dF)^T / J
    Tensor2_3d dJdF = this->compute_surface_cofactor( F, I );
    Tensor2_3d Finv;
    Finv( i_3, k_3 ) = dJdF( k_3, i_3 ) / this->compute_surface_determinant( F);
    return Finv;
};


Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_cofactor( const Tensor2_3d &F, const Tensor2_3d &I ) const
{
    // J^2 = ((tr C)^2 - tr C^2) / 2 gives 2 J dJ/dF = 2 F (tr C I - C); F I = F because F N = 0
    Tensor2_3d cofF, C, JcofF;
    C( i_3, j_3 )     = F( k_3, i_3 ) * F( k_3, j_3 );
    auto trC          = C( i_3, i_3 );
    JcofF( i_3, k_3 ) = F( i_3, j_3 ) * ( trC * I( j_3, k_3 ) - C( j_3, k_3 ) ); // J dJ/dF
    cofF( i_3, j_3 )  = 1. / this->compute_surface_determinant( F ) * JcofF( i_3, j_3 );

    return cofF;
};

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_dCof_dF( const Tensor2_3d &F, const Tensor2_3d &I ) const
{
    // d2J/dF2 = ( d(F (tr C I - C))/dF - dJ/dF x dJ/dF ) / J, where (B = F F^T)
    // d(F (tr C I - C))_ij / dF_kl = tr C delta_ik I_jl - delta_ik C_jl + 2 F_ij F_kl - F_il F_kj - B_ik I_jl.
    // The exact derivative has delta_jl instead of I_jl in the first and the last term; the difference
    // (tr C delta_ik - B_ik) N_j N_l does not act on variations with dF N = 0.
    Tensor4_3d A2;
    Tensor2_3d cofF          = this->compute_surface_cofactor( F, I );
    A2( i_3, j_3, k_3, l_3 ) = -cofF( i_3, j_3 ) * cofF( k_3, l_3 ); // -dJ/dF x dJ/dF

    Tensor4_3d A3, A4, A5, A6, d2J;
    Tensor2_3d C, B;
    C( i_3, j_3 )         = F( k_3, i_3 ) * F( k_3, j_3 );
    B( i_3, j_3 )         = F( i_3, k_3 ) * F( j_3, k_3 );
    FloatArrayF<9> vI     = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Ib( vI ); // 3D identity
    Tensor2_3d Is( I );
    auto trC = C( i_3, i_3 );

    A3( i_3, j_3, k_3, l_3 ) = trC * Ib( i_3, k_3 ) * Is( j_3, l_3 ) - Ib( i_3, k_3 ) * C( j_3, l_3 ); // from dF (tr C I - C)
    A4( i_3, j_3, k_3, l_3 ) = 2 * F( i_3, j_3 ) * F( k_3, l_3 ); // from d tr C = 2 F : dF
    A5( i_3, j_3, k_3, l_3 ) = -F( i_3, l_3 ) * F( k_3, j_3 ); // A5 + A6 from -F dC, dC = dF^T F + F^T dF
    A6( i_3, j_3, k_3, l_3 ) = -B( i_3, k_3 ) * Is( j_3, l_3 );

    d2J( i_3, j_3, k_3, l_3 ) = 1. / this->compute_surface_determinant( F ) * ( A2( i_3, j_3, k_3, l_3 ) + A3( i_3, j_3, k_3, l_3 ) + A4( i_3, j_3, k_3, l_3 ) + A5( i_3, j_3, k_3, l_3 ) + A6( i_3, j_3, k_3, l_3 ) );

    return d2J;
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

Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_d2_normF_dF2( const Tensor2_3d &F, const Tensor2_3d &I ) const
{
    Tensor4_3d answer, part1, part2;
    double normF = this->compute_surface_NormF( F );

    FloatArrayF<9> vI = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Ib( vI ); // 3D identity

    // derivative of F / |F|: delta_ik I_jl / |F| (second index projected by I) - F_ij F_kl / |F|^3
    part1( i_3, j_3, k_3, l_3 ) = 1. / normF * Ib( i_3, k_3 ) * I( j_3, l_3 );
    part2( i_3, j_3, k_3, l_3 ) = -1. / pow( normF, 3 ) * F( i_3, j_3 ) * F( k_3, l_3 );

    answer( i_3, j_3, k_3, l_3 ) = part1( i_3, j_3, k_3, l_3 ) + part2( i_3, j_3, k_3, l_3 );

    return answer;
};

Tensor2_3d BaseHyperElasticSurfaceMaterial::compute_surface_d_normFa_dF( const Tensor2_3d &F, const Tensor1_3d &a ) const
{
    // |F a| = sqrt(F a . F a), so d|F a|/dF = (F a) x a / |F a|
    Tensor1_3d Fa;
    Tensor2_3d answer;
    Fa( i_3 ) = F( i_3, j_3 ) * a( j_3 );
    double normFa = sqrt( Fa( i_3 ) * Fa( i_3 ) );
    answer( i_3, j_3 ) = 1 / normFa * F( i_3, k_3 ) * a( k_3 ) * a( j_3 );
    return answer;
}

 Tensor4_3d BaseHyperElasticSurfaceMaterial::compute_surface_d2_normFa_dF2( const Tensor2_3d &F, const Tensor1_3d &a ) const
{
     // derivative of (F a) x a / |F a|: delta_ik a_j a_l / |F a| - (F a)_i a_j (F a)_k a_l / |F a|^3
     FloatArrayF<9> vI = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
     Tensor2_3d Identity( vI ), aoa;
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
    // d|C a|/dF = F (C a x a + a x C a) / |C a| (from d(C a)_j / dF_il = delta_jl (F a)_i + F_ij a_l),
    // so d sqrt(|C a|)/dF = F (C a x a + a x C a) / (2 |C a|^(3/2))
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
    // Derivative with respect to F_op of 0.5 |C a|^(-3/2) F_ij ((C a)_j a_l + a_j (C a)_l), in three terms:
    // _a from the factor |C a|^(-3/2), _b from C a in the bracket (d(C a)_j / dF_op = delta_jp (F a)_o + F_oj a_p),
    // _c from F_ij itself.
    Tensor4_3d answer, Aaniso2_a, Aaniso2_b, Aaniso2_c;
    Tensor1_3d Ca, Fa;
    FloatArrayF<9> vI = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Identity( vI ), aoa, C, aFa, Faa;

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
