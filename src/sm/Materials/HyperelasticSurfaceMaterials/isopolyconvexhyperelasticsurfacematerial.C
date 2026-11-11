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

#include "isopolyconvexhyperelasticsurfacematerial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"
#include "domain.h"
#include "function.h"
#include "gausspoint.h"
#include "element.h"

#include <algorithm>
#include <cmath>

namespace oofem {
REGISTER_Material( IsotropicPolyconvexHyperelasticSurfaceMaterial );

IsotropicPolyconvexHyperelasticSurfaceMaterial::IsotropicPolyconvexHyperelasticSurfaceMaterial( int n, Domain *d ) :
    HyperElasticSurfaceMaterial( n, d ),
    ArcLengthMaterialInterface()
{}


void IsotropicPolyconvexHyperelasticSurfaceMaterial::readParameter( InputRecord &ir, SurfaceMaterialParameter &p, const char *name, const char *ltfName )
{
    IR_GIVE_OPTIONAL_FIELD( ir, p.value, name );
    IR_GIVE_OPTIONAL_FIELD( ir, p.ltf, ltfName );
}


void IsotropicPolyconvexHyperelasticSurfaceMaterial::initializeFrom( InputRecord &ir )
{
    HyperElasticSurfaceMaterial::initializeFrom( ir );

    // Keywords of the earlier versions whose meaning changed: stop instead of silently computing something else.
    // There, gamma_ltf scaled all parameters; now every parameter has its own time function.
    for ( const char *removed : { "fixgamma", "fixomega", "fixalphas", "aniso_ltf" } ) {
        if ( ir.hasField( removed ) ) {
            throw ValueInputException( ir, removed, "removed; every parameter has its own time function <name>_ltf, without it the value is constant" );
        }
    }

    readParameter( ir, this->gamma, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gamma1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gammaLTF );
    readParameter( ir, this->alpha, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alphaLTF );
    readParameter( ir, this->alphaS, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha_s1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha_sLTF );
    readParameter( ir, this->lambda, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_lambda, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_lambdaLTF );
    readParameter( ir, this->mu, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_mu, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_muLTF );
    IR_GIVE_OPTIONAL_FIELD( ir, this->gaminit, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gaminit );

    // older inputs with a single alpha and gamma
    double alphaTemp = 0., gammaTemp = 0.;
    IR_GIVE_OPTIONAL_FIELD( ir, alphaTemp, "alpha" );
    IR_GIVE_OPTIONAL_FIELD( ir, gammaTemp, "gamma" );
    if ( alphaTemp != 0. && this->alpha.value == 0. ) {
        this->alpha.value = alphaTemp;
    }
    if ( gammaTemp != 0. && this->gamma.value == 0. ) {
        this->gamma.value = gammaTemp;
    }

    if ( !this->supportsSecondGradient() ) {
        double omega1 = 0.;
        IR_GIVE_OPTIONAL_FIELD( ir, omega1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_omega1 );
        if ( omega1 != 0. ) {
            throw ValueInputException( ir, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_omega1, "this material has no bending stiffness, use surfacetensionmat2" );
        }
    }

    int par = calm_par_no;
    IR_GIVE_OPTIONAL_FIELD( ir, par, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_ArcLengthParameter );
    if ( par < calm_par_no || par > calm_par_omega || ( par == calm_par_omega && !this->supportsSecondGradient() ) ) {
        throw ValueInputException( ir, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_ArcLengthParameter,
            "must be 0 (none), 1 (gamma) or 2 (omega, only surfacetensionmat2)" );
    }
    this->arcLengthParameter = (CALM_parameter)par;
    this->hasALP             = par != calm_par_no;

    IR_GIVE_OPTIONAL_FIELD( ir, this->checkTangentTolerance, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_checkTangent );
}


Interface *IsotropicPolyconvexHyperelasticSurfaceMaterial::giveInterface( InterfaceType t )
{
    if ( t == ArcLengthMaterialInterfaceType ) {
        return static_cast<ArcLengthMaterialInterface *>( this );
    }
    return HyperElasticSurfaceMaterial::giveInterface( t );
}


MaterialStatus *IsotropicPolyconvexHyperelasticSurfaceMaterial::CreateStatus( GaussPoint *gp ) const
{
    return new StructuralSurfaceMaterialStatus( gp );
}


double IsotropicPolyconvexHyperelasticSurfaceMaterial::giveParameter( const SurfaceMaterialParameter &p, CALM_parameter role, TimeStep *tStep ) const
{
    if ( this->hasArcLengthParameter() ) {
        // The arc-length method assembles first the part of the internal forces which multiplies its load
        // parameter (returnAll false), then all of them; the time functions are not used.
        if ( role != calm_par_no && role == this->arcLengthParameter ) {
            return this->returnAll ? this->lambdaCALM * p.value : p.value;
        }
        return this->returnAll ? p.value : 0.;
    }
    if ( p.ltf ) {
        return p.value * this->domain->giveFunction( p.ltf )->evaluateAtTime( tStep->giveIntrinsicTime() );
    }
    return p.value;
}


double IsotropicPolyconvexHyperelasticSurfaceMaterial::giveGamma( TimeStep *tStep ) const
{
    double g = this->giveParameter( this->gamma, calm_par_gamma, tStep );
    if ( !this->hasArcLengthParameter() && this->gamma.ltf ) {
        g += this->gaminit;
    }
    return g;
}


Tensor2_3d IsotropicPolyconvexHyperelasticSurfaceMaterial::giveSurfaceIdentity( const FloatArray &normal )
{
    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Identity3d( indent3d ), Is;
    FloatArrayF<3> vN( normal );
    Tensor1_3d N( vN );
    Is( i_3, j_3 ) = Identity3d( i_3, j_3 ) - N( i_3 ) * N( j_3 );
    return Is;
}


Tensor2_3d IsotropicPolyconvexHyperelasticSurfaceMaterial::computeMembraneStress( const Tensor2_3d &F, const Tensor2_3d &Is, TimeStep *tStep ) const
{
    double gamma_t  = this->giveGamma( tStep );
    double alpha_t  = this->giveParameter( this->alpha, calm_par_no, tStep );
    double alphaS_t = this->giveParameter( this->alphaS, calm_par_no, tStep );
    double lambda_t = this->giveParameter( this->lambda, calm_par_no, tStep );
    double mu_t     = this->giveParameter( this->mu, calm_par_no, tStep );

    Tensor2_3d dJdF, Pgamma, Palpha, Pnh, Palpha_s, P;
    dJdF( i_3, j_3 )   = this->compute_surface_cofactor( F, Is )( i_3, j_3 );
    Pgamma( i_3, j_3 ) = gamma_t * dJdF( i_3, j_3 );
    Palpha( i_3, j_3 ) = alpha_t * this->compute_surface_d_normF_dF( F )( i_3, j_3 );

    double J    = this->compute_surface_determinant( F );
    double logJ = log( J );
    // neo-Hookean part mu (|F|^2 / 2 - 1 - ln J) + lambda / 2 (ln J)^2
    Pnh( i_3, j_3 ) = dJdF( i_3, j_3 ) / J * ( lambda_t * logJ - mu_t ) + mu_t * F( i_3, j_3 );
    // area part alphas / 2 (J - 1)^2
    Palpha_s( i_3, j_3 ) = alphaS_t * ( J - 1. ) * dJdF( i_3, j_3 );

    P( i_3, j_3 ) = Pgamma( i_3, j_3 ) + Palpha( i_3, j_3 ) + Pnh( i_3, j_3 ) + Palpha_s( i_3, j_3 );
    return P;
}


Tensor4_3d IsotropicPolyconvexHyperelasticSurfaceMaterial::computeMembraneTangent( const Tensor2_3d &F, const Tensor2_3d &Is, TimeStep *tStep ) const
{
    double gamma_t  = this->giveGamma( tStep );
    double alpha_t  = this->giveParameter( this->alpha, calm_par_no, tStep );
    double alphaS_t = this->giveParameter( this->alphaS, calm_par_no, tStep );
    double lambda_t = this->giveParameter( this->lambda, calm_par_no, tStep );
    double mu_t     = this->giveParameter( this->mu, calm_par_no, tStep );

    FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
    Tensor2_3d Identity3d( indent3d ), dJdF, IdentityUD( Is );
    Tensor4_3d d2JdFdF, Anh, A_alpha_s, A;

    d2JdFdF( i_3, j_3, k_3, l_3 ) = this->compute_surface_dCof_dF( F, Is )( i_3, j_3, k_3, l_3 );
    double J                      = this->compute_surface_determinant( F );
    double logJ                   = log( J );
    dJdF( i_3, j_3 )              = this->compute_surface_cofactor( F, Is )( i_3, j_3 );

    Anh( i_3, j_3, k_3, l_3 ) = lambda_t * ( logJ / J * d2JdFdF( i_3, j_3, k_3, l_3 ) + ( 1 - logJ ) / J / J * dJdF( i_3, j_3 ) * dJdF( k_3, l_3 ) ) +
        mu_t * ( Identity3d( i_3, k_3 ) * IdentityUD( j_3, l_3 ) - 1 / J * d2JdFdF( i_3, j_3, k_3, l_3 ) + 1 / J / J * dJdF( i_3, j_3 ) * dJdF( k_3, l_3 ) );

    A_alpha_s( i_3, j_3, k_3, l_3 ) = alphaS_t * ( J - 1. ) * d2JdFdF( i_3, j_3, k_3, l_3 ) + alphaS_t * dJdF( i_3, j_3 ) * dJdF( k_3, l_3 );

    A( i_3, j_3, k_3, l_3 ) = gamma_t * d2JdFdF( i_3, j_3, k_3, l_3 ) + alpha_t * this->compute_surface_d2_normF_dF2( F, Is )( i_3, j_3, k_3, l_3 ) +
        Anh( i_3, j_3, k_3, l_3 ) + A_alpha_s( i_3, j_3, k_3, l_3 );
    return A;
}


void IsotropicPolyconvexHyperelasticSurfaceMaterial::storeStress( GaussPoint *gp, const FloatArrayF<9> &vF, const Tensor2_3d &F, const Tensor2_3d &P ) const
{
    auto status = static_cast<StructuralMaterialStatus *>( this->giveStatus( gp ) );
    Tensor2_3d S, Pc( P ); // surface Cauchy stress; to_voigt_form is not const
    S( i_3, j_3 ) = ( 1 / this->compute_surface_determinant( F ) ) * P( i_3, k_3 ) * F( j_3, k_3 );
    auto vP       = Pc.to_voigt_form();
    auto vS       = S.to_voigt_form();
    status->letTempFVectorBe( vF );
    status->letTempPVectorBe( vP );
    status->letTempCVectorBe( vS );
}


FloatArrayF<9>
IsotropicPolyconvexHyperelasticSurfaceMaterial::giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    Tensor2_3d F( vF );
    Tensor2_3d P = this->computeMembraneStress( F, giveSurfaceIdentity( normal ), tStep );
    this->storeStress( gp, vF, F, P );
    return P.to_voigt_form();
}


FloatMatrixF<9, 9>
IsotropicPolyconvexHyperelasticSurfaceMaterial::give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    Tensor2_3d F( vF );
    auto answer = this->computeMembraneTangent( F, giveSurfaceIdentity( normal ), tStep ).to_voigt_form();
    if ( this->checkTangentTolerance != 0. ) {
        this->checkTangentsByFiniteDifferences( answer, normal, gp, tStep );
    }
    return answer;
}


void IsotropicPolyconvexHyperelasticSurfaceMaterial::checkTangentsByFiniteDifferences( const FloatMatrixF<9, 9> &dPdF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    // the perturbed evaluations overwrite the temporary status, which is restored at the end
    FloatArray sF = status->giveTempFVector(), sP = status->giveTempPVector(), sC = status->giveTempCVector();
    FloatArray sG = status->giveTempGVector(), sdI = status->giveTempdIVector(), sM = status->giveTempMVector(), sK = status->giveTempKVector();

    const bool bending = this->supportsSecondGradient();
    if ( bending && ( sG.giveSize() != 27 || sdI.giveSize() != 27 ) ) {
        OOFEM_ERROR( "the second order stress has to be evaluated before the tangent check" );
    }
    FloatArrayF<9> vF( sF );
    FloatArrayF<27> vG, vdI;
    FloatMatrixF<27, 27> dAdG;
    FloatMatrixF<27, 9> dAdF;
    if ( bending ) {
        vG   = FloatArrayF<27>( sG );
        vdI  = FloatArrayF<27>( sdI );
        dAdG = this->give3dSurfaceMaterialStiffnessMatrix_dAddF( TangentStiffness, normal, gp, tStep );
        dAdF = this->give3dSurfaceMaterialStiffnessMatrix_dAdF( TangentStiffness, normal, gp, tStep );
    }

    auto evaluate = [&]( const FloatArrayF<9> &F, const FloatArrayF<27> &G, FloatArrayF<9> &P, FloatArrayF<27> &A ) {
        if ( bending ) {
            A = this->giveSecondOrderSurfaceStressVector_3d( F, G, vdI, normal, gp, tStep );
        }
        P = this->giveFirstPKSurfaceStressVector_3d( F, normal, gp, tStep );
    };

    // tangent basis of the plane normal to N; the perturbations keep F N = 0
    FloatArrayF<3> n( normal );
    FloatArrayF<3> a = std::fabs( n[0] ) < 0.9 ? FloatArrayF<3>{ 1., 0., 0. } : FloatArrayF<3>{ 0., 1., 0. };
    FloatArrayF<3> t1 = cross( n, a );
    t1 = t1 / norm( t1 );
    FloatArrayF<3> t2 = cross( n, t1 );
    FloatArrayF<3> t[2] = { t1, t2 };

    double hF = 1.e-6 * std::max( 1., norm( vF ) );
    double hG = 1.e-6 * std::max( 1., norm( vG ) );
    // maximal difference and magnitude for the blocks dP/dF, dA/dF, dA/dgradF, dP/dgradF
    double diff[4] = { 0., 0., 0., 0. }, size[4] = { 0., 0., 0., 0. };
    auto compare = [&]( int block, const FloatArrayF<9> &fd, const FloatArrayF<9> &an ) {
        for ( int i = 0; i < 9; i++ ) {
            diff[block] = std::max( diff[block], std::fabs( fd[i] - an[i] ) );
            size[block] = std::max( { size[block], std::fabs( fd[i] ), std::fabs( an[i] ) } );
        }
    };
    auto compare27 = [&]( int block, const FloatArrayF<27> &fd, const FloatArrayF<27> &an ) {
        for ( int i = 0; i < 27; i++ ) {
            diff[block] = std::max( diff[block], std::fabs( fd[i] - an[i] ) );
            size[block] = std::max( { size[block], std::fabs( fd[i] ), std::fabs( an[i] ) } );
        }
    };

    FloatArrayF<9> Pp, Pm;
    FloatArrayF<27> Ap, Am;
    for ( int i = 0; i < 3; i++ ) {
        for ( int al = 0; al < 2; al++ ) {
            Tensor2_3d dF;
            for ( int j = 0; j < 3; j++ ) {
                dF( i, j ) = t[al][j];
            }
            FloatArrayF<9> vdF = dF.to_voigt_form();
            evaluate( vF + hF * vdF, vG, Pp, Ap );
            evaluate( vF - hF * vdF, vG, Pm, Am );
            compare( 0, ( Pp - Pm ) / ( 2. * hF ), dot( dPdF, vdF ) );
            if ( bending ) {
                compare27( 1, ( Ap - Am ) / ( 2. * hF ), dot( dAdF, vdF ) );
            }
        }
    }
    if ( bending ) {
        // grad F is symmetric in the last two indices, and so are the perturbations
        for ( int i = 0; i < 3; i++ ) {
            for ( int al = 0; al < 2; al++ ) {
                for ( int be = al; be < 2; be++ ) {
                    Tensor3_3d dG;
                    for ( int j = 0; j < 3; j++ ) {
                        for ( int k = 0; k < 3; k++ ) {
                            dG( i, j, k ) = 0.5 * ( t[al][j] * t[be][k] + t[be][j] * t[al][k] );
                        }
                    }
                    FloatArrayF<27> vdG = dG.to_voigt_form_27();
                    evaluate( vF, vG + hG * vdG, Pp, Ap );
                    evaluate( vF, vG - hG * vdG, Pm, Am );
                    compare27( 2, ( Ap - Am ) / ( 2. * hG ), dot( dAdG, vdG ) );
                    compare( 3, ( Pp - Pm ) / ( 2. * hG ), Tdot( dAdF, vdG ) ); // dP/dgradF = (dA/dF)^T
                }
            }
        }
    }

    // restore the temporary status
    status->letTempFVectorBe( sF );
    status->letTempPVectorBe( sP );
    status->letTempCVectorBe( sC );
    status->letTempGVectorBe( sG );
    status->letTempdIVectorBe( sdI );
    status->letTempMVectorBe( sM );
    status->letTempKVectorBe( sK );

    double ref = std::max( { size[0], size[1], size[2], size[3], 1.e-300 } );
    double err[4];
    for ( int b = 0; b < 4; b++ ) {
        err[b] = diff[b] / std::max( size[b], 1.e-8 * ref );
    }
    double worst = *std::max_element( err, err + 4 );
    if ( this->checkTangentTolerance < 0. ) {
        OOFEM_LOG_INFO( "tangent check, element %d gp %d: dP/dF %.2e dA/dF %.2e dA/dgradF %.2e dP/dgradF %.2e\n",
            gp->giveElement()->giveNumber(), gp->giveNumber(), err[0], err[1], err[2], err[3] );
    } else if ( worst > this->checkTangentTolerance ) {
        OOFEM_ERROR( "tangent check failed, element %d gp %d: relative differences dP/dF %.2e dA/dF %.2e dA/dgradF %.2e dP/dgradF %.2e (tolerance %.2e)",
            gp->giveElement()->giveNumber(), gp->giveNumber(), err[0], err[1], err[2], err[3], this->checkTangentTolerance );
    }
}


/////////////////////////////////////
// Membrane with bending stiffness

REGISTER_Material( IsotropicPolyconvexHyperelasticSurfaceMaterial2 );

IsotropicPolyconvexHyperelasticSurfaceMaterial2::IsotropicPolyconvexHyperelasticSurfaceMaterial2( int n, Domain *d ) :
    IsotropicPolyconvexHyperelasticSurfaceMaterial( n, d )
{}


void IsotropicPolyconvexHyperelasticSurfaceMaterial2::initializeFrom( InputRecord &ir )
{
    IsotropicPolyconvexHyperelasticSurfaceMaterial::initializeFrom( ir );
    readParameter( ir, this->omega, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_omega1, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_omegaLTF );
    IR_GIVE_OPTIONAL_FIELD( ir, this->use_H0, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_use_H0 );
    IR_GIVE_OPTIONAL_FIELD( ir, this->H0_val, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_H0 );
    if ( this->use_H0 < 0 || this->use_H0 > 2 ) {
        throw ValueInputException( ir, _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_use_H0, "must be 0 (H0 = 0), 1 (reference curvature) or 2 (value h0)" );
    }
}


double IsotropicPolyconvexHyperelasticSurfaceMaterial2::giveOmega( TimeStep *tStep ) const
{
    return this->giveParameter( this->omega, calm_par_omega, tStep );
}


double IsotropicPolyconvexHyperelasticSurfaceMaterial2::giveH0( const Tensor1_3d &N, const Tensor3_3d &gradIt ) const
{
    if ( this->use_H0 == 1 ) {
        Tensor2_3d K0, b0;
        K0( j_3, k_3 ) = -N( i_3 ) * gradIt( i_3, j_3, k_3 ); // gradient of the reference normal
        b0( j_3, l_3 ) = -K0( j_3, l_3 );
        return 0.5 * b0( j_3, j_3 );
    } else if ( this->use_H0 == 2 ) {
        return this->H0_val;
    }
    return 0.;
}


void IsotropicPolyconvexHyperelasticSurfaceMaterial2::giveSecondGradients( StructuralSurfaceMaterialStatus *status, FloatArrayF<27> &vG, FloatArrayF<27> &vdI )
{
    if ( status->giveTempGVector().giveSize() != 27 || status->giveTempdIVector().giveSize() != 27 ) {
        OOFEM_ERROR( "grad F is missing in the status, the second order stress has to be evaluated first" );
    }
    vG  = FloatArrayF<27>( status->giveTempGVector() );
    vdI = FloatArrayF<27>( status->giveTempdIVector() );
}


IsotropicPolyconvexHyperelasticSurfaceMaterial2::BendingState
IsotropicPolyconvexHyperelasticSurfaceMaterial2::computeBendingState( const Tensor2_3d &F, const Tensor2_3d &Is, const Tensor1_3d &N, const Tensor3_3d &G, const Tensor3_3d &gradIt, TimeStep *tStep ) const
{
    BendingState s;
    s.omega = this->giveOmega( tStep );
    s.J     = this->compute_surface_determinant( F );

    Tensor2_3d Finv = this->compute_surface_pseudoinverse( F, Is );
    Tensor1_3d n    = compute_surface_deformed_normal( F, N );
    s.Finv( i_3, j_3 ) = Finv( i_3, j_3 );
    s.n( i_3 )         = n( i_3 );
    s.B( j_3, m_3 )    = Finv( j_3, k_3 ) * Finv( m_3, k_3 );

    // curvature, formulation with the bending strain (Dortdivanlioglu, Tomassetti)
    Tensor2_3d b;
    s.K( j_3, k_3 ) = -n( i_3 ) * G( i_3, j_3, k_3 ); // current curvature in the referential frame
    b( j_3, l_3 )   = -Finv( k_3, j_3 ) * s.K( k_3, m_3 ) * Finv( m_3, l_3 );
    s.H             = 0.5 * b( j_3, j_3 );
    s.H0            = this->giveH0( N, gradIt );

    // second order stress conjugate to K, M = de/dK
    s.M( i_3, j_3 ) = -0.5 * s.omega * s.J * ( s.H - s.H0 ) * s.B( i_3, j_3 );
    return s;
}


Tensor2_3d IsotropicPolyconvexHyperelasticSurfaceMaterial2::computeBendingStress( const Tensor2_3d &F, const Tensor2_3d &Is, const Tensor1_3d &N, StructuralSurfaceMaterialStatus *status, TimeStep *tStep ) const
{
    FloatArrayF<27> vG, vdI;
    giveSecondGradients( status, vG, vdI );
    Tensor3_3d G( vG ), gradIt( vdI );
    BendingState s = this->computeBendingState( F, Is, N, G, gradIt, tStep );

    // e = J omega / 2 (H - H0)^2; through J and F^+ at fixed K
    Tensor2_3d dedF1, dedF2, P2, P;
    double e          = s.J * ( 0.5 * s.omega * ( s.H - s.H0 ) * ( s.H - s.H0 ) );
    dedF1( i_3, j_3 ) = e * s.Finv( j_3, i_3 );
    dedF2( i_3, l_3 ) = s.J * s.omega * ( s.H - s.H0 ) * s.Finv( j_3, i_3 ) * s.K( j_3, k_3 ) * s.B( k_3, l_3 );
    // through the dependence of the current normal in K on F
    P2( p_3, q_3 ) = s.n( p_3 ) * s.Finv( q_3, k_3 ) * G( k_3, i_3, j_3 ) * s.M( i_3, j_3 );

    P( i_3, j_3 ) = dedF1( i_3, j_3 ) + dedF2( i_3, j_3 ) + P2( i_3, j_3 );
    return P;
}


FloatArrayF<9>
IsotropicPolyconvexHyperelasticSurfaceMaterial2::giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    Tensor2_3d F( vF ), Is = giveSurfaceIdentity( normal );
    FloatArrayF<3> vN( normal );
    Tensor1_3d N( vN );

    Tensor2_3d Pm = this->computeMembraneStress( F, Is, tStep );
    Tensor2_3d Pb = this->computeBendingStress( F, Is, N, status, tStep );
    Tensor2_3d P;
    P( i_3, j_3 ) = Pm( i_3, j_3 ) + Pb( i_3, j_3 );

    this->storeStress( gp, vF, F, P );
    return P.to_voigt_form();
}


FloatArrayF<27>
IsotropicPolyconvexHyperelasticSurfaceMaterial2::giveSecondOrderSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    Tensor2_3d Is = giveSurfaceIdentity( normal );
    FloatArrayF<3> vN( normal );
    Tensor1_3d N( vN );
    Tensor2_3d F( vF );
    Tensor3_3d A, G( vG ), gradIt( gradI );
    BendingState s = this->computeBendingState( F, Is, N, G, gradIt, tStep );

    // the third order stress conjugate to grad F
    A( i_3, j_3, k_3 ) = -s.n( i_3 ) * s.M( j_3, k_3 );
    auto vA            = A.to_voigt_form_27();

    // P and the tangents take F, grad F and grad I from the temporary status; M and K are kept for output
    auto status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    status->letTempFVectorBe( vF );
    status->letTempGVectorBe( vG );
    status->letTempdIVectorBe( gradI );
    status->letTempMVectorBe( s.M.to_voigt_form() );
    status->letTempKVectorBe( s.K.to_voigt_form() );

    return vA;
}


Tensor4_3d IsotropicPolyconvexHyperelasticSurfaceMaterial2::computeBendingTangent( const Tensor2_3d &F, const Tensor2_3d &Is, const Tensor1_3d &N, StructuralSurfaceMaterialStatus *status, TimeStep *tStep ) const
{
    FloatArrayF<27> vG, vdI;
    giveSecondGradients( status, vG, vdI );
    Tensor3_3d G( vG ), gradIt( vdI );
    BendingState s = this->computeBendingState( F, Is, N, G, gradIt, tStep );
    Tensor2_3d &Finv = s.Finv, &B = s.B, &K = s.K, &Mk = s.M;
    Tensor1_3d &n_def = s.n;
    const double J = s.J, H = s.H, H0 = s.H0, omega_t = s.omega;

    Tensor4_3d d2e_dKdK;
    d2e_dKdK( a_3, b_3, c_3, d_3 ) = 0.25 * J * omega_t * B( a_3, b_3 ) * B( c_3, d_3 );

    double e = J * ( 0.5 * omega_t * ( H - H0 ) * ( H - H0 ) );

    Tensor2_3d dHdF;
    dHdF( i_3, l_3 ) = Finv( j_3, i_3 ) * K( j_3, k_3 ) * B( k_3, l_3 );

    Tensor2_3d dedF1, dedF2, dedF;
    dedF1( i_3, j_3 ) = e * Finv( j_3, i_3 );
    dedF2( i_3, l_3 ) = J * omega_t * ( H - H0 ) * dHdF( i_3, l_3 );
    dedF( i_3, j_3 )  = dedF1( i_3, j_3 ) + dedF2( i_3, j_3 );

    // d2e/dFdF at fixed K
    Tensor4_3d d2H_dFdF, d2e_dFdF;
    d2H_dFdF( a_3, b_3, c_3, d_3 ) = ( -Finv( o_3, c_3 ) * Finv( d_3, a_3 ) + B( o_3, d_3 ) * n_def( a_3 ) * n_def( c_3 ) ) * K( o_3, p_3 ) * B( p_3, b_3 ) +
        Finv( o_3, a_3 ) * K( o_3, p_3 ) * ( -Finv( p_3, c_3 ) * B( b_3, d_3 ) - B( p_3, d_3 ) * Finv( b_3, c_3 ) );

    d2e_dFdF( a_3, b_3, c_3, d_3 ) = Finv( b_3, a_3 ) * dedF( c_3, d_3 ) +
        e * ( -Finv( b_3, c_3 ) * Finv( d_3, a_3 ) + B( b_3, d_3 ) * n_def( a_3 ) * n_def( c_3 ) ) +
        J * omega_t * dHdF( a_3, b_3 ) * dHdF( c_3, d_3 ) +
        J * omega_t * ( H - H0 ) * dHdF( a_3, b_3 ) * Finv( d_3, c_3 ) +
        J * omega_t * ( H - H0 ) * d2H_dFdF( a_3, b_3, c_3, d_3 );

    // mixed derivatives
    Tensor4_3d d2e_dKdF, d2e_dFdK;
    d2e_dKdF( a_3, b_3, c_3, d_3 ) = -0.5 * J * omega_t * ( H - H0 ) * B( a_3, b_3 ) * Finv( d_3, c_3 ) +
        0.5 * J * omega_t * ( H - H0 ) * ( Finv( a_3, c_3 ) * B( b_3, d_3 ) + B( a_3, d_3 ) * Finv( b_3, c_3 ) ) -
        0.5 * J * omega_t * B( a_3, b_3 ) * dHdF( c_3, d_3 );
    d2e_dFdK( a_3, b_3, c_3, d_3 ) = d2e_dKdF( c_3, d_3, a_3, b_3 );

    // K depends on F through the current normal: dK/dF = T
    Tensor4_3d T;
    T( i_3, j_3, p_3, q_3 ) = ( n_def( p_3 ) * Finv( q_3, k_3 ) ) * G( k_3, i_3, j_3 );
    Tensor6_3d dTdF;
    dTdF( i_3, j_3, k_3, l_3, m_3, n_3 ) = ( -Finv( l_3, m_3 ) * Finv( n_3, a_3 ) + B( l_3, n_3 ) * n_def( a_3 ) * n_def( m_3 ) ) * ( G( a_3, i_3, j_3 ) * n_def( k_3 ) ) -
        ( Finv( l_3, a_3 ) * Finv( n_3, k_3 ) ) * ( G( a_3, i_3, j_3 ) * n_def( m_3 ) );

    Tensor4_3d A;
    A( k_3, l_3, m_3, n_3 ) = d2e_dFdF( k_3, l_3, m_3, n_3 ) + d2e_dFdK( k_3, l_3, o_3, p_3 ) * T( o_3, p_3, m_3, n_3 ) +
        ( d2e_dKdF( i_3, j_3, m_3, n_3 ) + d2e_dKdK( i_3, j_3, o_3, p_3 ) * T( o_3, p_3, m_3, n_3 ) ) * T( i_3, j_3, k_3, l_3 ) +
        Mk( i_3, j_3 ) * dTdF( i_3, j_3, k_3, l_3, m_3, n_3 );
    return A;
}


FloatMatrixF<9, 9>
IsotropicPolyconvexHyperelasticSurfaceMaterial2::give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    Tensor2_3d F( vF ), Is = giveSurfaceIdentity( normal );
    FloatArrayF<3> vN( normal );
    Tensor1_3d N( vN );

    Tensor4_3d Am = this->computeMembraneTangent( F, Is, tStep );
    Tensor4_3d Ab = this->computeBendingTangent( F, Is, N, status, tStep );
    Tensor4_3d A;
    A( i_3, j_3, k_3, l_3 ) = Am( i_3, j_3, k_3, l_3 ) + Ab( i_3, j_3, k_3, l_3 );

    auto answer = A.to_voigt_form();
    if ( this->checkTangentTolerance != 0. ) {
        this->checkTangentsByFiniteDifferences( answer, normal, gp, tStep );
    }
    return answer;
}


FloatMatrixF<27, 27>
IsotropicPolyconvexHyperelasticSurfaceMaterial2::give3dSurfaceMaterialStiffnessMatrix_dAddF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    FloatArrayF<27> vG, vdI;
    giveSecondGradients( status, vG, vdI );
    Tensor2_3d F( vF ), Is = giveSurfaceIdentity( normal );
    FloatArrayF<3> vN( normal );
    Tensor1_3d N( vN );
    Tensor3_3d G( vG ), gradIt( vdI );
    BendingState s = this->computeBendingState( F, Is, N, G, gradIt, tStep );

    Tensor4_3d d2e_dKdK;
    d2e_dKdK( a_3, b_3, c_3, d_3 ) = 0.25 * s.J * s.omega * s.B( a_3, b_3 ) * s.B( c_3, d_3 );

    Tensor6_3d dA;
    dA( i_3, j_3, k_3, l_3, m_3, n_3 ) = ( s.n( i_3 ) * s.n( l_3 ) ) * d2e_dKdK( j_3, k_3, m_3, n_3 );
    return dA.to_voigt_form();
}


FloatMatrixF<27, 9>
IsotropicPolyconvexHyperelasticSurfaceMaterial2::give3dSurfaceMaterialStiffnessMatrix_dAdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const
{
    auto status = static_cast<StructuralSurfaceMaterialStatus *>( this->giveStatus( gp ) );
    FloatArrayF<9> vF( status->giveTempFVector() );
    FloatArrayF<27> vG, vdI;
    giveSecondGradients( status, vG, vdI );
    Tensor2_3d F( vF ), Is = giveSurfaceIdentity( normal );
    FloatArrayF<3> vN( normal );
    Tensor1_3d N( vN );
    Tensor3_3d G( vG ), gradIt( vdI );
    BendingState s = this->computeBendingState( F, Is, N, G, gradIt, tStep );
    Tensor2_3d &Finv = s.Finv, &B = s.B, &K = s.K, &Mk = s.M;
    Tensor1_3d &n_def = s.n;
    const double J = s.J, H = s.H, H0 = s.H0, omega_t = s.omega;

    Tensor4_3d dAK; // d2e/dK2
    dAK( a_3, b_3, c_3, d_3 ) = 0.25 * J * omega_t * B( a_3, b_3 ) * B( c_3, d_3 );

    Tensor2_3d dHdF;
    dHdF( i_3, l_3 ) = Finv( j_3, i_3 ) * K( j_3, k_3 ) * B( k_3, l_3 );

    Tensor4_3d T;
    T( i_3, j_3, p_3, q_3 ) = ( n_def( p_3 ) * Finv( q_3, k_3 ) ) * G( k_3, i_3, j_3 );

    Tensor4_3d d2e_dKdF;
    d2e_dKdF( a_3, b_3, c_3, d_3 ) = -0.5 * J * omega_t * ( H - H0 ) * B( a_3, b_3 ) * Finv( d_3, c_3 ) +
        0.5 * J * omega_t * ( H - H0 ) * ( Finv( a_3, c_3 ) * B( b_3, d_3 ) + B( a_3, d_3 ) * Finv( b_3, c_3 ) ) -
        0.5 * J * omega_t * B( a_3, b_3 ) * dHdF( c_3, d_3 );

    Tensor5_3d dA;
    dA( i_3, j_3, k_3, a_3, b_3 ) = ( Mk( j_3, k_3 ) * n_def( a_3 ) ) * Finv( b_3, i_3 ) -
        ( d2e_dKdF( j_3, k_3, a_3, b_3 ) + dAK( j_3, k_3, o_3, p_3 ) * T( o_3, p_3, a_3, b_3 ) ) * n_def( i_3 );

    return dA.to_voigt_form_27x9();
}

} // end namespace oofem
