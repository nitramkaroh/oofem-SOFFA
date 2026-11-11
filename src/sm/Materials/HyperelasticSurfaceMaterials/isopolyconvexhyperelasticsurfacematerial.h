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

#ifndef isopolyconvexhyperelasticsurfacematerial_h
#define isopolyconvexhyperelasticsurfacematerial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralsurfacems.h"
#include "sm/Materials/HyperelasticSurfaceMaterials/hyperelasticsurfacematerial.h"
#include "sm/Materials/arclengthmaterialinterface.h"

///@name Input fields for IsotropicPolyconvexHyperelasticSurfaceMaterial (membrane)
//@{
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_Name "surfacetensionmat"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gamma1 "gamma1"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha1 "alpha1"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha_s1 "alphas"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_lambda "lambda"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_mu "mu"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gammaLTF "gamma_ltf"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alphaLTF "alpha_ltf"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_alpha_sLTF "alphas_ltf"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_lambdaLTF "lambda_ltf"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_muLTF "mu_ltf"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_gaminit "gaminit"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_ArcLengthParameter "arclenghtpar"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_checkTangent "checktangent"
//@}

///@name Input fields for IsotropicPolyconvexHyperelasticSurfaceMaterial2 (membrane with bending)
//@{
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_Name "surfacetensionmat2"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_omega1 "omega1"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_omegaLTF "omega_ltf"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_use_H0 "useh0"
#define _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_H0 "h0"
//@}

namespace oofem {
/// A material parameter which follows its own time function; without it (ltf = 0) the value is constant.
struct SurfaceMaterialParameter {
    double value = 0.;
    int ltf      = 0;
};

/**
 * Hyperelastic surface material without bending stiffness (membrane). Energy per unit reference area
 * @f[
 * W(F) = \gamma J + \alpha |F| + \frac{1}{2}\alpha_s (J - 1)^2 + \mu \left( \frac{1}{2}|F|^2 - 1 - \ln J \right) + \frac{\lambda}{2} (\ln J)^2,
 * @f]
 * where F is the surface deformation gradient and J = |cof F N| the area ratio. The terms gamma J and alpha |F|
 * form the surface-polyconvex model of Horak, Smejkal, Kruzik, Surface-polyconvex models for soft elastic solids,
 * JMPS 2025 (eq. 90): they are positively one-homogeneous, which the existence theory needs. The terms with alphas,
 * mu and lambda are not, so they are an extension outside that theory.
 *
 * Kinematic contract: F (9 components) and the unit reference normal N are given in the same orthonormal frame with
 * F N = 0, and the surface identity is I_s = I - N x N. The elements ensure this (IGA surfaces: global frame,
 * Lagrange surfaces: the frame of the element, with prestrain the normal of the original configuration).
 *
 * Parameters (keyword, time function): gamma1 (or gamma) with gamma_ltf, alpha1 (or alpha) with alpha_ltf, alphas with
 * alphas_ltf, lambda with lambda_ltf, mu with mu_ltf. A parameter without its time function is constant; with
 * gamma_ltf, gamma(t) = gaminit + gamma1 f(t). With arclenghtpar 1, gamma is the load parameter of the material
 * arc-length method (matcalm); the time functions are not used then.
 * checktangent tol > 0 compares the tangents with central finite differences at every evaluation and stops when the
 * relative difference exceeds tol; tol < 0 only prints the differences.
 */
class IsotropicPolyconvexHyperelasticSurfaceMaterial : public HyperElasticSurfaceMaterial, public ArcLengthMaterialInterface
{
protected:
    SurfaceMaterialParameter gamma, alpha, alphaS, lambda, mu;
    double gaminit = 0.;

    /// Material parameter which is the load parameter of the material arc-length method, if any.
    enum CALM_parameter {
        calm_par_no,
        calm_par_gamma,
        calm_par_omega,
    };
    CALM_parameter arcLengthParameter = calm_par_no;

    double checkTangentTolerance = 0.;

public:
    IsotropicPolyconvexHyperelasticSurfaceMaterial( int n, Domain *d );

    void initializeFrom( InputRecord &ir ) override;

    FloatMatrixF<6, 6> give3dSurfaceMaterialStiffnessMatrix( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for large strains only" ); }
    FloatArrayF<6> giveRealSurfaceStressVector_3d( const FloatArrayF<6> &strain, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override { OOFEM_ERROR( "not implemented, this material is designed for large strains only" ); }

    FloatArrayF<9> giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatMatrixF<9, 9> give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;

    MaterialStatus *CreateStatus( GaussPoint *gp ) const override;
    Interface *giveInterface( InterfaceType t ) override;

    const char *giveInputRecordName() const override { return _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial_Name; }
    const char *giveClassName() const override { return "IsotropicPolyconvexHyperelasticSurfaceMaterial"; }

protected:
    /// Surface identity I - N x N of the unit reference normal N.
    static Tensor2_3d giveSurfaceIdentity( const FloatArray &normal );
    /// Value of a parameter at the time step, or its arc-length share when the arc-length method drives the material.
    double giveParameter( const SurfaceMaterialParameter &p, CALM_parameter role, TimeStep *tStep ) const;
    double giveGamma( TimeStep *tStep ) const;
    static void readParameter( InputRecord &ir, SurfaceMaterialParameter &p, const char *name, const char *ltfName );

    /// First Piola-Kirchhoff stress of the membrane energy.
    Tensor2_3d computeMembraneStress( const Tensor2_3d &F, const Tensor2_3d &Is, TimeStep *tStep ) const;
    /// dP/dF of the membrane energy.
    Tensor4_3d computeMembraneTangent( const Tensor2_3d &F, const Tensor2_3d &Is, TimeStep *tStep ) const;
    /// Stores F, P and the Cauchy stress in the temporary status.
    void storeStress( GaussPoint *gp, const FloatArrayF<9> &vF, const Tensor2_3d &F, const Tensor2_3d &P ) const;

    /**
     * Compares dP/dF (given) and, for materials with bending, dA/dF, dA/dgradF and dP/dgradF = (dA/dF)^T with central
     * finite differences of P and A in the directions which keep F N = 0. The temporary status is restored.
     */
    void checkTangentsByFiniteDifferences( const FloatMatrixF<9, 9> &dPdF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const;
};


/**
 * Membrane material with bending stiffness. Energy per unit reference area
 * @f[
 * W(F, \nabla F) = W_{membrane}(F) + \frac{1}{2} \omega J (H - H_0)^2,
 * @f]
 * where H = tr b / 2 is the current mean curvature, b = -F^{+T} K F^{+}, K = -n . grad F, n the current normal.
 * The bending term needs the second gradient, so only IGA surface elements of degree >= 2 accept this material.
 *
 * Parameters as the membrane, plus omega1 with omega_ltf and useh0 (0: H0 = 0, 1: the curvature of the reference
 * geometry, 2: the value h0). With arclenghtpar 2, omega is the load parameter of the material arc-length method.
 * P and the tangents take grad F and the gradient of the reference surface identity from the temporary status,
 * where the second order stress A leaves them, so A has to be evaluated first. Everything else (the moment M,
 * d2e/dK2) is computed again from these and the current omega, so it never comes from an older evaluation.
 */
class IsotropicPolyconvexHyperelasticSurfaceMaterial2 : public IsotropicPolyconvexHyperelasticSurfaceMaterial
{
protected:
    SurfaceMaterialParameter omega;
    int use_H0    = 0;
    double H0_val = 0.;

public:
    IsotropicPolyconvexHyperelasticSurfaceMaterial2( int n, Domain *d );

    void initializeFrom( InputRecord &ir ) override;

    FloatArrayF<9> giveFirstPKSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatArrayF<27> giveSecondOrderSurfaceStressVector_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vG, const FloatArrayF<27> &gradI, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;

    FloatMatrixF<9, 9> give3dSurfaceMaterialStiffnessMatrix_dPdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatMatrixF<27, 27> give3dSurfaceMaterialStiffnessMatrix_dAddF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;
    FloatMatrixF<27, 9> give3dSurfaceMaterialStiffnessMatrix_dAdF( MatResponseMode mode, const FloatArray &normal, GaussPoint *gp, TimeStep *tStep ) const override;

    const char *giveInputRecordName() const override { return _IFT_IsotropicPolyconvexHyperelasticSurfaceMaterial2_Name; }
    const char *giveClassName() const override { return "IsotropicPolyconvexHyperelasticSurfaceMaterial2"; }

    bool supportsSecondGradient() const override { return true; }

protected:
    /// Quantities of the bending energy e = omega J (H - H0)^2 / 2 at one point.
    struct BendingState {
        Tensor2_3d Finv, B, K, M; ///< pseudoinverse F^+, B = F^+ F^+T, K = -n . grad F, M = de/dK
        Tensor1_3d n = Tensor1_3d( FloatArrayF<3>{ 0., 0., 0. } ); ///< current normal
        double J = 0., H = 0., H0 = 0., omega = 0.;
    };

    double giveOmega( TimeStep *tStep ) const;
    /// Reference mean curvature H0 according to useh0; gradIt is the gradient of the reference surface identity.
    double giveH0( const Tensor1_3d &N, const Tensor3_3d &gradIt ) const;
    /// grad F and the gradient of the reference surface identity, which the second order stress left in the temporary status.
    static void giveSecondGradients( StructuralSurfaceMaterialStatus *status, FloatArrayF<27> &vG, FloatArrayF<27> &vdI );
    BendingState computeBendingState( const Tensor2_3d &F, const Tensor2_3d &Is, const Tensor1_3d &N, const Tensor3_3d &G, const Tensor3_3d &gradIt, TimeStep *tStep ) const;
    /// First Piola-Kirchhoff stress of the bending energy.
    Tensor2_3d computeBendingStress( const Tensor2_3d &F, const Tensor2_3d &Is, const Tensor1_3d &N, StructuralSurfaceMaterialStatus *status, TimeStep *tStep ) const;
    /// dP/dF of the bending energy.
    Tensor4_3d computeBendingTangent( const Tensor2_3d &F, const Tensor2_3d &Is, const Tensor1_3d &N, StructuralSurfaceMaterialStatus *status, TimeStep *tStep ) const;
};
} // end namespace oofem
#endif
