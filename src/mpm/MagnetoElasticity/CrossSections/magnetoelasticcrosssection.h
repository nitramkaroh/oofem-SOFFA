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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#ifndef magnetoelasticcrosssection_h
#define magnetoelasticcrosssection_h

#include "crosssection.h"

#define _IFT_MagnetoElasticCrossSection_Name "magnetoelasticcrosssection"
#define _IFT_MagnetoElasticCrossSection_MaterialNumber "material"

namespace oofem {
class GaussPoint;
typedef GaussPoint IntegrationPoint;
class Element;
class FloatArray;
class FloatMatrix;

/**
 * Abstract base class for all structural cross section models. It declares commons services provided by all
 * structural cross section models. The implementation of this services is left on derived classes,
 * which will implement cross section model dependent part. However, some general services are
 * implemented here.
 * For information, how to introduce integration points in cross section volume for
 * macro integration point, see @ref CrossSection reference manual.
 *
 * At structural level of cross section or constitutive models are introduced several stress/strain modes.
 * Full and reduced formats of stress/strain vectors are also introduced for convenience.
 * The full format includes all components, even if they are zero due to stress/strain mode nature,
 * but in the reduced format, only generally nonzero components are stored.
 * (full format must used only if absolutely necessary, to avoid wasting of space. It is used
 * by output routines to print results in general form). Methods for converting vectors between
 * full and reduced format are provided.
 * General full strain vector has one of the following forms:
 * -# strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
 * -# For integrated cross section models (2d and 3d beams, plates and general shells)
 *    strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
 */
class OOFEM_EXPORT MagnetoElasticCrossSection : public CrossSection
{
public:
      int materialNumber = 0;   ///< Material number
    /**
     * Constructor. Creates cross section with given number, belonging to given domain.
     * @param n Cross section number.
     * @param d Domain to which new cross section will belong.
     */
     MagnetoElasticCrossSection(int n, Domain *d) : CrossSection(n, d)  { }

    /**
     * Computes the First Piola-Kirchoff stress vector for a given deformation gradient and integration point.
     * The service should use previously reached equilibrium history variables. Also
     * it should update temporary history variables in status according to newly reached state.
     * The temporary history variables are moved into equilibrium ones after global structure
     * equilibrium has been reached by iteration process.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains the First Piola-Kirchoff stresses.
     * @param gp Integration point.
     * @param reducedFIncrement Increment of the deformation gradient vector in reduced form. @todo should this then be in a multiplicative way? /JB
     * @param tStep Current time step (most models are able to respond only when tStep is current time step).
     */
    virtual std::vector<FloatArray> give_FluxVectors(const std::vector<FloatArray> &vGrads, GaussPoint *gp, TimeStep *tStep) const;
    //@{
    virtual std::tuple<FloatArrayF< 9 >,FloatArrayF< 3 >> give_FluxVectors_P_B_3d(const FloatArrayF< 9 > &reducedF, const FloatArrayF< 3 > &reducedH, GaussPoint *gp, TimeStep *tStep) const;
    virtual std::tuple<FloatArrayF< 5 >,FloatArrayF< 2 >> give_FluxVectors_P_B_PlaneStrain(const FloatArrayF< 5 > &reducedF, const FloatArrayF< 2 > &reducedH, GaussPoint *gp, TimeStep *tStep) const;
    //@}
    /**
     * Computes the large strain stiffness matrix of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equilibrium  history variables stored in integration point status
     * to compute and return required result.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
	virtual void give_dFlux_dGrad_Matrices(std::vector<FloatMatrix> &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the material stiffness matrix dPdF of receiver in a given integration point, respecting its history.
     * The algorithm should use temporary or equilibrium  history variables stored in integration point status
     * to compute and return required result.
     * Elements should always pass their requests to their cross section model, which
     * performs necessary integration over its volume and invokes necessary material
     * services for corresponding material model defined for given integration point.
     * @param answer Contains result.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    //@{
	
	virtual std::tuple<FloatMatrixF< 9, 9 >,FloatMatrixF< 9, 3 >,FloatMatrixF< 3, 9 >,FloatMatrixF< 3, 3 >> give_ConstitutiveMatrices_dPdF_dBdH_dPdH_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;
	virtual std::tuple<FloatMatrixF< 5, 5 >,FloatMatrixF< 5, 2 >,FloatMatrixF< 2, 5 >, FloatMatrixF< 2, 2 >> give_ConstitutiveMatrices_dPdF_dBdH_dPdH_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;


    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override{;}

    virtual void createMaterialStatus(GaussPoint &iGP);


    const char *giveClassName() const override { return "MagnetoElasticCrossSection"; }
    const char *giveInputRecordName() const override { return _IFT_MagnetoElasticCrossSection_Name; }

    int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)  override;
    int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *gp)  override;
    int estimatePackSize(DataStream &buff, GaussPoint *gp) override;

    Material *giveMaterial(IntegrationPoint *ip) const override;

    int giveMaterialNumber() const { return this->materialNumber; }
    void setMaterialNumber(int matNum) { this->materialNumber = matNum; }


    FloatArray give_SecondGradient_FluxVector(const FloatArray &vGrad, GaussPoint *gp, TimeStep *tStep) const;
    FloatArrayF< 27 > give_SecondGradient_FluxVector_3d(const FloatArrayF< 27 > &vG, GaussPoint *gp, TimeStep *tStep) const;
    FloatArrayF< 8 > give_SecondGradient_FluxVector_PlaneStrain(const FloatArrayF< 8 > &vG, GaussPoint *gp, TimeStep *tStep) const;
    void give_SecondGradient_dFlux_dGrad(FloatMatrix &answer, MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep);
    FloatMatrixF< 1, 1 > give_SecondGradient_ConstitutiveMatrix_3d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;
    FloatMatrixF< 1, 1 > give_SecondGradient_ConstitutiveMatrix_PlaneStrain(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;

    std::tuple<FloatArray, FloatArray> give_JacobianGradient_FluxVectors( const FloatArray &vF, const FloatArray &vGradF, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatArrayF<9>, FloatArrayF<27> > give_JacobianGradient_FluxVectors_3d( const FloatArrayF<9> &vF, const FloatArrayF<27> &vGradF, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatArrayF<5>, FloatArrayF<8> > give_JacobianGradient_FluxVectors_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatMatrix, FloatMatrix, FloatMatrix, FloatMatrix> give_JacobianGradient_dFluxes_dGrads( MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep );
    std::tuple<FloatMatrixF<9, 9>, FloatMatrixF<9, 27>, FloatMatrixF<27, 9>, FloatMatrixF<27, 27> > give_JacobianGradient_ConstitutiveMatrices_3d( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> > give_JacobianGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const;

    std::tuple<FloatArray, FloatArray> give_RotationGradient_FluxVectors( const FloatArray &vF, const FloatArray &vGradF, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatArrayF<5>, FloatArrayF<8> > give_RotationGradient_FluxVectors_PlaneStrain( const FloatArrayF<5> &vF, const FloatArrayF<8> &vGradF, GaussPoint *gp, TimeStep *tStep ) const;
    std::tuple<FloatMatrix, FloatMatrix, FloatMatrix, FloatMatrix> give_RotationGradient_dFluxes_dGrads( MatResponseMode rmode, GaussPoint *gp, TimeStep *tStep );
    std::tuple<FloatMatrixF<5, 5>, FloatMatrixF<5, 8>, FloatMatrixF<8, 5>, FloatMatrixF<8, 8> > give_RotationGradient_ConstitutiveMatrices_PlaneStrain( MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) const;
    

    //@}
};
  
} // end namespace oofem
#endif // magnetoelasticcrosssection_h
