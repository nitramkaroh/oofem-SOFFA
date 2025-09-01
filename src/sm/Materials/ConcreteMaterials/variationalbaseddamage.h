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

#ifndef variationalbaseddamage_h
#define variationalbaseddamage_h

#include "Materials/ConcreteMaterials/idm1.h"
#include "../sm/Materials/graddamagematerialextensioninterface.h"

#define _IFT_VarBasedDamageMaterial_Name "varbaseddamagematerial"
#define _IFT_VarBasedDamageMaterial_initDamage "initdamage"
#define _IFT_VarBasedDamageMaterial_beta "beta"
#define _IFT_VarBasedDamageMaterial_p "p"
#define _IFT_VarBasedDamageMaterial_phaseFieldModelType "phasefieldmodeltype"

#define _IFT_VarBasedDamageMaterial_a2 "a2"
#define _IFT_VarBasedDamageMaterial_a3 "a3"

#define _IFT_VarBasedDamageMaterial_Gf "fractenergy"
#define _IFT_VarBasedDamageMaterial_LdInf "ldinf"
#define _IFT_VarBasedDamageMaterial_ft "ft"
#define _IFT_VarBasedDamageMaterial_youngsModulus "youngsmodulus"

#define _IFT_VarBasedDamageMaterial_wu_softening_law "wusoftlaw"

#define _IFT_VarBasedDamageMaterial_equivstraintype "equivstraintype"
#define _IFT_VarBasedDamageMaterial_damageLaw "damlaw"
#define _IFT_VarBasedDamageMaterial_gf "gf"

namespace oofem {
/**
 * Variationally-based Gradient Isotropic Damage models,
 * Should it be inherited form IDM1?
 */
class VarBasedDamageMaterial : public IsotropicDamageMaterial1, public GradientDamageMaterialExtensionInterface
{
protected:
  double initialDamage;
  double beta;
  double p;
  double penalty; //???
  int pf; //???
  double a1; // Wu model - coefficient for approximation polynomial Q(gamma) in damage function 
  double a2;  // Wu model - coefficient for approximation polynomial Q(gamma) in damage function
  double a3;  // Wu model - coefficient for approximation polynomial Q(gamma) in damage function
  double Gf; // Wu model - fracture energy
  double LdInf; //Wu model - the size of the damage zone at complete failure
  double ft; //Wu model - tensile strength of the material
  double gfInf;
  double internalLengthInf;
  double maxDamage = 0.995;
  //
  enum PhaseFieldModelType {
			    phaseFieldModel_JZ=0,
			    phaseFieldModel_Miehe=1,
			    phaseFieldModel_Wu=2,
  }; // type of phase-field model

  PhaseFieldModelType phaseFieldModelType;

  enum WuSofteningLaw {
		     linear_softening=0,
		     exponential_softening=1,
		     user_specified_softening=2,
  }; // softening law for Wu phase-field model, user_specified_softening requires setting parameters a1, a2, a3

  WuSofteningLaw wuSofteningLaw;


public:
    /// Constructor
    VarBasedDamageMaterial(int n, Domain * d);
    /// Destructor
    virtual ~VarBasedDamageMaterial();

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
    // identification and auxiliary functions
    const char *giveClassName() const override { return "VarBasedDamageMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_VarBasedDamageMaterial_Name; }
    void initializeFrom(InputRecord &ir) override;


    Interface *giveInterface(InterfaceType t) override {
        if ( t == GradientDamageMaterialExtensionInterfaceType ) {
            return static_cast< GradientDamageMaterialExtensionInterface * >(this);
        } else {
            return NULL;
        }
    }
    virtual int hasMaterialModeCapability(MaterialMode mode);

    void giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
   void giveGradientDamageStiffnessMatrix_dd_NN(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
     void giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override;
    void giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep) override;

    void giveStiffnessMatrix(FloatMatrix &answer,  MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    
    void computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep) override;

    void giveNonlocalInternalForces_N_factor(double &answer, double nlddv, GaussPoint *gp, TimeStep *tStep)  override;
    void giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlddv, GaussPoint *gp, TimeStep *tStep)  override;
    
    
    virtual void computeDamage(double &answer, double damageDrivingVariable, GaussPoint *gp);
    virtual void computeDamagePrime(double &answer,  double damageDrivingVariable, GaussPoint *gp) const;
    virtual void computeDamagePrime2(double &answer, double damageDrivingVariable, GaussPoint *gp);
    double solveExpLaw(double dam, double c) const;
    virtual void computeDissipationFunctionPrime(double &answer, double damage, double damageDrivingVariable, GaussPoint *gp) const;
    virtual void computeDissipationFunctionPrime2(double &answer, double damage, double damageDrivingVariable, GaussPoint *gp);
#ifdef keep_track_of_dissipated_energy
    virtual void computeRegulirizingWork(GaussPoint *gp,const FloatArray &nonlocalDamageDrivingVariableGrad);
#endif

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    virtual double computeQ_Wu(double damageDrivingVariable) const;

    virtual double computeQ_Wu_prime1(double damageDrivingVariable) const;

    virtual double computeQ_Wu_prime2(double damageDrivingVariable);

    virtual double compute_dissipation_Wu_prime1_in_gamma(double damageDrivingVariable) const;

    virtual double compute_dissipation_Wu_prime2_in_gamma(double damageDrivingVariable);
     

    
};


/**
 *  * Variatinally-based Gradient Isotropic Damage models,
 */
class VarBasedDamageMaterialStatus : public IsotropicDamageMaterial1Status, public GradientDamageMaterialStatusExtensionInterface
{
 protected:
    FloatArray effectiveStressVector;
    FloatArray tempEffectiveStressVector;
    double strainEnergy, tempStrainEnergy;
    double tempRegularizingEnergy, regularizingEnergy;

 public:
  VarBasedDamageMaterialStatus(GaussPoint * g, double initialDamage, double localDamageDrivingVariable);
    virtual ~VarBasedDamageMaterialStatus();

    const char *giveClassName() const override { return "VarBasedDamageMaterialStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    virtual void letTempEffectiveStressVectorBe(const FloatArray &effectiveStress){this->tempEffectiveStressVector = effectiveStress;}
    const FloatArray &giveTempEffectiveStressVector() const { return tempEffectiveStressVector; }
	
    double giveStrainEnergy(){return strainEnergy;}
    void setTempStrainEnergy(double sE){strainEnergy = sE;}
    void setTempRegularizingEnergy(double tempReg){tempRegularizingEnergy = tempReg;}
    double giveRegularizingEnergy(){return regularizingEnergy;}
};
} // end namespace oofem
#endif // variationalbaseddamage_h
