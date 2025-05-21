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

#ifndef nlbeamcs_h
#define nlbeamcs_h

#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/Materials/structuralmaterial.h"
#include "element.h"
#include "dictionary.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "interface.h"

///@name Input fields for NlBeamElasticCS
//@{

/////////////////////////////////////////////////////////////////
#define _IFT_NlBeamSimpleElasticCrossSection_Name "nlbeamelasticcs"
//////////////////////////////////////////////////////////////////
#define _IFT_NlBeamSimpleElastoPlasticCrossSection_Name "nlbeamelastoplasticcs"
//@}

namespace oofem {
  class NlBeamCrossSection;
/**
 * This class implements a .... 
 * Tasks:
 * - Returning material ....
 */

  
class NlBeamCrossSection : public StructuralCrossSection
{
protected:
  double E = 0.;
public:
    NlBeamCrossSection(int n, Domain *d) : StructuralCrossSection(n, d) { ; }
    virtual ~NlBeamCrossSection(){;}


  Material *giveMaterial(IntegrationPoint *ip) const override {return 0;}


  FloatArray *imposeStressConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStressVector3d) override{
      OOFEM_ERROR("not implemented");
      FloatArray *A;
      return A;
  }
  FloatArray *imposeStrainConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStressVector3d) override {
      OOFEM_ERROR("not implemented");
      FloatArray *A;
      return A;
  }

    // identification and auxiliary functions
  //    void initializeFrom(InputRecord &ir) override;
    // pure virtual functions 
    virtual double computeMomentFromCurvature(double kappa) = 0;
    virtual double computeDerivativeOfMomentWrtCurvature(double kappa) = 0;
    virtual double computeCurvatureFromMoment(double M) = 0;

  /************************************************************************************************************/
    FloatArrayF< 6 >giveRealStress_3d(const FloatArrayF< 6 > &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return zeros< 6 >();
    }

    FloatArrayF< 4 >giveRealStress_PlaneStrain(const FloatArrayF< 4 > &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return zeros< 4 >();
    }

    FloatArrayF< 3 >giveRealStress_PlaneStress(const FloatArrayF< 3 > &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return zeros< 3 >();
    }

    FloatArrayF< 1 >giveRealStress_1d(const FloatArrayF< 1 > &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return zeros< 1 >();
    }

    FloatArrayF< 2 >giveRealStress_Warping(const FloatArrayF< 2 > &reducedStrain, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return zeros< 2 >();
    }



    FloatMatrixF< 6, 6 >giveStiffnessMatrix_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 6 ,6 >();
    }

    FloatMatrixF< 3, 3 >giveStiffnessMatrix_PlaneStress(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 3 ,3 >();

    }

    FloatMatrixF< 4, 4 >giveStiffnessMatrix_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 4 ,4 >();
    }

    FloatMatrixF< 1, 1 >giveStiffnessMatrix_1d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 1 ,1 >();
    }


    FloatArrayF< 3 >giveGeneralizedStress_Beam2d(const FloatArrayF< 3 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return zeros< 3 >();
    }

    FloatArrayF< 6 >giveGeneralizedStress_Beam3d(const FloatArrayF< 6 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return zeros< 6 >();
    }

    FloatArrayF< 5 >giveGeneralizedStress_Plate(const FloatArrayF< 5 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return zeros< 5 >();
    }

    FloatArrayF< 8 >giveGeneralizedStress_Shell(const FloatArrayF< 8 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return zeros< 8 >();
    }

    FloatArrayF< 9 >giveGeneralizedStress_ShellRot(const FloatArrayF< 9 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return zeros< 9 >();
    }

    FloatArrayF< 4 >giveGeneralizedStress_MembraneRot(const FloatArrayF< 4 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return zeros< 4 >();
    }

    FloatArrayF< 3 >giveGeneralizedStress_PlateSubSoil(const FloatArrayF< 3 > &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return zeros< 3 >();

    }


    void giveCharMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) override
          {
        OOFEM_ERROR("not implemented");
    }


    FloatMatrixF< 3, 3 >give2dBeamStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 3, 3 >();
    }

    FloatMatrixF< 6, 6 >give3dBeamStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 6, 6 >();
    }

    FloatMatrixF< 5, 5 >give2dPlateStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 5, 5 >();
    }

    FloatMatrixF< 8, 8 >give3dShellStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 8, 8 >();
    }

    FloatMatrixF< 9, 9 >give3dShellRotStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 9, 9 >();
    }

    FloatMatrixF< 4, 4 >giveMembraneRotStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 4, 4 >();
    }

    FloatMatrixF< 3, 3 >give2dPlateSubSoilStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 3, 3 >();
    }


    bool isCharacteristicMtrxSymmetric(MatResponseMode mode) const override
          {
        OOFEM_ERROR("not implemented");
        return false;
    }

    double give(int aProperty, GaussPoint *gp) const override
    {
        OOFEM_ERROR("not implemented yet");
        return 0;
    }

    FloatArrayF< 9 >giveFirstPKStress_3d(const FloatArrayF< 9 > &reducedvF, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return zeros< 9 >();
    }

    FloatArrayF< 5 >giveFirstPKStress_PlaneStrain(const FloatArrayF< 5 > &reducedvF, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return zeros< 5 >();
    }

    FloatArrayF< 4 >giveFirstPKStress_PlaneStress(const FloatArrayF< 4 > &reducedvF, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return zeros< 4 >();
    }

    FloatArrayF< 1 >giveFirstPKStress_1d(const FloatArrayF< 1 > &reducedvF, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return zeros< 1 >();
    }




    void giveCharMaterialStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override
          {
        OOFEM_ERROR("not implemented");
        //return 0;
    }

    FloatMatrixF< 9, 9 >giveStiffnessMatrix_dPdF_3d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 9, 9 >();
    }

    FloatMatrixF< 5, 5 >giveStiffnessMatrix_dPdF_PlaneStrain(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 5, 5 >();
    }

    FloatMatrixF< 4, 4 >giveStiffnessMatrix_dPdF_PlaneStress(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override
          {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 4, 4 >();
    }

    FloatMatrixF< 1, 1 >giveStiffnessMatrix_dPdF_1d(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override    {
        OOFEM_ERROR("not implemented");
	return FloatMatrixF< 1, 1 >();
    }
    void giveCauchyStresses(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedFIncrement, TimeStep *tStep) override
    { OOFEM_ERROR("not implemented"); }
    void giveStiffnessMatrix_dCde(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override
    { OOFEM_ERROR("not implemented"); }
  // what should this function do?
  int checkConsistency() override {return 1;}
  int packUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override {return 0;}
  int estimatePackSize(DataStream &buff, GaussPoint *ip) override {return  0;}
  int unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep, GaussPoint *ip) override {return 0;}
  void createMaterialStatus(GaussPoint &iGP) override{;}
  int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override{return 0;}

  void setE(double E){this->E = E;}
  virtual double giveEA(){ return 0;}
  virtual double giveGAs(){return 0;}
  virtual double giveEI(){ return 0;}


  


};



  
class NlBeamSimpleElasticCrossSection : public NlBeamCrossSection
{
protected:
    double Iy = 0.; ///< description of the paramter.
    double A = 0.; ///< description of the paramter.
    double gamma = 5./6.;
  
public:
    NlBeamSimpleElasticCrossSection(int n, Domain *d) : NlBeamCrossSection(n, d)
    { }
    // identification and auxiliary functions
    const char *giveInputRecordName() const override { return _IFT_NlBeamSimpleElasticCrossSection_Name; }
    const char *giveClassName() const override { return "NlBeamElasticCrossSection"; }
    void initializeFrom(InputRecord &ir) override;
    ///
    double computeMomentFromCurvature(double kappa) override;
    double computeDerivativeOfMomentWrtCurvature(double kappa)  override;
    double computeCurvatureFromMoment(double M)  override;
    void postInitialize();
    double giveEA() override {return E*A;}
    double giveEI() override {return E*Iy;}
    double giveGAs() override {return gamma * E*A;}


};

class NlBeamSimpleElastoPlasticCrossSection : public NlBeamSimpleElasticCrossSection
{
protected:
  double kappa_elastic; ///< description of the paramter.
  double delta_kappa;   ///< description of the paramter.
  double kappa_nl;     ///< description of the paramter. 
  double M_nl;        ///< description of the paramter. 
  double M_elastic;   ///< description of the paramter.
  double M_plastic;   ///< description of the paramter.
  
  double delta = 0.1;

public:
    NlBeamSimpleElastoPlasticCrossSection(int n, Domain *d) : NlBeamSimpleElasticCrossSection(n, d)
    { }

    // identification and auxiliary functions
    const char *giveInputRecordName() const override { return _IFT_NlBeamSimpleElastoPlasticCrossSection_Name; }
    const char *giveClassName() const override { return "NlBeamElastoPlasticCrossSection"; }
    void initializeFrom(InputRecord &ir) override;
    ///
    double computeMomentFromCurvature(double kappa) override;
    double computeDerivativeOfMomentWrtCurvature(double kappa)  override;
    double computeCurvatureFromMoment(double M)  override;


  void postInitialize();

  
};

  



} // end namespace oofem
#endif // nlbeamelasticcs_h
