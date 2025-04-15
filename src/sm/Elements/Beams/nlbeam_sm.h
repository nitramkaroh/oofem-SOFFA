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

#ifndef nlbeam_sm_h
#define nlbeam_sm_h


#include "dofmanager.h"
#include "Elements/nlstructuralelement.h"
#include "scalarfunction.h"
#include "vtkxmlexportmodule.h"
#include "elementmatlabexportmoduleinterface.h"
#include "floatarrayf.h"


///@name Input fields for NlBeam_SM
//@{
#define _IFT_NlBeam_SM_NIP "nip"
#define _IFT_NlBeam_SM_Material "materialtype"
#define _IFT_NlBeam_SM_Beam_Tolerance "btol"
#define _IFT_NlBeam_SM_Beam_MaxIteration "bmaxit"
#define _IFT_NlBeam_SMBeam_NumberMaxSubsteps "nsubsteps"
#define _IFT_NlBeam_SM_Section_Tolerance "stol"
#define _IFT_NlBeam_SM_Section_MaxIteration "smaxit"

#define _IFT_NlBeam_SM_fx "fx"
#define _IFT_NlBeam_SM_fxLTF "fx_ltf"

#define _IFT_NlBeam_SM_fz "fz"
#define _IFT_NlBeam_SM_fzLTF "fz_ltf"

#define _IFT_NlBeam_SM_m "m"
#define _IFT_NlBeam_SM_mLTF "m_ltf"


#define _IFT_NlBeam_SM_pressure "pressure"
#define _IFT_NlBeam_SM_pressureLTF "pressure_ltf"


#define _IFT_NlBeam_SM_num "num"
#define _IFT_NlBeam_SM_sym "sym"



#define _IFT_NlBeam_Reissner_Name "nlbeam_reissner"
#define _IFT_NlBeam_Ziegler_Name  "nlbeam_ziegler"





//@}

namespace oofem {
  class NlBeamCrossSection;

/**
 * This class implements a 2-dimensional large strain beam element
 * The shooting method is used to calculate internal forces and stiffness matrix
 * Add more description 
 */
  class NlBeam_SM : public NLStructuralElement, public VTKXMLExportModuleElementInterface, public ElementMatlabExportModuleInterface
{
protected:
  int NIP = 100;
  double pitch = 10, beamLength = 0;
  FloatArray internalForces, fb_init;
  FloatArrayF<6> X;
  //FloatArray s, ds, u, w, phi, kappa;
  //FloatMatrix jacobi;
  double EA, EI, GAs;
  double beam_tol = 1.e-14, beam_maxit = 100;
  double section_tol = 1.e-6,section_maxit = 20;
  int nsubsteps_init = 4;
  //
  double Btildax, Btildaz;
  FloatArray FBrx, FBrz;
   // loadtime function associated with the given loading
  int Ba_ltf = 0;
  //
  FloatArray x, z, phi, vKappa, s;
  // angle of inclincation
  double alpha;
  // loading
  FloatArray Fx, Fz;
  double M;
  // loadtime function associated with the given loading
  int fx_ltf = 0, fz_ltf = 0, m_ltf = 0;
  // testing parameters
  bool num = false;
  bool sym = false;
  // storing internal forces for postrocessing purposes
  FloatArray vM, vQ, vN;
  // cross-section of the beam
  NlBeamCrossSection *crossSection;
 
public:
    NlBeam_SM(int n, Domain *aDomain);
    virtual ~NlBeam_SM(){;}

    void giveDofManDofIDMask(int inode, IntArray &) const override;

 
   void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
   void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;
   int computeNumberOfDofs() override { return 6; }
   void  printOutputAt(FILE *file, TimeStep *tStep) override{;}

  void initializeFrom(InputRecord &ir) override;
  void postInitialize() override;
  // composite type - so we can do the postprocessing in the next function
  Element_Geometry_Type giveGeometryType() const override { return EGT_Composite; }
  // divide the beam into some small subelements
    void giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;
  Interface *giveInterface(InterfaceType it) override;
  
  void updateYourself(TimeStep *tStep) override;
  void initForNewStep() override;

  /*******************************************************************
   Pure virtual functions of structural element class not needed here
  *********************************************************************/
  void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override {OOFEM_ERROR("Shouldn't be called");}
  void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override {OOFEM_ERROR("Shouldn't be called");}
  void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override {OOFEM_ERROR("Shouldn't be called");}
  void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override {OOFEM_ERROR("Shouldn't be called");}
  /////////////////////////////////////////////////////////////////////////
  void doMatlabOutput(std::string fileName) override;

protected:
  /**
     Material mode is beam 2d
   **/
  MaterialMode giveMaterialMode() override { return _2dBeam; }

  /**
   **/
  double computeLength() override;

  /**
   **/
  double computeAlpha();

  /**
     Finite diference integration of the three first order differential equations
   **/
  virtual std::tuple<FloatArrayF<3>,FloatMatrixF<3,3>, FloatArrayF<3>, double, FloatArrayF<4>>
  integrateAlongBeam(const FloatArray &fab, const FloatArray &ua, TimeStep *tStep) = 0;
  /**
     
   **/
  std::tuple<bool, FloatArrayF<3>, double> findLeftEndForces(const FloatArrayF<6> &x_target, TimeStep *tStep);

  /**
     Line search function
   **/
  double performLineSearch(const FloatArray &fa, const FloatArray &dforces, const FloatArray &xa, const  FloatArray &xb_target, TimeStep *tStep);
  
  /**
     Function related to constitutive behavior of the beam -- consider moving to another class
  **/
  double computeMomentFromCurvature(double kappa);
  double computeDerivativeOfMomentWrtCurvature(double kappa);
  double computeCurvatureFromMoment(double M);
  double computeCenterlineStrainFromInternalForces(double M, double N, double curvature);
  double computeDerStrainMoment(double M, double N, double curvature);
  double computeDerStrainNormalForce(double M, double N);

  /**
     Distributed load
   **/
  void computeLoadingVariables(double fx, double fz, double Brx, double Brz);
  double computeLoadAt(FloatArray &Px, FloatArray &Pz, double &tBax, double &tBaz, TimeStep *tStep);
  /**
     Write output file
   **/
  FILE * giveOutputStream(TimeStep *tStep);

  /**
     Debugging staff
   **/
  void computeStiffnessMatrix_num(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
  void giveInternalForcesVector_from_u(FloatArray &answer, TimeStep *tStep, const FloatArray &u);


};

  class NlBeam_Reissner : public NlBeam_SM
  {
  public:
    NlBeam_Reissner(int n, Domain *aDomain) : NlBeam_SM(n, aDomain)
    {;}
    virtual ~NlBeam_Reissner(){;}
 
  protected:
    const char *giveClassName() const override { return "NlBeam_Reissner"; }
    const char *giveInputRecordName() const override { return _IFT_NlBeam_Reissner_Name; }
    /**
       Finite diference integration of the three first order differential equations
    **/
    std::tuple<FloatArrayF<3>,FloatMatrixF<3,3>, FloatArrayF<3>, double, FloatArrayF<4>>
    integrateAlongBeam(const FloatArray &fab, const FloatArray &ua, TimeStep *tStep) override;
  };


  
  class NlBeam_Ziegler : public NlBeam_SM
  {
  public:
     NlBeam_Ziegler(int n, Domain *aDomain) : NlBeam_SM(n, aDomain)
    {;}
    virtual ~NlBeam_Ziegler(){;}
 
  private:
    //@todo: to be added to initialization
    double chi_tol = 1.e-10;
    int chi_iter = 30;
  protected:
    const char *giveClassName() const override { return "NlBeam_Ziegler"; }
    const char *giveInputRecordName() const override { return _IFT_NlBeam_Ziegler_Name; }

    /**
       Finite diference integration of the three first order differential equations
    **/
    std::tuple<FloatArrayF<3>,FloatMatrixF<3,3>, FloatArrayF<3>, double, FloatArrayF<4>>
    integrateAlongBeam(const FloatArray &fab, const FloatArray &ua, TimeStep *tStep) override;
  };
  
} // end namespace oofem
#endif // beam2d_h
