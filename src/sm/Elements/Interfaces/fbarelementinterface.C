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
 *               Copyright (C) 1993 - 2023   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */


#include "fbarelementinterface.h"
#include "domain.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "gausspoint.h"		  
#include "engngm.h"
#include "../sm/Materials/structuralms.h"
#include "mathfem.h"
#include "tensor/tensor2.h"
#include "tensor/tensor4.h"

namespace oofem {


// constructor
FbarElementExtensionInterface :: FbarElementExtensionInterface()  : Interface()
{
}


  std::tuple<FloatArrayF<9>,FloatArrayF<9>,FloatArrayF<9>, double, double> 
FbarElementExtensionInterface :: computeFbarDeformationGradientVector_3d(GaussPoint *gp, TimeStep *tStep, NLStructuralElement *nlElem)
{
  // create a fictitious integration point in the element centroid
  FloatArray centroidCoords;
  auto matMode = nlElem->giveMaterialMode();
  this->giveElementParametricCentroid(centroidCoords);
  auto gp0 = std::make_unique< GaussPoint >( nullptr, 1, centroidCoords, 1.0, matMode);

  // compute compatible deformation gradient
  FloatArray  vF;
  FloatMatrix mF;
  nlElem->computeDeformationGradientVector(vF, gp, tStep);
  mF.beMatrixForm(vF);
  // compute deformation gradient in the element centroid
  FloatArray vF0;
  FloatMatrix mF0;
  nlElem->computeDeformationGradientVector(vF0, gp0.get(), tStep);
  mF0.beMatrixForm(vF0);
  // compute Fbar
  auto J = mF.giveDeterminant();
  auto J0 = mF0.giveDeterminant();
  FloatMatrixF<3,3> Fbar(mF), F(mF), F0(mF0);
  Fbar *= pow( J0 / J, 1. / 3. );
  
  return {to_voigt_form(Fbar), to_voigt_form(F), to_voigt_form(F0), J, J0};
}


void
FbarElementExtensionInterface :: computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, NLStructuralElement *elem)
{
    // Computes the first Piola-Kirchhoff stress vector containing the stresses at the  Gauss point
    // gp of the receiver at time step tStep. The deformation gradient F is computed and sent as
    // input to the material model.
    StructuralCrossSection *cs = elem->giveStructuralCrossSection();

    auto [vFbar, vF, vF0, J, J0] = this->computeFbarDeformationGradientVector_3d(gp, tStep, elem);
    answer = cs->giveFirstPKStresses(vFbar, gp, tStep);
    answer.times( pow( J/J0, 2. / 3. ) );
}



  
void
FbarElementExtensionInterface :: computeFbarStiffnessMatrix_3d(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep, NLStructuralElement *elem)
{
    FloatMatrix B, B0, A1B, A2B0, dPdF;
    answer.clear();

    // create a fictitious integration point in the element centroid
    GaussPoint *gp0;
    FloatArray centroidCoords, vP;
    MaterialMode matMode = elem->giveMaterialMode();
    this->giveElementParametricCentroid(centroidCoords);
    gp0 = new GaussPoint( NULL, 1, centroidCoords, 1.0, matMode);
    // loop over gauss points
    for ( auto &gp : *elem->giveDefaultIntegrationRulePtr() ) {
      if ( elem->giveGeometryMode() == 1 ) {
	elem->computeBHmatrixAt(gp, B);
	elem->computeBHmatrixAt(gp0,B0);
	if ( elem->giveDomain()->giveEngngModel()->giveFormulation() == TL ) {
	  // Material stiffness dP/dF
	  elem->computeConstitutiveMatrix_dPdF_At(dPdF, rMode, gp, tStep);
	  auto [vFbar, vF, vF0, J, J0] = this->computeFbarDeformationGradientVector_3d(gp, tStep, elem);
	  
	  StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
		  //	  if ( useUpdatedGpRecord == 1 ) {
	  vP = matStat->giveTempPVector();
	    /*	  } else {
	    
		  }*/
	  //this->computeFirstPKStressVector(vP, gp, tStep, elem);
	  //	  
	  auto [A1, A2] = comptuteFbarMaterialStiffnesses_3d(dPdF, vF, vF0, vP);
	  double dV = elem->computeVolumeAround(gp);
	  A1B.beProductOf(A1, B);
	  A2B0.beProductOf(A2, B0);
	  // Stiffness matrix for F-bar formulation is unsymmetric !!!
	  answer.plusProductUnsym(B, A1B, dV);
	  answer.plusProductUnsym(B, A2B0, dV); 
	}
      }
    }
}

std::pair<FloatMatrixF<9,9>,FloatMatrixF<9,9>>
FbarElementExtensionInterface :: comptuteFbarMaterialStiffnesses_3d(FloatMatrixF<9,9> mA, FloatArrayF<9> vF, FloatArrayF<9> vF0, FloatArrayF<9> vP)
{
  Tensor4_3d A(mA), t_answer_1, t_answer_2;
  Tensor2_3d F(vF), F0(vF0), P(vP);
  auto [J,iF] = F.compute_determinant_and_inverse();
  auto [J0,iF0] = F0.compute_determinant_and_inverse();
  //
  t_answer_1(i_3, j_3, k_3, l_3) = pow( J / J0, 1. /3. ) * A( i_3, j_3, k_3, l_3 ) - 1. / 3. * pow( J / J0, 1. / 3.) * A(i_3, j_3, m_3, n_3) * F(m_3,n_3) * iF(l_3,k_3) + 2. / 3. * pow( J / J0, 2. / 3.) * P(i_3,j_3) * iF( l_3, k_3 );
  t_answer_2(i_3, j_3, k_3, l_3) =  1. / 3. * pow( J / J0, 1. / 3.) * A(i_3, j_3, m_3, n_3) * F(m_3,n_3) * iF0(l_3,k_3) - 2. / 3. * pow( J / J0, 2. / 3.) * P(i_3,j_3) * iF0( l_3, k_3 );
  /*
  t_answer_1(i_3, j_3, k_3, l_3) = pow( J0 / J, 1. / 3. ) * ( A( i_3, j_3, k_3, l_3 ) - 1. / 3. * A(i_3, j_3, m_3, n_3) * F(m_3,n_3) * iF(l_3,k_3) );
  t_answer_2(i_3, j_3, k_3, l_3) =  1. / 3. * pow( J0 / J, 1. / 3.) * A(i_3, j_3, m_3, n_3) * F(m_3,n_3) * iF0(l_3,k_3);
  */
  //
  auto m_answer_1 = t_answer_1.to_voigt_form();
  auto m_answer_2 = t_answer_2.to_voigt_form();
  return { m_answer_1, m_answer_2 };
}


 

} // end namespace oofem
