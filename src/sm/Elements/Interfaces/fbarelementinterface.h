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
 *               Copyright (C) 1993 - 2014    Borek Patzak
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

#ifndef fbarelementinterface_h
#define fbarelementinterface_h

#include "interface.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matresponsemode.h"
#include "nlstructuralelement.h"



namespace oofem {

class GaussPoint;
class TimeStep;

/**
 * F-bar element extension interface
 * @author Martin Horak, nitramkaroh@seznam.cz
 * @note Reference: EA de Souze Neto, D. Peric, M. Dutko, DRJ Owen: Design of simple low order finite elements for large strain analysis of nearly incompressible solids
 * The method is based on replacing of deformation gradient by its modified counterpart
 * \bar{F} = (\frac{J_0}{J})^{/frac{1}{3}}F, where J_0 is determinant of the deformation gradient evaluated at the element centroid
 * Appropriate to tackle incompressibility locking in large strains
 */

class FbarElementExtensionInterface : public Interface
{


public:
    /**
     * Constructor. Creates element interface belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    FbarElementExtensionInterface();
    /// Destructor.
    virtual ~FbarElementExtensionInterface() { }
  
  /**
   * Provides coordinates of element centroid. 
   */
  virtual void giveElementParametricCentroid(FloatArray &centroidCoords) = 0;


  /**
   * Computes the first Piola-Kirchhoff stress tensor in Voigt format for Fbar formulation.
   * It is slight modification of the original method. In the calculation of the 1PK stress,the deformation gradient is relaced by its Fbar counterpart and its it hen multiplied by (J/J0)^2/3
   * @param answer Computed stress vector in Voigt form.
   * @param gp Gauss point at which the stress is evaluated.
   * @param tStep Time step.
   */
  void computeFirstPKStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *ts, NLStructuralElement *e);               

  /**
   * Computes the modified deformation gradient Fbar in Voigt form at integration point ip and at time step tStep. 
  * @param gp Gauss point.
  * @param tStep Time step.
  * @param nlElement nlstructural element
  * @return tuple containing deformation gradient Fbar, F, and F0; Moreover, det(F), and det(F0)
  */
  std::tuple<FloatArrayF<9>,FloatArrayF<9>,FloatArrayF<9>, double, double> computeFbarDeformationGradientVector_3d(GaussPoint *gp, TimeStep *tStep, NLStructuralElement *nlElem);
  
  /**
   * Computes the stiffness matrix corrected due to F-bar formulation.
   * The response is evaluated using @f$ \int B_{\mathrm{H}}^{\mathrm{T}} A1 (B_{\mathrm{H}}) + B_{\mathrm{H}}^{\mathrm{T}} A2 (B_{\mathrm{H}})_0 \;\mathrm{d}v @f$, where
   * @f$ B_{\mathrm{H}} @f$ is the B-matrix which produces the displacement gradient vector @f$ H_{\mathrm{V}} @f$ when multiplied with
   * the solution vector a.
   * @f $ (B_{\mathrm{H}})_0 @f$ is the B-matrix which produces the displacement gradient vector in the element parametric centroid
   * @param answer Computed stiffness matrix correction.
   * @param rMode Response mode.
   * @param tStep Time step.
   */
  void computeFbarStiffnessMatrix_3d(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep, NLStructuralElement *elem);
  
   
  /**
   * Computes two material stiffness terms due to the F-bar formulation
   * @f$ \bar{A}^1_{ijmn} = (\frac{J_0}{J})^{1/3}(A_{ijmn} - \frac{1}{3}A_ijkl F_{nm}^{-1}F_{kl}) @f
   * @f$ \bar{A}^2_{ijmn} = \frac{1}{3}(\frac{J_0}{J})^{1/3}(A_ijkl (F_{nm})^{-1}_0F_{kl}) @f
   * @param answer Computed tranformed first elasticty stiffness.
   * @param rMode Response mode.
   * @param tStep Time step.
   */
  std::pair<FloatMatrixF<9,9>,FloatMatrixF<9,9>> comptuteFbarMaterialStiffnesses_3d(FloatMatrixF<9,9> mA, FloatArrayF<9> vF, FloatArrayF<9> vF0, FloatArrayF<9> vP);
  
};
} 
#endif //fbarelementinterface_h
