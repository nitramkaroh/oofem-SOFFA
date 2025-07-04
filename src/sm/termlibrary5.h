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
 *               Copyright (C) 1993 - 2024   Borek Patzak
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
#ifndef termlibrary4_h
#define termlibrary4_h

#include "boundaryload.h"
#include "termlibrary.h"
#include "structuralelement.h"

// file containing various term definitions

namespace oofem {


/**
 * @brief Class representing a weak form expression to be evaluated (integrated).
 * It defines two key methods:
 * - evaluate method to evaluate the term value, all unknowns known
 * - evaluate_lin to evaluate the consistent linearization of the term, so if Term is T(u), depending on unknown u,
 *   this term evaluates dT/du, which typically contributes to the LHS.
 */
class StructuralTerm
{
  public:
  const Variable &field;
  const Variable &testField;
  MaterialMode mode;

  public:
  StructuralTerm( const Variable &testField, const Variable &unknownField, MaterialMode m = MaterialMode::_Unknown ) :
      field( unknownField ), testField( testField ) { mode = m; }

  // evaluate linearized term contribution to weak form on given cell at given point
  virtual void evaluate_lin( FloatMatrix &, StructuralElement &cell, GaussPoint *gp, TimeStep *tStep ) const = 0;
  // evaluate contribution (all vars known) on given cell
  virtual void evaluate( FloatArray &, StructuralElement &cell, GaussPoint *gp, TimeStep *tStep ) const = 0;
  virtual void getDimensions( Element &cell ) const = 0;
  virtual void initializeCell( Element &cell ) const = 0;
  virtual IntegrationRule *giveElementIntegrationRule( Element *e ) const { return NULL; };
};


/**
 * @brief A Linear momentum balance equation term ($B^T\P(F, H)$)
 *
 */
class ThirdMedium_GradGrad_SecondGradientTerm : public StructuralTerm
{
  protected:
  public:
  ThirdMedium_GradGrad_SecondGradientTerm( const Variable &testField, const Variable &displacementField );
  /**
   * @brief Evaluates Internal forces vector, i.e. $B^TP(F, H)$
   *
   * @param cell
   * @param coords
   */
  void evaluate( FloatArray &, StructuralElement &cell, GaussPoint *gp, TimeStep *tstep ) const override;
  void evaluate_lin( FloatMatrix &answer, StructuralElement &e, GaussPoint *gp, TimeStep *tstep ) const override;

  void getDimensions( Element &cell ) const override { ; }
  void initializeCell( Element &cell ) const override { ; }
  int computeSecondGradientField( FloatArray &grad, FloatMatrix &B, StructuralElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const;

  protected:
  void computeGmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const;
};

/**
 * @brief A Linear momentum balance equation term ($B^T\P(F, H)$)
 *
 */
class ThirdMedium_GradGrad_JacobianGradientTerm : public StructuralTerm
{
  protected:
  public:
  ThirdMedium_GradGrad_JacobianGradientTerm( const Variable &testField, const Variable &displacementField );
  /**
   * @brief Evaluates Internal forces vector, i.e. $B^TP(F, H)$
   *
   * @param cell
   * @param coords
   */
  void evaluate( FloatArray &, StructuralElement &cell, GaussPoint *gp, TimeStep *tstep ) const override;
  void evaluate_lin( FloatMatrix &answer, StructuralElement &e, GaussPoint *gp, TimeStep *tstep ) const override;

  void getDimensions( Element &cell ) const override { ; }
  void initializeCell( Element &cell ) const override { ; }
  int computeGradientField( FloatArray &grad, FloatMatrix &B, StructuralElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const;
  int computeSecondGradientField( FloatArray &grad, FloatMatrix &B, StructuralElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const;

  protected:
  void computeBHmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const;
  void computeGmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const;
};

} // end namespace oofem
#endif // termlibrary4_h
