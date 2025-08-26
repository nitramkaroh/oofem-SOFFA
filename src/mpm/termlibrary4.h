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

#include "mpm.h"
#include "boundaryload.h"
#include "termlibrary.h"

// file containing various term definitions 

namespace oofem {

 /**
 * @brief A Linear momentum balance equation term ($B^T\P(F, H)$)
 * 
 */
  class MagnetoElasticity_GradGrad_Term : public Term {
    protected:
    const Variable& magneticPotentialField;
    public:
    MagnetoElasticity_GradGrad_Term (const Variable &testField, const Variable& displacementField, const Variable& magneticPotentialField) ;
    /**
     * @brief Evaluates Internal forces vector, i.e. $B^TP(F, H)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
  void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;

    void getDimensions(Element& cell) const override{;}
    void initializeCell(Element& cell) const override{;}
    int computeGradientFields (std::vector<FloatArray> &grad, std::vector<FloatMatrix> &B, MPElement& cell, const FloatArray& lcoords, MaterialMode mmode, TimeStep* tstep) const;

protected:
  void computeBHmatrixAt(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) const;

    //@todo: delete
    void compute_lin_num(FloatMatrix &D1, FloatMatrix &D2, FloatMatrix &D3, FloatMatrix &D4, MPElement &cell, GaussPoint *gp, TimeStep* tstep) const;
    int computeGradientFields_from_u(std::vector<FloatArray>& vGrad, std::vector<FloatMatrix>& vB, MPElement& cell, const FloatArray& lcoords, MaterialMode mmode, TimeStep* tstep, const FloatArray &u, const FloatArray &phi) const;
};



 /**
 * @brief A Linear momentum balance equation term ($B^T\P(F, H)$)
 * 
 */
  class MagnetoElasticity_GradGrad_SecondGradientTerm : public Term {
    
    public:
    MagnetoElasticity_GradGrad_SecondGradientTerm (const Variable &testField, const Variable& displacementField) ;
    /**
     * @brief Evaluates Internal forces vector, i.e. $B^TP(F, H)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;

    void getDimensions(Element& cell) const override{;}
    void initializeCell(Element& cell) const override{;}
    int computeSecondGradientField (FloatArray &grad, FloatMatrix &B, MPElement& cell, const FloatArray& lcoords, MaterialMode mmode, TimeStep* tstep) const;

protected:
  void computeGmatrixAt(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) const;
};

/**
 * @brief 
 *
 */
class MagnetoElasticity_GradGrad_FirstSecondGradientTerm : public Term
{

  public:
  MagnetoElasticity_GradGrad_FirstSecondGradientTerm( const Variable &testField, const Variable &displacementField );
  /**
   * @brief Evaluates Internal forces vector, i.e. $B^TP(F, H)$
   *
   * @param cell
   * @param coords
   */
  void evaluate( FloatArray &, MPElement &cell, GaussPoint *gp, TimeStep *tstep ) const override;
  void evaluate_lin( FloatMatrix &answer, MPElement &e, GaussPoint *gp, TimeStep *tstep ) const override;

  void getDimensions( Element &cell ) const override { ; }
  void initializeCell( Element &cell ) const override { ; }
  int computeGradientField( FloatArray &grad, FloatMatrix &B, MPElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const;
  int computeSecondGradientField( FloatArray &grad, FloatMatrix &B, MPElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const;

  protected:
  void computeBHmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const;
  void computeGmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const;
};

/**
 * @brief Evaluates term dependent on F and Fbar (= F at cell center)
 *
 */
class MagnetoElasticity_GradGrad_FbarTerm : public Term
{
  protected:
  public:
  MagnetoElasticity_GradGrad_FbarTerm( const Variable &testField, const Variable &displacementField );
  /**
   * @brief Evaluates Internal forces vector, i.e. $B^TP(F, H)$
   *
   * @param cell
   * @param coords
   */
  void evaluate( FloatArray &, MPElement &cell, GaussPoint *gp, TimeStep *tstep ) const override;
  void evaluate_lin( FloatMatrix &answer, MPElement &e, GaussPoint *gp, TimeStep *tstep ) const override;

  void getDimensions( Element &cell ) const override { ; }
  void initializeCell( Element &cell ) const override { ; }
  int computeGradientField( FloatArray &grad, FloatMatrix &B, MPElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const;
  
  protected:
  void computeBHmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const;

};

} // end namespace oofem
#endif // termlibrary4_h
