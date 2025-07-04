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


#include "termlibrary5.h"
#include "termlibrary.h"
#include "element.h"
#include "../nlstructuralelement.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "intarray.h"
#include "classfactory.h"
#include "gaussintegrationrule.h"
#include "crosssection.h"
#include "sm/CrossSections/structuralcrosssection.h"

#include "fei2dquadquad.h"
#include "mathfem.h"
#include "floatarrayf.h"
#include "lobattoir.h"

namespace oofem {

/**
 * @brief Base class for third medium elements (with local second gradient terms)
 *
 */
class ThirdMediumElement : public NLStructuralElement
{

  private:
  virtual int giveNumberOfDofs() const = 0;
  virtual const Variable &getU() const = 0;
  FEICellGeometry *cellGeometryWrapper = nullptr;

  public:
  ThirdMediumElement( int n, Domain *d ) :
          NLStructuralElement( n, d )
  {
    cellGeometryWrapper == NULL;
  }


  // Note: performance can be probably improved once it will be possible
  // to directly assemble multiple term contributions to the system matrix.
  // template metaprogramming?
  virtual void giveCharacteristicMatrix( FloatMatrix &answer, CharType type, TimeStep *tStep ) override
  {
    if ( type == TangentStiffnessMatrix ) {
      FloatMatrix standardAnswer, gradFanswer, gradJanswer;
      this->computeStiffnessMatrix( standardAnswer, TangentStiffness, tStep );

      //add contribution from gradient terms
      int udofs = this->giveNumberOfDofs();
      IntegrationRule *ir = this->giveDefaultIntegrationRulePtr();

      gradFanswer.resize( udofs, udofs );
      gradFanswer.zero();
      this->integrateTerm_dw( gradFanswer, ThirdMedium_GradGrad_SecondGradientTerm( getU(), getU() ), ir, tStep );

      gradJanswer.resize( udofs, udofs );
      gradJanswer.zero();
      this->integrateTerm_dw( gradJanswer, ThirdMedium_GradGrad_JacobianGradientTerm( getU(), getU() ), ir, tStep );

      answer.clear();
      answer.add( standardAnswer );
      answer.add( gradFanswer );
      answer.add( gradJanswer );
    } else {
      NLStructuralElement::giveCharacteristicMatrix( answer, type, tStep );
    }
  }

  virtual void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp,TimeStep *tStep ) override 
  {
    OOFEM_ERROR( "Small strain unsupported" );
  }

  virtual void computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override
  {
    MaterialMode mmode = gp->giveMaterialMode();
    if ( mmode == _PlaneStrain ) {
      answer = this->giveStructuralCrossSection()->giveStiffnessMatrix_dPdF_PlaneStrain( rMode, gp, tStep );
    } else if ( mmode == _3dMat ) {
      answer = this->giveStructuralCrossSection()->giveStiffnessMatrix_3d( rMode, gp, tStep );
    } else{
      OOFEM_ERROR("Unsupported material mode");
    }
  }

  virtual void giveCharacteristicVector( FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep ) override
  {
    if ( ( type == InternalForcesVector ) && ( mode == VM_Total ) ) {
      
      FloatArray standardAnswer, gradFanswer, gradJanswer;

      this->giveInternalForcesVector( standardAnswer, tStep );

      IntegrationRule *ir = this->giveDefaultIntegrationRulePtr();
      
      gradFanswer.resize( this->giveNumberOfDofs() );
      gradFanswer.zero();
      this->integrateTerm_c( gradFanswer, ThirdMedium_GradGrad_SecondGradientTerm( getU(), getU() ), ir, tStep );

      gradJanswer.resize( this->giveNumberOfDofs() );
      gradJanswer.zero();
      this->integrateTerm_c( gradJanswer, ThirdMedium_GradGrad_JacobianGradientTerm( getU(), getU() ), ir, tStep );

      answer.clear();
      answer.add( standardAnswer );
      answer.add( gradFanswer );
      answer.add( gradJanswer );

    } else {
      NLStructuralElement::giveCharacteristicVector( answer, type, mode, tStep );
    }
  }

  virtual void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override
  {
    OOFEM_ERROR("Small strain unsupported");
  }

  virtual void computeBmatrixAt( GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS ) override
  {
    OOFEM_ERROR( "Small strain unsupported" );
  }

  void integrateTerm_dw( FloatMatrix &answer, const StructuralTerm &term, IntegrationRule *iRule, TimeStep *tstep )
  {
    // need for integration domain and rule.
    // who should determine integration domain? Element or term? Term is just integrand, not integral
    // so integral type (surface, volume, etc) defined by element ---
    FloatMatrix dw;
    for ( GaussPoint *igp : *iRule ) {
      term.evaluate_lin( dw, *this, igp, tstep );
      dw.times( this->computeVolumeAround( igp ) );
      answer.add( dw );
    }
  }

  void integrateTerm_c( FloatArray &answer, const StructuralTerm &term, IntegrationRule *iRule, TimeStep *tstep )
  {
    // need for integration domain and rule.
    // who should determine integration domain? Element or term? Term is just integrand, not integral
    // so integral type (surface, volume, etc) defined by element ---
    FloatArray dw;
    for ( GaussPoint *igp : *iRule ) {
      term.evaluate( dw, *this, igp, tstep );
      dw.times( this->computeVolumeAround( igp ) );
      answer.add( dw );
    }
  }

  virtual void computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer ) override
  {
    FloatMatrix dNdx;
    MaterialMode mmode = gp->giveMaterialMode();
    FEInterpolation* interpol = this->giveInterpolation();
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    interpol->evaldNdx( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper( this ) );

    if ( mmode == _3dMat ) {
      answer.resize( 9, dNdx.giveNumberOfRows() * 3 );
      answer.zero();
      for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at( 1, 3 * i - 2 ) = dNdx.at( i, 1 ); // du/dx
        answer.at( 2, 3 * i - 1 ) = dNdx.at( i, 2 ); // dv/dy
        answer.at( 3, 3 * i - 0 ) = dNdx.at( i, 3 ); // dw/dz
        answer.at( 4, 3 * i - 1 ) = dNdx.at( i, 3 ); // dv/dz
        answer.at( 7, 3 * i - 0 ) = dNdx.at( i, 2 ); // dw/dy
        answer.at( 5, 3 * i - 2 ) = dNdx.at( i, 3 ); // du/dz
        answer.at( 8, 3 * i - 0 ) = dNdx.at( i, 1 ); // dw/dx
        answer.at( 6, 3 * i - 2 ) = dNdx.at( i, 2 ); // du/dy
        answer.at( 9, 3 * i - 1 ) = dNdx.at( i, 1 ); // dv/dx
      }
    } else if ( mmode == _PlaneStrain ) {
      answer.resize( 5, dNdx.giveNumberOfRows() * 2 );
      answer.zero();
      for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at( 1, 2 * i - 1 ) = dNdx.at( i, 1 ); // du/dx -1
        answer.at( 2, 2 * i - 0 ) = dNdx.at( i, 2 ); // dv/dy -2
        answer.at( 4, 2 * i - 1 ) = dNdx.at( i, 2 ); // du/dy -6
        answer.at( 5, 2 * i - 0 ) = dNdx.at( i, 1 ); // dv/dx -9
      }
    } else {
      OOFEM_ERROR( "Unsupported material mode" );
    }
  }

  double computeVolumeAround( GaussPoint *gp ) override
  {
    MaterialMode mmode = gp->giveMaterialMode();

    if ( mmode == _PlaneStrain ) {
      double weight = gp->giveWeight();
      const FloatArray &lCoords = gp->giveNaturalCoordinates(); // local/natural coords of the gp (parent domain)
      double detJ = fabs( this->giveInterpolation()->giveTransformationJacobian( lCoords, *this->giveCellGeometryWrapper() ) );
      double thickness = this->giveCrossSection()->give( CS_Thickness, gp ); // the cross section keeps track of the thickness

      return detJ * thickness * weight; // dV
    } else if ( mmode == _3dMat ) {
      double determinant, weight, volume;
      determinant = fabs( this->giveInterpolation()->giveTransformationJacobian( gp->giveNaturalCoordinates(),
          FEIElementGeometryWrapper( this ) ) );

      weight = gp->giveWeight();
      volume = determinant * weight;

      return volume;
    }
  }

  FEICellGeometry * giveCellGeometryWrapper()
  {
    if ( !cellGeometryWrapper ) {
      cellGeometryWrapper = new FEIElementGeometryWrapper( this );
    }

    return cellGeometryWrapper;
  }

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief 2D 3M elastic element with quadratic interpolation of displacements
 *
 */
class ThirdMediumQuad_q : public ThirdMediumElement
{
  protected:
  static FEInterpolation &interpol;
  const static Variable &u;

  public:
  ThirdMediumQuad_q( int n, Domain *d ) :
          ThirdMediumElement( n, d )
  {
    numberOfDofMans = 8;
    numberOfGaussPoints = 8;
    this->computeGaussPoints();
  }

  //int getNumberOfSurfaceDOFs() const override { return 0; }
  //int getNumberOfEdgeDOFs() const override { return 0; }
  //void getSurfaceLocalCodeNumbers( oofem::IntArray &, oofem::Variable::VariableQuantity ) const override { ; }

  //void getDofManLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q, int num ) const override
  //{
  //  /* dof ordering: u1 v1 u2 v2 u3 v3 u4 v4 u5 v5 u6 v6 */
  //  if ( q == Variable::VariableQuantity::Displacement ) {
  //    // answer={1,2,3 4,5,6,7,8,9,10,11,12,13,14,15,16};
  //    int o = ( num - 1 ) * 2 + 1;
  //    answer = { o, o + 1 };
  //  }
  //}
  /*void getInternalDofManLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q, int num ) const override
  {
    answer = {};
  }*/

  virtual FEInterpolation* giveInterpolation() const override { return &interpol; }

  void giveDofManDofIDMask( int inode, IntArray &answer ) const override
  {

    answer = { 1, 2 };
  }
  int giveNumberOfDofs() override { return 16; }
  const char *giveInputRecordName() const override { return "ThirdMediumQuad_q"; }
  const char *giveClassName() const override { return "ThirdMediumQuad_q"; }

  //const FEInterpolation &getGeometryInterpolation() const override { return this->interpol; }

  Element_Geometry_Type giveGeometryType() const override
  {
    return EGT_quad_2;
  }
  //void getEdgeLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q ) const override {}


  private:
  virtual int giveNumberOfDofs() const override { return 16; }
  virtual const Variable &getU() const override { return u; }
  void computeGaussPoints() override
  {
    if ( integrationRulesArray.size() == 0 ) {
      integrationRulesArray.resize( 1 );
      integrationRulesArray[0] = std::make_unique<LobattoIntegrationRule>( 1, this );
      integrationRulesArray[0]->SetUpPointsOnSquare( numberOfGaussPoints, _PlaneStrain );
    }
  }
};

FEInterpolation &ThirdMediumQuad_q ::interpol = FEI2dQuadQuad( 1, 2 );
const Variable &ThirdMediumQuad_q::u = Variable( ThirdMediumQuad_q::interpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 2, NULL, { 1, 2 } );

#define _IFT_ThirdMediumQuad_q_Name "thirdmediumquad_q"
REGISTER_Element( ThirdMediumQuad_q )

} //namespace oofem

