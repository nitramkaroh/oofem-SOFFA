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
#include "fei2dquadlin.h"
#include "fei3dhexalin.h"
#include "fei3dwedgelin.h"
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
      answer.clear();
      FloatMatrix standardAnswer, secondGradAnswer, firstSecondGradAnswer, FbarAnswer;
      int udofs = this->giveNumberOfDofs();
      IntegrationRule *ir = this->giveDefaultIntegrationRulePtr();

      //standard answer
      this->computeStiffnessMatrix( standardAnswer, TangentStiffness, tStep );
      answer.add( standardAnswer );

      //terms using Fbar
      FbarAnswer.resize( udofs, udofs );
      FbarAnswer.zero();
      this->integrateTerm_dw( FbarAnswer, ThirdMedium_Grad_FbarTerm( getU(), getU() ), ir, tStep );
      answer.add( FbarAnswer );

      //gradient terms
      if ( this->giveInterpolation()->giveInterpolationOrder() > 1 ) {
        secondGradAnswer.resize( udofs, udofs );
        secondGradAnswer.zero();
        this->integrateTerm_dw( secondGradAnswer, ThirdMedium_GradGrad_SecondGradientTerm( getU(), getU() ), ir, tStep );

        firstSecondGradAnswer.resize( udofs, udofs );
        firstSecondGradAnswer.zero();
        this->integrateTerm_dw( firstSecondGradAnswer, ThirdMedium_GradGrad_FirstSecondGradientTerm( getU(), getU() ), ir, tStep );

        answer.add( secondGradAnswer );
        answer.add( firstSecondGradAnswer );
      }

    } else {
      NLStructuralElement::giveCharacteristicMatrix( answer, type, tStep );
    }
  }

  virtual void computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep ) override
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
    } else {
      OOFEM_ERROR( "Unsupported material mode" );
    }
  }

  virtual void giveCharacteristicVector( FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep ) override
  {
    if ( ( type == InternalForcesVector ) && ( mode == VM_Total ) ) {
      answer.clear();
      FloatArray standardAnswer, firstGradAnswer, firstSecondGradAnswer, FbarAnswer;
      IntegrationRule *ir = this->giveDefaultIntegrationRulePtr();

      //standard answer
      this->giveInternalForcesVector( standardAnswer, tStep );
      answer.add( standardAnswer );

      //Fbar answer
      FbarAnswer.resize( this->giveNumberOfDofs() );
      FbarAnswer.zero();
      this->integrateTerm_c( FbarAnswer, ThirdMedium_Grad_FbarTerm( getU(), getU() ), ir, tStep );
      answer.add( FbarAnswer );

      // gradient answers
      if ( this->giveInterpolation()->giveInterpolationOrder() > 1 ) {
        firstGradAnswer.resize( this->giveNumberOfDofs() );
        firstGradAnswer.zero();
        this->integrateTerm_c( firstGradAnswer, ThirdMedium_GradGrad_SecondGradientTerm( getU(), getU() ), ir, tStep );
        answer.add( firstGradAnswer );

        firstSecondGradAnswer.resize( this->giveNumberOfDofs() );
        firstSecondGradAnswer.zero();
        this->integrateTerm_c( firstSecondGradAnswer, ThirdMedium_GradGrad_FirstSecondGradientTerm( getU(), getU() ), ir, tStep );
        answer.add( firstSecondGradAnswer );
      }
      
    } else {
      NLStructuralElement::giveCharacteristicVector( answer, type, mode, tStep );
    }
  }

  virtual void computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep ) override
  {
    OOFEM_ERROR( "Small strain unsupported" );
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
    FEInterpolation *interpol = this->giveInterpolation();
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

  FEICellGeometry *giveCellGeometryWrapper()
  {
    if ( !cellGeometryWrapper ) {
      cellGeometryWrapper = new FEIElementGeometryWrapper( this );
    }

    return cellGeometryWrapper;
  }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define _IFT_ThirdMediumQuad_q_lobatto "lobatto"


/**
 * @brief 2D 3M elastic element with quadratic interpolation of displacements
 *
 */
class ThirdMediumQuad_q : public ThirdMediumElement
{
  protected:
  static FEI2dQuadQuad interpol;
  const static Variable &u;
  bool lobatto = false;

  public:
  ThirdMediumQuad_q( int n, Domain *d ) :
      ThirdMediumElement( n, d )
  {
    numberOfDofMans = 8;
    numberOfGaussPoints = 8;
    //this->computeGaussPoints();
  }

  virtual void initializeFrom( InputRecord &ir ) override
  {
    ThirdMediumElement::initializeFrom( ir );
    lobatto = ir.hasField( _IFT_ThirdMediumQuad_q_lobatto );
  }

  virtual FEInterpolation *giveInterpolation() const override { return &interpol; }

  void giveDofManDofIDMask( int inode, IntArray &answer ) const override
  {

    answer = { 1, 2 };
  }
  int giveNumberOfDofs() override { return 16; }
  const char *giveInputRecordName() const override { return "ThirdMediumQuad_q"; }
  const char *giveClassName() const override { return "ThirdMediumQuad_q"; }

  // const FEInterpolation &getGeometryInterpolation() const override { return this->interpol; }

  Element_Geometry_Type giveGeometryType() const override
  {
    return EGT_quad_2;
  }
  // void getEdgeLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q ) const override {}


  private:
  virtual int giveNumberOfDofs() const override { return 16; }
  virtual const Variable &getU() const override { return u; }
  void computeGaussPoints() override
  {
    if ( integrationRulesArray.size() == 0 ) {
      integrationRulesArray.resize( 1 );
      if (lobatto) {
        integrationRulesArray[0] = std::make_unique<LobattoIntegrationRule>( 1, this );
      } else{
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this );
      }
      integrationRulesArray[0]->SetUpPointsOnSquare( numberOfGaussPoints, _PlaneStrain );
    }
  }
};

FEI2dQuadQuad ThirdMediumQuad_q ::interpol( 1, 2 );
const Variable &ThirdMediumQuad_q::u = Variable( ThirdMediumQuad_q::interpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 2, NULL, { 1, 2 } );

#define _IFT_ThirdMediumQuad_q_Name "thirdmediumquad_q"
REGISTER_Element( ThirdMediumQuad_q )

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define _IFT_ThirdMediumQuad_l_lobatto "lobatto"

/**
 * @brief 2D 3M elastic element with linear interpolation of displacements
 *
 */
class ThirdMediumQuad_l : public ThirdMediumElement
{
  protected:
  static FEI2dQuadLin interpol;
  const static Variable &u;
  bool lobatto = false;

  public:
  ThirdMediumQuad_l( int n, Domain *d ) :
      ThirdMediumElement( n, d )
  {
    numberOfDofMans = 4;
    numberOfGaussPoints = 4;
    //this->computeGaussPoints();
  }

  virtual void initializeFrom( InputRecord &ir ) override
  {
    ThirdMediumElement::initializeFrom( ir );
    lobatto = ir.hasField( _IFT_ThirdMediumQuad_l_lobatto );
  }

  virtual FEInterpolation *giveInterpolation() const override { return &interpol; }

  void giveDofManDofIDMask( int inode, IntArray &answer ) const override
  {

    answer = { 1, 2 };
  }
  int giveNumberOfDofs() override { return 8; }
  const char *giveInputRecordName() const override { return "ThirdMediumQuad_l"; }
  const char *giveClassName() const override { return "ThirdMediumQuad_l"; }

  // const FEInterpolation &getGeometryInterpolation() const override { return this->interpol; }

  Element_Geometry_Type giveGeometryType() const override
  {
    return EGT_quad_1;
  }
  // void getEdgeLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q ) const override {}


  private:
  virtual int giveNumberOfDofs() const override { return 8; }
  virtual const Variable &getU() const override { return u; }
  void computeGaussPoints() override
  {
    if ( integrationRulesArray.size() == 0 ) {
      integrationRulesArray.resize( 1 );
      if ( lobatto ) {
        integrationRulesArray[0] = std::make_unique<LobattoIntegrationRule>( 1, this );
      } else {
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this );
      }
      integrationRulesArray[0]->SetUpPointsOnSquare( numberOfGaussPoints, _PlaneStrain );
    }
  }
};

FEI2dQuadLin ThirdMediumQuad_l ::interpol( 1, 2 );
const Variable &ThirdMediumQuad_l::u = Variable( ThirdMediumQuad_l::interpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 2, NULL, { 1, 2 } );

#define _IFT_ThirdMediumQuad_l_Name "thirdmediumquad_l"
REGISTER_Element( ThirdMediumQuad_l )

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define _IFT_ThirdMediumBrick_l_lobatto "lobatto"

/**
 * @brief 3D 3M elastic brick element with linear interpolation of displacements
 *
 */
class ThirdMediumBrick_l : public ThirdMediumElement
{
  protected:
  static FEI3dHexaLin interpol;
  const static Variable &u;
  bool lobatto = false;

  public:
  ThirdMediumBrick_l( int n, Domain *d ) :
      ThirdMediumElement( n, d )
  {
    numberOfDofMans = 8;
    numberOfGaussPoints = 8;
  }

  virtual void initializeFrom( InputRecord &ir ) override
  {
    ThirdMediumElement::initializeFrom( ir );
    lobatto = ir.hasField( _IFT_ThirdMediumBrick_l_lobatto );
  }

  virtual FEInterpolation *giveInterpolation() const override { return &interpol; }

  void giveDofManDofIDMask( int inode, IntArray &answer ) const override
  {

    answer = { 1, 2, 3 };
  }
  int giveNumberOfDofs() override { return 8; }
  const char *giveInputRecordName() const override { return "ThirdMediumBrick_l"; }
  const char *giveClassName() const override { return "ThirdMediumBrick_l"; }

  // const FEInterpolation &getGeometryInterpolation() const override { return this->interpol; }

  Element_Geometry_Type giveGeometryType() const override
  {
    return EGT_hexa_1;
  }
  // void getEdgeLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q ) const override {}


  private:
  virtual int giveNumberOfDofs() const override { return 3*8; }
  virtual const Variable &getU() const override { return u; }
  void computeGaussPoints() override
  {
    if ( integrationRulesArray.size() == 0 ) {
      integrationRulesArray.resize( 1 );
      if ( lobatto ) {
        integrationRulesArray[0] = std::make_unique<LobattoIntegrationRule>( 1, this );
      } else {
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this );
      }
      integrationRulesArray[0]->SetUpPointsOnCube( numberOfGaussPoints, _3dMat );
    }
  }
};

FEI3dHexaLin ThirdMediumBrick_l ::interpol;
const Variable &ThirdMediumBrick_l::u = Variable( ThirdMediumBrick_l::interpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 3, NULL, { 1, 2, 3 } );

#define _IFT_ThirdMediumBrick_l_Name "thirdmediumbrick_l"
REGISTER_Element( ThirdMediumBrick_l )

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define _IFT_ThirdMediumWedge_l_lobatto "lobatto"

    /**
     * @brief 3D 3M elastic wedge element with linear interpolation of displacements
     *
     */
    class ThirdMediumWedge_l : public ThirdMediumElement
{
  protected:
  static FEI3dWedgeLin interpol;
  const static Variable &u;
  bool lobatto = false;

  public:
  ThirdMediumWedge_l( int n, Domain *d ) :
      ThirdMediumElement( n, d )
  {
    numberOfDofMans = 6;
    numberOfGaussPoints = 6;
  }

  virtual void initializeFrom( InputRecord &ir ) override
  {
    ThirdMediumElement::initializeFrom( ir );
    lobatto = ir.hasField( _IFT_ThirdMediumWedge_l_lobatto );
  }

  virtual FEInterpolation *giveInterpolation() const override { return &interpol; }

  void giveDofManDofIDMask( int inode, IntArray &answer ) const override
  {

    answer = { 1, 2, 3 };
  }
  int giveNumberOfDofs() override { return 8; }
  const char *giveInputRecordName() const override { return "ThirdMediumWedge_l"; }
  const char *giveClassName() const override { return "ThirdMediumWedge_l"; }

  // const FEInterpolation &getGeometryInterpolation() const override { return this->interpol; }

  Element_Geometry_Type giveGeometryType() const override
  {
    return EGT_wedge_1;
  }
  // void getEdgeLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q ) const override {}


  private:
  virtual int giveNumberOfDofs() const override { return 3 * numberOfDofMans; }
  virtual const Variable &getU() const override { return u; }
  void computeGaussPoints() override
  {
    if ( integrationRulesArray.size() == 0 ) {
      integrationRulesArray.resize( 1 );
      if ( lobatto ) {
        integrationRulesArray[0] = std::make_unique<LobattoIntegrationRule>( 1, this );
      } else {
        integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this );
      }
      integrationRulesArray[0]->SetUpPointsOnWedge( 3, 2, _3dMat );
    }
  }
};

FEI3dWedgeLin ThirdMediumWedge_l ::interpol;
const Variable &ThirdMediumWedge_l::u = Variable( ThirdMediumWedge_l::interpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 3, NULL, { 1, 2, 3 } );

#define _IFT_ThirdMediumWedge_l_Name "thirdmediumwedge_l"
REGISTER_Element( ThirdMediumWedge_l )


} // namespace oofem
