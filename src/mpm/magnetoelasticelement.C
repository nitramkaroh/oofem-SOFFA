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


#include "mpm.h"
#include "termlibrary4.h"
#include "termlibrary.h"
#include "element.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "intarray.h"
#include "classfactory.h"
#include "gaussintegrationrule.h"
#include "crosssection.h"

#include "fei2dquadquad.h"
#include "fei2dquadlin.h"
#include "mathfem.h"
#include "floatarrayf.h"


namespace oofem {

#define _IFT_MagnetoElasticElement_F1 "f01"
#define _IFT_MagnetoElasticElement_F2 "f02"
#define _IFT_MagnetoElasticElement_F3 "f03"
#define _IFT_MagnetoElasticElement_F4 "f04"

/**
 * @brief Base class for fully coupled, nonlinear thermo mechanical elements
 *
 */
class MagnetoElasticElement : public MPElement
{

  private:
  virtual int giveNumberOfUDofs() const = 0;
  virtual int giveNumberOfPhiDofs() const = 0;
  virtual const Variable &getU() const = 0;
  virtual const Variable &getPhi() const = 0;

  FloatArrayF<9> F0_1, F0_2, F0_3, F0_4; // prestrain for each GP

  public:
  MagnetoElasticElement( int n, Domain *d ) :
      MPElement( n, d ) {}

  void initializeFrom( InputRecord &ir ) override
  {

    MPElement::initializeFrom( ir );

    FloatArray F1_temp( 9 ), F2_temp( 9 ), F3_temp( 9 ), F4_temp( 9 );

    IR_GIVE_OPTIONAL_FIELD( ir, F1_temp, _IFT_MagnetoElasticElement_F1 );
    IR_GIVE_OPTIONAL_FIELD( ir, F2_temp, _IFT_MagnetoElasticElement_F2 );
    IR_GIVE_OPTIONAL_FIELD( ir, F3_temp, _IFT_MagnetoElasticElement_F3 );
    IR_GIVE_OPTIONAL_FIELD( ir, F4_temp, _IFT_MagnetoElasticElement_F4 );

    if ( F1_temp.giveSize() == 4 ) {
      F0_1 = assemble<9>( FloatArrayF<4>( F1_temp ), { 0, 1, 5, 8 } );
      F0_1.at( 3 ) = 1.0; // 1-based, unforunately
    } else {
      F0_1 = FloatArrayF<9>( F1_temp );
    }
    if ( F2_temp.giveSize() == 4 ) {
      F0_2 = assemble<9>( FloatArrayF<4>( F2_temp ), { 0, 1, 5, 8 } );
      F0_2.at( 3 ) = 1.0; // 1-based, unforunately
    } else {
      F0_2 = FloatArrayF<9>( F2_temp );
    }
    if ( F3_temp.giveSize() == 4 ) {
      F0_3 = assemble<9>( FloatArrayF<4>( F3_temp ), { 0, 1, 5, 8 } );
      F0_3.at( 3 ) = 1.0; // 1-based, unforunately
    } else {
      F0_3 = FloatArrayF<9>( F3_temp );
    }
    if ( F4_temp.giveSize() == 4 ) {
      F0_4 = assemble<9>( FloatArrayF<4>( F4_temp ), { 0, 1, 5, 8 } );
      F0_4.at( 3 ) = 1.0; // 1-based, unforunately
    } else {
      F0_4 = FloatArrayF<9>( F4_temp );
    }
  }

  // Note: performance can be probably improved once it will be possible
  // to directly assemble multiple term contributions to the system matrix.
  // template metaprogramming?
  virtual void giveCharacteristicMatrix( FloatMatrix &answer, CharType type, TimeStep *tStep ) override
  {
    IntegrationRule *ir = this->giveDefaultIntegrationRulePtr();

    if ( type == MagnetoElasticity_GradGrad_dFluxdGrad ) {
      int udofs = this->giveNumberOfUDofs();
      int phidofs = this->giveNumberOfPhiDofs();
      answer.resize( udofs + phidofs, udofs + phidofs );
      answer.zero();
      this->integrateTerm_dw( answer, MagnetoElasticity_GradGrad_Term( getU(), getU(), getPhi() ), ir, tStep );
    } else if ( type == MagnetoElasticity_GradGrad_SecondGradient_dFluxdGrad ) {
      answer.resize(this->giveNumberOfUDofs(),this->giveNumberOfUDofs());
    } else {
      OOFEM_ERROR( "Unknown characteristic matrix type" );
    }
  }

  virtual void giveCharacteristicVector( FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep ) override
  {
    IntegrationRule *ir = this->giveDefaultIntegrationRulePtr();
    if ( type == MagnetoElasticity_GradGrad_Flux ) {
      answer.resize( this->giveNumberOfUDofs() + this->giveNumberOfPhiDofs() );
      answer.zero();
      this->integrateTerm_c( answer, MagnetoElasticity_GradGrad_Term( getU(), getU(), getPhi() ), ir, tStep );
    } else if ( type == MagnetoElasticity_GradGrad_SecondGradient_Flux ) {
      answer.resize( this->giveNumberOfUDofs() );
    } else if ( type == ExternalForcesVector ) {
      ;
    } else {
      OOFEM_ERROR( "Unknown characteristic vector type" );
    }
  }

  virtual int giveIPValue( FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep ) override
  {
    if ( type == IST_PrestrainDeformationGradient ) {
      switch ( gp->giveNumber() ) {
      case 1:
        answer = FloatArray( F0_1 );
        break;
      case 2:
        answer = FloatArray( F0_2 );
        break;
      case 3:
        answer = FloatArray( F0_3 );
        break;
      case 4:
        answer = FloatArray( F0_4 );
        break;
      default:
        OOFEM_ERROR( "Prestrain requested for wrong GP number %i", gp->giveNumber() );
      }
    } else {
      return MPElement::giveIPValue( answer, gp, type, tStep );
    }
    return 1;
  }
};

/**
 * @brief 2D Magneto elatic element with quadratic interpolation of displacements and linear interpolation of magnetic potential
 *
 */
class MagnetoElasticQuad_qq : public MagnetoElasticElement
{
  protected:
  const static FEInterpolation &phiInterpol;
  const static FEInterpolation &uInterpol;
  const static Variable &phi;
  const static Variable &u;

  public:
  MagnetoElasticQuad_qq( int n, Domain *d ) :
      MagnetoElasticElement( n, d )
  {
    numberOfDofMans = 8;
    numberOfGaussPoints = 4;
    this->computeGaussPoints();
  }

  int getNumberOfSurfaceDOFs() const override { return 0; }
  int getNumberOfEdgeDOFs() const override { return 0; }
  void getSurfaceLocalCodeNumbers( oofem::IntArray &, oofem::Variable::VariableQuantity ) const override { ; }

  void getDofManLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q, int num ) const override
  {
    /* dof ordering: u1 v1 phi1  u2 v2 phi2  u3 v3 phi3  u4 v4 phi4   u5 v5 u6 v6 */
    if ( q == Variable::VariableQuantity::Displacement ) {
      // answer={1,2, 4,5,  7,8, 10,11,  13,14  15,16, 17,18  19,20};
      int o = ( num - 1 ) * 3 + 1;
      answer = { o, o + 1 };
    } else if ( q == Variable::VariableQuantity::MagneticPotential ) {
      answer = { num * 3 };
    }
  }
  void getInternalDofManLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q, int num ) const override
  {
    answer = {};
  }

  void giveDofManDofIDMask( int inode, IntArray &answer ) const override
  {

    answer = { 1, 2, 55 };
  }
  int giveNumberOfDofs() override { return 24; }
  const char *giveInputRecordName() const override { return "MagnetoElasticQuad_qq"; }
  const char *giveClassName() const override { return "MagnetoElasticQuad_qq"; }

  void giveCharacteristicMatrix( FloatMatrix &answer, CharType type, TimeStep *tStep ) override
  {
    IntegrationRule *ir = this->giveDefaultIntegrationRulePtr();

    if ( type == MagnetoElasticity_GradGrad_SecondGradient_dFluxdGrad ) {
      int udofs = this->giveNumberOfUDofs();
      answer.resize( udofs, udofs );
      answer.zero();
      this->integrateTerm_dw( answer, MagnetoElasticity_GradGrad_SecondGradientTerm( getU(), getU(), getPhi() ), ir, tStep );
    } else {
      MagnetoElasticElement::giveCharacteristicMatrix(answer, type, tStep);
    }
  }

  void giveCharacteristicVector( FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep ) override
  {
    IntegrationRule *ir = this->giveDefaultIntegrationRulePtr();
    if ( type == MagnetoElasticity_GradGrad_SecondGradient_Flux ) {
      answer.resize( this->giveNumberOfUDofs() );
      answer.zero();
      this->integrateTerm_c( answer, MagnetoElasticity_GradGrad_SecondGradientTerm( getU(), getU(), getPhi() ), ir, tStep );
    } else {
      MagnetoElasticElement::giveCharacteristicVector(answer,type,mode,tStep);
    }
  }


  const FEInterpolation &getGeometryInterpolation() const override { return this->uInterpol; }

  Element_Geometry_Type giveGeometryType() const override
  {
    return EGT_quad_2;
  }
  void getEdgeLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q ) const override {}


  private:
  virtual int giveNumberOfUDofs() const override { return 16; }
  virtual int giveNumberOfPhiDofs() const override { return 8; }
  virtual const Variable &getU() const override { return u; }
  virtual const Variable &getPhi() const override { return phi; }
  void computeGaussPoints() override
  {
    if ( integrationRulesArray.size() == 0 ) {
      integrationRulesArray.resize( 1 );
      integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this );
      integrationRulesArray[0]->SetUpPointsOnSquare( numberOfGaussPoints, _PlaneStrain );
    }
  }
};

const FEInterpolation &MagnetoElasticQuad_qq ::uInterpol = FEI2dQuadQuad( 1, 2 );
const FEInterpolation &MagnetoElasticQuad_qq ::phiInterpol = FEI2dQuadQuad( 1, 2 );
const Variable &MagnetoElasticQuad_qq::phi = Variable( MagnetoElasticQuad_qq::phiInterpol, Variable::VariableQuantity::MagneticPotential, Variable::VariableType::scalar, 1, NULL, { 55 } );
const Variable &MagnetoElasticQuad_qq::u = Variable( MagnetoElasticQuad_qq::uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 2, NULL, { 1, 2 } );

#define _IFT_MagnetoElasticQuad_qq_Name "magnetoelasticquad_qq"
REGISTER_Element( MagnetoElasticQuad_qq )


    /**
     * @brief 2D Magneto elatic element with quadratic interpolation of displacements and linear interpolation of magnetic potential
     *
     */
    class MagnetoElasticQuad_ql : public MagnetoElasticElement
{
  protected:
  const static FEInterpolation &phiInterpol;
  const static FEInterpolation &uInterpol;
  const static Variable &phi;
  const static Variable &u;

  public:
  MagnetoElasticQuad_ql( int n, Domain *d ) :
      MagnetoElasticElement( n, d )
  {
    numberOfDofMans = 8;
    numberOfGaussPoints = 8;
    this->computeGaussPoints();
  }

  int getNumberOfSurfaceDOFs() const override { return 0; }
  int getNumberOfEdgeDOFs() const override { return 0; }
  void getSurfaceLocalCodeNumbers( oofem::IntArray &, oofem::Variable::VariableQuantity ) const override { ; }

  void getDofManLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q, int num ) const override
  {
    /* dof ordering: u1 v1 phi1  u2 v2 phi2  u3 v3 phi3  u4 v4 phi4   u5 v5 u6 v6 */
    if ( q == Variable::VariableQuantity::Displacement ) {
      // answer={1,2, 4,5,  7,8, 10,11,  13,14  15,16, 17,18  19,20};
      int o = ( num - 1 ) * 3 + 1 - ( num > 4 ) * ( num - 5 );
      answer = { o, o + 1 };
    } else if ( q == Variable::VariableQuantity::MagneticPotential ) {
      if ( num <= 4 ) {
        answer = { num * 3 };
      } else {
        answer = {};
      }
    }
  }
  void getInternalDofManLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q, int num ) const override
  {
    answer = {};
  }

  void giveDofManDofIDMask( int inode, IntArray &answer ) const override
  {
    if ( inode <= 4 ) {
      answer = { 1, 2, 55 };
    } else {
      answer = { 1, 2 };
    }
  }
  int giveNumberOfDofs() override { return 20; }
  const char *giveInputRecordName() const override { return "MagnetoElasticQuad_ql"; }
  const char *giveClassName() const override { return "MagnetoElasticQuad_ql"; }


  const FEInterpolation &getGeometryInterpolation() const override { return this->uInterpol; }

  Element_Geometry_Type giveGeometryType() const override
  {
    return EGT_quad_2;
  }
  void getEdgeLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q ) const override {}


  private:
  virtual int giveNumberOfUDofs() const override { return 16; }
  virtual int giveNumberOfPhiDofs() const override { return 4; }
  virtual const Variable &getU() const override { return u; }
  virtual const Variable &getPhi() const override { return phi; }
  void computeGaussPoints() override
  {
    if ( integrationRulesArray.size() == 0 ) {
      integrationRulesArray.resize( 1 );
      integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this );
      integrationRulesArray[0]->SetUpPointsOnSquare( numberOfGaussPoints, _PlaneStrain );
    }
  }
};

const FEInterpolation &MagnetoElasticQuad_ql ::uInterpol = FEI2dQuadQuad( 1, 2 );
const FEInterpolation &MagnetoElasticQuad_ql ::phiInterpol = FEI2dQuadLin( 1, 2 );
const Variable &MagnetoElasticQuad_ql::phi = Variable( MagnetoElasticQuad_ql::phiInterpol, Variable::VariableQuantity::MagneticPotential, Variable::VariableType::scalar, 1, NULL, { 55 } );
const Variable &MagnetoElasticQuad_ql::u = Variable( MagnetoElasticQuad_ql::uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 2, NULL, { 1, 2 } );

#define _IFT_MagnetoElasticQuad_ql_Name "magnetoelasticquad_ql"
REGISTER_Element( MagnetoElasticQuad_ql )


/**
 * @brief 2D Magneto elatic element with linear interpolation of displacements and magnetic potential
 *
 */

class MagnetoElasticQuad_ll : public MagnetoElasticElement
{

  protected:
  const static FEInterpolation &phiInterpol;
  const static FEInterpolation &uInterpol;
  const static Variable &phi;
  const static Variable &u;

  public:
  MagnetoElasticQuad_ll( int n, Domain *d ) :
      MagnetoElasticElement( n, d )
  {
    numberOfDofMans = 4;
    numberOfGaussPoints = 4;
    this->computeGaussPoints();
  }

  int getNumberOfSurfaceDOFs() const override { return 0; }
  int getNumberOfEdgeDOFs() const override { return 0; }
  void getSurfaceLocalCodeNumbers( oofem::IntArray &, oofem::Variable::VariableQuantity ) const override { ; }

  void getDofManLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q, int num ) const override
  {
    /* dof ordering: u1 v1 phi1  u2 v2 phi2  u3 v3 phi3  u4 v4 phi4   u5 v5 u6 v6 */
    if ( q == Variable::VariableQuantity::Displacement ) {
      // answer={1,2, 4,5,  7,8, 10,11,  13,14  15,16, 17,18  19,20};
      int o = ( num - 1 ) * 3 + 1;
      answer = { o, o + 1 };
    } else if ( q == Variable::VariableQuantity::MagneticPotential ) {
      answer = { num * 3 };
    }
  }
  void getInternalDofManLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q, int num ) const override
  {
    answer = {};
  }

  void giveDofManDofIDMask( int inode, IntArray &answer ) const override
  {
    answer = { 1, 2, 55 };
  }
  int giveNumberOfDofs() override { return 12; }
  const char *giveInputRecordName() const override { return "MagnetoElasticQuad_ll"; }
  const char *giveClassName() const override { return "MagnetoElasticQuad_ll"; }


  const FEInterpolation &getGeometryInterpolation() const override { return this->uInterpol; }

  Element_Geometry_Type giveGeometryType() const override
  {
    return EGT_quad_1;
  }
  void getEdgeLocalCodeNumbers( IntArray &answer, const Variable::VariableQuantity q ) const override {}

  private:
  virtual int giveNumberOfUDofs() const override { return 8; }
  virtual int giveNumberOfPhiDofs() const override { return 4; }
  virtual const Variable &getU() const override { return u; }
  virtual const Variable &getPhi() const override { return phi; }
  void computeGaussPoints() override
  {
    if ( integrationRulesArray.size() == 0 ) {
      integrationRulesArray.resize( 1 );
      integrationRulesArray[0] = std::make_unique<GaussIntegrationRule>( 1, this );
      integrationRulesArray[0]->SetUpPointsOnSquare( numberOfGaussPoints, _PlaneStrain );
    }
  }
};

const FEInterpolation &MagnetoElasticQuad_ll ::uInterpol = FEI2dQuadLin( 1, 2 );
const FEInterpolation &MagnetoElasticQuad_ll ::phiInterpol = FEI2dQuadLin( 1, 2 );
const Variable &MagnetoElasticQuad_ll::phi = Variable( MagnetoElasticQuad_ll::phiInterpol, Variable::VariableQuantity::MagneticPotential, Variable::VariableType::scalar, 1, NULL, { 55 } );
const Variable &MagnetoElasticQuad_ll::u = Variable( MagnetoElasticQuad_ll::uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 2, NULL, { 1, 2 } );

#define _IFT_MagnetoElasticQuad_ll_Name "magnetoelasticquad_ll"
REGISTER_Element( MagnetoElasticQuad_ll )


} // end namespace oofem
