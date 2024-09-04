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

#include "fei2dquadquad.h"
#include "fei2dquadlin.h"
#include "mathfem.h"


namespace oofem {


/**
 * @brief Base class for fully coupled, nonlinear thermo mechanical elements
 * 
 */
class MagnetoElasticElement : public MPElement {
        
    private:
        virtual int  giveNumberOfUDofs() const = 0;
        virtual int  giveNumberOfPhiDofs() const = 0;
        virtual const Variable& getU() const = 0;
        virtual const Variable& getPhi() const = 0;

    public:
    MagnetoElasticElement(int n, Domain* d): 
        MPElement(n,d) { }

    // Note: performance can be probably improved once it will be possible 
    // to directly assemble multiple term contributions to the system matrix.
    // template metaprogramming?
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) override {
        IntegrationRule* ir = this->giveDefaultIntegrationRulePtr();

        if (type == MagnetoElasticity_GradGrad_dFluxdGrad) {
            int udofs = this->giveNumberOfUDofs();
	    int phidofs = this->giveNumberOfPhiDofs(); 
            answer.resize(udofs+phidofs,udofs+phidofs);
            answer.zero();
            this->integrateTerm_dw (answer, MagnetoElasticity_GradGrad_Term(getU(),getU(), getPhi()), ir, tStep) ;
        } else {
	        OOFEM_ERROR("Unknown characteristic matrix type");
	    }
    }

    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) override {
        IntegrationRule* ir = this->giveDefaultIntegrationRulePtr();
        if (type == MagnetoElasticity_GradGrad_Flux) {
            answer.resize(this->giveNumberOfUDofs() + this->giveNumberOfPhiDofs());
            answer.zero();
            this->integrateTerm_c (answer, MagnetoElasticity_GradGrad_Term(getU(),getU(),getPhi()), ir, tStep);
        } else if(type == ExternalForcesVector) {
	  ;
	} else {
	  OOFEM_ERROR("Unknown characteristic vector type");
	}
    }    
};



/**
 * @brief 2D Magneto elatic element with quadratic interpolation of displacements and linear interpolation of magnetic potential
 * 
 */
class MagnetoElasticQuad_ql : public MagnetoElasticElement {
    protected:
        const static FEInterpolation & phiInterpol;
        const static FEInterpolation & uInterpol;
        const static Variable& phi;
        const static Variable& u;
      
    public:
  MagnetoElasticQuad_ql(int n, Domain* d):  MagnetoElasticElement(n,d)
    {
        numberOfDofMans  = 8;
        numberOfGaussPoints = 8;
        this->computeGaussPoints();
    }

  int getNumberOfSurfaceDOFs() const override {return 0;}
  int getNumberOfEdgeDOFs() const override {return 0;}
  void getSurfaceLocalCodeNumbers(oofem::IntArray&, oofem::Variable::VariableQuantity) const override {;}

  void getDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        /* dof ordering: u1 v1 phi1  u2 v2 phi2  u3 v3 phi3  u4 v4 phi4   u5 v5 u6 v6 */
        if (q == Variable::VariableQuantity::Displacement) {
          //answer={1,2, 4,5,  7,8, 10,11,  13,14  15,16, 17,18  19,20};
	  int o = (num-1)*3+1 - (num>4)*(num-5);
          answer={o, o+1};
        } else if (q == Variable::VariableQuantity::MagneticPotential) {
	  if(num <= 4) {
	    answer={num*3};
	  } else {
	    answer = {};
	  }
        }
    }
    void getInternalDofManLocalCodeNumbers (IntArray& answer, const Variable::VariableQuantity q, int num ) const  override {
        answer={};
    }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override {
      if(inode <= 4) {
            answer = {1,2,55};
      } else {
	answer = {1,2};
      }
    }
    int giveNumberOfDofs() override { return 20; }
    const char *giveInputRecordName() const override {return "MagnetoElasticQuad_ql";}
    const char *giveClassName() const override { return "MagnetoElasticQuad_ql"; }

    
    const FEInterpolation& getGeometryInterpolation() const override {return this->uInterpol;}
  
    Element_Geometry_Type giveGeometryType() const override {
        return EGT_quad_2;
    }
    void getEdgeLocalCodeNumbers(IntArray& answer, const Variable::VariableQuantity q) const override {}


private:
        virtual int  giveNumberOfUDofs() const override {return 16;} 
        virtual int  giveNumberOfPhiDofs() const override {return 4;}
        virtual const Variable& getU() const override {return u;}
        virtual const Variable& getPhi() const override {return phi;}
        void computeGaussPoints() override {
            if ( integrationRulesArray.size() == 0 ) {
                integrationRulesArray.resize( 1 );
                integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this);
                integrationRulesArray [ 0 ]->SetUpPointsOnSquare(numberOfGaussPoints, _PlaneStrain);
            }
        }
};

  const FEInterpolation &  MagnetoElasticQuad_ql :: uInterpol = FEI2dQuadQuad(1,2);
  const FEInterpolation &  MagnetoElasticQuad_ql :: phiInterpol = FEI2dQuadLin(1,2);
const Variable&   MagnetoElasticQuad_ql::phi = Variable(  MagnetoElasticQuad_ql::phiInterpol, Variable::VariableQuantity::MagneticPotential, Variable::VariableType::scalar, 1, NULL, {55});
const Variable&   MagnetoElasticQuad_ql::u = Variable(MagnetoElasticQuad_ql::uInterpol, Variable::VariableQuantity::Displacement, Variable::VariableType::vector, 2, NULL, {1,2});

#define _IFT_MagnetoElasticQuad_ql_Name "magnetoelasticquad_ql"
REGISTER_Element(MagnetoElasticQuad_ql)




} // end namespace oofem
