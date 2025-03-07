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
 *
 *               Copyright (C) 1993 - 2015   Borek Patzak
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


#include "../sm/Elements/MixedPressure/basemixedpressureelement.h"
#include "../sm/Materials/MixedPressure/mixedpressurematerialextensioninterface.h"


#include "../sm/CrossSections/simplemixedupcrosssection.h"

#include "../sm/Materials/structuralms.h"

#include "material.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "domain.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "unknownnumberingscheme.h"


#include <cstdio>

namespace oofem {
BaseMixedPressureElement :: BaseMixedPressureElement()
{}



void
BaseMixedPressureElement :: giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_p, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u, const IntArray &dofIdArray_p)
{
    // Routine to extract the location array of an element for given dofid array.
    locationArray_u.clear();
    locationArray_p.clear();
    NLStructuralElement *el = this->giveElement();
    int k = 0;
    IntArray nodalArray;
    for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = el->giveDofManager(i);
        int itt = 1;
        for ( int j = 1; j <= dofIdArray_u.giveSize(); j++ ) {
            if ( dMan->hasDofID( ( DofIDItem ) dofIdArray_u.at(j) ) ) {
                //  Dof *d = dMan->giveDofWithID( dofIdArray_u.at( j ) );
                locationArray_u.followedBy(k + itt);
		itt++;
            }

        }
        for ( int j = 1; j <= dofIdArray_p.giveSize(); j++ ) {
            if ( dMan->hasDofID( ( DofIDItem ) dofIdArray_p.at(j) ) ) {
                //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
                locationArray_p.followedBy(k + itt);
		itt++;
            }

        }
        k += dMan->giveNumberOfDofs();
    }

    for ( int i = 1; i <= el->giveNumberOfInternalDofManagers(); i++ ) {
        DofManager *dMan = el->giveInternalDofManager(i);
        int itt = 1;
        for ( int j = 1; j <= dofIdArray_u.giveSize(); j++ ) {
            if ( dMan->hasDofID( ( DofIDItem ) dofIdArray_u.at(j) ) ) {
                //  Dof *d = dMan->giveDofWithID( dofIdArray_u.at( j ) );
                locationArray_u.followedBy(k + itt);
		itt++;
            }

        }
        for ( int j = 1; j <= dofIdArray_p.giveSize(); j++ ) {
            if ( dMan->hasDofID( ( DofIDItem ) dofIdArray_p.at(j) ) ) {
                //Dof *d = dMan->giveDofWithID( dofIdArray_m.at( j ) );
                locationArray_p.followedBy(k + itt);
		itt++;
            }

        }
        k += dMan->giveNumberOfDofs();
    }
    
}


void
BaseMixedPressureElement :: computeStressVector(FloatArray &stress, double &volumeChange , GaussPoint *gp, TimeStep *tStep)
{
    NLStructuralElement *elem = this->giveElement();
    double pressure;

    if ( elem->giveGeometryMode() == 0 ) {
        FloatArray strain;
        this->computeStrainVector(strain, gp, tStep);
        this->computePressure(pressure, gp, tStep);
        this->computeSmallStrainStressVector(stress, volumeChange, strain, pressure, gp, tStep);
    } else {
         FloatArray F;
	 elem->computeDeformationGradientVector(F, gp, tStep);
	 this->computePressure(pressure, gp, tStep);
	 this->computeLargeStrainStressVector(stress, volumeChange,F, pressure, gp, tStep);
    }
}





  
void
BaseMixedPressureElement :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray d_u;
    FloatMatrix B;
    NLStructuralElement *elem = this->giveElement();
    IntArray IdMask_u;
    this->giveDofManDofIDMask_u(IdMask_u);
    this->giveElement()->computeVectorOf(IdMask_u, VM_Total, tStep, d_u);
    elem->computeBmatrixAt(gp, B);
    answer.beProductOf(B, d_u);
}


void
BaseMixedPressureElement :: computePressure(double &answer, GaussPoint *gp, TimeStep *tStep)
{
    IntArray IdMask_p;
    FloatArray d_p, N_p;

    this->giveDofManDofIDMask_p(IdMask_p);
    this->giveElement()->computeVectorOf(IdMask_p, VM_Total, tStep, d_p);
    this->computePressureNMatrixAt(gp, N_p);
    answer = N_p.dotProduct(d_p);
}


void
BaseMixedPressureElement :: computeSmallStrainStressVector(FloatArray &answer, double &volChange, const FloatArray &strain, double p, GaussPoint *gp, TimeStep *tStep)
{
  
  SimpleMixedUPCrossSection *cs = static_cast<  SimpleMixedUPCrossSection * >( this->giveElement()->giveCrossSection());
  if ( !cs ) {
    OOFEM_ERROR("Wrong Corss-section!");
  }

  std::tie(answer, volChange) = cs->giveRealStresses(strain, p, gp, tStep);
}


void
BaseMixedPressureElement :: computeLargeStrainStressVector(FloatArray &answer, double &volChange, const FloatArray &F, double p, GaussPoint *gp, TimeStep *tStep)
{
  SimpleMixedUPCrossSection *cs = static_cast<  SimpleMixedUPCrossSection * >( this->giveElement()->giveCrossSection());
  if ( !cs ) {
    OOFEM_ERROR("Wrong Corss-section!");
  }

  std::tie(answer, volChange) = cs->giveFirstPKStresses(F, p, gp, tStep);
  
}




void
BaseMixedPressureElement :: computeSmallStrainMaterialStiffnessMatrix(FloatMatrix &Kuu, FloatArray &Kup, double &K, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  
  SimpleMixedUPCrossSection *cs = static_cast<  SimpleMixedUPCrossSection * >( this->giveElement()->giveCrossSection());
  if ( !cs ) {
    OOFEM_ERROR("Wrong Corss-section!");
  }

  cs->giveSmallStrainStiffnessMatrix(Kuu, Kup, K, rMode, gp, tStep);
}


void
BaseMixedPressureElement :: computeLargeStrainMaterialStiffnessMatrix(FloatMatrix &Kuu, FloatArray &Kup, double &K, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
  SimpleMixedUPCrossSection *cs = static_cast<  SimpleMixedUPCrossSection * >( this->giveElement()->giveCrossSection());
  if ( !cs ) {
    OOFEM_ERROR("Wrong Corss-section!");
  }
  double pressure;
  this->computePressure(pressure, gp, tStep);
  cs->giveLargeStrainStiffnessMatrix(Kuu, Kup, K, pressure, rMode, gp, tStep);
    
}


  

  

void
BaseMixedPressureElement :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    NLStructuralElement *elem = this->giveElement();
    SimpleMixedUPCrossSection *cs = static_cast<  SimpleMixedUPCrossSection * >( elem->giveCrossSection());
  if ( !cs ) {
    OOFEM_ERROR("Wrong Corss-section!");
  }
    double volumeChange;
    FloatArray BS, vStress, Np;
    FloatMatrix B, Nv, mNp;
    //
    answer.resize( this->giveNumberOfDofs() );
    answer.zero();
    //
    FloatArray answer_u( this->giveNumberOfDisplacementDofs() );
    answer_u.zero();
    FloatArray answer_p( this->giveNumberOfPressureDofs() );
    answer_p.zero();
    //
    for ( GaussPoint *gp : *elem->giveIntegrationRule(0) ) {
      this->computeStressVector(vStress, volumeChange, gp, tStep);
      // Compute nodal internal forces at nodes as f = B^T*Stress dV
      double dV  = elem->computeVolumeAround(gp);
      //
      if ( elem->giveGeometryMode() == 0 ) {
 	elem->computeBmatrixAt(gp, B);
      } else {
	elem->computeBHmatrixAt(gp, B);
      }
      BS.beTProductOf(B, vStress);
      answer_u.add(dV, BS);
      //
      this->computePressureNMatrixAt(gp, Np);
      answer_p.add(dV*volumeChange, Np);
      
    }
    answer.assemble(answer_u, locationArray_u);
    answer.assemble(answer_p, locationArray_p);  
}


void
BaseMixedPressureElement :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    FloatArray localForces( this->giveNumberOfDisplacementDofs() );
    answer.resize( this->giveNumberOfDofs() );
    this->computeLocForceLoadVector(localForces, tStep, mode);
    answer.assemble(localForces, locationArray_u);
}


/************************************************************************/
void
BaseMixedPressureElement :: computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting
// on element volume (surface).
// When reactions forces are computed, they are computed from element::GiveRealStressVector
// in this vector a real forces are stored (temperature part is subtracted).
// so we need further sobstract part corresponding to non-nodeal loading.
{
    FloatMatrix T;
    NLStructuralElement *elem = this->giveElement();
    //@todo check this
    //    elem->computeLocalForceLoadVector(answer, tStep, mode);

    // transform result from global cs to nodal cs. if necessary
    if ( answer.isNotEmpty() ) {
        if ( elem->computeGtoLRotationMatrix(T) ) {
            // first back to global cs from element local
            answer.rotatedWith(T, 't');
        }
    } else {
        answer.resize( this->giveNumberOfDisplacementDofs() );
        answer.zero();
    }
}


  
void
BaseMixedPressureElement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //set displacement and nonlocal location array
    answer.resize( this->giveNumberOfDofs(), this->giveNumberOfDofs() );
    answer.zero();

    FloatMatrix Kuu, Kup, Kpu, Kpp;
    //
    NLStructuralElement *elem = this->giveElement();
    double minusK;
    FloatArray Np, Dup;
    FloatMatrix Bu, DuuBu, DupNp, Duu;
    FloatMatrix answer_uu, answer_up, answer_pu, answer_pp;
    answer_uu.clear();
    answer_up.clear();
    answer_pu.clear();
    answer_pp.clear();
    for ( GaussPoint *gp : *elem->giveIntegrationRule(0) ) {
      if ( elem->giveGeometryMode() == 0 ) {
	this->computeSmallStrainMaterialStiffnessMatrix(Duu, Dup, minusK, rMode, gp, tStep);
	elem->computeBmatrixAt(gp, Bu);
      } else {
	this->computeLargeStrainMaterialStiffnessMatrix(Duu, Dup, minusK, rMode, gp, tStep);
	elem->computeBHmatrixAt(gp, Bu);
      }
      double dV = elem->computeVolumeAround(gp);
      DuuBu.beProductOf(Duu, Bu);
      answer_uu.plusProductUnsym(Bu, DuuBu, dV);
      //
      this->computePressureNMatrixAt(gp, Np);
      FloatMatrix mNp(Np, true);
      DupNp.beProductOf(Dup, mNp);
      answer_up.plusProductUnsym(Bu, DupNp, dV);
      //
      answer_pp.plusProductUnsym(mNp, mNp, dV / minusK);
    }   
    answer_pu.beTranspositionOf(answer_up);
    //
    answer.assemble(answer_uu, locationArray_u);
    answer.assemble(answer_up, locationArray_u, locationArray_p);
    answer.assemble(answer_pu, locationArray_p, locationArray_u);
    answer.assemble(answer_pp, locationArray_p);
}


void
BaseMixedPressureElement :: initializeFrom(InputRecord &ir)
{
    // @todo Is this function necessary???
}

void
BaseMixedPressureElement :: updateInternalState(TimeStep *tStep)
// Updates the receiver at end of step.
{
    FloatArray stress, strain;
    /*
     * // force updating strains & stresses
     * for ( auto &iRule: integrationRulesArray ) {
     *  for ( GaussPoint *gp: *iRule ) {
     *      this->computeStrainVector(strain, gp, tStep);
     *      this->computeStressVector(stress, strain, gp, tStep);
     *  }
     * }
     */
}


void
BaseMixedPressureElement :: postInitialize()
{
    IntArray IdMask_u, IdMask_p;
    this->giveDofManDofIDMask_u(IdMask_u);
    this->giveDofManDofIDMask_p(IdMask_p);
    this->giveLocationArrayOfDofIDs(locationArray_u, locationArray_p, EModelDefaultEquationNumbering(), IdMask_u, IdMask_p);
}
} // end namespace oofem
