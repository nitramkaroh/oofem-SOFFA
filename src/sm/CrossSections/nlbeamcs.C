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

#include "sm/CrossSections/nlbeamcs.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "material.h"
#include "floatarray.h"
#include "verbose.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "mathfem.h"

namespace oofem {
REGISTER_CrossSection(NlBeamSimpleElasticCrossSection);
REGISTER_CrossSection(NlBeamSimpleElastoPlasticCrossSection);

void
NlBeamSimpleElasticCrossSection::initializeFrom(InputRecord &ir)
{
    CrossSection::initializeFrom(ir);
    IR_GIVE_FIELD(ir, Iy, "iy");
    IR_GIVE_FIELD(ir, A, "area");
    IR_GIVE_OPTIONAL_FIELD(ir, gamma, "gamma");
}


double
NlBeamSimpleElasticCrossSection :: computeMomentFromCurvature(double kappa)
{
  return E*Iy*kappa;

}


double
NlBeamSimpleElasticCrossSection :: computeDerivativeOfMomentWrtCurvature(double kappa)  
{
  return E*Iy;
}


  
double
NlBeamSimpleElasticCrossSection :: computeCurvatureFromMoment(double M)
{
  return M/(E*Iy);
}
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////



void
NlBeamSimpleElastoPlasticCrossSection::initializeFrom(InputRecord &ir)
{
    NlBeamSimpleElasticCrossSection::initializeFrom(ir);
    IR_GIVE_FIELD(ir, M_elastic, "mel");
    IR_GIVE_FIELD(ir, M_plastic, "mpl");
    IR_GIVE_OPTIONAL_FIELD(ir, delta, "delta");

    
}

void 
NlBeamSimpleElastoPlasticCrossSection::postInitialize()
{
    kappa_elastic = M_elastic/E/Iy;
    delta_kappa = 2. * (M_plastic-M_elastic)/E/Iy/M_PI;
    kappa_nl = kappa_elastic + delta_kappa * sqrt(1./delta-1.);
    M_nl = M_elastic + E *Iy * delta_kappa * atan(sqrt(1./delta-1.));  
 
}



  
double
NlBeamSimpleElastoPlasticCrossSection :: computeMomentFromCurvature(double kappa)
{
  double M;
  if(fabs(kappa) <= kappa_elastic) {
    M = E * Iy * kappa;
  } else if(fabs(kappa) <= kappa_nl) {
    M = signum(kappa) * (M_elastic + E * Iy * delta_kappa * atan(fabs(kappa - kappa_elastic)/delta_kappa));
  } else {
    M = signum(kappa) * (M_nl + delta * E * Iy * (fabs(kappa) - kappa_nl));
  }
  return M;
}


double
NlBeamSimpleElastoPlasticCrossSection :: computeDerivativeOfMomentWrtCurvature(double kappa)  
{
    kappa_elastic = M_elastic/E/Iy;
    delta_kappa = 2. * (M_plastic-M_elastic)/E/Iy/M_PI;
    kappa_nl = kappa_elastic + delta_kappa * sqrt(1./delta-1.);
    M_nl = M_elastic + E *Iy * delta_kappa * atan(sqrt(1./delta-1.));  
    //
    double dMdk;
  
  if(fabs(kappa) <= kappa_elastic) {
    dMdk = E * Iy;
  } else if(fabs(kappa) <= kappa_nl) {
    dMdk = E * Iy / ( 1 + (fabs(kappa)-kappa_elastic) / delta_kappa * (fabs(kappa)-kappa_elastic)/delta_kappa );
  } else {
    dMdk = delta * E * Iy;
  }
  return dMdk;
}


  
double
NlBeamSimpleElastoPlasticCrossSection :: computeCurvatureFromMoment(double M)
{
    kappa_elastic = M_elastic/E/Iy;
    delta_kappa = 2. * (M_plastic-M_elastic)/E/Iy/M_PI;
    kappa_nl = kappa_elastic + delta_kappa * sqrt(1./delta-1.);
    M_nl = M_elastic + E *Iy * delta_kappa * atan(sqrt(1./delta-1.));  
    //  
  double kappa;
  if(fabs(M) <= M_elastic) {
    kappa = M / (E * Iy);
  } else if(fabs(M) <= M_nl) {
    kappa = signum(M) * (kappa_elastic + delta_kappa * tan((fabs(M) - M_elastic) / (E*Iy*delta_kappa)));
  } else {
    kappa = signum(M) * (kappa_nl + (fabs(M) - M_elastic) / (E*Iy*delta));
  }
  return kappa;

}


  


} // end namespace oofem
