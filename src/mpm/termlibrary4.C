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
#include "termlibrary2.h"
#include "termlibrary4.h"
#include "element.h"
#include "material.h"
#include "mathfem.h"
#include "crosssection.h"
#include "mpm/MagnetoElasticity/CrossSections/magnetoelasticcrosssection.h"
#include <ranges>
#include <vector>
namespace oofem {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MagnetoElasticity_GradGrad_Term :: MagnetoElasticity_GradGrad_Term (const Variable &testField, const Variable& unknownField, const Variable& magneticPotentialField) : Term(testField, unknownField), magneticPotentialField(magneticPotentialField) {}

int
MagnetoElasticity_GradGrad_Term :: computeGradientFields(std::vector<FloatArray>& vGrad, std::vector<FloatMatrix>& vB, MPElement& cell, const FloatArray& lcoords, MaterialMode mmode, TimeStep* tstep) const
{
    FloatArray u, vF;
    FloatMatrix Bu;
    //
    cell.getUnknownVector(u, this->field, VM_TotalIntrinsic, tstep);
    this->computeBHmatrixAt(Bu, this->field, this->field.interpolation, cell, lcoords, mmode);
    vF.beProductOf(Bu, u);
    // Deformation gradient F = gradU + I
    if ( mmode == _3dMat || mmode == _PlaneStrain ) {
      vF.at(1) += 1.0;
      vF.at(2) += 1.0;
      vF.at(3) += 1.0;
    } else if ( mmode == _PlaneStress ) {
      vF.at(1) += 1.0;
      vF.at(2) += 1.0;
    } else if ( mmode == _1dMat ) {
       vF.at(1) += 1.0;
    } else {
      OOFEM_ERROR("MaterialMode is not supported yet (%s)", __MaterialModeToString(mmode) );
    }
    /////////////////////////
    FloatArray phi, vH;
    FloatMatrix Bphi;
    //
    cell.getUnknownVector(phi, magneticPotentialField, VM_TotalIntrinsic, tstep);
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    FloatMatrix dNdx;
    this->magneticPotentialField.interpolation.evaldNdx(dNdx, lcoords, FEIElementGeometryWrapper(&cell));
    Bphi.beTranspositionOf(dNdx);
    // evaluate magnetic field at given point
    vH.beProductOf(Bphi, phi);
    vH.times(-1.);
    vGrad.emplace_back(vF);
    vGrad.emplace_back(vH);
    vB.emplace_back(Bu);
    vB.emplace_back(Bphi);
    return (u.giveSize() + phi.giveSize());
}




int
MagnetoElasticity_GradGrad_Term :: computeGradientFields_from_u(std::vector<FloatArray>& vGrad, std::vector<FloatMatrix>& vB, MPElement& cell, const FloatArray& lcoords, MaterialMode mmode, TimeStep* tstep, const FloatArray &u, const FloatArray &phi) const
{
    FloatArray  vF;
    FloatMatrix Bu;
    //
    this->computeBHmatrixAt(Bu, this->field, this->field.interpolation, cell, lcoords, mmode);
    vF.beProductOf(Bu, u);
    // Deformation gradient F = gradU + I
    if ( mmode == _3dMat || mmode == _PlaneStrain ) {
      vF.at(1) += 1.0;
      vF.at(2) += 1.0;
      vF.at(3) += 1.0;
    } else if ( mmode == _PlaneStress ) {
      vF.at(1) += 1.0;
      vF.at(2) += 1.0;
    } else if ( mmode == _1dMat ) {
       vF.at(1) += 1.0;
    } else {
      OOFEM_ERROR("MaterialMode is not supported yet (%s)", __MaterialModeToString(mmode) );
    }
    /////////////////////////
    FloatArray vH;
    FloatMatrix Bphi;
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    FloatMatrix dNdx;
    this->magneticPotentialField.interpolation.evaldNdx(dNdx, lcoords, FEIElementGeometryWrapper(&cell));
    Bphi.beTranspositionOf(dNdx);
    // evaluate magnetic field at given point
    vH.beProductOf(Bphi, phi);
    vH.times(-1.);    
    vGrad.emplace_back(vF);
    vGrad.emplace_back(vH);
    vB.emplace_back(Bu);
    vB.emplace_back(Bphi);
    return (u.giveSize() + phi.giveSize());
}




  
void
MagnetoElasticity_GradGrad_Term :: evaluate(FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const
{
  //
  std::vector<FloatArray> vGrad;
  std::vector<FloatMatrix> vB;
  // Grad contains deformation gradient and magnetic field
  auto size = this->computeGradientFields(vGrad, vB, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep);
  //
  auto cs = cell.giveCrossSection();
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *> (cs);
  auto vFlux = mcs->give_FluxVectors(vGrad, gp, tstep);
  //  for(auto& [B,Flux] : std::zip (vB, vFlux)) {
  answer.reserve(size);
  for(int i = 0; i < vB.size(); i++ ) {
	//auto& [B,Flux] : std::zip (vB, vFlux)) {
    FloatArray BF;
    BF.beTProductOf(vB.at(i), vFlux.at(i));
    answer.append(BF);
  }
}

void
MagnetoElasticity_GradGrad_Term :: compute_lin_num(FloatMatrix &D1, FloatMatrix &D2, FloatMatrix &D3, FloatMatrix &D4, MPElement &cell, GaussPoint *gp, TimeStep* tstep) const
{
  std::vector<FloatArray> vGrad;
  std::vector<FloatMatrix> vB;
  auto cs = cell.giveCrossSection();
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *> (cs);
  //
  auto size = this->computeGradientFields(vGrad, vB, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep);
  auto vFlux = mcs->give_FluxVectors(vGrad, gp, tstep);

  ////////////////////////////////////
  FloatArray u, phi;
  cell.getUnknownVector(u, this->field, VM_TotalIntrinsic, tstep);
  cell.getUnknownVector(phi, magneticPotentialField, VM_TotalIntrinsic, tstep);
  auto pert = 1.e-6;
  FloatArray u_pert, phi_pert;
  std::vector<FloatArray> vGrad_pert;
  std::vector<FloatMatrix> vB_pert;
  D1.resize(16,16);
  D2.resize(4,16);
  D3.resize(16,4);
  D4.resize(4,4);
  for(int i = 1; i <= u.giveSize(); i++) {
    // Grad contains deformation gradient and magnetic field
    u_pert = u;
    u_pert.at(i) += pert;
    auto size = this->computeGradientFields_from_u(vGrad_pert, vB_pert, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep, u_pert, phi);  //

    auto vFlux_pert = mcs->give_FluxVectors(vGrad_pert, gp, tstep);
    FloatArray BF1, BF2;
    auto dF1 = vFlux_pert.at(0) - vFlux.at(0);
    dF1.times(1./pert);
    BF1.beTProductOf(vB.at(0), dF1);
    auto dF2 = vFlux_pert.at(1) - vFlux.at(1);
    dF2.times(1./pert);
    BF2.beTProductOf(vB.at(1), dF2);
    D1.setColumn(BF1, i);
    D2.setColumn(BF2, i);
  }


  for(int i = 1; i <= phi.giveSize(); i++) {
    // Grad contains deformation gradient and magnetic field
    phi_pert = phi;
    phi_pert.at(i) += pert;
    auto size = this->computeGradientFields_from_u(vGrad_pert, vB, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep, u, phi_pert);  //
    auto cs = cell.giveCrossSection();
    auto mcs = dynamic_cast<MagnetoElasticCrossSection *> (cs);
    auto vFlux_pert = mcs->give_FluxVectors(vGrad_pert, gp, tstep);
    FloatArray BF1, BF2;
    auto dF1 = vFlux_pert.at(0) - vFlux.at(0);
    dF1.times(1./pert);
    BF1.beTProductOf(vB.at(0), dF1);
    auto dF2 = vFlux_pert.at(1) - vFlux.at(1);
    dF2.times(1./pert);
    BF2.beTProductOf(vB.at(1), dF2);
    //    
    D3.setColumn(BF1, i);
    D4.setColumn(BF2, i);
  }

  

}

void
MagnetoElasticity_GradGrad_Term :: evaluate_lin (FloatMatrix& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const
{
  std::vector<FloatArray> vGrad;
  std::vector<FloatMatrix> vB, dFdG;
  // Grad contains deformation gradient and magnetic field
  auto size = this->computeGradientFields(vGrad, vB, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep);
  FloatMatrix D1,D2,D3,D4;
  this->compute_lin_num(D1,D2,D3,D4, cell, gp, tstep);
  //
  auto cs = cell.giveCrossSection();
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *> (cs);
  mcs->give_dFlux_dGrad_Matrices(dFdG, TangentStiffness,gp, tstep);
  FloatMatrix BF, BFB;
  //
  answer.resize(size, size);
  int nr = 1, nc = 1;
  int nFields = 2;
  int index = 0;
  for(int i = 0; i < nFields; i++) {
    for(int j = 0; j < nFields; j++) {
      BF.beTProductOf(vB.at(i), dFdG.at(index));
      BFB.beProductOf(BF,vB.at(j));
      answer.setSubMatrix(BFB, nr, nc);
      nc += BFB.giveNumberOfColumns();
      index++;
    }
    nc = 1;
    nr += BFB.giveNumberOfRows();
  }

}



void
MagnetoElasticity_GradGrad_Term :: computeBHmatrixAt(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) const
{
    FloatMatrix dNdx;
    // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
    interpol.evaldNdx(dNdx, coords, FEIElementGeometryWrapper(&cell));
    
    if (mmode == _3dMat) {
      answer.resize(9, dNdx.giveNumberOfRows() * 3);
      answer.zero();
      for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);     // dv/dy
        answer.at(3, 3 * i - 0) = dNdx.at(i, 3);     // dw/dz
        answer.at(4, 3 * i - 1) = dNdx.at(i, 3);     // dv/dz
        answer.at(7, 3 * i - 0) = dNdx.at(i, 2);     // dw/dy
        answer.at(5, 3 * i - 2) = dNdx.at(i, 3);     // du/dz
        answer.at(8, 3 * i - 0) = dNdx.at(i, 1);     // dw/dx
        answer.at(6, 3 * i - 2) = dNdx.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dNdx.at(i, 1);     // dv/dx
      }
    } else if (mmode == _PlaneStrain) {
      answer.resize(5, dNdx.giveNumberOfRows() * 2);
      answer.zero();
      for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 2 * i - 1) = dNdx.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = dNdx.at(i, 2);     // dv/dy -2
        answer.at(4, 2 * i - 1) = dNdx.at(i, 2);     // du/dy -6
        answer.at(5, 2 * i - 0) = dNdx.at(i, 1);     // dv/dx -9
      }
    } else {
      OOFEM_ERROR("Unsupported material mode");
    }
}



  

} // end namespace oofem
