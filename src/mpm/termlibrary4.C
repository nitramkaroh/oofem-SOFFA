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
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *>( cell.giveCrossSection() );
  if ( mcs == nullptr ) {
    OOFEM_ERROR( "MagnetoElastic terms require magnetoelastic cross sections." );
  }
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
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *>( cell.giveCrossSection() );
  if ( mcs == nullptr ) {
    OOFEM_ERROR( "MagnetoElastic terms require magnetoelastic cross sections." );
  }
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
    auto mcs = dynamic_cast<MagnetoElasticCrossSection *>( cell.giveCrossSection() );
    if ( mcs == nullptr ) {
      OOFEM_ERROR( "MagnetoElastic terms require magnetoelastic cross sections." );
    }
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
  //  this->compute_lin_num(D1,D2,D3,D4, cell, gp, tstep);
  //
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *>( cell.giveCrossSection() );
  if ( mcs == nullptr ) {
    OOFEM_ERROR( "MagnetoElastic terms require magnetoelastic cross sections." );
  }
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MagnetoElasticity_GradGrad_SecondGradientTerm :: MagnetoElasticity_GradGrad_SecondGradientTerm (const Variable &testField, const Variable& unknownField) : Term(testField, unknownField) {}




int MagnetoElasticity_GradGrad_SecondGradientTerm ::computeSecondGradientField( FloatArray &vGradF, FloatMatrix &G, MPElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const
{
    FloatArray u;
    //
    cell.getUnknownVector(u, this->field, VM_TotalIntrinsic, tstep);
    this->computeGmatrixAt(G, this->field, this->field.interpolation, cell, lcoords, mmode);
    vGradF.beProductOf( G, u );
    /////////////////////////
    return u.giveSize();
}





  
void
MagnetoElasticity_GradGrad_SecondGradientTerm :: evaluate(FloatArray& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const
{
  //
  FloatArray vGrad;
  FloatMatrix B;
  // Grad contains deformation gradient and magnetic field
  auto size = this->computeSecondGradientField(vGrad, B, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep);
  //
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *>( cell.giveCrossSection() );
  if ( mcs == nullptr ) {
    OOFEM_ERROR( "MagnetoElastic terms require magnetoelastic cross sections." );
  }
  auto vFlux = mcs->give_SecondGradient_FluxVector(vGrad, gp, tstep);

  answer.beTProductOf(B, vFlux);  
}


void
MagnetoElasticity_GradGrad_SecondGradientTerm :: evaluate_lin (FloatMatrix& answer, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const
{
  FloatArray vGrad;
  FloatMatrix B, dFdG;
  // Grad contains deformation gradient and magnetic field
  auto size = this->computeSecondGradientField(vGrad, B, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep);
  //
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *>( cell.giveCrossSection() );
  if ( mcs == nullptr ) {
    OOFEM_ERROR( "MagnetoElastic terms require magnetoelastic cross sections." );
  }
  mcs->give_SecondGradient_dFlux_dGrad(dFdG, TangentStiffness,gp, tstep);
  FloatMatrix BD;
  //
  answer.beTProductOf(B, B);
  answer.times(dFdG.at(1,1));
}



void
MagnetoElasticity_GradGrad_SecondGradientTerm :: computeGmatrixAt(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) const
{

    FloatMatrix d2Ndx2; 
    interpol.evald2Ndx2( d2Ndx2, coords, FEIElementGeometryWrapper(&cell) );
    //
    if (mmode == _PlaneStrain) {
      answer.resize(8,   d2Ndx2.giveNumberOfRows() * 2);
      answer.zero();
      for ( int i = 1; i <= d2Ndx2.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = d2Ndx2.at(i, 1);
        answer.at(2, i * 2 - 1) = d2Ndx2.at(i, 3);
	
     	answer.at(3, i * 2 - 0) = d2Ndx2.at(i, 3);
        answer.at(4, i * 2 - 0) = d2Ndx2.at(i, 2);

        answer.at(5, i * 2 - 1) = d2Ndx2.at(i, 3);
        answer.at(6, i * 2 - 1) = d2Ndx2.at(i, 2);

        answer.at(7, i * 2 - 0) = d2Ndx2.at(i, 1);
        answer.at(8, i * 2 - 0) = d2Ndx2.at(i, 3);
      }     
    } else {
      OOFEM_ERROR("Unsupported material mode");
    }  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MagnetoElasticity_GradGrad_FirstSecondGradientTerm ::MagnetoElasticity_GradGrad_FirstSecondGradientTerm( const Variable &testField, const Variable &unknownField) :
    Term( testField, unknownField ) {}


int MagnetoElasticity_GradGrad_FirstSecondGradientTerm ::computeGradientField( FloatArray &vF, FloatMatrix &BH, MPElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const
{
  FloatArray u;
  //
  cell.getUnknownVector( u, this->field, VM_TotalIntrinsic, tstep );
  this->computeBHmatrixAt( BH, this->field, this->field.interpolation, cell, lcoords, mmode );
  vF.beProductOf( BH, u );
  // Deformation gradient F = gradU + I
  if ( mmode == _3dMat || mmode == _PlaneStrain ) {
    vF.at( 1 ) += 1.0;
    vF.at( 2 ) += 1.0;
    vF.at( 3 ) += 1.0;
  } else if ( mmode == _PlaneStress ) {
    vF.at( 1 ) += 1.0;
    vF.at( 2 ) += 1.0;
  } else if ( mmode == _1dMat ) {
    vF.at( 1 ) += 1.0;
  } else {
    OOFEM_ERROR( "MaterialMode is not supported yet (%s)", __MaterialModeToString( mmode ) );
  }

  return u.giveSize();
}

int MagnetoElasticity_GradGrad_FirstSecondGradientTerm ::computeSecondGradientField( FloatArray &vGradF, FloatMatrix &G, MPElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const
{
  FloatArray u;
  //
  cell.getUnknownVector( u, this->field, VM_TotalIntrinsic, tstep );
  this->computeGmatrixAt( G, this->field, this->field.interpolation, cell, lcoords, mmode );
  vGradF.beProductOf( G, u );
  /////////////////////////
  return u.giveSize();
}


void MagnetoElasticity_GradGrad_FirstSecondGradientTerm ::evaluate( FloatArray &answer, MPElement &cell, GaussPoint *gp, TimeStep *tstep ) const
{
  // Working with both deformation gradient (F) and its gradient (gradF);
  FloatArray vF, vGradF, GtimesT;
  FloatMatrix BH, G;
  //
  auto size1 = this->computeGradientField( vF, BH, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep );
  auto size2 = this->computeSecondGradientField(vGradF, G, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep);
  //
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *>( cell.giveCrossSection() );
  if ( mcs == nullptr ) {
    OOFEM_ERROR( "MagnetoElastic terms require magnetoelastic cross sections." );
  }
  auto [vP, vT] = mcs->give_FirstSecondGradient_FluxVectors( vF, vGradF, gp, tstep );
  answer.beTProductOf( BH, vP );
  GtimesT.beTProductOf( G, vT );
  answer.add(GtimesT);
}


void MagnetoElasticity_GradGrad_FirstSecondGradientTerm ::evaluate_lin( FloatMatrix &answer, MPElement &cell, GaussPoint *gp, TimeStep *tstep ) const
{
  FloatArray vGradF, vF; //useless since they are not passed on
  FloatMatrix BH, G; //necessary, on the other hand - base function derivatives
  //FloatMatrix dPdF, dPdGradF, dTdF, dTdGradF; //constitutive matrices
  FloatMatrix B_dPdF_B, B_dPdGradF_G, G_dTdF_B, G_dTdGradF_G; // constitutive matrices multiplied by base function derivatives
  FloatMatrix B_dPdF, B_dPdGradF, G_dTdF, G_dTdGradF; //intermediate results
  //
  auto size1 = this->computeGradientField( vF, BH, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep );
  auto size2 = this->computeSecondGradientField( vGradF, G, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep );
  // get constitutive matrices
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *>( cell.giveCrossSection() );
  if ( mcs == nullptr ) {
    OOFEM_ERROR( "MagnetoElastic terms require magnetoelastic cross sections." );
  }
  //mcs->give_FirstSecondGradient_dFluxes_dGrads( std::make_tuple( dPdF, dPdGradF, dTdF, dTdGradF ), TangentStiffness, gp, tstep );
  auto [dPdF, dPdGradF, dTdF, dTdGradF] = mcs->give_FirstSecondGradient_dFluxes_dGrads( TangentStiffness, gp, tstep );
  // construct result
  answer.resize(0,0);
  B_dPdF.beTProductOf(BH, dPdF);
  B_dPdF_B.beProductOf(B_dPdF, BH);
  answer.add(B_dPdF_B);

  B_dPdGradF.beTProductOf( BH, dPdGradF );
  B_dPdGradF_G.beProductOf( B_dPdGradF, G );
  answer.add( B_dPdGradF_G );

  G_dTdF.beTProductOf( G, dTdF );
  G_dTdF_B.beProductOf( G_dTdF, BH );
  answer.add( G_dTdF_B );


  G_dTdGradF.beTProductOf( G, dTdGradF );
  G_dTdGradF_G.beProductOf( G_dTdGradF, G );
  answer.add( G_dTdGradF_G );
}


void MagnetoElasticity_GradGrad_FirstSecondGradientTerm ::computeGmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const
{

  FloatMatrix d2Ndx2;
  interpol.evald2Ndx2( d2Ndx2, coords, FEIElementGeometryWrapper( &cell ) );
  //
  if ( mmode == _PlaneStrain ) {
    answer.resize( 8, d2Ndx2.giveNumberOfRows() * 2 );
    answer.zero();
    for ( int i = 1; i <= d2Ndx2.giveNumberOfRows(); i++ ) {
      answer.at( 1, i * 2 - 1 ) = d2Ndx2.at( i, 1 );
      answer.at( 2, i * 2 - 1 ) = d2Ndx2.at( i, 3 );

      answer.at( 3, i * 2 - 0 ) = d2Ndx2.at( i, 3 );
      answer.at( 4, i * 2 - 0 ) = d2Ndx2.at( i, 2 );

      answer.at( 5, i * 2 - 1 ) = d2Ndx2.at( i, 3 );
      answer.at( 6, i * 2 - 1 ) = d2Ndx2.at( i, 2 );

      answer.at( 7, i * 2 - 0 ) = d2Ndx2.at( i, 1 );
      answer.at( 8, i * 2 - 0 ) = d2Ndx2.at( i, 3 );
    }
  } else {
    OOFEM_ERROR( "Unsupported material mode" );
  }
}

void MagnetoElasticity_GradGrad_FirstSecondGradientTerm ::computeBHmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const
{

  FloatMatrix dNdx;
  // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
  interpol.evaldNdx( dNdx, coords, FEIElementGeometryWrapper( &cell ) );

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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
MagnetoElasticity_GradGrad_FbarTerm::MagnetoElasticity_GradGrad_FbarTerm( const Variable &testField, const Variable &unknownField ) :
    Term( testField, unknownField ) {}


int MagnetoElasticity_GradGrad_FbarTerm ::computeGradientField( FloatArray &vF, FloatMatrix &BH, MPElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const
{
  FloatArray u;
  //
  //cell.computeVectorOf( VM_Total, tstep, u );
  cell.getUnknownVector( u, this->field, VM_TotalIntrinsic, tstep );
  // cell.getUnknownVector(u, this->field, VM_TotalIntrinsic, tstep);
  this->computeBHmatrixAt( BH, this->field, this->field.interpolation, cell, lcoords, mmode );
  vF.beProductOf( BH, u );
  // Deformation gradient F = gradU + I
  if ( mmode == _3dMat || mmode == _PlaneStrain ) {
    vF.at( 1 ) += 1.0;
    vF.at( 2 ) += 1.0;
    vF.at( 3 ) += 1.0;
  } else if ( mmode == _PlaneStress ) {
    vF.at( 1 ) += 1.0;
    vF.at( 2 ) += 1.0;
  } else if ( mmode == _1dMat ) {
    vF.at( 1 ) += 1.0;
  } else {
    OOFEM_ERROR( "MaterialMode is not supported yet (%s)", __MaterialModeToString( mmode ) );
  }

  return u.giveSize();
}


void MagnetoElasticity_GradGrad_FbarTerm ::evaluate( FloatArray &answer, MPElement &cell, GaussPoint *gp, TimeStep *tstep ) const
{
  // Working with both deformation gradient (F) and its gradient (gradF);
  FloatArray vF, vFbar, BbarTimesPbar;
  FloatMatrix BH, Bbar;
  //
  auto size1 = this->computeGradientField( vF, BH, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep );
  FloatArray zeroCoords( gp->giveNaturalCoordinates() );
  zeroCoords.times( 0.0 );
  auto size2 = this->computeGradientField( vFbar, Bbar, cell, zeroCoords, gp->giveMaterialMode(), tstep );
  //
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *>( cell.giveCrossSection() );
  if ( mcs == nullptr ) {
    OOFEM_ERROR( "MagnetoElastic terms require magnetoelastic cross sections." );
  }
  auto [vP, vPbar] = mcs->give_Fbar_FluxVectors( vF, vFbar, gp, tstep );
  answer.beTProductOf( BH, vP );
  BbarTimesPbar.beTProductOf( Bbar, vPbar );
  answer.add( BbarTimesPbar );
}


void MagnetoElasticity_GradGrad_FbarTerm ::evaluate_lin( FloatMatrix &answer, MPElement &cell, GaussPoint *gp, TimeStep *tstep ) const
{
  FloatArray vFbar, vF; // useless since they are not passed on
  FloatMatrix BH, Bbar; // necessary, on the other hand - base function derivatives
  // FloatMatrix dPdF, dPdGradF, dTdF, dTdGradF; //constitutive matrices
  FloatMatrix B_dPdF_B, B_dPdFbar_Bbar, Bbar_dPbardF_B, Bbar_dPbardFbar_Bbar; // constitutive matrices multiplied by base function derivatives
  FloatMatrix B_dPdF, B_dPdFbar, Bbar_dPbardF, Bbar_dPbardFbar; // intermediate results
  //
  auto size1 = this->computeGradientField( vF, BH, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep );
  auto size2 = this->computeGradientField( vFbar, Bbar, cell, { 0, 0 }, gp->giveMaterialMode(), tstep );
  // get constitutive matrices
  auto mcs = dynamic_cast<MagnetoElasticCrossSection *>( cell.giveCrossSection() );
  if ( mcs == nullptr ) {
    OOFEM_ERROR( "MagnetoElastic terms require magnetoelastic cross sections." );
  }
  if ( mcs == nullptr ) {
    OOFEM_ERROR( "MagnetoElastic terms require magnetoelastic cross sections." );
  }
  auto [dPdF, dPdFbar, dPbardF, dPbardFbar] = mcs->give_Fbar_dFluxes_dGrads( TangentStiffness, gp, tstep );
  // construct result
  answer.resize( 0, 0 );
  B_dPdF.beTProductOf( BH, dPdF );
  B_dPdF_B.beProductOf( B_dPdF, BH );
  answer.add( B_dPdF_B );

  B_dPdFbar.beTProductOf( BH, dPdFbar );
  B_dPdFbar_Bbar.beProductOf( B_dPdFbar, Bbar );
  answer.add( B_dPdFbar_Bbar );

  Bbar_dPbardF.beTProductOf( Bbar, dPbardF );
  Bbar_dPbardF_B.beProductOf( Bbar_dPbardF, BH );
  answer.add( Bbar_dPbardF_B );

  Bbar_dPbardFbar.beTProductOf( Bbar, dPbardFbar );
  Bbar_dPbardFbar_Bbar.beProductOf( Bbar_dPbardFbar, Bbar );
  answer.add( Bbar_dPbardFbar_Bbar );
}

void MagnetoElasticity_GradGrad_FbarTerm ::computeBHmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const
{

  FloatMatrix dNdx;
  // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
  interpol.evaldNdx( dNdx, coords, FEIElementGeometryWrapper( &cell ) );

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

  

} // end namespace oofem
