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
#include "material.h"
#include "mathfem.h"
#include "crosssection.h"
#include <ranges>
#include <vector>
#include "sm/CrossSections/thirdmediumcrosssection.h"

namespace oofem {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ThirdMedium_GradGrad_SecondGradientTerm ::ThirdMedium_GradGrad_SecondGradientTerm( const Variable &testField, const Variable &unknownField ) :
    StructuralTerm( testField, unknownField ) {}




int ThirdMedium_GradGrad_SecondGradientTerm ::computeSecondGradientField( FloatArray &vGradF, FloatMatrix &G, StructuralElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const
{
    FloatArray u;
    //
    cell.computeVectorOf(VM_Total, tstep, u);
    //cell.getUnknownVector(u, this->field, VM_TotalIntrinsic, tstep);
    this->computeGmatrixAt(G, this->field, this->field.interpolation, cell, lcoords, mmode);
    vGradF.beProductOf( G, u );
    /////////////////////////
    return u.giveSize();
}





  
void ThirdMedium_GradGrad_SecondGradientTerm ::evaluate( FloatArray &answer, StructuralElement &cell, GaussPoint *gp, TimeStep *tstep ) const
{
  //
  FloatArray vGrad;
  FloatMatrix B;
  // Grad contains deformation gradient and magnetic field
  auto size = this->computeSecondGradientField(vGrad, B, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep);
  //
  auto time = tstep->giveIntrinsicTime();

  auto cs = cell.giveCrossSection();
  auto tcs = dynamic_cast<ThirdMediumCrossSection *> (cs);
  if ( tcs == nullptr ) {
    OOFEM_ERROR( "Third medium terms require third medium cross sections." );
  }
  auto vFlux = tcs->give_SecondGradient_FluxVector(vGrad, gp, tstep);
  //  for(auto& [B,Flux] : std::zip (vB, vFlux)) {
  FloatArray BF;
  answer.beTProductOf(B, vFlux);
  
}


void ThirdMedium_GradGrad_SecondGradientTerm ::evaluate_lin( FloatMatrix &answer, StructuralElement &cell, GaussPoint *gp, TimeStep *tstep ) const
{
  FloatArray vGrad;
  FloatMatrix B, dFdG;
  // Grad contains deformation gradient and magnetic field
  auto size = this->computeSecondGradientField(vGrad, B, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep);
  //
  auto cs = cell.giveCrossSection();
  auto tcs = dynamic_cast<ThirdMediumCrossSection *>( cs );
  if ( tcs == nullptr ) {
    OOFEM_ERROR( "Third medium terms require third medium cross sections." );
  }
  tcs->give_SecondGradient_dFlux_dGrad(dFdG, TangentStiffness,gp, tstep);
  FloatMatrix BD;
  //
  answer.beTProductOf(B, B);
  answer.times(dFdG.at(1,1));
}



void ThirdMedium_GradGrad_SecondGradientTerm ::computeGmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const
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
	
     	  answer.at( 3, i * 2 - 0 ) = d2Ndx2.at( i, 3 );
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
ThirdMedium_GradGrad_FirstSecondGradientTerm::ThirdMedium_GradGrad_FirstSecondGradientTerm( const Variable &testField, const Variable &unknownField ) :
    StructuralTerm( testField, unknownField ) {}


int ThirdMedium_GradGrad_FirstSecondGradientTerm ::computeGradientField( FloatArray &vF, FloatMatrix &BH, StructuralElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const
{
  FloatArray u;
  //
  cell.computeVectorOf( VM_Total, tstep, u );
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

int ThirdMedium_GradGrad_FirstSecondGradientTerm ::computeSecondGradientField( FloatArray &vGradF, FloatMatrix &G, StructuralElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const
{
  FloatArray u;
  //
  cell.computeVectorOf( VM_Total, tstep, u );
  // cell.getUnknownVector(u, this->field, VM_TotalIntrinsic, tstep);
  this->computeGmatrixAt( G, this->field, this->field.interpolation, cell, lcoords, mmode );
  vGradF.beProductOf( G, u );
  /////////////////////////
  return u.giveSize();
}


void ThirdMedium_GradGrad_FirstSecondGradientTerm ::evaluate( FloatArray &answer, StructuralElement &cell, GaussPoint *gp, TimeStep *tstep ) const
{
  // Working with both deformation gradient (F) and its gradient (gradF);
  FloatArray vF, vGradF, GtimesT;
  FloatMatrix BH, G;
  //
  auto size1 = this->computeGradientField( vF, BH, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep );
  auto size2 = this->computeSecondGradientField(vGradF, G, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep);
  //
  auto cs = cell.giveCrossSection();
  auto tcs = dynamic_cast<ThirdMediumCrossSection *>( cs );
  if ( tcs == nullptr ) {
    OOFEM_ERROR( "Third medium terms require third medium cross sections." );
  }
  auto [vP, vT] = tcs->give_FirstSecondGradient_FluxVectors( vF, vGradF, gp, tstep );
  answer.beTProductOf( BH, vP );
  GtimesT.beTProductOf( G, vT );
  answer.add(GtimesT);
}


void ThirdMedium_GradGrad_FirstSecondGradientTerm ::evaluate_lin( FloatMatrix &answer, StructuralElement &cell, GaussPoint *gp, TimeStep *tstep ) const
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
  auto cs = cell.giveCrossSection();
  auto tcs = dynamic_cast<ThirdMediumCrossSection *>( cs );
  if ( tcs == nullptr ) {
    OOFEM_ERROR( "Third medium terms require third medium cross sections." );
  }
  auto [dPdF, dPdGradF, dTdF, dTdGradF] =  tcs->give_FirstSecondGradient_dFluxes_dGrads(TangentStiffness, gp, tstep );
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


void ThirdMedium_GradGrad_FirstSecondGradientTerm ::computeGmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const
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

void ThirdMedium_GradGrad_FirstSecondGradientTerm ::computeBHmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const
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
ThirdMedium_Grad_FbarTerm::ThirdMedium_Grad_FbarTerm( const Variable &testField, const Variable &unknownField ) :
    StructuralTerm( testField, unknownField ) {}


int ThirdMedium_Grad_FbarTerm ::computeGradientField( FloatArray &vF, FloatMatrix &BH, StructuralElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const
{
  FloatArray u;
  //
  cell.computeVectorOf( VM_Total, tstep, u );
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


void ThirdMedium_Grad_FbarTerm ::evaluate( FloatArray &answer, StructuralElement &cell, GaussPoint *gp, TimeStep *tstep ) const
{
  // Working with both deformation gradient (F) and its gradient (gradF);
  FloatArray vF, vFbar, GtimesT;
  FloatMatrix BH, Bbar;
  //
  auto size1 = this->computeGradientField( vF, BH, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep );
  FloatArray zeroCoords(gp->giveNaturalCoordinates());
  zeroCoords.times(0.0);
  auto size2 = this->computeGradientField(vFbar, Bbar, cell, zeroCoords, gp->giveMaterialMode(), tstep);
  //
  auto cs = cell.giveCrossSection();
  auto tcs = dynamic_cast<ThirdMediumCrossSection *>( cs );
  if ( tcs == nullptr ) {
    OOFEM_ERROR( "Third medium terms require third medium cross sections." );
  }
  auto [vP, vPbar] = tcs->give_Fbar_FluxVectors( vF, vFbar, gp, tstep );
  answer.beTProductOf( BH, vP );
  GtimesT.beTProductOf( Bbar, vPbar );
  answer.add(GtimesT);
}


void ThirdMedium_Grad_FbarTerm ::evaluate_lin( FloatMatrix &answer, StructuralElement &cell, GaussPoint *gp, TimeStep *tstep ) const
{
  FloatArray vFbar, vF; //useless since they are not passed on
  FloatMatrix BH, Bbar; //necessary, on the other hand - base function derivatives
  //FloatMatrix dPdF, dPdGradF, dTdF, dTdGradF; //constitutive matrices
  FloatMatrix B_dPdF_B, B_dPdFbar_Bbar, Bbar_dPbardF_B, Bbar_dPbardFbar_Bbar; // constitutive matrices multiplied by base function derivatives
  FloatMatrix B_dPdF, B_dPdFbar, Bbar_dPbardF, Bbar_dPbardFbar; //intermediate results
  //
  auto size1 = this->computeGradientField( vF, BH, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep );
  FloatArray zeroCoords(gp->giveNaturalCoordinates());
  zeroCoords.times(0.0);
  auto size2 = this->computeGradientField( vFbar, Bbar, cell, zeroCoords, gp->giveMaterialMode(), tstep );
  // get constitutive matrices
  auto cs = cell.giveCrossSection();
  auto tcs = dynamic_cast<ThirdMediumCrossSection *>( cs );
  if ( tcs == nullptr ) {
    OOFEM_ERROR( "Third medium terms require third medium cross sections." );
  }
  auto [dPdF, dPdFbar, dPbardF, dPbardFbar] =  tcs->give_Fbar_dFluxes_dGrads(TangentStiffness, gp, tstep );
  // construct result
  answer.resize(0,0);
  B_dPdF.beTProductOf(BH, dPdF);
  B_dPdF_B.beProductOf(B_dPdF, BH);
  answer.add(B_dPdF_B);

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

void ThirdMedium_Grad_FbarTerm ::computeBHmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const
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
ThirdMedium_Grad_LinearElasticTerm::ThirdMedium_Grad_LinearElasticTerm( const Variable &testField, const Variable &unknownField ) :
    StructuralTerm( testField, unknownField ) {}


int ThirdMedium_Grad_LinearElasticTerm ::computeDeformationField( FloatArray &eps, FloatMatrix &B, StructuralElement &cell, const FloatArray &lcoords, MaterialMode mmode, TimeStep *tstep ) const
{
  FloatArray u;
  //
  cell.computeVectorOf( VM_Total, tstep, u );

  this->computeBmatrixAt( B, this->field, this->field.interpolation, cell, lcoords, mmode );
  eps.beProductOf( B, u );
  
  return u.giveSize();
}


void ThirdMedium_Grad_LinearElasticTerm ::evaluate( FloatArray &answer, StructuralElement &cell, GaussPoint *gp, TimeStep *tstep ) const
{
  // Working with both deformation gradient (F) and its gradient (gradF);
  FloatArray eps;
  FloatMatrix B;
  //
  auto size = this->computeDeformationField( eps, B, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep );
  //
  auto cs = cell.giveCrossSection();
  auto tcs = dynamic_cast<ThirdMediumCrossSection *>( cs );
  if ( tcs == nullptr ) {
    OOFEM_ERROR( "Third medium terms require third medium cross sections." );
  }
  auto sigma = tcs->give_LinearElasticity_FluxVector( eps, gp, tstep );
  answer.beTProductOf( B, sigma );
}


void ThirdMedium_Grad_LinearElasticTerm ::evaluate_lin( FloatMatrix &answer, StructuralElement &cell, GaussPoint *gp, TimeStep *tstep ) const
{
  FloatArray eps; // useless since they are not passed on
  FloatMatrix B; // necessary, on the other hand - base function derivatives

  FloatMatrix B_dSigdEps;

  auto size = this->computeDeformationField( eps, B, cell, gp->giveNaturalCoordinates(), gp->giveMaterialMode(), tstep );

  // get constitutive matrices
  auto cs = cell.giveCrossSection();
  auto tcs = dynamic_cast<ThirdMediumCrossSection *>( cs );
  if ( tcs == nullptr ) {
    OOFEM_ERROR( "Third medium terms require third medium cross sections." );
  }
  auto dSigdEps = tcs->give_LinearElasticity_dFlux_dGrad( TangentStiffness, gp, tstep );
  // construct result
  answer.resize( 0, 0 );
  B_dSigdEps.beTProductOf( B, dSigdEps );
  answer.beProductOf( B_dSigdEps, B );
}

void ThirdMedium_Grad_LinearElasticTerm ::computeBmatrixAt( FloatMatrix &answer, const Variable &v, const FEInterpolation &interpol, const Element &cell, const FloatArray &coords, const MaterialMode mmode ) const
{

  FloatMatrix dNdx;
  // evaluate matrix of derivatives, the member at i,j position contains value of dNi/dxj
  interpol.evaldNdx( dNdx, coords, FEIElementGeometryWrapper( &cell ) );

  if ( mmode == _3dMat ) {
    answer.resize( 6, dNdx.giveNumberOfRows() * 3 );
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
      answer.at( 1, 3 * i - 2 ) = dNdx.at( i, 1 );
      answer.at( 2, 3 * i - 1 ) = dNdx.at( i, 2 );
      answer.at( 3, 3 * i - 0 ) = dNdx.at( i, 3 );

      answer.at( 5, 3 * i - 2 ) = answer.at( 4, 3 * i - 1 ) = dNdx.at( i, 3 );
      answer.at( 6, 3 * i - 2 ) = answer.at( 4, 3 * i - 0 ) = dNdx.at( i, 2 );
      answer.at( 6, 3 * i - 1 ) = answer.at( 5, 3 * i - 0 ) = dNdx.at( i, 1 );
    }
  } else if ( mmode == _PlaneStrain ) {

    answer.resize( 4, dNdx.giveNumberOfRows() * 2 );
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
      answer.at( 1, i * 2 - 1 ) = dNdx.at( i, 1 );
      answer.at( 2, i * 2 - 0 ) = dNdx.at( i, 2 );

      answer.at( 4, 2 * i - 1 ) = dNdx.at( i, 2 );
      answer.at( 4, 2 * i - 0 ) = dNdx.at( i, 1 );
    }

  } else {
    OOFEM_ERROR( "Unsupported material mode" );
  }
}  

} // end namespace oofem
