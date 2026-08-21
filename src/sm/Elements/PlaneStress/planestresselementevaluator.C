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

#include "sm/Elements/PlaneStress/planestresselementevaluator.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "domain.h"
#include "node.h"
#include "element.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "matresponsemode.h"
#include "crosssection.h"
#include "mathfem.h"
#include "iga/iga.h"

namespace oofem {
void PlaneStressStructuralElementEvaluator :: computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray N;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    answer.beNMatrixOf(N, 2);
}

void PlaneStressStructuralElementEvaluator :: computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix d;

    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evaldNdx( d, gp->giveNaturalCoordinates(),
                     FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );

    answer.resize(3, d.giveNumberOfRows() * 2);
    answer.zero();

    for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
        answer.at(1, i * 2 - 1) = d.at(i, 1);
        answer.at(2, i * 2 - 0) = d.at(i, 2);

        answer.at(3, 2 * i - 1) = d.at(i, 2);
        answer.at(3, 2 * i - 0) = d.at(i, 1);
    }
}

void PlaneStressStructuralElementEvaluator::computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    // Returns the [ 4 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
    // evaluated at gp.
    /// @todo not checked if correct

    FloatMatrix dNdx;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(),
    FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );


    answer.resize( 4, dNdx.giveNumberOfRows() * 2 );
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at( 1, 2 * i - 1 ) = dNdx.at( i, 1 ); // du/dx -1
        answer.at( 2, 2 * i - 0 ) = dNdx.at( i, 2 ); // dv/dy -2
        answer.at( 3, 2 * i - 1 ) = dNdx.at( i, 2 ); // du/dy -6
        answer.at( 4, 2 * i - 0 ) = dNdx.at( i, 1 ); // dv/dx -9
    }
}


double PlaneStressStructuralElementEvaluator :: computeVolumeAround(GaussPoint *gp)
{
    double determinant, weight, thickness, volume;
    determinant = fabs( this->giveElement()->giveInterpolation()
                       ->giveTransformationJacobian( gp->giveNaturalCoordinates(),
                                                    FEIIGAElementGeometryWrapper( this->giveElement(),
                                                                                 gp->giveIntegrationRule()->giveKnotSpan() ) ) );
    weight      = gp->giveWeight();
    thickness   = this->giveElement()->giveCrossSection()->give(CS_Thickness, gp);
    volume      = determinant * weight * thickness;

    return volume;
}


void PlaneStressStructuralElementEvaluator :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->giveRealStress_PlaneStress(strain, gp, tStep);
}

void PlaneStressStructuralElementEvaluator :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = static_cast< StructuralCrossSection * >( this->giveElement()->giveCrossSection() )->giveStiffnessMatrix_PlaneStress(rMode, gp, tStep);
}

void PlaneStressStructuralElementEvaluator ::initializeFrom( InputRecord &ir )
{
    StructuralElementEvaluator::initializeFrom( ir );
}

void PlaneStressStructuralElementEvaluator::computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    answer = static_cast<StructuralCrossSection *>( this->giveElement()->giveCrossSection() )->giveStiffnessMatrix_dPdF_PlaneStress( rMode, gp, tStep );
}


/////  Plane Strain
void PlaneStrainStructuralElementEvaluator ::computeNMatrixAt( FloatMatrix &answer, GaussPoint *gp )
{
    FloatArray N;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    answer.beNMatrixOf( N, 2 );
}

void PlaneStrainStructuralElementEvaluator ::computeBMatrixAt( FloatMatrix &answer, GaussPoint *gp )
{
    FloatMatrix d;

    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evaldNdx( d, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    
    answer.resize( 4, d.giveNumberOfRows() * 2 );
    answer.zero();

    for ( int i = 1; i <= d.giveNumberOfRows(); i++ ) {
        answer.at( 1, i * 2 - 1 ) = d.at( i, 1 );
        answer.at( 2, i * 2 - 0 ) = d.at( i, 2 );

        answer.at( 4, 2 * i - 1 ) = d.at( i, 2 );
        answer.at( 4, 2 * i - 0 ) = d.at( i, 1 );
    }
}

void PlaneStrainStructuralElementEvaluator::computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    // Returns the [ 4 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
    // evaluated at gp.
    /// @todo not checked if correct

    FloatMatrix dNdx;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );


    answer.resize( 5, dNdx.giveNumberOfRows() * 2 );
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at( 1, 2 * i - 1 ) = dNdx.at( i, 1 ); // du/dx -1
        answer.at( 2, 2 * i - 0 ) = dNdx.at( i, 2 ); // dv/dy -2
        answer.at( 4, 2 * i - 1 ) = dNdx.at( i, 2 ); // du/dy -6
        answer.at( 5, 2 * i - 0 ) = dNdx.at( i, 1 ); // dv/dx -9
    }
}


double PlaneStrainStructuralElementEvaluator ::computeVolumeAround( GaussPoint *gp )
{
    double determinant, weight, thickness, volume;
    determinant = fabs( this->giveElement()->giveInterpolation()->giveTransformationJacobian( gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( this->giveElement(),
            gp->giveIntegrationRule()->giveKnotSpan() ) ) );
    weight      = gp->giveWeight();
    thickness   = this->giveElement()->giveCrossSection()->give( CS_Thickness, gp );
    volume      = determinant * weight * thickness;

    return volume;
}


void PlaneStrainStructuralElementEvaluator ::computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep )
{
    answer = static_cast<StructuralCrossSection *>( this->giveElement()->giveCrossSection() )->giveRealStress_PlaneStrain( strain, gp, tStep );
}

void PlaneStrainStructuralElementEvaluator ::computeConstitutiveMatrixAt( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    answer = static_cast<StructuralCrossSection *>( this->giveElement()->giveCrossSection() )->giveStiffnessMatrix_PlaneStrain( rMode, gp, tStep );
}

void PlaneStrainStructuralElementEvaluator ::initializeFrom( InputRecord &ir )
{
    StructuralElementEvaluator::initializeFrom( ir );
}

void PlaneStrainStructuralElementEvaluator::computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    answer = static_cast<StructuralCrossSection *>( this->giveElement()->giveCrossSection() )->giveStiffnessMatrix_dPdF_PlaneStrain( rMode, gp, tStep );
}


////////////////
// second order stuff

void PlaneStrainStructuralElementEvaluator::computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    FloatMatrix ders;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();

    // returns index where nonzero basis functions start (1-based)
    // ders is a matrix of size [nno x 6] where nno is the number of nonzero basis functions at the gauss point, the columns are: [N, dN/drho, dN/dxi, d2N/drho2, d2N/dxi2, d2N/drhodxi]
    auto knotSpan  = gp->giveIntegrationRule()->giveKnotSpan();
    int indexStart = interp->evalSecondDerivatives( ders, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), knotSpan ) );

    // Get nodal coordinates
    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );
    FloatArray nodesCoords;
    elem->giveNodalCoordinates( nodesCoords );

    //double R     = 0.; // undeformed radius (X1)
    double R_rho = 0., R_xi = 0., Y_rho = 0., Y_xi = 0.; // 1st derivatives
    double R_rho_rho = 0., R_xi_xi = 0., R_rho_xi = 0., Y_rho_rho = 0., Y_xi_xi = 0., Y_rho_xi = 0.; // 2nd derivatives

    IntArray mask;
    interp->giveKnotSpanBasisFuncMask( *knotSpan, mask );

    // Evaluate reference coordinate derivatives
    for ( int k = 0; k < ders.giveNumberOfRows(); k++ ) {
        // int indCords         = indexStart + k;
        int indCords         = mask.at( k + 1 );
        double vertexCoord_1 = nodesCoords.at( 2 * indCords - 1 );
        double vertexCoord_2 = nodesCoords.at( 2 * indCords );
        // position
        //R += vertexCoord_1 * ders( k, 0 ); // compute radial coordinate at the gauss point
        // first derivatives
        R_rho += ders( k, 1 ) * vertexCoord_1;
        Y_rho += ders( k, 1 ) * vertexCoord_2;

        R_xi += ders( k, 2 ) * vertexCoord_1;
        Y_xi += ders( k, 2 ) * vertexCoord_2;
        // second derivatives
        R_rho_rho += ders( k, 3 ) * vertexCoord_1;
        Y_rho_rho += ders( k, 3 ) * vertexCoord_2;

        R_xi_xi += ders( k, 4 ) * vertexCoord_1;
        Y_xi_xi += ders( k, 4 ) * vertexCoord_2;

        R_rho_xi += ders( k, 5 ) * vertexCoord_1;
        Y_rho_xi += ders( k, 5 ) * vertexCoord_2;
    }
    double R_xi_rho = R_rho_xi;
    double Y_xi_rho = Y_rho_xi;

    // Create part proportional to individual coordinates in cylindrical system
    Tensor3_3d dFcyl_r_rho_rho;
    Tensor3_3d dFcyl_r_xi_rho;
    Tensor3_3d dFcyl_r_xi_xi;
    Tensor3_3d dFcyl_r_rho;
    Tensor3_3d dFcyl_r_xi;
    Tensor3_3d dFcyl_r;
    Tensor3_3d dFcyl_y_rho_rho;
    Tensor3_3d dFcyl_y_xi_rho;
    Tensor3_3d dFcyl_y_xi_xi;
    Tensor3_3d dFcyl_y_rho;
    Tensor3_3d dFcyl_y_xi;

    const double J  = R_rho * Y_xi - R_xi * Y_rho;
    const double J3 = J * J * J;
    //const double RJ = R * J;

    // r_rho_rho
    dFcyl_r_rho_rho( 0, 0, 0 ) = ( R_rho * Y_xi * Y_xi * Y_xi - R_xi * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_r_rho_rho( 0, 1, 0 ) = -( -Y_rho * R_xi * R_xi * Y_xi + R_rho * R_xi * Y_xi * Y_xi ) / J3;
    dFcyl_r_rho_rho( 0, 0, 1 ) = dFcyl_r_rho_rho( 0, 1, 0 );
    dFcyl_r_rho_rho( 0, 1, 1 ) = ( R_rho * R_xi * R_xi * Y_xi - R_xi * R_xi * R_xi * Y_rho ) / J3;

    // r_xi_rho
    dFcyl_r_xi_rho( 0, 0, 0 ) = -( -2.0 * R_xi * Y_rho * Y_rho * Y_xi + 2.0 * R_rho * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_r_xi_rho( 0, 1, 0 ) = ( R_rho * R_rho * Y_xi * Y_xi - R_xi * R_xi * Y_rho * Y_rho ) / J3;
    dFcyl_r_xi_rho( 0, 0, 1 ) = dFcyl_r_xi_rho( 0, 1, 0 );
    dFcyl_r_xi_rho( 0, 1, 1 ) = ( 2.0 * R_rho * R_xi * R_xi * Y_rho - 2.0 * R_rho * R_rho * R_xi * Y_xi ) / J3;

    // r_xi_xi
    dFcyl_r_xi_xi( 0, 0, 0 ) = -( R_xi * Y_rho * Y_rho * Y_rho - R_rho * Y_xi * Y_rho * Y_rho ) / J3;
    dFcyl_r_xi_xi( 0, 1, 0 ) = ( -Y_xi * R_rho * R_rho * Y_rho + R_xi * R_rho * Y_rho * Y_rho ) / J3;
    dFcyl_r_xi_xi( 0, 0, 1 ) = dFcyl_r_xi_xi( 0, 1, 0 );
    dFcyl_r_xi_xi( 0, 1, 1 ) = ( R_rho * R_rho * R_rho * Y_xi - R_rho * R_rho * R_xi * Y_rho ) / J3;

    // r_rho
    dFcyl_r_rho( 0, 0, 0 ) = ( R_xi * ( Y_xi_xi * Y_rho * Y_rho - 2.0 * Y_xi_rho * Y_rho * Y_xi + Y_rho_rho * Y_xi * Y_xi ) - R_rho_rho * Y_xi * Y_xi * Y_xi - R_xi_xi * Y_rho * Y_rho * Y_xi + 2.0 * R_xi_rho * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_r_rho( 0, 1, 0 ) = ( R_xi * ( R_rho_rho * Y_xi * Y_xi + ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) * Y_xi - R_rho * Y_rho * Y_xi_xi ) + R_xi * R_xi * ( Y_rho * Y_xi_rho - Y_rho_rho * Y_xi ) - R_rho * R_xi_rho * Y_xi * Y_xi + R_rho * R_xi_xi * Y_rho * Y_xi ) / J3;
    //dFcyl_r_rho( 2, 2, 0 ) = Y_xi / RJ;
    dFcyl_r_rho( 0, 0, 1 ) = dFcyl_r_rho( 0, 1, 0 );
    dFcyl_r_rho( 0, 1, 1 ) = ( Y_xi_xi * R_rho * R_rho * R_xi - 2.0 * Y_xi_rho * R_rho * R_xi * R_xi + Y_rho_rho * R_xi * R_xi * R_xi ) / J3 - ( Y_xi * ( R_xi_xi * R_rho * R_rho - 2.0 * R_xi_rho * R_rho * R_xi + R_rho_rho * R_xi * R_xi ) ) / J3;
    //dFcyl_r_rho( 2, 2, 1 ) = -R_xi / RJ;
    //dFcyl_r_rho( 0, 2, 2 ) = dFcyl_r_rho( 2, 2, 0 );
    //dFcyl_r_rho( 2, 0, 2 ) = dFcyl_r_rho( 2, 2, 0 );
    //dFcyl_r_rho( 2, 1, 2 ) = dFcyl_r_rho( 2, 2, 1 );

    // r_xi
    dFcyl_r_xi( 0, 0, 0 ) = -( R_rho * ( Y_xi_xi * Y_rho * Y_rho - 2.0 * Y_xi_rho * Y_rho * Y_xi + Y_rho_rho * Y_xi * Y_xi ) - R_xi_xi * Y_rho * Y_rho * Y_rho - R_rho_rho * Y_rho * Y_xi * Y_xi + 2.0 * R_xi_rho * Y_rho * Y_rho * Y_xi ) / J3;
    dFcyl_r_xi( 0, 1, 0 ) = ( R_xi * ( Y_xi * ( R_rho * Y_rho_rho - R_rho_rho * Y_rho ) + R_xi_rho * Y_rho * Y_rho - R_rho * Y_rho * Y_xi_rho ) - Y_xi * ( Y_xi_rho * R_rho * R_rho - R_xi_rho * Y_rho * R_rho ) - R_rho * R_xi_xi * Y_rho * Y_rho + R_rho * R_rho * Y_rho * Y_xi_xi ) / J3;
    //dFcyl_r_xi( 2, 2, 0 ) = -Y_rho / RJ;
    dFcyl_r_xi( 0, 0, 1 ) = dFcyl_r_xi( 0, 1, 0 );
    dFcyl_r_xi( 0, 1, 1 ) = -( Y_xi_xi * R_rho * R_rho * R_rho - 2.0 * Y_xi_rho * R_rho * R_rho * R_xi - R_xi_xi * Y_rho * R_rho * R_rho + Y_rho_rho * R_rho * R_xi * R_xi + 2.0 * R_xi_rho * Y_rho * R_rho * R_xi - R_rho_rho * Y_rho * R_xi * R_xi ) / J3;
    //dFcyl_r_xi( 2, 2, 1 ) = R_rho / RJ;
    //dFcyl_r_xi( 0, 2, 2 ) = dFcyl_r_xi( 2, 2, 0 );
    //dFcyl_r_xi( 2, 0, 2 ) = dFcyl_r_xi( 2, 2, 0 );
    //dFcyl_r_xi( 2, 1, 2 ) = dFcyl_r_xi( 2, 2, 1 );

    // r
    //dFcyl_r( 2, 2, 0 ) = -1.0 / ( R * R );
    //dFcyl_r( 0, 2, 2 ) = dFcyl_r( 2, 2, 0 );
    //dFcyl_r( 2, 0, 2 ) = dFcyl_r( 2, 2, 0 );

    // y_rho_rho
    dFcyl_y_rho_rho( 1, 0, 0 ) = ( R_rho * Y_xi * Y_xi * Y_xi - R_xi * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_y_rho_rho( 1, 1, 0 ) = -( -Y_rho * R_xi * R_xi * Y_xi + R_rho * R_xi * Y_xi * Y_xi ) / J3;
    dFcyl_y_rho_rho( 1, 0, 1 ) = dFcyl_y_rho_rho( 1, 1, 0 );
    dFcyl_y_rho_rho( 1, 1, 1 ) = ( R_rho * R_xi * R_xi * Y_xi - R_xi * R_xi * R_xi * Y_rho ) / J3;

    // y_xi_rho
    dFcyl_y_xi_rho( 1, 0, 0 ) = -( -2.0 * R_xi * Y_rho * Y_rho * Y_xi + 2.0 * R_rho * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_y_xi_rho( 1, 1, 0 ) = ( R_rho * R_rho * Y_xi * Y_xi - R_xi * R_xi * Y_rho * Y_rho ) / J3;
    dFcyl_y_xi_rho( 1, 0, 1 ) = dFcyl_y_xi_rho( 1, 1, 0 );
    dFcyl_y_xi_rho( 1, 1, 1 ) = ( 2.0 * R_rho * R_xi * R_xi * Y_rho - 2.0 * R_rho * R_rho * R_xi * Y_xi ) / J3;

    // y_xi_xi
    dFcyl_y_xi_xi( 1, 0, 0 ) = -( R_xi * Y_rho * Y_rho * Y_rho - R_rho * Y_xi * Y_rho * Y_rho ) / J3;
    dFcyl_y_xi_xi( 1, 1, 0 ) = ( -Y_xi * R_rho * R_rho * Y_rho + R_xi * R_rho * Y_rho * Y_rho ) / J3;
    dFcyl_y_xi_xi( 1, 0, 1 ) = dFcyl_y_xi_xi( 1, 1, 0 );
    dFcyl_y_xi_xi( 1, 1, 1 ) = ( R_rho * R_rho * R_rho * Y_xi - R_rho * R_rho * R_xi * Y_rho ) / J3;

    // y_rho
    dFcyl_y_rho( 1, 0, 0 ) = ( R_xi * ( Y_xi_xi * Y_rho * Y_rho - 2.0 * Y_xi_rho * Y_rho * Y_xi + Y_rho_rho * Y_xi * Y_xi ) - R_rho_rho * Y_xi * Y_xi * Y_xi - R_xi_xi * Y_rho * Y_rho * Y_xi + 2.0 * R_xi_rho * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_y_rho( 1, 1, 0 ) = ( R_xi * ( R_rho_rho * Y_xi * Y_xi + ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) * Y_xi - R_rho * Y_rho * Y_xi_xi ) + R_xi * R_xi * ( Y_rho * Y_xi_rho - Y_rho_rho * Y_xi ) - R_rho * R_xi_rho * Y_xi * Y_xi + R_rho * R_xi_xi * Y_rho * Y_xi ) / J3;
    dFcyl_y_rho( 1, 0, 1 ) = dFcyl_y_rho( 1, 1, 0 );
    dFcyl_y_rho( 1, 1, 1 ) = ( Y_xi_xi * R_rho * R_rho * R_xi - 2.0 * Y_xi_rho * R_rho * R_xi * R_xi + Y_rho_rho * R_xi * R_xi * R_xi ) / J3 - ( Y_xi * ( R_xi_xi * R_rho * R_rho - 2.0 * R_xi_rho * R_rho * R_xi + R_rho_rho * R_xi * R_xi ) ) / J3;
    //dFcyl_y_rho( 1, 2, 2 ) = Y_xi / RJ;

    // y_xi
    dFcyl_y_xi( 1, 0, 0 ) = -( R_rho * ( Y_xi_xi * Y_rho * Y_rho - 2.0 * Y_xi_rho * Y_rho * Y_xi + Y_rho_rho * Y_xi * Y_xi ) - R_xi_xi * Y_rho * Y_rho * Y_rho - R_rho_rho * Y_rho * Y_xi * Y_xi + 2.0 * R_xi_rho * Y_rho * Y_rho * Y_xi ) / J3;
    dFcyl_y_xi( 1, 1, 0 ) = ( R_xi * ( Y_xi * ( R_rho * Y_rho_rho - R_rho_rho * Y_rho ) + R_xi_rho * Y_rho * Y_rho - R_rho * Y_rho * Y_xi_rho ) - Y_xi * ( Y_xi_rho * R_rho * R_rho - R_xi_rho * Y_rho * R_rho ) - R_rho * R_xi_xi * Y_rho * Y_rho + R_rho * R_rho * Y_rho * Y_xi_xi ) / J3;
    dFcyl_y_xi( 1, 0, 1 ) = dFcyl_y_xi( 1, 1, 0 );
    dFcyl_y_xi( 1, 1, 1 ) = -( Y_xi_xi * R_rho * R_rho * R_rho - 2.0 * Y_xi_rho * R_rho * R_rho * R_xi - R_xi_xi * Y_rho * R_rho * R_rho + Y_rho_rho * R_rho * R_xi * R_xi + 2.0 * R_xi_rho * Y_rho * R_rho * R_xi - R_rho_rho * Y_rho * R_xi * R_xi ) / J3;
    //dFcyl_y_xi( 1, 2, 2 ) = -Y_rho / RJ;
    //////////////////////

    // change to voigt
    FloatArray dFcyl_r_rho_rho_v = dFcyl_r_rho_rho.to_voigt_form_27();
    FloatArray dFcyl_r_xi_xi_v   = dFcyl_r_xi_xi.to_voigt_form_27();
    FloatArray dFcyl_r_xi_rho_v  = dFcyl_r_xi_rho.to_voigt_form_27();
    FloatArray dFcyl_r_rho_v     = dFcyl_r_rho.to_voigt_form_27();
    FloatArray dFcyl_r_xi_v      = dFcyl_r_xi.to_voigt_form_27();
    FloatArray dFcyl_r_v         = dFcyl_r.to_voigt_form_27();

    FloatArray dFcyl_y_rho_rho_v = dFcyl_y_rho_rho.to_voigt_form_27();
    FloatArray dFcyl_y_xi_xi_v   = dFcyl_y_xi_xi.to_voigt_form_27();
    FloatArray dFcyl_y_xi_rho_v  = dFcyl_y_xi_rho.to_voigt_form_27();
    FloatArray dFcyl_y_rho_v     = dFcyl_y_rho.to_voigt_form_27();
    FloatArray dFcyl_y_xi_v      = dFcyl_y_xi.to_voigt_form_27();

    // Transform to global coordinates!!
    int nRows = ders.giveNumberOfRows();
    answer.resize( 27, nRows * 2 );
    answer.zero();

    // assume
    double N, N_rho, N_xi, N_rho_rho, N_xi_xi, N_rho_xi, N_xi_rho;
    FloatArray col1, col2;
    for ( int i = 1; i <= nRows; i++ ) {
        // basis and its derivatives
        //N         = ders( i - 1, 0 );
        N_rho     = ders( i - 1, 1 );
        N_xi      = ders( i - 1, 2 );
        N_rho_rho = ders( i - 1, 3 );
        N_xi_xi   = ders( i - 1, 4 );
        N_rho_xi  = ders( i - 1, 5 );
        N_xi_rho  = N_rho_xi;

        col1 = N_rho_rho * dFcyl_r_rho_rho_v
            + N_xi_xi * dFcyl_r_xi_xi_v
            + N_xi_rho * dFcyl_r_xi_rho_v
            + N_rho * dFcyl_r_rho_v
            + N_xi * dFcyl_r_xi_v;
            //+ N * dFcyl_r_v;

        col2 = N_rho_rho * dFcyl_y_rho_rho_v
            + N_xi_xi * dFcyl_y_xi_xi_v
            + N_xi_rho * dFcyl_y_xi_rho_v
            + N_rho * dFcyl_y_rho_v
            + N_xi * dFcyl_y_xi_v;

        answer.setColumn( col1, 2 * i - 1 );
        answer.setColumn( col2, 2 * i );
    }
}

void PlaneStrainStructuralElementEvaluator::computeGradientOfDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );
    elem->computeVectorOf( VM_Total, tStep, u );

    //// Get nodal coordinates
    FloatArray nodesCoords;
    elem->giveNodalCoordinates( nodesCoords );

    FloatArray xDef = nodesCoords + u;

    // Displacement gradient H = du/dX
    FloatMatrix G;
    this->computeGMatrixAt( gp, G );

    // get local code numbers corresponding to ir
    IntArray lc;
    FloatArray xr;
    this->giveIntegrationElementLocalCodeNumbers( lc, elem, gp->giveIntegrationRule() );
    xr.resize( G.giveNumberOfColumns() );
    for ( int i = 1; i <= lc.giveSize(); i++ ) {
        xr.at( i ) = xDef.at( lc.at( i ) );
    }

    answer.beProductOf( G, xr );
}

//void PlaneStrainStructuralElementEvaluator::computeSecondOrderStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
//{
//    FloatArray vF, vG;
//    this->computeDeformationGradientVector( vF, gp, tStep );
//    this->computeGradientOfDeformationGradientVector( vG, gp, tStep );
//    answer = static_cast<StructuralCrossSection *>( this->giveElement()->giveCrossSection() )->giveSecondOrderStresses( vF, vG, gp, tStep );
//}
//
//
//void PlaneStrainStructuralElementEvaluator ::computeConstitutiveMatrix_dAddF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
//{
//    answer = static_cast<StructuralCrossSection *>( this->giveElement()->giveCrossSection() )->giveStiffnessMatrix_dAddF_3d( rMode, gp, tStep );
//}
//
//void PlaneStrainStructuralElementEvaluator ::computeConstitutiveMatrix_dAdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
//{
//    answer = static_cast<StructuralCrossSection *>( this->giveElement()->giveCrossSection() )->giveStiffnessMatrix_dAdF_3d( rMode, gp, tStep );
//}


///////
////Axisymetric element
void AxisymStructuralElementEvaluator ::computeNMatrixAt( FloatMatrix &answer, GaussPoint *gp )
{
    FloatArray N;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    interp->evalN( N, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    answer.beNMatrixOf( N, 2 );
}

void AxisymStructuralElementEvaluator::computeBHmatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    // Returns the [ 4 x (nno*2) ] displacement gradient matrix {BH} of the receiver,
    // evaluated at gp.
    /// @todo not checked if correct
    FloatArray n;
    FloatMatrix dNdx;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    // this uses FEInterpolation::nodes2coords - quite inefficient in this case (large num of dofmans)
    interp->evaldNdx( dNdx, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    interp->evalN( n, gp->giveNaturalCoordinates(), 
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );


    FloatArray gcoords;
    interp->local2global( gcoords, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    double r = gcoords.at( 1 ); 

    
    answer.resize( 9, dNdx.giveNumberOfRows() * 2 );
    answer.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at( 1, 2 * i - 1 ) = dNdx.at( i, 1 ); // du/dx -1
        answer.at( 2, 2 * i - 0 ) = dNdx.at( i, 2 ); // dv/dy -2
        answer.at( 6, 2 * i - 1 ) = dNdx.at( i, 2 ); // du/dy -6
        answer.at( 9, 2 * i - 0 ) = dNdx.at( i, 1 ); // dv/dx -9
    }

    for ( int i = 0; i < dNdx.giveNumberOfRows(); i++ ) {
        answer.at( 3, 2 * i + 1 ) = n.at( i + 1 ) / r;
    }

}


double AxisymStructuralElementEvaluator ::computeVolumeAround( GaussPoint *gp )
{
    double determinant, weight, thickness, volume;
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();
    determinant = fabs( interp->giveTransformationJacobian( gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( this->giveElement(),gp->giveIntegrationRule()->giveKnotSpan() ) ) );
    weight      = gp->giveWeight();
    thickness   = this->giveElement()->giveCrossSection()->give( CS_Thickness, gp );

    FloatArray gcoords;
    interp->local2global( gcoords, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    double r = gcoords.at( 1 ); 

    volume      = determinant * weight * thickness*r;

    return volume;
}



void AxisymStructuralElementEvaluator ::initializeFrom( InputRecord &ir )
{
    StructuralElementEvaluator::initializeFrom( ir );
}

void AxisymStructuralElementEvaluator::computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
{
    answer = static_cast<StructuralCrossSection *>( this->giveElement()->giveCrossSection() )->giveStiffnessMatrix_dPdF_3d( rMode, gp, tStep );
}


////////////////
// second order stuff

void AxisymStructuralElementEvaluator::computeGMatrixAt( GaussPoint *gp, FloatMatrix &answer )
{
    FloatMatrix ders; 
    FEInterpolation *interp = gp->giveElement()->giveInterpolation();

    // returns index where nonzero basis functions start (1-based)
    // ders is a matrix of size [nno x 6] where nno is the number of nonzero basis functions at the gauss point, the columns are: [N, dN/drho, dN/dxi, d2N/drho2, d2N/dxi2, d2N/drhodxi]
    auto knotSpan  = gp->giveIntegrationRule()->giveKnotSpan();
    int indexStart = interp->evalSecondDerivatives( ders, gp->giveNaturalCoordinates(),
        FEIIGAElementGeometryWrapper( gp->giveElement(), knotSpan ) );

    // Get nodal coordinates
    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );
    FloatArray nodesCoords;
    elem->giveNodalCoordinates( nodesCoords );

    double R = 0.; // undeformed radius (X1)
    double R_rho = 0., R_xi = 0., Y_rho = 0., Y_xi = 0.; // 1st derivatives
    double R_rho_rho = 0., R_xi_xi = 0., R_rho_xi = 0., Y_rho_rho = 0., Y_xi_xi = 0., Y_rho_xi = 0.; // 2nd derivatives

    IntArray mask;
    interp->giveKnotSpanBasisFuncMask( *knotSpan, mask );
    
    // Evaluate reference coordinate derivatives
    for ( int k = 0; k < ders.giveNumberOfRows(); k++ ) {
        //int indCords         = indexStart + k;
        int indCords         = mask.at(k+1);
        double vertexCoord_1 = nodesCoords.at( 2 * indCords - 1 );
        double vertexCoord_2 = nodesCoords.at( 2 * indCords );
        // position
        R += vertexCoord_1 * ders( k, 0 ); // compute radial coordinate at the gauss point
        // first derivatives
        R_rho += ders( k, 1 ) * vertexCoord_1;
        Y_rho += ders( k, 1 ) * vertexCoord_2;

        R_xi += ders( k, 2) * vertexCoord_1;
        Y_xi += ders( k, 2) * vertexCoord_2;
        // second derivatives
        R_rho_rho += ders( k, 3 ) * vertexCoord_1;
        Y_rho_rho += ders( k, 3 ) * vertexCoord_2;

        R_xi_xi += ders( k, 4 ) * vertexCoord_1;
        Y_xi_xi += ders( k, 4 ) * vertexCoord_2;

        R_rho_xi += ders( k, 5 ) * vertexCoord_1;
        Y_rho_xi += ders( k, 5 ) * vertexCoord_2; 
    }
    double R_xi_rho = R_rho_xi;
    double Y_xi_rho = Y_rho_xi;
    
    ////////////////////////////////////////////////////////////
    //// sqrt of determinant of the reference metric tensor
    //double Gb = R * ( R_rho * Y_xi - R_xi * Y_rho ); // sqrt of determinant of the reference metric tensor
    //double Gb2 = Gb * Gb;
    //double R2  = R * R;
    //double R3  = R * R * R;

    //// Components greadF^i_jk in basis {G_i} \otimes {G^j} \otimes {G^k}
    //Tensor3_3d dF_r_rho_rho;
    //Tensor3_3d dF_r_xi_rho;
    //Tensor3_3d dF_r_xi_xi;
    //Tensor3_3d dF_r_rho;
    //Tensor3_3d dF_r_xi;
    //Tensor3_3d dF_r;
    //Tensor3_3d dF_y_rho_rho;
    //Tensor3_3d dF_y_xi_rho;
    //Tensor3_3d dF_y_xi_xi;
    //Tensor3_3d dF_y_rho;
    //Tensor3_3d dF_y_xi;

    //// Tensor 1: dF_r_rho_rho
    //dF_r_rho_rho( 0, 0, 0 ) = ( R * Y_xi ) / Gb;
    //dF_r_rho_rho( 1, 0, 0 ) = -( R * Y_rho ) / Gb;

    //// Tensor 3: dF_r_xi_xi
    //dF_r_xi_xi( 0, 1, 1 ) = ( R * Y_xi ) / Gb;
    //dF_r_xi_xi( 1, 1, 1 ) = -( R * Y_rho ) / Gb;

    //// Tensor 2: dF_r_xi_rho
    //dF_r_xi_rho( 0, 1, 0 ) = ( R * Y_xi ) / Gb;
    //dF_r_xi_rho( 1, 1, 0 ) = -( R * Y_rho ) / Gb;
    //dF_r_xi_rho( 0, 0, 1 ) = ( R * Y_xi ) / Gb;
    //dF_r_xi_rho( 1, 0, 1 ) = -( R * Y_rho ) / Gb;

    //// Tensor 4: dF_r_rho
    //dF_r_rho( 0, 0, 0 ) = -( R2 * Y_xi * ( R_rho_rho * Y_xi - R_xi * Y_rho_rho ) ) / Gb2;
    //dF_r_rho( 0, 1, 0 ) = ( R2 * Y_xi * ( R_xi * Y_xi_rho - R_xi_rho * Y_xi ) ) / Gb2;
    //dF_r_rho( 1, 0, 0 ) = ( R2 * Y_rho * ( R_rho_rho * Y_xi - R_xi * Y_rho_rho ) ) / Gb2;
    //dF_r_rho( 1, 1, 0 ) = -( R2 * Y_rho * ( R_xi * Y_xi_rho - R_xi_rho * Y_xi ) ) / Gb2;
    //dF_r_rho( 2, 2, 0 ) = 1.0 / R;

    //dF_r_rho( 0, 0, 1 ) = ( R2 * Y_xi * ( R_xi * Y_xi_rho - R_xi_rho * Y_xi ) ) / Gb2;
    //dF_r_rho( 0, 1, 1 ) = ( R2 * Y_xi * ( R_xi * Y_xi_xi - R_xi_xi * Y_xi ) ) / Gb2;
    //dF_r_rho( 1, 0, 1 ) = -( R2 * Y_rho * ( R_xi * Y_xi_rho - R_xi_rho * Y_xi ) ) / Gb2;
    //dF_r_rho( 1, 1, 1 ) = -( R2 * Y_rho * ( R_xi * Y_xi_xi - R_xi_xi * Y_xi ) ) / Gb2;

    //dF_r_rho( 0, 2, 2 ) = ( R3 * Y_xi * Y_xi ) / Gb2;
    //dF_r_rho( 1, 2, 2 ) = -( R3 * Y_rho * Y_xi ) / Gb2;
    //dF_r_rho( 2, 0, 2 ) = 1.0 / R;

    //// Tensor 5: dF_r_xi
    //dF_r_xi( 0, 0, 0 ) = -( R2 * Y_xi * ( R_rho * Y_rho_rho - R_rho_rho * Y_rho ) ) / Gb2;
    //dF_r_xi( 0, 1, 0 ) = -( R2 * Y_xi * ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) ) / Gb2;
    //dF_r_xi( 1, 0, 0 ) = ( R2 * Y_rho * ( R_rho * Y_rho_rho - R_rho_rho * Y_rho ) ) / Gb2;
    //dF_r_xi( 1, 1, 0 ) = ( R2 * Y_rho * ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) ) / Gb2;

    //dF_r_xi( 0, 0, 1 ) = -( R2 * Y_xi * ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) ) / Gb2;
    //dF_r_xi( 0, 1, 1 ) = -( R2 * Y_xi * ( R_rho * Y_xi_xi - R_xi_xi * Y_rho ) ) / Gb2;
    //dF_r_xi( 1, 0, 1 ) = ( R2 * Y_rho * ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) ) / Gb2;
    //dF_r_xi( 1, 1, 1 ) = ( R2 * Y_rho * ( R_rho * Y_xi_xi - R_xi_xi * Y_rho ) ) / Gb2;
    //dF_r_xi( 2, 2, 1 ) = 1.0 / R;

    //dF_r_xi( 0, 2, 2 ) = -( R3 * Y_rho * Y_xi ) / Gb2;
    //dF_r_xi( 1, 2, 2 ) = ( R3 * Y_rho * Y_rho ) / Gb2;
    //dF_r_xi( 2, 1, 2 ) = 1.0 / R;

    //// Tensor 6: dF_r
    //dF_r( 2, 2, 0 ) = -R_rho / R2;
    //dF_r( 2, 2, 1 ) = -R_xi / R2;
    //dF_r( 0, 2, 2 ) = -( R * Y_xi ) / Gb;
    //dF_r( 1, 2, 2 ) = -( R2 * ( R_xi * Y_rho * Y_rho - R_rho * Y_xi * Y_rho ) ) / Gb2;
    //dF_r( 2, 0, 2 ) = -R_rho / R2;
    //dF_r( 2, 1, 2 ) = -R_xi / R2;

    //// Tensor 7: dF_y_rho_rho
    //dF_y_rho_rho( 0, 0, 0 ) = -( R * R_xi ) / Gb;
    //dF_y_rho_rho( 1, 0, 0 ) = ( R * R_rho ) / Gb;

    //// Tensor 9: dF_y_xi_xi
    //dF_y_xi_xi( 0, 1, 1 ) = -( R * R_xi ) / Gb;
    //dF_y_xi_xi( 1, 1, 1 ) = ( R * R_rho ) / Gb;

    //// Tensor 8: dF_y_xi_rho
    //dF_y_xi_rho( 0, 1, 0 ) = -( R * R_xi ) / Gb;
    //dF_y_xi_rho( 2, 2, 1 ) = ( R * R_rho ) / Gb;
    //dF_y_xi_rho( 0, 0, 1 ) = -( R * R_xi ) / Gb;
    //dF_y_xi_rho( 1, 0, 1 ) = ( R * R_rho ) / Gb;

    //// Tensor 10: dF_y_rho
    //dF_y_rho( 0, 0, 0 ) = ( R2 * R_xi * ( R_rho_rho * Y_xi - R_xi * Y_rho_rho ) ) / Gb2;
    //dF_y_rho( 0, 1, 0 ) = -( R2 * R_xi * ( R_xi * Y_xi_rho - R_xi_rho * Y_xi ) ) / Gb2;
    //dF_y_rho( 1, 0, 0 ) = -( R2 * R_rho * ( R_rho_rho * Y_xi - R_xi * Y_rho_rho ) ) / Gb2;
    //dF_y_rho( 1, 1, 0 ) = ( R2 * R_rho * ( R_xi * Y_xi_rho - R_xi_rho * Y_xi ) ) / Gb2;

    //dF_y_rho( 0, 0, 1 ) = -( R2 * R_xi * ( R_xi * Y_xi_rho - R_xi_rho * Y_xi ) ) / Gb2;
    //dF_y_rho( 0, 1, 1 ) = -( R2 * R_xi * ( R_xi * Y_xi_xi - R_xi_xi * Y_xi ) ) / Gb2;
    //dF_y_rho( 1, 0, 1 ) = ( R2 * R_rho * ( R_xi * Y_xi_rho - R_xi_rho * Y_xi ) ) / Gb2;
    //dF_y_rho( 1, 1, 1 ) = ( R2 * R_rho * ( R_xi * Y_xi_xi - R_xi_xi * Y_xi ) ) / Gb2;

    //dF_y_rho( 0, 2, 2 ) = -( R3 * R_xi * Y_xi ) / Gb2;
    //dF_y_rho( 1, 2, 2 ) = ( R3 * R_rho * Y_xi ) / Gb2;

    //// Tensor 11: dF_y_xi
    //dF_y_xi( 0, 0, 0 ) = ( R2 * R_xi * ( R_rho * Y_rho_rho - R_rho_rho * Y_rho ) ) / Gb2;
    //dF_y_xi( 0, 1, 0 ) = ( R2 * R_xi * ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) ) / Gb2;
    //dF_y_xi( 1, 0, 0 ) = -( R2 * R_rho * ( R_rho * Y_rho_rho - R_rho_rho * Y_rho ) ) / Gb2;
    //dF_y_xi( 1, 1, 0 ) = -( R2 * R_rho * ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) ) / Gb2;

    //dF_y_xi( 0, 0, 1 ) = ( R2 * R_xi * ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) ) / Gb2;
    //dF_y_xi( 0, 1, 1 ) = ( R2 * R_xi * ( R_rho * Y_xi_xi - R_xi_xi * Y_rho ) ) / Gb2;
    //dF_y_xi( 1, 0, 1 ) = -( R2 * R_rho * ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) ) / Gb2;
    //dF_y_xi( 1, 1, 1 ) = -( R2 * R_rho * ( R_rho * Y_xi_xi - R_xi_xi * Y_rho ) ) / Gb2;

    //dF_y_xi( 0, 2, 2 ) = ( R3 * R_xi * Y_rho ) / Gb2;
    //dF_y_xi( 1, 2, 2 ) = -( R3 * R_rho * Y_rho ) / Gb2;


    //// Reference base vectors in terms of cartesian basis ordered as {G_rho, G_xi, G_theta} and their duals {G^rho, G^xi, G^theta}
    //FloatMatrixF<3, 3> G_covMF( R_rho, Y_rho, 0., R_xi, Y_xi, 0., 0., 0., R ); // Collumn wise
    //FloatMatrixF<3, 3> G_contraMF( R * Y_xi / Gb,  -R*R_xi / Gb, 0.,  -R * Y_rho / Gb,  R*R_rho / Gb,  0., 0., 0., 1. / R ); // 3rd col.

    //FloatMatrix G_covM( G_covMF ), G_contraM( G_contraMF );
    //Tensor2_3d G_cov( G_covM ), G_contra( G_contraM ); // create tensors
    ////////////////////////////////////////////////////////////

    // Create part proportional to individual coordinates in cylindrical system
    Tensor3_3d dFcyl_r_rho_rho;
    Tensor3_3d dFcyl_r_xi_rho;
    Tensor3_3d dFcyl_r_xi_xi;
    Tensor3_3d dFcyl_r_rho;
    Tensor3_3d dFcyl_r_xi;
    Tensor3_3d dFcyl_r;
    Tensor3_3d dFcyl_y_rho_rho;
    Tensor3_3d dFcyl_y_xi_rho;
    Tensor3_3d dFcyl_y_xi_xi;
    Tensor3_3d dFcyl_y_rho;
    Tensor3_3d dFcyl_y_xi;

    ////////////////////////////////////////////////////////////
    //dFcyl_r_rho_rho( i_3, j_3, k_3 ) = dF_r_rho_rho( a_3, b_3, c_3 ) * G_cov( i_3, a_3 ) * G_contra( j_3, b_3 ) * G_contra( k_3, c_3 ); 
    //dFcyl_r_xi_xi( i_3, j_3, k_3 )   = dF_r_xi_xi( a_3, b_3, c_3 ) * G_cov( i_3, a_3 ) * G_contra( j_3, b_3 ) * G_contra( k_3, c_3 );
    //dFcyl_r_xi_rho( i_3, j_3, k_3 )  = dF_r_xi_rho( a_3, b_3, c_3 ) * G_cov( i_3, a_3 ) * G_contra( j_3, b_3 ) * G_contra( k_3, c_3 );
    //dFcyl_r_rho( i_3, j_3, k_3 )     = dF_r_rho( a_3, b_3, c_3 ) * G_cov( i_3, a_3 ) * G_contra( j_3, b_3 ) * G_contra( k_3, c_3 );
    //dFcyl_r_xi( i_3, j_3, k_3 )      = dF_r_xi( a_3, b_3, c_3 ) * G_cov( i_3, a_3 ) * G_contra( j_3, b_3 ) * G_contra( k_3, c_3 );
    //dFcyl_r( i_3, j_3, k_3 )         = dF_r( a_3, b_3, c_3 ) * G_cov( i_3, a_3 ) * G_contra( j_3, b_3 ) * G_contra( k_3, c_3 );

    //dFcyl_y_rho_rho( i_3, j_3, k_3 ) = dF_y_rho_rho( a_3, b_3, c_3 ) * G_cov( i_3, a_3 ) * G_contra( j_3, b_3 ) * G_contra( k_3, c_3 );
    //dFcyl_y_xi_xi( i_3, j_3, k_3 )   = dF_y_xi_xi( a_3, b_3, c_3 ) * G_cov( i_3, a_3 ) * G_contra( j_3, b_3 ) * G_contra( k_3, c_3 );
    //dFcyl_y_xi_rho( i_3, j_3, k_3 )  = dF_y_xi_rho( a_3, b_3, c_3 ) * G_cov( i_3, a_3 ) * G_contra( j_3, b_3 ) * G_contra( k_3, c_3 );
    //dFcyl_y_rho( i_3, j_3, k_3 )     = dF_y_rho( a_3, b_3, c_3 ) * G_cov( i_3, a_3 ) * G_contra( j_3, b_3 ) * G_contra( k_3, c_3 );
    //dFcyl_y_xi( i_3, j_3, k_3 )      = dF_y_xi( a_3, b_3, c_3 ) * G_cov( i_3, a_3 ) * G_contra( j_3, b_3 ) * G_contra( k_3, c_3 );
    ////////////////////////////////////////////////////////////

    ////////////////////
    const double J  = R_rho * Y_xi - R_xi * Y_rho;
    const double J3 = J * J * J;
    const double RJ = R * J;

    // r_rho_rho
    dFcyl_r_rho_rho( 0, 0, 0 ) = ( R_rho * Y_xi * Y_xi * Y_xi - R_xi * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_r_rho_rho( 0, 1, 0 ) = -( -Y_rho * R_xi * R_xi * Y_xi + R_rho * R_xi * Y_xi * Y_xi ) / J3;
    dFcyl_r_rho_rho( 0, 0, 1 ) = dFcyl_r_rho_rho( 0, 1, 0 );
    dFcyl_r_rho_rho( 0, 1, 1 ) = ( R_rho * R_xi * R_xi * Y_xi - R_xi * R_xi * R_xi * Y_rho ) / J3;

    // r_xi_rho
    dFcyl_r_xi_rho( 0, 0, 0 ) = -( -2.0 * R_xi * Y_rho * Y_rho * Y_xi + 2.0 * R_rho * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_r_xi_rho( 0, 1, 0 ) = ( R_rho * R_rho * Y_xi * Y_xi - R_xi * R_xi * Y_rho * Y_rho ) / J3;
    dFcyl_r_xi_rho( 0, 0, 1 ) = dFcyl_r_xi_rho( 0, 1, 0 );
    dFcyl_r_xi_rho( 0, 1, 1 ) = ( 2.0 * R_rho * R_xi * R_xi * Y_rho - 2.0 * R_rho * R_rho * R_xi * Y_xi ) / J3;

    // r_xi_xi
    dFcyl_r_xi_xi( 0, 0, 0 ) = -( R_xi * Y_rho * Y_rho * Y_rho - R_rho * Y_xi * Y_rho * Y_rho ) / J3;
    dFcyl_r_xi_xi( 0, 1, 0 ) = ( -Y_xi * R_rho * R_rho * Y_rho + R_xi * R_rho * Y_rho * Y_rho ) / J3;
    dFcyl_r_xi_xi( 0, 0, 1 ) = dFcyl_r_xi_xi( 0, 1, 0 );
    dFcyl_r_xi_xi( 0, 1, 1 ) = ( R_rho * R_rho * R_rho * Y_xi - R_rho * R_rho * R_xi * Y_rho ) / J3;

    // r_rho
    dFcyl_r_rho( 0, 0, 0 ) = ( R_xi * ( Y_xi_xi * Y_rho * Y_rho - 2.0 * Y_xi_rho * Y_rho * Y_xi + Y_rho_rho * Y_xi * Y_xi ) - R_rho_rho * Y_xi * Y_xi * Y_xi - R_xi_xi * Y_rho * Y_rho * Y_xi + 2.0 * R_xi_rho * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_r_rho( 0, 1, 0 ) = ( R_xi * ( R_rho_rho * Y_xi * Y_xi + ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) * Y_xi - R_rho * Y_rho * Y_xi_xi ) + R_xi * R_xi * ( Y_rho * Y_xi_rho - Y_rho_rho * Y_xi ) - R_rho * R_xi_rho * Y_xi * Y_xi + R_rho * R_xi_xi * Y_rho * Y_xi ) / J3;
    dFcyl_r_rho( 2, 2, 0 ) = Y_xi / RJ;
    dFcyl_r_rho( 0, 0, 1 ) = dFcyl_r_rho( 0, 1, 0 );
    dFcyl_r_rho( 0, 1, 1 ) = ( Y_xi_xi * R_rho * R_rho * R_xi - 2.0 * Y_xi_rho * R_rho * R_xi * R_xi + Y_rho_rho * R_xi * R_xi * R_xi ) / J3 - ( Y_xi * ( R_xi_xi * R_rho * R_rho - 2.0 * R_xi_rho * R_rho * R_xi + R_rho_rho * R_xi * R_xi ) ) / J3;
    dFcyl_r_rho( 2, 2, 1 ) = -R_xi / RJ;
    dFcyl_r_rho( 0, 2, 2 ) = dFcyl_r_rho( 2, 2, 0 );
    dFcyl_r_rho( 2, 0, 2 ) = dFcyl_r_rho( 2, 2, 0 );
    dFcyl_r_rho( 2, 1, 2 ) = dFcyl_r_rho( 2, 2, 1 );

    // r_xi
    dFcyl_r_xi( 0, 0, 0 ) = -( R_rho * ( Y_xi_xi * Y_rho * Y_rho - 2.0 * Y_xi_rho * Y_rho * Y_xi + Y_rho_rho * Y_xi * Y_xi ) - R_xi_xi * Y_rho * Y_rho * Y_rho - R_rho_rho * Y_rho * Y_xi * Y_xi + 2.0 * R_xi_rho * Y_rho * Y_rho * Y_xi ) / J3;
    dFcyl_r_xi( 0, 1, 0 ) = ( R_xi * ( Y_xi * ( R_rho * Y_rho_rho - R_rho_rho * Y_rho ) + R_xi_rho * Y_rho * Y_rho - R_rho * Y_rho * Y_xi_rho ) - Y_xi * ( Y_xi_rho * R_rho * R_rho - R_xi_rho * Y_rho * R_rho ) - R_rho * R_xi_xi * Y_rho * Y_rho + R_rho * R_rho * Y_rho * Y_xi_xi ) / J3;
    dFcyl_r_xi( 2, 2, 0 ) = -Y_rho / RJ;
    dFcyl_r_xi( 0, 0, 1 ) = dFcyl_r_xi( 0, 1, 0 );
    dFcyl_r_xi( 0, 1, 1 ) = -( Y_xi_xi * R_rho * R_rho * R_rho - 2.0 * Y_xi_rho * R_rho * R_rho * R_xi - R_xi_xi * Y_rho * R_rho * R_rho + Y_rho_rho * R_rho * R_xi * R_xi + 2.0 * R_xi_rho * Y_rho * R_rho * R_xi - R_rho_rho * Y_rho * R_xi * R_xi ) / J3;
    dFcyl_r_xi( 2, 2, 1 ) = R_rho / RJ;
    dFcyl_r_xi( 0, 2, 2 ) = dFcyl_r_xi( 2, 2, 0 );
    dFcyl_r_xi( 2, 0, 2 ) = dFcyl_r_xi( 2, 2, 0 );
    dFcyl_r_xi( 2, 1, 2 ) = dFcyl_r_xi( 2, 2, 1 );

    // r
    dFcyl_r( 2, 2, 0 ) = -1.0 / ( R * R );
    dFcyl_r( 0, 2, 2 ) = dFcyl_r( 2, 2, 0 );
    dFcyl_r( 2, 0, 2 ) = dFcyl_r( 2, 2, 0 );

    // y_rho_rho
    dFcyl_y_rho_rho( 1, 0, 0 ) = ( R_rho * Y_xi * Y_xi * Y_xi - R_xi * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_y_rho_rho( 1, 1, 0 ) = -( -Y_rho * R_xi * R_xi * Y_xi + R_rho * R_xi * Y_xi * Y_xi ) / J3;
    dFcyl_y_rho_rho( 1, 0, 1 ) = dFcyl_y_rho_rho( 1, 1, 0 );
    dFcyl_y_rho_rho( 1, 1, 1 ) = ( R_rho * R_xi * R_xi * Y_xi - R_xi * R_xi * R_xi * Y_rho ) / J3;

    // y_xi_rho
    dFcyl_y_xi_rho( 1, 0, 0 ) = -( -2.0 * R_xi * Y_rho * Y_rho * Y_xi + 2.0 * R_rho * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_y_xi_rho( 1, 1, 0 ) = ( R_rho * R_rho * Y_xi * Y_xi - R_xi * R_xi * Y_rho * Y_rho ) / J3;
    dFcyl_y_xi_rho( 1, 0, 1 ) = dFcyl_y_xi_rho( 1, 1, 0 );
    dFcyl_y_xi_rho( 1, 1, 1 ) = ( 2.0 * R_rho * R_xi * R_xi * Y_rho - 2.0 * R_rho * R_rho * R_xi * Y_xi ) / J3;

    // y_xi_xi
    dFcyl_y_xi_xi( 1, 0, 0 ) = -( R_xi * Y_rho * Y_rho * Y_rho - R_rho * Y_xi * Y_rho * Y_rho ) / J3;
    dFcyl_y_xi_xi( 1, 1, 0 ) = ( -Y_xi * R_rho * R_rho * Y_rho + R_xi * R_rho * Y_rho * Y_rho ) / J3;
    dFcyl_y_xi_xi( 1, 0, 1 ) = dFcyl_y_xi_xi( 1, 1, 0 );
    dFcyl_y_xi_xi( 1, 1, 1 ) = ( R_rho * R_rho * R_rho * Y_xi - R_rho * R_rho * R_xi * Y_rho ) / J3;

    // y_rho
    dFcyl_y_rho( 1, 0, 0 ) = ( R_xi * ( Y_xi_xi * Y_rho * Y_rho - 2.0 * Y_xi_rho * Y_rho * Y_xi + Y_rho_rho * Y_xi * Y_xi ) - R_rho_rho * Y_xi * Y_xi * Y_xi - R_xi_xi * Y_rho * Y_rho * Y_xi + 2.0 * R_xi_rho * Y_rho * Y_xi * Y_xi ) / J3;
    dFcyl_y_rho( 1, 1, 0 ) = ( R_xi * ( R_rho_rho * Y_xi * Y_xi + ( R_rho * Y_xi_rho - R_xi_rho * Y_rho ) * Y_xi - R_rho * Y_rho * Y_xi_xi ) + R_xi * R_xi * ( Y_rho * Y_xi_rho - Y_rho_rho * Y_xi ) - R_rho * R_xi_rho * Y_xi * Y_xi + R_rho * R_xi_xi * Y_rho * Y_xi ) / J3;
    dFcyl_y_rho( 1, 0, 1 ) = dFcyl_y_rho( 1, 1, 0 );
    dFcyl_y_rho( 1, 1, 1 ) = ( Y_xi_xi * R_rho * R_rho * R_xi - 2.0 * Y_xi_rho * R_rho * R_xi * R_xi + Y_rho_rho * R_xi * R_xi * R_xi ) / J3 - ( Y_xi * ( R_xi_xi * R_rho * R_rho - 2.0 * R_xi_rho * R_rho * R_xi + R_rho_rho * R_xi * R_xi ) ) / J3;
    dFcyl_y_rho( 1, 2, 2 ) = Y_xi / RJ;

    // y_xi
    dFcyl_y_xi( 1, 0, 0 ) = -( R_rho * ( Y_xi_xi * Y_rho * Y_rho - 2.0 * Y_xi_rho * Y_rho * Y_xi + Y_rho_rho * Y_xi * Y_xi ) - R_xi_xi * Y_rho * Y_rho * Y_rho - R_rho_rho * Y_rho * Y_xi * Y_xi + 2.0 * R_xi_rho * Y_rho * Y_rho * Y_xi ) / J3;
    dFcyl_y_xi( 1, 1, 0 ) = ( R_xi * ( Y_xi * ( R_rho * Y_rho_rho - R_rho_rho * Y_rho ) + R_xi_rho * Y_rho * Y_rho - R_rho * Y_rho * Y_xi_rho ) - Y_xi * ( Y_xi_rho * R_rho * R_rho - R_xi_rho * Y_rho * R_rho ) - R_rho * R_xi_xi * Y_rho * Y_rho + R_rho * R_rho * Y_rho * Y_xi_xi ) / J3;
    dFcyl_y_xi( 1, 0, 1 ) = dFcyl_y_xi( 1, 1, 0 );
    dFcyl_y_xi( 1, 1, 1 ) = -( Y_xi_xi * R_rho * R_rho * R_rho - 2.0 * Y_xi_rho * R_rho * R_rho * R_xi - R_xi_xi * Y_rho * R_rho * R_rho + Y_rho_rho * R_rho * R_xi * R_xi + 2.0 * R_xi_rho * Y_rho * R_rho * R_xi - R_rho_rho * Y_rho * R_xi * R_xi ) / J3;
    dFcyl_y_xi( 1, 2, 2 ) = -Y_rho / RJ;
    //////////////////////

    // change to voigt
    FloatArray dFcyl_r_rho_rho_v = dFcyl_r_rho_rho.to_voigt_form_27();
    FloatArray dFcyl_r_xi_xi_v   = dFcyl_r_xi_xi.to_voigt_form_27();
    FloatArray dFcyl_r_xi_rho_v  = dFcyl_r_xi_rho.to_voigt_form_27();
    FloatArray dFcyl_r_rho_v     = dFcyl_r_rho.to_voigt_form_27();
    FloatArray dFcyl_r_xi_v      = dFcyl_r_xi.to_voigt_form_27();
    FloatArray dFcyl_r_v         = dFcyl_r.to_voigt_form_27();

    FloatArray dFcyl_y_rho_rho_v = dFcyl_y_rho_rho.to_voigt_form_27();
    FloatArray dFcyl_y_xi_xi_v   = dFcyl_y_xi_xi.to_voigt_form_27();
    FloatArray dFcyl_y_xi_rho_v  = dFcyl_y_xi_rho.to_voigt_form_27();
    FloatArray dFcyl_y_rho_v     = dFcyl_y_rho.to_voigt_form_27();
    FloatArray dFcyl_y_xi_v      = dFcyl_y_xi.to_voigt_form_27();

    // Transform to global coordinates!!
    int nRows = ders.giveNumberOfRows();
    answer.resize( 27, nRows * 2 );
    answer.zero();

    // assume
    double N, N_rho, N_xi, N_rho_rho, N_xi_xi, N_rho_xi, N_xi_rho;
    FloatArray col1, col2;
    for ( int i = 1; i <= nRows; i++ ) {
        // basis and its derivatives
        N          = ders( i - 1, 0 );
        N_rho      = ders( i - 1, 1 );
        N_xi       = ders( i - 1, 2 );
        N_rho_rho  = ders( i - 1, 3 );
        N_xi_xi    = ders( i - 1, 4 );
        N_rho_xi   = ders( i - 1, 5 );
        N_xi_rho   = N_rho_xi;

        col1 = N_rho_rho * dFcyl_r_rho_rho_v
            + N_xi_xi * dFcyl_r_xi_xi_v
            + N_xi_rho * dFcyl_r_xi_rho_v
            + N_rho * dFcyl_r_rho_v
            + N_xi * dFcyl_r_xi_v
            + N * dFcyl_r_v;

        col2 = N_rho_rho * dFcyl_y_rho_rho_v
            + N_xi_xi * dFcyl_y_xi_xi_v
            + N_xi_rho * dFcyl_y_xi_rho_v
            + N_rho * dFcyl_y_rho_v
            + N_xi * dFcyl_y_xi_v;

        answer.setColumn( col1, 2 * i - 1 );
        answer.setColumn( col2, 2 * i );
    }
}

void AxisymStructuralElementEvaluator::computeGradientOfDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );
    elem->computeVectorOf( VM_Total, tStep, u );

    //// Get nodal coordinates
    FloatArray nodesCoords;
    elem->giveNodalCoordinates( nodesCoords );

    FloatArray xDef = nodesCoords + u;

    // Displacement gradient H = du/dX
    FloatMatrix G;
    this->computeGMatrixAt( gp, G );

    // get local code numbers corresponding to ir
    IntArray lc;
    FloatArray xr;
    this->giveIntegrationElementLocalCodeNumbers( lc, elem, gp->giveIntegrationRule() );
    xr.resize( G.giveNumberOfColumns() );
    for ( int i = 1; i <= lc.giveSize(); i++ ) {
        xr.at( i ) = xDef.at( lc.at( i ) );
    }

    answer.beProductOf( G, xr );
}

//void AxisymStructuralElementEvaluator::computeSecondOrderStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
//{
//
//    FloatArray vF, vG;
//    this->computeDeformationGradientVector( vF, gp, tStep );
//    this->computeGradientOfDeformationGradientVector( vG, gp, tStep );
//
//    answer = static_cast<StructuralCrossSection *>( this->giveElement()->giveCrossSection() )->giveSecondOrderStresses( vF, vG, gp, tStep );
//}
//
//
//void AxisymStructuralElementEvaluator ::computeConstitutiveMatrix_dAddF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
//{
//    //IGAElement *elem              = static_cast<IGAElement *>( this->giveElement() );
//    //SimpleSurfaceCrossSection *cs = static_cast<SimpleSurfaceCrossSection *>( elem->giveCrossSection() );
//    //FloatArray normal             = elem->giveNormal( gp );
//    //answer                        = static_cast<SimpleSurfaceCrossSection *>( this->giveElement()->giveCrossSection() )->giveSurfaceStiffnessMatrix_dAddF_3d( rMode, normal, gp, tStep );
//    answer = static_cast<StructuralCrossSection *>( this->giveElement()->giveCrossSection() )->giveStiffnessMatrix_dAddF_3d( rMode, gp, tStep );
//}
//
//void AxisymStructuralElementEvaluator ::computeConstitutiveMatrix_dAdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
//{
//    //IGAElement *elem              = static_cast<IGAElement *>( this->giveElement() );
//    //SimpleSurfaceCrossSection *cs = static_cast<SimpleSurfaceCrossSection *>( elem->giveCrossSection() );
//    //FloatArray normal             = elem->giveNormal( gp );
//    //answer                        = static_cast<SimpleSurfaceCrossSection *>( this->giveElement()->giveCrossSection() )->giveSurfaceStiffnessMatrix_dAdF_3d( rMode, normal, gp, tStep );
//    answer = static_cast<StructuralCrossSection *>( this->giveElement()->giveCrossSection() )->giveStiffnessMatrix_dAdF_3d( rMode, gp, tStep );
//}


} // end namespace oofem
