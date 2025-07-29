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

} // end namespace oofem
