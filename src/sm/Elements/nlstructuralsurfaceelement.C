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

#include "sm/Elements/nlstructuralsurfaceelement.h"
#include "sm/Materials/structuralms.h"
#include "sm/CrossSections/structuralcrosssection.h"
#include "feinterpol.h"
#include "domain.h"
#include "material.h"
#include "sm/CrossSections/simplesurfacecrosssection.h"
#include "integrationrule.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "dynamicinputrecord.h"
#include "gausspoint.h"
#include "engngm.h"
#include "mathfem.h"

namespace oofem {
NLStructuralSurfaceElement::NLStructuralSurfaceElement( int n, Domain *aDomain ) :
    NLStructuralElement( n, aDomain )
// Constructor. Creates an element with number n, belonging to aDomain.
{}

void NLStructuralSurfaceElement::computeFirstPKStressVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );

    FloatArray normal = this->giveNormal( gp );
    FloatArray vF;
    this->computeDeformationGradientVector( vF, gp, tStep );
    answer = cs->giveFirstPKSurfaceStresses( vF, normal, gp, tStep );
}


void NLStructuralSurfaceElement ::computeStressVector( FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep )
{
    auto cs = dynamic_cast<SimpleSurfaceCrossSection *>( this->giveStructuralCrossSection() );
    FloatArray normal = this->giveNormal( gp );   
    answer = cs->giveRealSurfaceStresses( strain, normal, gp, tStep );
}

void NLStructuralSurfaceElement::computeDeformationGradientVector( FloatArray &answer, GaussPoint *gp, TimeStep *tStep )
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    FloatArray u;
    this->computeVectorOf( { D_u, D_v, D_w }, VM_Total, tStep, u ); // solution vector
    if ( initialDisplacements ) {
        u.subtract( *initialDisplacements );
    }

    // Displacement gradient H = du/dX
    FloatMatrix B;
    this->computeBHmatrixAt( gp, B );
    answer.beProductOf( B, u );

    // Deformation gradient F = H + I
    answer.at( 1 ) += 1.0;
    answer.at( 2 ) += 1.0;

}



} // end namespace oofem
