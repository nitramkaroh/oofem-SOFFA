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

//void NLStructuralSurfaceElement ::computeConstitutiveMatrix_dPdF_At( FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep )
//{
//
//    answer = this->giveStructuralCrossSection()->giveStiffnessMatrix_dPdF_3d( rMode, gp, tStep );
//    if ( this->matRotation ) {
//        FloatArray x, y, z;
//        FloatMatrix Q;
//        this->giveMaterialOrientationAt( x, y, z, gp->giveNaturalCoordinates() );
//        Q = {
//            { x( 0 ) * x( 0 ), x( 1 ) * x( 1 ), x( 2 ) * x( 2 ), x( 1 ) * x( 2 ), x( 0 ) * x( 2 ), x( 0 ) * x( 1 ), x( 2 ) * x( 1 ), x( 2 ) * x( 0 ), x( 1 ) * x( 0 ) },
//            { y( 0 ) * y( 0 ), y( 1 ) * y( 1 ), y( 2 ) * y( 2 ), y( 1 ) * y( 2 ), y( 0 ) * y( 2 ), y( 0 ) * y( 1 ), y( 2 ) * y( 1 ), y( 2 ) * y( 0 ), y( 1 ) * y( 0 ) },
//            { z( 0 ) * z( 0 ), z( 1 ) * z( 1 ), z( 2 ) * z( 2 ), z( 1 ) * z( 2 ), z( 0 ) * z( 2 ), z( 0 ) * z( 1 ), z( 2 ) * z( 1 ), z( 2 ) * z( 0 ), z( 1 ) * z( 0 ) },
//            { y( 0 ) * z( 0 ), y( 1 ) * z( 1 ), y( 2 ) * z( 2 ), y( 1 ) * z( 2 ), y( 0 ) * z( 2 ), y( 0 ) * z( 1 ), y( 2 ) * z( 1 ), y( 2 ) * z( 0 ), y( 1 ) * z( 0 ) },
//            { x( 0 ) * z( 0 ), x( 1 ) * z( 1 ), x( 2 ) * z( 2 ), x( 1 ) * z( 2 ), x( 0 ) * z( 2 ), x( 0 ) * z( 1 ), x( 2 ) * z( 1 ), x( 2 ) * z( 0 ), x( 1 ) * z( 0 ) },
//            { x( 0 ) * y( 0 ), x( 1 ) * y( 1 ), x( 2 ) * y( 2 ), x( 1 ) * y( 2 ), x( 0 ) * y( 2 ), x( 0 ) * y( 1 ), x( 2 ) * y( 1 ), x( 2 ) * y( 0 ), x( 1 ) * y( 0 ) },
//            { z( 0 ) * y( 0 ), z( 1 ) * y( 1 ), z( 2 ) * y( 2 ), z( 1 ) * y( 2 ), z( 0 ) * y( 2 ), z( 0 ) * y( 1 ), z( 2 ) * y( 1 ), z( 2 ) * y( 0 ), z( 1 ) * y( 0 ) },
//            { z( 0 ) * x( 0 ), z( 1 ) * x( 1 ), z( 2 ) * x( 2 ), z( 1 ) * x( 2 ), z( 0 ) * x( 2 ), z( 0 ) * x( 1 ), z( 2 ) * x( 1 ), z( 2 ) * x( 0 ), z( 1 ) * x( 0 ) },
//            { y( 0 ) * x( 0 ), y( 1 ) * x( 1 ), y( 2 ) * x( 2 ), y( 1 ) * x( 2 ), y( 0 ) * x( 2 ), y( 0 ) * x( 1 ), y( 2 ) * x( 1 ), y( 2 ) * x( 0 ), y( 1 ) * x( 0 ) },
//        };
//        answer.rotatedWith( Q, 't' );
//    }
//
//}



} // end namespace oofem
