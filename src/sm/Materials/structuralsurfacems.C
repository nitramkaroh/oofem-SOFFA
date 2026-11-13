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

#include "sm/Materials/structuralsurfacems.h"
#include "sm/Materials/structuralmaterial.h"
#include "gausspoint.h"
#include "iga/iga.h"
#include "tensor/tensor1.h"
#include "tensor/tensor2.h"

namespace oofem {
StructuralSurfaceMaterialStatus ::StructuralSurfaceMaterialStatus( GaussPoint *g ) :
    StructuralMaterialStatus( g )
{
    // Surface identity in the frame with the normal as the third axis, which all Lagrange
    // surface elements (NLStructuralSurfaceElement, 2D and 3D) use for F
    // (StructuralMaterialStatus has already set F = diag(1, 1, 1) for nlgeo 1)
    FVector.resize( 9 );
    FVector.zero();
    FVector.at( 1 ) = FVector.at( 2 ) = 1.;

    if ( IGAElement *el = dynamic_cast<IGAElement *>( gp->giveElement() ) ) {
        // IGA surfaces work in the global frame: I - N x N of the reference normal
        FloatArray normal = el->giveNormal( g );
        FloatArrayF<3> UDn( normal );
        Tensor1_3d Nt( UDn );
        Tensor2_3d id2D;
        FloatArrayF<9> indent3d = { 1., 1., 1., 0., 0., 0., 0., 0., 0. };
        Tensor2_3d Identity3d( indent3d );

        id2D( i_3, j_3 ) = Identity3d( i_3, j_3 ) - Nt( i_3 ) * Nt( j_3 );
        FVector          = id2D.to_voigt_form();
    }
    tempFVector = FVector;
    KVector.resize( 9 );
    tempKVector = KVector;

    // Empty gradients mark a fresh status. The surface evaluator supplies the reference
    // gradients before the first tangent; on a curved surface they are not zero.
}


void StructuralSurfaceMaterialStatus ::updateYourself( TimeStep *tStep )
// Performs end-of-step updates.
{
    StructuralMaterialStatus ::updateYourself( tStep );
    CVector = tempCVector;
    KVector = tempKVector;
    GVector = tempGVector;
    dIVector = tempdIVector;
}


void StructuralSurfaceMaterialStatus ::initTempStatus()
//
// initialize record at the begining of new load step
//
{
    StructuralMaterialStatus ::initTempStatus();

    // reset temp vars.
    tempCVector = CVector;
    tempKVector = KVector;
    tempGVector = GVector;
    tempdIVector = dIVector;
}


void StructuralSurfaceMaterialStatus::printOutputAt( FILE *file, TimeStep *tStep ) const
{
    // Surface elements store all nine components in their constitutive frame, including IGA
    // elements, which do not derive from NLStructuralElement.
    MaterialStatus::printOutputAt( file, tStep );
    fprintf( file, "  F " );
    for ( double value : FVector ) {
        fprintf( file, " %+.4e", value );
    }
    fprintf( file, "\n  P" );
    for ( double value : PVector ) {
        fprintf( file, " %+.4e", value );
    }
    fprintf( file, "\n" );
}


void StructuralSurfaceMaterialStatus ::copyStateVariables( const MaterialStatus &iStatus )
{
    StructuralMaterialStatus ::copyStateVariables(iStatus );

    const StructuralSurfaceMaterialStatus &structStatus = static_cast<const StructuralSurfaceMaterialStatus &>( iStatus );

    KVector     = structStatus.giveKVector();
    tempKVector = structStatus.giveTempKVector();

    GVector     = structStatus.giveGVector();
    tempGVector = structStatus.giveTempGVector();

    dIVector     = structStatus.givedIVector();
    tempdIVector = structStatus.giveTempdIVector();
}

} // end namespace oofem
