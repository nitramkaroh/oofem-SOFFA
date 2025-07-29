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
#include "contextioerr.h"
#include "sm/Elements/nlstructuralelement.h"
#include "sm/Elements/nlstructuralsurfaceelement.h"
#include "sm/Elements/structural2dsurfaceelement.h"
#include "gausspoint.h"
#include "oofemlib/iga/iga.h"

namespace oofem {
StructuralSurfaceMaterialStatus ::StructuralSurfaceMaterialStatus( GaussPoint *g ) :
    StructuralMaterialStatus( g ), MVector(), tempMVector(), GVector(), tempGVector(), dIVector(), tempdIVector(), D2eDK2_Matrix(), tempD2eDK2_Matrix()
{
    FVector.resize( 9 );
    FVector.at( 1 ) = FVector.at( 2 ) = 1.; // True only if coordinate system with outward z is used

    if ( Structural2DSurfaceElement *el = dynamic_cast<Structural2DSurfaceElement *>( gp->giveElement() ) ) {
        FVector.zero();
        FVector.at( 1 ) = FVector.at( 3 ) = 1.; // For 2d surface element z is out of plane
    }

    if ( IGAElement *el = dynamic_cast<IGAElement *>( gp->giveElement() ) ) {
        // Construct 2d identity
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
    MVector.resize( 9 );
    tempMVector = MVector;

    GVector.resize( 27 );
    tempGVector = GVector;

    dIVector.resize( 27 );
    tempdIVector = dIVector;
    
    D2eDK2_Matrix.resize( 9, 9 );
    tempD2eDK2_Matrix = D2eDK2_Matrix;
}


void StructuralSurfaceMaterialStatus ::updateYourself( TimeStep *tStep )
// Performs end-of-step updates.
{
    StructuralMaterialStatus ::updateYourself( tStep );
    MVector      = tempMVector;
    GVector = tempGVector;
    dIVector    = tempdIVector;
    D2eDK2_Matrix = tempD2eDK2_Matrix;
}


void StructuralSurfaceMaterialStatus ::initTempStatus()
//
// initialize record at the begining of new load step
//
{
    StructuralMaterialStatus ::initTempStatus();

    // reset temp vars.
    tempMVector = MVector;
    tempGVector = GVector;
    tempdIVector = dIVector;
    tempD2eDK2_Matrix = D2eDK2_Matrix;
}


void StructuralSurfaceMaterialStatus ::copyStateVariables( const MaterialStatus &iStatus )
{
    StructuralMaterialStatus ::copyStateVariables(iStatus );

    const StructuralSurfaceMaterialStatus &structStatus = static_cast<const StructuralSurfaceMaterialStatus &>( iStatus );

    MVector     = structStatus.giveMVector();
    tempMVector = structStatus.giveTempMVector();

    GVector     = structStatus.giveGVector();
    tempGVector = structStatus.giveTempGVector();

    dIVector     = structStatus.givedIVector();
    tempdIVector = structStatus.giveTempdIVector();

    D2eDK2_Matrix     = structStatus.giveD2eDK2_Matrix();
    tempD2eDK2_Matrix = structStatus.giveTempD2eDK2_Matrix();
}

} // end namespace oofem
