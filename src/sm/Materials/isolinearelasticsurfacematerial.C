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

#include "linearelasticsurfacematerial.h"
#include "isolinearelasticsurfacematerial.h"
//#include "sm/CrossSections/simplecrosssection.h"
#include "sm/Materials/structuralms.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "contextioerr.h"

namespace oofem {
REGISTER_Material( IsotropicLinearElasticSurfaceMaterial );

IsotropicLinearElasticSurfaceMaterial ::IsotropicLinearElasticSurfaceMaterial( int n, Domain *d ) :
    LinearElasticSurfaceMaterial( n, d )
{
}

IsotropicLinearElasticSurfaceMaterial ::IsotropicLinearElasticSurfaceMaterial( int n, Domain *d,double _E, double _nu ) :
    LinearElasticSurfaceMaterial( n, d ),
    E( _E ),
    nu( _nu ),
    G( E / ( 2.0 * ( 1. + nu ) ) ),
    a( 0. )
{
    this->initTangents();
}


void IsotropicLinearElasticSurfaceMaterial ::initializeFrom( InputRecord &ir )
{
    LinearElasticSurfaceMaterial ::initializeFrom( ir );

    IR_GIVE_FIELD( ir, E, _IFT_IsotropicLinearElasticSurfaceMaterial_e );
    IR_GIVE_FIELD( ir, nu, _IFT_IsotropicLinearElasticSurfaceMaterial_n );
    IR_GIVE_FIELD( ir, a, _IFT_IsotropicLinearElasticSurfaceMaterial_talpha );
    //    propertyDictionary.add(tAlpha, a);
    // compute  value of shear modulus
    G = E / ( 2.0 * ( 1. + nu ) );
    this->initTangents();
}


void IsotropicLinearElasticSurfaceMaterial ::initTangents()
{
    double K      = E / ( 3.0 * ( 1. - 2. * nu ) );
    this->tangent = 2 * G * I_dev6 + K * I6_I6;
    this->computesSubTangents();
    this->alpha = {
        a, a, 0., 0., 0., 0.
    }; // CHANGED
    //this->alpha = {
    //    a, a, a, 0., 0., 0.
    //};
}

void IsotropicLinearElasticSurfaceMaterial ::giveInputRecord( DynamicInputRecord &input )
{
    this->LinearElasticSurfaceMaterial ::giveInputRecord( input );
    StructuralMaterial ::giveInputRecord( input );

    input.setField( this->E, _IFT_IsotropicLinearElasticSurfaceMaterial_e );
    input.setField( this->nu, _IFT_IsotropicLinearElasticSurfaceMaterial_n );
    input.setField( this->a, _IFT_IsotropicLinearElasticSurfaceMaterial_talpha );
}


void IsotropicLinearElasticSurfaceMaterial ::saveContext( DataStream &stream, ContextMode mode )
{
    LinearElasticSurfaceMaterial ::saveContext( stream, mode );

    if ( ( mode & CM_Definition ) ) {
        if ( !stream.write( E ) ) {
            THROW_CIOERR( CIO_IOERR );
        }
        if ( !stream.write( nu ) ) {
            THROW_CIOERR( CIO_IOERR );
        }
        if ( !stream.write( a ) ) {
            THROW_CIOERR( CIO_IOERR );
        }
    }
}


void IsotropicLinearElasticSurfaceMaterial ::restoreContext( DataStream &stream, ContextMode mode )
{
    LinearElasticSurfaceMaterial ::restoreContext( stream, mode );

    if ( mode & CM_Definition ) {
        if ( !stream.read( E ) ) {
            THROW_CIOERR( CIO_IOERR );
        }
        if ( !stream.read( nu ) ) {
            THROW_CIOERR( CIO_IOERR );
        }
        if ( !stream.read( a ) ) {
            THROW_CIOERR( CIO_IOERR );
        }
    }
    G = E / ( 2.0 * ( 1. + nu ) );
    this->initTangents();
}


double
IsotropicLinearElasticSurfaceMaterial ::give( int aProperty, GaussPoint *gp ) const
{
    if ( ( aProperty == NYxy ) || ( aProperty == NYxz ) || ( aProperty == NYyz ) ) {
        return nu;
    }

    if ( ( aProperty == 'G' ) || ( aProperty == Gyz ) || ( aProperty == Gxz ) || ( aProperty == Gxy ) ) {
        return G;
    }

    if ( ( aProperty == 'E' ) || ( aProperty == Ex ) || ( aProperty == Ey ) || ( aProperty == Ez ) ) {
        return E;
    }

    if ( ( aProperty == 'n' ) || ( aProperty == NYzx ) || ( aProperty == NYzy ) || ( aProperty == NYyx ) ) {
        return nu;
    }

    return this->Material ::give( aProperty, gp );
}


} // end namespace oofem
