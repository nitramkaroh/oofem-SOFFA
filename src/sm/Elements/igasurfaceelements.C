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

#include "sm/Elements/igasurfaceelements.h"
//#include "sm/CrossSections/structuralcrosssection.h"
#include "sm/CrossSections/simplesurfacecrosssection.h"
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
#include "classfactory.h"

namespace oofem {
REGISTER_Element( BsplinePlaneStrainSurfaceElement );


BsplinePlaneStrainSurfaceElement ::BsplinePlaneStrainSurfaceElement( int n, Domain *aDomain ) :
    IGAElement( n, aDomain ), Space2dStructuralSurfaceElementEvaluator(), interpolation( 1 ) {
    this->setSecondGradient( true );
}


void BsplinePlaneStrainSurfaceElement ::initializeFrom( InputRecord &ir )
{
    Space2dStructuralSurfaceElementEvaluator::initializeFrom( ir );
    IGAElement ::initializeFrom( ir );
}


int BsplinePlaneStrainSurfaceElement ::checkConsistency()
{
    BSplineInterpolation *interpol = static_cast<BSplineInterpolation *>( this->giveInterpolation() );
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints( 1 ) ) {
        OOFEM_WARNING( "number of control points mismatch" );
        return 0;
    }
    return 1;
}

FloatArray BsplinePlaneStrainSurfaceElement::giveNormal( GaussPoint *gp ) const
{
    FloatArray normal;
    this->interpolation.boundaryEvalNormal( normal, 1, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    return normal;
};


void BsplinePlaneStrainSurfaceElement ::giveCompositeExportData( std::vector<VTKPiece> &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{
    
    int numSubCells = this->ExportIntegrationRulesArray.size();
    int numPointsInCell = this->ExportIntegrationRulesArray[0]->giveNumberOfIntegrationPoints();

    vtkPieces.resize( 1 );
    int numCells           = (numPointsInCell-1) * numSubCells;
    const int numCellNodes = 2; // linear line
    int nNodes             = numPointsInCell * numSubCells;
    //
    vtkPieces.at( 0 ).setNumberOfCells( numCells );
    vtkPieces.at( 0 ).setNumberOfNodes( nNodes );
    //
    IntArray nodes( numCellNodes );

    FloatArray nodeCoords( 2 );
    FloatMatrix Nmat;
    IntArray irlocnum;
    int numberOfIntegrationRules = numSubCells;
    int iNode = 1;
    // loop over individual integration rules
    for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        auto iRule = this->ExportIntegrationRulesArray[ir].get();
        for ( GaussPoint *gp : *iRule ) {
            this->computeNMatrixAt( Nmat,gp );

            FloatArray X_vert, X;
            this->giveNodalCoordinates( X_vert );

            // get positions of nonzero SFs
            IntArray lc;
            this->giveIntegrationElementLocalCodeNumbers( lc, this, gp->giveIntegrationRule() );
            X.resize( Nmat.giveNumberOfColumns() );
            for ( int i = 1; i <= lc.giveSize(); i++ ) {
                X.at( i ) = X_vert.at( lc.at( i ) );
            }

            // Compute coordinates of the points
            nodeCoords.beProductOf( Nmat, X );
            //nodeCoords.printYourself();

            vtkPieces.at( 0 ).setNodeCoords( iNode, nodeCoords );
            iNode++;
        }
    } // end loop over irules
    ///////////////////



    // write connectivity
    IntArray connectivity( 2 );
    int offset = 0;
    int iN, jN;
    int count = 0;
    for ( int iElement = 1; iElement <= numCells; iElement++ ) {
        iN                   = iElement + count;
        jN                   = iElement + 1 + count;
        connectivity.at( 1 ) = iN;
        connectivity.at( 2 ) = jN;
        //connectivity.printYourself();

        vtkPieces.at( 0 ).setConnectivity( iElement, connectivity );
        offset += 2;
        vtkPieces.at( 0 ).setOffset( iElement, offset );
        vtkPieces.at( 0 ).setCellType( iElement, 3 );
        if ( iElement % (numPointsInCell-1) == 0 ) {
            count++;
        }
    }
    
    // export the displacements
    int n = primaryVarsToExport.giveSize();
    vtkPieces[0].setNumberOfPrimaryVarsToExport( primaryVarsToExport, nNodes );
    for ( int i = 1; i <= n; i++ ) {
        UnknownType utype = (UnknownType)primaryVarsToExport.at( i );
        if ( utype == DisplacementVector ) {
            FloatMatrix Nmat;
            int numberOfIntegrationRules = numSubCells;
            int iNode                    = 1;
            // loop over individual integration rules
            for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
                auto iRule = this->ExportIntegrationRulesArray[ir].get();
                for ( GaussPoint *gp : *iRule ) {
                    this->computeNMatrixAt( Nmat, gp );

                    // Get nodal coordinates to the corresponding nodes
                    FloatArray u;
                    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );
                    elem->computeVectorOf( VM_Total, tStep, u );

                    // get positions of nonzero SFs
                    IntArray lc;
                    FloatArray ur, uval, uval_3d(3);
                    uval_3d.zero();
                    this->giveIntegrationElementLocalCodeNumbers( lc, this, gp->giveIntegrationRule() );

                    ur.resize( Nmat.giveNumberOfColumns() );
                    for ( int i = 1; i <= lc.giveSize(); i++ ) {
                        ur.at( i ) = u.at( lc.at( i ) );
                    }
                    uval.beProductOf( Nmat, ur );
                    uval_3d.at( 1 ) = uval.at( 1 );
                    uval_3d.at( 2 ) = uval.at( 2 );

                    vtkPieces.at( 0 ).setPrimaryVarInNode( utype, iNode, uval_3d );
                    iNode++;
                }
            } // end loop over irules
        }
    }
}


//////////////////
// axisym

REGISTER_Element( BsplineAxisymSurfaceElement );


BsplineAxisymSurfaceElement ::BsplineAxisymSurfaceElement( int n, Domain *aDomain ) :
    IGAElement( n, aDomain ), AxisymSurfaceElementEvaluator(), interpolation( 1 )
{
    this->setSecondGradient( true );
}


void BsplineAxisymSurfaceElement ::initializeFrom( InputRecord &ir )
{
    AxisymSurfaceElementEvaluator::initializeFrom( ir );
    IGAElement ::initializeFrom( ir );
}


int BsplineAxisymSurfaceElement ::checkConsistency()
{
    BSplineInterpolation *interpol = static_cast<BSplineInterpolation *>( this->giveInterpolation() );
    if ( giveNumberOfDofManagers() != interpol->giveNumberOfControlPoints( 1 ) ) {
        OOFEM_WARNING( "number of control points mismatch" );
        return 0;
    }
    return 1;
}

FloatArray BsplineAxisymSurfaceElement::giveNormal( GaussPoint *gp ) const
{
    FloatArray normal;
    this->interpolation.boundaryEvalNormal( normal, 1, gp->giveNaturalCoordinates(), FEIIGAElementGeometryWrapper( gp->giveElement(), gp->giveIntegrationRule()->giveKnotSpan() ) );
    return normal;
};


void BsplineAxisymSurfaceElement ::giveCompositeExportData( std::vector<VTKPiece> &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{

    int numSubCells     = this->ExportIntegrationRulesArray.size();
    int numPointsInCell = this->ExportIntegrationRulesArray[0]->giveNumberOfIntegrationPoints();

    vtkPieces.resize( 1 );
    int numCells           = ( numPointsInCell - 1 ) * numSubCells;
    const int numCellNodes = 2; // linear line
    int nNodes             = numPointsInCell * numSubCells;
    //
    vtkPieces.at( 0 ).setNumberOfCells( numCells );
    vtkPieces.at( 0 ).setNumberOfNodes( nNodes );
    //
    IntArray nodes( numCellNodes );

    FloatArray nodeCoords( 2 );
    FloatMatrix Nmat;
    IntArray irlocnum;
    int numberOfIntegrationRules = numSubCells;
    int iNode                    = 1;
    // loop over individual integration rules
    for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
        auto iRule = this->ExportIntegrationRulesArray[ir].get();
        for ( GaussPoint *gp : *iRule ) {
            this->computeNMatrixAt( Nmat, gp );

            FloatArray X_vert, X;
            this->giveNodalCoordinates( X_vert );

            // get positions of nonzero SFs
            IntArray lc;
            this->giveIntegrationElementLocalCodeNumbers( lc, this, gp->giveIntegrationRule() );
            X.resize( Nmat.giveNumberOfColumns() );
            for ( int i = 1; i <= lc.giveSize(); i++ ) {
                X.at( i ) = X_vert.at( lc.at( i ) );
            }

            // Compute coordinates of the points
            nodeCoords.beProductOf( Nmat, X );
            // nodeCoords.printYourself();

            vtkPieces.at( 0 ).setNodeCoords( iNode, nodeCoords );
            iNode++;
        }
    } // end loop over irules
    ///////////////////


    // write connectivity
    IntArray connectivity( 2 );
    int offset = 0;
    int iN, jN;
    int count = 0;
    for ( int iElement = 1; iElement <= numCells; iElement++ ) {
        iN                   = iElement + count;
        jN                   = iElement + 1 + count;
        connectivity.at( 1 ) = iN;
        connectivity.at( 2 ) = jN;
        // connectivity.printYourself();

        vtkPieces.at( 0 ).setConnectivity( iElement, connectivity );
        offset += 2;
        vtkPieces.at( 0 ).setOffset( iElement, offset );
        vtkPieces.at( 0 ).setCellType( iElement, 3 );
        if ( iElement % ( numPointsInCell - 1 ) == 0 ) {
            count++;
        }
    }

    // export the displacements
    int n = primaryVarsToExport.giveSize();
    vtkPieces[0].setNumberOfPrimaryVarsToExport( primaryVarsToExport, nNodes );
    for ( int i = 1; i <= n; i++ ) {
        UnknownType utype = (UnknownType)primaryVarsToExport.at( i );
        if ( utype == DisplacementVector ) {
            FloatMatrix Nmat;
            int numberOfIntegrationRules = numSubCells;
            int iNode                    = 1;
            // loop over individual integration rules
            for ( int ir = 0; ir < numberOfIntegrationRules; ir++ ) {
                auto iRule = this->ExportIntegrationRulesArray[ir].get();
                for ( GaussPoint *gp : *iRule ) {
                    this->computeNMatrixAt( Nmat, gp );

                    // Get nodal coordinates to the corresponding nodes
                    FloatArray u;
                    IGAElement *elem = static_cast<IGAElement *>( this->giveElement() );
                    elem->computeVectorOf( VM_Total, tStep, u );

                    // get positions of nonzero SFs
                    IntArray lc;
                    FloatArray ur, uval, uval_3d( 3 );
                    uval_3d.zero();
                    this->giveIntegrationElementLocalCodeNumbers( lc, this, gp->giveIntegrationRule() );

                    ur.resize( Nmat.giveNumberOfColumns() );
                    for ( int i = 1; i <= lc.giveSize(); i++ ) {
                        ur.at( i ) = u.at( lc.at( i ) );
                    }
                    uval.beProductOf( Nmat, ur );
                    uval_3d.at( 1 ) = uval.at( 1 );
                    uval_3d.at( 2 ) = uval.at( 2 );

                    vtkPieces.at( 0 ).setPrimaryVarInNode( utype, iNode, uval_3d );
                    iNode++;
                }
            } // end loop over irules
        }
    }
}

} // end namespace oofem
