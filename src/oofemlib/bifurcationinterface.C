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

#include "bifurcationinterface.h"
#include "mathfem.h"
namespace oofem {

void BifurcationInterface::printPDStatus( FILE *outputStream ) const
{
    if ( !this->pdCheckPerformed ) {
        return; // no converged solution has been checked in this step
    }

    // The number of negative pivots of the LDL^T factorization equals the number of
    // negative eigenvalues by Sylvester's law of inertia. The smallest pivot is reported
    // as a pivot on purpose, it is not the smallest eigenvalue and does not scale like one.
    fprintf( outputStream, "Stability of converged solution : %s\n",
        this->pdLastResult ? "positive definite" : "INDEFINITE" );
    fprintf( outputStream, "  number of negative eigenvalues : %d\n", this->pdNumNegPivots );
    fprintf( outputStream, "  smallest diagonal pivot of D   : %+11.5e\n", this->pdMinPivot );
    fprintf( outputStream, "  size of the checked matrix     : %d\n\n", this->pdMatrixSize );
}


void BifurcationInterface::setEigenValuesVectors( SparseMtrx&A, FloatArray &Xeigs )
{
    FloatArray evaluesFA;
    FloatMatrix evectorsFM;
    int tmp = this->computeEigenValuesVectors( A, evaluesFA, evectorsFM );
    this->storeEigenValuesVectors( evaluesFA, evectorsFM, Xeigs );
}

void BifurcationInterface::performBifurcation( SparseMtrx &Ae, FloatArray &b, FloatArray &x )
{
    if ( this->bifurcation ) { // if bifurcatuion should be performed
        if ( this->choleskyBif ) { // using cholesky
            this->CholeskyUpdate( Ae, b, x ); // Update cholesky at the derived class

        } else if ( this->deflationBifurcation ) { // Do deflation analysis
            //int p             = 2;
            //int pp             = 0.5; // now set to member var
            double alph_defl  = 1.;
            double dx_norm    = this->dx_Defl.computeNorm();
            //double gamma      = this->p / ( dx_norm * dx_norm + alph_defl * pow( dx_norm, this->p + 2 ) );
            //double gamma =pp / ( dx_norm * dx_norm + alph_defl * pow( dx_norm,pp + 2 ) );

            // normalize gamma 
            double xd_norm = this->giveX0Defl().computeNorm();
            double xd_norm_p = pow( xd_norm, this->p );
            double gamma     = this->p / ( dx_norm * dx_norm + pow( dx_norm, this->p + 2 ) / xd_norm_p );

            // apply Sherman-Morrison
            double vt_xold = this->dx_Defl.dotProduct( x );
            x              = x * ( 1 - gamma * vt_xold / ( 1 + gamma * vt_xold ) ); // Modify the classical solution

            //// alternatively
            //double s = dx_norm * dx_norm;
            //double fact = ( s * s + s ) / ( s * s + s + 2. * vt_xold );
            //x           = x * fact;

        } else { // Spectra eigenvector bifurcation
            FloatArray evectorFAmax;
            this->getEigenVectors().copyColumn( evectorFAmax, this->getEigenVectors().giveNumberOfColumns() ); // Chosen eigenvector in array

            // perturbate the solution by the eigenvector
            x   = this->alpha * evectorFAmax;
            this->setBifurcation( false ); // In the next iterations dont do bifurcation
        }
    }
}

void BifurcationInterface::storeEigenValuesVectors( FloatArray &evaluesFA, FloatMatrix &evectorsFM, FloatArray &Xeigs )
{
    this->eigenvalues = evaluesFA;
    this->evectors    = evectorsFM;
    this->Xeigs       = Xeigs;

    // check size
    int nx  = Xeigs.giveSize();
    int nev = evectorsFM.giveNumberOfColumns();
    if ( evaluesFA.giveSize() != nx ) {
        //this->eigenvalues.resizeWithValues( nx );
        this->evectors.resizeWithData( nx, nev );
    }


}

} // namespace oofem
