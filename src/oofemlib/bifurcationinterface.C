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

namespace oofem {

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
            int p             = 2;
            double alph_defl  = 1.;
            double dx_norm    = this->dx_Defl.computeNorm();
            double gamma      = p / ( dx_norm * dx_norm + alph_defl * pow( dx_norm, p + 2 ) );

            // apply Sherman-Morrison
            double vt_xold = this->dx_Defl.dotProduct( x );
            x              = x * ( 1 - gamma * vt_xold / ( 1 + gamma * vt_xold ) ); // Modify the classical solution

        } else { // Spectra eigenvector bifurcation
            FloatArray evectorFAmax;
            this->getEigenVectors().copyColumn( evectorFAmax, this->getEigenVectors().giveNumberOfColumns() ); // Chosen eigenvector in array

            // perturbate the solution by the eigenvector
            x   = this->alpha * evectorFAmax;
            this->setBifurcation( false ); // In the next iterations dont do bifurcation
        }
    }
}

} // namespace oofem